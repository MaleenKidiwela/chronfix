#!/usr/bin/env python
"""Standalone uncertainty diagnostic for the HYS14 Δt clock-correction track.

Quantifies three sources of uncertainty in Δt(t) and combines them:

  1. σ_lag(t)        per-hour CCF peak-lag pick uncertainty, from a parabolic
                     fit to the envelope-of-squared-CC around the picked peak.
                     σ_x0 = σ_noise / (|a| · √2), converted to seconds via Δlag.
                     σ_noise is the robust off-peak envelope scatter (MAD).

  2. σ_model(t)      per-segment residual scatter of (raw cleaned Δt) vs
                     (deployed segment-smoothed model). One MAD per inter-
                     trigger segment, broadcast to every hour in the segment.
                     Captures how well the smoother fits the picker series.

  3. σ_nonlin(t)     per-hour |modeled(t) − linear_fit(t)| within each
                     segment. Reported for reference only — this is the
                     error a *worse* method would incur if it assumed drift
                     were linear over a whole inter-trigger segment.
                     Chronfix does not do that: it linearly interpolates
                     between hourly samples of the 24 h rolling-median +
                     6 h MA smoothed model, which tracks curvature. Excluded
                     from σ_total.

  σ_total(t) = √(σ_lag² + σ_model²)

Reads existing pipeline outputs only — does not modify the correction.

Inputs
    data/peak_lag_hourly/<pair>/cc_hourly.npy
    data/peak_lag_hourly/<pair>/peak_lag_hourly_global.npy
    data/peak_lag_hourly/<pair>/hour_times.npy
    data/ccf/<pair>/lags.npy
    data/clock_estimate/<target>/delta_t_hourly_clean.npy
    data/clock_estimate/<target>/delta_t_hourly_filtered_raw.npy
    data/clock_estimate/<target>/trigger_periods.csv

Outputs (data/uncertainty/<target>/)
    sigma_lag_hourly.npy
    sigma_model_hourly.npy
    sigma_nonlin_hourly.npy
    sigma_total_hourly.npy
    segment_summary.csv
    uncertainty.png

Usage
    python -m chronos.scripts.uncertainty \
        --pair HYS12-HYS14 --target HYS14
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import numpy as np
from scipy.signal import hilbert

CHRONOS_ROOT = Path("/home/seismic/chronos")
LOG = logging.getLogger("uncertainty")

# --- Off-peak window (in seconds) used to estimate envelope noise --------
# Anything farther than this from zero is treated as "noise" for σ_n.
OFFPEAK_LAG_S = 30.0


def parabolic_peak_uncertainty(
    cc_hourly: np.ndarray,
    peak_lag: np.ndarray,
    lags: np.ndarray,
) -> np.ndarray:
    """Per-hour σ_lag in seconds from parabolic fit + off-peak noise.

    For each hour h:
      env = |Hilbert(cc**2)|
      k   = argmin |lags - peak_lag[h]|
      a   = (env[k-1] - 2*env[k] + env[k+1]) / 2     (in units of env / idx²)
      σ_n = 1.4826 * MAD(env at |lag| > OFFPEAK_LAG_S)
      σ_x0 = σ_n / (|a| * √2)            (peak index uncertainty)
      σ_lag = Δlag * σ_x0                (seconds)

    Hours where peak is at the lag-axis edge, the curvature is non-positive,
    or env is all-zero/NaN return NaN.
    """
    n_hours, n_lags = cc_hourly.shape
    dlag = float(lags[1] - lags[0])
    offpeak_mask = np.abs(lags) > OFFPEAK_LAG_S

    sigma_lag = np.full(n_hours, np.nan, dtype=np.float64)

    # Process in chunks to keep memory bounded — Hilbert on full array works
    # but doubles memory. Chunk size tuned for ~ 100 MB working set.
    chunk = 1024
    for i0 in range(0, n_hours, chunk):
        i1 = min(i0 + chunk, n_hours)
        block = cc_hourly[i0:i1].astype(np.float64)
        sq = block ** 2
        env = np.abs(hilbert(sq, axis=-1))

        for j in range(i1 - i0):
            h = i0 + j
            pl = peak_lag[h]
            if not np.isfinite(pl):
                continue
            k = int(np.argmin(np.abs(lags - pl)))
            if k <= 0 or k >= n_lags - 1:
                continue
            e_m, e_0, e_p = env[j, k - 1], env[j, k], env[j, k + 1]
            curvature = 0.5 * (e_m - 2.0 * e_0 + e_p)  # negative at a max
            if not np.isfinite(curvature) or curvature >= 0:
                continue
            offpeak = env[j, offpeak_mask]
            offpeak = offpeak[np.isfinite(offpeak)]
            if offpeak.size < 32:
                continue
            sigma_n = 1.4826 * np.median(np.abs(offpeak - np.median(offpeak)))
            if sigma_n <= 0:
                continue
            # σ_x0 in index units, then convert to seconds.
            sigma_x0 = sigma_n / (abs(curvature) * np.sqrt(2.0))
            sigma_lag[h] = dlag * sigma_x0

    return sigma_lag


def load_trigger_indices(csv_path: Path, n_hours: int) -> np.ndarray:
    """Boolean mask of length n_hours: True where the hour falls in a trigger."""
    mask = np.zeros(n_hours, dtype=bool)
    if not csv_path.exists():
        return mask
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            i0 = int(float(row["start_index"]))
            i1 = int(float(row["end_index"]))
            mask[i0:i1 + 1] = True
    return mask


def segment_bounds(trigger_mask: np.ndarray) -> list[tuple[int, int]]:
    """Inter-trigger segments as (start, stop) half-open index ranges."""
    n = trigger_mask.size
    out: list[tuple[int, int]] = []
    i = 0
    while i < n:
        if trigger_mask[i]:
            i += 1
            continue
        j = i
        while j < n and not trigger_mask[j]:
            j += 1
        out.append((i, j))
        i = j
    return out


def per_segment_uncertainty(
    raw: np.ndarray,
    modeled: np.ndarray,
    segs: list[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """σ_model (segment MAD of raw-modeled) and σ_nonlin (|modeled - linear|).

    Returns (sigma_model, sigma_nonlin, summary_rows).
    """
    n = raw.size
    sigma_model = np.full(n, np.nan, dtype=np.float64)
    sigma_nonlin = np.full(n, np.nan, dtype=np.float64)
    rows: list[dict] = []

    for seg_idx, (a, b) in enumerate(segs):
        idx = np.arange(a, b)
        m = modeled[a:b]
        r = raw[a:b]
        finite_m = np.isfinite(m)
        finite_both = finite_m & np.isfinite(r)

        n_finite = int(finite_both.sum())
        if n_finite < 5:
            rows.append(dict(
                segment=seg_idx, start_index=a, end_index=b,
                n_hours=b - a, n_finite=n_finite,
                sigma_model_mad=np.nan, sigma_nonlin_max=np.nan,
                drift_total=np.nan, drift_rate_s_per_day=np.nan,
            ))
            continue

        # σ_model: robust MAD of (raw - modeled) residuals across the segment.
        resid = r[finite_both] - m[finite_both]
        seg_mad = 1.4826 * np.median(np.abs(resid - np.median(resid)))
        sigma_model[a:b] = seg_mad

        # σ_nonlin: deviation of the modeled curve from a linear fit through it.
        # x is hours (0..n-1), y is modeled Δt. Linear fit on finite-only.
        x_full = np.arange(b - a, dtype=np.float64)
        x = x_full[finite_m]
        y = m[finite_m]
        # Ordinary least-squares is fine here — modeled series is already
        # smooth and outlier-free.
        slope, intercept = np.polyfit(x, y, 1)
        linear = slope * x_full + intercept
        nonlin = np.full(b - a, np.nan, dtype=np.float64)
        nonlin[finite_m] = np.abs(m[finite_m] - linear[finite_m])
        sigma_nonlin[a:b] = nonlin

        rows.append(dict(
            segment=seg_idx,
            start_index=a,
            end_index=b,
            n_hours=b - a,
            n_finite=n_finite,
            sigma_model_mad=float(seg_mad),
            sigma_nonlin_max=float(np.nanmax(nonlin)) if np.any(finite_m) else np.nan,
            drift_total=float(y[-1] - y[0]) if y.size else np.nan,
            drift_rate_s_per_day=float(slope * 24.0),
        ))

    return sigma_model, sigma_nonlin, rows


def plot_diagnostic(
    out_path: Path,
    hour_times: np.ndarray,
    modeled: np.ndarray,
    sigma_lag: np.ndarray,
    sigma_model: np.ndarray,
    sigma_nonlin: np.ndarray,
    sigma_total: np.ndarray,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = hour_times.astype("datetime64[h]").astype("datetime64[s]").astype(object)

    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)

    ax = axes[0]
    ax.plot(t, modeled, lw=0.7, color="k", label="modeled Δt")
    ax.fill_between(t, modeled - sigma_total, modeled + sigma_total,
                    alpha=0.3, color="C0", label="±σ_total")
    ax.set_ylabel("Δt (s)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(t, sigma_lag, lw=0.5, color="C1", label="σ_lag (CCF pick)")
    ax.plot(t, sigma_model, lw=0.5, color="C2", label="σ_model (segment MAD)")
    ax.plot(t, sigma_nonlin, lw=0.4, color="C3", alpha=0.5,
            label="σ_nonlin (reference only — not in total)")
    ax.set_ylabel("component σ (s)")
    ax.set_yscale("log")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    ax.plot(t, sigma_total, lw=0.6, color="C0")
    ax.set_ylabel("σ_total (s)")
    ax.set_xlabel("UTC")
    ax.grid(alpha=0.3)

    fig.suptitle("HYS14 Δt uncertainty decomposition", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def run(pair: str, target: str) -> None:
    pair_dir = CHRONOS_ROOT / "data" / "peak_lag_hourly" / pair
    ccf_dir = CHRONOS_ROOT / "data" / "ccf" / pair
    clock_dir = CHRONOS_ROOT / "data" / "clock_estimate" / target
    out_dir = CHRONOS_ROOT / "data" / "uncertainty" / target
    out_dir.mkdir(parents=True, exist_ok=True)

    cc_hourly = np.load(pair_dir / "cc_hourly.npy")
    peak_lag = np.load(pair_dir / "peak_lag_hourly_global.npy")
    hour_times = np.load(pair_dir / "hour_times.npy")
    lags = np.load(ccf_dir / "lags.npy")

    raw = np.load(clock_dir / "delta_t_hourly_filtered_raw.npy")
    modeled = np.load(clock_dir / "delta_t_hourly_clean.npy")
    trig_csv = clock_dir / "trigger_periods.csv"

    n_hours = hour_times.size
    assert cc_hourly.shape[0] == n_hours == raw.size == modeled.size, (
        "input arrays disagree on hour count"
    )

    LOG.info("σ_lag: parabolic fit on %d hourly CCs", n_hours)
    sigma_lag = parabolic_peak_uncertainty(cc_hourly, peak_lag, lags)

    trigger_mask = load_trigger_indices(trig_csv, n_hours)
    segs = segment_bounds(trigger_mask)
    LOG.info("σ_model + σ_nonlin: %d inter-trigger segments", len(segs))
    sigma_model, sigma_nonlin, rows = per_segment_uncertainty(raw, modeled, segs)

    # σ_total reflects the actual chronfix method (linear interp between
    # smoothed hourly samples), so it combines only σ_lag and σ_model.
    # σ_nonlin is kept as a side metric — see module docstring.
    components = np.stack([sigma_lag, sigma_model], axis=0) ** 2
    sigma_total = np.sqrt(np.nansum(components, axis=0))
    all_nan = np.all(~np.isfinite(np.stack([sigma_lag, sigma_model])), axis=0)
    sigma_total[all_nan] = np.nan

    np.save(out_dir / "sigma_lag_hourly.npy", sigma_lag)
    np.save(out_dir / "sigma_model_hourly.npy", sigma_model)
    np.save(out_dir / "sigma_nonlin_hourly.npy", sigma_nonlin)
    np.save(out_dir / "sigma_total_hourly.npy", sigma_total)

    with open(out_dir / "segment_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    plot_diagnostic(out_dir / "uncertainty.png", hour_times, modeled,
                    sigma_lag, sigma_model, sigma_nonlin, sigma_total)

    def _q(x: np.ndarray, q: float) -> float:
        x = x[np.isfinite(x)]
        return float(np.quantile(x, q)) if x.size else float("nan")

    LOG.info(
        "median σ (s): lag=%.3f  model=%.3f  total=%.3f  [nonlin(ref)=%.3f]",
        _q(sigma_lag, 0.5), _q(sigma_model, 0.5),
        _q(sigma_total, 0.5), _q(sigma_nonlin, 0.5),
    )
    LOG.info(
        "p90    σ (s): lag=%.3f  model=%.3f  total=%.3f  [nonlin(ref)=%.3f]",
        _q(sigma_lag, 0.9), _q(sigma_model, 0.9),
        _q(sigma_total, 0.9), _q(sigma_nonlin, 0.9),
    )
    LOG.info("wrote %s", out_dir)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pair", default="HYS12-HYS14")
    p.add_argument("--target", default="HYS14")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    run(args.pair, args.target)


if __name__ == "__main__":
    main()
