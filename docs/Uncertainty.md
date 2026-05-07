## Uncertainty estimation for [[HYS14 — Correction Method V2]]

Yes, I agree that the uncertainty should be included, but I think it needs to be separated into several terms rather than reported as one global number.

In this workflow, the drift is measured at an **hourly cadence**, not daily. Each hourly $\Delta t$ estimate is obtained by median-stacking the approximately eight 30-minute cross-correlation windows whose midpoints fall within that UTC hour, and then picking the peak lag of that hourly stack. So each hourly drift value should be interpreted as representative of that hour, approximately centered on the measurement window, not as an instantaneous value at the exact start or end of the hour. The correction itself is then applied per sample by interpolating the segment-smoothed hourly $\Delta t(t)$ model, so the working assumption is that the clock drift varies smoothly between adjacent hourly estimates. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

For uncertainty, I would separate four pieces:

### 1. Cross-correlation lag-pick uncertainty

The native sampling is 8 Hz, so the discrete lag spacing is:

$$\Delta t_{\mathrm{sample}} = \frac{1}{8} = 0.125 \ \mathrm{s}$$

That is the **lag-pick sample interval**, not the full uncertainty. If the only issue were nearest-sample rounding, the RMS uncertainty would be:

$$\sigma_{\mathrm{quant}} = \frac{0.125}{\sqrt{12}} \approx 0.036 \ \mathrm{s}$$

The picker in `peak_lag_hourly.py` is nearest-sample on the envelope of CC², so the chronos hourly Δt fed to chronfix lives on the 0.125 s grid. On top of that, picks are also affected by cross-correlation coherence, SNR, and whether the picker locks onto the correct lobe.

Two independent ways to quantify this directly from the data:

**(a) Parabolic fit to the envelope of each hourly CC²** (already implemented in `chronos/scripts/uncertainty.py`). Fit a parabola through the three envelope samples around the picked peak; the peak-position uncertainty is

$$\sigma_{x_0} = \frac{\sigma_{\mathrm{noise}}}{|a|\sqrt{2}}$$

where $a$ is the parabolic curvature coefficient and $\sigma_{\mathrm{noise}}$ is the robust off-peak envelope MAD. Converted to seconds via $\Delta t_{\mathrm{lag}}$. For HYS12-HYS14 this gives **median σ_lag ≈ 0.002 s, p90 ≈ 0.005 s** — well below the 0.036 s nearest-sample floor, because high-SNR CCFs let the parabolic fit recover sub-sample information that the nearest-sample picker discarded.

**(b) Spread of 30-minute window lag picks within each hour** — for each hour, calculate the peak lag for each individual 30-minute CCF window, then compute

$$\sigma_{\mathrm{lag,hour}} = 1.4826 \times \mathrm{MAD}(\delta t_{\mathrm{30min}})$$

or a block-bootstrap uncertainty on the hourly median (the 30-min windows overlap by 75% and are not independent, so block bootstrap is the right move).

(b) is a useful cross-check on (a) and would also catch lobe-jump events that the parabolic fit can miss. It is not yet implemented.

### 2. Model / smoothing uncertainty

This is the uncertainty in how well the smoothed clock model represents the cleaned hourly lag estimates. This is where my earlier “few tenths of a second” came from, but I should have explained it more precisely.

The methods file says that the deployed model residual is:

$$r(t) = \Delta t_{\mathrm{filtered/raw}}(t) - \Delta t_{\mathrm{modeled}}(t)$$

and that this residual has median $\approx 0.000$ s and MAD $\leq 0.115$ s on every inter-trigger segment. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

Converting that **worst-case** MAD into a robust Gaussian-equivalent 1-sigma value:

$$\sigma_{\mathrm{model,worst}} \approx 1.4826 \times 0.115 = 0.170 \ \mathrm{s}$$

with an approximate 95 % interval $1.96 \times 0.170 \approx 0.33$ s. That is the upper bound; **typical** behaviour across the full HYS14 record is smaller. The diagnostic script (`chronos/scripts/uncertainty.py`, see §"Code" below) computes σ_model = 1.4826 × MAD per inter-trigger segment and broadcasts it to every hour in that segment. Across all 33 segments of the HYS14 record:

- **median σ_model ≈ 0.077 s**
- **p90 σ_model ≈ 0.108 s**
- **worst-segment σ_model ≈ 0.170 s** (matches the bound above)

The script also writes per-segment statistics (median, MAD, robust σ, drift total, drift rate, max residual). Adding RMS and 68/95/99 percentile residuals to that CSV would be a one-liner if needed.

### 3. Drift magnitude

I would **not** treat the magnitude of the drift itself as an uncertainty term.

A 50 s drift is not uncertain just because it is large. The uncertainty is how well we estimate the drift curve $\Delta t(t)$. The magnitude matters because it tells us how consequential the correction is and how large the signal is relative to the uncertainty.

So instead of saying “uncertainty from the magnitude of drift,” I would report a relative error metric per segment:

$$\mathrm{relative\ uncertainty} = \frac{\sigma_{\mathrm{model}}}{\Delta t_{\mathrm{range}}}$$

For example, the methods file mentions a long curved segment with 56.7 s total drift and residual MAD $\leq 0.115$ s. That means the robust 1-sigma model error is at most about 0.17 s, which is only about:

$$\frac{0.17}{56.7} \times 100 \approx 0.3\%$$

of the total drift amplitude. That is a much better way to discuss drift magnitude: not as an uncertainty source, but as a signal-to-error or correction-size metric. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

### 4. Non-linearity in the drift

The current correction does not assume each drift segment is linear. Stable segments are replaced by a robust median, while drifting segments are modeled with a 24-hour centered rolling robust median followed by a 6-hour moving average. That is specifically meant to track curved drift rather than force a straight-line fit. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

However, there is still uncertainty from unresolved curvature between hourly estimates. The current method assumes the drift is smooth between hourly samples and applies linear interpolation between the segment-smoothed hourly values. The methods file explicitly notes that this cannot resolve sub-hour events. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

Note that this is **not** the same as the σ_nonlin term that an earlier draft of the diagnostic script computed (deviation of the smoothed model from a whole-segment linear fit). That figure of ≈ 0.3 s median, ≈ 2.7 s p90 is the error a *worse* method would incur if it forced linearity over an entire inter-trigger segment, and it was correctly removed from σ_total. Chronfix tracks curvature on multi-hour timescales via the 24 h rolling-median + 6 h MA, so the only residual non-linearity uncertainty is the sub-hour piece.

To quantify the sub-hour piece, I think we should add a leave-one-out or even/odd-hour test:

1. Build the smoothed drift model using only even hours.
2. Predict the withheld odd-hour $\Delta t$ values.
3. Compute the residuals between predicted and observed odd-hour values.
4. Repeat in the opposite direction.
5. Report MAD, robust (\sigma), and 95th percentile prediction error.

That would directly estimate the uncertainty from unresolved non-linearity and interpolation.

### Trigger / resync intervals

Trigger intervals should be reported separately from normal drift uncertainty. The trigger time is only localized to approximately one hour because the drift estimates are hourly and are built from 30-minute windows. The method therefore treats trigger intervals as undefined, keeps them as NaN, and splits the corrected output at those boundaries rather than interpolating through the clock jump. ([GitHub](https://github.com/coszo-hub/chronfix/blob/main/examples/HYS14%20%E2%80%94%20Correction%20Method.md "chronfix/examples/HYS14 — Correction Method.md at main · coszo-hub/chronfix · GitHub"))

So the uncertainty near triggers is not “0.2 s.” It is a different category: the exact resync time is only known at roughly hour-level precision, and no continuous correction should be assigned across that interval.

---

### Code

**The correction itself does not change — only diagnostics/reporting.**

A standalone uncertainty diagnostic now exists at `src/chronos/scripts/uncertainty.py`. Run with:

```bash
python -m chronos.scripts.uncertainty --pair HYS12-HYS14 --target HYS14
```

It writes to `data/uncertainty/HYS14/`:

```text
sigma_lag_hourly.npy        # parabolic-fit CCF pick uncertainty (§1)
sigma_model_hourly.npy      # per-segment robust σ of (raw − modeled) (§2)
sigma_nonlin_hourly.npy     # whole-segment linear-fit deviation, reference only
sigma_total_hourly.npy      # √(σ_lag² + σ_model²) — the deployed-method uncertainty
segment_summary.csv         # per-segment stats (see fields below)
uncertainty.png             # 3-panel diagnostic
```

`segment_summary.csv` currently includes:

```text
segment, start_index, end_index, n_hours, n_finite,
sigma_model_mad, sigma_nonlin_max, drift_total, drift_rate_s_per_day
```

Recommended additions (one-liners on top of the existing residual array):

```text
median_residual_s, rms_residual_s,
p68_abs_residual_s, p95_abs_residual_s, p99_abs_residual_s,
max_abs_residual_s,
relative_sigma_percent = robust_sigma_s / drift_range_s * 100
```

Two stronger diagnostics that are **not yet implemented** and would meaningfully tighten the answer:

```text
lag_pick_uncertainty_from_30min_windows    # cross-check on σ_lag
leave_one_out_prediction_error_for_nonlinearity   # sub-hour smoothness test
```

### Headline numbers (HYS12-HYS14, full record)

| component | median | p90 | worst |
|---|---|---|---|
| σ_lag (parabolic CCF pick) | 0.002 s | 0.005 s | — |
| σ_model (raw − smoothed, per-segment robust σ) | 0.077 s | 0.108 s | ≈ 0.170 s (one segment) |
| **σ_total = √(σ_lag² + σ_model²)** | **0.077 s** | **0.108 s** | — |

For the worst long curved segment (126 days, 56.7 s of total drift), σ_model ≈ 0.170 s ⇒ relative uncertainty ≈ 0.30 %.
