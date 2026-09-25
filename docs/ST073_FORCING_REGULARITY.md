# ST073 forcing near the critical time

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) permits an external force, but Theorem 1.1 requires it to be smooth and compactly supported in space-time. Defining the force as the current momentum residual on a truncated time slab is only an identity there; it does not establish a smooth extension through the critical time. The paper addresses this in Sections 9–10: its improved residual is flat to arbitrary order near the singularity (equation 9.20), and the localized force extends smoothly (equation 10.5).

The current robust ST073 candidate remains unaccepted. `forcing_extension_screen.py` evaluates the full physical momentum residual at eight relative axial-collar points per time. These are selected points, not a spatial maximum or L2 norm. The registered time interval starts at backward time `tau = 0.5/64 = 0.0078125`.

| Backward time `tau` | Largest residual at sampled collar points |
| ---: | ---: |
| 0.024 | 5.0063e4 |
| 0.016 | 1.1619e5 |
| 0.012 | 2.1741e5 |
| 0.010 | 3.2658e5 |
| 0.0084 | 4.8484e5 |

The first five registered points fit a log-log slope of about `-2.16`. This is evidence that the present residual grows over the registered interval, not an asymptotic proof. The numerical report is `experiments/root_st073/compact_potential/forcing_extension_screen.json`.

With the explicit `--extrapolate` option, the diagnostic bypasses the candidate's finite-slab gate and uses central time differences. This is **outside the registered model domain**: the old fit coefficients and field formulas have not been validated at these times. It gives residuals `2.4390e6` at `tau=0.0042`, `1.2835e7` at `0.0021`, `6.8758e7` at `0.00105`, and `3.6441e8` at `0.000525`. The first five extrapolated points fit a log-log slope of about `-2.39`. At the smallest time and largest sampled point, halving the finite-difference step twice changes the residual norm by less than `1e-8` relative; that excludes an obvious step-size artifact at that one point, but it does not validate the extrapolated field. See `forcing_extension_extrapolation.json` and `forcing_extrapolation_fd_check.json` in the same report directory.

The next construction task is to derive and cancel the leading collar residual across scales, then build a time-valid inner/transition/outer match and verify smooth compact forcing. Retain the existing finite-slab gate by default. The experimental bypass is a diagnostic only. The required full momentum max and physical-volume L2 thresholds remain `1e-3`; neither is achieved by this screen.
