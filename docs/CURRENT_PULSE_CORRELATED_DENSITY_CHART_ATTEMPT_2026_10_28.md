# Correlated pulse density and exact-chart attempt — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new [correlated density kernel](../src/openai_ns_reconstruction/pulse_correlated_density_delta.py) computes `Delta S` and ordinary `Delta S_Z` directly from the frame jets and separately supplied physical-minus-frame jets. It does not form two large `S` values and subtract them. Fraction and mpmath input types are preserved; caller-supplied difference jets must themselves come from a valid shared source. The kernel does not establish that condition and does not turn interval caps into values.

With `dQ=Q_P-Q_F`, `dU=U_P-U_F`, and frame variables without a prefix, the source formula uses

```
dA = dQ + (1-delta)*Z*dQ_Z/2 + dQ_y
dB = dQ*(Q+Q_y) + Q*dQ_y + dQ^2/2 + dQ*dQ_y
dC = -2*Z*dQ + (1-Z^2)*dQ_Z - 2*Z*dQ_y
d(UC) = dU*C + U*dC + dU*dC
dS = (dA + d(UC))/(1-delta*Z^2) + dB.
```

The implementation differentiates every term in ordinary Z, including `Q_ZZ`, `Q_yZ` and `L_Z=-2*delta*Z`; it requires the caller to declare `delta_Z=0` and rejects `L=0`. That declaration is a precondition, not a source certificate. It accepts point arithmetic only; a real interval containing zero in `L`, precision loss in upstream jets, source correlation, physical units and signed log-scale integration still require separate handling.

The frozen density fixture ran once after exclusive entry and its complete output was saved before review. All 13 synthetic controls passed: exact Fraction direct and independently differentiated references, separate quadratic/product components, swap sign, zero supplied difference, required assumptions, singular denominator rejection and a `10^-100` mpmath difference. For that synthetic difference, the 180-digit result differed from a separate 240-digit direct reference by approximately `2.81e-242` in absolute magnitude. This is a point-arithmetic check, not a certified physical error or an evaluated physical pulse. Result SHA-256: `054b4749ecfdbcd3b2501812e4befcd81022fdd279d9f0ac9365f418b496666f`.

## Exact chart attempt: failure preserved

A separate, frozen exact-coordinate draft represents the real source recipe `mu=exp(-4*(exp(40)+11))/1000` symbolically, keeps `13-4*mu` exact, and proposes the source transformation `xi=13+mu*s`, `s in [-2/mu,-4]`, `y=13/mu+s` and `dxi/mu=ds`. The end chart is `s in [-4,0]`. This form avoids rounding `13-4*mu` to 13 and does not materialize `exp(13/mu)`. It remains an **unaccepted draft**: its one fixture invocation failed at `gap_end_measure_map_same`, because the check required the literal label `dy = ds` for the gap while the draft emits `d_xi / mu = ds`. The full exception, traceback, frozen inputs and entered record are preserved; no result or chart certificate is accepted from that invocation. Do not re-enter the failed chart fixture or its entered source family. Failure-result SHA-256: `b9204f6d322bfaccfefacf950dc0cacc55c5018f4db6fcdb336da26863a97eb5`.

Existing saved own-pulse internal joins are within one original/native source; they are not a proof of equality between the newer physical and frame-only sources. The common `mu`, `Pstar`, coordinate, measure, full field jets, finite endpoint source owners and Z1 still need an actual cross-source connection. In particular, nonzero gap mean and both signed sides must remain. The previous M/Z0 failure and all error packets remain unchanged; no pulse correction, old error, or reference offset is cancelled. The real `mu` is too small for the float/SciPy pulse kernel, and the new density kernel alone does not solve nested exponential scaling or integrate actual fields.

Retained: 212 cumulative / 213 raw / 183 reduced terms, six signed channels, one reserved contribution, the unresolved independent reference offset and all old errors. No actual source pulse value, global five-moment closure, heat/stress correction, time recursion or Cartesian NS validation is claimed. Runtime model metadata was unavailable (`null`).

Next: construct a new, non-replaying source-owned gap/end coordinate and finite-endpoint witness, then bind real frame and physical jets plus their same-source difference to the correlated kernel. Preserve structured `13-4mu`, log radii and separate linear/quadratic scales; implement arbitrary-precision/log-scaled integration and certified error transport before judging the terminal correction. An old error or signed pulse term may be removed only with its corresponding same-source identity.

The report `experiments/root_st073/lei_ren_part1_paper_compliant_current_pulse_correlated_density_chart_attempt.json` embeds the exact code, fixture, runner, source, capsule, entry and result bytes for both leaves.
