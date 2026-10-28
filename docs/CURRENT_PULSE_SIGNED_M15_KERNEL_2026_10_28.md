# Signed pulse m15 integration kernel — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new numerical kernel, `src/openai_ns_reconstruction/pulse_signed_m15_correction.py`, computes a signed physical-minus-frame pulse correction from caller-supplied physical field jets and positive radii. It retains each chart's actual oriented interval, Jacobian, upper and lower endpoint contributions, and ordinary Z1. It does not call an ancestor constructor or equate the two sources.

Given common chart query `x`, common source `mu` and `delta`, Z-independent radii/endpoints and the source law `dR/dy=R`, it computes each side's physical source

```
S=(Q+(1-delta)*Z*Q_Z/2+Q_y)/(1-delta*Z^2)
  +Q*(Q/2+Q_y)
  +U*(-2*Z*Q+(1-Z^2)*Q_Z-2*Z*Q_y)/(1-delta*Z^2)
m15 = (1/2)*integral R^2*S dy - [R*Q_y]_(lower)^(upper).
```

The implementation applies `dy=dxi/mu` on main/gap and `dy=ds` on end, differentiates the density analytically including `Q_ZZ`, `Q_yZ` and `L_Z=-2*delta*Z`, and forms physical-minus-frame differences before summing charts. The integration-by-parts expression requires the stated radius law. A boolean input declaration is not a proof of that law, common coordinates, common `mu`, endpoint identity, physical units or the `Pstar` source mapping. The kernel does not itself verify equal main/gap `mu` values or identify the current source callbacks. It rejects undeclared Z-independent radius or endpoints, missing radius law, nonpositive radii and unconverged `quad_vec` runs. `quad_vec` reports estimates, not certified error bounds.

The frozen independent synthetic fixture was invoked once, with its entry saved before invocation and its complete result saved before review. All 13 recorded controls passed. It uses two unequal nonzero smooth fields and positive exponential radii, distinct physical chart intervals and synthetic `mu=2.5`; this reverses the gap interval to exercise orientation and is not the paper parameter. A separate direct implementation of the original `N1p=-S/2+(Q_yy+Q_y)/R` integral checks each side on each chart, the signed pair and finite endpoint reduction. Ordinary Z1 agrees with an independent centered difference; the maximum absolute chart discrepancy was `6.096306606195867e-06` on this fixture. Equal-side exact zero and rejection of a missing radius-independence declaration also passed. These results are synthetic algebra/numerics checks, not source-function equality or an evaluated physical pulse. The fixture does not exercise a zero-width chart.

**Important numerical limit:** the real source describes `mu=(1/1000)*exp(-4*(exp(40)+11))`, whose magnitude underflows ordinary binary64. This initial float/SciPy kernel therefore cannot yet evaluate the actual source parameter as a nonzero float, and the synthetic result must not be substituted for it. A future source-bound path needs an appropriate logarithmic or arbitrary-precision representation, verified physical field callbacks and common-coordinate/scale premises. It must not materialize a source handle as a number.

The prior M/Z0 attempt remains a recorded failure; neither it nor the earlier generic M family was re-entered. Existing 212 cumulative / 213 raw / 183 reduced terms, six signed channels, both sides' incoming history, one reserved contribution, independent reference offset, all old errors and the `physical_pulse-frame_only_pulse` term remain. No new source value, equality, certified physical error, global five-moment closure, heat/stress correction or NS validation is claimed. Runtime model metadata was unavailable (`null`).

Next: design a source-bound arbitrary-precision or log-scaled version for the actual tiny `mu`, bind genuine physical/frame field jets, source-owned radii, both sides' physical chart coordinates and finite endpoints, and establish the same-source conditions for the signed correction. Keep ordinary Z1 and all endpoint terms; never repeat this entered fixture or any ancestor computation. Only a proved common-source correction and error-packet transport can support removing old errors.

The published report `experiments/root_st073/lei_ren_part1_paper_compliant_current_pulse_signed_m15_fixture.json` embeds raw frozen source/capsule/entry/result, kernel, fixture and runner bytes. Result SHA-256: `40c9f2eb2f6ce794134b0b41016c69febf3e5595bb108973955280576a99ca19`.
