# Three-time derivative-cardinal patch collocation

Run `python experiments/root_st073/curl_wave_time_cardinal.py` from the repository root. It reads the three-stage evolution and preceding midpoint fit, then writes `experiments/root_st073/compact_potential/curl_wave_time_cardinal.json` with all fitted coefficients and direct residual screens.

The previous single midpoint time bubble reduced the residual there but shifted the error to neighboring times. This experiment fits the second Hermite interval at normalized times `s=0.25, 0.5, 0.75`. The compact potential time polynomials have zero value at all three fit times, zero value and first derivative at both element endpoints, and **cardinal time derivatives** at the fit times. Pressure polynomials are value-cardinal at those times and vanish with first derivative at the endpoints. Their derivative and pressure interpolation matrices are the identity to the reported numerical precision. Thus each fitted time derivative and pressure affects its intended time node without changing the velocity at any fit time or the existing element interfaces.

The full nonlinear momentum operator was recomputed directly on 16 held-out spatial nodes and eight angles per node. The temporal finite-difference step is `dt/64`, small relative to the `2.5e-6` element. Results for maximum residual:

| Normalized time | Original field | Base Hermite element | Three-time correction |
| ---: | ---: | ---: | ---: |
| `0.25` | `727,348` | `1,748,616` | `422,364` |
| `0.50` | `727,061` | `1,748,646` | `462,672` |
| `0.75` | `726,774` | `832,393` | `444,479` |

This is a genuine improvement at **all three fitted times**, and the compact exact-curl potential retains divergence freedom. The recorded endpoint coefficient mismatch is `3.1e-12`; the potential mismatch at fit times is `5.9e-15`.

The high-degree cardinal time basis oscillates between the fit times. On a four-node spatial subset not used for fitting, the maxima at `s=0.125, 0.375, 0.625, 0.875` are `2,685,385`, `2,861,865`, `2,415,328`, and `1,461,770`, whereas the original field is about `579,000` at each of those points. The candidate is therefore rejected. The next solver must include interstitial time samples in the optimization and control temporal polynomial norms, rather than enforcing exact collocation at a few nodes. Even the best recorded fitted-time residual remains roughly eight orders of magnitude above the `1e-3` target; there is no global maximum or physical-volume L2 bound.
