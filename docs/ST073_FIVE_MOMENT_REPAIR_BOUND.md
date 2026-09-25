# Necessary radius for repairing the current remote patch

The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), equation (4.15), requires the cumulative kinetic/axial moment `S` and centrifugal-pressure moment `Cp` to join consistently. The present two-moment `[4,16]` correction changes both relative to the unchanged pure-heat exterior. Before implementing a five-coefficient Newton solve on a reserved outer interval, we can test whether any such correction could satisfy even these two moment targets.

Suppose a correction is confined to `[16,R]`, with the old heat exterior restored beyond `R`. The base has `U_base=0` there. Let `M0=integral_16^R E_base^2 dX` and `C0=integral_16^R E_base^2/(2X) dX`. Restoring the existing outgoing differences `delta S>0`, `delta Cp>0` requires

```text
integral E_new^2 dX >= M0 + 2 delta S,
integral E_new^2/(2X) dX = C0 - delta Cp.
```

The first inequality also allows any new axial flow, whose square only increases the required `E_new^2` mass. Because `X<=R` and `E_new^2>=0`, a necessary condition is

```text
2 R (C0 - delta Cp) >= M0 + 2 delta S.
```

Using the analytic heat-tail moments, this gives the following lower bounds at the three audited eta slices:

| `eta` | `delta S` | `delta Cp` | Minimum `R` from `Cp` alone | Minimum `R` from `S` and `Cp` |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 3.58165 | 0.00167078 | 82.64 | **9711.48** |
| 0.2 | 3.46342 | 0.00159568 | 69.50 | **7887.62** |
| 0.35 | 3.24686 | 0.00144302 | 52.56 | **5582.19** |

For `[16,64]`, even the available nonnegative `Cp` after repair is negative at `eta=0` and `.2`; at `.35` the combined mass inequality still fails by a wide margin. Thus no profile correction confined to that reserved interval can restore both moments while retaining the unchanged heat exterior. At `eta=0`, the lower bound `R=9711.48` corresponds to physical radius about `1.277` at `tau=.0084` and `4.986` at `tau=.128` for the current `nu=.01` similarity map. This scale is much larger than the compact interface region and exceeds the radius-2 support of the repository's separate canonical constrained delivery at the earlier time; these are different candidate protocols.

The bound is **necessary, not sufficient**. Matching `M`, `I`, and `J`, keeping `E>0`, obtaining a relaxed-cone margin, preserving smoothness, and meeting the full Navier–Stokes residual can only add requirements. It applies to repairing this specific two-moment field back to the unchanged heat exterior on `[16,R]`. It does not rule out redesigning the inner profile, heat exterior, pressure datum, or moment target together. That co-design is now the better route than attempting a five-bump solve on `[16,64]`.

Run `python experiments/root_st073/five_moment_repair_bound.py`; the JSON result is `experiments/root_st073/five_moment_repair_bound.json`, marked `accepted: false`.
