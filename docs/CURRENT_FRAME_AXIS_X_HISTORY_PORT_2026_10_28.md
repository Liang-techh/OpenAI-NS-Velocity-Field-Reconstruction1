# Axis evaluation and saved X-history ports — 2026-10-28

Goal ACTIVE / INCOMPLETE. Two new independently frozen leaves were each entered once. They authenticated the raw Rp archive before decoding; the second also authenticated the raw Gamma archive and the preceding saved results. They read only saved function-binding records, port lists and the previously recorded axis descriptors. They did not traverse expression edges or a protected geometry closure, call an ancestor builder, evaluate source functions or integrals, or remove old errors.

The observed replacement `31655 -> 265772` is the **value port of native angular history `X`** in both saved pulse inputs. Its first four saved rows have the following direct inlet-map correspondence:

| Ordinary Z order | Original X reference | Current X reference |
| --- | ---: | ---: |
| 0 | 31655 | 265772 |
| 1 | 31678 | 265776 |
| 2 | 32564 | 265782 |
| 3 | 32590 | 265789 |

These are saved reference/port matches, not a proof that the actual current functions are ordinary derivatives of one another or that inlet rewriting commutes with differentiation. In the pressure-C3 input, the original fourth and fifth `X` references are the global zero node `0`, while the current frame's respective references are `265797` and `265806`; the direct saved mapping does not replace global constants. The original source `native_C5_functions` sets the two higher `X` rows to zero and records that they require an actual constant-native-X identity. This is a structural compatibility gate for later high-order work: a nonzero node ID does not prove its function is nonzero, and neither equality nor inequality of those functions has been established here. The separate C3 input only supplies orders 0–3.

The authenticated Gamma binding stores `Xt=[31764,31765]` and `Xf=[31725,31731]`. The prior axis receipt records `31776` as substitution of `Z=0` in expression `31764`, and `31944` as substitution of `Z=0` in expression `31725`; hence the saved history value references agree with the corresponding axis-expression references. The Rp input record pins that Gamma archive hash, and the source-family metadata matches. The leaf deliberately leaves the stronger common-namespace semantics and actual composition identity unproved. Gamma's saved imported X C1 handles are `[31655,31678]`, matching the original `X` value and first-derivative port references. In source definitions, `Xt(0)` feeds the waiting scalar and `Xf(0)-Xf(Z)` enters the future angular control; this is not a new m15 or pressure segment.

Seven exact symbolic checks use arbitrary quadratic coefficients to distinguish axis evaluation `E0` from ordinary differentiation: `E0(f')=f'(0)` and `E0(f'')=f''(0)` need not vanish, whereas differentiating the fixed value `E0(f)` with respect to the evaluated variable gives zero. A separate polynomial example retains the first and second derivatives of an outer parameter that still depends on Z. These checks illustrate the required distinction; they are not actual derivative evaluations of the source graph, nor do they prove interchange of an improper integral and differentiation.

No actual `X` jet transport, capture freedom, same-source frame endpoint identity, m15 cancellation, numerical value, certified error reduction, global five-moment closure or subsequent physical construction is claimed. Existing errors, nonzero histories, the independent reference and the once-only reserved contribution remain. `P1` and `m15` remain distinct. Runtime model metadata was unavailable (`null`).

Next: establish the same-source ordinary-Z jet contract for the current `X` port through order two and the actual `Xt(0)`/`Xf(0)` composition, including explicit outer parameters and Gamma differentiation. Independently resolve the fourth/fifth-row compatibility before relying on selected C5 output. Keep unknown source semantics and protected geometry as gates; do not infer functional nonzero from node IDs or remove errors.
