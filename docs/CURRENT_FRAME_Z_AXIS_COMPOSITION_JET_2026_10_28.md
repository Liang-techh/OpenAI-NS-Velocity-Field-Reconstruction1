# Saved Z-axis substitutions and conditional C2 chain rule — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new source-and-symbolic leaf was frozen and entered once. It authenticated the raw Rp archive before decoding and parsing it, and authenticated the previous report/result; after parsing it selected the two saved `Z` substitution rows and their immediate variable and value descriptors. It did not traverse either expression subtree, value subtree or any protected geometry closure, invoke an ancestor, evaluate a source function or integral, or alter a numerical error.

The previously recorded 52 binder stacks contain exactly two `Z` substitution node IDs. Their saved definitions are:

| Substitution node | Variable | Value | Opaque expression reference |
| --- | --- | --- | --- |
| `31776` | node `20`, `bound_variable("Z")` | node `0`, exact rational `0/1` | `31764` |
| `31944` | node `20`, `bound_variable("Z")` | node `0`, exact rational `0/1` | `31725` |

Both source rows also carry `Z_independent_substitution=true`. The leaf records the value descriptors without an explicit equality assertion for zero; the table reads their `0/1` fields from the saved output. The exact zero is a source definition, not an interval midpoint or evaluated solution. Under the source's substitution-as-composition interpretation these are axis evaluations at `Z=0`; their argument has zero first and second ordinary Z derivatives. This is consistent with deliberate composition of an inlet function before axis evaluation, but the descriptors do not by themselves prove the provenance or correctness of the complete composed expressions, or that opaque source atoms have no other Z dependence.

For an ordinary C2 function `F(z,q)` and argument `a(z)`, the leaf used arbitrary quadratic Taylor jets to check exactly

```
d/dz F(z,a(z)) = F_z + F_q a'
d²/dz² F(z,a(z)) = F_zz + 2 F_zq a' + F_qq (a')² + F_q a''.
```

All five recorded symbolic checks passed, including the fixed-argument specializations. The computed checks use arbitrary C2 Taylor coefficients at `(z,q)=(0,a0)`, where `a1=a'(0)` and `a2=a''(0)`; the displayed formulas hold at any point under the stated C2 assumptions, with partial derivatives evaluated at the same `(z,a(z))`. In the fixed-zero case the terms involving `a'` and `a''` vanish, leaving the explicit `F_z` and `F_zz` terms. These remaining terms vanish only if there is no other explicit outer-Z dependence; that condition was not established for the actual saved source. This generic jet algebra is not an independent differentiation of the archived source graph.

The saved result deliberately leaves `actual_argument_derivatives_proved=false`: it records the rational descriptor and generic specialization, but does not contain a source-graph derivative proof. Likewise it does not prove that the saved replacement jets are genuine derivatives, that inlet rewriting commutes with either axis substitution, that all relevant binders are capture-free, or that a frame/m15 endpoint identity holds. The prior replacement free-name audit still has 1,790 unresolved states. No numerical source value, certified bound, global five-moment closure or later physical layer is claimed. The independent reference, once-only reserved contribution, nonzero histories and old errors remain; `P1` and `m15` remain distinct. Runtime model metadata was unavailable (`null`).

Next: trace the two saved axis-evaluation bodies and the `31655 -> 265772` inlet jet ports through an authenticated, bounded composition contract, distinguishing explicit outer-Z dependence from the evaluated Z argument. Prove the corresponding Z1/Z2 transport for the actual source functions, with opaque operations and protected geometry left as explicit gates; then return to the frame endpoint identity.
