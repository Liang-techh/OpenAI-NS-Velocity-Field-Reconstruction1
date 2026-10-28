# Compact profile function and M source-role port — 2026-10-28

Goal ACTIVE / INCOMPLETE. This turn adds an independently callable numerical compact profile adapter and separately selects the four previously deferred M source-role rows. Both leaves were frozen before entry and invoked once; neither reran an ancestor, traversed the source graph, or evaluated an original-source function or integral. The original user changes and old error ledger remain untouched.

## Numerical profile interface

`src/openai_ns_reconstruction/compact_profile_function_port.py` accepts a fixed dimensionless radius `x=r/Rc` in `[1,2]`, a caller-supplied positive original bump normalization `J0`, positive integer `N`, and ordinary Z0..Z5 derivative jets of the amplitude `A` and the three source controls `h2,h3,h4`. It computes the source-form profile

```
ell=log(2)/40;  c=(1/5,1/2,4/5) log(2)
b_i=beta((log(x)-c_i)/ell)/(ell*J0*x)
F=(h2*b_0+h3*b_1+h4*b_2)/N;  delta_E=A*F.
```

Product rows use ordinary binomial Leibniz coefficients; `x` is held fixed. `beta(t)=exp(-1/(1-t^2))` for `|t|<1`, zero outside. The implementation uses an isolated mpmath context at caller-selected precision (at least 50 decimal digits), accepts arbitrarily large Python integers without binary64 conversion, and does not obtain source inputs or J0 itself. Its numerical results are uncertified, and support classification away from the exact fixed endpoints is at working precision. It does not compute `delta_V`, transported correction histories, physical velocity, a moment target, an error bound or an all-Z source identity.

One newly frozen synthetic fixture returned with **19/19 checks** passing: independent direct values at each bump center and the expected one-active-bump pattern, ordinary-Z product rows at an interior point, both fixed endpoints, and a `10^400` integer N scaling check. Fixture values, including J0, are synthetic; no actual selected N, source input, normalization integral or physical field was evaluated. Exact beta/provider correspondence and capture-free source substitution remain separate obligations.

## M history references

A separate new source-only leaf reads the saved source-role publication and selects the named field before ordinary Z order. It returned four C3/C5 M endpoint, terminal m-history and amplitude role triples for Z0..Z3. These are saved integer references in distinct source namespaces, not evaluated values or a newly proven cross-packet binding. They do not establish either endpoint as a history/amplitude quotient.

The conditional algebra retains `delta_m=T_N,m-T_C,m` and `a=A_N-A_C`; if the denominators are nonzero and the same-source normalizations and source interpretation hold, the M gap is `(A_C*delta_m-a*T_C,m)/(A_C*A_N)`. The C5 history is the saved N-scaled terminal history; this leaf adds no N division. It neither assumes `a=0` nor sets incoming histories to zero. The failed joint-family `KeyError: 0` remains recorded, and its consumer was not reentered. The 20 previously defined signed m/h/k/e/p history gaps and the M numerator are distinct objects; the new selection does not prove their common-domain relation or cancel any historical error.

## Remaining gates

The previously logged 212 cumulative / 213 raw / 183 reduced terms, six signed channels, once-only reserved contribution, independent reference offset and four physical-minus-frame pulse correction roots remain. Counts are inherited, not re-audited here; P1 and M are not m15. The physical pulse/frame equality remains open, as do actual compact source-provider binding, full-Z endpoint and derivative identities, common physical chart and units, and certified signed-gap transport.

Static producer review distinguishes `J_sigma(x)=integral_0^x sigma(t)dt` from the physical Jacobian `dy/dx`; their names do not make them interchangeable. Source text traces `sigma` to a flat exact-switch node and `mu=exp(log(1/1000)-4*(exp(40)+11))`, but actual graph-to-function semantics, capture-free substitution and finite-parameter admission have not been proved by this turn.

Next: bind the selected M references to the previously defined signed-history and numerator packets on the same physical chart, including exact N normalization and ordinary Z derivatives, using a fresh frozen leaf; separately bind actual compact beta/provider semantics and the physical pulse correction's radius/measure/endpoints. Do not rerun entered leaves or delete old errors before same-source identities. Global five-moment closure, exact heat exterior, stress/remainder, true order-n recursion and independent Cartesian NS validation remain incomplete. Runtime model metadata was unavailable (`null`).
