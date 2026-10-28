# Compact profiles and correction densities with formal N — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new numerical interface covers both compact correction profiles and their radial y derivatives while keeping N symbolic; its single frozen synthetic fixture passed 79 checks. No source graph was traversed, no original solver or ancestor was rerun, and no actual selected-N source input or original integral was evaluated.

`src/openai_ns_reconstruction/compact_scaled_profiles.py` takes caller-provided ordinary Z0..Z5 jets of the Rc amplitude A and five controls, a positive J0, and fixed `x=r/Rc` in `[1,2]`. It returns the numerical coefficient jets `X=N*delta_E` and `Y=N*delta_V`, plus their first and second `y=log(x)` derivatives. The bump and its radial derivatives follow the original compact source formula; support classification is at working precision and is not certified. The interface does not differentiate a composed `x(Z)`.

For caller-supplied background E, it separately exposes these formal source-density terms without materializing N or discarding the quadratic corrections:

```
m = Y/N;                     h = X/N
k = E*Y/N + X*Y/N^2
e = -E*X/N + (Y^2-X^2/2)/N^2
p = E*X/N + X^2/(2*N^2).
```

The inverse powers and their signs are exact algebraic bookkeeping; the coefficient jets are rounded numerical outputs, not certified values. The selected original integer is saved as `N=2^(2^1358356628656378313)` in the prior repair-band handoff. It is not constructed or replaced here, and its source admissibility is not checked by this interface. The synthetic fixture tests distinct controls at all three bump centers, both fixed endpoints, ordinary derivative products, radial derivatives, density coefficients at Z0/Z1, and rejection of a radius just outside `[1,2]` before rounding.

A separate static audit found no complete actual input packet for A, the five control jets and J0. Existing exact function graph references, finite-N upper bounds, diagnostic N candidates, sampled values and interval midpoints are not interchangeable with actual numerical source values. True source evaluation and provenance, J0, full-Z provider equality, error bounds and physical transport remain open.

For the m15 comparison, the original own-field producer explicitly defines gap radius `exp(logRp+xi/mu)` and the pressure producer constructs finite endpoint tensors from that field radius. The gap between original F endpoints and their post-inlet-substitution frame-only counterparts is still unbound. In particular, the separate current-geometry metadata radius cannot be substituted for F's field radius without a same-source proof. Existing physical pulse correction roots, the independent reference offset, the 212/213/183 ledger, six signed channels, once-only reserved term and all old error terms remain; P1 and M are distinct from m15.

Next: bind F's original field radius and endpoint tensors through the saved inlet substitution, including capture-free ordinary Z1 transport, physical measure, orientation and Pstar; retain `P-F=(P-G)+(G-F)` and do not drop any correction or error without the needed same-source identity. In parallel, obtain actual compact A/control/J0 source values and a rigorous error account before using these numerical interfaces for physical conclusions. Global five-moment closure, exact heat exterior, stress/remainder, order-n recursion and independent Cartesian NS validation remain incomplete. Runtime model metadata was not available (`null`).

## Input-guard correction

Static review found that the previously published numerical adapter converted radius to working precision before checking `[1,2]` or endpoint identity. That could accept an exact input just above 2 or label one just below 2 an exact endpoint. The adapter now checks the original comparable value before rounding, accepts Decimal and Fraction inputs, and labels near-support classification as working-precision uncertainty rather than source-certified zero. It also records that the original native route requires `N>=160`; a general positive integer may still be used for synthetic algebra, without claiming source admissibility. A separately frozen new guard fixture passed 7 checks, including `2+10^-100` rejected and `2-10^-100` not called the endpoint at 50 digits. These are synthetic checks, not a proof of actual field admissibility.

The prior version and its one-shot receipt remain recoverable byte-for-byte from commit `439357ac72bb229c5ff11ab8ee4e30dbe0803044` and `current_compact_profile_function_port.json`; its old fixture and runner were not reentered. Its frozen checkpoint refers to the historical source bytes, not the updated file. The new guard source/checkpoint/entry/result are embedded in `current_compact_profile_input_guards.json`.
