# Correlated Gamma scale difference — 2026-10-29

Goal ACTIVE / INCOMPLETE. A new Gamma difference kernel derives an exact source-shaped expression for the change in `Dhat` when the inverse-tail scale changes. It factors out `DeltaS` before numerical evaluation instead of subtracting two Gamma expectations. The complete original Gamma source bytes, accepted receipt and raw archive are pinned and authenticated before decoding. Direct saved roles are `schedule.S=31794`, `functions.Dhat=[31833,31840]`, `functions.energy=[31915,31916]`; these are archive references, not numerical values or proof of P/F source equivalence. No old producer, source-value reducer, ancestor quadrature or graph traversal was called.

For `a=delta/2`, `d=1-Z²`, `k=2d exp(-t)`, `w(v)=exp(-v)v^(a+1)/Gamma(1+a)`, `S_theta=S_F+theta*DeltaS`, the exact uncomputed functions are

```
Delta_Dhat = DeltaS*K0
K0 = -k²(a+1) integral_theta integral_s integral_v
     w(v)*v*s*(1+s*k*v*S_theta)^(-a-2)

Delta_Dhat_Z = DeltaS*K1
K1 = 4Z exp(-t)(a+1)k integral_theta integral_v
     w(v)*v*(1+k*v*S_theta)^(-a-2).
```

Here theta,s run over `[0,1]` and v over `[0,infinity)`. They follow from the actual original `Dhat` integral using `f_a(x)=integral_0^1(1+s*x)^(-a-1) ds` and the fundamental theorem of calculus. No old interval cap or midpoint is selected. The Z1 formula requires `a`, `S_F`, and `S_P` to be Z-independent; otherwise the two branch-dependent `Dhat_S*S_Z` chain terms must be included. Positive branch, shared source parameters, P/F ownership, input identities and differentiation of the complete future integral are not proved by this leaf.

The point API evaluates the normalized kernels using an analytic theta average expressed with `log1p`/`expm1` relative functions and a direct positive integrand. It does not compute a P/F numerical difference quotient or request precision proportional to `-log10(abs(DeltaS))`. The output is a rounded estimate, not a certificate. For nonnegative branch scales, the analytic inequalities `|K0| <= k²(a+1)²(a+2)/2` and `|K1| <= 4|Z|exp(-t)k(a+1)²(a+2)` follow by a Gamma moment; the API's rounded bound endpoints are explicitly **not certified**. Neither inequalities nor estimates are used as actual source values or to delete an old error.

The frozen one-shot fixture compared moderate synthetic inputs against independent integrals and checked `DeltaS -> 0` and `Z=±1`. All 11 reported diagnostic errors were below `1e-9`; the largest was `3.46e-13` in the ordinary-Z finite difference. At `Z=±1` both difference kernels vanish while the original ordinary-Z boundary derivative is retained and nonzero in the fixture. Actual `mu`, `S_F`, `S_P` and the full future-energy difference were not numerically evaluated; true `S=exp(-logRtail)` may have a double-exponential scale that cannot be represented as an ordinary mpf. Runtime model metadata was unavailable (`null`).

The result directly authenticates archive/receipt hashes and source-family metadata but does not independently certify each archived reference's graph semantics or bind P/F owners. The normalized point gate is not application-level source admissibility. The non-directed mpmath bound endpoints are estimates. The runner's failure path would label numerical evaluation as true even if it failed before that stage; the saved successful run did perform synthetic numerical evaluation. Entry/result writes are exclusive but not fsync-protected, and the runner does not claim protection from concurrent file modification.

The 212/213/183 inherited entries, six signed channels, once-only reserved contribution, four physical pulse corrections, `(P-G)+(G-F)`, independent reference offset and all old errors remain. P1 is not m15. Full Gamma energy also contains `Dhat*(2pre-a*S*Dhat)`; its explicit S difference, angular repair, waiting and all other future components must be retained before constructing correlated `Delta_future` and `Delta_C`. No global five-moment closure, exact heat exterior, admissible stress, true temporal recursion or corrected Cartesian NS validation is claimed.

Next: source-bind P/F `S`, common `a`, waiting, `pre` and the full Gamma energy factors; carry the explicit S difference in the quadratic heat integrand, derive ordinary Z1 including any nonzero S_Z, and combine all five complete-future terms before using the selector's positive-root difference. Independently bind F's actual radius/endpoints and prove common coordinates and certified error transport. Do not reenter this leaf or older computations.

Frozen source SHA-256: `519462985650302c1710d21f41ef37c0fcd3474ecfbb6cb98c2f09dcfe0ca43e`. Complete result SHA-256: `d08d019bd7e1d7903e3c710b8113b094885e791d444b1ca027b534067722fb51`.
