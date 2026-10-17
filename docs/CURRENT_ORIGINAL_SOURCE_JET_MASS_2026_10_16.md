# Original source jets, transported lower endpoint and interior pulse mass

Tasks 14/15 remain **IN PROGRESS**; the long-term goal remains **ACTIVE / INCOMPLETE**. This stage adds new ordinary source-jet arithmetic, three transported lower-dy jets through ordinary Z2, local lower endpoint/mean-source C0/C1 operators, seven bounded pulse masses and 19 signed raw46 main-component C0 evaluations. It does not replace any active m15 term or complete terminal five-moment matching.

## Source-jet and lower-endpoint implementation

The new kernel implements ordinary derivatives through order two for products, positive-denominator quotients, exponentials and transported lower-dy. Order-one and order-two inputs are explicit; missing derivatives and nonpositive source denominators raise. This is spatial source differentiation, not the order-dependent temporal n recursion.

At second-micro native coordinate 2 (graph query node 48), the original source formula is

```
lower_dy = actual.Phi / comparison.Phi[0] * comparison.Psi[1] * rho.
```

The implementation hydrates same-family saved second-micro Phi and first-collar core radial jets. Their saved coefficients are Taylor coefficients; both approximations and enclosures are multiplied by k! before entering the ordinary-derivative kernel. Core comparison continuation uses the fixed radial derivative index in its Cauchy norm, then the separate axial weight and factorial conversion. Positive continuation caps widen source enclosures only. They are never signed values. The inherited symmetric finite-N envelope is conservative input uncertainty, not a selected signed finite-N correction.

The actual source rho=4*exp(2hb) remains enclosed with the width factored. Approximation rho=4 is explicitly a numerical proxy. The rho remainder bound retains the factor s=2. Exact interval endpoints are copied between contexts; no macro J66987 factor is applied. Independent leading graph ports and the native coordinate/radius/J/phase bindings remain pinned.

Three Z queries (-.3, 0, .3) now have ordinary Z0/Z1/Z2 lower-dy approximations and directed enclosures. These feed signed local E=2Z*B/(1-delta Z^2) and C(B)=[(1-delta)Z*B+(1-Z^2)B_Z]/(1-delta Z^2) C0/C1 operators. The lower endpoint keeps the 4*epsilon^2*hbB physical reference factor; actual radius discrepancies, mean histories and raw24 remain active. This does not evaluate the whole macro viscous source integral.

## Original pulse mass and signed C0 component

The pinned original flat pulse is gp(t)=exp(-1/t^2) for t>0, zero for t<=0, and sigma(t)=gp(t)/(gp(t)+gp(1-t)). The new independent mass is

```
M(s)=integral_0^s (1-sigma(t)) dt.
```

Source symmetry gives exact M(0)=0 and M(1)=1/2. Interior displayed approximations use a fresh scalar adaptive quadrature of this explicit formula; they are not rectangle midpoints or exact source values. Independently computed monotone endpoint rectangles with 128 and 256 panels certify directed boxes. Since 1-sigma decreases, their analytic width bound is s/N. The finite quadrature nodes are machinery, not field parameters. Exact MPF endpoint tuples and serialized-approximation error tuples are authoritative; displayed decimal bounds are labels.

| s | Numerical mass approximation | Certified absolute error upper (display) |
| --- | --- | --- |
| 0 | 0.0 | 0.0 |
| 0.125 | 0.1249999999999999999999999999994363936234 | 1.5659015e-31 |
| 0.25 | 0.2499999953658073349167873386053756825758 | 3.321377e-10 |
| 0.5 | 0.4793661246122048874190013670773094388455 | 0.00049082446 |
| 0.75 | 0.4999999953658073349167873386053756825758 | 0.0014648428 |
| 0.875 | 0.4999999999999999999999999999994363936234 | 0.0017089844 |
| 1 | 0.5 | 0.0 |

The table rounds both columns for display. Certified error tuples apply to the full serialized JSON approximation, rather than a further-rounded table label.

Mass is independent of Z. Symmetry comparisons between numerical approximations and their error radii are labeled numerical consistency; the functional symmetry identity comes from the source formula. No interval overlap is used to claim a source equality.

At all 19 inherited first-switch queries, the new interface computes C(B)*M(s) from the saved signed B/B_Z approximations and their exact source enclosures. The physical factor remains R100*epsilon*hbS^2, with hbB and hbS kept distinct and R100 separate from actual R(s). Raw47 remains the weighted-primitive remainder. The result is the raw46 main component with a certified enclosure; it is not the full weighted source primitive.

All 16 independent saved C(B) C0 comparisons and 6 exact-mass endpoint raw46 comparisons pass fixed relative tolerance 1e-70, with no source-error radius or interval-overlap tolerance. Interior comparisons report differences between two numerical mass approximations; the maximum observed relative difference is 1.00482151238e-8. They do not claim exact equality with a legacy quadrature result or functional closure.

## Actual verification and preservation

New production exits 0 in 20.794 s. The saved report is 677,531 bytes, SHA256 `c3b2e2a974ebca402e2d2762ab2c508150fdd72f80aa2882152e944bdad84efe`, with 1958 input pins and 216 implementation pins. Serialized inspection checks pins and packet invariants without constructing an owner, fixture or ancestor. All 61 independent analytic controls pass; maximum supporting difference is 5.714936956411374911078918e-101.

The controls differentiate nonconstant exp, polynomial and sin functions independently, test directed transport enclosures and reject missing rows, Taylor/ordinary confusion and an extra macro J factor. The production runs only the new bounded kernels over hydrated data. Accepted ancestor constructors, producers, checkers, fixtures and calculations are not rerun.

The raw213/reduced183 ledger witnesses, six signed inner/lower contributions, independent negative complete-pulse m15 target with reference*A*Rv^2, reserved n1 own-swirl exactly once and all prior error budgets remain unchanged by parent SHA/digest. P1/m13 and m15 remain separate. No positive source error is deleted, centered or canceled by this source export.

## Remaining construction and next task

The first-switch factored B_ZZ source row and weighted primitive/remainder remain the immediate dependencies for signed source-integral C1. The local transported lower B_ZZ now exists as a numerical approximation/enclosure; it is a different source and cannot substitute for first-switch B_ZZ. Continue from saved same-source C5 leaves, retaining ordinary derivatives and physical factors, then compare coupled signed contributions with the independent negative full-pulse m15 target.

A local untracked `original_factored_source_jet.py` draft was produced separately. Its worker smoke computes enclosure-only modal B0/B1/B2 rows with approximation=None and retains a formal swirl scale. Root has not accepted or integrated it; it is excluded from this report and commit. Review its actual source bindings and scale algebra next, and do not treat its existence or interval bounds as recovered source point values.

Total signed m15, functional all-five/common-N terminal closure, exact heat and finite-energy tail, admissible stress and flat remainder, true temporal n>=2 recursion/summation, realizable mean/oscillatory corrections and independent corrected Cartesian NS remain incomplete. No leading-background residual is asserted to satisfy the final 1e-3 corrected-NS threshold.
