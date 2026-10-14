# Actual current pulse second radial derivatives and pressure FTC, 2026-10-14

Scope: actual current Rp inputs now define leading main/gap/end `U_yy` C5, `Q_yy` C4, `N1p` C3 and signed partial pressure primitives. This is a source-function milestone. Subsequent Rp/Rv radial-pressure source joins and canonical/current partial-provider equivalence are accepted in [CURRENT_ORIGINAL_N1_OUTER_RADIAL_PRESSURE_SOURCE_JOINS_2026_10_14.md](CURRENT_ORIGINAL_N1_OUTER_RADIAL_PRESSURE_SOURCE_JOINS_2026_10_14.md); tasks6..9 below are complete within that scope. Higher mixed matching, combined outer acceptance, numerical bounds and global n1 remain open. The full goal is **ACTIVE / INCOMPLETE**.

## Implemented construction

- Consume the accepted eight-function current Rp C5 callback and its actual selected `ap/C1/C2` functions. Preserve the finite original mu, gp defining integral, beta normalization, original end scale and nonzero gap mean.
- Main/gap physical-y differentiation is `mu*d/dxi`; end differentiation is `d/ds`. With `bp=1/2+mu`, construct `U_yy=Utheta*(B_yy-2*bp*B_y+bp^2*B)` through ordinary Z5. Field rows are already ordinary derivatives. Only sigma's second Taylor row needs the factor two.
- Use `m_yy=U_y-m_y` and the leading `1-delta` recovery to obtain `Q_yy` through Z4. Keep the distinct `1+delta` own recovery unchanged.
- Construct the original reduced density `P1_y=R*N1p=-R*S/2+Q_y+Q_yy`, with the complete signed linear, quadratic and axial terms.
- Build actual current-query partial primitives for all ten I0/I1 products on all three charts. The integral binder is the field query's name. Preserve the original chart domains and exactly one physical Jacobian: `R/mu`, `R^2/mu` on main/gap and `R`, `R^2` on end.
- Retain all complete preceding charts and the finite global lower endpoint tensors. Return signed pressure increments rather than resetting incoming own or absolute leading pressure.
- A local globals adapter runs the unchanged, pinned pressure-map body with the restored owner's current FunctionJet/Algebra classes. It changes no accepted source file and no derivative row. This fixes the restore-time Python class-identity failure.

The saved graph has **285239 nodes**. The acceptance contains three scoped proof groups:

| Evidence | Scope |
| --- | --- |
| 54 actual second-radial source rows | B_y, B_yy and U_yy through Z5 on three charts |
| 36 actual flat pulse rows | B/B_y/B_yy at main start Rp and end exit Rv |
| 120 actual integral FTC rows | Ten products, Z0..3, three charts, with actual named binders and measures |
| 36 actual mass recovery rows | m_y and m_yy through Z5 |
| 33 directly checked Q recovery rows | Existing Q/Q_y Z0..2 and new Q_yy Z0..4 |
| 12 existing higher Q rows | Consumed from pinned genuine ordinary C5/Q4 and actual Rp substitution receipts |
| 6 pressure FTC whole-function identities | m13_y=P1_y and m15_y=-R*P1_y on each chart |
| 18 higher pressure FTC rows | Follow by differentiating those identities with the pinned ordinary source algebra; not a separate expansion of every rational derivative |

No numerical ancestor replay or old 779-second Rp check is required. Numerical oracle availability is false. Source gates for this milestone are accepted; combined outer seam/accessor, temporal recursion and corrected Navier-Stokes gates remain false.

## Files and recovery

Source and checker:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check.py`
- Matching `.json.gz` source archive and `_check.json` acceptance receipt.

Inspect the pinned receipt, or reconstruct the checked source on any host:

```powershell
python scripts/run_original_pulse_radial_pressure_C3_functions.py
python scripts/run_original_pulse_radial_pressure_C3_functions.py --restore-source
```

The public constructor is `CurrentOriginalPulseRadialPressureC3Functions()`. `at_fraction('main', '1/2')` returns actual source expressions for fields, density, FTC and signed pressure increments. It supplies no interval-valued numerical callback and makes no global n1 claim.

## Agent queue

Take the first unchecked construction dependency; retain the accepted scopes and publish evidence with each milestone. Do not repeat broad historical validation.

The immediate Rp target is quiet `true_quiet_Rp_C5_C3_bindings.Rp_endpoint.physical.{U,U_y,m,m_y,Q,Q_y,Q_yy}` and top-level `N1p` evaluated at `s=length`. Recover quiet `P1_y` as `R*N1p`; the endpoint stores no separate P1_y. Compare against new `endpoints.main.left.density.{N1p,P1_y}`. The quiet linear-mass theorem retains nonzero `m=K*Z/R`, `Q=-K/R`, `Q_y=-Q`. The native eight-function frame lacks Q/Q_y/Q_yy/N1p, so the exact physical recovery/source projection is still required before pressure matching. Absolute leading pressure is a separate quantity.

1. [x] Construct actual gp/beta B_yy with original finite mu and normalization; keep the gp integral at the entrance.
2. [x] Construct U_yy C5, actual mass m_yy C5 and leading Q_yy C4.
3. [x] Construct leading N1p/P1_y C3 with all signed reduced-density terms.
4. [x] Construct actual current-query signed partial pressure primitives and physical-y FTC, preserving named binders, finite lower tensors and complete preceding charts.
5. [x] Save scoped source acceptance and a source-only recovery entry point.
6. [x] Bind Rp main-start Q/Q_y/Q_yy to the actual quiet Rp C5 frame. Use the accepted 48-row inlet and 146-row consumer substitution identity; retain the actual nonzero mean and the proved physical radius identity.
7. [x] Prove all four total own pressure-y Z0..3 joins at Rp: compose the same incoming own pressure history with the new fixed leading N1p. Direct source pressure matching and recorded FTC/endpoint composition are accepted in the newer scoped receipt.
8. [x] Prove all four total own pressure-y Z0..3 joins at Rv. Consume the actual flat axial pulse, empty backward leading histories and specified own tail; preserve absolute leading/swirl pressure and m13/m15 endpoint continuity.
9. [x] Compare the new current-input partial pressure provider with the existing canonical partial pressure body using the accepted current Rp source substitution and old/current radius identity. Actual integrands and finite tensors are bound before pressure outputs are compared. A checked combined outer accessor still requires tasks10..12.
10. [ ] Establish main/gap and gap/end leading radial-pressure/mixed joins with the true original domains and full gp/beta flatness. Keep the nonzero gap mean.
11. [ ] Determine and construct the required additional physical-y/Z mixed rows for the complete outer own seam theorem. Newly available U_yy/Q_yy does not itself finish every mixed row.
12. [ ] Publish the combined 2Rc/Rp/Rv own seam receipt only after all ordinary, radial-pressure and required mixed joins pass. Then enable its checked accessor; keep require_checked=False until that receipt exists.
13. [ ] Complete actual own Rin/rminus/Rc joins from the admitted nonlinear core, real lower collar, finite-N modification and repaired own histories.
14. [ ] Assemble a checked global n1 source provider using those common-source joins. Preserve separate leading `1-delta`, own `1+delta`, fixed leading N1p and all signed five-moment histories.
15. [ ] Add original-scale quantitative containment for the new radial sources and partial pressure functions. Restore a numerical callback from genuine source data; never use caps, comparison covers or interval midpoints as derivatives.
16. [ ] Complete the common global N, axis/energy, analytic preheat/exact heat and admissible stress/flat remainder dependencies.
17. [ ] Construct genuine n>=2 coefficient recursion, independent moment repairs, finite-order remainders and smooth summation; coordinate rescaling is insufficient.
18. [ ] Add oscillatory stress cancellation and corrected Cartesian [u,v,w], then independently check the fixed-forcing NS Linf and physical volume-L2 targets below 1e-3.

## Stream recovery

A response-stream error does not classify any mathematical gate. Inspect saved archives, receipts and the current branch before retrying work. The supplied request ID is a platform diagnostic handle; local session logs identify the failed turn but do not reveal its server-side root cause. Keep the goal active and reuse accepted results.
