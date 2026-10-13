# Full original n1 defect C1 functions and selected linear repair

The five true signed moment functions are now complete. The inner cutoff
endpoints supply m11/m12/m14, and both pressure moments add all three
exterior source intervals exactly once. These genuine C1 functions select
the two axial and three angular repair coefficients using the original
reserved physical bump matrix. C1 repair velocities also have actual
function bodies.

The five algebraic physical-response identities and their first Z
derivatives are checked. This does not install the repaired own M/V/P,
transfer coefficient norms or higher derivatives, or certify global
five-moment closure. The full goal stays **ACTIVE / INCOMPLETE**.

## Restore and source ownership

```powershell
python scripts/run_original_n1_full_defect_linear_functions.py
```

The checked quartet under `experiments/root_st073` is
`lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`).
`CurrentOriginalN1FullDefectLinearFunctions()` restores the accepted graph
of 53,398 nodes. Its interface exports `.moments`, `.coefficients`,
`.response_matrix`, `.linear_response_residuals`, `.repair_velocity`,
`.segments`, `.partition` and `.reserved_background`.

The exact retained owner chain is

```
full.cutoff.inner.collar.core.pulse
    .own.jets.complete.Gamma.pulse.preRp.operator
```

The lower collar, preRp, pulse, cutoff and reserved operator use the same
graph object and source family. Each owner's accepted prefix and source
hashes is checked. No graph import, numerical pressure grid or serialized
signed-sector selection is needed. Reuse these saved ancestors.

Local installed gates are
`n1_complete_signed_defect_C1_function_adapter_installed`,
`n1_selected_linear_coefficient_functions_installed`, and
`current_original_n1_linear_repair_velocity_C1_functions_installed`.
Earlier receipts remain historical and retain their original narrower scope.

## Full signed pressure partition

For j=3,5,

```
m1j_full = m1j_inner(0,Rin)
          + m1j_lower(Rin,r_minus)
          + m1j_preRp(r_minus,Rp)
          + m1j_pulse(Rp,Rv)
```

The preRp interval contains all 17 active charts, the compact Rc-to-2Rc
band and the quiet 2Rc-to-Rp interval. The pulse includes its complete
main/gap/end integrals and finite endpoints. Every m15 segment already
contains its negative weighted-pressure sign; all four functions are
added without another negation.

The accepted original meridional support theorem gives exact zero N1p
and R*N1p integrals from Rv to infinity. It retains angular velocity and
pressure/heat fields there. Gamma/future functions are upstream pulse
inputs, not another pressure segment. The complete velocity-support
m11/m12/m14 are the cutoff's true Rin endpoint functions.

## Exact reserved background and bump normalization

The new centers are j=(12,13,14,15,16), measured in log(R/Rc). Their
half-width is log(2)/400, with center j*log(2)/40. Axial bumps use 12,16;
angular bumps use 13,14,15. The five supports lie in the reserved gap and
are disjoint from the old compact corrections at j=8,20,32.

The saved actual band V rows are the full leading axial profile in its
physical common unit. Zero propagation through their arithmetic, using
the exact disjoint flat supports, proves all three y rows and four Z rows
of V and delta_E vanish on the new supports. Thus the leading axial
velocity is zero and the swirl is its original power there. Different
preserved C2 arithmetic node IDs are not treated as different functions.
Nonlocal mean, pressure and the fixed leading Omega/N1p remain unchanged;
they are not assigned zero in this reserved band.

Each physical bump b satisfies integral b*dR/sqrt(2R)=1. Write
Utheta0=Astar*(2R)^(-1/2-mu), with the original positive mu retained.
The physical matrices are

```
Bz rows:     integral (2R)^0 b_z dR,
             integral (2R)^(-mu) b_z dR
Btheta rows: integral (2R)^(1/2) b_theta dR,
             integral (2R)^(-3/2-mu) b_theta dR,
             integral (2R)^(-1/2-mu) b_theta dR
```

The selected coefficient functions solve

```
Bz * alpha = -(m11, m14/Astar)
Btheta * beta = -(m12, m13/(2Astar), -m15/Astar)
```

The stable divided axial inverse keeps the nonzero original mu. Its
positive denominator equals the negative normalized axial determinant.
The angular inverse retains the original positive Vandermonde determinant.
The accepted matrix/admitted amplitude guards cover the original source;
they do not yet bound the resulting repair coefficients.

The actual moment changes are

```
Delta m11 = Bz_row0 * alpha
Delta m12 = Btheta_row0 * beta
Delta m13 = 2 Astar * Btheta_row1 * beta
Delta m14 = Astar * Bz_row1 * alpha
Delta m15 = -Astar * Btheta_row2 * beta
```

The factor two and minus sign are physical identities, not conventions
to change during composition. Astar_Z is retained in all C1 quotient and
product rules. The C1 velocity patches are the corresponding physical
bump sums; the angular F1 patch is Utheta1/sqrt(2R) once. Their supports
are away from the axis, with the canonical zero axis extension.

## Evidence and remaining limits

An independent signed polynomial partition checks both value and first Z
rows and the single m15 sign. The universal proof reads the actual selected
coefficient/physical response expressions. It substitutes the previously
accepted exact determinant identities and proves all ten residual rows
zero for arbitrary C1 defects and Astar/Astar_Z, keeping mu and Rc.
This is an algebraic response theorem; it is not a numerical original
point-field calculation or a repaired own-field terminal theorem.

Independent physical dR integrals with Rc=3.7, finite mu=.08 and variable
Astar check the actual normalization, response rows, first derivative
cancellations and angular F1 units. Its finite-fixture difference is below
4.53e-71. Source/archive guards and checked CLI restoration complete the
handoff. These fixture errors are not uniform original parameter accuracy.

The coefficient C1 functions and the formal linear response are installed.
Original numerical n1 point evaluation, coefficient/C5 norm transfer,
genuine C3 repair velocity functions, repaired own M/V/P, terminal five
zero identities, global heat/stress and temporal recursion remain open.

## Next executable tasks

- [x] Compose all five true signed defect C1 functions, with the full
  lower/preRp/pulse pressure partition and beyond-Rv support theorem.
- [x] Select actual two axial/three angular C1 coefficient functions.
- [x] Build physical compact C1 repair velocity functions with one
  Utheta-to-F conversion and the exact disjoint-support source guard.
- [x] Verify all five physical linear-response identities through Z1.
- [ ] `N1-DIRECTED-FUNCTION-EVALUATOR`: attach the same-source complete
  regular-series encloser, proving the exact physical conjugation and
  parameter/amplitude/P0/chart identities first. Publish physical
  factored/log outputs and full error tails; do not select centers as fields.
- [ ] `N1-DEFECT-HIGHER-DERIVATIVES`: extend genuine pressure defect
  functions beyond C1 as required for repair C3 and original C5 bounds.
  Differentiate their actual field/integral/implicit definitions, retaining
  ordinary/Taylor factors and endpoints; do not differentiate norm envelopes.
- [ ] `N1-COEFFICIENT-BOUND-TRANSFER`: identify the old five-defect covers
  with these exact same-source functions, including nonzero source errors
  and the complete pressure intervals. Apply the stable physical inverse
  to obtain directed coefficient and derivative ranges on the original Z
  cells/tubes. Matrix positivity alone is not a coefficient bound.
- [ ] `N1-REPAIR-VELOCITY-C3`: generate true coefficient Z2/Z3 functions
  with the amplitude derivatives retained; keep all compact support and
  physical-scale identities through radial/axial differentiation.
- [ ] `N1-REPAIRED-OWN-FIELDS`: add the bumps to the true cutoff velocities,
  reconstruct own M1 by its full signed integral, then K1/Q1/V1 and P1
  with the n1 mean weight 1+delta and the unchanged leading N1p. Preserve
  lower/core/pressure histories; do not multiply M/V/P by the velocity cutoff.
- [ ] Prove all repaired terminal moment zero identities as functions,
  verify their own radial FTC/axis data and required smooth joins, and
  transfer the original quantitative bounds before closing the n1 repair gate.
- [ ] Complete global matching, energy, exact heat exterior and admissible
  stress/flat remainder; then construct actual n>=2 temporal recursion.
- [ ] Add oscillatory stress cancellation and independently validate the
  corrected Cartesian NS residual with fixed/restricted forcing.

The full objective remains in [PROJECT_GOAL.md](PROJECT_GOAL.md).
