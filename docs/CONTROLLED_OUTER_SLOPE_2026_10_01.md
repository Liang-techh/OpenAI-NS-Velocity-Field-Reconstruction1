# C120-F16: fresh reference continuation and first outer slope transition

The fresh local axial family now has a callable O.1 reference continuation
after the five-bump repair and an O.2 angular slope transition on
`y=log(R/Rref)` in `[0,1]`. This is an additional physical field stage, not
temporal coefficient recursion or a complete exterior construction.

Implementation: `experiments/root_st073/lei_ren_part1_paper_interval_outer_slope_field.py`.
Production receipt and independent check have the same basename with `.json`
and `_check.py` / `_check.json` suffixes.

## Field and inherited data

`IntervalOuterSlopeField.evaluate_reference_x(x)` uses `x=R/Rm`,
`2 <= x <= exp(6)`, on the exact reference branch. The uniform implicit
five-moment inverse proves the terminal moment identities at `x=2`.
The actual solution satisfying those identities supplies this continuation;
midpoint controls and residual interval boxes are not treated as exact roots.
The terminal identities propagate by the reference radial moment equations.

`evaluate_y(y,cells=256)` implements Eq. (4.6) for `0 <= y <= 1`:

\[
U^\theta=\frac{e^{14}}{1+Z^2}
  \exp\{y/10-(3/5)J(y)\},\quad
J(y)=\int_0^y\sigma(s)\,ds,\quad U^z=4Z.
\]

The angular slope changes from `1/10` to `-1/2`. All five cumulative moments
are inherited from the reference inlet and advanced by their defining radial
integrals. The pressure remains `P0+Mp`, using the same analytic axis datum.
The radial velocity and radial derivative are recovered from the axial moment.
No velocities, moments or pressure are independently fitted.

Directed monotone rectangles bound `J` and the three angular moment integrals.
On every rectangle, the complete range of `J` bounds the exponential integrand.
No finite quadrature result is presented as a certificate. Symmetry gives the
exact terminal identity `J(1)=1/2`; the exact terminal normalized shear is `-2`.
The field has the existing local `C_Z^1` scope, `Z in [0.49,0.51]`; `Ur_Z`
and whole-axis higher smoothness remain unavailable.

## Evidence and limits

The independent moderate-scale fixture integrates a coupled primitive ODE
with SciPy DOP853 and checks 71 enclosed quantities, including all five moment
values and axial derivatives, pressure, and intermediate angular fields.
This numerical fixture checks formulas and is not a rigorous numerical oracle.
Additional checks cover inlet C1 moment overlap, exact terminal J/shear,
endpoint angular amplitude, and nested rectangle refinement.

Six production point packets certify the relaxed cone. These are point
diagnostics, not complete radial coverage or strong admissibility. No heat
exterior, temporal recursion, oscillatory correction, or full residual claim
is made.

## Discovered parameter mismatch and next work

The accepted pressure schedule SHA is
`736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4`,
validated through `lei_ren_part1_paper_reference_endpoint_targets.py` and
`lei_ren_part1_paper_coherent_pressure_source_alignment.json`.
It fixes **Md=0.5**. The cached paper text
`work_paper_cache/lei_ren_part1.txt`, lines 3418–3428, requires **Md>1**,
a sufficiently large auxiliary constant, and `Pstar>exp(Td)` with
`Td=exp(Md)+10`. The new initial slope segment does not use Md. The old
schedule is therefore not evidence that the full paper parameter hypotheses
hold. Do not silently continue Md=0.5 and label the result fully paper-faithful.

- [x] Extend the repaired terminal reference field to Rref.
- [x] Install the first O.2 angular slope transition with inherited C1 moments.
- [x] Preserve the same analytic pressure datum and provide directed integrals.
- [ ] Resolve the outer parameter regime against the actual paper inequalities,
  including how a changed Md/Pstar propagates through the fourteen pressure
  masses, core/source data, defects, and implicit repair controls. Preserve
  the existing receipts; generate a coherent new companion chain where needed.
- [ ] Implement the axial cutoff `B(y-1)=1-sigma(log(y)/Md)` for y>1,
  retain accumulated Mz after Uz becomes zero, and transport all five moments
  and pressure from this slope endpoint.
- [ ] Implement the next angular turn from -1/2 to -1/2-mu, and the power buffer
  reaching the reserved Section 11 interval with Uz=0 and kappa>2.
- [ ] Certify complete radial cells; resolve angular normalization near kappa=2
  before asserting strong cone prerequisites.
- [ ] Implement pulse, flatten, angular repair, waiting and heat layers using
  this same functional family. Then proceed to genuine temporal recursion.

Reproduce with Python 3.11 by running the field module and its check module.
The earlier core, connecting and repaired-annulus receipts remain frozen.
