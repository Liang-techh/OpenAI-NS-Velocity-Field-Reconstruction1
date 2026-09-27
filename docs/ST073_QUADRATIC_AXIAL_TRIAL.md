# Quadratic axial shape restores sampled moments at k=11, at excessive cost

The existing streamfunction construction now accepts `axial_powers`,
defaulting to `(0,1)` for compatibility. The trial uses `(0,1,2)` for
swirl and poloidal factors, with the axial derivative included in the
streamfunction-derived radial velocity. Old amplitudes embed with zeros
in the new quadratic entries. Default-field samples match the previous
implementation exactly; a sampled new-mode finite-difference divergence
check has relative error 7.33e-10. This is not a global divergence bound.

At k=11, split 24-point quadrature and an analytic quadratic-moment
Jacobian yield a moment-only fit after 50 evaluations. All 18 local
amplitude increments stay within +/-40; max absolute increment is
22.56334. The four fitted moments have maximum 2.79e-11.

Direct independent field replay gives:

| Points per radial panel | Maximum absolute sampled moment |
| --- | ---: |
| 48 | 0.0004185183 |
| 96 | 0.0004190285 |

Thus this added axial freedom removes the old bounded-basis obstruction
for these four sampled moments at this scale. It is not complete-field
momentum acceptance, continuous moment closure, or scale recursion.

The quadrature-node full momentum peak increases from 1.65332e6 to
2.12913e7, a factor 12.8779. Only 3 of 9 physical support nodes pass the
48-point stress-cone check; maximum reported cone ratio is 1.86375.
The candidate is rejected. The polynomial axial factor is a local
exploratory ansatz; finite energy and global axial closure are unproved.

Next optimize momentum cost within the resolved moment-feasible family,
then restore the physical cone. Compare a minimum-norm constrained
solution and a full-momentum penalized solution before extending the
same axial modes to k=15/19 and intermediate times. Refine training
quadrature if small sampled moments are being used as a constraint.
Neither the 1e-3 full momentum maximum nor spatial-volume L2 goal has
been approached by this candidate.

Reproduce: `python experiments/root_st073/midplane_quadratic_axial_trial.py`.
