# ST050 in-progress checkpoint — NOT scientific acceptance

Task #408, parent #390 at `97695a86f85ce68fb4ae70c41fc81c904d655183`.
The original ST048-S mathematical field was reconstructed and matched all eight published velocity references exactly. Its local JSON metadata differs from the original raw artifact.

## Actual completed fits so far

| Run | Variables | Iterations / evaluations | Stop | Minimum normalized training constraint |
|---|---:|---:|---|---:|
| C (same-configuration retry) | 398 | 163 / 190 | optimizer success | -2.52e-10 |
| P (pressure-oriented reconfiguration) | 398 | 173 / 229 | optimizer success | -4.19e-12 |
| E (edge/morphology continuation of C) | 542 | 213 / 268 | 480-second budget, NOT convergence | -4.67e-11 |
| PE (sequential recovery, parent P) | 542 | 159 / 203 | 480-second budget, NOT convergence | -8.45e-10 |
| PC (stronger pressure-Poisson penalty, parent P) | 398 | 159 / 196 | optimizer success | -3.29e-12 |

These are training checkpoints, not held-out residual reports. At this checkpoint neither planned full independent seed9175001 nor9175002 has been evaluated. Original nu, force family/bounds, support, energy normalization and scientific gates are unchanged. P deliberately changes the prior autonomous core-anchor allowance from8% to60% and removes a previous training-only no-worsening cap; this is documented, not called preservation of every prior auxiliary constraint. E/PE add actual vorticity distribution constraints.

## Other actual attempts and resource limits

An initial foreground C run was interrupted by a tool timeout; the completed retry is separate. The first parallel PE attempt was OOM-killed under the4GiB cgroup limit. An initial PC attempt was deliberately terminated to avoid another OOM; its sequential retry completed. All logs/checkpoints remain retained, not relabeled successful.

The temporal-expansion/reference-morphology experiment T was deliberately stopped after50saved iterations, with no feasible final candidate and repeated original parameter-bound failures. The saved minimum constraint was approximately-4269.97. This is a rejected numerical attempt, not an infeasibility theorem or proof that more time modes are useless.

A further HM run is currently testing necessary degree2-through12 harmonic pressure moments jointly with the already pressure-oriented field. Its numerical lower-bound target is not a replacement for full momentum acceptance.

## New checked identities

For smooth axisymmetric divergence-free u=(xA-yB,yA+xB,C),s=r² and the original divergence-free force:
`div R = tr((grad u)^2)+Delta p`; on the axis, `div R=6A²-2B²+4p_s+p_zz`.
At exact R=0, positive radial and axial pressure curvatures imply B²>3A². Small sampled |R| does not itself bound divR. This is NOT a whole-family no-go result.

On a finite cylinder, if the pressure has nonnegative outward normal derivative on its ENTIRE boundary, the frozen velocity obeys the conditional bound
`||R||inf >= max(integral tr((grad u)^2),0)/boundary_area`.
For frozen ST048-S on the stated small cylinder at t=.75, the quadrature estimate is .00434673. Exact conditional inequality, floating quadrature value, not interval-certified; its pressure-orientation hypothesis is not true of the current S field.

Eight focused new tests actually passed in12.11s, including symbolic divergence/axis-pressure identities, analytic derivative calibration, independent Cartesian residual-divergence refinement, vorticity moment derivatives and the cylinder flux identity. Full inherited suite and Lean were not run.

Original scientific1e-3 target and full source correspondence remain UNMET. No default, old artifact, other agent branch or schedule changed. A final source/recipe/independent-report delivery will replace this provisional handoff; no future result is claimed here.
