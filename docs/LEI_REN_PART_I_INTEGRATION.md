# Lei–Ren Part I integration — 2026-09-29

Source: Zhen Lei and Xiao Ren, *Finite-Time Blowup for Navier–Stokes with Smooth
Forcing, Part I: Construction of Self-Similar Solutions with Admissible Stress
and Flat Remainder*, [arXiv:2609.35406v1](https://arxiv.org/abs/2609.35406v1),
submitted 2026-09-28. [Full text](https://arxiv.org/html/2609.35406v1).
This document records implementation decisions after reading the paper, not
an independent verification of its proof. Source version is pinned to v1.

## What changes in our work

Part I constructs the axisymmetric background and its admissible stress.
It does not supply the oscillatory cancellation deferred to Part II.
Our final requirement remains a nonzero, divergence-free, finite-energy
time-dependent field with the requested morphology, repeated scale dynamics,
and full momentum maximum and volume L2 below `1e-3` under declared forcing.

We now distinguish three objects in every relevant report:

1. Background momentum residual `R_B`.
2. Profile-derived stress contribution `-D T_B` and background remainder
   `E_B = R_B + D T_B`.
3. Full corrected-field residual after actual oscillatory velocities,
   interactions, pressure, cutoffs and prescribed forcing are included.

A large first quantity need not reject a Part I background. A small second
quantity does not pass our final PDE gate. A freely chosen stress or force
cannot be used to manufacture acceptance. Stress must come from the same
profile, its moments and its derivatives, with support and cone evidence.
Existing nonaxisymmetric reference fits are not automatically Part I backgrounds.

The current ST073 comparisons remain structural experiments on their declared
domains. They do not replace the retained ST006 repository PDE baseline
(maximum 0.1082289305112118, volume L2 0.10758432876230622) without the same
validation contract. Neither route currently passes `1e-3`.

## Coordinate and parameter bridge

Equations (2.2)–(2.7) correspond to the existing implicit-similarity convention:

| Part I | Repository source coordinates |
| --- | --- |
| `lambda^2` | `q = tau/(1-eta^2)` |
| `R = r_source^2/(2 lambda^2)` | `X` |
| `Z` | `eta` |
| `delta` | `2h` |
| `lambda^(-1-delta)` | `q^(-A)`, `A=1/2+h` |
| `lambda^(1-delta)` | `q^D`, `D=1/2-h` |

Use `tau=T-t`; the paper uses `T=1`. For the current viscosity convention,
`x_source=x_physical/sqrt(nu)`, velocity is multiplied by `sqrt(nu)` and
pressure by `nu` when returning to physical coordinates. This bridge applies
to the implicit profile coordinates, not to an assertion that our global
Piola scale map is an exact Navier–Stokes symmetry.

Section 1.2 states `0<delta<1/200`. Therefore this source route requires
`0<h<1/400`. Existing `h=0.005` lies outside that stated range. No old
candidate or baseline parameter is silently changed. A separate versioned
seed must be built and its actual profile/constraints recomputed; substituting
the exponent into a frozen field is not a valid transfer of the paper's result.

## Construction order and compatibility

The useful implementation order in Sections 4–12 is: prepare the outer
connection and heat exterior; compute its actual axis pressure; construct
the regular core with that same pressure; connect velocities; restore five
moments; modify shear while restoring the moments again.

The five data in Section 2.5 are functions of `Z`, not five numbers at the
midplane: angular, axial, two mixed/quadratic moments, and the pressure
moment. The exterior angular moment is renormalized (1.2), so our compact
global axial angular-momentum cancellation is not a substitute for it.
Pressure normalization and the inner core cannot be chosen independently.

Theorem 1.1 uses compatible inputs. Theorem 12.4 supplies those inputs for
the paper's construction; we must demonstrate compatibility for our own
parameters and profiles. We must not describe the paper's construction as
leaving Assumption 12.1 unproved, nor inherit its conclusion for our candidates.

Section 8.2's linear model is useful for core initialization and exit-direction
checks. It matches the nonlinear regular core's axis values and first radial
slopes. It is not itself a general stress-free nonlinear solution. The new
`linear_core_axis_slopes` helper implements (8.7) with an externally supplied
consistent `P0(Z)`; it does not choose a replacement pressure.

## Stress, cone, and radial equation

Sections 2.4 and 17.2 use

`D T = (0, d_r Ttheta + 2 Ttheta/r, d_r Tz + Tz/r)`.

The Cartesian symmetric tensor realizing this operator needs an extra
theta-theta entry `r*d_z Tz`, in addition to its r-theta and r-z entries.
Without it, its divergence introduces a spurious radial term `d_z Tz`.
The helper `completed_stress_tensor` implements (17.24); physical stress
jets must already include cutoff derivatives. It is distinct from the
Cauchy stress `-p I + nu(grad u + grad u^T)` in `interface_stress.py`.

Since `D T` has zero radial component, `E_B^r=R_B^r`. We cannot discard radial
momentum while evaluating the two stress components. Nor does completing the
tensor assert a separate cone condition for its added diagonal component.

The admissible cone uses the actual shear and stress, including the positive
swirl and negative angular-shear signs. Check `kappa>2`, the negative stress–shear
dot product and the strict angular-width inequality (1.3). Section 17 requires
directional margins near vanishing-stress edges, not merely a finite list of
interior passes. Reuse current cone code after verifying its variable conventions.

## What flatness and recursion mean here

Sections 13–16 require higher-order coefficients, moment restoration and
cutoff-aware smooth summation. A few decreasing sampled residuals do not
establish their all-orders conclusion. Record formal order separately from
finite numerical estimates for each fixed Cartesian space/time derivative.

Section 17.2 converts lambda-flatness to tau-flatness on fixed interior
sectors `R<=K`, `|Z|<=1-epsilon`. It does not make that conversion uniformly
on the closed endpoints `Z=±1`. Our reports must retain the sector and
derivative order rather than call three normalized scalar norms “flatness”.
Existing near-constant pulled-back residuals provide no such evidence.

After a compatible background exists, the next route is to realize its stress
through actual waves and control every remaining interaction. Existing
wave-mean flux, analytic curl jets and cached full residual machinery can be
reused. Part I alone does not justify marking that phase complete.

## Reuse map from the current repository

Read-only source inspection identified these reusable implementations. Their
existence is not evidence that the new background has been assembled.

| Need | Existing implementation | Integration gap |
| --- | --- | --- |
| Regular nonlinear core | `src/openai_ns_reconstruction/paper_core_series.py`: `build_paper_core_series`, `paper_core_profile`; `experiments/root_st073/NS_ST073_Full_Local_Recurrence/full_radial.py` | Recompute under the compatible axis pressure and parameter manifest; local series do not give an outer connection |
| Five-moment repair | `src/openai_ns_reconstruction/background_moment_repair.py`: `Lemma52MomentRepair.solve`, `corrected_u_value`, `corrected_e_value`; `moments.py` | Existing finite supplied-moment repairs need actual moment functions from one profile, including exterior pressure normalization |
| Cone / wave stress | `src/openai_ns_reconstruction/stress_cone.py`: `stress_target`, `PositiveSignedStressDecomposition`; `stress_cone_perturbation.py` | Bind actual profile stress/shear, pulse integrals and edge margins; supplied scalar error bounds are not measured field evidence |
| Full residual and truncation | `src/openai_ns_reconstruction/background_truncation_residual.py`: `recurrence_truncation_breakdown`; existing ST073 cached momentum evaluators | Bind the physical profile-derived stress and all cutoff derivatives to `R_B`, `D T_B`, `E_B` |
| Formal summation / flatness | `endpoint_borel.py`: `BorelRightExtension`; `endpoint_scale_schedule.py`: `SpatialBorelScaleSchedule`; `section9_flat_remainder.py`: `certify_flat_remainder_power_ladder` in the same package | Finite supplied jets/witnesses and caller-provided bounds do not establish all-order closure; materialized coefficient corrections and measured finite-order remainders come first |

Do not restart these modules, or make every historical theorem wrapper a new
prerequisite. The immediate missing connection is an actual shared
profile–pressure–moments–stress data path. The new Part I helper is an explicit
interface for that path, not a claim that this missing connection is closed.

## Execution queue

These are project tasks, not claims that every original proof must be reproduced.
Use one candidate identity and keep paper facts, choices and observations separate.

| ID | Task and concrete output | Depends on | State / acceptance |
| --- | --- | --- | --- |
| LR1-01 | Pin source, map variables, record source/goal boundaries | none | DONE: this document; no scientific promotion |
| LR1-02 | Implement stress completion, radial remainder, axis-slope and sector interfaces | LR1-01 | DONE: `lei_ren_part1.py`; manufactured Cartesian-divergence check |
| LR1-03 | Create an explicit parameter manifest in the paper's stated delta range, with nu, tau interval, fixed interior sectors and nontriviality limits | LR1-01 | OPEN; retain old seed and ST006 baseline |
| LR1-04 | Bind an existing outer/heat profile to its full `P0(Z)` and five moment functions; identify all missing moments and tail normalizations | LR1-03 | OPEN; one source identity, no independent pressure refit |
| LR1-05 | Use the Section 8.2 model to initialize a regular nonlinear core under that pressure; compare axis slopes and exit signs | LR1-04 | OPEN; model agreement alone is not core acceptance |
| LR1-06 | Connect core and exterior and restore all five moment functions over the declared Z range | LR1-05 | OPEN; verify pressure and exterior field restoration, not just midplane scalars |
| LR1-07 | Recompute actual stress and shear; check interior cones and both edge directional margins | LR1-06 | OPEN; distinguish relaxed from admissible cone |
| LR1-08 | Evaluate paired `R_B`, `D T_B`, `E_B` including radial momentum and all cutoff terms on matched grids | LR1-07 | OPEN; never infer E_B by assigning arbitrary T_B |
| LR1-09 | Implement a first lower-order coefficient correction with consistent axial viscosity, moments and cutoff derivatives | LR1-08 | OPEN; demonstrate a measured remainder-order gain, preserve support/cone |
| LR1-10 | Repeat on at least three scales and declared Cartesian derivative orders; distinguish finite-order evidence from infinite flatness | LR1-09 | OPEN; fixed interior sector and independent discretization |
| LR1-11 | Fit actual nonaxisymmetric wave flux to the profile-derived stress and recompute full interactions | LR1-07/09 | OPEN; Part II/OpenAI correction mechanisms still required |
| LR1-12 | Validate repeated dynamics, fixed exterior/forcing regularity, nonzero field, energy and three morphology observables | LR1-11 | OPEN; no free residual-defined forcing or amplitude collapse |
| LR1-13 | Run original full max/L2 `1e-3` gate; reuse Python/MATLAB export and visualization only after field selection | LR1-12 | OPEN; same-contract ST006 comparison and independent holdout |

Reproduce the completed algebraic checks:

```powershell
python experiments/root_st073/lei_ren_part1_checks.py
```

The generated JSON keeps `pde_validated=false` and
`scale_recursion_established=false`. It does not certify an existing profile.
