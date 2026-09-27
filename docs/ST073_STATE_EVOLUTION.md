# Nonlinear state evolution diagnostic

The 165-point derivative fit now advances the 25 velocity coefficients using
new slopes and pressure coefficients at every changed scale. The full
quadratic advection is retained. Two explicit Euler steps reach each tested
endpoint, where an independently sampled 176-point Cartesian residual is
evaluated using the endpoint's local tangent field.

| Delta k | Fixed-slope maximum | Refreshed maximum | Fixed volume L2 | Refreshed volume L2 |
|---|---:|---:|---:|---:|
| 0.001 | 4.38230e9 | 4.20736e9 | 189537.55 | 180899.30 |
| 0.01 | 6.00640e10 | 2.27427e9 | 1473368.61 | 84464.34 |

These results in `meridional_state_evolution.json` support state-dependent
updates instead of frozen initial slopes. They do not establish a converged
time integrator, a differentiable field across step boundaries, maintained
moment/cone compatibility, or the required residual threshold. No actual
nonaxisymmetric waves are included. The baseline broad-shear mean already
differs substantially from the earlier low-residual core continuation.

The earlier 44-point fit is retained as
`meridional_state_evolution_coarse.json`. With only two axial nodes it failed
to distinguish some quadratic axial modes, and reported rank after ridge
augmentation. Its poor extrapolation should not be interpreted as physical
infeasibility. The dense initial matrix has unregularized rank 44, condition
4.1312, and smallest scaled singular value 0.3944. Both unregularized and
augmented spectra are now recorded explicitly. Spatial holdouts are common
to the comparisons, and both integral L2 and volume-normalized RMS are saved.

Four steps instead of two at delta k=0.01 give maximum 2.00952e9 and volume
L2 91731.66 (RMS 2.60911e8), compared with 2.27427e9 and 84464.34 for two
steps. State norms differ substantially, 215.913 versus 280.624. These
results in `meridional_state_evolution_refined.json` demonstrate sensitivity,
not convergence. The fixed-cylinder refined observables still broaden and
weaken: radial RMS 0.00140017, aspect ratio 0.445366, weighted angular speed
1.16941e6. Thus the wrong target trend persists in both tested step counts.

## First dynamically constrained step

`meridional_constrained_evolution.py/.json` advances a zero correction state
by delta k=0.0001, rebuilding the nonlinear residual, growth geometry,
moment equations and cone inequalities at the endpoint. A fresh feasibility
phase precedes the scaled constrained least-squares fit. The original basis
knots remain fixed. Four moment rows use the order-96 integrated identities;
cone assembly uses order 48 and independent replay order 64.

| Quantity | Initial | Endpoint |
|---|---:|---:|
| Independent momentum maximum | 4.44494e9 | 4.43477e9 |
| Physical-volume L2 | 188619.22 | 188215.12 |
| Integrated moment maximum | 2.91806e-8 | 8.88285e-8 |
| Sampled cone passes | 27/27 | 27/27 |

This is one explicit state step with compatible initial/endpoint tangent
fields. Interior compatibility and derivative continuity at step boundaries
are not established. SLSQP did not certify convergence; the selected points
were checked for feasibility and objective improvement. Full residuals still
miss the requested 1e-3 thresholds by many orders of magnitude.

`vortex_constrained_observables.json` reports the same fixed-cylinder
diagnostic at this endpoint: radial RMS 0.000559536, aspect ratio 0.934013,
weighted angular speed 2.11821e6 and maximum swirl 2020.66. The shape changes
are tiny and lack a spatial convergence study; do not call them successful
shrinkage or slenderization. The sampled spin measures do not increase.
This supplies dynamic compatibility, not the missing amplification mechanism.

## Target-shape diagnostic

`vortex_state_observables.py/.json` measures the actual velocity curl on a
fixed physical cylinder: radius 0.00441896, axial half-length 0.000938970,
1944 quadrature points. The domain does not shrink with the similarity
coordinates. At delta k=0.01, the refreshed path changes:

| Observable | Initial | Refreshed endpoint |
|---|---:|---:|
| Enstrophy-weighted radial RMS | 0.000559513 | 0.00167749 |
| Enstrophy axial/radial RMS ratio | 0.933673 | 0.390416 |
| Enstrophy-weighted signed angular speed | 2.12518e6 | 967990 |
| Sampled maximum swirl speed | 2027.25 | 1153.82 |
| Cylinder kinetic energy | 0.00656263 | 0.00262920 |

The observed trend is broadening and weakening, not the requested narrowing,
relative axial slenderization and spin-up. Reject this unconstrained path as
evidence of target amplification even though its residual improves. These
are whole-cylinder weighted measures, not a uniquely identified vortex core;
the cylinder energy is not global energy. No quadrature-convergence claim is
made. An analytic solid-body rotation check verifies the curl, radial/axial
moments, angular speed and enstrophy formulas.

Next actions are step-size sensitivity and state-dependent constrained
evolution using the integrated conservation rows. The construction must
also couple actual nonaxisymmetric waves and monitor target growth. Pure
mean residual minimization can lower the error by weakening the swirl and
therefore cannot by itself be the recurrence objective.
