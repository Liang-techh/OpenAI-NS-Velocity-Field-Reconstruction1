# Coupled meridional momentum correction

## Completed finite-basis experiment

`experiments/root_st073/broad_meridional_momentum.py` and its JSON implement
a reconstructable 44-direction field: 15 meridional slopes, 10 swirl slopes,
18 compact pressure modes and one free compact centrifugal-pressure direction.
The last direction allows the fit to undo the initial pressure primitive.
Fixing that primitive at unit strength reduced peak error but increased volume
L2; the free coefficient is -0.780356, leaving net primitive strength 0.219644.

At k=11.0003, nu=0.01, the domain is the full angular annulus with
eta in [-0.5,0.5] and y in [0.01,0.99]. Training uses 165 nodes and replay
uses 176 nodes with different axial/radial Gauss orders and a rotated angle.
Radial quadrature is split at the declared transition/basis breakpoints.
The L2 values below use the physical volume integral, without normalization.

| Field on the same holdout | Sampled maximum | Volume L2 estimate |
| --- | ---: | ---: |
| Original 19-control broad mean | 7.01009e9 | 237497.46 |
| Coupled 44-direction correction | 4.38510e9 | 186820.70 |

The decreases are 37.4% and 21.3%, respectively. The corrected volume RMS
is 5.2862e8. These are local, instantaneous numerical improvements; they
remain far from the 1e-3 gate and do not replace the repository's retained
PDE baseline. Quadrature/differentiation convergence is not established.
Finite-difference divergence max remains 0.154809 for both fields because
the added velocity is zero at this instant; no numerical divergence gate
is claimed. The new slope increments are structurally solenoidal.

Moment and stress-cone constraints were not imposed in this fit. Their earlier
pass counts do not transfer to this candidate. Next, assemble the actual
coupled pressure/meridional columns in those constraints and solve the full
momentum objective subject to them. Only a compatible spatial improvement
should proceed to nonlinear time evolution and inter-scale contraction.

Run `python experiments/root_st073/broad_meridional_momentum.py` to regenerate
the report; `load_saved_field()` reconstructs the candidate. No finite-time
evolution was run. Initial sparse-grid overfitting and an omitted mapped-weight
factor were corrected before the saved final report; no sparse-grid result
is used in the table above.

## Construction and source

### Off-time replay: local benefit does not extend across one scale

`broad_meridional_time_audit.py` replays the saved callable with zero forcing
and no refit. `broad_meridional_time_audit.json` covers delta k in
{0, 1e-5, 1e-4, 1e-3}; `broad_meridional_scale_audit.json` covers {0.1,1}.
The same similarity annulus and 176-node quadrature are used at each time.
These are evaluations of an explicit affine-k field, not NS time integration.

| Delta k | Corrected momentum max | Corrected volume L2 | Both below the dynamic reference? |
| --- | ---: | ---: | --- |
| 0 | 4.38510e9 | 186820.70 | yes |
| 1e-5 | 4.38114e9 | 186772.48 | yes |
| 1e-4 | 4.34574e9 | 186363.32 | yes |
| 1e-3 | 4.01805e9 | 184723.73 | yes |
| 0.1 | 1.52546e12 | 3.69057e7 | no |
| 1 | 3.82489e14 | 5.95543e9 | no |

At delta k=1, remaining time has halved. The large residual rejects naive
fixed-slope continuation across that scale. The combined broad-amplitude
slope is -1316.5884 with initial amplitude 44.8774, so that coefficient
crosses zero already at delta k=0.034086. Sampled peak speed and local
kinetic energy decrease over the short tested interval and then grow sharply
under extrapolation. This is not evidence of the requested recursive vortex
growth. Pulse decay by itself is not a contradiction of the paper's mechanism.

The next evolution must update slopes and pressure from the changed state,
retaining nonlinear interactions and rechecking compatibility. The pressure
primitive in this saved candidate has a fixed amplitude parameter; this audit
does not silently reconstruct a different pressure away from the initial time.
The small-domain kinetic-energy values in the JSON are not global finite-energy
evidence. The unchanged 1e-3 gates remain unmet at every tested time.

Regenerate the larger-step record with:

```text
python experiments/root_st073/broad_meridional_time_audit.py --delta-k 0.1 1 --output experiments/root_st073/broad_meridional_scale_audit.json
```

The broad-shear field has a positive instantaneous wave-growth direction but
its full momentum residual is large. The previous compact centrifugal-pressure
primitive transfers error between radial and axial equations. A coupled
velocity-time and pressure update is therefore required before evolving it.

[The reference paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Section 8, equations (8.1)–(8.3), retains radial/axial mean corrections and the
complete wave covariance. Its compact primitive keeps the cutoff remainder.
Section 9 treats pressure as part of the evolving correction state and
recomputes the full residual. These motivate the coupled experiment here;
the finite basis and least-squares solver are autonomous numerical choices,
not implementations of the paper's supported inverse or recursion theorem.

For a compact axisymmetric streamfunction Psi, define

```text
B_r = -(1/r) partial_z Psi
B_z =  (1/r) partial_r Psi
delta u = (k(t)-k0) B,       k(t) = -log2(2 tau), t = -tau.
```

Mixed derivatives cancel in div B. At the reference time delta u and all
its spatial derivatives vanish, while partial_t delta u = B/(tau log 2).
Thus the instantaneous full momentum update is exactly linear in these
time-slope coefficients and in an added pressure gradient. The meridional
basis also preserves integral(r B_z dr)=0 where Psi vanishes at both radial
ends. This identity does not ensure every momentum or stress-cone constraint.

Pressure gradients must be solved together with the solenoidal acceleration:
a solenoidal time correction alone cannot generally cancel the gradient part
of a residual. Full radial, azimuthal and axial components enter the fit.
No unrestricted residual-dependent force is introduced.

## Scope and acceptance

`instantaneous_control_columns.py` supplies a limited initial-time assembly
operator: analytic physical-time slope derivatives, analytic compact bump
pressure gradients, and the required finite differences of the centrifugal
primitive. It rejects velocity slopes away from their registered k0 and
requires zero-background direction types. It does not replace the nonlinear
jets of an evolved state. A four-point off-axis comparison of the 44 columns
against full Cartesian FD gives maximum absolute difference 4.48013e-5 and
maximum column-scaled difference 8.67373e-10. This is derivative consistency
evidence, not the 1e-3 PDE gate. The measured cold analytic-path runtime was
17.00 s versus 11.18 s for the subsequent FD path; that small ordered test
does not establish a speedup. Tangential moment assembly can request only
the axial pressure derivative to avoid unnecessary primitive evaluations.

The comparison must keep the instantaneous nonzero broad-shear velocity,
viscosity and spatial domain fixed. Report both integral volume L2 and volume
normalized RMS distinctly, alongside the maximum Euclidean residual norm.
In coordinates (y,eta,theta), the volume element is

```text
q = tau/(1-eta^2), D = 1/2-h,
r_i = sqrt(2 nu q join_X), width = (ratio-1) r_i,
r = r_i + width*y,
z_eta = sqrt(nu) q^D [1 + 2 D eta^2/(1-eta^2)],
dV = r * width * z_eta dy deta dtheta.
```

An independent grid must evaluate the saved callable with Cartesian
derivatives, not just reuse the fitted linear matrix. A correction that loses
the previous moment/cone compatibility remains exploratory even if its full
residual decreases. A single-time improvement is not time integration, finite
global energy, or scale recursion. Those remain required by PROJECT_GOAL.md.
