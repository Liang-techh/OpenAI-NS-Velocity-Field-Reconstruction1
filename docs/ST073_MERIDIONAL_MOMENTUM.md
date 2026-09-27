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
