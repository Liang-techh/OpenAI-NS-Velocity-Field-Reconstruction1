# ST056 completed results and route diagnosis

Task #746. Same full physical problem as ST054; no publication-default change. All fields were frozen at 2026-09-20T02:07:58.725762+00:00 before both fresh samples. No tuning or candidate selection after the holdouts.

## Six actually completed original-protocol validations

Each seed uses 4096 Cartesian points, six original times, the complete separate spatial/time refinement ladders and energy-quadrature ladder. Full-vector Euclidean maximum and fixed-time spatial volume L2=sqrt(64*mean(|R|^2)); worst time at spatial h=.005, temporal h=.0025. All six scientific CLI commands returned exit 1, failing momentum_max and momentum_L2. All other original sampled gates pass. Sampling is not continuum certification.

| Field | Seed | Full maximum | Spatial volume L2 |
|---|---|---:|---:|
| ST054-Q2 | 9205691 | 0.03776742629978789 | 0.047003537390649726 |
| ST056-J2 | 9205691 | 0.037759414075332606 | 0.04653784312879544 |
| ST056-M1 | 9205691 | 0.038534893979973166 | 0.05998556746363937 |
| ST054-Q2 | 9205692 | 0.03764253081582378 | 0.04854318006629012 |
| ST056-J2 | 9205692 | 0.03760476164201231 | 0.0479179815831483 |
| ST056-M1 | 9205692 | 0.03834656385586908 | 0.06149837003789147 |

J2 lowers L2 by 0.990764% / 1.287922%, and sampled maximum by only 0.021215% / 0.100336%. This is modest paired improvement, not a peak breakthrough or global superiority. M1 worsens L2 by 27.619262% / 26.687971% and max by 2.032089% / 1.870313%. It is a negative structural control, not promoted.

Additional fixed-sample peak checks select the eight largest analytic-residual locations at each original time, then independently vary only spatial finite-difference step. Worst norms at h=.005/.0025/.00125:

- Seed 9205691, Q2: .037767426300 / .037776548139 / .037777128420.
- Seed 9205691, J2: .037759414075 / .037768544304 / .037769125121.
- Seed 9205692, Q2: .037642530816 / .037646798801 / .037647064215.
- Seed 9205692, J2: .037604761642 / .037609034549 / .037609300268.

The tiny J2 advantage persists; it remains tiny. These are refinements at selected sample locations, not new continuous maximizations. Original reports are unchanged.

## New spatial and structure audits

Independent 55x89 spatial grid including the axis, thirteen times: max Q2 .038696100090; J2 .038477615686; M1 .039121435713. Q2 peaks on the axis at t=.708333; J2 and M1 peak at t=.25, r=.29925, z=-1.9025. Each time's grid peak has an independent Cartesian-FD vector cross-check. No sampled maximum is an upper bound.

On 800 new core points at 17 times (15 random interior times plus endpoints, seed9205693), all three fields retain all five measured inward/swirl/bipolar/radial-pressure/axial-pressure directions. The range is R[.035,.215], |Z|[.025,.215], not the entire support. J2 changes core swirl by at most .303918% on this sample; poloidal velocity is unchanged. In 1025 fresh midplane probes its axial velocity and radial shear differences are exactly zero. Shape is essentially unchanged: comparable new-grid final drift .196026044553 -> .196020940822. This is NOT a substantial geometry improvement. Fine enstrophy-moment quadrature also shows nearly unchanged global extents.

M1 changes physical poloidal amplitude by 2.080673%, despite retaining initial total energy. Its shear amplitude changes accordingly. It is not an unchanged-Q2 initial field and is not an improvement merely because it retains sign tests.

## Global compatibility diagnostic: numerical evidence, exact identity

The weak-moment idea already existed in ST053; this round applies it specifically to the locked initial field. It is not claimed as a newly invented general identity.

For smooth compact incompressible u, compact p, and divergence-free compact independently prescribed force, phi=(x/2,y/2,-z) is the gradient of the harmonic quadratic H=(x^2+y^2)/4-z^2/2. Integration by parts gives

```
Integral phi dot R = D(t) = Integral [uz^2 - (ux^2+uy^2)/2]
||phi||_L2(r<2,|z|<2)^2 = 88*pi/3
||R||_L2 >= |D(t)| / sqrt(88*pi/3)
```

At t=.25 the Q2/J2 96-point Gauss result is D=-.04034981786822833, yielding a lower-bound ESTIMATE .004203256891236408. Orders48/72/96 give .004203258047/.004203256877/.004203256891. This strongly indicates that preserving the entire Q2 initial velocity cannot reach true L2<.001. The integral value is not interval-certified; this is not a family-wide impossibility or a proof from Monte Carlo samples. The observed residual near .047 is still much larger than this lower estimate, so the identity does not explain the whole current plateau.

With initial energy one, D=(3/2)*Integral uz^2-1. M1's initial redistribution makes D approximately -1.33e-15 at the same quadrature, but its complete residual worsens. At t=.5 and .75 its corresponding lower-bound estimates are .00112446 and .00171583. Repairing the initial moment alone does not repair all-time compatibility or momentum balance.

The 96-point independent integral of phi dot analytic R agrees with the velocity-moment route within 7.44e-8 across the recorded cases. This is a quadrature cross-check, not interval proof. The harmonicity, convection contraction and exact cylinder norm are checked symbolically.

## Implemented solver versus proposed next step

Implemented: joint swirl-time and pressure correction, exact centrifugal nonlinearity, initial-flat swirl, frozen poloidal field/force, 25-time probe nullspace, bounded training-feasibility selection, independent full validation and global-moment diagnostics. The preconditioner retains initial swirl-pressure cross blocks. It is NOT exact variable projection, nor a full three-velocity-block optimizer.

Next construction should allow initial energy redistribution and poloidal dynamics while keeping E0=1, the original force/support and geometry constraints. Weak-moment compatibility must be controlled through time; otherwise locking initial data or fixing only one moment may give misleading progress. A pressure Schur-complement/variable-projection implementation with active pressure constraints is a proposed acceleration, not a completed feature. No fixed completion date follows from this small residual change.

## Persistence and checks

The complete user package contains correct raw Q2, J, J2, M0 and M1 fields, both J2/M1 reconstruction arrays, all fit registrations/source snapshots/checkpoints/histories, all six full reports, independent structure/dense/moment audits and peak refinements. The large deterministic preconditioner arrays are included in the package, not in GitHub.

GitHub tracks the executed solver, J2 frozen recipe, replay/audit/tests and this report. A malformed M1 transport copy was removed before PR publication; it was never used for fitting, validation or any reported result. The correct M1 recipe remains in the complete user package; GitHub does not claim to contain that control array. Its original raw SHA is 4bca9aa7c1ce6ccee0673a6e10f2936a667ae2e0d309d30beb6fb8ff6ed0a8ac. J2 raw SHA is 31961e6c7ff58358954967308fdaa3dcda3d52dfe05c796e35d898203bf0375d; modifier NPY SHA384af675a0b136660876a4056e2e3ba8e70f1ec8778b1327d97ae892b02fe1b4. The executed solver and uploaded J2 recipe Git blobs match local bytes exactly.

Source/package test receipts record the actual focused tests. Full historical tests, native MATLAB and Lean were not run for this round. Cloud status, if any, must be read separately; no cloud success is inferred. Original1e-3 target remains unmet. No default promotion or source-field/blow-up assertion.
