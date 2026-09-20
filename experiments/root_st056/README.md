# ST056 — joint momentum correction and initial-data compatibility

Task #746. This round is a new bounded experiment from the complete ST054-Q2 checkpoint, not a recovery of missing ST055 arrays. The publication repository and defaults are unchanged.

## Completed fits, before independent validation

The main new field is **ST056-J2**. Swirl-time and compact pressure are optimized together, including the exact nonlinear centrifugal change in radial momentum. The poloidal field and original restricted force remain fixed. Swirl corrections vanish at t=0.25, preserving the entire original initial velocity. The original core probe is additionally held fixed at 25 training times using a rank-12 numerical nullspace. That is a finite-probe preservation device, not a continuum constraint.

ST056-J was an earlier control with weaker penalty settings. Its unscaled final iterate was infeasible; a training-only scaled direction supplied a small feasible step. It was not chosen as the main field.

ST056-J2 completed 650 iterations / 966 solver calls in 152.03 seconds, reaching the iteration limit, NOT optimizer convergence. The selected feasible training direction reduces the broad training-pool maximum from 0.03799680163 to 0.03670154081 and weighted training MSE from 3.341090492e-5 to 3.303192771e-5. Every training-time MSE remains below its parent value. These are training observations, not independent NS acceptance.

A separately registered **ST056-M1** control starts from a moment-balanced initial field. The parent is not unchanged Q2: its initial poloidal amplitude is 1.02080673083 times Q2, and relative swirl factor is changed by 0.958132801471448 before energy renormalization. Initial energy remains one. M1 then performs the same joint correction with an additional weak-moment penalty. It stopped at the registered 210-second budget after 508 iterations; numerical convergence is not claimed. Its training L2 is worse than Q2. It must not be promoted merely because an initial integral condition is repaired.

All fits finished before freeze **2026-09-20T02:07:58.725762+00:00** and before either fresh seed **9205691 / 9205692**. No subsequent fitting or coefficient selection is permitted from these holdouts.

Frozen original JSON SHA256 identities:

- Parent Q2: `d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5`.
- ST056-J2: `31961e6c7ff58358954967308fdaa3dcda3d52dfe05c796e35d898203bf0375d`.
- ST056-M1: `4bca9aa7c1ce6ccee0673a6e10f2936a667ae2e0d309d30beb6fb8ff6ed0a8ac`.

## Original contract

nu=0.01; time [0.25,0.75]; R3 with evaluation box [-2,2]^3; smooth compact u AND p in r<2, |z|<2; original two-parameter divergence-free curl forcing; initial energy one; original coefficient, energy, support, divergence, core and full-momentum gates. Full-vector sampled maximum and fixed-time spatial volume L2 must BOTH be <0.001. The latter is sqrt(64*mean(|R|^2)), not RMS or a time average.

New training guards and weights are autonomous, not source-paper constants. They preserve pressure directions, core swirl within a 3% band and the existing single-probe drift cap. They do not replace new off-grid audits or the original independent Cartesian validator.

## Why initial data must eventually change

For smooth compact incompressible u, compact p, and the independently prescribed divergence-free force, phi=(x/2,y/2,-z) gives the exact weak identity

```
Integral phi dot R = Integral [uz^2 - (ux^2+uy^2)/2] = D(t).
||R||_L2 >= |D(t)| / sqrt(88*pi/3).
```

The denominator is the exact L2 norm of phi on the support cylinder. The pressure term cancels because div(phi)=0; time derivative, viscosity and force cancel because phi is a gradient of a harmonic quadratic and u/f are divergence-free and compactly supported.

Q2 at t=0.25 gives a converged quadrature estimate D=-0.0403498178682 and L2 lower-bound estimate **0.00420325689**. The 48/72/96 quadrature values agree closely. These numerical integrals are NOT interval-certified bounds, and this is not an impossibility theorem for the candidate family. It does identify a structural limitation of preserving this entire initial field forever. J2 is a valid intermediate experiment, not an end-to-end route to true L2<0.001 with unchanged Q2 initial data.

At E(0.25)=1, the exact zero-defect condition would require Integral uz^2=2/3. The M1 control explores this necessary initial compatibility without changing viscosity, support, force or normalization. Necessary compatibility does not guarantee small momentum residual.

## Verification and publication boundary

Seven focused local tests passed, covering original-evaluator parity, the exact joint nonlinear update, parameter derivatives, initial/poloidal/force preservation, independent Cartesian residual, energy/support, and four symbolic structure identities. Full legacy tests and Lean were not run.

Independent full-protocol comparisons and fresh spatial/structure audits are being completed after the freeze. Their actual results will be recorded separately; no cloud run or PDE success is claimed here. Complete raw fields, source snapshots, training histories and reports are retained for the user research package. This README alone is not a recovery recipe or proof certificate.

`pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `blowup_proved=false`.
