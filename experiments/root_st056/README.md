# ST056 — joint residual correction and compatibility audit

Task #746. Start from complete ST054-Q2, source commit c77492a48e9c0f13d4d51244987c57c28519409b. No missing ST055 arrays are treated as recovered. The publication repository, defaults, physical thresholds and other agents are unchanged.

**Completed:** a frozen joint swirl/pressure candidate ST056-J2, a negative moment-balanced control ST056-M1, six full independent comparisons and new spatial/structural/weak-moment diagnostics. **Original full NS1e-3 target remains unmet.**

Read [RESULTS.md](RESULTS.md) for paired numbers and [HARMONIC_DIAGNOSTIC.md](HARMONIC_DIAGNOSTIC.md) for the stronger initial-data compatibility diagnostic. The latter is a post-freeze check, not a new optimization round or interval proof.

## Actual fit scope

J2 jointly adjusts 756 initial-flat swirl-time and 864 compact pressure variables. It retains the exact nonlinear centrifugal term, fixes poloidal velocity and original forcing, and locks the old core probe at 25 training times through a rank-12 numerical nullspace. It completed650iterations/966calls, reaching the iteration limit in152.03seconds, NOT optimizer convergence. The chosen feasible training step preserves sampled peak and every training-time MSE relative to parent.

M1 is a deliberately different control: initial raw swirl is multiplied by .958132801471448 before initial-energy renormalization; physical poloidal amplitude grows2.080673%. It then uses the same joint correction with an additional weak-moment penalty. It stopped at the210second budget after508iterations, not convergence. Its full independent L2 worsens; it is not promoted. The earlier J fit is retained as a training-only control in the archive.

All fields were frozen before both fresh seeds9205691/9205692, at2026-09-20T02:07:58.725762+00:00. No subsequent tuning or coefficient selection. Original nu=.01,t[.25,.75],compact u AND p,solenoidal bounded two-parameter force,E0=1 and the original full max/volume-L2 gates are unchanged. New optimizer penalties are autonomous, not source-paper constants.

## Reproduce without fitting

From repository root:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st056/test_joint.py experiments/root_st056/test_harmonic.py
python experiments/root_st056/replay_st056.py --id ST056-J2 --out outputs/ST056-J2 --seed 9205691 --validate --structure
```

The last command currently exits1 for the two failed momentum gates. It refuses a nonempty output directory. In the research repository, the parent is reconstructed through the existing ST054 recipes; the complete user package contains the original SHA-pinned raw parent. Frozen references verify u/p at1e-8. Reconstructed JSON metadata and file hashes differ; raw byte identity is not asserted for regenerated files.

`audit.py` computes new core/shear, dense residual and degree-2 weak diagnostics. `harmonic_audit.py` computes nested degree2-8 harmonic-gradient witnesses. Finite samples and converged quadrature estimates are not continuum acceptance certificates.

## Persistence

GitHub contains the executed solver, exact checksum-verified J2 correction, replay/audit/tests and completed reports. The complete user archive additionally contains correct M1 data/recipe, all original raw fields, fit registrations/source snapshots/checkpoints/history, deterministic preconditioner arrays, six full reports and supplementary diagnostics. The malformed M1 transport copy was removed before PR publication and never entered a computation; the correct local control remains available in the archive.

J2 original raw SHA31961e6c7ff58358954967308fdaa3dcda3d52dfe05c796e35d898203bf0375d.
M1 original raw SHA4bca9aa7c1ce6ccee0673a6e10f2936a667ae2e0d309d30beb6fb8ff6ed0a8ac.

Ten focused local tests passed6.12s. Full historical tests, native MATLAB and Lean were not run. Cloud checks have not been verified for this round. Software tests do not promote the scientific status. All pde/source-correspondence/paper-exact/blow-up flags remain false.
