# ST050R — recovered, reproducible pressure/morphology experiment

This is a NEW bounded experiment, not a claimed recovery of the unavailable ST050 C/P/E/PE/PC arrays. The earlier WORK_IN_PROGRESS note is preserved separately. Neither a lower training loss nor a smaller pressure-Poisson residual establishes full NS acceptance.

## Starting point and identity

Task #408. The branch starts from `71b5a6cfff866fbca9c2127c47e39b54cf8ae969`, which contains the ST048 source stack and the original ST050 WIP note. The actual starting mathematical field is ST048-S, reconstructed from immutable ST047-E and the original SHA-checked 398-modifier recipe. All eight published parent velocity AND pressure references match locally with maximum difference 0.0.

Reconstructed parent raw JSON SHA256:
`13127759debcf40b3f3d6c3627747ace67f1471b825b2b0e8c53ae79a9911344`.
Its metadata differs from original ST048-S raw SHA `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`; byte identity is not claimed.

## Unchanged original physical contract

nu=.01; t=[.25,.75]; physical R3, evaluation box[-2,2]^3; smooth compact velocity AND pressure in r<2, |z|<2; original two-parameter divergence-free curl force with a,c in[0,10]; E(.25)=1 and the original energy range, parameter bounds, support, divergence, core signs and scaled-core gate. Both original full-vector momentum max and spatial volume L2 must be <=.001. Volume L2 is sqrt(64*mean(|R|^2)) at each time, not training MSE or a space-time average. No arbitrary residual-defined force or zero-amplitude solution.

Both fields use the existing 2594-stored-coefficient asymmetric axisymmetric family and 398 reduced modifiers. No added global reflection restriction and no mathematical formula changes in inherited evaluators.

## New mechanisms, separate from scientific acceptance

`pressure_morph.py` adds the exact pressure-Poisson consistency expression

```
div(R) = tr((grad u)^2) + Delta(p)
```

for this divergence-free velocity and force. With u=(xA-yB,yA+xB,C), s=r^2:

```
tr((grad u)^2)=(A+2sAs)^2+A^2+Cz^2-2B(B+2sBs)+4sAzCs
on the axis: div(R)=6A^2-2B^2+4p_s+p_zz
```

At exact R=0, strictly positive radial and axial pressure curvatures require B^2>3A^2. This is a pointwise necessary condition, not family-wide infeasibility, not a bound on div(R) from a sampled max of R, and not a replacement acceptance metric.

The actual vorticity squared is used to form quadratic moment matrices: normalized axial second/fourth moments must remain at least .999 of the parent's values; the normalized radial second moment is bounded by 1.01 times parent. Initial-energy normalization cancels in these ratios. These finite quadrature safeguards are autonomous, not source-calibrated shape matching or a certificate for all vorticity tails.

Original signed-shear preservation .995, the 3423 support-edge training probes and inherited core/pressure/nonpressure-dynamics constraints are reused. All derivatives account for energy normalization. A fixed numerical whitening map changes solver coordinates, not physical acceptance norms.

## Two deliberately different fits

**ST050R-C:** Poisson weight .003; inherited autonomous core-anchor tolerance .08, per-training-time L2 cap1.0, pressure ratio.99, profile ratio.995, nonpressure axial ratio.98. Completed163 iterations/190 evaluations, numerical optimizer success. Best feasible training constraint minimum -1.00e-12. Training Poisson RMS .03707165 -> .03245941. This is not independent NS validation.

**ST050R-P:** stronger Poisson weight .03, nonnegative signed axial-pressure target on the training core, and axis swirl/strain ratio>=1.75. The previous AUTONOMOUS anchor is explicitly relaxed from .08 to .60 and the previous TRAINING-only per-time L2 cap is removed. Original scientific/physical gates are unchanged, but this is NOT preservation of every auxiliary constraint. Completed192 iterations before480-second wall budget, NOT optimizer convergence; best feasible training minimum -1.36e-11. The reported271 calls at budget stop count the inherited objective history, not a recovered SciPy final nfev. Training Poisson RMS .03398060.

Training:32x48 Gauss space,13 Gauss time nodes+endpoints;17 structure constraint times; morphology32x48 at five times; Poisson18x26 at seven times. Morphology/shear/edge parameters and all budgets were recorded before each fit. C's initial launch was interrupted before a candidate; the separately logged retry completed. C ran the same default numerical formulas before the optional axis/P switches were added; its executed source snapshot is retained in the archive. Checkpoints are written every five iterations and the feasible selection tolerance is1e-7. No claim of global optimality.

Both chosen fields were frozen BEFORE fresh independent seeds9175091 and9175092. No further parameter tuning is allowed from those observations in this experiment. See `results.json` for executed comparisons and the user archive for the original full reports, additional diagnostics, freeze receipt and command exit records.

## Reproduction without fitting

From repository root:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st050r/test_recovery.py
python experiments/root_st050r/replay_recovery.py --id ST050R-C --out outputs/ST050R-C --seed 9175091 --validate --structure
```

Use `--id ST050R-P` to inspect the stronger structural control. The script refuses a nonempty output directory and returns the original scientific gate result. A scientific rejection is not a software crash. The recipe checks modifier SHA256, parent, finite bounds, unsupported scientific flags, and eight frozen u/p references at1e-8. GitHub uses the ancestor ST048 recipe when raw parent is absent; the user archive contains the SHA-pinned reconstructed parent. Regenerated metadata and file SHA may differ, while mathematical fields must match.

To repeat the bounded fits with an actual saved ST048-S parent:

```bash
python experiments/root_st050r/pressure_morph.py --parent artifacts/research/ST048-S/candidate.json --out outputs/refit_C
python experiments/root_st050r/pressure_morph.py --parent artifacts/research/ST048-S/candidate.json --out outputs/refit_P --id ST050R-P --pp .03 --anchor .6 --cap -1 --axis-ratio 1.75 --pressure-target 0 --maxiter 240
```

Wall budgets and platform changes can alter retrained parameters; frozen replay is the primary reproducibility path. No external numerical solver or source-paper coefficient recovery is claimed.

## Delivered evidence boundary

Local final seven tests passed43.67s, warnings-as-errors, including finite differences of the complete Cartesian residual, Poisson/moment modifier derivatives, symbolic trace/axis/divergence identities, actual energy/support, both recipes and claim/checksum mutation rejection. Five uploaded source/recipe Git blob identities were checked against actual local files. Full inherited local suite and Lean were not run.

GitHub contains readable new source, both frozen modifier arrays/reference values, result summary and replay CI. The accompanying archive additionally contains raw parent/C/P JSONs, actual registrations, interrupted-start and completed-run logs/checkpoints, all six original validations and separate structure/dense audits. It is not asserted that every old historical experiment has been recovered or uploaded. No original source/default/config/threshold/other-agent schedule is modified.

Original local raw child SHA256:
- C: `eed7321a3d44e7a6aff5b99432cfa73215546a36df263ea8968cdc95663fc8b9`
- P: `f5b5a2869991387e503a522b599af86ea814ba02bdae2f8b56fd16674f7dab19`

Scientific status is determined by the independent reports, not this file's existence. `pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `blowup_proved=false`. No automatic default promotion.
