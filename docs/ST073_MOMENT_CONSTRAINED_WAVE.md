# Restore four mean moments during wave correction

The earlier momentum-aware candidate reduced complete momentum but failed
renewed compatibility: maximum actual integrated moment error was
583600.74, with 0/27 cone locations passing. It was not a valid recursive
step. Its complete diagnostic is preserved in
`wave_dynamics_mean_compatibility.json`.

The new solve enforces four moment equalities while eliminating linear
time-derivative/pressure controls. The right-hand sides retain the actual
wave's quadratic force; they vary with wave shape. The remaining linear
least-squares problem has rank 175. The outer search retains the five-node
flux targets and original growth/norm constraints. Reports are separate
from the previous unconstrained candidate.

| Check | Earlier candidate | Moment-constrained candidate |
|---|---:|---:|
| Training complete volume L2 | 2673580.61 | 2694421.28 |
| Independent actual complete volume L2 | 2645336.89 | 2661414.30 |
| Independent actual complete maximum | 2.00979e11 | 2.03734e11 |
| Actual order-96 moment max | 583600.74 | 2.54e-8 |
| Cone64 locations passing | 0/27 | 12/27 |
| Cone64 individual inequalities passing | 2/81 | 50/81 |

Restoring the moments increases independent residual L2 by about 0.61%,
retaining most of the earlier improvement. The assembled moment error is
5.73e-9; actual integrated replay gives 2.54e-8. Cone constraints remain
unsatisfied: their minimum margin is -12508.24. More passing locations
does not mean every margin improved. Full momentum remains far above
1e-3, and no time integration or scale recursion is accepted.

Applying the saved higher z20/r11 energy operators to this candidate gives
growth lambda +997.0958. Those operators apply because the instantaneous
mean velocity and wave geometry are unchanged; mode-0 additions here are
time derivatives and pressure. Positive finite-grid growth is supported,
but the imposed 1000 floor is not independently attained.

## Reproduction

After generating the original projection cache, run:

```powershell
python experiments/root_st073/wave_dynamics_codesign.py --moment-report experiments/root_st073/wave_dynamics_mean_compatibility.json --seed experiments/root_st073/wave_dynamics_codesign.json --output experiments/root_st073/wave_dynamics_moment_codesign.json --maxiter 80
python experiments/root_st073/wave_dynamics_replay.py --source experiments/root_st073/wave_dynamics_moment_codesign.json --output experiments/root_st073/wave_dynamics_moment_replay.json
python experiments/root_st073/wave_dynamics_mean_compatibility.py --candidate experiments/root_st073/wave_dynamics_moment_codesign.json --replay experiments/root_st073/wave_dynamics_moment_replay.json --output experiments/root_st073/wave_dynamics_moment_compatibility.json
python experiments/root_st073/cached_wave_growth.py experiments/root_st073/wave_dynamics_moment_codesign.json experiments/root_st073/wave_dynamics_moment_growth.json
```

The next inner solve must include the remaining cone inequalities as well
as moments. Its needed affine mean-control rows and quadratic wave-force
forms are being assembled separately; no completed cone-constrained solve
is claimed here.

## Derivative precision work

`supported_fourier_analytic_jets.py` computes compact Fourier velocity,
Cartesian gradient and vector Laplacian analytically, including cylindrical
connection terms. Factored cutoff derivatives avoid cancellation near the
support edge. Existing basis values agree to about 1e-11 or better; sampled
analytic divergence is at most 7.47e-9. These are basis checks, not full-field
accuracy or continuum-divergence certificates.

`wave_momentum_projection_analytic.py` reuses the stored mean jets and
replaces only the wave derivatives. Its difference from the older frozen
FD cache has maximum 6467.67 and volume L2 0.005137. This exposes a numerical
difference already exceeding the ultimate absolute tolerance. The new
cache is available for later optimization through `--projection-cache`;
it is not silently substituted into the reported moment candidate. Mean
derivatives still use finite differences, so full analytic accuracy remains
unfinished. Regenerable NPZ caches remain local and ignored by Git.
