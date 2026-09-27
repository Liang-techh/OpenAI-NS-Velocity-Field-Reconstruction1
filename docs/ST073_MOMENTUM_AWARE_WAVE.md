# Wave shape optimized against complete momentum

The new construction retains the ten scalar flux targets at the existing
five nodes while minimizing complete instantaneous momentum over the
physical patch. It includes all wave self-interaction terms. The fixed
180 potential time-derivative/pressure columns are eliminated by weighted
SVD projection; 54 real wave-shape coefficients remain to optimize.

The cached representation is exactly quadratic in those coefficients:
`R(x) = R0 + L x + (B x)(A x)`. It reproduces the existing frozen physical
finite-difference residual to relative volume L2 1.38e-10, but its absolute
L2 difference is 0.001384 and maximum difference is 98.0. This is useful
relative consistency at large residuals, not evidence for the final
absolute 1e-3 gate. Two analytic Jacobian checks give relative differences
2.95e-10 and 4.69e-12.

The first SLSQP attempt left the feasible set and retained the initial
candidate; it is archived separately. A subsequent feasible descent search
retracts each step onto the ten flux equalities and the original matrix
growth-lambda=1000 surface, checking the original coefficient norm bound.
It accepted 17 steps and stopped when none of the tested step lengths
improved the objective. This is a feasible local improvement, not a
converged or global optimum.

## Physical replay result

The order-9, 9072-point training grid is the former holdout, explicitly
reused as training. The new wave is replayed as an actual Cartesian field
on the separate order-13, 18720-point grid, with all time, diffusion,
advection and pressure terms and zero external forcing in this diagnostic.
The fixed physical patch volume is 2.13613448e-8.

| Complete momentum metric | Previous dense tangent | New wave shape |
|---|---:|---:|
| Training volume L2 | 7454759.29 | 2673580.61 |
| Independent actual volume L2 | 7376639.77 | 2645336.89 |
| Independent actual maximum | 5.64090e11 | 2.00979e11 |

The independent volume L2 and maximum both fall by about 64%. The new
residual squared L2 remains dominated by mode 0 (63.76%) and mode 2
(35.53%). These values remain vastly above the requested tolerance.

The normalized five-node flux mismatch is at most 1.79e-12. The original
growth matrix still gives lambda approximately 1000. However, direct
production/dissipation on the order-9 training grid gives lambda -8936 for
the new wave, and -9971 for the original wave. The source's positive
higher-quadrature result and this negative result initially disagreed.
`wave_growth_consistency.py/.json` resolves the discrepancy: the source
z16/r9 grid reproduces its stored operators to approximately 1e-17 relative
error, and the higher z20/r11 grid gives lambda +997.125 for the original
wave and +997.096 for the new wave. The coarse z9 integration changes the
large production and dissipation terms by about 1%, reversing the sign
of their small difference. Higher quadrature supports positive growth,
but still does not independently meet the imposed 1000 floor or establish
continuum convergence. Higher-grid operators are saved for reuse.
The fitted wave time derivative has positive energy rate
419.92, which does not itself prove the unforced wave energy balance.

Mode-0 time and pressure corrections also change mean compatibility. The
new field has not passed renewed moment/cone checks, continuous-time
matching, global-domain residual checks, or a scale transition. All
acceptance and recursion flags remain false.

`constrained_tangent_projection.py` now supplies an equality-constrained
variable projection for the next update: the linear time/pressure solve
can enforce four physical moment equations with a quadratic wave-dependent
right-hand side. Its algebra was compared with a direct KKT solution and
its variable-target derivative with finite differences. This is operator
verification, not evidence that the physical candidate satisfies moments.
The optimizer accepts optional `--moment-report`, `--seed`, and `--output`
arguments; replay accepts `--source` and `--output`, so subsequent
candidates can retain prior reports unchanged.

## Reproduction

From the repository root, with NumPy and SciPy available:

```powershell
python experiments/root_st073/wave_momentum_projection.py
python experiments/root_st073/wave_dynamics_codesign.py --method retracted --maxiter 80
python experiments/root_st073/wave_dynamics_replay.py
```

The adjacent JSON reports and frozen source coefficient snapshots are
tracked. The 67.6 MB NPZ cache is generated locally by the first command
and ignored by Git; later optimizations reuse it without assembling jets
again. The independent grid is saved in `full_wave_dense_tangent.json`.
The method follows the paper's need to retain complete wave covariance
and nonlinear remainders and to revisit mean corrections after wave
updates; this finite numerical search does not implement or certify the
paper's full construction.

Reference: [Sections 7–9 of the supplied OpenAI-hosted manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
