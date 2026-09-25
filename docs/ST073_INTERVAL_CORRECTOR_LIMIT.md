# The current coefficient march is not a pulse amplitude inverse

The late-bridge two-harmonic curl patch at `k=11` loses direct interior-time
momentum control when its explicit coefficient step grows from `1e-9` to
`1e-8` (`ST073_WAVE_STEP_HORIZON.md`). We tested two local repairs at
`dt=1e-8` on the same three-knot moment-closed base:

1. `adaptive_bridge_wave_interval_corrector.py` forms a damped trapezoid
   endpoint equation for the spatial harmonic and mean-potential
   coefficients. With damping `0.25` and four updates, its relative
   fixed-point defects are `2.10`, `10.45`, `13.90`, and `9.33`; the final
   defect is `7.80`. Its first-interval direct midpoint momentum maximum
   is `5.93e11`, compared with `1.04e10` for the frozen wave. This is
   neither a converged implicit solve nor a useful corrected interval.
2. `adaptive_bridge_wave_two_stage.py` now accepts an explicit ridge
   parameter for the harmonic/mean slope fits. Increasing it from `1e-4`
   to `1e-2` leaves the second-interval direct midpoint maximum essentially
   unchanged: `3.28e11` becomes `3.27e11`. The first-interval Hermite
   midpoint remains about `7.50e10`.

These are finite collocation results, with a temporal difference stencil of
`dt/8`; they do not prove impossibility of another stable discretization.
They do show that a trapezoid correction and modest ridge adjustment do not
extend this ansatz to an appreciable pulse interval. Endpoint projected
residuals still understate the directly evaluated interior error.

The [OpenAI paper, Proposition 7.2, Equation (7.13)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
solves the transverse amplitude ODE along the pulse path with zero initial
data and reconstructs pressure from the normal constraint. Its hypotheses
include slow, transverse, and shell support with smooth zero extension;
Proposition 9.6 then recomputes the full residual after wave, stress, mean,
and moment corrections. Our polynomial-in-space physical-time Taylor fit
does not implement that inverse or those support estimates.

The next constructive step is a pulse-coordinate amplitude/pressure solver
for at least one nonzero harmonic, using the actual background/phase data.
It must verify the transverse constraint, pressure identity, support and
cutoff budget, then feed the resulting exact-curl field into direct
interior-time momentum checks. Only after a full correction cycle improves
the residual on one band should it be transferred to another dyadic scale.
