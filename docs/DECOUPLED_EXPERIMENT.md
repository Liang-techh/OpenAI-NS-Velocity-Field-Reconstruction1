# Decoupled coefficient experiment

For each velocity shape, solve a five-column bounded linear least-squares
problem for three pressure coefficients and two force coefficients. Then
optimize only nonlinear velocity parameters. The force family and all bounds
are unchanged. This is variable projection using the existing candidate,
not a reproduction of the neural-basis Decoupled-DFNN algorithm.

Reproduce (installed package or PYTHONPATH=src):

```python
from openai_ns_reconstruction.constrained_optimize import run
run('configs/constraints_v3.json', 'artifacts/constrained/decoupled_v3', decoupled=True)
```

Then run the independent validator with the saved candidate, training artifact
and constraints_v3.json. Results: 660 outer function evaluations, ftol
termination; training loss 0.255627. Independent maximum sampled residual at
finest step: 4.036601, versus 4.115165 for the joint v3 run. Threshold 0.001
still fails. Inner solves consume additional work beyond the outer call count;
no equal-wall-time speed comparison is claimed. The test checks that fitting
pressure/force decreases a fixed velocity's residual without changing velocity
or escaping coefficient bounds (2 related tests passed in 0.57 s).

## Relevant external work examined

- [nsblowup](https://github.com/CokieMiner/nsblowup): finite-scale surrogate,
  pseudo-spectral solver and frozen forcing. Its forcing is manufactured from
  a reference-field residual; that cannot satisfy our restricted-force objective
  by itself. Potential reuse: numerical operators and resolution diagnostics.
- [Decoupled-DFNN](https://arxiv.org/abs/2603.17906): vector-potential velocity,
  separated pressure/velocity subproblems and Gauss-Newton linearization.
  Motivates investigating separation, but this experiment does not implement
  its neural basis or reproduce its reported results.
- [Axisymmetric Euler profiles](https://arxiv.org/abs/2201.06780): neural
  self-similar profile search; different equations from our viscous NS problem.

These sources were read, not independently reproduced. This experiment is
original local implementation; no external code was copied. The modest result
supports investing in a richer velocity representation rather than expecting
pressure fitting alone to solve the remaining residual.
