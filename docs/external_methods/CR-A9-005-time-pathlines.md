# CR-A9-005 — time-dependent pathline integration

## Classification

**Directly reusable public API**.

## External source

- repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- public API: `scipy.integrate.solve_ivp`
- relevant source: `scipy/integrate/_ivp/ivp.py`
- license: SciPy BSD-3-Clause

SciPy's public `solve_ivp` interface solves first-order initial-value ODE systems and supports explicitly requested output times. This repository already depends on `scipy>=1.10,<2`, so this increment adds no dependency and copies no SciPy implementation code.

## Why this lowers the current blocker

The final deliverable is a time-dependent callable `velocity(x,y,z,t)->[u,v,w]`. Frozen-time streamlines do not show how material trajectories move through the evolving field. A pathline solves

`dX/dt = velocity(X,t)`

and therefore gives a direct time-evolution visualization channel from the same public velocity callable, without re-entering the field formula in MATLAB/Python or changing the candidate.

This is complementary to the open frozen-time streamline fingerprint work: streamline geometry samples one fixed time, while pathlines integrate the explicitly time-dependent velocity.

## Migration scope and differences

Only SciPy's installed public `solve_ivp` API is called. No integrator source, tableaux, dense-output implementation, or adaptive-step logic is copied.

The repository wrapper deliberately narrows the interface to:

- seed positions shaped `(n,3)`;
- strictly increasing requested sample times;
- one batched state containing all seed particles;
- optional declared allowed time interval;
- fail-closed finite shape checks on every velocity evaluation;
- returned positions and speeds with explicit visualization-only truth labels.

## Truth boundary

Pathline integration is visualization/kinematics evidence only. A stable, converged, or visually similar pathline is not evidence that the velocity satisfies Navier--Stokes, is divergence-free, matches any hidden OpenAI numerical field, reproduces paper coefficients, or proves blow-up. This increment changes no domain, viscosity, forcing convention, nontriviality normalization, PDE residual definition, training/validation split, or acceptance threshold.
