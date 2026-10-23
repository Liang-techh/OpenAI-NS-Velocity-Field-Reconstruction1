# Constrained Agent 9 external-method record

## CR-A9-002 — rank-revealing basis-capacity screening

Current active constrained work has grown to a large velocity/pressure coefficient family while independent held-out momentum residual remains O(1). Before adding another large correction block, a cheap rank/capacity screen can detect redundant residual-Jacobian directions and estimate whether the current linearized span can even address the observed residual channel.

### External candidate

- Source: `scipy/scipy`
- Screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- Public API: `scipy.linalg.qr(..., pivoting=True)`
- Method: QR with column pivoting / rank-revealing QR
- License: BSD-3-Clause
- Classification: **可直接迁移（public API only）**
- Migration scope: only call SciPy's existing public QR API, already permitted by this repository's `scipy>=1.10,<2` dependency. No SciPy implementation, LAPACK wrapper, or source code is copied.

### Repository adaptation

`constrained_basis_rank.py` normalizes residual-Jacobian columns before pivoted QR so the diagnostic highlights geometric redundancy rather than coefficient-unit scale. It reports numerical rank, pivot order, independent/dependent basis labels, normalized condition number, and an unbounded first-order least-squares projection of the residual onto the selected span.

The projection is intentionally diagnostic only. It does **not** enforce coefficient bounds, nonlinear feasibility, nontriviality, support/boundary constraints, the restricted forcing contract, or held-out validation. Its projected residual must never be reported as a PDE validation result or acceptance evidence.

### Why this method is useful now

The active constrained branch records a current `coupled_joint` family with 117 velocity and 27 pressure coefficients while the finest held-out momentum maximum is still about 1.0 against the preregistered 1e-3 threshold. A rank-revealing screen can answer a narrower question before another optimization cycle: are proposed correction directions genuinely independent in residual space, or are we paying nonlinear optimization cost for nearly duplicate directions?

This does not replace Agent 1's representation work, Agent 2's bounded pressure projection, Agent 3's independent refinement audit, or Agent 7's channel-capacity diagnostics. It is a generic pre-fit filter that can consume their residual Jacobians without modifying candidate state.

### Verification and truth boundary

Focused regression tests cover duplicated/zero directions, an exactly representable manufactured residual, an orthogonal residual floor, and malformed/non-finite input rejection. Local isolated execution before upload: `3 passed in 0.25s`; `py_compile` completed with exit 0.

No domain, viscosity, time interval, support/boundary condition, forcing convention, nontriviality normalization, residual norm, validation seed, or acceptance threshold is changed by this increment. No `u -> 0` path and no `f = R(u,p)` path are introduced. A green unit test or CI run is not PDE evidence.
