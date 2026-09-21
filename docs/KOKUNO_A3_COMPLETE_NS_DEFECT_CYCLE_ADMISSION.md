# Kokuno Agent 3: complete-NS-defect / finite-cycle admission

This increment is a fail-closed bridge between the current correction-side machinery and any future real finite correction cycle. It is stacked on exact A3 #1045 (`3ec20a77b551be819a71e3308e9ecab2a2226406`).

The finite-stage engineering iteration is `u^(k+1) = u^k + delta_u_k`. The corrected Kokuno 2026-09-09 reader is structural provenance only; it is not independent PDE validation.

A real cycle is **not admitted** until one semantic candidate identity carries all of the following evidence: global velocity; matched global pressure; preregistered restricted forcing with evidence that it was not defined from the measured residual; a complete identity-bound NS defect; disjoint held-in and held-out sample sets; an independent second covariance column and bounded-inverse preflight; a Cartesian correction velocity; and an accepted current correction-cycle gain receipt.

Historical A3 #293 is retained only as the prior gain-guard mechanism. Its rejected one-column experiment is not rerun or widened. Historical A3 #282/#302 are retained as the second-column coverage/bounded-inverse prerequisites; their existence is not treated as evidence that the current candidate already has the required second column.

The current exact #1045 lineage has only the scoped path `nonlinear m=0 mean -> compact stress -> radial force partial_z sigma_1` through the unmodulated RF40 power-law endpoint X4. Its own truth boundary keeps global completion, matched pressure, restricted forcing, complete NS defect, Cartesian correction velocity, finite cycle, held-out normalized residual, same-protocol ST006 comparison, and PDE validation false. Therefore this admission gate must return `finite_cycle_admitted=false`.

The public materializer intentionally accepts no caller-supplied residual/defect arrays, pressure, forcing, correction, gain, held-in/out sample payload, or scientific threshold. Final normalized momentum and divergence gates remain fixed at `1e-3` and `1e-5` respectively. Scoped mean/stress/radial-force evidence cannot be relabeled as a complete NS defect, and residual-defined free forcing remains forbidden.
