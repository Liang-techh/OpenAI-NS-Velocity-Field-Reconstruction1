# Active project instructions

The user explicitly replaced the exact-reconstruction objective with constrained independent reconstruction. Read `docs/PROJECT_GOAL.md`, `docs/AGENT_TASKS.md`, and `docs/CURRENT_CHECKPOINT.md` first.

- Work in `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`. Old repository task IDs/claims are historical and do not authorize duplicate work here.
- Deliver a nontrivial candidate velocity family, constrained optimization, independent numerical validation, and 3–5 selected symbolic/formal structure checks. Exact original coefficients and the complete original proof are not acceptance requirements.
- Reuse existing working modules. Do not spend the project closing every old theorem/provenance adapter. Preserve their truthful labels without making them the new critical path.
- Define viscosity, domain, boundary/support conditions, forcing model, scales and validation thresholds before optimization. Prevent zero-field or amplitude-collapse solutions with explicit nontriviality constraints.
- Never choose an unconstrained force equal to the computed residual and call that an independent NS validation. Forcing must be fixed or restricted in advance, and its regularity/support checked separately.
- Keep optimization samples and independent validation samples/operators separate. Vary discretization, quadrature, truncation and differentiation error independently as applicable.
- Record paper-sourced constraints, autonomous choices, numerical observations and proved identities separately. A candidate is not an exact reproduction or a blow-up theorem.
- Claim one unowned task, work on an isolated branch, and record commit/PR, actual checks, limits and completion. DONE, accepted and merged are distinct states.
- Favor concrete implementation over repeated inventories. Run focused tests for changes; the full inherited legacy suite is optional unless changes affect that scope. Never fabricate test, CI or Lean success.
- Use Luna/max for bounded implementation where useful; root handles mathematical choices and integration. Do not spawn copies of Astra. Do not let a stalled worker delay delivery indefinitely.

Historical instructions under `docs/legacy_exact_reconstruction/` and historical manifests describe the previous project only. They do not override this goal.
