# Simultaneous swirl, poloidal and pressure fit

Exposed the existing90 swirl coefficients and27 poloidal coefficients in one117-column velocity basis, fitted jointly with27 pressure coefficients. The previous anchored candidate is exactly embedded by moving its90 inner-swirl coefficients into the leading block and zeroing the nested base coefficients. A regression checks velocity/pressure equivalence and serialization. Three focused coupled/poloidal tests passed.

The same48x48x9 plus2048 uniform training set and fourth-power residual objective were used. All144 coefficients remain in[-1,1]; force, initial/core preservation and acceptance thresholds remain unchanged. Analytic parameter derivatives and cached energy are reused.

29 function calls reached ftol. Training objective decreased to .01332709304. Independent standard-step momentum maximum changed1.00633693 ->1.00221830, a small gain. Core/energy probes pass with energy[.68503645,1.0] and unchanged core drift. Numerical divergence at the standard .005 step is .0036918883 and fails1e-5. No fine-step results have yet been generated for this exact candidate; the previous candidate's refinement must not be inherited. Full PDE acceptance remains failed.

Reproduce with constrained_poloidal_optimize.run(output='artifacts/constrained/coupled_joint',coupled=True), then the independent validation CLI with configs/constraints_coupled.json. Candidate and training/validation artifacts are committed. The scheduled-agent queue in docs/SCHEDULED_AGENT_TASKS.md defines follow-up work and completion evidence.
