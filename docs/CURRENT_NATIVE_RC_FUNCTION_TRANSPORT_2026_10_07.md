# Original Rc source-function integral representation and target interface

Implementation: [665263bc](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/665263bc64dfea1d61f0f2b60be0f858d8299a34). This follows the checked centered-phase bound stage in [CURRENT_NATIVE_CENTERED_PHASE_CONDITIONING_2026_10_07.md](CURRENT_NATIVE_CENTERED_PHASE_CONDITIONING_2026_10_07.md).

The new `NativeRcFunctionTransport.build()` binds the accepted original full signed density graphs to **24 true cells / 17 charts**, constructs their exact Duhamel integral representations, carries all incoming histories, and forms all five Rc target functions and their first Z derivatives. These are source-bound function definitions. Original numerical integral values and five control functions have not yet been computed.

## Implemented function interface

Files use the `experiments/root_st073/lei_ren_part1_paper_compliant_` prefix:

- `current_native_Rc_function_transport.py`: typed function graph and optional explicit-oracle evaluator.
- `current_native_Rc_function_transport.json`: 1590 nodes, original geometry/source references, histories and five target roots.
- `current_native_Rc_function_transport_check.py` / `.json`: focused independent geometry, signed source, derivative and evaluator checks.

`FunctionRef` and `C1Function(value,Z)` are separate from logarithmic caps and signed enclosures. The new graph contains no density-cover coefficients chosen as defining values. Original source parameters are bound to the same native seed, selected constants, bridge/switch widths, sc and original mu. The accepted source graph hash and source family are carried through the evaluator contract.

For each native coordinate s, use the original offset `y(s)=log(R(s)/r_minus)` and `phi=frac(N*y(s))`, with one shared integer N>=160. The integrals are

`I_j(Z,N)=integral_left^right exp(-lambda_j*(y_right-y(s)))*delta_density_j(s,Z,N)*y'(s) ds`.

The full signed original density and its first Z derivative are used. The native Jacobian is applied once. The exact width is collected symbolically before converting to a graph, so tiny bridge/switch widths are not obtained by subtracting rounded huge endpoints. Each step is

`D_out=exp(-lambda_j*width)*D_in+I_j`.

The initial correction is the actual zero inlet. On the proved initial collar and original O3 power cells, the incremental source is exactly zero; incoming correction continues to propagate. Pressure has lambda=0, so its memory is retained. Original underlying velocity, absolute pressure P0 and P0_Z are not zeroed or redefined.

All integration endpoints, weights and global phase are Z independent. Thus `I_j_Z` integrates the full source derivative without extra boundary terms. The graph carries these derivatives through every incoming/outgoing step.

At the original Rc endpoint, obtain A and A_Z from the original O3 power E/E_Z function roots. With the same original positive mu, use

`r=(D_m/A,(D_k-A*D_m)/(mu*A^2),D_h/A,D_e/A^2,D_p/A^2)`.

For `C=D_k-A*D_m`, retain `C_Z=D_k_Z-A_Z*D_m-A*D_m_Z` and

`r_D_Z=(C_Z-2*(A_Z/A)*C)/(mu*A^2)`.

Both r and N*r have function roots. N dependence is inside the original inverse/modulation/source functions; the N^-1/N^-2 coefficient covers from earlier stages are not used as polynomial function coefficients.

## Evaluation contract and evidence

The optional `evaluate` API requires an injected oracle with `parameter`, `source` and `integrate` methods, matching the source family. Its mode must explicitly distinguish `original_source` from `synthetic_reference`; original mode also binds the accepted source graph hash. Missing oracle, a cap record or a cap coefficient fails closed. Numerical quadrature remains approximate unless the injected backend certifies it.

The present implementation **does not install that original numerical oracle**. The included signed variable-Jacobian quadrature reference uses a manufactured smooth source and is explicitly labelled synthetic; its error is approximately 5.74e-42. It tests integral/Jacobian/sign/derivative wiring and is not a physical-field result.

Focused checks passed:

- 17 radius maps and 17 Jacobians independently compared with the accepted original radius binder; 23 exact route joins and the original zero offset inlet.
- 210 full signed C0/Z density-root bindings, 240 inherited C0/Z rows and 15 exact quiet contribution rows.
- 105 independent integral C1 functions; 16 symbolic first-Z history/target/numerator identities and five exact target transformations.
- Explicit no-oracle and cap-as-function rejection; the labelled signed variable-Jacobian reference above.
- Producer 0.375s; focused checker 6.063s. Working and staged dependency audit: 1089 hashes. Read-only review: **GPT-5.6 Luna / max**, scoped PASS.

## Next agent tasks

- [x] **SOURCE1a function-integral representation:** original full signed source graph roots on the complete 24-cell route, exact native coordinate maps/Jacobians, common global phase, zero inlet and complete incoming memory. The checked artifact now provides the function roots; numerical original integration remains open.
- [x] **SOURCE3a exact target representation:** bind same-original Rc E/E_Z and mu, form the joint numerator before division, and emit r/r_Z/N*r roots. This is a function representation, not evaluated target numbers.
- [ ] **SOURCE2a original source oracle:** implement an oracle resolving `original_function_graph` nodes from live original analytic source providers. Resolve derivative leaves from their actual function recipes and source coordinates; never use stored interval coefficient midpoints. Bind source hash/family and one N. Preserve the actual nonzero E/V/b/p2, original log scales and exact phase origin.
- [ ] **SOURCE2b inverse and precision:** evaluate the original monotone phase inverse and A/B C0/Z source functions using the correlated Poisson denominator, stable small increments and exact flat branches. Reuse accepted actual point/cell backend where available. Cover signed r, r=0, active body/transition, flat support, phase endpoints and all 17 charts. Return directed enclosures or explicit unresolved status when precision is inadequate; never silently return zero or substitute a different family.
- [ ] **SOURCE1b numerical transport:** attach source-range/point/quad capabilities to the 105 integral pairs. Control radial and Z refinement, quadrature and precision independently. Carry microscopic widths and huge scales in the accepted factored representation; apply original dy/dcoordinate once. Evaluate at least one full original route and report rigorous error/remaining unresolved cells. A function AST alone does not complete this task.
- [ ] **SOURCE3b evaluated original targets:** produce same-original r/r_Z from the numerical transport, retaining correlations in C/C_Z and original A/A_Z/mu. Record values separately from range covers and error bounds. Do not infer a small divided row from separate small k/m bounds.
- [ ] **CONTROL1a exact repair function map:** consume the typed N*r target functions and original exact bump-integral matrix/Gram functions. Build `h_next=-B_exact(mu)^-1*(N*r+Q_exact(mu,h)/N)` and its complete first Z map; the interval inverse enclosures remain sidecar bounds, not selected matrix entries. Preserve the inverse-mu cross row. Distinguish control h from repair increment target -r.
- [ ] **CONTROL1b actual controls:** evaluate or construct the five control functions with a convergent source-bound iteration and quantified C1 error. Show the actual repair equations, function-domain derivative consistency and nonlinear cross terms. A control-ball bound or sampled vector alone is insufficient.
- [ ] **CONTROL2 original closure and N:** justify a compatible finite integer N without materializing an astronomical integer unnecessarily. Repair-local log conditions alone do not admit global N. Verify all five terminal functions and their first Z derivatives on the required domain, then Rc->2Rc and post-repair2Rc->Rb admission using the same changed field.
- [ ] **HIGH/OUTER/ENERGY/REC/WAVE/PHYS:** actual higher derivatives/seams, whole global cone and admissible stress, pressure/heat/physical tail energy, true n=1/n>=2 coefficient recursion and independent repairs, oscillatory stress cancellation, corrected Cartesian NS and measured scale recursion/material winding remain separate unfinished layers.

Scoped gate: `current_original_exact_N_dependent_Rc_source_integral_and_target_functions_defined`. Original numerical source oracle/integrals, actual controls, terminal closure, global finite N/cone, higher jets/heat/energy, genuine coefficient recursion and corrected NS remain **false/open**. The long-term goal remains active.
