# Current native O3 explicit cylindrical spatial jet - passed conditionally (2026-10-27)

This is a source-only metadata record for `native_O3_explicit_cylindrical_spatial_jet`. It records the explicit streamfunction/cylindrical formula packet, six physical cylindrical direction-partial entries per branch, nine Cartesian gradient entries per branch, pressure-gradient, curl and divergence expressions, and native-x/physical-volume speed and vorticity integrands. It does not claim evaluated numbers, a global Navier-Stokes closure, a phase/axis extension, or a physical norm.

The capsule is a literal source-data wrapper assembled by static extraction. Its `status` records static literal extraction only, `metadata_only=true`, `mathematical_functions_called=0`, and `new_recipe_executed=false`. The kernel consumes the frozen accepted exterior, radial, flux, and closed-parameterization data without replaying an ancestor helper family. Physical azimuth remains distinct from native fast phase. The source-jet domain records `off_phase_seams=True` and `conditional_C2_regular_source_functions=True`; this is not a seam-extension claim.

The conditional closed-parameterization frame flag is `True`; frame proof count is `1` and overall bundle proof count is `1`. Actual source-Z jet admission is `False` and actual source-parameter admission is `False`. The source-Z1/Z2 contract is consumed only conditionally: `True`.

Formula packet:
- `d=1-Z^2`
- `L=1-delta*Z^2`
- `a=1-delta`
- `lambda=sqrt(tau/d)`
- `r=lambda*sqrt(2R)`
- `kappa=lambda^a`
- `C=Pstar*lambda^(-1-delta)`
- `B=Pstar*R/(rL)`
- `S=d*m_Z+a*Z*m-2*Z*V`
- `u_r=-B*S`
- `u_theta=C*E`
- `u_z=C*V`

Physical direction packet: `D_r=(2R/r)*D_R`, `D_z=d/(kappa*L)*(D_Z-(2RZ/d)*D_R)`, fixed physical time `True`. The pressure receipt is `p=C^2*(P0+p_hist)` with P0 entering once; `∇p` is the gradient of `C^2*(P0+p_hist)` with the proper C^2 scaling. Curl receipt: `["omega_r=-u_theta_z", "omega_theta=u_r_z-u_z_r", "omega_z=u_theta_r+u_theta/r"]`. Divergence receipt: `u_r_r+u_r/r+u_z_z`. The Cartesian basis derivative receipt is `{"D_theta_u": "-v", "D_theta_v": "u", "D_theta_w": "0"}`.

Both background and modified branches are represented by saved scalar-program references. Physical compositions are recursively bound inside all native and saved source-function references. Source graph IDs are references only, not numeric source values. The closed frame keeps `core_logC` as an unselected opaque fixed input and does not materialize the opaque N atom. Source-provided D_y, m_y, and other source-port derivatives remain data-only expressions, not values; unknown partials are rejected. The source ports and radial targets are `["native_source_function_port", "saved_radial_derivative_port", "saved_flux_source_differential_port"]` and `["E_background", "V_background", "E_N", "V_N"]`.

The publication scope is conditional on the frozen capsule and saved result. Conditional flags are `{"conditional_cartesian_spatial_jet": true, "conditional_closed_parameterization_frame": true, "conditional_curl_and_divergence_formula": true, "conditional_cylindrical_velocity_formula": true, "conditional_pressure_gradient_formula": true, "conditional_source_Z1_Z2_contract_consumed": true, "conditional_volume_integrands_defined": true}`. Shape summary is `{"branches": ["background", "modified"], "cartesian_gradient_entries_per_branch": 9, "cylindrical_partial_entries_per_branch": 6, "cylindrical_velocity_components_per_branch": 3, "physical_Cartesian_function_wrappers_per_branch": {"background": 19, "modified": 19}, "physical_direction_operators": 2, "proof_rows": 9}`. Scalar-program and callback counts are `{"callback_rows": 9, "journal_formula_proof_rows": 9, "scalar_program_nodes": 577, "scalar_program_roots": 14, "scalar_program_source_bindings": 1}`. The primary volume receipt is `saved r*J2` from `accepted_exterior_payload.physical_volume.profile_R_Z_theta_absolute_measure`, with profile measure `physical_expr:201`; native-x contributes the extra R exactly once. No numerical norms, global seam/axis admission, time recursion, five-moment closure, or NS acceptance follows from this leaf.

Geometry data are reused from the accepted payload: `geometry_reused=True`, `geometry_dual_calls=0`, `old_math_calls=0`, and `source_function_port_calls=0`. These are metadata counters, not additional proof calls.

Actual saved counts: `{"N_materialization_calls": 0, "failed_family_calls": 0, "new_cartesian_composition_family_calls": 1, "new_cylindrical_jet_family_calls": 1, "new_explicit_formula_family_calls": 1, "new_mathematical_proofs": 1, "new_source_identity_family_calls": 1, "numeric_evaluation_calls": 0, "old_math_calls": 0, "old_provider_calls": 0, "saved_flux_source_differential_port_calls": 0, "saved_radial_derivative_port_calls": 0, "source_function_port_calls": 0}`.

Actual saved flags: `{"actual_incoming_x_independence_proved": false, "actual_parameter_admission": false, "actual_source_Z_jets_admitted": false, "conditional_cartesian_spatial_jet": true, "conditional_closed_parameterization_frame": true, "conditional_curl_and_divergence_formula": true, "conditional_cylindrical_velocity_formula": true, "conditional_physical_direction_operators": true, "conditional_pressure_gradient_formula": true, "conditional_source_Z1_Z2_contract_consumed": true, "conditional_volume_integrands_defined": true, "full_goal_complete": false, "global_axis_regularity_proved": false, "global_phase_seam_regularity_proved": false, "native_field_derivative_identity_proved": false, "numeric_evaluation": false, "physical_layer_completed": false}`.

Attempt status: `passed`. Nested failure, if any, is preserved exactly: `null`.

Independent review result SHA: `a9b74a4cea64e05fe79e2cfc971798576e25f9fc60d5785a488e12799aff9991`.

Journal format: `jsonl_compact_one_object_per_line`; complete objects: `11`; callback objects: `9`. The compact JSONL stream is parsed as a durable object stream; no post-save parser failure is fabricated.

The next frozen task is `source_native_O3_source_value_and_axial_jet_admission`. Its unexecuted plan is retained as `C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work/native_O3_source_value_and_axial_jet_admission_plan_2026_10_27.json` (SHA `8091f1a1b016a40d9340db0ba1688cc15a3424472b8278144ff19f3fb505a7a9`). This leaf does not execute or complete that next task.

The additional static value-anchor trace is retained losslessly in this report as `native_O3_core_logC_value_anchor_trace_2026_10_27.json`. It distinguishes the degenerate encoded `selected_logCstar` interval, the unmaterialized exact source anchor and the provisional root-interval center. `attach_oracle` supplies identity metadata, not point values. The next consumer must record a legitimate source-parameter selection and a new independent value-side admission; no bound, midpoint or provisional center is promoted to an actual function value by this leaf.
