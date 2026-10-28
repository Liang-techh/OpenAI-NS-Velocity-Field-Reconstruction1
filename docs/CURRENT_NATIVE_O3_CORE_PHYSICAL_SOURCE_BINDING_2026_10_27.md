# Current native O3 core physical source binding - saved metadata (2026-10-27)

This is a publication metadata record for native_O3_core_physical_source_binding. The saved consumer is metadata_only=false, source_only=true, numeric_evaluation=true, numeric_values_computed=false and current_leaf_numeric_evaluation=true. It evaluates only the Z=0 core slice: four rho values, four tau values and two quarter-turn theta cases (32 packets). physical_layer_completed=false, full_physical_field_materialized=false and full_goal_complete=false remain explicit.

Conditional on the authenticated saved result, source/operator identity and mathematical proof acceptance are true. Review receipts are metadata-only and independent; they do not add mathematical calls.

The saved coefficient data define analytic functions; they are not selected point substitutions. With lambda^2=tau/(1-Z^2), kappa=lambda^(1-delta), R=rho*epsilon_core, r=lambda*sqrt(2R), z=kappa*Z and psi=kappa*R*M, the local records use ur=lambda^-1*sqrt(R/2)*Q, utheta=lambda^(-1-delta)*sqrt(2R)*F0*Phi, uz=lambda^(-1-delta)*U and p=lambda^(-2-2delta)*P_profile, where P_scaled=epsilon_core*P_profile. The rotation is u=ur*cos(theta)-utheta*sin(theta), v=ur*sin(theta)+utheta*cos(theta), w=uz.

The global Q identity Q=-[d M_Z+(1-delta)Z M-2Z U]/L is retained as scope context only; this leaf consumes saved Z=0 Q, Qrho, QZ and ordinary source jets, then newly evaluates physical cylindrical jets and common-source growth without rerunning an ancestor symbolic pullback. Axis, increment and source relations remain admitted data, and no inherited global error term is removed. F0=exp(logF0_Z0) is a strict-positive formal atom: formal utheta, omega and logabs records remain distinct from ordinary numeric ur, uz and p, and the exponential is not materialized. The divergence interval is a numerical diagnostic, not an overlap proof. Common-source ratios and positive delta-excess certificates are not material trajectories, observed growth, or n-time recursion.

The source capsule carries frozen literal saved source DATA with preserved guards and source-result hash. Capsule scope metadata is {"grid": {"Z": "0", "rho": ["1/16", "1/4", "1", "3"], "tau": ["1", "1/4", "1/16", "1/64"], "theta": ["0", "pi/2"]}, "physical_contract": {"F0_positive_exponential_materialized": false, "core_U_has_extra_Pstar": false, "epsilon_core": "1/Lambda", "evaluated_scope": "same-source Z=0 core physical values and first spatial jets; not global construction", "global_Q": "-[d*M_Z+(1-delta)*Z*M-2*Z*U]/L", "map": "lambda=sqrt(tau/(1-Z^2)); kappa=lambda^(1-delta); r=lambda*sqrt(2R); z=kappa*Z", "ordinary_coordinate": "rho=Lambda*R", "pressure_unit": "P_scaled=epsilon_core*P", "source_comoving_points_are_material_trajectories": false, "streamfunction": "psi=kappa*R*M", "tau": "1-t", "velocity_source": "u_r=-psi_z/r;u_z=psi_r/r;u_theta=lambda^(-1-delta)*sqrt(2R)*exp(logF0)*Phi", "viscosity": "1"}, "runtime_model_metadata": null, "source_projection_keys": ["ordinary_jets", "local_core_recovery", "mean_minus_axial", "majorant_records", "source_identity", "row_units", "preservation", "packet_identity", "accepted_function_definition_records", "core_amplitude_Z0.parameter_intervals", "core_amplitude_Z0.exact_parameter_recipes", "core_amplitude_Z0.amplitude", "core_amplitude_Z0.ordinary_Z0_gauge", "core_amplitude_Z0.source_data"]}.

Saved callback metadata: {"callback_rows": 11, "jet_value_keys": 15, "new_math_family_calls": 4, "ordinary_value_keys": 10, "packet_count": 16, "payload_keys": 25, "proof_rows": 1, "rho_values": 4, "tau_values": 4}. Proof summary: {"formal_F0_atom_only": true, "overall_proof_count": 1, "physical_geometry_binding_proof_present": true, "positive_excess_growth_records_present": true, "source_binding_present": true}. Selected rows (approximate display only; exact endpoints remain in the raw result/report):

Rows labelled ut/F0, omega_r/F0 and omega_z/F0 show coefficients of the same strictly positive formal F0 atom. The full quantities equal F0 times these coefficients; a large coefficient is not a large evaluated vorticity. The aspect_excess row is exp(x)-1, rather than the aspect factor exp(x).

| Packet | Location | Lower (approx.) | Upper (approx.) |
|---|---|---:|---:|
| physical_values | quantity=ur, rho=1, tau=1 | `-3.204413665550565142612e-204453045017284784...` | `-3.204413665550565142612e-204453045017284784...` |
| physical_values | quantity=uz, rho=1, tau=1 | `1.747835043960599111830e-22...` | `1.747835043960599111830e-22...` |
| physical_values | quantity=p, rho=1, tau=1 | `-2.084440827025636323978e+204453045017284567...` | `-2.084440827025636323978e+204453045017284567...` |
| physical_values | quantity=ut/F0 (coefficient only), rho=1, tau=1 | `1.233679177816785394171e-204453045017284784...` | `1.233679177816785394171e-204453045017284784...` |
| cylindrical_jets | quantity=omega_r/F0 (coefficient only), rho=1, tau=1 | `5.499116134595799691325e+204453045017284805...` | `5.499116134595799691325e+204453045017284805...` |
| cylindrical_jets | quantity=omega_theta (numeric interval), rho=1, tau=1 | `-1.8713163478976546207428e-248...` | `1.87131634789765462074281e-248...` |
| cylindrical_jets | quantity=omega_z/F0 (coefficient only), rho=1, tau=1 | `1.118271368773813883037e0...` | `1.118271368773813883037e0...` |
| cylindrical_jets | quantity=divergence_diagnostic, rho=1, tau=1 | `-2.213161865127226132340e-220...` | `2.213161865127226132334e-220...` |
| common_source_growth | quantity=aspect_excess=exp(x)-1, rho=1, tau=1 -> 1/16 | `3.280273580956959656052e-408906090034569147...` | `3.280273580956959656052e-408906090034569147...` |

Selected scope flags: {"F0_exp_materialized": false, "N_materialized": false, "common_source_growth_evaluated": true, "full_goal_complete": false, "old_math_reentered": false, "old_provider_called": false, "physical_core_packets_evaluated": true, "physical_field_materialized": false, "physical_geometry_evaluated": true, "physical_jets_evaluated": true, "physical_values_evaluated": true}. Complete exact-real and directed records remain in the saved result/report artifacts rather than this concise document.

Global source-Z jet and parameter admission remain false; no full finite-N/O3 field or global physical admission is claimed. The saved count contract is retained in the report, including 4 new math-family calls, 1 proof, 11 callback rows and 13 journal objects. Attempt status is passed; any nested failure is preserved as null.

No old family/provider/checker was reentered, no F0 exponential or N was materialized, and preserved error/cap metadata remains inherited rather than removed.

Independent review result SHA: 8d825638974a9fcb6430b34ee9a4c7b78376df42ee17a6bcb6494135e518ead8. Journal format: jsonl_compact_one_object_per_line; complete objects: 13; callback objects: 11.

The next frozen task is source_native_O3_off_equator_core_source_value_jet. Its unexecuted plan is retained at C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work/native_O3_off_equator_core_source_value_jet_plan_2026_10_27.json (SHA 1919e1599c98f1e57da8b069b463af9de79c209b5c519f1312c7bef117292fc6); this leaf does not execute that task.
