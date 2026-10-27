# Current native O3 exterior differential spacetime pullback — passed conditionally (2026-10-27)

This is a source-only metadata record for native_O3_exterior_differential_spacetime_pullback. The saved result and independent review accept only the fresh conditional geometry/potential differentials and exterior-form closure. It is not a numerical field construction, a numeric inverse oracle, or a Navier–Stokes proof.

Readable formula packet for the declared leading n = 0 regional O3 route:

```text
d = 1 - Z^2                         L = 1 - delta*Z^2
lambda_phys = sqrt(tau/d)             r = lambda_phys*sqrt(2*R)
z = lambda_phys^(1-delta)*Z            t = 1 - tau
psi = lambda_phys^(1-delta)*Pstar*R*m
A_one_form = psi*dtheta
m_y = V - m,   d_profile(m) = (m_y/R)dR + m_Z*dZ

Fresh dual/one-form algebra in (R,Z):
    r_R, r_Z, z_R, z_Z, psi_R, psi_Z
    J2 = r_R*z_Z - r_Z*z_R
    u_r = -(-psi_R*r_Z + psi_Z*r_R)/(r*J2)
    u_z =  (psi_R*z_Z - psi_Z*r_R)/(r*J2)
    u_theta = lambda_phys^(-1-delta)*Pstar*E

u = u_r*cos(theta) - u_theta*sin(theta)
v = u_r*sin(theta) + u_theta*cos(theta)
w = u_z
p_phys = lambda_phys^(-2-2*delta)*Pstar^2*(P0 + p_hist)

q_inverse > 0, q_inverse - z_phys^2*q_inverse^delta = 1 - t_phys
R = (X^2 + Y^2)/(2*q_inverse)
Z = z_phys / q_inverse^((1-delta)/2)
x_native = log(R/Rm) - 17 - exp(40)

physical volume primary: r*J2 dR dZ dtheta (positive orientation (R,theta,Z))
positive Jacobian magnitude in order (R,Z,theta):
    lambda_phys^(3-delta)*L/d
native O3 x measure carries the additional profile-R factor.
The paper closed-form expression is retained as a string-level reference only;
closed_form_symbolic_comparison_performed = False and this is not a numeric L2 norm.

The meridional identity is i_u(r dr wedge dtheta wedge dz) =
dpsi wedge dtheta. It is conditional on C2 regular phase cells and fixed
parameters. Positivity premises are Pstar>0, Rm>0, R>0, tau>0, |Z|<1,
and 0<delta<1/200. It is not a numeric inverse, a global physical field, an NS
validation, an axis/phase-seam extension, an unconditional Z4 admission, or
higher-order time recursion.
```

The two-component dual computes the six R/Z derivatives shown above from the new geometry and potential ansatz. The exterior identity is conditional on C2 mixed-partial commutation and antisymmetry in regular interior phase cells; fixed Pstar, delta, and Rm are premises. Physical theta is used in the Cartesian rotation; it is not the source fast phase.

The background and actual modified branches use this same source frame; incoming histories remain supplied graph data and m_y = V - m.

In psi, the profile coordinate R multiplies the mass-mean m. Rm is a separate inlet-radius port and appears only in the native x inverse; R*m is never merged into Rm.

Public source function definitions are limited to leading n = 0, regular O3 phase cells, and fixed parameters. There is no axis extension, unconditional Z4/source-value admission, n >= 2 recursion, global physical closure, NS validation, or cap deletion. Native identifiers Pstar = 33572, delta = 25220, R = 241787, Z = 20, and Rm = 85743 are graph references, not values; the N graph reference is node 18.

Actual saved counts: {"Cartesian_composition": 2, "N_materialization_calls": 0, "ancestor_calls": 0, "api_refusal": 4, "exterior_closure": 2, "failed_family_calls": 0, "flux_form": 2, "geometry_dual": 1, "new_differential_family_calls": 1, "new_mathematical_proofs": 1, "numeric_evaluation_calls": 0, "old_math_calls": 0, "potential": 2, "provider_calls": 0, "source_family_calls": 1, "source_join": 1, "source_value_reducer_calls": 0, "temporal_recursion_calls": 0, "velocity_pressure": 2, "volume_form": 1}.
Actual saved flags: {"Cartesian_velocity_pressure_functions_defined": true, "actual_incoming_x_independence_proved": false, "axis_extension_proved": false, "cap_deletion": false, "conditional_exterior_divergence_identity_proved": true, "conditional_spacetime_pullback_defined": true, "full_goal_complete": false, "global_phase_seam_regularity_proved": false, "implicit_Cartesian_function_composition_defined": true, "native_field_derivative_identity_proved": false, "numeric_closure": false, "numerical_oracle_available": false, "physical_layer_completed": false, "physical_volume_expression_defined": true, "same_source_exterior_potential_defined": true, "temporal_recursion": false, "unconditional_Z4_source_admission": false, "unconditional_parameter_admission": false}.
Review result SHA: 1706f3092b0f5e7f2996c22852c23831038d1cab64c5daf3c97f0e848253de20.
Attempt status: passed. Nested failure, if any, is preserved exactly: null.

Equation references: [Lei-Ren Part I](https://arxiv.org/abs/2609.35406v2) and [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

Frozen prospective next task: source_native_O3_physical_spatial_jet_and_volume_norm_route. Next task after this successful attempt: source_native_O3_physical_spatial_jet_and_volume_norm_route. A failed attempt must not mutate the frozen prospective manifest task; any retry is a different independent method.
