from __future__ import annotations

import inspect

import numpy as np

from openai_ns_reconstruction.kokuno_strict_inner_transport_axis_completed_radial_stress import (
    INDEPENDENT_SCOPE_AUDIT_HEAD,
    PARENT_AGENT3_HEAD,
    _axis_completed_compact_radial_stress,
    _materialize_axis_completed_from_mean_witness,
    build_mechanics_report,
    materialize_axis_completed_strict_inner_transport_radial_stress,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_mean import (
    _materialize_from_transport_provider,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
    _mechanics_transport_provider,
)


def test_axis_completed_discrete_primitive_includes_zero_to_first_radius_interval():
    radii = np.linspace(0.02, 0.80, 79)
    source = np.ones_like(radii)

    for exponent in (1, 2):
        result = _axis_completed_compact_radial_stress(
            radii,
            source,
            exponent=exponent,
            bump_center=0.55,
            bump_halfwidth=0.20,
        )
        # The compact bump is zero on [0,r_min], so P_e(F)=1 there for this
        # manufactured source.  The discrete axis-completed trapezoid therefore
        # contributes 0.5*r_min^(e+1), whereas the legacy positive-radius
        # primitive starts at r_min and reports zero by construction.
        expected_first_primitive = 0.5 * radii[0] ** (exponent + 1)
        np.testing.assert_allclose(
            result["axis_to_first_radius_weighted_primitive"],
            expected_first_primitive,
            atol=2.0e-18,
            rtol=2.0e-14,
        )
        np.testing.assert_allclose(
            result["stress_inner_positive_edge"],
            -0.5 * radii[0],
            atol=2.0e-16,
            rtol=2.0e-14,
        )
        assert result["legacy_truncated_stress_inner_edge"] == 0.0
        assert result["stress_delta_vs_legacy_truncated_max_abs"] > 0.0
        assert abs(float(result["bump_weighted_integral"]) - 1.0) < 3.0e-14
        assert abs(float(result["moment_complement_weighted_moment"])) < 4.0e-14
        assert abs(float(result["stress_outer_edge"])) < 4.0e-14


def _mechanics_case():
    radii = np.linspace(0.02, 0.80, 79)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.50,
        bump_halfwidth=0.28,
    )
    mean_witness = _materialize_from_transport_provider(
        _mechanics_transport_provider(),
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    result = _materialize_axis_completed_from_mean_witness(
        mean_witness,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    return radii, geometry, mean_witness, result


def test_scoped_typed_mean_is_routed_without_reintroducing_surrogate_source():
    radii, _, mean_witness, result = _mechanics_case()
    np.testing.assert_allclose(
        result.legacy_witness.theta_mean_profile,
        mean_witness.mean_transport_cylindrical[:, 1],
        atol=0.0,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        result.legacy_witness.axial_mean_profile,
        mean_witness.mean_transport_cylindrical[:, 2],
        atol=0.0,
        rtol=0.0,
    )

    theta = result.theta_axis_completed_stress
    axial = result.axial_axis_completed_stress
    assert theta["quadrature_lower_limit"] == 0.0
    assert axial["quadrature_lower_limit"] == 0.0
    assert theta["first_positive_radius"] == radii[0]
    assert axial["first_positive_radius"] == radii[0]
    assert theta["legacy_truncated_stress_inner_edge"] == 0.0
    assert axial["legacy_truncated_stress_inner_edge"] == 0.0
    assert abs(float(theta["stress_inner_positive_edge"])) > 0.0
    assert abs(float(axial["stress_inner_positive_edge"])) > 0.0
    assert float(theta["stress_delta_vs_legacy_truncated_rms"]) > 0.0
    assert float(axial["stress_delta_vs_legacy_truncated_rms"]) > 0.0
    assert abs(float(theta["moment_complement_weighted_moment"])) < 4.0e-14
    assert abs(float(axial["moment_complement_weighted_moment"])) < 4.0e-14
    assert abs(float(theta["stress_outer_edge"])) < 4.0e-14
    assert abs(float(axial["stress_outer_edge"])) < 4.0e-14


def test_public_materializer_has_no_defect_mean_stress_or_retuning_inputs():
    parameters = inspect.signature(
        materialize_axis_completed_strict_inner_transport_radial_stress
    ).parameters
    assert list(parameters) == ["transport_backend", "inner_backend", "geometry"]
    forbidden = {
        "residual",
        "defect",
        "mean",
        "source",
        "stress",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "alpha",
        "damping",
        "angular_order",
        "spatial_step",
        "time_step",
        "viscosity",
        "nu",
        "delta_y",
        "delta_a",
        "normalized_score",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_answers_cr002_lower_limit_audit_without_overclaiming():
    boundary = truth_boundary()
    assert PARENT_AGENT3_HEAD == "9140b0bab1bfc34d496460675e92e63fc6dd721b"
    assert INDEPENDENT_SCOPE_AUDIT_HEAD == "7dddb17c9a853f0aa328c251050901104628ba9c"
    assert boundary["axis_completed_discrete_lower_limit_node_included"] is True
    assert boundary["virtual_axis_weighted_integrand_zero_for_e_ge_1"] is True
    assert boundary["axis_to_first_positive_radius_primitive_recorded"] is True
    assert boundary["legacy_positive_radius_truncation_difference_recorded"] is True
    assert boundary["scoped_discrete_lower_limit_realization_improved"] is True
    assert boundary["formal_axis_based_inverse_verified"] is False
    assert boundary["formal_full_domain_moment_verified"] is False
    assert boundary["axis_regularity_verified"] is False
    assert boundary["bounded_source_to_axis_independently_certified"] is False
    assert boundary["global_radial_domain_materialized"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["global_corrected_leading_join_materialized"] is False
    assert boundary["scoped_transport_stress_authorized_as_correction_target"] is False
    assert boundary["finite_head_debt_authorized_from_scoped_transport"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["same_protocol_comparable_to_st006"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_is_explicitly_not_formal_or_candidate_evidence():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["formal_axis_inverse_certificate"] is False
    assert report["theta_e2_axis_completed"]["stress_delta_vs_legacy_truncated_rms"] > 0.0
    assert report["axial_e1_axis_completed"]["stress_delta_vs_legacy_truncated_rms"] > 0.0
    boundary = report["truth_boundary"]
    assert boundary["formal_axis_based_inverse_verified"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
