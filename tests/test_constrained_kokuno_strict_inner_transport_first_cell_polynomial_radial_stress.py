from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_strict_inner_transport_first_cell_polynomial_radial_stress import (
    INDEPENDENT_SCOPE_AUDIT_HEAD,
    PARENT_AGENT3_HEAD,
    _cumulative_with_polynomial_first_cell,
    _first_cell_polynomial_integral,
    _formal_polynomial_first_cell,
    _polynomial_first_cell_compact_radial_stress,
    build_mechanics_report,
    materialize_polynomial_first_cell_strict_inner_transport_radial_stress,
    truth_boundary,
)


def test_quadratic_unweighted_first_cell_is_integrated_exactly_for_e1_e2():
    radii = np.array([0.2, 0.31, 0.47, 0.62, 0.79], dtype=float)
    coefficients = (1.2, -0.7, 0.4)
    values = coefficients[0] + coefficients[1] * radii + coefficients[2] * radii**2

    for exponent in (1, 2):
        actual = _first_cell_polynomial_integral(
            values, radii, exponent=exponent, degree=2
        )
        expected = _formal_polynomial_first_cell(
            coefficients, float(radii[0]), exponent
        )
        np.testing.assert_allclose(actual, expected, rtol=3.0e-14, atol=3.0e-17)

        cumulative = _cumulative_with_polynomial_first_cell(
            values, radii, exponent=exponent, degree=2
        )
        np.testing.assert_allclose(
            cumulative[0], expected, rtol=3.0e-14, atol=3.0e-17
        )
        assert np.all(np.isfinite(cumulative))


def test_constant_source_removes_cr002_e2_first_cell_overestimate():
    radii = np.linspace(0.2, 1.0, 9)
    source = np.ones_like(radii)

    e1 = _polynomial_first_cell_compact_radial_stress(
        radii,
        source,
        exponent=1,
        bump_center=0.7,
        bump_halfwidth=0.2,
    )
    e2 = _polynomial_first_cell_compact_radial_stress(
        radii,
        source,
        exponent=2,
        bump_center=0.7,
        bump_halfwidth=0.2,
    )

    r0 = float(radii[0])
    np.testing.assert_allclose(
        e1["axis_to_first_radius_weighted_primitive"],
        r0**2 / 2.0,
        rtol=3.0e-14,
        atol=3.0e-17,
    )
    np.testing.assert_allclose(
        e2["axis_to_first_radius_weighted_primitive"],
        r0**3 / 3.0,
        rtol=3.0e-14,
        atol=3.0e-17,
    )
    np.testing.assert_allclose(
        e2["parent_868_first_cell_primitive"],
        0.5 * r0**3,
        rtol=3.0e-14,
        atol=3.0e-17,
    )
    np.testing.assert_allclose(
        e2["stress_inner_positive_edge"],
        -r0 / 3.0,
        rtol=3.0e-14,
        atol=3.0e-16,
    )
    np.testing.assert_allclose(
        e2["parent_868_stress_inner_positive_edge"],
        -0.5 * r0,
        rtol=3.0e-14,
        atol=3.0e-16,
    )
    assert float(e2["stress_delta_vs_parent_868_rms"]) > 0.0

    for result in (e1, e2):
        assert abs(float(result["bump_weighted_integral"]) - 1.0) < 3.0e-13
        assert abs(float(result["moment_complement_weighted_moment"])) < 3.0e-12
        assert abs(float(result["stress_outer_edge"])) < 3.0e-12
        assert result["first_cell_polynomial_degree"] == 2
        assert np.isfinite(
            float(result["source_first_cell_linear_quadratic_relative_difference"])
        )


def test_polynomial_first_cell_rejects_bad_grids_and_shapes():
    radii = np.array([0.2, 0.3, 0.4], dtype=float)
    values = np.ones(3)
    with pytest.raises(ValueError, match="matching 1D arrays"):
        _cumulative_with_polynomial_first_cell(
            values[:, None], radii, exponent=2, degree=2
        )
    with pytest.raises(ValueError, match="positive and strictly increasing"):
        _cumulative_with_polynomial_first_cell(
            values, np.array([0.2, 0.2, 0.4]), exponent=2, degree=2
        )
    with pytest.raises(ValueError, match="only e=1 and e=2"):
        _cumulative_with_polynomial_first_cell(
            values, radii, exponent=3, degree=2
        )


def test_public_materializer_keeps_surrogate_and_retuning_inputs_closed():
    parameters = inspect.signature(
        materialize_polynomial_first_cell_strict_inner_transport_radial_stress
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


def test_truth_boundary_answers_cr002_first_cell_scope_without_overclaiming():
    boundary = truth_boundary()
    assert PARENT_AGENT3_HEAD == "c8268100a36942de9fa966d76f22ad6148ceaa02"
    assert INDEPENDENT_SCOPE_AUDIT_HEAD == "594af2273c2c27a7a64a9a4340da43472343b83a"
    assert boundary["constant_source_cr002_e2_mechanics_defect_removed"] is True
    assert boundary["quadratic_unweighted_first_cell_mechanics_exact"] is True
    assert boundary["real_source_first_cell_linear_quadratic_disagreement_recorded"] is True
    assert boundary["real_source_first_cell_convergence_verified"] is False
    assert boundary["formal_first_cell_accuracy_verified"] is False
    assert boundary["formal_axis_based_inverse_verified"] is False
    assert boundary["formal_full_domain_moment_verified"] is False
    assert boundary["axis_regularity_verified"] is False
    assert boundary["bounded_source_to_axis_independently_certified"] is False
    assert boundary["global_radial_domain_materialized"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["scoped_transport_stress_authorized_as_correction_target"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["same_protocol_comparable_to_st006"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_records_exactness_but_not_candidate_or_formal_evidence():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["formal_axis_inverse_certificate"] is False
    for exponent in (1, 2):
        check = report["polynomial_first_cell_checks"][f"e{exponent}"]
        assert check["quadratic_rule_absolute_error"] < 5.0e-17
        assert check["constant_source_absolute_error"] < 5.0e-17
    typed = report["typed_routing_receipt"]
    assert abs(
        float(
            typed["theta_e2_polynomial_first_cell"][
                "moment_complement_weighted_moment"
            ]
        )
    ) < 1.0e-11
    assert abs(
        float(
            typed["axial_e1_polynomial_first_cell"][
                "moment_complement_weighted_moment"
            ]
        )
    ) < 1.0e-11
    assert report["truth_boundary"]["pde_validated"] is False
