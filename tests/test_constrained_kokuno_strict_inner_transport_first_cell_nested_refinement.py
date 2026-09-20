from __future__ import annotations

import inspect

import numpy as np

from openai_ns_reconstruction.kokuno_strict_inner_transport_first_cell_nested_refinement import (
    ABSOLUTE_STABILITY_FLOOR,
    INDEPENDENT_OBSERVABILITY_AUDIT_HEAD,
    PARENT_AGENT3_HEAD,
    REFINEMENT_RATIO_GATE,
    _nested_geometries,
    build_mechanics_report,
    materialize_nested_first_cell_refinement,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)


def _geometry() -> StrictInnerTransportRadialGeometry:
    radii = np.linspace(0.20, 1.00, 17)
    return StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.65,
        bump_halfwidth=0.18,
    )


def test_nested_geometry_only_adds_frozen_near_axis_positive_radii():
    base = _geometry()
    coarse, medium, fine = _nested_geometries(base)
    r = np.asarray(base.radii, dtype=float)
    np.testing.assert_array_equal(np.asarray(coarse.radii), r)
    np.testing.assert_array_equal(np.asarray(medium.radii)[1:], r)
    np.testing.assert_array_equal(np.asarray(fine.radii)[2:], r)
    assert medium.radii[0] == 0.5 * r[0]
    assert fine.radii[0] == 0.25 * r[0]
    assert fine.radii[1] == 0.5 * r[0]
    for level in (coarse, medium, fine):
        assert level.time == base.time
        assert level.axial_z == base.axial_z
        assert level.bump_center == base.bump_center
        assert level.bump_halfwidth == base.bump_halfwidth


def test_mechanics_nested_refinement_improves_on_common_radii_without_formal_promotion():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["formal_first_cell_certificate"] is False
    witness = report["witness"]
    assert witness["common_mean_replay_max_abs"] < 2.0e-13
    assert witness["engineering_stability_passed"] is True
    for channel in ("theta_e2_refinement", "axial_e1_refinement"):
        metrics = witness[channel]
        assert metrics["engineering_stability_gate_passed"] is True
        assert metrics["medium_fine_common_radii_rms"] < metrics["coarse_medium_common_radii_rms"]
        if not metrics["absolute_stability_floor_reached"]:
            assert metrics["raw_refinement_ratio"] >= REFINEMENT_RATIO_GATE
        assert metrics["finite_positive_radius_refinement_is_formal_error_bound"] is False
    assert report["truth_boundary"]["real_source_first_cell_convergence_verified"] is False


def test_public_api_keeps_surrogate_refinement_retuning_inputs_closed():
    parameters = inspect.signature(materialize_nested_first_cell_refinement).parameters
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
        "refinement_ratio_gate",
        "stability_floor",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_answers_cr002_observability_without_overclaiming():
    boundary = truth_boundary()
    assert PARENT_AGENT3_HEAD == "4b302412fe5733690394002a947eeedaa61a1835"
    assert INDEPENDENT_OBSERVABILITY_AUDIT_HEAD == "1763282ab49c72a8998628805b53ed69f5e3890b"
    assert boundary["nested_first_radius_factors"] == [1.0, 0.5, 0.25]
    assert boundary["refinement_ratio_gate"] == REFINEMENT_RATIO_GATE == 1.5
    assert boundary["absolute_stability_floor"] == ABSOLUTE_STABILITY_FLOOR == 1.0e-10
    assert boundary["engineering_gate_frozen_before_exact_head_output"] is True
    assert boundary["real_source_positive_radius_nested_refinement_assessed"] is True
    assert boundary["finite_positive_radius_refinement_stability_is_formal_error_bound"] is False
    assert boundary["finite_nested_refinement_closes_hidden_cell_observability"] is False
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
