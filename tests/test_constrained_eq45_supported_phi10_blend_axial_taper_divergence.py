import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_divergence import (
    AXIAL_PLATEAU_VALUES,
    REGIONS,
    REGION_SIZE,
    SPATIAL_STEPS,
    TIMES,
    audit_axial_taper_divergence,
    fresh_stratified_points,
)


def test_fresh_axial_taper_divergence_probe_contract_is_deterministic_and_buffered():
    first_points, first_labels = fresh_stratified_points()
    second_points, second_labels = fresh_stratified_points()
    np.testing.assert_array_equal(first_points, second_points)
    np.testing.assert_array_equal(first_labels, second_labels)
    assert first_points.shape == (REGION_SIZE * len(REGIONS), 3)
    assert set(first_labels.tolist()) == set(REGIONS)
    for region in REGIONS:
        assert int(np.sum(first_labels == region)) == REGION_SIZE

    radii = np.hypot(first_points[:, 0], first_points[:, 1])
    assert float(np.max(radii)) < 1.90
    assert float(np.max(np.abs(first_points[:, 2]))) < 1.93


def test_axial_taper_family_public_divergence_converges_and_detects_mutation():
    report = audit_axial_taper_divergence()
    assert report["velocity_access"] == "serialized_reloaded_public_at_points_only"
    assert report["operator"] == "independent_centered_second_order_cartesian_space"
    assert report["axial_taper_value_selected"] is False
    assert report["registered_full_domain_divergence_gate_assessed"] is False
    assert report["full_momentum_or_vorticity_residual_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["pde_validated"] is False

    members = report["members"]
    assert [member["axial_plateau_q"] for member in members] == list(AXIAL_PLATEAU_VALUES)
    for member in members:
        assert member["artifact_reloaded"] is True
        assert len(member["candidate_sha256"]) == 64
        assert len(member["time_reports"]) == len(TIMES)
        for time_report in member["time_reports"]:
            levels = time_report["levels"]
            assert [level["step"] for level in levels] == list(SPATIAL_STEPS)
            assert levels[-1]["sampled_rms"] < levels[0]["sampled_rms"]
            assert levels[-1]["sampled_max_abs"] < levels[0]["sampled_max_abs"]
            assert levels[-1]["gradient_normalized_rms"] < levels[0]["gradient_normalized_rms"]
            assert all(order is not None and order > 0.5 for order in time_report["observed_rms_orders"])
            for level in levels:
                assert set(level["regions"]) == set(REGIONS)
                assert np.isfinite(level["sampled_rms"])
                assert np.isfinite(level["sampled_max_abs"])
                assert np.isfinite(level["gradient_normalized_rms"])

    baseline = members[0]
    for time_report in baseline["time_reports"]:
        response = time_report["velocity_response_vs_q064"]
        assert response["rms_vector_change"] == 0.0
        assert response["relative_rms_vector_change"] == 0.0
        assert response["sampled_max_vector_change"] == 0.0

    extended = members[-1]
    assert any(
        item["velocity_response_vs_q064"]["rms_vector_change"] > 1.0e-8
        for item in extended["time_reports"]
    )

    mutation = report["mutation_calibration"]
    assert mutation["mutation"] == "u <- u + 0.03*x"
    assert mutation["expected_added_divergence"] == 0.03
    assert mutation["mutant_finest_rms"] > 0.02
    assert mutation["mutant_finest_rms"] > 5.0 * mutation["baseline_finest_rms"]
    assert mutation["mutant_finest_max_abs"] > 0.02
