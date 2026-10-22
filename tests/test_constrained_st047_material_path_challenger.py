import numpy as np
import pytest

from openai_ns_reconstruction import constrained_st047_material_path_challenger as m


def _linear_swirl(points, _time):
    p = np.asarray(points, dtype=float)
    x, y, z = p.T
    alpha, omega, beta = 0.2, 1.1, 0.25
    return np.column_stack((-alpha * x - omega * y, omega * x - alpha * y, beta * z))


def test_frozen_seed_contract_is_48_paths_and_24_pairs():
    seeds, metadata = m._fixed_seed_table()
    assert seeds.shape == (48, 3)
    assert len(metadata) == 48
    assert m.SEED_RADII == (0.6, 0.9, 1.2)
    assert m.SEED_Z == (-0.3, 0.3)
    assert m.REGISTERED_TIME_INTERVAL == (0.25, 0.75)
    assert m.ST047_HEAD_SHA == "2c51cd20e036ab29954cf77a70814bc918d5c6a5"


def test_analytic_contracting_swirl_and_axial_stretch():
    report = m.measure_material_paths(_linear_swirl)
    assert report["path_count"] == 48
    assert report["paired_material_line_count"] == 24
    assert report["inward_path_count"] == 48
    assert report["outward_path_count"] == 0
    assert report["mean_radius_change"] < 0.0
    assert report["mean_absolute_turns"] > 0.08
    assert report["pair_axial_separation_growth_count"] == 24
    assert report["pair_axial_separation_shrink_count"] == 0
    assert report["mean_pair_axial_separation_ratio"] > 1.0


def test_zero_velocity_fails_closed():
    with pytest.raises(ValueError, match="numerically zero"):
        m.measure_material_paths(lambda points, _time: np.zeros_like(points))


def test_receipt_binds_both_references_without_promoting_truth():
    summary = m.measure_material_paths(_linear_swirl)
    metadata = {
        "pde_validated": False,
        "source_correspondence_verified": False,
        "source_pr": 379,
    }
    receipt = m.build_receipt(summary, metadata)
    assert receipt["references"]["st046a"]["source_report_sha256"] == m.ST046_REFERENCE_REPORT_SHA256
    assert receipt["references"]["st006"]["source_report_sha256"] == m.ST006_REFERENCE_REPORT_SHA256
    assert set(receipt["descriptive_comparison"]) == {"vs_st046a", "vs_st006"}
    truth = receipt["truth_boundary"]
    assert truth["comparison_is_descriptive_not_acceptance"] is True
    assert truth["visual_acceptance_threshold_defined"] is False
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
    assert "pass" not in receipt["descriptive_comparison"]
    assert "threshold" not in receipt["descriptive_comparison"]


def test_claim_laundering_and_malformed_shapes_fail_closed():
    summary = m.measure_material_paths(_linear_swirl)
    with pytest.raises(ValueError, match="PDE truth state"):
        m.build_receipt(summary, {"pde_validated": True, "source_correspondence_verified": False})
    with pytest.raises(ValueError, match="source-correspondence truth state"):
        m.build_receipt(summary, {"pde_validated": False, "source_correspondence_verified": True})
    with pytest.raises(ValueError, match=r"shape \(n,3\)"):
        m.measure_material_paths(lambda points, _time: np.zeros((len(points), 2)))
