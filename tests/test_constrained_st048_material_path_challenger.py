import numpy as np
import pytest

from openai_ns_reconstruction import constrained_st048_material_path_challenger as m


def _linear_swirl(points, _time):
    p = np.asarray(points, dtype=float)
    x, y, z = p.T
    alpha, omega, beta = 0.2, 1.1, 0.25
    return np.column_stack((-alpha * x - omega * y, omega * x - alpha * y, beta * z))


def _synthetic_receipt_inputs():
    summary = m.measure_material_paths(_linear_swirl)
    measurements = {"ST048-S": summary, "ST048-B": dict(summary)}
    metadata = {
        candidate: {
            "candidate_id": candidate,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "source_pr": 390,
        }
        for candidate in measurements
    }
    return measurements, metadata


def test_frozen_contract_and_source_identity():
    seeds, metadata = m._fixed_seed_table()
    assert seeds.shape == (48, 3)
    assert len(metadata) == 48
    assert m.SEED_RADII == (0.6, 0.9, 1.2)
    assert m.SEED_Z == (-0.3, 0.3)
    assert m.REGISTERED_TIME_INTERVAL == (0.25, 0.75)
    assert m.ST048_HEAD_SHA == "97695a86f85ce68fb4ae70c41fc81c904d655183"
    assert set(m.ST048_RAW_SHA256) == {"ST048-S", "ST048-B"}


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


def test_zero_and_malformed_velocity_fail_closed():
    with pytest.raises(ValueError, match="numerically zero"):
        m.measure_material_paths(lambda points, _time: np.zeros_like(points))
    with pytest.raises(ValueError, match=r"shape \(n,3\)"):
        m.measure_material_paths(lambda points, _time: np.zeros((len(points), 2)))


def test_receipt_binds_two_children_and_frozen_references_without_promotion():
    measurements, metadata = _synthetic_receipt_inputs()
    receipt = m.build_receipt(measurements, metadata)
    assert receipt["references"]["st047e"]["source_report_sha256"] == m.ST047_REFERENCE_REPORT_SHA256
    assert receipt["references"]["st006"]["source_report_sha256"] == m.ST006_REFERENCE_REPORT_SHA256
    assert set(receipt["descriptive_comparison"]) == {"ST048-S", "ST048-B"}
    assert set(receipt["descriptive_comparison"]["ST048-S"]) == {"vs_st047e", "vs_st006"}
    assert "pass" not in receipt["within_st048"]["status"]
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


def test_identity_and_claim_laundering_fail_closed():
    measurements, metadata = _synthetic_receipt_inputs()
    bad = {k: dict(v) for k, v in metadata.items()}
    bad["ST048-S"]["pde_validated"] = True
    with pytest.raises(ValueError, match="PDE truth state"):
        m.build_receipt(measurements, bad)

    bad = {k: dict(v) for k, v in metadata.items()}
    bad["ST048-B"]["source_correspondence_verified"] = True
    with pytest.raises(ValueError, match="source-correspondence truth state"):
        m.build_receipt(measurements, bad)

    bad_measurements = {"ST048-S": measurements["ST048-S"]}
    with pytest.raises(ValueError, match="exactly ST048-S and ST048-B"):
        m.build_receipt(bad_measurements, {"ST048-S": metadata["ST048-S"]})
