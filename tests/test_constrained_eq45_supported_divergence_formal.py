import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_divergence_formal import (
    CLAIM_SCOPE,
    _load_contract,
    audit_supported_candidate_registered_divergence,
    registered_held_out_points,
)


ROOT = Path(__file__).resolve().parents[1]
PARENT_SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"


def _roundtrip_child(tmp_path):
    parent = Eq45VelocityCandidate.load_json(PARENT_SEED)
    child = Eq45SupportedVelocityCandidate(parent=parent)
    path = tmp_path / "eq45_supported_candidate.json"
    child.save_json(path)
    loaded = Eq45SupportedVelocityCandidate.load_json(path)
    assert loaded.sha256 == child.sha256
    return path, loaded


def test_registered_held_out_point_cloud_is_exact_and_reproducible():
    contract = _load_contract(CONSTRAINTS)
    assert contract["seed"] == 914027
    assert contract["point_count"] == 4096
    assert contract["times"] == (0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75)
    assert contract["steps"] == (0.02, 0.01, 0.005)
    assert np.array_equal(contract["lower"], np.array([-2.0, -2.0, -2.0]))
    assert np.array_equal(contract["upper"], np.array([2.0, 2.0, 2.0]))
    assert contract["volume"] == 64.0
    assert contract["divergence_max"] == 1e-5
    assert contract["divergence_L2"] == 1e-5

    first = registered_held_out_points(contract)
    second = registered_held_out_points(contract)
    assert np.array_equal(first, second)
    assert first.shape == (4096, 3)
    np.testing.assert_allclose(
        first[:3],
        np.array(
            [
                [-0.67590863, 0.59701732, 0.57435127],
                [0.74912670, -1.40066948, -0.31510757],
                [0.08397452, 0.34735862, -1.30126462],
            ]
        ),
        rtol=0.0,
        atol=5e-9,
    )
    assert np.all(first >= -2.0)
    assert np.all(first <= 2.0)


def test_supported_child_registered_divergence_metrics(tmp_path, capsys):
    artifact, child = _roundtrip_child(tmp_path)
    report = audit_supported_candidate_registered_divergence(artifact, CONSTRAINTS)

    assert report["claim_scope"] == CLAIM_SCOPE
    assert report["artifact_reloaded"] is True
    assert report["candidate_sha256"] == child.sha256
    assert report["parent_sha256"] == child.parent_sha256
    assert report["velocity_access"] == "at_points_only"
    assert report["sampling"] == {
        "distribution": "uniform_cartesian_box",
        "seed": 914027,
        "point_count": 4096,
        "same_spatial_cloud_reused_at_each_time": True,
        "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "box_volume": 64.0,
        "times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
    }
    assert report["registered_thresholds"] == {
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    assert [level["step"] for level in report["levels"]] == [0.02, 0.01, 0.005]

    for level in report["levels"]:
        assert len(level["per_time"]) == 6
        assert np.isfinite(level["worst_sampled_max_abs"])
        assert np.isfinite(level["worst_volume_weighted_L2"])
        for item in level["per_time"]:
            assert np.isfinite(item["sampled_max_abs"])
            assert np.isfinite(item["sampled_rms"])
            assert np.isfinite(item["volume_weighted_L2"])
            # The registered evaluation box has volume 64, so this catches an
            # accidental substitution of sampled RMS for spatial L2.
            np.testing.assert_allclose(
                item["volume_weighted_L2"],
                8.0 * item["sampled_rms"],
                rtol=2e-13,
                atol=1e-15,
            )

    assert len(report["observed_volume_L2_orders_by_time"]) == 6
    for time_report in report["observed_volume_L2_orders_by_time"]:
        assert len(time_report["volume_L2_orders"]) == 2
        for order in time_report["volume_L2_orders"]:
            assert order is None or np.isfinite(order)

    finest = report["levels"][-1]
    assert report["finest_level"] == {
        "step": 0.005,
        "worst_sampled_max_abs": finest["worst_sampled_max_abs"],
        "worst_volume_weighted_L2": finest["worst_volume_weighted_L2"],
        "registered_sample_max_pass": finest["worst_sampled_max_abs"] <= 1e-5,
        "registered_volume_L2_pass": finest["worst_volume_weighted_L2"] <= 1e-5,
    }
    assert report["registered_finite_divergence_gate_assessed"] is True
    assert report["registered_finite_divergence_gate_passed"] is (
        report["finest_level"]["registered_sample_max_pass"]
        and report["finest_level"]["registered_volume_L2_pass"]
    )

    # Even a numerical divergence pass would not establish the momentum PDE or
    # a continuum theorem, and it must not promote visualization/source claims.
    assert report["continuum_supremum_proved"] is False
    assert report["full_momentum_residual_assessed"] is False
    assert report["pde_validated"] is False
    assert report["physical_support_validated"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False

    # Preserve deterministic production numerics in exact-head Actions logs.
    with capsys.disabled():
        print("CR006_SUPPORTED_FORMAL_DIVERGENCE=" + json.dumps(report, sort_keys=True))
