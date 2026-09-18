import pytest

from openai_ns_reconstruction.constrained_public_visual_target import (
    audit_visual_evidence_receipt,
    public_visual_target,
)


CANDIDATE = "a" * 64


def test_target_records_only_public_qualitative_observables():
    target = public_visual_target()
    assert target.publisher == "OpenAI"
    assert target.published_date == "2026-09-08"
    assert target.source_url == "https://openai.com/index/navier-stokes-solution/"
    assert set(target.observable_ids) == {
        "inward_spiraling",
        "axial_stretching",
        "vortex_swirl",
        "increasing_elongation",
        "central_region_shrinks",
        "speed_increases",
        "relative_angular_rotation_color_encoding",
        "circulating_speed_depends_on_radius",
    }
    assert "numerical_velocity_samples" in target.explicitly_unobserved
    assert "numerical_frame_time" in target.explicitly_unobserved
    assert "camera_pose_or_projection" in target.explicitly_unobserved
    assert "quantitative_visual_acceptance_thresholds" in target.explicitly_unobserved


def test_receipt_reports_coverage_without_visual_or_pde_promotion():
    receipt = audit_visual_evidence_receipt(
        candidate_sha256=CANDIDATE,
        evidence=[
            {
                "observable_id": "inward_spiraling",
                "status": "measured",
                "method": "true-3d streamline review",
                "provenance": "fixed candidate/time/seed receipt",
                "value": {"qualitative_direction": "inward"},
            },
            {
                "observable_id": "axial_stretching",
                "status": "inconclusive",
                "method": "pathline axial-span trend",
                "provenance": "fixed candidate/time window receipt",
            },
            {
                "observable_id": "vortex_swirl",
                "status": "not_measured",
            },
        ],
    )
    assert receipt["measured_observables"] == ("inward_spiraling",)
    assert receipt["inconclusive_observables"] == ("axial_stretching",)
    assert receipt["declared_not_measured"] == ("vortex_swirl",)
    assert receipt["measured_coverage_fraction"] == pytest.approx(1 / 8)
    assert len(receipt["missing_observables"]) == 5
    assert receipt["posthoc_visual_threshold_defined"] is False
    assert receipt["visualization_ready"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False
    with pytest.raises(TypeError):
        receipt["pde_validated"] = True


@pytest.mark.parametrize(
    "bad_key",
    [
        "target_value",
        "target_threshold",
        "acceptance_threshold",
        "hidden_frame_time",
        "camera_pose",
        "camera_registration",
        "streamline_seed_target",
        "recovered_velocity",
        "recovered_parameters",
    ],
)
def test_receipt_rejects_hidden_or_posthoc_target_injection(bad_key):
    with pytest.raises(ValueError):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[
                {
                    "observable_id": "inward_spiraling",
                    "status": "measured",
                    "method": "diagnostic",
                    "provenance": "receipt",
                    bad_key: 0.5,
                }
            ],
        )


def test_receipt_rejects_unknown_duplicate_and_truth_laundering():
    with pytest.raises(ValueError, match="unknown public observable"):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[{"observable_id": "hidden_camera_angle", "status": "not_measured"}],
        )
    with pytest.raises(ValueError, match="duplicate"):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[
                {"observable_id": "vortex_swirl", "status": "not_measured"},
                {"observable_id": "vortex_swirl", "status": "not_measured"},
            ],
        )
    with pytest.raises(ValueError, match="may not promote pde_validated"):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[
                {
                    "observable_id": "vortex_swirl",
                    "status": "measured",
                    "method": "lambda-ci",
                    "provenance": "fixed grid",
                    "pde_validated": True,
                }
            ],
        )


def test_receipt_fails_closed_on_identity_and_provenance_errors():
    with pytest.raises(ValueError, match="candidate_sha256"):
        audit_visual_evidence_receipt(candidate_sha256="deadbeef", evidence=[])
    with pytest.raises(ValueError, match="method"):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[
                {
                    "observable_id": "speed_increases",
                    "status": "measured",
                    "method": "",
                    "provenance": "fixed times",
                }
            ],
        )
    with pytest.raises(ValueError, match="provenance"):
        audit_visual_evidence_receipt(
            candidate_sha256=CANDIDATE,
            evidence=[
                {
                    "observable_id": "speed_increases",
                    "status": "inconclusive",
                    "method": "speed trend",
                    "provenance": "",
                }
            ],
        )
