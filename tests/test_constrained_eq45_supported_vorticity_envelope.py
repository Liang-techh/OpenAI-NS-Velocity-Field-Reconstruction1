import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_vorticity_envelope import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIMES,
    audit_supported_vorticity_envelope,
    structured_vorticity_magnitude,
)


def test_structured_vorticity_calibrates_rigid_rotation():
    class RigidRotation:
        @staticmethod
        def at_points(points, time):
            points = np.asarray(points, dtype=float)
            out = np.zeros_like(points)
            out[..., 0] = -points[..., 1]
            out[..., 1] = points[..., 0]
            return out

    axis = np.linspace(-1.0, 1.0, 17)
    magnitude, radius, abs_z = structured_vorticity_magnitude(RigidRotation(), axis, 0.5)
    np.testing.assert_allclose(magnitude, 2.0, rtol=0.0, atol=2e-14)
    assert magnitude.shape == radius.shape == abs_z.shape == (17, 17, 17)


def test_supported_vorticity_envelope_detects_resolved_early_radial_collar():
    report = audit_supported_vorticity_envelope()

    assert report["schema"] == "eq45_supported_vorticity_envelope_audit_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["times"] == list(DEFAULT_TIMES)
    assert len(report["rows"]) == len(DEFAULT_RESOLUTIONS) * len(DEFAULT_TIMES)

    summaries = {row["time"]: row for row in report["time_summaries"]}
    early = summaries[0.25]
    middle = summaries[0.50]
    late = summaries[0.75]

    # Whole-domain view adds information not covered by the pointwise collar audit:
    # at the early slice the support connection creates a resolved radial
    # superlevel branch even though the axial 99%-extent is unchanged.
    assert early["child_parent_radial_q99_ratio"] > 1.35
    assert abs(early["child_parent_axial_q99_ratio"] - 1.0) < 0.02
    assert early["child_parent_aspect_ratio"] < 0.80
    assert early["child_radial_q99_mid_to_fine_relative_change"] < 0.02
    assert early["child_minus_parent_collar_vorticity2_fraction"] > 0.03

    # Once the similarity core has contracted into the taper plateau, the same
    # 25%-of-core envelope is essentially unchanged by the support transform.
    for summary in (middle, late):
        assert abs(summary["child_parent_radial_q99_ratio"] - 1.0) < 0.02
        assert abs(summary["child_parent_axial_q99_ratio"] - 1.0) < 0.02
        assert abs(summary["child_parent_aspect_ratio"] - 1.0) < 0.02
        assert abs(summary["child_minus_parent_collar_vorticity2_fraction"]) < 0.002

    truth = report["truth_boundary"]
    assert truth["physical_support_connection_implemented"] is True
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
