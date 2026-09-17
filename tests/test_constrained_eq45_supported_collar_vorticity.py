import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    DEFAULT_SPATIAL_STEPS,
    PARENT_SHA256,
    SUPPORTED_SHA256,
    audit_supported_collar_vorticity,
    centered_point_vorticity,
    governed_supported_seed,
)


def test_centered_point_vorticity_calibrates_rigid_rotation():
    def rigid_rotation(points, time):
        points = np.asarray(points, dtype=float)
        out = np.zeros_like(points)
        out[:, 0] = -points[:, 1]
        out[:, 1] = points[:, 0]
        return out

    points = np.asarray(
        [[0.2, -0.3, 0.4], [-0.8, 0.1, -0.2], [1.1, 0.7, 0.0]],
        dtype=float,
    )
    omega = centered_point_vorticity(rigid_rotation, points, 0.5, 0.01)
    np.testing.assert_allclose(omega, np.asarray([[0.0, 0.0, 2.0]] * 3), atol=2e-14)


def test_supported_child_collar_vorticity_is_resolved_and_truth_bounded():
    report = audit_supported_collar_vorticity(governed_supported_seed())

    assert report["schema"] == "eq45_supported_collar_vorticity_audit_v1"
    assert report["parent_sha256"] == PARENT_SHA256
    assert report["supported_child_sha256"] == SUPPORTED_SHA256
    assert report["steps"] == list(DEFAULT_SPATIAL_STEPS)
    assert len(report["rows"]) == 3 * 3 * 4
    assert len(report["time_summaries"]) == 3

    early_taper_signal = 0.0
    for summary in report["time_summaries"]:
        assert summary["plateau_parent_child_relative_mismatch"] < 1e-10
        for region in ("radial_collar", "axial_collar", "corner_collar"):
            collar_ratio = summary[f"{region}_child_to_plateau_rms_ratio"]
            taper_ratio = summary[f"{region}_taper_delta_to_plateau_rms_ratio"]
            refinement = summary[f"{region}_mid_to_fine_relative_delta"]
            assert np.isfinite(collar_ratio)
            assert np.isfinite(taper_ratio)
            assert np.isfinite(refinement)
            # The support connection must not create a visualization-dominating
            # vorticity shell tens of times stronger than the resolved core.
            assert collar_ratio < 25.0
            assert taper_ratio < 25.0
            # Ignore relative refinement only when the collar is effectively
            # zero compared with the core; otherwise require a resolved signal.
            if collar_ratio > 1e-8:
                assert refinement < 0.10
        if summary["time"] == 0.25:
            early_taper_signal = max(
                summary["radial_collar_taper_delta_to_plateau_rms_ratio"],
                summary["axial_collar_taper_delta_to_plateau_rms_ratio"],
                summary["corner_collar_taper_delta_to_plateau_rms_ratio"],
            )

    # The child is genuinely different in the early support collar; this is not
    # merely a replay of the untapered parent.
    assert early_taper_signal > 1e-5

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


def test_supported_collar_audit_rejects_malformed_step_ladders():
    child = governed_supported_seed()
    with pytest.raises(ValueError, match="three"):
        audit_supported_collar_vorticity(child, steps=(0.02, 0.01))
    with pytest.raises(ValueError, match="decreasing"):
        audit_supported_collar_vorticity(child, steps=(0.01, 0.02, 0.005))
