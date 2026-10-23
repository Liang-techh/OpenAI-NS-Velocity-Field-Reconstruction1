from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_vorticity import (
    audit_supported_eq45,
    vorticity_equation_terms,
)


class _RigidRotation:
    def at_points(self, points, time):
        pts = np.asarray(points, dtype=float)
        return np.column_stack((-pts[:, 1], pts[:, 0], np.zeros(len(pts))))


class _SinShear:
    def at_points(self, points, time):
        pts = np.asarray(points, dtype=float)
        return np.column_stack((np.sin(pts[:, 1]), np.zeros(len(pts)), np.zeros(len(pts))))


def test_independent_operator_calibrates_rigid_rotation_and_viscous_shear():
    points = np.asarray(
        [
            [-0.21, -0.37, 0.11],
            [0.17, -0.05, -0.09],
            [0.08, 0.31, 0.14],
            [0.26, 0.48, -0.18],
        ],
        dtype=float,
    )
    times = np.asarray([0.36, 0.44, 0.56, 0.64], dtype=float)

    rigid = vorticity_equation_terms(
        _RigidRotation(),
        points,
        times,
        spatial_step=0.01,
        time_step=0.0025,
        nu=0.01,
    )
    assert np.max(np.abs(rigid["residual"])) < 1.0e-8

    shear = vorticity_equation_terms(
        _SinShear(),
        points,
        times,
        spatial_step=0.01,
        time_step=0.0025,
        nu=0.01,
    )
    expected = np.column_stack(
        (
            np.zeros(len(points)),
            np.zeros(len(points)),
            -0.01 * np.cos(points[:, 1]),
        )
    )
    np.testing.assert_allclose(shear["residual"], expected, rtol=0.0, atol=4.0e-6)


def test_supported_child_revalidation_uses_public_velocity_and_preserves_truth_boundary():
    report = audit_supported_eq45()

    assert report["schema"] == "eq45_supported_vorticity_revalidation_v1"
    assert report["serialized_reloaded_child"] is True
    assert report["velocity_access"] == "public_at_points_only"
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["probe_count"] == 16
    assert report["regions"] == [
        "plateau",
        "radial_collar",
        "axial_collar",
        "corner_collar",
    ]
    assert report["parent_sha256"] != report["supported_child_sha256"]
    assert report["forcing_family"] == "restricted_two_parameter_family"
    assert report["supported_force_coefficients_status"] == (
        "pending_refit_and_independent_revalidation"
    )
    assert report["untapered_fitted_force_transferred"] is False
    assert report["formal_pde_gate_assessed"] is False

    rows = report["rows"]
    assert len(rows) == 3
    for row in rows:
        for field_name in ("parent", "supported_child"):
            metrics = row[field_name]
            assert metrics["max"] > 0.0
            assert metrics["rms"] > 0.0
            assert metrics["term_normalized_rms"] > 0.0
            assert np.isfinite(list(metrics.values())).all()
        assert set(row["by_region"]) == set(report["regions"])

    # The support transform is exactly identity on the interior plateau, including
    # all nested finite-difference offsets used by this probe set.
    assert report["finest_plateau_parent_child_rms_relative_difference"] < 1.0e-10

    # At least one collar must feel the transform; otherwise this would not be a
    # child-specific revalidation.
    finest = rows[-1]
    collar_changes = []
    for region in ("radial_collar", "axial_collar", "corner_collar"):
        parent_rms = finest["by_region"][region]["parent"]["rms"]
        child_rms = finest["by_region"][region]["supported_child"]["rms"]
        collar_changes.append(abs(child_rms - parent_rms) / max(parent_rms, 1.0e-300))
    assert max(collar_changes) > 1.0e-4

    assert report["physical_support_validated"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


def test_vorticity_operator_fails_closed_on_invalid_steps():
    points = np.asarray([[0.1, 0.2, 0.3]], dtype=float)
    with pytest.raises(ValueError, match="finite-difference step"):
        vorticity_equation_terms(
            _RigidRotation(),
            points,
            0.5,
            spatial_step=0.0,
            time_step=0.0025,
            nu=0.01,
        )
    with pytest.raises(ValueError, match="time_step"):
        vorticity_equation_terms(
            _RigidRotation(),
            points,
            0.5,
            spatial_step=0.01,
            time_step=0.0,
            nu=0.01,
        )
