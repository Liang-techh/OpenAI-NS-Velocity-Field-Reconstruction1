from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.st054_cylindrical_morphology import (
    EXTERNAL_METHOD_SCREEN,
    TRUTH_BOUNDARY,
    fingerprint_velocity_field,
    validate_fingerprint_receipt,
)
from openai_ns_reconstruction.st054_snapshot import load_st054


class ManufacturedSpiralStretchField:
    """Analytic axisymmetric field with inward, swirl and axial-stretch kinematics."""

    def velocity(self, x, y, z, t):
        x, y, z = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        radius = np.hypot(x, y)
        inward_rate = 0.20 + 0.10 * float(t)
        angular_rate = 0.80 + 0.25 * radius
        axial_rate = 0.35
        u = -inward_rate * x - angular_rate * y
        v = -inward_rate * y + angular_rate * x
        w = axial_rate * z
        return np.stack((u, v, w), axis=-1)


class AcceleratingConcentratingField(ManufacturedSpiralStretchField):
    """Adds a time-dependent mid-plane speed profile for the core proxy regression."""

    def velocity(self, x, y, z, t):
        base = super().velocity(x, y, z, t)
        x, y, z = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        radius = np.hypot(x, y)
        scale = 1.0 + 2.0 * float(t)
        width = 1.25 - 0.55 * float(t)
        gain = scale * np.exp(-(radius / width) ** 2)
        return base * gain[..., None]


def test_cylindrical_projection_recovers_manufactured_kinematics():
    field = ManufacturedSpiralStretchField()
    receipt = fingerprint_velocity_field(
        field,
        radii=(0.5, 1.0),
        z_magnitudes=(0.4,),
        times=(0.5,),
        azimuth_count=16,
        core_radii=(0.25, 0.5, 0.75, 1.0),
    )
    rows = receipt["measurements"]["ring_rows"]
    assert len(rows) == 4
    expected_inward_rate = 0.25
    for row in rows:
        radius = row["radius"]
        expected_omega = 0.80 + 0.25 * radius
        assert row["mean_u_r"] == pytest.approx(-expected_inward_rate * radius, abs=2e-15)
        assert row["std_u_r"] < 2e-15
        assert row["mean_u_theta"] == pytest.approx(expected_omega * radius, abs=2e-15)
        assert row["mean_angular_rate"] == pytest.approx(expected_omega, abs=2e-15)
        assert row["mean_away_from_midplane_u_z"] == pytest.approx(0.35 * 0.4, abs=2e-15)
        assert row["inward_fraction"] == 1.0
        assert row["swirl_nonzero_fraction"] == 1.0
        assert row["inward_swirl_fraction"] == 1.0
        assert row["away_from_midplane_fraction"] == 1.0
        assert row["mean_inward_to_circulation_ratio"] == pytest.approx(
            expected_inward_rate / expected_omega,
            abs=2e-15,
        )

    variations = receipt["measurements"]["radial_variation"]
    assert len(variations) == 2
    assert all(item["circulation_speed_span_across_radii"] > 0.0 for item in variations)
    assert all(item["angular_rate_span_across_radii"] > 0.0 for item in variations)
    assert all(item["all_rings_inward"] for item in variations)


def test_core_proxy_reports_autonomous_trend_without_source_acceptance():
    receipt = fingerprint_velocity_field(
        AcceleratingConcentratingField(),
        radii=(0.5,),
        z_magnitudes=(0.4,),
        times=(0.25, 0.75),
        azimuth_count=16,
        core_radii=tuple(np.linspace(0.15, 1.6, 24)),
    )
    trend = receipt["measurements"]["time_evolution_proxy"]
    assert trend["speed_weighted_radius_decreased"] is True
    assert trend["peak_mean_speed_increased"] is True
    assert receipt["protocol"]["source_numeric_targets_used"] is False
    assert receipt["protocol"]["renderer_or_camera_used"] is False
    assert receipt["protocol"]["pixel_loss_used"] is False
    for key, expected in TRUTH_BOUNDARY.items():
        assert receipt[key] is expected


def test_actual_st054_continuous_callable_can_be_fingerprinted_on_small_protocol():
    field = load_st054("ST054-Q2")
    receipt = fingerprint_velocity_field(
        field,
        model_id="ST054-Q2",
        radii=(0.6, 1.0),
        z_magnitudes=(0.3,),
        times=(0.5,),
        azimuth_count=8,
        core_radii=(0.3, 0.6, 0.9, 1.2),
    )
    assert receipt["source_commit"] == field.source_commit
    assert receipt["source_mat_sha256"] == field.mat_sha256
    assert len(receipt["measurement_sha256"]) == 64
    assert len(receipt["measurements"]["ring_rows"]) == 4
    assert all(
        np.isfinite(row["mean_speed"])
        for row in receipt["measurements"]["ring_rows"]
    )
    assert receipt["cylindrical_morphology_diagnostic_ready"] is True
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False


def test_protocol_rejects_degenerate_or_ambiguous_sampling_inputs():
    field = ManufacturedSpiralStretchField()
    with pytest.raises(ValueError, match="radii"):
        fingerprint_velocity_field(field, radii=(0.0,))
    with pytest.raises(ValueError, match="z_magnitudes"):
        fingerprint_velocity_field(field, z_magnitudes=(0.0,))
    with pytest.raises(ValueError, match="azimuth_count"):
        fingerprint_velocity_field(field, azimuth_count=7)
    with pytest.raises(ValueError, match="times"):
        fingerprint_velocity_field(field, times=(np.nan,))
    with pytest.raises(ValueError, match="core_radii"):
        fingerprint_velocity_field(field, core_radii=(0.0,))


def test_receipt_is_deterministic_and_external_screen_is_non_migration():
    kwargs = dict(
        radii=(0.5, 1.0),
        z_magnitudes=(0.4,),
        times=(0.25, 0.5),
        azimuth_count=8,
        core_radii=(0.25, 0.5, 0.75),
    )
    first = fingerprint_velocity_field(ManufacturedSpiralStretchField(), **kwargs)
    second = fingerprint_velocity_field(ManufacturedSpiralStretchField(), **kwargs)
    assert first["measurement_sha256"] == second["measurement_sha256"]
    assert first["measurements"] == second["measurements"]
    assert EXTERNAL_METHOD_SCREEN[0]["repo"] == "pyvista/pyvista"
    assert EXTERNAL_METHOD_SCREEN[0]["screened_commit"] == "f749a1b0a10a5a4c3c5ca3eedbc8f7860c6f9d9e"
    assert EXTERNAL_METHOD_SCREEN[0]["license"] == "MIT"
    assert EXTERNAL_METHOD_SCREEN[0]["classification"] == "screened_not_adopted"


def test_receipt_verifier_rejects_truth_promotion_target_invention_and_tampering():
    receipt = fingerprint_velocity_field(
        ManufacturedSpiralStretchField(),
        radii=(0.5,),
        z_magnitudes=(0.4,),
        times=(0.5,),
        azimuth_count=8,
        core_radii=(0.25, 0.5),
    )
    validate_fingerprint_receipt(receipt)

    promoted = copy.deepcopy(receipt)
    promoted["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="truth boundary promoted"):
        validate_fingerprint_receipt(promoted)

    invented_target = copy.deepcopy(receipt)
    invented_target["protocol"]["source_numeric_targets_used"] = True
    with pytest.raises(ValueError, match="source_numeric_targets_used"):
        validate_fingerprint_receipt(invented_target)

    pixelized = copy.deepcopy(receipt)
    pixelized["protocol"]["pixel_loss_used"] = True
    with pytest.raises(ValueError, match="pixel_loss_used"):
        validate_fingerprint_receipt(pixelized)

    migrated = copy.deepcopy(receipt)
    migrated["external_method_screen"][0]["classification"] = "direct_migration"
    with pytest.raises(ValueError, match="external method screen drift"):
        validate_fingerprint_receipt(migrated)

    tampered = copy.deepcopy(receipt)
    tampered["measurements"]["ring_rows"][0]["mean_u_r"] += 1.0
    with pytest.raises(ValueError, match="digest mismatch"):
        validate_fingerprint_receipt(tampered)
