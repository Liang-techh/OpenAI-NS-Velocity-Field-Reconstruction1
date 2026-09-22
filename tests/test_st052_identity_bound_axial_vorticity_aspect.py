from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction import st052_identity_bound_axial_vorticity_aspect as mod


class _SolidBodyRotation:
    """Manufactured field with constant omega_z=1 on the whole cube."""

    def velocity(self, x, y, z, t):
        del z, t
        xx = np.asarray(x, dtype=float)
        yy = np.asarray(y, dtype=float)
        return np.stack((-0.5 * yy, 0.5 * xx, np.zeros_like(xx)), axis=-1)


def _receipt() -> dict:
    measurement = mod.measure_field(_SolidBodyRotation())
    identity = {
        "candidate_id": mod.CANDIDATE_ID,
        "candidate_sha256": "a" * 64,
        "velocity_identity_sha256": "b" * 64,
    }
    measurement_sha = mod._canonical_sha256(measurement)
    receipt = {
        "schema": mod.SCHEMA,
        "task_id": mod.TASK_ID,
        "candidate_identity": identity,
        "identity_ancestry": {},
        "protocol": {
            "time": mod.TIME,
            "grid_resolution": mod.GRID_RESOLUTION,
            "box": list(mod.BOX),
            "cartesian_curl_order": 2,
            "integration_rule": "tensor_product_trapezoid",
            "source_numeric_targets_used": False,
            "renderer_or_camera_used": False,
            "pixel_loss_used": False,
        },
        "public_observable": dict(mod.PUBLIC_OBSERVABLE),
        "metric_provenance": dict(mod.METRIC_PROVENANCE),
        "external_method_screen": dict(mod.EXTERNAL_METHOD_SCREEN),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": mod._measurement_binding_sha256(
            identity, measurement_sha
        ),
        "truth_boundary": dict(mod.TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = mod._canonical_sha256(receipt)
    return receipt


def _rehash(receipt: dict) -> dict:
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = mod._canonical_sha256(receipt)
    return receipt


def test_constant_axial_vorticity_calibrates_aspect_metric() -> None:
    measurement = mod.measure_field(_SolidBodyRotation())
    metrics = measurement["metrics"]
    assert metrics["full_aspect_ratio"] == pytest.approx(1.0 / math.sqrt(2.0), abs=2e-14)
    assert metrics["full_axial_rms"] > 0.0
    assert metrics["full_radial_rms"] > metrics["full_axial_rms"]
    assert metrics["axially_elongated_proxy"] is False
    assert metrics["full_enstrophy_trapezoid"] > 0.0


def test_curl_detects_manufactured_constant_omega_z() -> None:
    axis, velocity = mod.sample_velocity_grid(_SolidBodyRotation())
    omega, magnitude = mod.cartesian_vorticity(velocity, float(axis[1] - axis[0]))
    assert np.max(np.abs(omega[..., 0])) < 2e-14
    assert np.max(np.abs(omega[..., 1])) < 2e-14
    assert np.max(np.abs(omega[..., 2] - 1.0)) < 2e-14
    assert np.max(np.abs(magnitude - 1.0)) < 2e-14


def test_receipt_binds_candidate_metric_and_keeps_truth_false() -> None:
    receipt = _receipt()
    mod.validate_receipt(receipt)
    assert receipt["public_observable"]["numerical_target"] is None
    assert receipt["external_method_screen"]["classification"] == "screened_not_adopted"
    assert receipt["external_method_screen"]["new_external_method_migrated"] is False
    assert receipt["truth_boundary"]["axial_vorticity_aspect_measured"] is True
    assert receipt["truth_boundary"]["axial_stretching_correspondence_verified"] is False
    assert receipt["truth_boundary"]["visualization_ready"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_candidate_or_measurement_drift_fails_closed() -> None:
    bad = copy.deepcopy(_receipt())
    bad["candidate_identity"]["candidate_sha256"] = "c" * 64
    _rehash(bad)
    with pytest.raises(ValueError, match="binding mismatch"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["measurement"]["metrics"]["full_aspect_ratio"] *= 1.01
    bad["measurement_sha256"] = mod._canonical_sha256(bad["measurement"])
    bad["candidate_measurement_binding_sha256"] = mod._measurement_binding_sha256(
        bad["candidate_identity"], bad["measurement_sha256"]
    )
    _rehash(bad)
    with pytest.raises(ValueError, match="aspect ratio arithmetic drift"):
        mod.validate_receipt(bad)


def test_source_target_or_truth_promotion_fails_closed() -> None:
    bad = copy.deepcopy(_receipt())
    bad["public_observable"]["numerical_target"] = 1.5
    _rehash(bad)
    with pytest.raises(ValueError, match="public qualitative observable drift"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["protocol"]["source_numeric_targets_used"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="protocol drift"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["truth_boundary"]["visualization_ready"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="truth boundary drift"):
        mod.validate_receipt(bad)


def test_metric_and_external_migration_provenance_are_frozen() -> None:
    assert mod.METRIC_PROVENANCE["source_pr"] == 1069
    assert mod.METRIC_PROVENANCE["source_head"] == "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7"
    assert mod.METRIC_PROVENANCE["source_blob_sha"] == "862f5b25554b1a6e7c038418092e0b20d9435bef"
    assert mod.METRIC_PROVENANCE["metric_origin_blob_sha"] == "c072712ced553eb6aef09a5789f53d5cea1b7aa3"
    assert mod.EXTERNAL_METHOD_SCREEN == {
        "source_repo": "pyvista/pyvista",
        "source_commit": "f749a1b0a10a5a4c3c5ca3eedbc8f7860c6f9d9e",
        "license": "MIT",
        "classification": "screened_not_adopted",
        "migration_scope": "none",
        "difference": "direct NumPy callable/grid diagnostics are sufficient; no renderer dependency is needed",
        "new_external_method_migrated": False,
    }
