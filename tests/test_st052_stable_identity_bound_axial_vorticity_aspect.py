from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_axial_vorticity_aspect as mod


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
        "candidate_semantic_identity_sha256": mod.EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": mod.EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    measurement_sha = mod._canonical_sha256(measurement)
    receipt = {
        "schema": mod.SCHEMA,
        "task_id": mod.TASK_ID,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_main_merge": mod.STABLE_IDENTITY_MAIN_MERGE,
            "stable_identity_receipt_sha256": "a" * 64,
            "legacy_whole_candidate_identity_sha256": "b" * 64,
            "materialization_evidence_sha256": "c" * 64,
            "included_in_measurement_binding": False,
        },
        "protocol": mod._protocol(),
        "public_observable": dict(mod.PUBLIC_OBSERVABLE),
        "metric_provenance": dict(mod.METRIC_PROVENANCE),
        "external_method_screen": dict(mod.EXTERNAL_METHOD_SCREEN),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": mod._measurement_binding_sha256(identity, measurement_sha),
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


def test_receipt_requires_admitted_stable_semantic_identity() -> None:
    receipt = _receipt()
    mod.validate_receipt(receipt)
    assert receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"] == mod.EXPECTED_STABLE_CANDIDATE_IDENTITY
    assert receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"] == mod.EXPECTED_STABLE_VELOCITY_IDENTITY
    assert receipt["materialization_join"]["included_in_measurement_binding"] is False
    assert receipt["truth_boundary"]["historical_cr_a9_104_receipt_reused"] is False
    assert receipt["truth_boundary"]["visualization_ready"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_stable_identity_drift_fails_even_after_rehash() -> None:
    bad = copy.deepcopy(_receipt())
    bad["stable_candidate_identity"]["velocity_semantic_identity_sha256"] = "d" * 64
    bad["candidate_measurement_binding_sha256"] = mod._measurement_binding_sha256(
        bad["stable_candidate_identity"], bad["measurement_sha256"]
    )
    _rehash(bad)
    with pytest.raises(ValueError, match="stable candidate/callable identity drift"):
        mod.validate_receipt(bad)


def test_materialization_rehash_does_not_change_semantic_measurement_binding() -> None:
    first = _receipt()
    second = copy.deepcopy(first)
    second["materialization_join"]["stable_identity_receipt_sha256"] = "d" * 64
    second["materialization_join"]["legacy_whole_candidate_identity_sha256"] = "e" * 64
    second["materialization_join"]["materialization_evidence_sha256"] = "f" * 64
    _rehash(second)
    mod.validate_receipt(second)
    assert first["candidate_measurement_binding_sha256"] == second["candidate_measurement_binding_sha256"]
    assert first["receipt_sha256"] != second["receipt_sha256"]


def test_materialization_cannot_be_promoted_into_measurement_identity() -> None:
    bad = copy.deepcopy(_receipt())
    bad["materialization_join"]["included_in_measurement_binding"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="materialization evidence leaked"):
        mod.validate_receipt(bad)


def test_measurement_or_public_target_or_truth_promotion_fails_closed() -> None:
    bad = copy.deepcopy(_receipt())
    bad["measurement"]["metrics"]["full_aspect_ratio"] *= 1.01
    bad["measurement_sha256"] = mod._canonical_sha256(bad["measurement"])
    bad["candidate_measurement_binding_sha256"] = mod._measurement_binding_sha256(
        bad["stable_candidate_identity"], bad["measurement_sha256"]
    )
    _rehash(bad)
    with pytest.raises(ValueError, match="aspect ratio arithmetic drift"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["public_observable"]["numerical_target"] = 1.5
    _rehash(bad)
    with pytest.raises(ValueError, match="public qualitative observable drift"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["truth_boundary"]["visualization_ready"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="truth boundary drift"):
        mod.validate_receipt(bad)


def test_external_screen_is_frozen_and_no_new_method_is_migrated() -> None:
    assert mod.EXTERNAL_METHOD_SCREEN == {
        "source_repo": "pyvista/pyvista",
        "source_commit": "f749a1b0a10a5a4c3c5ca3eedbc8f7860c6f9d9e",
        "license": "MIT",
        "classification": "screened_not_adopted",
        "migration_scope": "none",
        "difference": "direct NumPy callable/grid diagnostics are sufficient; no renderer dependency is needed",
        "new_external_method_migrated": False,
    }
