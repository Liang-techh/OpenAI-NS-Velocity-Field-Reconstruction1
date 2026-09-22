from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_radial_rotation_variation as diag


class RigidSwirl:
    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        x, y, z = np.broadcast_arrays(x, y, z)
        omega = 1.25
        return np.stack((-omega * y, omega * x, np.zeros_like(z)), axis=-1)


class DifferentialSwirl:
    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        x, y, z = np.broadcast_arrays(x, y, z)
        r = np.hypot(x, y)
        omega = 0.9 + 0.35 * r + 0.12 * np.abs(z) + 0.05 * float(t)
        return np.stack((-omega * y, omega * x, np.zeros_like(z)), axis=-1)


def _synthetic_receipt(metrics: dict) -> dict:
    identity = {
        "candidate_id": diag.CANDIDATE_ID,
        "candidate_semantic_identity_sha256": diag.EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": diag.EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    measurement = {
        "time": diag.TIME,
        "metrics": metrics,
        "proxy_interpretation": (
            "nonzero variation means only that the frozen candidate is not radius/spatially rigid under this autonomous ring diagnostic; it is not a numerical OpenAI match or adequacy threshold"
        ),
    }
    measurement_sha = diag._canonical_sha256(measurement)
    receipt = {
        "schema": diag.SCHEMA,
        "task_id": diag.TASK_ID,
        "preregister_issue": diag.PREREG_ISSUE,
        "base_main": diag.BASE_MAIN,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_receipt_sha256": "1" * 64,
            "legacy_whole_candidate_identity_sha256": "2" * 64,
            "materialization_evidence_sha256": "3" * 64,
            "stable_identity_main_merge": diag.STABLE_IDENTITY_MAIN_MERGE,
            "stable_axial_main_merge": diag.STABLE_AXIAL_MAIN_MERGE,
            "included_in_measurement_binding": False,
        },
        "protocol": diag._protocol(),
        "public_observable": diag.PUBLIC_OBSERVABLE,
        "metric_provenance": diag.METRIC_PROVENANCE,
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": diag._measurement_binding_sha256(identity, measurement_sha),
        "decision": diag._decision(metrics),
        "truth_boundary": diag.TRUTH_BOUNDARY,
    }
    receipt["receipt_sha256"] = diag._canonical_sha256(receipt)
    return receipt


def test_rigid_rotation_has_zero_variation_proxy() -> None:
    metrics = diag.rotation_variation_metrics(RigidSwirl())
    assert metrics["angular_rotation_radial_variation_t050"] < 1.0e-12
    assert metrics["nonzero_radius_dependent_rotation_proxy"] is False
    assert len(metrics["samples"]) == 16
    omegas = np.asarray([row["omega_star"] for row in metrics["samples"]], dtype=float)
    assert np.max(np.abs(omegas - omegas[0])) < 1.0e-12


def test_differential_rotation_is_detected() -> None:
    metrics = diag.rotation_variation_metrics(DifferentialSwirl())
    assert metrics["angular_rotation_radial_variation_t050"] > 1.0e-3
    assert metrics["nonzero_radius_dependent_rotation_proxy"] is True
    assert metrics["omega_star_range"] > 0.0
    decision = diag._decision(metrics)
    assert decision["radius_invariant_rotation_deficit_established"] is False
    assert decision["new_basis_authorized"] is False
    assert decision["coefficient_change_authorized"] is False
    assert decision["candidate_mutation_authorized"] is False


def test_cylindrical_transform_recovers_signed_swirl() -> None:
    pts = diag._ring_points(0.8, 0.3)
    omega = 2.0
    vel = np.column_stack((-omega * pts[:, 1], omega * pts[:, 0], np.zeros(len(pts))))
    cyl = diag._cylindrical_components(pts, vel)
    np.testing.assert_allclose(cyl[:, 0], 0.0, atol=2.0e-15, rtol=0.0)
    np.testing.assert_allclose(cyl[:, 1], omega * 0.8, atol=2.0e-15, rtol=0.0)
    np.testing.assert_allclose(cyl[:, 2], 0.0, atol=0.0, rtol=0.0)


def test_protocol_and_truth_boundary_are_frozen() -> None:
    protocol = diag._protocol()
    assert protocol["time"] == 0.50
    assert protocol["times_for_baseline_speed"] == [0.25, 0.50, 0.75]
    assert protocol["ring_radii"] == [0.40, 0.70, 1.00, 1.30]
    assert protocol["ring_abs_z"] == [0.30, 0.80]
    assert protocol["azimuth_count"] == 16
    assert protocol["nonzero_variation_floor"] == 1.0e-10
    assert protocol["source_numeric_targets_used"] is False
    assert diag.TRUTH_BOUNDARY["candidate_changed"] is False
    assert diag.TRUTH_BOUNDARY["basis_dimension_changed"] is False
    assert diag.TRUTH_BOUNDARY["direct_visualization_fingerprint_improvement"] == 0.0
    assert diag.TRUTH_BOUNDARY["pde_validated"] is False


def test_validate_receipt_accepts_consistent_synthetic_receipt_and_rejects_truth_promotion() -> None:
    metrics = diag.rotation_variation_metrics(DifferentialSwirl())
    receipt = _synthetic_receipt(metrics)
    diag.validate_receipt(receipt)

    bad = copy.deepcopy(receipt)
    bad["truth_boundary"]["pde_validated"] = True
    bad["receipt_sha256"] = diag._canonical_sha256({k: v for k, v in bad.items() if k != "receipt_sha256"})
    with pytest.raises(ValueError, match="truth boundary drift"):
        diag.validate_receipt(bad)


def test_validate_receipt_rejects_metric_arithmetic_drift() -> None:
    metrics = diag.rotation_variation_metrics(DifferentialSwirl())
    receipt = _synthetic_receipt(metrics)
    bad = copy.deepcopy(receipt)
    bad_metrics = bad["measurement"]["metrics"]
    bad_metrics["omega_star_range"] += 0.5
    bad["measurement_sha256"] = diag._canonical_sha256(bad["measurement"])
    bad["candidate_measurement_binding_sha256"] = diag._measurement_binding_sha256(
        bad["stable_candidate_identity"], bad["measurement_sha256"]
    )
    bad["decision"] = diag._decision(bad_metrics)
    bad["receipt_sha256"] = diag._canonical_sha256({k: v for k, v in bad.items() if k != "receipt_sha256"})
    with pytest.raises(ValueError, match="range arithmetic drift"):
        diag.validate_receipt(bad)
