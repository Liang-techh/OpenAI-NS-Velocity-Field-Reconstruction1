from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_tip_taper_proxy as mod


def _synthetic_magnitude(*, tapered: bool) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(-2.0, 2.0, mod.GRID_RESOLUTION, dtype=float)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    rr = np.hypot(xx, yy)
    az = np.abs(zz) / mod.SUPPORT_HALF_HEIGHT
    if tapered:
        sigma = np.where(az >= 0.50, 0.34, np.where(az <= 0.35, 0.92, 0.62))
    else:
        sigma = np.where(az >= 0.50, 1.02, np.where(az <= 0.35, 0.42, 0.66))
    magnitude = np.exp(-np.square(rr / sigma)) * (0.8 + 0.2 * np.cos(0.4 * zz) ** 2)
    return axis, magnitude


def _measurement(*, tapered: bool) -> dict:
    axis, magnitude = _synthetic_magnitude(tapered=tapered)
    return {
        "time": mod.TIME,
        "grid_resolution": mod.GRID_RESOLUTION,
        "box": list(mod.BOX),
        "spacing": float(axis[1] - axis[0]),
        "metrics": mod.tip_taper_metrics(magnitude, axis),
        "proxy_interpretation": "synthetic regression fixture only",
    }


def _receipt(*, tapered: bool = True) -> dict:
    measurement = _measurement(tapered=tapered)
    identity = {
        "candidate_id": mod.CANDIDATE_ID,
        "candidate_semantic_identity_sha256": mod.EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": mod.EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    measurement_sha = mod._canonical_sha256(measurement)
    receipt = {
        "schema": mod.SCHEMA,
        "task_id": mod.TASK_ID,
        "prereg_issue": mod.PREREG_ISSUE,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_receipt_sha256": "a" * 64,
            "legacy_whole_candidate_identity_sha256": "b" * 64,
            "materialization_evidence_sha256": "c" * 64,
            "included_in_measurement_binding": False,
        },
        "protocol": mod._protocol(),
        "metric_provenance": dict(mod.METRIC_PROVENANCE),
        "public_context": dict(mod.PUBLIC_CONTEXT),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": mod._measurement_binding_sha256(identity, measurement_sha),
        "routing": mod.routing_from_metrics(measurement["metrics"]),
        "truth_boundary": dict(mod.TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = mod._canonical_sha256(receipt)
    return receipt


def _rehash(receipt: dict) -> dict:
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = mod._canonical_sha256(receipt)
    return receipt


def _remeasure(receipt: dict) -> None:
    receipt["measurement_sha256"] = mod._canonical_sha256(receipt["measurement"])
    receipt["candidate_measurement_binding_sha256"] = mod._measurement_binding_sha256(
        receipt["stable_candidate_identity"], receipt["measurement_sha256"]
    )
    receipt["routing"] = mod.routing_from_metrics(receipt["measurement"]["metrics"])
    _rehash(receipt)


def test_symmetric_narrow_tip_fixture_is_detected_as_tapered() -> None:
    axis, magnitude = _synthetic_magnitude(tapered=True)
    metrics = mod.tip_taper_metrics(magnitude, axis)
    assert metrics["tip_radial_rms_upper"] == pytest.approx(
        metrics["tip_radial_rms_lower"], rel=0.0, abs=2e-14
    )
    assert metrics["worst_tip_to_central_ratio"] < 1.0
    assert metrics["both_tips_radially_tapered_proxy"] is True
    assert metrics["upper_lower_tip_relative_mismatch"] < 1e-13


def test_wide_tip_fixture_is_detected_without_promoting_basis_growth() -> None:
    axis, magnitude = _synthetic_magnitude(tapered=False)
    metrics = mod.tip_taper_metrics(magnitude, axis)
    assert metrics["worst_tip_to_central_ratio"] > 1.0
    assert metrics["both_tips_radially_tapered_proxy"] is False
    routing = mod.routing_from_metrics(metrics)
    assert routing["classification"] == "candidate_side_tip_taper_deficit_established"
    assert routing["tip_specific_basis_growth_authorized"] is False
    assert routing["axial_turnover_change_authorized"] is False
    assert routing["next_minimal_action"] == "route_to_existing_tip_thickness_control_sensitivity_before_basis_growth"


def test_tapered_route_closes_blunt_tip_trigger_only() -> None:
    metrics = _measurement(tapered=True)["metrics"]
    routing = mod.routing_from_metrics(metrics)
    assert routing["classification"] == "no_candidate_side_blunt_tip_deficit_established"
    assert routing["tip_specific_basis_growth_authorized"] is False
    assert routing["axial_turnover_change_authorized"] is False
    assert routing["next_minimal_action"] == "require_distinct_stable_identity_bound_trajectory_or_return_flow_discrepancy"


def test_receipt_requires_stable_semantic_identity_and_truth_boundary() -> None:
    receipt = _receipt(tapered=True)
    mod.validate_receipt(receipt)
    assert receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"] == mod.EXPECTED_STABLE_CANDIDATE_IDENTITY
    assert receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"] == mod.EXPECTED_STABLE_VELOCITY_IDENTITY
    assert receipt["public_context"]["public_numerical_target"] is None
    assert receipt["truth_boundary"]["candidate_velocity_changed"] is False
    assert receipt["truth_boundary"]["basis_dimension_changed"] is False
    assert receipt["truth_boundary"]["direct_visualization_fingerprint_improvement"] == 0.0
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


def test_materialization_evidence_cannot_enter_semantic_measurement_binding() -> None:
    bad = copy.deepcopy(_receipt())
    bad["materialization_join"]["included_in_measurement_binding"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="materialization evidence leaked"):
        mod.validate_receipt(bad)


def test_metric_arithmetic_mutation_fails_closed() -> None:
    bad = copy.deepcopy(_receipt())
    bad["measurement"]["metrics"]["worst_tip_to_central_ratio"] *= 1.01
    _remeasure(bad)
    with pytest.raises(ValueError, match="tip/central ratio arithmetic drift"):
        mod.validate_receipt(bad)


def test_truth_or_public_target_promotion_fails_closed() -> None:
    bad = copy.deepcopy(_receipt())
    bad["truth_boundary"]["visualization_ready"] = True
    _rehash(bad)
    with pytest.raises(ValueError, match="truth boundary drift"):
        mod.validate_receipt(bad)

    bad = copy.deepcopy(_receipt())
    bad["public_context"]["public_numerical_target"] = 0.75
    _rehash(bad)
    with pytest.raises(ValueError, match="public-context boundary drift"):
        mod.validate_receipt(bad)


def test_empty_or_invalid_weight_regions_fail_closed() -> None:
    axis = np.linspace(-2.0, 2.0, mod.GRID_RESOLUTION)
    zero = np.zeros((mod.GRID_RESOLUTION,) * 3)
    with pytest.raises(ValueError, match="zero full-grid enstrophy"):
        mod.tip_taper_metrics(zero, axis)

    _, magnitude = _synthetic_magnitude(tapered=True)
    magnitude = magnitude.copy()
    magnitude[0, 0, 0] = -1.0
    with pytest.raises(ValueError, match="nonnegative magnitude"):
        mod.tip_taper_metrics(magnitude, axis)
