from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_core_contraction_speed as target


class _PrescribedSwirl:
    """Axisymmetric manufactured field with time-prescribed radial concentration."""

    def __init__(self, sigma_by_time: dict[float, float], amplitude_by_time: dict[float, float]):
        self._sigma = {float(k): float(v) for k, v in sigma_by_time.items()}
        self._amplitude = {float(k): float(v) for k, v in amplitude_by_time.items()}

    def velocity(self, x, y, z, t):
        del z
        time = float(t)
        sigma = self._sigma[time]
        amplitude = self._amplitude[time]
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        r = np.hypot(x, y)
        phi = np.arctan2(y, x)
        u_theta = amplitude * r * np.exp(-np.square(r / sigma))
        return np.stack(
            (-u_theta * np.sin(phi), u_theta * np.cos(phi), np.zeros_like(r)),
            axis=-1,
        )


def _monotone_field() -> _PrescribedSwirl:
    return _PrescribedSwirl(
        {0.25: 1.20, 0.50: 0.80, 0.75: 0.45},
        {0.25: 1.0, 0.50: 10.0, 0.75: 100.0},
    )


def _endpoint_only_field() -> _PrescribedSwirl:
    return _PrescribedSwirl(
        {0.25: 1.00, 0.50: 1.45, 0.75: 0.50},
        {0.25: 1.0, 0.50: 4.0, 0.75: 40.0},
    )


def test_frozen_protocol_and_truth_boundary() -> None:
    p = target.protocol()
    assert p["times"] == [0.25, 0.5, 0.75]
    assert p["z"] == 0.0
    assert p["core_radii"] == list(target.CORE_RADII)
    assert len(p["core_radii"]) == 15
    assert p["core_radii"][0] == pytest.approx(0.2)
    assert p["core_radii"][-1] == pytest.approx(1.6)
    assert p["azimuth_count"] == 16
    assert p["source_numeric_targets_used"] is False
    assert p["renderer_or_camera_used"] is False
    assert p["pixel_loss_used"] is False
    assert target.PUBLIC_OBSERVABLE["numerical_target"] is None
    assert target.TRUTH_BOUNDARY["candidate_changed"] is False
    assert target.TRUTH_BOUNDARY["basis_dimension_changed"] is False
    assert target.TRUTH_BOUNDARY["coefficient_selected"] is False
    assert target.TRUTH_BOUNDARY["visualization_fingerprint_direct_improvement"] == 0.0
    assert target.TRUTH_BOUNDARY["visualization_ready"] is False
    assert target.TRUTH_BOUNDARY["pde_validated"] is False


def test_monotone_manufactured_core_closes_only_coarse_absence_trigger() -> None:
    measurement = target.measure_field(_monotone_field())
    target.validate_measurement(measurement)
    d = measurement["derived"]
    assert d["endpoint_width_fraction"] < 0.0
    assert d["endpoint_peak_speed_fraction"] > 0.0
    assert d["endpoint_shrink_and_speedup_proxy"] is True
    assert d["three_time_width_monotone_decrease"] is True
    assert d["three_time_peak_speed_monotone_increase"] is True
    assert d["three_time_monotone_proxy"] is True
    routing = measurement["routing"]
    assert routing["decision"] == "close_coarse_absence_trigger_no_temporal_basis_growth_from_this_presence_proxy"
    assert routing["candidate_mutation_authorized"] is False
    assert routing["sixth_basis_authorized"] is False
    assert routing["coefficient_selection_authorized"] is False
    assert routing["failures"] == []


def test_midpoint_shape_defect_routes_to_existing_temporal_capacity() -> None:
    measurement = target.measure_field(_endpoint_only_field())
    target.validate_measurement(measurement)
    d = measurement["derived"]
    assert d["endpoint_shrink_and_speedup_proxy"] is True
    assert d["three_time_monotone_proxy"] is False
    routing = measurement["routing"]
    assert routing["decision"] == "route_temporal_shape_defect_to_existing_temporal_curvature_and_offgrid_time_skewed_capacity"
    assert "three_time_midpoint_monotonicity_failed" in routing["failures"]
    assert routing["candidate_mutation_authorized"] is False
    assert routing["sixth_basis_authorized"] is False


def test_endpoint_sign_failure_does_not_authorize_sixth_basis() -> None:
    routing = target._routing(
        endpoint_proxy=False,
        monotone_proxy=False,
        width_decreased=False,
        speed_increased=True,
    )
    assert routing["decision"] == "route_endpoint_sign_defect_to_existing_five_control_temporal_sensitivity"
    assert "speed_weighted_core_radius_did_not_decrease" in routing["failures"]
    assert routing["sixth_basis_authorized"] is False
    assert routing["coefficient_selection_authorized"] is False


def test_measurement_arithmetic_and_routing_tampering_fail_closed() -> None:
    measurement = target.measure_field(_monotone_field())

    bad_fraction = copy.deepcopy(measurement)
    bad_fraction["derived"]["endpoint_width_fraction"] += 1.0e-4
    with pytest.raises(ValueError, match="width-fraction arithmetic drift"):
        target.validate_measurement(bad_fraction)

    bad_route = copy.deepcopy(measurement)
    bad_route["routing"]["sixth_basis_authorized"] = True
    with pytest.raises(ValueError, match="routing drift"):
        target.validate_measurement(bad_route)


def test_semantic_measurement_binding_excludes_materialization_evidence() -> None:
    identity = {
        "candidate_id": target.CANDIDATE_ID,
        "candidate_semantic_identity_sha256": target.EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": target.EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    measurement_sha = "a" * 64
    digest = target._measurement_binding_sha256(identity, measurement_sha)
    assert isinstance(digest, str) and len(digest) == 64

    # Materialization hashes are intentionally not arguments to the semantic binding.
    assert digest == target._measurement_binding_sha256(dict(identity), measurement_sha)
