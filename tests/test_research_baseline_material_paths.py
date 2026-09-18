from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from research_baseline.material_paths import (
    CANDIDATE_SHA256,
    OUTPUT_SAMPLES,
    SEED_ANGLES,
    SEED_RADII,
    SEED_Z,
    _canonical_hash,
    measure_material_paths,
)


def _linear_helical_velocity(a: float, omega: float, b: float):
    def velocity(points, time):
        del time
        p = np.asarray(points, dtype=float)
        out = np.empty_like(p)
        out[:, 0] = -a * p[:, 0] - omega * p[:, 1]
        out[:, 1] = omega * p[:, 0] - a * p[:, 1]
        out[:, 2] = b * p[:, 2]
        return out

    return velocity


def _measure(velocity):
    return measure_material_paths(
        velocity,
        candidate="analytic-test-field",
        candidate_sha256="1" * 64,
        provenance="analytic regression",
    )


def test_contract_recovers_inward_helical_motion_and_axial_stretching():
    a, omega, b = 0.4, 4.0, 0.6
    report = _measure(_linear_helical_velocity(a, omega, b))
    dt = 0.5
    expected_radius_factor = math.exp(-a * dt)
    expected_turns = omega * dt / (2.0 * math.pi)
    expected_separation_ratio = math.exp(b * dt)

    assert report["registered_contract"]["seed_radii"] == list(SEED_RADII)
    assert report["registered_contract"]["seed_z"] == list(SEED_Z)
    assert report["registered_contract"]["seed_angles"] == SEED_ANGLES
    assert report["registered_contract"]["output_samples"] == OUTPUT_SAMPLES
    assert report["summary"]["path_count"] == len(SEED_RADII) * SEED_ANGLES * len(SEED_Z)
    assert report["summary"]["inward_path_count"] == report["summary"]["path_count"]
    assert report["summary"]["pair_axial_separation_growth_count"] == len(SEED_RADII) * SEED_ANGLES
    assert report["summary"]["mean_absolute_turns"] == pytest.approx(expected_turns, rel=2e-8, abs=2e-10)
    assert report["summary"]["mean_pair_axial_separation_ratio"] == pytest.approx(
        expected_separation_ratio, rel=2e-8, abs=2e-10
    )
    for path in report["per_path"]:
        assert path["radius_end"] / path["radius_start"] == pytest.approx(
            expected_radius_factor, rel=2e-8, abs=2e-10
        )
    assert report["interpretation_boundary"]["pde_validated"] is False
    assert report["interpretation_boundary"]["visual_correspondence_verified"] is False
    assert report["interpretation_boundary"]["openai_seed_locations_used"] is False


def test_rigid_rotation_does_not_fake_contraction_or_axial_stretching():
    omega = 3.0
    report = _measure(_linear_helical_velocity(0.0, omega, 0.0))
    summary = report["summary"]
    assert summary["mean_radius_change"] == pytest.approx(0.0, abs=2e-10)
    assert summary["mean_pair_axial_separation_change"] == pytest.approx(0.0, abs=2e-10)
    assert summary["mean_pair_axial_separation_ratio"] == pytest.approx(1.0, abs=2e-10)
    assert summary["mean_absolute_turns"] == pytest.approx(omega * 0.5 / (2.0 * math.pi), rel=2e-8)


def test_report_hash_is_deterministic_and_binds_payload():
    report1 = _measure(_linear_helical_velocity(0.2, 2.0, 0.1))
    report2 = _measure(_linear_helical_velocity(0.2, 2.0, 0.1))
    assert report1["report_sha256"] == report2["report_sha256"]
    payload = copy.deepcopy(report1)
    digest = payload.pop("report_sha256")
    assert digest == _canonical_hash(payload)
    payload["summary"]["mean_radius_change"] += 1e-6
    assert digest != _canonical_hash(payload)


def test_fail_closed_for_malformed_nonfinite_and_zero_velocity():
    with pytest.raises(ValueError, match="shape"):
        _measure(lambda points, time: np.zeros((len(points), 2)))

    def nonfinite(points, time):
        out = np.ones_like(points, dtype=float)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        _measure(nonfinite)
    with pytest.raises(ValueError, match="numerically zero"):
        _measure(lambda points, time: np.zeros_like(points, dtype=float))


def test_fail_closed_for_bad_identity_and_provenance():
    velocity = _linear_helical_velocity(0.1, 1.0, 0.1)
    with pytest.raises(ValueError, match="candidate_sha256"):
        measure_material_paths(velocity, candidate="x", candidate_sha256="bad", provenance="test")
    with pytest.raises(ValueError, match="hexadecimal"):
        measure_material_paths(velocity, candidate="x", candidate_sha256="z" * 64, provenance="test")
    with pytest.raises(ValueError, match="candidate"):
        measure_material_paths(velocity, candidate="", candidate_sha256=CANDIDATE_SHA256, provenance="test")
    with pytest.raises(ValueError, match="provenance"):
        measure_material_paths(velocity, candidate="x", candidate_sha256=CANDIDATE_SHA256, provenance="")
