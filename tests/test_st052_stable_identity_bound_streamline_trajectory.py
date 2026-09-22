from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_streamline_trajectory as target


class LinearHelixField:
    def __init__(self, *, radial_rate: float = -0.30, swirl: float = 1.0, axial: float = 0.20):
        self.radial_rate = float(radial_rate)
        self.swirl = float(swirl)
        self.axial = float(axial)

    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        u = self.radial_rate * x - self.swirl * y
        v = self.radial_rate * y + self.swirl * x
        w = self.axial * z
        return np.column_stack((u, v, w))


def _synthetic_receipt(metrics):
    identity = {
        "candidate_id": target.CANDIDATE_ID,
        "candidate_semantic_identity_sha256": target.EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": target.EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    measurement = {
        "time": target.TIME,
        "metrics": metrics,
        "interpretation": "target-free candidate-side fixed-seed trajectory structure; not an OpenAI magnitude/profile match",
    }
    measurement_sha = target._canonical_sha256(measurement)
    receipt = {
        "schema": target.SCHEMA,
        "task_id": target.TASK_ID,
        "preregister_issue": target.PREREG_ISSUE,
        "base_main": target.BASE_MAIN,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_receipt_sha256": "0" * 64,
            "legacy_whole_candidate_identity_sha256": "1" * 64,
            "materialization_evidence_sha256": "2" * 64,
            "included_in_measurement_binding": False,
            "stable_identity_main_merge": target.STABLE_IDENTITY_MAIN_MERGE,
            "stable_axial_main_merge": target.STABLE_AXIAL_MAIN_MERGE,
        },
        "protocol": target._protocol(),
        "public_observables": list(target.PUBLIC_OBSERVABLES),
        "metric_provenance": target.METRIC_PROVENANCE,
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": target._measurement_binding_sha256(identity, measurement_sha),
        "decision": target._decision(metrics),
        "truth_boundary": target.TRUTH_BOUNDARY,
        "receipt_sha256": None,
    }
    unsigned = dict(receipt)
    unsigned.pop("receipt_sha256")
    receipt["receipt_sha256"] = target._canonical_sha256(unsigned)
    return receipt


def test_frozen_seed_and_integrator_contract():
    seeds = target.seed_points()
    assert seeds.shape == (48, 3)
    assert set(np.round(np.hypot(seeds[:, 0], seeds[:, 1]), 12)) == {0.6, 0.9, 1.2}
    assert set(np.round(seeds[:, 2], 12)) == {-0.3, 0.3}

    axis, grid = target.sample_velocity_grid(LinearHelixField())
    full, forward = target.integrate_streamlines(target.grid_interpolator(axis, grid), seeds)
    assert full.shape == (48, 181, 3)
    assert forward.shape == (48, 91, 3)
    assert np.isfinite(full).all()


def test_inward_helical_field_closes_only_coarse_structure_triggers():
    metrics = target.trajectory_metrics(LinearHelixField())
    proxies = metrics["proxies"]
    assert proxies["inward_forward_trajectory_proxy"] is True
    assert proxies["nonzero_spiral_winding_proxy"] is True
    assert proxies["nonplanar_axial_trajectory_proxy"] is True
    assert proxies["coarse_inward_spiraling_axial_structure_proxy"] is True
    assert metrics["forward_overall"]["mean_endpoint_radial_displacement"] < 0.0
    assert metrics["forward_upper_seeds"]["mean_endpoint_radial_displacement"] < 0.0
    assert metrics["forward_lower_seeds"]["mean_endpoint_radial_displacement"] < 0.0

    decision = target._decision(metrics)
    assert decision["candidate_side_trajectory_structure_defect_established"] is False
    assert decision["new_basis_authorized"] is False
    assert decision["candidate_mutation_authorized"] is False
    assert decision["direct_visualization_fingerprint_improvement"] == 0.0


def test_planar_or_outward_fields_expose_specific_structure_defects():
    planar = target.trajectory_metrics(LinearHelixField(axial=0.0))
    assert planar["proxies"]["nonplanar_axial_trajectory_proxy"] is False
    assert "nonplanar_axial_trajectory_proxy" in target._decision(planar)["failed_structural_proxies"]

    outward = target.trajectory_metrics(LinearHelixField(radial_rate=+0.30))
    assert outward["proxies"]["inward_forward_trajectory_proxy"] is False
    assert "inward_forward_trajectory_proxy" in target._decision(outward)["failed_structural_proxies"]
    assert target._decision(outward)["new_basis_authorized"] is False


def test_receipt_validation_freezes_identity_protocol_and_truth_boundary():
    metrics = target.trajectory_metrics(LinearHelixField())
    receipt = _synthetic_receipt(metrics)
    target.validate_receipt(receipt)

    identity_drift = copy.deepcopy(receipt)
    identity_drift["stable_candidate_identity"]["velocity_semantic_identity_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="identity drift"):
        target.validate_receipt(identity_drift)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth boundary drift"):
        target.validate_receipt(promoted)


def test_receipt_validation_rejects_unbound_metric_mutation():
    metrics = target.trajectory_metrics(LinearHelixField())
    receipt = _synthetic_receipt(metrics)
    mutated = copy.deepcopy(receipt)
    mutated["measurement"]["metrics"]["forward_overall"]["mean_endpoint_radial_displacement"] += 0.01
    with pytest.raises(ValueError, match="routing decision drift|measurement checksum mismatch"):
        target.validate_receipt(mutated)
