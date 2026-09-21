from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_a4_current_rf40_lambda_turn_leading_oscillatory_divergence_independent_audit as mod,
)


def _stat(weighted: bool = True) -> dict:
    out = {
        "sampled_max": 1.0e-8,
        "rms": 5.0e-9,
        "normalized_sampled_max": 1.0e-8,
        "normalized_rms": 5.0e-9,
        "worst_index": 0,
        "worst_divergence": 1.0e-8,
        "worst_normalized_divergence": 1.0e-8,
    }
    if weighted:
        out.update({
            "estimated_volume": 1.0,
            "weighted_rms": 5.0e-9,
            "volume_l2_estimate": 5.0e-9,
            "normalized_weighted_rms": 5.0e-9,
        })
    return out


def _region(weighted: bool = True) -> dict:
    return {
        "total": _stat(weighted),
        "leading": _stat(weighted),
        "oscillatory_increment": _stat(weighted),
        "total_speed_rms": 1.0,
        "oscillatory_speed_rms": 1.0e-4,
        "oscillatory_speed_abs_max": 2.0e-4,
    }


def _synthetic_receipt() -> dict:
    ladder = []
    for step in mod.STEPS:
        ladder.append({
            "step": step,
            "inner": _region(True),
            "lambda_turn": _region(True),
            "X2_seam": _region(False),
            "axis": _region(False),
        })
    return {
        "schema": mod.SCHEMA,
        "task": mod.TASK,
        "seed": mod.SEED,
        "steps": list(mod.STEPS),
        "times": list(mod.TIMES),
        "upstream": {
            "pr": mod.UPSTREAM_PR,
            "head": mod.UPSTREAM_HEAD,
            "source_blob_sha1": mod.UPSTREAM_SOURCE_BLOB,
            "semantic_sha256": "1" * 64,
        },
        "agent1": {
            "pr": mod.AGENT1_PR,
            "head": mod.AGENT1_HEAD,
            "source_blob_sha1": mod.AGENT1_SOURCE_BLOB,
            "semantic_sha256": "2" * 64,
        },
        "probe_counts": {"inner": 32, "lambda_turn": 72, "X2_seam": 16, "axis": 5},
        "resolution_ladder": ladder,
        "firewalls": {
            "save_load_velocity_replay_max_abs": 0.0,
            "divergence_mutation_detected": 1.0e-3,
            "ordering_aggregate_max_abs": 0.0,
            "manufactured_solenoidal_divergence_max_abs": 0.0,
            "agent1_semantic_mutation_detected": True,
            "oscillatory_runtime_mutation_rejected": True,
            "beyond_X3_fail_closed": True,
        },
        "frozen_gates": {},
        "truth_boundary": dict(mod.TRUTH_BOUNDARY),
        "receipt_sha256": "0" * 64,
    }


def test_exact_target_and_truth_boundary_are_lambda_turn_only() -> None:
    assert mod.TASK == "K4-VAL-101"
    assert mod.UPSTREAM_PR == 1004
    assert mod.UPSTREAM_HEAD == "e9697cdfe2249e1b8bc2f9241bd16d456c60418f"
    assert mod.AGENT1_PR == 998
    assert mod.PARENT_A4_PR == 1001
    truth = mod.TRUTH_BOUNDARY
    assert truth["RF40_lambda_turn_scoped_divergence_assessed"] is True
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["RF40_power_law_composite_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["agent3_correction_velocity_materialized"] is False
    assert truth["pde_validated"] is False


def test_frozen_protocol_has_three_resolutions_and_fixed_gate() -> None:
    assert mod.STEPS == (0.02, 0.01, 0.005)
    assert mod.SEED == 9173731
    assert mod.DIVERGENCE_GATE == 1.0e-5
    assert mod.LAMBDA_TURN_ZONES == (
        ("lambda_turn_early", 0.08, 0.28),
        ("lambda_turn_middle", 0.38, 0.62),
        ("lambda_turn_late", 0.72, 0.92),
    )


def test_probe_generator_is_strictly_inside_lambda_turn() -> None:
    class Dummy:
        X_2 = 3.0
        X_3 = 3.0 * math.e
        D = 0.31

    p = mod.make_probes(Dummy())
    assert len(p.inner.points) == 32
    assert len(p.lambda_turn.points) == 72
    assert len(p.seam_points) == 16
    assert len(p.axis_points) == 5
    assert np.all(np.isfinite(p.lambda_turn.points))
    assert np.all(p.lambda_turn.weights > 0.0)

    # Invert the public source map enough to recover X for the held-out points.
    r2 = np.sum(p.lambda_turn.points[:, :2] ** 2, axis=1)
    z = p.lambda_turn.points[:, 2]
    t = p.lambda_turn.times
    # eta is not directly stored. Recover it monotonically by scalar bisection.
    recovered = []
    for rr, zz, tt in zip(r2, z, t):
        lo, hi = -0.9, 0.9
        for _ in range(80):
            eta = 0.5 * (lo + hi)
            q = (1.0 - tt) / (1.0 - eta * eta)
            f = q ** Dummy.D * eta - zz
            if f < 0.0:
                lo = eta
            else:
                hi = eta
        eta = 0.5 * (lo + hi)
        q = (1.0 - tt) / (1.0 - eta * eta)
        recovered.append(rr / (2.0 * q))
    recovered = np.asarray(recovered)
    assert np.all(recovered > Dummy.X_2)
    assert np.all(recovered < Dummy.X_3)


def test_fd2_calibrates_on_time_dependent_solenoidal_linear_field() -> None:
    field = mod._ManufacturedSolenoidal()
    points = np.asarray([[0.2, -0.1, 0.3], [-0.4, 0.2, -0.2], [0.0, 0.0, 0.1]])
    times = np.asarray([0.31, 0.47, 0.71])
    for h in mod.STEPS:
        J = mod.fd2_jacobian(field, points, times, h)
        assert np.max(np.abs(np.trace(J, axis1=1, axis2=2))) < 1.0e-12


def test_divergence_mutation_is_detected_analytically() -> None:
    field = mod._ManufacturedSolenoidal()
    mutated = mod._DivergenceMutation(field)
    points = np.asarray([[0.2, -0.1, 0.3], [-0.4, 0.2, -0.2]])
    times = np.asarray([0.31, 0.71])
    base = mod.fd2_jacobian(field, points, times, mod.STEPS[-1])
    child = mod.fd2_jacobian(mutated, points, times, mod.STEPS[-1])
    delta = np.trace(child, axis1=1, axis2=2) - np.trace(base, axis1=1, axis2=2)
    assert np.allclose(delta, mod.MUTATION_EPSILON, rtol=0.0, atol=1.0e-12)


def test_recursive_eta_fd_mutation_changes_only_requested_leaf() -> None:
    payload = {"a": {"b": [{"eta_fd_step": 2.0e-5, "keep": 7}]}}
    assert mod._recursive_eta_fd_perturb(payload) is True
    assert payload["a"]["b"][0]["eta_fd_step"] == pytest.approx(2.002e-5)
    assert payload["a"]["b"][0]["keep"] == 7


def test_synthetic_receipt_accepts_frozen_scope_and_rejects_laundering() -> None:
    receipt = _synthetic_receipt()
    mod.enforce_receipt(receipt)

    bad = copy.deepcopy(receipt)
    bad["resolution_ladder"][-1]["lambda_turn"]["total"]["sampled_max"] = 2.0e-5
    with pytest.raises(AssertionError, match="lambda_turn/total"):
        mod.enforce_receipt(bad)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="truth boundary"):
        mod.enforce_receipt(promoted)

    moved = copy.deepcopy(receipt)
    moved["upstream"]["head"] = "0" * 40
    with pytest.raises(AssertionError, match="A2 #1004"):
        mod.enforce_receipt(moved)

    unstable = copy.deepcopy(receipt)
    unstable["resolution_ladder"][-2]["inner"]["leading"]["weighted_rms"] = 1.0e-7
    unstable["resolution_ladder"][-1]["inner"]["leading"]["weighted_rms"] = 2.0e-7
    with pytest.raises(AssertionError, match="resolution instability"):
        mod.enforce_receipt(unstable)


def test_receipt_checksum_helper_is_deterministic() -> None:
    payload = {"b": 2, "a": {"x": 1}}
    assert mod._sha256(payload) == mod._sha256({"a": {"x": 1}, "b": 2})
