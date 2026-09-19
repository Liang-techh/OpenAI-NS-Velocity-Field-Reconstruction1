from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_post_axial_oscillatory_routing_checkpoint import (
    AGENT1,
    AGENT2,
    AGENT3,
    AGENT4,
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_v38_closes_support_only_and_keeps_oscillatory_rejected() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    state = payload["states"]
    assert state["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"]
    assert state["public_oscillatory_independent_audit_assessed"]
    assert state["public_oscillatory_project_support_passed"]
    assert not state["public_oscillatory_local_divergence_passed"]
    assert not state["public_oscillatory_covariance_rank_passed"]
    assert not state["oscillatory_ready"]
    assert not state["correction_ingest_allowed"]
    assert not state["pde_validated"]


def test_v38_routes_latest_exact_heads_and_green_software_ci() -> None:
    upstream = build_checkpoint()["upstream"]
    assert upstream["agent1_latest"] == AGENT1
    assert upstream["agent2_support_repair"] == AGENT2
    assert upstream["agent3_admission"] == AGENT3
    assert upstream["agent4_frozen_rerun"] == AGENT4

    assert AGENT1["head"] == "7375fc961bbbc3ee469936104fb1eea3e742a979"
    assert AGENT1["l0"] == pytest.approx(66.40155058292201)
    assert AGENT1["dedicated_status"] == AGENT1["standard_status"] == "success"

    assert AGENT2["head"] == "985d1c3fc43081d2d87feb1f0c27643cc4ebe75b"
    assert AGENT2["vector_potential_axial_support_localization"]
    assert AGENT2["analytic_D_z_coefficient_derivative_supplied"]
    assert AGENT2["dedicated_status"] == AGENT2["standard_status"] == "success"

    assert AGENT3["head"] == "e59cdd518460983bea2073f59524d4082769d94b"
    assert AGENT3["correction_admission_fail_closed"]
    assert not AGENT3["actual_audited_candidate_identity_bound_in_admission"]

    assert AGENT4["head"] == "9498209376fd9dc5abf228248967958759e9f64a"
    assert AGENT4["actual_audited_agent2_head"] == AGENT2["head"]
    assert AGENT4["artifact_id"] == 10576031714
    assert AGENT4["dedicated_status"] == AGENT4["standard_status"] == "success"


def test_v38_freezes_support_pass_and_remaining_numeric_failures() -> None:
    a4 = build_checkpoint()["upstream"]["agent4_frozen_rerun"]
    metrics = a4["metrics"]
    guards = a4["guards"]

    assert a4["project_support_preflight_passed"]
    assert a4["project_axial_support_passed"]
    assert metrics["project_axial_exterior_absolute_max"] == 0.0
    assert metrics["project_axial_exterior_absolute_max"] <= guards["project_axial_exterior_absolute_max"]

    assert metrics["relative_divergence_rms_by_step"][-1] == pytest.approx(0.15167269397404395)
    assert metrics["finest_relative_divergence_max"] == pytest.approx(0.19755153174056717)
    assert metrics["divergence_refinement_ratios"] == pytest.approx([1.3747692035625707, 3.488563215878927])
    assert metrics["relative_divergence_rms_by_step"][-1] > guards["finest_relative_divergence_rms"]
    assert metrics["finest_relative_divergence_max"] > guards["finest_relative_divergence_max"]
    assert min(metrics["divergence_refinement_ratios"]) < guards["minimum_divergence_refinement_ratio"]

    assert metrics["minimum_covariance_rank_ratio"] == pytest.approx(0.01873858406014804)
    assert metrics["minimum_covariance_rank_ratio"] < guards["minimum_covariance_rank_ratio"]
    assert metrics["maximum_covariance_resolution_drift"] <= guards["maximum_covariance_resolution_drift"]
    assert not a4["public_oscillatory_preflight_passed"]


def test_v38_distinguishes_protocol_origin_from_actual_audited_candidate() -> None:
    payload = build_checkpoint()
    a3 = payload["upstream"]["agent3_admission"]
    a4 = payload["upstream"]["agent4_frozen_rerun"]
    correction = payload["typed_handoffs"]["correction"]

    assert a4["protocol_origin_parent_agent2_head"] == a3["embedded_protocol_origin_parent_agent2_head"]
    assert a4["actual_audited_agent2_head"] != a4["protocol_origin_parent_agent2_head"]
    assert correction["required_actual_identity"]["agent2_head"] == a4["actual_audited_agent2_head"]
    assert correction["required_actual_identity"]["agent4_head"] == a4["head"]
    assert correction["required_actual_identity"]["agent4_artifact_id"] == a4["artifact_id"]
    assert not a3["actual_audited_candidate_identity_bound_in_admission"]
    assert not a3["actual_agent4_artifact_identity_bound_in_admission"]
    assert not payload["states"]["correction_ingest_allowed"]


def test_v38_preserves_fixed_final_gate_and_st006_comparison_boundary() -> None:
    payload = build_checkpoint()
    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1e-3
    assert payload["formal_gates"]["held_out_normalized_momentum_l2"] == 1e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1e-5
    assert payload["formal_gates"]["held_out_divergence_l2"] == 1e-5
    baseline = payload["baseline_vs_kokuno"]
    assert baseline["st006"]["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["st006"]["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert baseline["kokuno_current_comparable_full_domain_receipt"] is None


def test_v38_fails_closed_on_gate_state_receipt_or_identity_laundering() -> None:
    payload = copy.deepcopy(build_checkpoint())
    payload["formal_gates"]["held_out_normalized_momentum_max"] = 2e-3
    with pytest.raises(ValueError, match="sha256 mismatch"):
        validate_checkpoint(payload)

    payload = copy.deepcopy(build_checkpoint())
    payload["states"]["oscillatory_ready"] = True
    with pytest.raises(ValueError, match="routing states changed"):
        validate_checkpoint(_rehash(payload))

    payload = copy.deepcopy(build_checkpoint())
    payload["upstream"]["agent4_frozen_rerun"]["public_oscillatory_preflight_passed"] = True
    with pytest.raises(ValueError, match="upstream receipt changed"):
        validate_checkpoint(_rehash(payload))

    payload = copy.deepcopy(build_checkpoint())
    payload["upstream"]["agent4_frozen_rerun"]["actual_audited_agent2_head"] = payload["upstream"]["agent4_frozen_rerun"]["protocol_origin_parent_agent2_head"]
    with pytest.raises(ValueError, match="upstream receipt changed"):
        validate_checkpoint(_rehash(payload))


def test_v38_roundtrip_is_deterministic(tmp_path) -> None:
    expected = build_checkpoint()
    path = tmp_path / "v38.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert written == expected
    assert loaded == expected
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)
