from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_rejection_routing_checkpoint import (
    AGENT1,
    AGENT4,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_v37_records_public_provider_but_rejects_oscillatory_readiness() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]
    assert states["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"]
    assert states["public_oscillatory_independent_audit_assessed"]
    assert not states["public_oscillatory_local_divergence_passed"]
    assert not states["public_oscillatory_covariance_rank_passed"]
    assert not states["public_oscillatory_project_support_passed"]
    assert not states["oscillatory_ready"]
    assert not states["correction_ingest_allowed"]
    assert not states["pde_validated"]


def test_agent4_rejection_metrics_and_frozen_guards_are_exactly_routed() -> None:
    a4 = build_checkpoint()["upstream"]["agent4_public_provider_audit"]
    metrics = a4["metrics"]
    guards = a4["guards"]
    assert a4 == AGENT4
    assert a4["dedicated_run"] == 35413803718
    assert a4["dedicated_status"] == "success"
    assert a4["artifact_id"] == 10575081799
    assert not a4["public_oscillatory_preflight_passed"]
    assert metrics["finest_relative_divergence_rms"] == pytest.approx(0.10460795155383215)
    assert metrics["finest_relative_divergence_max"] == pytest.approx(0.09811443995058594)
    assert metrics["minimum_covariance_rank_ratio"] == pytest.approx(0.01874819346354583)
    assert metrics["project_axial_exterior_absolute_max"] == pytest.approx(769255.3191046526)
    assert metrics["finest_relative_divergence_rms"] > guards["finest_relative_divergence_rms"]
    assert metrics["finest_relative_divergence_max"] > guards["finest_relative_divergence_max"]
    assert min(metrics["divergence_refinement_ratios"]) < guards["minimum_divergence_refinement_ratio"]
    assert metrics["minimum_covariance_rank_ratio"] < guards["minimum_covariance_rank_ratio"]
    assert metrics["project_axial_exterior_absolute_max"] > guards["project_axial_exterior_absolute_max"]


def test_v37_preserves_passed_subchecks_without_laundering_failure() -> None:
    a4 = build_checkpoint()["upstream"]["agent4_public_provider_audit"]
    metrics = a4["metrics"]
    guards = a4["guards"]
    assert a4["radial_axis_support_passed"]
    assert metrics["axis_near_absolute_max"] <= guards["axis_near_absolute_max"]
    assert metrics["radial_exterior_absolute_max"] <= guards["radial_exterior_absolute_max"]
    assert metrics["maximum_covariance_resolution_drift"] <= guards["maximum_covariance_resolution_drift"]
    assert metrics["heldout_velocity_rms"] >= guards["minimum_nontrivial_velocity_rms"]
    assert min(
        metrics["h_minus_10pct_relative_change"],
        metrics["h_plus_10pct_relative_change"],
    ) >= guards["minimum_parameter_perturbation_relative_change"]
    assert metrics["divergence_mutation_relative_rms"] >= guards["minimum_divergence_mutation_relative_rms"]
    assert metrics["duplicated_covariance_rank_ratio"] <= guards["maximum_duplicated_covariance_rank_ratio"]
    assert not a4["formal_full_domain_pde_gate_assessed"]
    assert not a4["heldout_ns_residual_assessed"]
    assert not a4["pde_validated"]


def test_v37_updates_leading_truth_boundary_to_source_choice_order_guard() -> None:
    a1 = build_checkpoint()["upstream"]["agent1_latest"]
    assert a1 == AGENT1
    assert a1["head"] == "61b24c1ced2b9246df26e1ee4148392d2ff521d9"
    assert a1["dedicated_status"] == "success"
    assert a1["standard_status"] == "success"
    assert a1["source_choice_order_guard_executable"]
    assert a1["fresh_coupled_C_pointwise_screen_passed"]
    assert not a1["fresh_coupled_C_pass_is_source_ordered_certificate"]
    assert not a1["source_B0_dependencies_machine_bound"]
    assert not a1["source_T_sh_lower_bound_verified"]
    assert not a1["selected_pa16_handoff_allowed"]


def test_v37_fail_closed_on_threshold_or_state_laundering() -> None:
    payload = copy.deepcopy(build_checkpoint())
    payload["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="sha256 mismatch"):
        validate_checkpoint(payload)

    payload = copy.deepcopy(build_checkpoint())
    payload["states"]["oscillatory_ready"] = True
    with pytest.raises(ValueError, match="routing states changed"):
        validate_checkpoint(_rehash(payload))

    payload = copy.deepcopy(build_checkpoint())
    payload["upstream"]["agent4_public_provider_audit"]["public_oscillatory_preflight_passed"] = True
    with pytest.raises(ValueError, match="upstream receipt changed"):
        validate_checkpoint(_rehash(payload))


def test_v37_roundtrip_is_deterministic(tmp_path) -> None:
    expected = build_checkpoint()
    path = tmp_path / "v37.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert written == expected
    assert loaded == expected
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)
