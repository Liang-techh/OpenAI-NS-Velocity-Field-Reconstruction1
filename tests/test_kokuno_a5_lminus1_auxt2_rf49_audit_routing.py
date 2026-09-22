from __future__ import annotations

import copy
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/openai_ns_reconstruction/kokuno_a5_lminus1_auxt2_rf49_audit_routing.py"
SPEC = importlib.util.spec_from_file_location("kokuno_a5_lminus1_auxt2_rf49_audit_routing", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_artifact_registers_fresh_frontier_without_truth_promotion() -> None:
    artifact = mod.build_artifact()
    routes = artifact["upstream_routes"]

    assert artifact["parent_a5"] == {
        "pr": 1192,
        "exact_head": "18f9df62154dd5230eede98eb3634a74560b9e97",
    }
    assert routes["agent1"]["pr"] == 1196
    assert routes["agent1"]["source_l_minus1_hold_materialized"] is True
    assert routes["agent1"]["source_l_minus1_to_minus_h_transition_materialized"] is False
    assert routes["agent1"]["outer_global_leading_velocity_materialized"] is False

    assert routes["agent2"]["pr"] == 1198
    assert routes["agent2"]["source_blob"] == "96168ac6583ad5e71aa044bd558c6feeadf051f4"
    assert routes["agent2"]["current_i4_raw_auxiliary_t2_provider_materialized"] is True
    assert routes["agent2"]["provider_checksum_pinned_by_agent3"] is False
    assert routes["agent2"]["rf30_repository_candidate_state_authorized"] is False

    assert routes["agent3"]["pr"] == 1197
    assert routes["agent3"]["rf44_rf49_postupdate_handoff_executable"] is True
    assert routes["agent3"]["pinned_same_identity_repository_provider_blob"] is None
    assert routes["agent3"]["cartesian_correction_velocity_materialized"] is False

    assert routes["agent4"]["pr"] == 1199
    assert routes["agent4"]["audits_agent1_pr"] == 1196
    assert routes["agent4"]["audits_agent1_exact_head"] == routes["agent1"]["exact_head"]
    assert routes["agent4"]["scoped_gate_passed"] is None
    assert routes["agent4"]["pde_validated"] is False

    assert artifact["latest_self_contained_candidate"]["agent2_pr"] == 1117
    assert artifact["latest_self_contained_candidate"]["agent1_pr"] == 1107
    assert artifact["latest_self_contained_candidate"]["stage"] == "xi=11"
    assert artifact["core_states"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert not any(artifact["truth_boundary"].values())


def test_closed_and_remaining_blockers_are_distinguished() -> None:
    blockers = mod.build_artifact()["remaining_blockers"]
    assert blockers["a1_l_minus1_hold"] is False
    assert blockers["raw_same_identity_auxiliary_t2_provider_missing"] is False
    assert blockers["a1_minus1_to_minus_h_transition"] is True
    assert blockers["agent3_checksum_pin_of_exact_agent2_1198_provider_blob"] is True
    assert blockers["scientific_candidate_specific_rf30_to_rf49_correction"] is True
    assert blockers["cartesian_correction_velocity"] is True
    assert blockers["complete_identity_bound_ns_defect"] is True
    assert blockers["agent4_heldout_canonical_complete_ns_gate"] is True


def test_validation_contract_is_frozen() -> None:
    protocol = mod.build_artifact()["fixed_validation_protocol"]
    assert protocol["st006_momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert protocol["st006_momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert protocol["normalized_momentum_sampled_max_gate"] == 1.0e-3
    assert protocol["normalized_momentum_volume_l2_gate"] == 1.0e-3
    assert protocol["normalized_divergence_sampled_max_gate"] == 1.0e-5
    assert protocol["normalized_divergence_volume_l2_gate"] == 1.0e-5
    assert protocol["canonical_quadrature"] == [24, 48, 96]
    assert protocol["residual_defined_free_forcing_forbidden"] is True
    assert protocol["posthoc_threshold_relaxation_forbidden"] is True


def test_provider_availability_cannot_be_mutated_into_agent3_authorization() -> None:
    artifact = mod.build_artifact()
    mutated = copy.deepcopy(artifact)
    mutated["upstream_routes"]["agent2"]["provider_checksum_pinned_by_agent3"] = True
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)

    mutated = copy.deepcopy(artifact)
    mutated["upstream_routes"]["agent3"]["pinned_same_identity_repository_provider_blob"] = mod.AGENT2_SOURCE_BLOB
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)

    mutated = copy.deepcopy(artifact)
    mutated["upstream_routes"]["agent3"]["repository_candidate_rf44_rf49_evidence"] = True
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)


def test_scoped_a4_evidence_cannot_promote_full_science() -> None:
    artifact = mod.build_artifact()
    for key, value in (
        ("scoped_gate_passed", True),
        ("scientifically_admitted", True),
        ("pde_validated", True),
    ):
        mutated = copy.deepcopy(artifact)
        mutated["upstream_routes"]["agent4"][key] = value
        with pytest.raises(mod.KokunoA5RoutingError):
            mod.validate_artifact(mutated)


def test_latest_self_contained_candidate_cannot_be_silently_relabelled() -> None:
    artifact = mod.build_artifact()
    mutated = copy.deepcopy(artifact)
    mutated["latest_self_contained_candidate"]["superseded_by_agent1_1196_as_unified_candidate"] = True
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)


def test_threshold_and_forcing_firewalls_fail_closed() -> None:
    artifact = mod.build_artifact()

    mutated = copy.deepcopy(artifact)
    mutated["fixed_validation_protocol"]["normalized_momentum_sampled_max_gate"] = 2.0e-3
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)

    mutated = copy.deepcopy(artifact)
    mutated["fixed_validation_protocol"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(mod.KokunoA5RoutingError):
        mod.validate_artifact(mutated)


def test_artifact_write_is_deterministic(tmp_path: Path) -> None:
    first = mod.write_artifact(tmp_path / "first.json").read_bytes()
    second = mod.write_artifact(tmp_path / "second.json").read_bytes()
    assert first == second
