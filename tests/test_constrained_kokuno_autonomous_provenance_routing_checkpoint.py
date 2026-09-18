from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_autonomous_provenance_routing_checkpoint import (
    AUTONOMOUS_REALIZATION_POLICY,
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    validate_component_provenance,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_checkpoint_opens_only_labelled_autonomous_candidate_route() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    assert states["displayed_source_phase_frame_independently_audited"]
    assert states["autonomous_signed_rectangle_geometry_ready"]
    assert states["formal_mean_factor_dependency_pinned"]
    assert not states["source_exact_numeric_realization_required_for_candidate"]
    assert states["candidate_autonomous_realization_route_open"]
    assert states["candidate_autonomous_mean_profile_allowed"]

    for name in (
        "candidate_autonomous_mean_profile_selected",
        "candidate_state_requestedStress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "actual_source_numeric_bump_profile_recovered",
        "actual_source_missing_weight_materialized",
        "actual_source_requested_stress_materialized",
        "actual_source_finite_head_mean_debt_materialized",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]


def test_autonomous_provenance_cannot_be_laundered_into_source_truth() -> None:
    validate_component_provenance(
        {
            "evidence_class": "repository_autonomous_source_compatible",
            "source_exact": False,
            "recovered_source_numeric": False,
            "discharges_opaque_formal_numeric_witness": False,
        }
    )

    with pytest.raises(ValueError, match="promoted to recovered source truth"):
        validate_component_provenance(
            {
                "evidence_class": "repository_autonomous_source_compatible",
                "source_exact": True,
            }
        )
    with pytest.raises(ValueError, match="opaque formal numeric witness"):
        validate_component_provenance(
            {
                "evidence_class": "repository_autonomous_source_compatible",
                "discharges_opaque_formal_numeric_witness": True,
            }
        )
    with pytest.raises(ValueError, match="unknown evidence_class"):
        validate_component_provenance({"evidence_class": "paper_exact_by_assertion"})


def test_checkpoint_records_latest_agent2_and_agent3_green_receipts() -> None:
    payload = build_checkpoint()
    a2 = payload["upstream"]["agent2_sibling"]
    a3 = payload["upstream"]["agent3_ancestry"]

    assert a2["head"] == "58b95eafb3429685f5d7e83c22fa7b0a8f9eeba8"
    assert a2["dedicated_run"] == 35402140447
    assert a2["standard_run"] == 35402140319
    assert a2["dedicated_status"] == "success"
    assert a2["standard_status"] == "success"
    assert a2["evidence_class"] == "repository_autonomous_source_compatible"
    assert a2["autonomous_signed_rectangle_geometry_ready"]
    assert not a2["source_rectangle_centers_recovered"]
    assert not a2["source_rectangle_radius_r0_recovered"]

    assert a3["head"] == "d4ce470abe655378746e8e3131e2891280476619"
    assert a3["dedicated_run"] == 35402618032
    assert a3["standard_run"] == 35402618004
    assert a3["dedicated_status"] == "success"
    assert a3["standard_status"] == "success"
    assert a3["source_algebra_pinned"]
    assert not a3["numeric_bump_profile_export_present"]
    assert not a3["finite_head_mean_debt_materialized"]
    assert not a3["arbitrary_python_bump_may_be_promoted_to_actual"]
    assert a3["consumed_in_executable_ancestry"]


def test_checkpoint_keeps_latest_agent1_pending_separate_from_last_green() -> None:
    payload = build_checkpoint()
    latest = payload["upstream"]["agent1_latest_sibling"]
    stable = payload["upstream"]["agent1_last_fully_green_sibling"]

    assert latest["pr"] == 500
    assert latest["head"] == "c360b7a130e007f160bb2680415cc8f020db2e04"
    assert latest["standard_run"] == 35403080310
    assert latest["standard_status_at_freeze"] == "in_progress"
    assert not latest["global_pressure_matched"]
    assert stable["pr"] == 491
    assert stable["standard_status"] == "success"
    assert stable["PA16_input_tuple_formable"]
    assert not payload["states"]["leading_ready"]


def test_checkpoint_preserves_agent4_independence_and_formal_gates() -> None:
    payload = build_checkpoint()
    a4 = payload["upstream"]["agent4_sibling"]

    assert a4["head"] == "31f7510545488818b3bcf979c5cac3cb5d826a49"
    assert a4["dedicated_run"] == 35398848549
    assert a4["standard_run"] == 35398848498
    assert a4["dedicated_status"] == "success"
    assert a4["standard_status"] == "success"
    assert a4["evidence_class"] == "independent_validation"
    assert a4["local_structural_preflight_passed"]
    assert a4["fd_finest_relative_rms"] == pytest.approx(4.3321017055747775e-12)
    assert not a4["formal_full_domain_pde_gate_assessed"]
    assert not a4["pde_validated"]

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_checkpoint_fails_closed_on_source_or_candidate_overpromotion() -> None:
    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["candidate_autonomous_mean_profile_selected"] = True
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["actual_source_numeric_bump_profile_recovered"] = True
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent3_ancestry"]["numeric_bump_profile_export_present"] = True
    with pytest.raises(ValueError, match="opaque formal mean factor"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent2_sibling"]["source_rectangle_centers_recovered"] = True
    with pytest.raises(ValueError, match="promoted to recovered source truth"):
        validate_checkpoint(_rehash(payload))


def test_policy_and_checkpoint_roundtrip_are_deterministic(tmp_path) -> None:
    payload = build_checkpoint()
    assert payload["autonomous_realization_policy"] == AUTONOMOUS_REALIZATION_POLICY
    assert payload["candidate_artifact_contract"]["instantiated"] is False

    path = tmp_path / "checkpoint.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == written
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
