from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_divergence_floor_routing_checkpoint import (
    AGENT2_HEAD,
    AGENT4_ARTIFACT_ID,
    AGENT4_HEAD,
    CURRENT_A4_METRICS,
    FORMAL_GATES,
    FROZEN_GUARDS,
    PRIOR_A4_543_METRICS,
    build_checkpoint,
    validate_checkpoint,
)


def test_v39_closes_covariance_but_keeps_divergence_and_correction_closed() -> None:
    checkpoint = build_checkpoint()
    states = checkpoint["states"]

    assert states["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"] is True
    assert states["public_oscillatory_project_support_passed"] is True
    assert states["public_oscillatory_covariance_rank_passed"] is True
    assert states["candidate_independent_second_covariance_column_ready"] is True
    assert states["public_oscillatory_local_divergence_max_passed"] is True
    assert states["public_oscillatory_local_divergence_rms_passed"] is False
    assert states["public_oscillatory_divergence_refinement_passed"] is False
    assert states["public_oscillatory_local_divergence_passed"] is False
    assert states["public_oscillatory_preflight_passed"] is False
    assert states["oscillatory_ready"] is False
    assert states["correction_ingest_allowed"] is False
    assert states["correction_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False


def test_v39_records_exact_independent_audit_identity_and_frozen_guards() -> None:
    checkpoint = build_checkpoint()
    a4 = checkpoint["upstream"]["agent4_carrier_resolved_independent_audit"]
    handoff = checkpoint["typed_handoffs"]["correction"]["required_actual_identity"]

    assert a4["audited_agent2_head"] == AGENT2_HEAD
    assert a4["head"] == AGENT4_HEAD
    assert a4["artifact_id"] == AGENT4_ARTIFACT_ID
    assert a4["dedicated_workflow_conclusion"] == "success"
    assert a4["standard_workflow_conclusion"] == "success"
    assert a4["scientific_guards_changed_from_agent4_543"] is False
    assert handoff["agent2_head"] == AGENT2_HEAD
    assert handoff["agent4_head"] == AGENT4_HEAD
    assert handoff["agent4_artifact_id"] == AGENT4_ARTIFACT_ID
    assert checkpoint["formal_gates"] == FORMAL_GATES


def test_v39_matches_scientific_pass_fail_partition() -> None:
    rms = CURRENT_A4_METRICS["relative_divergence_rms_by_step"]["0.005"]
    vmax = CURRENT_A4_METRICS["finest_relative_divergence_max"]
    refs = CURRENT_A4_METRICS["relative_rms_refinement_ratios"]
    ranks = CURRENT_A4_METRICS["covariance_rank_ratio_by_phase_resolution"].values()
    drift = CURRENT_A4_METRICS["covariance_relative_drift_to_finest"]

    assert rms > FROZEN_GUARDS["finest_relative_divergence_rms"]
    assert vmax <= FROZEN_GUARDS["finest_relative_divergence_max"]
    assert min(refs) < FROZEN_GUARDS["minimum_divergence_refinement_ratio"]
    assert min(ranks) >= FROZEN_GUARDS["minimum_covariance_rank_ratio"]
    assert max(drift) <= FROZEN_GUARDS["maximum_covariance_resolution_drift"]
    assert CURRENT_A4_METRICS["axis_near_absolute_max"] == 0.0
    assert CURRENT_A4_METRICS["radial_exterior_absolute_max"] == 0.0
    assert CURRENT_A4_METRICS["project_axial_exterior_absolute_max"] == 0.0


def test_v39_records_large_component_improvement_without_full_ns_claim() -> None:
    checkpoint = build_checkpoint()
    comparison = checkpoint["baseline_vs_kokuno"]
    component = comparison["component_preflight_only"]

    assert component["finest_rms_improvement_factor"] > 2.0e3
    assert component["finest_max_improvement_factor"] > 2.0e3
    assert PRIOR_A4_543_METRICS["minimum_covariance_rank_ratio"] < 0.02
    assert component["current_agent4_553"]["minimum_covariance_rank_ratio"] > 0.09
    assert component["is_full_ns_comparison"] is False
    assert comparison["kokuno_current_comparable_full_domain_receipt"] is None


def test_v39_hash_round_trip_is_deterministic() -> None:
    first = build_checkpoint()
    second = json.loads(json.dumps(build_checkpoint(), sort_keys=True))
    assert first == second
    validate_checkpoint(second)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("states", "oscillatory_ready"), True),
        (("states", "correction_ingest_allowed"), True),
        (("states", "pde_validated"), True),
        (("truth_boundary", "paper_exact"), True),
        (("truth_boundary", "thresholds_relaxed_after_result"), True),
    ],
)
def test_v39_tampering_fails_closed(path: tuple[str, str], value: bool) -> None:
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered[path[0]][path[1]] = value
    # Preserve the old digest deliberately: identity and science-state tampering
    # must never silently turn into a valid checkpoint.
    with pytest.raises(ValueError):
        validate_checkpoint(tampered)
