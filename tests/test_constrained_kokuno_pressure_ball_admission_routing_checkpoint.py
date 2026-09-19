import copy

import pytest

from openai_ns_reconstruction.kokuno_pressure_ball_admission_routing_checkpoint import (
    AGENT4_ARTIFACT_DIGEST,
    AGENT4_ARTIFACT_ID,
    AGENT4_AUDIT_HEAD,
    AGENT4_AUDIT_PR,
    AGENT4_DEDICATED_RUN,
    FORMAL_GATES,
    build_checkpoint,
    validate_checkpoint,
)


def test_v49_admits_only_independently_audited_pressure_ball_operator_algebra():
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)

    a4 = checkpoint["upstream"]["agent4_independent_pressure_ball_audit"]
    assert a4["pr"] == AGENT4_AUDIT_PR
    assert a4["head"] == AGENT4_AUDIT_HEAD
    assert a4["dedicated_workflow_run"] == AGENT4_DEDICATED_RUN
    assert a4["dedicated_workflow_conclusion"] == "success"
    assert a4["artifact_id"] == AGENT4_ARTIFACT_ID
    assert a4["artifact_zip_digest"] == AGENT4_ARTIFACT_DIGEST
    assert a4["fresh_replay_failed_guards"] == []
    assert a4["fresh_replay_passed"] is True

    metrics = checkpoint["independent_pressure_ball_metrics"]
    assert metrics["case_count"] == 37
    assert metrics["minimum_public_to_independent_upper_ratio"] >= 1.0
    assert metrics["public_underbounds"] == 0
    assert metrics["eta_factor_underbounds"] == 0
    assert metrics["downward_mutation_detected"] == 296
    assert metrics["downward_mutation_total"] == 296

    states = checkpoint["states"]
    assert states["pressure_ball_operator_algebra_independent_preflight_passed"] is True
    assert states["pressure_ball_operator_algebra_independently_admitted"] is True
    assert states["conditional_pressure_ball_operator_algebra_ready"] is True
    assert states["leading_ready"] is False
    assert states["source_axis_domain_independently_admitted"] is False
    assert states["source_pressure_radius_one_ball_bounds_ready"] is False
    assert states["source_R1_R2_MK_ready"] is False
    assert states["global_matched_pressure_ready"] is False
    assert states["correction_ready"] is False
    assert states["velocity_export_ready"] is False
    assert states["pde_validated"] is False
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v49_keeps_newer_source_axis_domain_and_correction_ledger_unadmitted():
    checkpoint = build_checkpoint()
    siblings = checkpoint["upstream"]["newer_sibling_work"]
    a1 = siblings["agent1_source_axis_domain"]
    assert a1["independent_agent4_audit_required"] is True
    assert a1["new_sigma_selected_center_pressure_preflight_required"] is True
    assert a1["admitted_by_this_checkpoint"] is False
    assert siblings["agent2_public_bundle"]["admitted_by_this_checkpoint"] is False
    assert siblings["agent3_finite_cycle_ledger"]["admitted_by_this_checkpoint"] is False
    assert checkpoint["truth_boundary"]["old_selected_center_pressure_pass_may_transfer_to_new_sigma"] is False


def test_v49_is_deterministic():
    a = build_checkpoint()
    b = build_checkpoint()
    assert a == b
    assert a["checkpoint_sha256"] == b["checkpoint_sha256"]


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("states", "leading_ready"), True),
        (("states", "source_pressure_radius_one_ball_bounds_ready"), True),
        (("states", "correction_ready"), True),
        (("states", "pde_validated"), True),
        (("truth_boundary", "old_selected_center_pressure_pass_may_transfer_to_new_sigma"), True),
        (("upstream", "newer_sibling_work", "agent1_source_axis_domain", "admitted_by_this_checkpoint"), True),
        (("upstream", "agent4_independent_pressure_ball_audit", "dedicated_workflow_conclusion"), "failure"),
        (("upstream", "agent4_independent_pressure_ball_audit", "artifact_zip_digest"), "sha256:mutated"),
    ],
)
def test_v49_fail_closed_mutations(path, value):
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    node = mutated
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    # Re-sign to prove validation checks semantics, not only the checksum.
    from openai_ns_reconstruction import kokuno_pressure_ball_admission_routing_checkpoint as routing

    mutated["checkpoint_sha256"] = routing._canonical_sha256(mutated)
    with pytest.raises(ValueError):
        validate_checkpoint(mutated)
