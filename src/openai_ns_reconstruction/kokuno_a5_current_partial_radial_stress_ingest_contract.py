"""Agent-5 integration contract for the current partial nonlinear radial-stress seam.

This module is integration/provenance glue only. It registers the newest current-
lineage correction-side artifact and its independent audit without reimplementing
Agent-3 mathematics or Agent-4 validation.

Registered scientific seam:
  A5 #973 partial leading+oscillatory/nonlinear-mean registration
  -> A3 #976 current partial nonlinear m=0 means routed through the inherited
     compact radial-stress inverse
  -> A4 #977 implementation-distinct serialized-receipt radial-stress audit.

Agent-2 #975 is recorded as the latest production composite differential sibling,
but #976/#977 still consume the authenticated #970/#960 surfaces. This contract
must not silently rewrite that lineage.

The radial stress remains a scoped partial-domain diagnostic through current X_h.
It is not a complete NS defect, is not an authorized correction target, and does
not establish a Cartesian correction velocity, matched pressure, restricted
forcing, same-protocol residual reduction, or PDE validation.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-partial-radial-stress-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-PARTIAL-RADIAL-STRESS-INGEST-094"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 973,
    "head": "f989bd99888e6263f9acb7008a932aa1db08d0aa",
    "branch": "codex/kokuno-a5-current-partial-composite-mean-ingest-093",
    "role": "current partial leading+oscillatory and nonlinear-mean registration through X_h",
}

AGENT1_LEADING = {
    "pr": 965,
    "head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "source_blob": "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c",
    "role": "current Cartesian leading velocity through current X_h",
    "velocity_beyond_xh_materialized": False,
}

AGENT2_COMPOSITE = {
    "pr": 970,
    "head": "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4",
    "source_blob": "532435705bd5db341e82371180da8c04c3dd6c14",
    "oscillatory_differentials_pr": 960,
    "oscillatory_differentials_head": "6d2fb1f701a34f783dca15a267ae2ce0734ba741",
    "oscillatory_differentials_source_blob": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "role": "authenticated partial u_lead + u_osc surface consumed by A3 #976",
}

AGENT2_LATEST_DIFFERENTIAL_SIBLING = {
    "pr": 975,
    "head": "b716dfe495763439ac51f6e866211337910017cb",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_partial_composite_batch_differentials.py",
    "source_blob": "c7772c481c3bebeffd6774ed29fbc37c46647ed8",
    "test_blob": "4cd3c3ff753326fdc676945d7c1b6cd9faa71ccb",
    "workflow_blob": "f0d48ca204bc6b6ae45efb0ef36a17680c45589a",
    "role": "production batch velocity/Jacobian/divergence/vorticity surface",
    "consumed_by_agent3_976": False,
    "replaces_independent_agent4_validation": False,
}

AGENT3_RADIAL_STRESS = {
    "pr": 976,
    "head": "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6",
    "branch": "codex/kokuno-a3-current-partial-nonlinear-radial-stress-103",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_partial_nonlinear_radial_stress.py",
    "source_blob": "9afa50c573186b86cd86fa7ffe44a5986a808b6d",
    "test_path": "tests/test_constrained_kokuno_current_partial_nonlinear_radial_stress.py",
    "test_blob": "f5fea08095ff14ac7cecc92ad0387d76b09a0c46",
    "workflow_path": ".github/workflows/kokuno-agent3-current-partial-nonlinear-radial-stress.yml",
    "workflow_blob": "9ea77ba1bd0e3d51e45044e1f1ec3a88109a7eb3",
    "parent_mean_pr": 971,
    "parent_mean_head": "bdb02ad487a175ab62748687c3172127d09aeb9c",
    "parent_mean_source_blob": "9c62ae3c3b4b3ca2d9082f822d57feef700ea6cf",
    "inherited_radial_operator_blob": "19f8cf9fa909fef7be2cf53c84df269e696b647c",
    "public_entry_point": "materialize_current_partial_nonlinear_radial_stress(backend, geometry)",
    "channels": {"theta": 2, "axial": 1},
    "radial_component_inverted": False,
    "mean_piece_closure_gate": 1.0e-12,
    "stress_piece_relative_closure_gate": 2.0e-11,
    "partial_domain_through_current_xh_only": True,
    "complete_ns_defect": False,
    "authorized_as_correction_target": False,
    "candidate_residual_evidence": False,
}

AGENT4_RADIAL_STRESS_AUDIT = {
    "pr": 977,
    "head": "2b015feccb5d55735ced8763fcb05a563a4cf12e",
    "branch": "codex/kokuno-a4-current-partial-radial-stress-audit-096",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_partial_radial_stress_independent_audit.py",
    "source_blob": "016215939f3cda805ae2dbbe421787ad0dfbee16",
    "test_blob": "815e5707a4acc4bac5692f478aebc8c5cfcc87c4",
    "workflow_blob": "066e2a38d83312f9f07feb20d50081d7b0a56b1c",
    "audited_agent3_pr": 976,
    "audited_agent3_head": "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6",
    "reference_path": "piecewise cubic interpolation plus order-8 Gauss-Legendre cell integration",
    "uses_production_radial_inverse_as_reference": False,
    "seed": 9173681,
    "time": 0.39,
    "axial_z": -0.08,
    "radial_interval": [0.005, 0.44],
    "nested_radial_counts": [43, 85, 169],
    "bump_center": 0.30,
    "bump_halfwidth": 0.10,
    "fine_stress_relative_rms_gate": 5.0e-2,
    "sampled_relative_max_gate": 1.5e-1,
    "weighted_moment_relative_error_gate": 3.0e-2,
    "positive_radius_axis_near_normalized_error_gate": 1.5e-1,
    "normalized_outer_edge_stress_gate": 1.0e-8,
    "exact_cartesian_axis_evidence": False,
    "candidate_residual_evidence": False,
    "scientific_admission": False,
}

FROZEN_SCIENCE = {
    "viscosity": 0.01,
    "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
    "time_interval": [0.25, 0.75],
    "forcing_family": "restricted two-parameter curl forcing",
    "residual_defined_free_forcing_forbidden": True,
}

FINAL_GATE = {
    "normalized_momentum_sampled_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_sampled_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
    "canonical_volume_quadrature_ladder": [24, 48, 96],
}

ST006_BASELINE = {
    "same_protocol_full_pde_baseline": True,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
}

READINESS = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}

TRUTH_BOUNDARY = {
    "current_cartesian_leading_velocity_materialized_through_xh": True,
    "current_partial_leading_plus_oscillatory_velocity_materialized_through_xh": True,
    "current_partial_nonlinear_m0_means_materialized": True,
    "current_partial_nonlinear_radial_inverse_performed": True,
    "current_partial_quadratic_radial_stress_materialized": True,
    "independent_a4_radial_stress_audit_registered": True,
    "independent_a4_radial_stress_audit_admitted": False,
    "velocity_beyond_xh_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "scoped_current_partial_radial_stress_authorized_as_correction_target": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "canonical_whole_domain_divergence_l2_assessed": False,
    "heldout_normalized_ns_residual_assessed": False,
    "same_protocol_st006_comparison_available_now": False,
    "residual_reduction_claimed": False,
    "scientific_admission": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "stage": "current partial nonlinear mean -> scoped radial inverse -> independent radial-stress audit registration",
    "input": "A5 #973 + A3 #976 + A4 #977; A2 #975 retained as latest differential sibling",
    "new_output": "checksum-bound registration of current partial nonlinear radial stress and its implementation-distinct public-receipt audit",
    "not_output": "complete NS defect, authorized correction target, correction velocity, global velocity, matched pressure, restricted forcing, full residual, or PDE validation",
    "next_shortest_blocker": "global completion beyond X_h plus matched pressure/restricted forcing are required before the scoped stress can be promoted into a complete-defect correction cycle",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA")


def build_contract(exact_head: str) -> dict[str, Any]:
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_leading": copy.deepcopy(AGENT1_LEADING),
        "agent2_composite": copy.deepcopy(AGENT2_COMPOSITE),
        "agent2_latest_differential_sibling": copy.deepcopy(AGENT2_LATEST_DIFFERENTIAL_SIBLING),
        "agent3_radial_stress": copy.deepcopy(AGENT3_RADIAL_STRESS),
        "agent4_radial_stress_audit": copy.deepcopy(AGENT4_RADIAL_STRESS_AUDIT),
        "frozen_science": copy.deepcopy(FROZEN_SCIENCE),
        "final_gate": copy.deepcopy(FINAL_GATE),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
        "registration_only": True,
    }
    payload["contract_sha256"] = _sha256(payload)
    return payload


def validate_contract(payload: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    observed = copy.deepcopy(dict(payload))
    supplied_sha = observed.pop("contract_sha256", None)
    if supplied_sha != _sha256(observed):
        errors.append("contract_sha256_mismatch")

    expected = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent1_leading": AGENT1_LEADING,
        "agent2_composite": AGENT2_COMPOSITE,
        "agent2_latest_differential_sibling": AGENT2_LATEST_DIFFERENTIAL_SIBLING,
        "agent3_radial_stress": AGENT3_RADIAL_STRESS,
        "agent4_radial_stress_audit": AGENT4_RADIAL_STRESS_AUDIT,
        "frozen_science": FROZEN_SCIENCE,
        "final_gate": FINAL_GATE,
        "st006_baseline": ST006_BASELINE,
        "readiness": READINESS,
        "truth_boundary": TRUTH_BOUNDARY,
        "pipeline_position": PIPELINE_POSITION,
        "registration_only": True,
    }
    for key, value in expected.items():
        if observed.get(key) != value:
            errors.append(f"{key}_drift")

    try:
        _require_hex40(observed.get("agent5_exact_head"), "agent5_exact_head")
    except ValueError:
        errors.append("agent5_exact_head_invalid")

    a2 = observed.get("agent2_latest_differential_sibling", {})
    a3 = observed.get("agent3_radial_stress", {})
    a4 = observed.get("agent4_radial_stress_audit", {})
    truth = observed.get("truth_boundary", {})
    readiness = observed.get("readiness", {})
    gate = observed.get("final_gate", {})
    science = observed.get("frozen_science", {})

    if a3.get("head") != a4.get("audited_agent3_head"):
        errors.append("a3_a4_lineage_mismatch")
    if a2.get("consumed_by_agent3_976") is not False:
        errors.append("a2_975_lineage_laundered")
    if a4.get("uses_production_radial_inverse_as_reference") is not False:
        errors.append("a4_independence_laundered")
    if a4.get("exact_cartesian_axis_evidence") is not False:
        errors.append("positive_radius_scope_laundered")

    required_true = (
        "current_cartesian_leading_velocity_materialized_through_xh",
        "current_partial_leading_plus_oscillatory_velocity_materialized_through_xh",
        "current_partial_nonlinear_m0_means_materialized",
        "current_partial_nonlinear_radial_inverse_performed",
        "current_partial_quadratic_radial_stress_materialized",
        "independent_a4_radial_stress_audit_registered",
    )
    required_false = (
        "independent_a4_radial_stress_audit_admitted",
        "velocity_beyond_xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "scoped_current_partial_radial_stress_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "canonical_whole_domain_divergence_l2_assessed",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in required_true:
        if truth.get(key) is not True:
            errors.append(f"truth_true_required:{key}")
    for key in required_false:
        if truth.get(key) is not False:
            errors.append(f"truth_false_required:{key}")

    if readiness != READINESS:
        errors.append("readiness_promoted")
    if gate != FINAL_GATE:
        errors.append("final_gate_drift")
    if science.get("residual_defined_free_forcing_forbidden") is not True:
        errors.append("free_forcing_firewall_removed")
    if a3.get("candidate_residual_evidence") is not False:
        errors.append("a3_stress_laundered_as_residual")
    if a4.get("candidate_residual_evidence") is not False:
        errors.append("a4_stress_audit_laundered_as_residual")
    return errors


def enforce_contract(payload: Mapping[str, Any]) -> None:
    errors = validate_contract(payload)
    if errors:
        raise AssertionError(";".join(errors))
