"""CR002 audit for current Kokuno partial radial-stress channel completeness.

The current A3/A5 lane materializes compact radial stresses only for the
tangential (e=2) and axial (e=1) cylindrical channels. The radial cylindrical
nonlinear mean is recorded but is not inverted by that operator. This module
fails closed if the two-channel result is promoted to a full three-component
stress, a complete NS defect, an authorized correction target, or a Cartesian
correction velocity.

This is representation/provenance governance only. It changes no candidate,
pressure, forcing, validation data, or scientific threshold.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-partial-radial-stress-channel-completeness-v1"
TASK_ID = "CR002-KOKUNO-PARTIAL-RADIAL-STRESS-CHANNEL-COMPLETENESS-099"
PARENT_A5_HEAD = "cb67fed42c459d9d753ad17cbc339744e51a9b84"
PARENT_A5_SOURCE_BLOB = "faf868c4b5b536cfcb46494c494697acda8e61b1"
A3_HEAD = "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6"
A3_SOURCE_BLOB = "9afa50c573186b86cd86fa7ffe44a5986a808b6d"
A4_HEAD = "2b015feccb5d55735ced8763fcb05a563a4cf12e"
A4_SOURCE_BLOB = "016215939f3cda805ae2dbbe421787ad0dfbee16"
CR001_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

CONFIG_PATH = (
    Path(__file__).resolve().parents[2]
    / "configs"
    / "kokuno_partial_radial_stress_channel_completeness_scope.json"
)

_REQUIRED_PENDING = {
    "materialize_or_mathematically_eliminate_the_missing_radial_cylindrical_stress_channel",
    "bind_all_required_correction_channels_to_one_cylindrical_and_cartesian_component_convention",
    "derive_the_correction_target_from_a_complete_ns_defect_on_the_same_global_candidate",
    "complete_velocity_beyond_current_xh_and_global_support_join",
    "materialize_matched_pressure_and_preregistered_restricted_forcing",
    "run_independent_whole_domain_heldout_validation_before_pde_promotion",
}


def load_scope(path: str | Path | None = None) -> dict[str, Any]:
    source = CONFIG_PATH if path is None else Path(path)
    return json.loads(source.read_text(encoding="utf-8"))


def validate_scope(scope: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []

    if scope.get("schema") != SCHEMA:
        errors.append("schema_drift")
    if scope.get("task_id") != TASK_ID:
        errors.append("task_id_drift")

    parent = scope.get("parent_a5", {})
    if parent.get("pr") != 978 or parent.get("head") != PARENT_A5_HEAD:
        errors.append("parent_a5_identity_drift")
    if parent.get("source_blob") != PARENT_A5_SOURCE_BLOB:
        errors.append("parent_a5_source_blob_drift")

    a3 = scope.get("upstream_agent3", {})
    if a3.get("pr") != 976 or a3.get("head") != A3_HEAD:
        errors.append("a3_identity_drift")
    if a3.get("source_blob") != A3_SOURCE_BLOB:
        errors.append("a3_source_blob_drift")
    if a3.get("materialized_stress_channels") != {
        "theta": {"exponent": 2},
        "axial": {"exponent": 1},
    }:
        errors.append("a3_materialized_channel_set_drift")
    if a3.get("radial_cylindrical_mean_recorded") is not True:
        errors.append("radial_mean_record_missing")
    if a3.get("radial_cylindrical_stress_inverted") is not False:
        errors.append("radial_stress_inversion_laundered")
    if a3.get("partial_domain_through_current_xh_only") is not True:
        errors.append("a3_partial_domain_scope_laundered")
    if a3.get("candidate_residual_evidence") is not False:
        errors.append("a3_candidate_residual_evidence_laundered")

    a4 = scope.get("upstream_agent4", {})
    if a4.get("pr") != 977 or a4.get("head") != A4_HEAD:
        errors.append("a4_identity_drift")
    if a4.get("source_blob") != A4_SOURCE_BLOB:
        errors.append("a4_source_blob_drift")
    if a4.get("audited_stress_channels") != ["theta", "axial"]:
        errors.append("a4_audited_channel_set_drift")
    if a4.get("radial_cylindrical_stress_independently_audited") is not False:
        errors.append("a4_radial_stress_audit_laundered")
    if a4.get("exact_cartesian_axis_evidence") is not False:
        errors.append("a4_axis_scope_laundered")
    if a4.get("candidate_residual_evidence") is not False:
        errors.append("a4_candidate_residual_evidence_laundered")
    if a4.get("scientific_admission") is not False:
        errors.append("a4_scientific_admission_laundered")

    distinctions = scope.get("governed_distinctions", {})
    for key in (
        "theta_axial_two_channel_stress_is_full_three_component_cylindrical_stress",
        "theta_axial_a4_audit_is_full_three_component_stress_validation",
        "recorded_radial_mean_is_inverted_radial_stress",
        "scoped_two_channel_stress_is_complete_ns_defect",
        "scoped_two_channel_stress_is_authorized_correction_target",
        "scoped_two_channel_stress_is_cartesian_correction_velocity",
        "scoped_two_channel_stress_implies_pde_validation",
    ):
        if distinctions.get(key) is not False:
            errors.append(f"{key}_laundered")

    pending = set(scope.get("pending_requirements", []))
    missing_pending = sorted(_REQUIRED_PENDING - pending)
    if missing_pending:
        errors.append("pending_requirements_dropped:" + ",".join(missing_pending))

    readiness = scope.get("readiness", {})
    if readiness.get("kokuno_two_channel_radial_stress_materialized") is not True:
        errors.append("two_channel_materialization_truth_dropped")
    for key in (
        "kokuno_full_three_component_cylindrical_stress_materialized",
        "kokuno_complete_ns_defect_materialized",
        "kokuno_complete_correction_target_materialized",
        "kokuno_cartesian_correction_velocity_materialized",
        "kokuno_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        if readiness.get(key) is not False:
            errors.append(f"{key}_premature_promotion")

    eq45 = scope.get("canonical_eq45_delivery", {})
    if eq45.get("velocity_export_ready") is not True:
        errors.append("eq45_delivery_incorrectly_downgraded")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        if eq45.get(key) is not False:
            errors.append(f"eq45_{key}_incorrectly_promoted")

    cr = scope.get("cr001", {})
    expected_cr = {
        "constraints_blob": CR001_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "residual_defined_pointwise_free_force_allowed": False,
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "collapsed_candidate_allowed": False,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max_gate": 0.001,
        "momentum_l2_gate": 0.001,
        "divergence_max_gate": 1.0e-5,
        "divergence_l2_gate": 1.0e-5,
        "post_hoc_threshold_relaxation_allowed": False,
    }
    if cr != expected_cr:
        errors.append("cr001_contract_drift")

    provenance = scope.get("four_way_provenance", {})
    if set(provenance) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        errors.append("four_way_provenance_partition_drift")

    witness = scope.get("mechanics_witness", {})
    if witness.get("classification") != "autonomous_design":
        errors.append("mechanics_witness_source_laundered")
    if witness.get("not_kokuno_or_openai_data") is not True:
        errors.append("mechanics_witness_disclaimer_dropped")
    va = witness.get("vector_a_radial_theta_axial")
    vb = witness.get("vector_b_radial_theta_axial")
    if not (
        isinstance(va, list)
        and isinstance(vb, list)
        and len(va) == len(vb) == 3
        and va[1:] == vb[1:]
        and va[0] != vb[0]
        and witness.get("same_observed_theta_axial") is True
        and witness.get("full_vectors_equal") is False
    ):
        errors.append("mechanics_witness_invalid")

    return errors


def enforce_scope(scope: Mapping[str, Any]) -> None:
    errors = validate_scope(scope)
    if errors:
        raise ValueError("CR002 radial-stress channel-completeness audit failed: " + "; ".join(errors))


def audit_scope(path: str | Path | None = None) -> dict[str, Any]:
    scope = load_scope(path)
    enforce_scope(scope)
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "status": "pass",
        "parent_a5_head": PARENT_A5_HEAD,
        "materialized_stress_channels": ["theta", "axial"],
        "missing_radial_stress_channel": True,
        "full_three_component_stress_materialized": False,
        "complete_ns_defect_materialized": False,
        "cartesian_correction_velocity_materialized": False,
        "kokuno_velocity_export_ready": False,
        "canonical_eq45_velocity_export_ready": True,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_scope(), indent=2, sort_keys=True))
