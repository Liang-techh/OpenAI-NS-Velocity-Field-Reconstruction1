"""Fail-closed CR002 audit for Kokuno final-bridge width/cone semantics.

This auditor governs a narrow representation/provenance boundary.  The current
candidate-side final bridge uses repository-autonomous logarithmic widths and
has an implementation-distinct numerical profile/derivative audit.  Neither
fact is evidence that those concrete widths satisfy the source's later
proof-dependent cone/global-smallness conditions, nor that hidden source width
values have been recovered.

No scientific field is evaluated or modified here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-final-bridge-width-cone-scope-v1"
PROVENANCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _load_json(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"{path} must contain a JSON object")
    return payload


def _assert_parent_contract(parent_source: str) -> None:
    required = (
        '"axial_shutdown_log_width": 0.02',
        '"angular_settle_log_width": 0.02',
        '"selected_kappa0_is_repository_autonomous": True',
        '"source_admitted_global_kappa0_smallness": False',
        '"final_bridge_cone_admissibility_independently_certified": False',
        '"actual_final_bridge_independent_a4_audit_admitted": False',
        '"source_coordinate_profiles_are_cartesian_velocity": False',
    )
    for snippet in required:
        _require(snippet in parent_source, f"parent A5 contract drifted: missing {snippet}")


def _assert_canonical_constraints(scope: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    lock = scope["cr001_lock"]
    _require(constraints.get("nu") == lock["nu"] == 0.01, "CR001 viscosity drift")
    domain = constraints["domain"]
    _require(domain["physical"] == lock["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drift")
    _require(lock["evaluation_box"] == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]], "locked evaluation box drift")
    _require(domain["support"] == lock["support"] == "r < 2 and abs(z) < 2", "support drift")
    _require(domain["time_interval"] == lock["time_interval"] == [0.25, 0.75], "time interval drift")

    forcing = constraints["forcing"]
    _require(forcing["mode"] == lock["forcing_mode"] == "restricted_two_parameter_family", "forcing-family drift")
    _require("No residual-dependent basis or pointwise free force" in forcing["restriction"], "canonical free-force prohibition missing")
    _require(lock["residual_defined_free_force_forbidden"] is True, "scope enabled residual-defined force")

    nontrivial = constraints["nontriviality"]
    _require(nontrivial["reference_energy"] == lock["reference_energy"] == 1.0, "reference energy drift")
    _require(nontrivial["reference_energy_abs_tolerance"] == lock["reference_energy_abs_tolerance"] == 0.001, "energy tolerance drift")
    _require("reject collapsed candidates" in nontrivial["enforcement"], "collapse rejection missing")

    validation = constraints["validation"]
    _require(validation["seed"] == lock["validation_seed"] == 914027, "validation seed drift")
    _require(validation["held_out_points"] == lock["held_out_points"] == 4096, "held-out count drift")
    _require(validation["derivative_steps"] == lock["derivative_steps"] == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(validation["quadrature_orders_per_axis"] == lock["quadrature_orders_per_axis"] == [24, 48, 96], "quadrature ladder drift")
    thresholds = validation["thresholds"]
    _require(thresholds["pde_residual_max"] == lock["momentum_max_gate"] == 1.0e-3, "momentum max gate drift")
    _require(thresholds["pde_residual_L2"] == lock["momentum_volume_l2_gate"] == 1.0e-3, "momentum L2 gate drift")
    _require(thresholds["divergence_max"] == lock["divergence_max_gate"] == 1.0e-5, "divergence max gate drift")
    _require(thresholds["divergence_L2"] == lock["divergence_volume_l2_gate"] == 1.0e-5, "divergence L2 gate drift")
    _require("changing thresholds requires a new experiment version" in validation["failure_policy"], "post-hoc threshold firewall missing")


def _assert_provenance(scope: Mapping[str, Any]) -> None:
    provenance = scope["provenance"]
    _require(set(provenance) == PROVENANCE_CLASSES, "four-way provenance vocabulary drift")
    for name in PROVENANCE_CLASSES:
        _require(isinstance(provenance[name], list) and provenance[name], f"empty provenance class: {name}")

    public_text = " ".join(provenance["public_source_fact"])
    autonomous_text = " ".join(provenance["autonomous_design"])
    pending_text = " ".join(provenance["pending_unknown"])
    _require("0.02" not in public_text, "autonomous numeric widths laundered into public-source facts")
    _require("0.02" in autonomous_text, "autonomous width provenance missing")
    _require("source-hidden numerical values" in pending_text, "hidden source widths no longer pending")
    _require("source-admitted global kappa_0" in pending_text, "global kappa/cone admission no longer pending")


def _assert_width_and_cone_scope(scope: Mapping[str, Any]) -> None:
    realization = scope["registered_realization"]
    _require(realization["X_b_ref"] == 100.0 and realization["X_i"] == 110.0, "final-bridge source-coordinate interval drift")
    _require(realization["axial_shutdown_log_width"] == 0.02, "axial autonomous width drift")
    _require(realization["angular_settle_log_width"] == 0.02, "angular autonomous width drift")
    _require(realization["selected_kappa_0"] == 0.1, "selected autonomous kappa drift")
    _require(realization["bridge_quadrature_order"] == 48, "bridge quadrature drift")
    _require(realization["widths_are_repository_autonomous"] is True, "width provenance promotion")
    _require(realization["widths_recovered_from_public_source"] is False, "hidden width recovery invented")
    _require(realization["widths_fixed_before_residual_evaluation"] is True, "width preregistration lost")
    _require(realization["widths_fit_strictly_inside_100_to_110"] is True, "registered geometric-fit fact lost")

    evidence = scope["evidence_scope"]
    _require(evidence["candidate_side_final_bridge_materialized"] is True, "actual final bridge registration lost")
    _require(evidence["candidate_side_G_i_ell_i_materialized"] is True, "Xi handoff registration lost")
    _require(evidence["independent_profile_derivative_audit_protocol_registered"] is True, "A4 audit registration lost")
    _require(evidence["independent_endpoint_derivative_method_registered"] is True, "A4 endpoint method registration lost")
    _require(evidence["agent4_exact_head_scientific_admission_observed"] is False, "queued/unresolved A4 evidence laundered into admission")
    for key in (
        "geometric_width_fit_implies_source_cone_admission",
        "profile_derivative_consistency_implies_source_cone_admission",
        "endpoint_handoff_consistency_implies_source_cone_admission",
        "autonomous_widths_identify_source_hidden_widths",
        "final_bridge_cone_admissibility_independently_certified",
        "source_admitted_global_kappa0_smallness",
    ):
        _require(evidence[key] is False, f"unsupported cone/source promotion: {key}")


def _assert_mechanics_witness(scope: Mapping[str, Any]) -> None:
    witness = scope["mechanics_witness"]
    _require("not a Kokuno/OpenAI bound" in witness["scope"], "mechanics witness provenance drift")
    lo, hi = witness["toy_allowed_interval"]
    bound = witness["toy_admissible_upper_bound"]
    existing = witness["toy_existing_admissible_width"]
    selected = witness["toy_autonomous_selected_width"]
    _require(lo < existing < bound < selected < hi, "toy witness inequalities no longer establish the distinction")
    _require(witness["existing_width_satisfies_toy_admissibility"] is True, "toy existential witness drift")
    _require(witness["autonomous_selected_width_is_geometrically_inside_interval"] is True, "toy geometric-fit witness drift")
    _require(witness["autonomous_selected_width_satisfies_toy_admissibility"] is False, "toy witness no longer fail-closes selected-width admission")


def _assert_state_boundary(scope: Mapping[str, Any]) -> None:
    state = scope["state_boundary"]
    _require(state["canonical_eq45_velocity_export_ready"] is True, "canonical callable delivery was incorrectly cancelled")
    _require(state["kokuno_velocity_export_ready"] is False, "source-coordinate bridge promoted to Cartesian velocity delivery")
    for key in (
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pending_kokuno_cone_or_pde_work_blocks_canonical_callable_velocity_delivery",
    ):
        _require(state[key] is False, f"unsupported state promotion/coupling: {key}")


def validate_scope(scope: Mapping[str, Any], parent_source: str, constraints: Mapping[str, Any]) -> dict[str, Any]:
    _require(scope.get("schema") == SCHEMA, "scope schema drift")
    _require(scope.get("task_id") == "CR002-KOKUNO-FINAL-BRIDGE-WIDTH-CONE-SCOPE-090", "task id drift")
    _require(scope["audited_parent"]["pr"] == 944, "wrong A5 parent PR")
    _require(scope["audited_parent"]["head"] == "a732424661a615f24b7034f733bf5f281ac9a2d2", "wrong A5 parent head")
    _require(scope["audited_parent"]["source_blob"] == "fc0f463a451fd26d0d7a5e90570a2c1698f725e6", "wrong A5 parent source blob")
    _require(scope["cr001_lock"]["constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814", "canonical constraints blob drift")
    _assert_parent_contract(parent_source)
    _assert_provenance(scope)
    _assert_width_and_cone_scope(scope)
    _assert_mechanics_witness(scope)
    _assert_canonical_constraints(scope, constraints)
    _assert_state_boundary(scope)
    return {
        "schema": SCHEMA,
        "status": "pass",
        "claim": "local/autonomous final-bridge realization consistency is not source cone/global-smallness admission",
        "scientific_promotion": False,
    }


def audit_files(scope_path: str | Path, parent_path: str | Path, constraints_path: str | Path) -> dict[str, Any]:
    scope = _load_json(scope_path)
    constraints = _load_json(constraints_path)
    parent_source = Path(parent_path).read_text(encoding="utf-8")
    return validate_scope(scope, parent_source, constraints)


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", default="configs/kokuno_final_bridge_width_cone_scope.json")
    parser.add_argument("--parent", default="src/openai_ns_reconstruction/kokuno_a5_actual_final_bridge_xi110_ingest_contract.py")
    parser.add_argument("--constraints", default="configs/constraints.json")
    args = parser.parse_args()
    print(json.dumps(audit_files(args.scope, args.parent, args.constraints), sort_keys=True))


if __name__ == "__main__":
    _main()
