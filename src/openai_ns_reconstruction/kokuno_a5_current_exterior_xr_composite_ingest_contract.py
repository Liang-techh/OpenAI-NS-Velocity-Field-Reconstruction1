"""Agent-5 integration contract for the current through-X_R composite seam.

This module is integration/provenance glue only.  It registers the first current
lineage in which all of the following coexist on the same scoped domain:

    A1 #980 Cartesian leading velocity through X_R
      -> A2 #987 identity-bound leading + complete-curl oscillatory velocity
      -> A3 #988 nonlinear m=0 mean attribution of that exact composite
      -> A4 #989 public-velocity-only independent staged divergence audit

The contract does not reimplement any mathematical lane.  It also records A1
#986 as the newer leading-only first post-X_R RF40 turn, but does not silently
retarget A2/A3/A4 onto that newer domain.  Therefore the integrated composite
registered here still stops at X_R.

Registration is not scientific admission.  In particular, A4 #989 is scoped
incompressibility evidence only while exact-head Actions remain unresolved; it
is not a momentum/full-NS residual, canonical [24,48,96] quadrature evidence,
or PDE validation.  The fixed 1e-3 momentum and 1e-5 divergence gates and the
ban on residual-defined free forcing remain unchanged.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-exterior-xr-composite-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-EXTERIOR-XR-COMPOSITE-INGEST-096"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 984,
    "head": "dae4797fbd270a4f676b113f5ca3f7801d71a47e",
    "branch": "codex/kokuno-a5-current-cartesian-exterior-xr-ingest-095",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_cartesian_exterior_xr_ingest_contract.py",
    "source_blob": "0a1c6492680c27ab98988afd64d27004a4ada0ef",
    "role": "register current Cartesian leading through X_R plus A4 leading-only divergence audit",
}

AGENT2_XR_COMPOSITE = {
    "pr": 987,
    "head": "4be2c9ee898c44dd1ad2217e90601b161fe81964",
    "branch": "codex/kokuno-a2-exterior-xr-composite-083",
    "base_pr": 981,
    "source_path": "src/openai_ns_reconstruction/kokuno_current_exterior_xr_leading_oscillatory_identity.py",
    "source_blob": "93a3da5c3aed4be8310cc510edfab6af940b64dc",
    "test_path": "tests/test_constrained_kokuno_current_exterior_xr_leading_oscillatory_identity.py",
    "test_blob": "ff9f6e419ebf4ff24dfaf5716e8a7b047c1f2893",
    "workflow_path": ".github/workflows/kokuno-agent2-current-exterior-xr-leading-oscillatory.yml",
    "workflow_blob": "c00bea38acfa1342ba7f77a17238029c58c157f3",
    "leading_pr": 980,
    "leading_head": "d3c971f2c62e272333e124e532212d23cca4908d",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "materialized_domain": "current leading + frozen complete-curl oscillation through X_R",
    "identity_preserving_save_load": True,
    "complete_oscillatory_runtime_digest_bound": True,
    "fails_closed_after_X_R": True,
    "post_XR_composite_materialized": False,
    "matched_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
}

AGENT3_XR_NONLINEAR_MEAN = {
    "pr": 988,
    "head": "fa95dee3709a82326167af1a2cb42b9b5f9a88a5",
    "branch": "codex/kokuno-a3-current-exterior-xr-nonlinear-mean-105",
    "base_pr": 982,
    "source_path": "src/openai_ns_reconstruction/kokuno_current_exterior_xr_nonlinear_mean_attribution.py",
    "source_blob": "fdb01b9f02942060e8c343533c91fd44431b908e",
    "test_path": "tests/test_constrained_kokuno_current_exterior_xr_nonlinear_mean_attribution.py",
    "test_blob": "71a47b51a9335d58f5b6809ff515d7dee80c6c62",
    "workflow_path": ".github/workflows/kokuno-agent3-current-exterior-xr-nonlinear-mean.yml",
    "workflow_blob": "b6d959a720c91ffa9b7e4d5d2b0e8eb1f5a94fac",
    "consumed_agent2_pr": 987,
    "consumed_agent2_head": "4be2c9ee898c44dd1ad2217e90601b161fe81964",
    "decomposition": "A=(u_lead·grad)u_osc; B=(u_osc·grad)u_lead; Q=(u_osc·grad)u_osc; N=A+B+Q",
    "projector": "existing rotating cylindrical m=0 projector",
    "materialized_domain": "through X_R",
    "radial_inverse_in_this_increment": False,
    "complete_ns_defect": False,
    "authorized_as_correction_target": False,
    "real_cartesian_correction_velocity": False,
}

AGENT4_XR_COMPOSITE_AUDIT = {
    "pr": 989,
    "head": "84b8aab2f5a56dbfc8ef02b6053d19bea910a959",
    "branch": "codex/kokuno-a4-exterior-xr-composite-divergence-audit-098",
    "audited_agent2_pr": 987,
    "audited_agent2_head": "4be2c9ee898c44dd1ad2217e90601b161fe81964",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_exterior_xr_leading_oscillatory_divergence_independent_audit.py",
    "source_blob": "be05390250ce21056e315a8a2b893c36cb1cc7e8",
    "test_path": "tests/test_constrained_kokuno_a4_current_exterior_xr_leading_oscillatory_divergence_independent_audit.py",
    "test_blob": "b738f3613684175701fac4a6e26b2e8762b8c908",
    "workflow_path": ".github/workflows/kokuno-agent4-current-exterior-xr-leading-oscillatory-divergence-audit.yml",
    "workflow_blob": "d17a0cc8067cdc5266e8b852dba5fc9b47258d37",
    "reference_path": "save/reloaded public velocity only; centered Cartesian FD2 Jacobians",
    "uses_agent2_production_jacobian_or_divergence": False,
    "uses_agent1_source_derivative_helpers": False,
    "uses_pressure_forcing_or_residual_helper": False,
    "seed": 9173701,
    "times": [0.31, 0.47, 0.63, 0.71],
    "fd2_step_ladder": [0.02, 0.01, 0.005],
    "scoped_divergence_gate": 1.0e-5,
    "nontrivial_total_speed_rms_floor": 1.0e-10,
    "nontrivial_oscillatory_speed_rms_floor": 1.0e-12,
    "medium_to_fine_degradation_ratio_cap": 1.25,
    "numerical_floor": 2.0e-8,
    "scoped_divergence_evidence_only": True,
    "canonical_whole_domain_admission": False,
    "momentum_or_full_ns_evidence": False,
    "covers_agent1_986_post_XR_first_turn": False,
    "scientific_admission": False,
}

AGENT1_POST_XR_FIRST_TURN_SIBLING = {
    "pr": 986,
    "head": "23c00b98526e187ff04d432745300a084ec859f2",
    "branch": "agent-kokuno-1/current-cartesian-rf40-first-turn-097",
    "base_pr": 980,
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_rf40_first_turn.py",
    "source_blob": "bf71948cbdb827d791b86ab8f00496bd1b231961",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_rf40_first_turn.py",
    "test_blob": "14703a4151ef86fecc2a3731f7e0bd3f81602ad2",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-rf40-first-turn.yml",
    "workflow_blob": "c723e37ec024e5bbde4a4c05e83401b6d0e9c6ab",
    "materialized_domain": "leading-only through X_1=e X_R",
    "preserves_current_incompressibility_memory": True,
    "complete_post_XR_RF40": False,
    "cone_i1_i2_i3_i4_completed": False,
    "consumed_by_agent2_987": False,
    "consumed_by_agent3_988": False,
    "covered_by_agent4_989": False,
    "consumed_by_this_increment": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_984": {
        "repository_tests_run": 35586185733,
        "dedicated_run": 35586185713,
        "status": "queued",
        "conclusion": None,
    },
    "agent1_986": {
        "repository_tests_run": 35589029470,
        "dedicated_run": 35589029596,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_987": {
        "repository_tests_run": 35589369692,
        "dedicated_run": 35589369729,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_988": {
        "repository_tests_run": 35590088520,
        "dedicated_run": 35590088517,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_989": {
        "repository_tests_run": 35591242652,
        "dedicated_run": 35591242635,
        "status": "queued",
        "conclusion": None,
    },
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
    "current_cartesian_leading_velocity_materialized_through_xr": True,
    "current_leading_plus_oscillatory_velocity_through_xr_materialized": True,
    "identity_preserving_composite_save_load_through_xr_available": True,
    "current_nonlinear_m0_mean_attribution_through_xr_materialized": True,
    "agent4_989_independent_composite_divergence_audit_registered": True,
    "agent4_989_independent_composite_divergence_audit_admitted": False,
    "agent1_986_leading_only_post_xr_first_turn_materialized_as_sibling": True,
    "current_leading_plus_oscillatory_velocity_post_xr_materialized": False,
    "current_nonlinear_mean_post_xr_materialized": False,
    "independent_post_xr_composite_audit_available": False,
    "complete_post_xr_rf40_current_lineage_materialized": False,
    "cone_i1_i2_i3_i4_outer_overlays_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "current_xr_nonlinear_mean_authorized_as_correction_target": False,
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
    "stage": "identity-bound current leading+oscillatory composite and nonlinear mean through X_R, with matching independent staged divergence audit registered",
    "input": "A5 #984 + A2 #987 + A3 #988 + A4 #989; A1 #986 retained as newer non-consumed leading-only post-X_R sibling",
    "new_output": "checksum-bound through-X_R composite/save-load -> nonlinear m=0 mean -> independent public-velocity staged divergence registration",
    "not_output": "post-X_R composite/mean/audit, remaining RF40 plus cone/I1-I4 global completion, pressure, forcing, complete defect, correction velocity, full residual, or PDE validation",
    "next_shortest_blocker": "recompose the exact oscillatory runtime onto A1 #986 or a later completed post-X_R/global leading identity, then carry A3 mean/radial plumbing and A4 independent validation forward before matched pressure/restricted forcing and complete-defect correction",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA")


def build_contract(exact_head: str) -> dict[str, Any]:
    """Build the immutable Agent-5 registration receipt for this increment."""
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent2_xr_composite": copy.deepcopy(AGENT2_XR_COMPOSITE),
        "agent3_xr_nonlinear_mean": copy.deepcopy(AGENT3_XR_NONLINEAR_MEAN),
        "agent4_xr_composite_audit": copy.deepcopy(AGENT4_XR_COMPOSITE_AUDIT),
        "agent1_post_xr_first_turn_sibling": copy.deepcopy(AGENT1_POST_XR_FIRST_TURN_SIBLING),
        "observed_ci_at_freeze": copy.deepcopy(OBSERVED_CI_AT_FREEZE),
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
    """Return fail-closed validation errors for a serialized receipt."""
    errors: list[str] = []
    observed = copy.deepcopy(dict(payload))
    supplied_sha = observed.pop("contract_sha256", None)
    if supplied_sha != _sha256(observed):
        errors.append("contract_sha256_mismatch")

    expected = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent2_xr_composite": AGENT2_XR_COMPOSITE,
        "agent3_xr_nonlinear_mean": AGENT3_XR_NONLINEAR_MEAN,
        "agent4_xr_composite_audit": AGENT4_XR_COMPOSITE_AUDIT,
        "agent1_post_xr_first_turn_sibling": AGENT1_POST_XR_FIRST_TURN_SIBLING,
        "observed_ci_at_freeze": OBSERVED_CI_AT_FREEZE,
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

    # Explicit scientific promotion firewalls, kept redundant with exact-payload
    # matching so accidental schema edits fail with useful diagnostics.
    readiness = observed.get("readiness", {})
    if readiness.get("leading_ready") is not False:
        errors.append("leading_ready_promoted")
    if readiness.get("correction_ready") is not False:
        errors.append("correction_ready_promoted")
    if readiness.get("velocity_export_ready") is not False:
        errors.append("velocity_export_ready_promoted")
    if readiness.get("pde_validated") is not False:
        errors.append("pde_validated_promoted")

    truth = observed.get("truth_boundary", {})
    for key in (
        "agent4_989_independent_composite_divergence_audit_admitted",
        "current_leading_plus_oscillatory_velocity_post_xr_materialized",
        "independent_post_xr_composite_audit_available",
        "outer_global_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "current_xr_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "heldout_normalized_ns_residual_assessed",
        "scientific_admission",
    ):
        if truth.get(key) is not False:
            errors.append(f"truth_promotion:{key}")

    if observed.get("agent4_xr_composite_audit", {}).get(
        "covers_agent1_986_post_XR_first_turn"
    ) is not False:
        errors.append("agent4_scope_laundering")
    if observed.get("agent1_post_xr_first_turn_sibling", {}).get(
        "consumed_by_agent2_987"
    ) is not False:
        errors.append("agent1_986_lineage_laundering")

    final_gate = observed.get("final_gate", {})
    if final_gate.get("normalized_momentum_sampled_max") != 1.0e-3:
        errors.append("momentum_gate_changed")
    if final_gate.get("normalized_momentum_volume_l2") != 1.0e-3:
        errors.append("momentum_l2_gate_changed")
    if final_gate.get("divergence_sampled_max") != 1.0e-5:
        errors.append("divergence_gate_changed")
    if final_gate.get("divergence_volume_l2") != 1.0e-5:
        errors.append("divergence_l2_gate_changed")
    if final_gate.get("canonical_volume_quadrature_ladder") != [24, 48, 96]:
        errors.append("canonical_quadrature_changed")
    if observed.get("frozen_science", {}).get(
        "residual_defined_free_forcing_forbidden"
    ) is not True:
        errors.append("free_forcing_firewall_relaxed")

    return errors


def write_contract(path: str | Path, exact_head: str) -> dict[str, Any]:
    payload = build_contract(exact_head)
    errors = validate_contract(payload)
    if errors:
        raise ValueError("invalid Agent-5 contract: " + ", ".join(errors))
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    payload = write_contract(args.output, args.exact_head)
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "task_id": payload["task_id"],
                "contract_sha256": payload["contract_sha256"],
                "output": str(args.output),
                "registration_only": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
