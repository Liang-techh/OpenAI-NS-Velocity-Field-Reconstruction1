"""Agent-5 registration seam for the current Kokuno numerical T_sh certificate.

This module is integration/provenance glue only. It checksum-binds the current
Agent-1 candidate-side numerical B_0 -> T_sh certificate, the matching Agent-4
implementation-distinct audit protocol, the existing Agent-3 current-Xi
five-moment handoff, and the latest Agent-2 axis-safe oscillatory evaluator.

Registration is deliberately weaker than scientific admission. Finite sampling
of ell_i does not prove the source analytic C^0 bound, does not apply PA.16, and
does not materialize a global Cartesian leading field, matched pressure,
restricted forcing, correction velocity, or complete Navier--Stokes residual.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-tsh-certificate-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-TSH-CERTIFICATE-INGEST-090"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 951,
    "head": "b3fb36d91239a2ce2669971b7b07bfe892cfdf1c",
    "branch": "codex/kokuno-a5-current-xi-moment-discrepancy-ingest-089",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_xi_moment_discrepancy_ingest_contract.py",
    "observed_ci": {
        "dedicated": {"run_id": 35559565416, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35559565469, "status": "queued", "conclusion": None},
    },
}

AGENT1_CURRENT_TSH = {
    "pr": 953,
    "head": "bf16b88db92d23923b3e049c00e000a5bfa09992",
    "branch": "agent-kokuno-1/current-tsh-certificate-090",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_tsh_certificate.py",
    "source_blob": "5fd610b417ccaa1c4673ea64985296b94dfe5946",
    "test_path": "tests/test_kokuno_pa16_current_tsh_certificate.py",
    "test_blob": "dbf43e26bb6b37f6167c923cc9880e2911a09177",
    "workflow_path": ".github/workflows/kokuno-agent1-current-tsh-certificate.yml",
    "workflow_blob": "e2c666e31cd124081b88102d09622f7255bac9fe",
    "upstream_xi_moment_pr": 947,
    "upstream_xi_moment_head": "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "public_source_tex_blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "corrected_reader_date": "2026-09-09",
    "source_formula": "T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)",
    "sigma_prime_sup": 8.0,
    "log_f_sup": "log(2)",
    "source_eta_interval": [-1.0, 1.0],
    "separation_condition": "log(X_i/X_R)+T_sh < -8",
    "X_i": 110.0,
    "candidate_side_numerical_B0_envelope_materialized": True,
    "candidate_side_numerical_T_sh_instantiated": True,
    "selected_join_route_executable_when_geometry_allows": True,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "pa16_repair_applied": False,
    "inner_to_outer_join_completed": False,
    "observed_ci": {
        "dedicated": {"run_id": 35561326177, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35561326100, "status": "queued", "conclusion": None},
    },
}

AGENT4_CURRENT_TSH_AUDIT = {
    "pr": 955,
    "head": "811fcd9bf5cc3d48264d2f7f7d3763d2167db83e",
    "branch": "agent-kokuno-4/current-tsh-independent-audit-092",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_tsh_independent_audit.py",
    "source_blob": "e517a400be283131d530e0b34f3d6c0712667a7b",
    "test_path": "tests/test_constrained_kokuno_a4_current_tsh_independent_audit.py",
    "test_blob": "d1441d0047deaa93d27622cedd2176e99c2e57de",
    "workflow_path": ".github/workflows/kokuno-agent4-current-tsh-independent-audit.yml",
    "workflow_blob": "5f5940b1f774640c4e2385d3be9ada2bca4b220c",
    "audited_agent1_pr": 953,
    "audited_agent1_head": "bf16b88db92d23923b3e049c00e000a5bfa09992",
    "saved_reloaded_public_ell_i_only": True,
    "agent1_envelope_cache_used_as_reference": False,
    "independent_reference": "nested uniform eta grids plus fresh stratified off-grid probes and separately recoded public algebra",
    "protocol": {
        "seed": 9173641,
        "nested_uniform_eta_grids": [129, 257, 513],
        "fresh_stratified_off_grid_eta_count": 257,
        "selected_B0_must_cover_all_independent_probes": True,
        "medium_to_fine_scale_normalized_envelope_change_gate": 2.0e-3,
        "ell_i_nontrivial_max_floor": 1.0e-12,
        "T_sh_formula_relative_agreement_gate": 2.0e-13,
        "separation_geometry_absolute_agreement_gate": 2.0e-11,
        "exact_save_reload_semantic_identity": True,
    },
    "audit_protocol_registered": True,
    "audit_admitted": False,
    "candidate_side_route_feasibility": "unresolved_pending_exact_head_receipt",
    "source_analytic_C0_bound_proved": False,
    "complete_ns_residual_audit": False,
    "observed_ci": {
        "dedicated": {"run_id": 35562409849, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35562409842, "status": "queued", "conclusion": None},
    },
}

AGENT3_CURRENT_XI_HANDOFF = {
    "pr": 949,
    "head": "87c1ef765d6bf3613c7c6c301669ffccf5260117",
    "role": "typed current-Xi five-moment discrepancy consumer and public PA.16 row transform",
    "public_pa16_row_transform": ["M", "J-4*eta*I", "I", "S-8*eta*M", "C_p"],
    "pa16_inverse_or_fixed_point_solved": False,
    "pa16_repair_applied": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "real_correction_velocity_materialized": False,
    "observed_ci": {
        "dedicated": {"run_id": 35558337953, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35558338012, "status": "queued", "conclusion": None},
    },
}

AGENT2_SIBLING_STATUS = {
    "pr": 954,
    "head": "e43c32c3d0258f635e8b1f999d89c38524355856",
    "role": "batch axis-safe adapter for the frozen public oscillatory velocity",
    "source_path": "src/openai_ns_reconstruction/kokuno_oscillatory_batch_axis_safety.py",
    "source_blob": "599baa190ec742d0e532e5517078534124e68d6b",
    "test_path": "tests/test_constrained_kokuno_oscillatory_batch_axis_safety.py",
    "test_blob": "1c08a6ce8f7ffc30cec9cd38e19e4b7dc156b9c1",
    "workflow_path": ".github/workflows/kokuno-agent2-oscillatory-batch-axis-safety.yml",
    "workflow_blob": "439686dbdc3a4d8b0a536e412585465de444cf29",
    "axis_and_support_masked_before_interior_cylindrical_eval": True,
    "creates_global_leading_velocity": False,
    "creates_correction_velocity": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35562205146, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35562205067, "status": "queued", "conclusion": None},
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
    "current_candidate_side_xi_discrepancy_registered": True,
    "candidate_side_numerical_tsh_certificate_registered": True,
    "agent4_current_tsh_independent_audit_registered": True,
    "agent4_current_tsh_independent_audit_admitted": False,
    "upstream_ci_admitted_as_pass": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "pa16_repair_applied": False,
    "five_moment_repair_closed": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "actual_inner_to_outer_global_join_materialized": False,
    "corrected_global_cartesian_leading_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "real_agent3_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "same_protocol_full_ns_residual_available": False,
    "current_candidate_eligible_for_full_ns_validation": False,
    "same_protocol_st006_comparison_available_now": False,
    "scientific_admission": False,
    "paper_exact": False,
}

PIPELINE_POSITION = {
    "stage": "current candidate-side numerical T_sh certificate registration",
    "input": "A5 #951 current Xi moment-discrepancy seam",
    "new_output": "A1 #953 numerical B0/T_sh certificate paired with A4 #955 implementation-distinct audit protocol",
    "not_output": "source analytic B0/T_sh theorem, PA.16 repair, global Cartesian leading velocity, matched pressure/forcing, correction velocity, or full NS validation",
    "next_shortest_blocker": "resolve exact-head A4 evidence, then actually solve/apply the current-lineage PA.16 repair and materialize the inner-to-outer/global Cartesian leading join",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: str, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA")


def build_contract(exact_head: str) -> dict[str, Any]:
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_current_tsh": copy.deepcopy(AGENT1_CURRENT_TSH),
        "agent4_current_tsh_audit": copy.deepcopy(AGENT4_CURRENT_TSH_AUDIT),
        "agent3_current_xi_handoff": copy.deepcopy(AGENT3_CURRENT_XI_HANDOFF),
        "agent2_sibling_status": copy.deepcopy(AGENT2_SIBLING_STATUS),
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

    expected_top = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent1_current_tsh": AGENT1_CURRENT_TSH,
        "agent4_current_tsh_audit": AGENT4_CURRENT_TSH_AUDIT,
        "agent3_current_xi_handoff": AGENT3_CURRENT_XI_HANDOFF,
        "agent2_sibling_status": AGENT2_SIBLING_STATUS,
        "frozen_science": FROZEN_SCIENCE,
        "final_gate": FINAL_GATE,
        "st006_baseline": ST006_BASELINE,
        "readiness": READINESS,
        "truth_boundary": TRUTH_BOUNDARY,
        "pipeline_position": PIPELINE_POSITION,
        "registration_only": True,
    }
    for key, expected in expected_top.items():
        if observed.get(key) != expected:
            errors.append(f"{key}_drift")

    exact_head = observed.get("agent5_exact_head")
    try:
        _require_hex40(exact_head, "agent5_exact_head")
    except (TypeError, ValueError):
        errors.append("agent5_exact_head_invalid")

    if observed.get("truth_boundary", {}).get("source_T_sh_lower_bound_verified") is not False:
        errors.append("source_tsh_laundering")
    if observed.get("truth_boundary", {}).get("pa16_repair_applied") is not False:
        errors.append("pa16_repair_laundering")
    if observed.get("truth_boundary", {}).get("scientific_admission") is not False:
        errors.append("scientific_admission_laundering")
    if observed.get("readiness", {}).get("pde_validated") is not False:
        errors.append("pde_validation_laundering")
    if not observed.get("frozen_science", {}).get("residual_defined_free_forcing_forbidden", False):
        errors.append("free_forcing_firewall_dropped")
    return sorted(set(errors))


def write_contract(path: str | Path, exact_head: str) -> dict[str, Any]:
    payload = build_contract(exact_head)
    errors = validate_contract(payload)
    if errors:
        raise ValueError(f"contract failed validation: {errors}")
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = write_contract(args.output, args.exact_head)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
