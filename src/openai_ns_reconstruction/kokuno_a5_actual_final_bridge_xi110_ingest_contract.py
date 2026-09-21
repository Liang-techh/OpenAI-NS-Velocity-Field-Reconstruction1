"""Agent-5 registration seam for the current actual Kokuno final bridge through Xi=110.

Integration glue only.  This module pins Agent 1 PR #940's executable
candidate-side source-coordinate bridge on ``100 <= X <= 110`` and Agent 4 PR
#943's implementation-distinct save/reload public-values-only audit of that
same bridge, including one-sided endpoint differentiation at Xi.

Registration is not scientific admission: the matching A4 exact-head CI is
still unresolved at the frozen snapshot.  The source-coordinate ``F/U/E`` +
``G_i/ell_i`` handoff is also not a Cartesian ``velocity(x,y,z,t)`` artifact.
The reconstruction remains provenance-labelled and is not paper-exact.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-actual-final-bridge-xi110-ingest-v1"
TASK_ID = "KOKUNO-A5-ACTUAL-FINAL-BRIDGE-XI110-INGEST-088"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 937,
    "head": "44b9349ef07d947d0da7eff053ae2409348374b8",
    "branch": "codex/kokuno-a5-reference-stress-xi110-ingest-087",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_reference_stress_xi110_ingest_contract.py",
    "source_blob": "531106d200b4163018f0331537b7bf6706042bd0",
    "test_path": "tests/test_kokuno_a5_reference_stress_xi110_ingest_contract.py",
    "test_blob": "32864416a1270cfacb5cd896a36ecec28a64dda8",
    "workflow_path": ".github/workflows/kokuno-agent5-reference-stress-xi110-ingest.yml",
    "workflow_blob": "7fee1034269a6dd7cb2c183c2d5e16ec4c3a285d",
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35552584976, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35552585005, "status": "queued", "conclusion": None},
    },
}

AGENT1_ACTUAL_FINAL_BRIDGE_XI110 = {
    "pr": 940,
    "head": "f430f8f97df3bfafce56bf094047c769ecec922d",
    "branch": "codex/kokuno-a1-actual-final-bridge-xi110-086",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_actual_final_bridge_xi110.py",
    "source_blob": "0cbd49ad0726c1573f2e649c96a3b516d8df053e",
    "test_path": "tests/test_kokuno_pa10_actual_final_bridge_xi110.py",
    "test_blob": "227d132a7ef1ea74323cfeb6c94f0deda490d3c1",
    "workflow_path": ".github/workflows/kokuno-agent1-actual-final-bridge-xi110.yml",
    "workflow_blob": "99c03716ac33a029931a67bc90d9511401963f89",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "public_source_tex_blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "corrected_reader_date": "2026-09-09",
    "source_coordinate_range": [100.0, 110.0],
    "X_i": 110.0,
    "kappa_0": 0.1,
    "axial_shutdown_log_width": 0.02,
    "angular_settle_log_width": 0.02,
    "bridge_quadrature_order": 48,
    "quadrature_replay_order": 64,
    "public_profiles": ["F_final_bridge", "U_final_bridge", "E_final_bridge"],
    "public_radial_derivatives_registered": True,
    "G_i_registered": True,
    "ell_i_registered": True,
    "semantic_identity_and_save_load_registered": True,
    "candidate_side_actual_final_interpolation_materialized": True,
    "autonomous_widths_fixed_before_residual": True,
    "selected_kappa0_is_repository_autonomous": True,
    "selected_kappa0_chosen_from_ns_residual": False,
    "source_prepared_appendixA_pressure": False,
    "source_prepared_full_reference_stress_pair": False,
    "source_admitted_global_kappa0_smallness": False,
    "final_bridge_cone_admissibility_independently_certified": False,
    "actual_inner_to_outer_global_join_materialized": False,
    "global_cartesian_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35554352985, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35554353013, "status": "queued", "conclusion": None},
    },
}

AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT = {
    "pr": 943,
    "head": "d090572fdc766626c417cb53b8acc8632b2bdb63",
    "branch": "codex/kokuno-a4-actual-final-bridge-xi110-independent-audit-088",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_actual_final_bridge_xi110_independent_audit.py",
    "source_blob": "afe00af9274408051bb01ec9924ac287c538fe5f",
    "test_path": "tests/test_constrained_kokuno_a4_actual_final_bridge_xi110_independent_audit.py",
    "test_blob": "0a5de622f718321c74de8733718daf989928a97b",
    "workflow_path": ".github/workflows/kokuno-agent4-actual-final-bridge-xi110-independent-audit.yml",
    "workflow_blob": "2e96ff193b603cd283ac8e342434a0336a9df931",
    "audited_agent1_pr": 940,
    "audited_agent1_head": "f430f8f97df3bfafce56bf094047c769ecec922d",
    "saved_reloaded_public_values_only": True,
    "production_derivatives_used_as_reference": False,
    "independent_operator": (
        "degree-six local polynomial derivatives reconstructed from public values only; "
        "asymmetric seven-node interior stencil plus distinct forward/backward one-sided endpoint stencils"
    ),
    "protocol": {
        "seed": 9173621,
        "random_per_region": 40,
        "interior_sample_count": 130,
        "endpoint_sample_count": 7,
        "interior_halfwidths": [0.12, 0.06, 0.03],
        "endpoint_steps": [0.08, 0.04, 0.02],
        "derivative_relative_rms_gate": 5.0e-3,
        "derivative_relative_max_gate": 2.0e-2,
        "endpoint_scale_normalized_rms_gate": 5.0e-3,
        "endpoint_scale_normalized_max_gate": 2.0e-2,
        "refinement_ratio_gate": 8.0,
        "refinement_floor": 2.0e-9,
        "value_identity_relative_gate": 5.0e-12,
        "derivative_chain_relative_gate": 5.0e-5,
        "endpoint_value_handoff_abs_gate": 2.0e-12,
        "xi_final_slope_abs_gate": 2.0e-4,
        "xi_zero_Ux_scale_normalized_gate": 5.0e-3,
        "value_rms_floor": 1.0e-10,
        "derivative_rms_floor": 1.0e-12,
    },
    "endpoint_safe_independent_method_registered": True,
    "negative_result_receipt_retained_before_gate": True,
    "actual_final_bridge_independent_a4_audit_registered": True,
    "actual_final_bridge_independent_a4_audit_admitted": False,
    "actual_Xi_endpoint_derivative_independently_verified": False,
    "complete_ns_residual_audit": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35555427211, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35555427197, "status": "queued", "conclusion": None},
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
    "candidate_side_actual_final_bridge_100_to_Xi_registered": True,
    "candidate_side_actual_G_i_at_Xi_registered": True,
    "candidate_side_actual_ell_i_at_Xi_registered": True,
    "semantic_identity_and_save_load_registered": True,
    "actual_final_bridge_independent_a4_audit_registered": True,
    "Xi_endpoint_independent_derivative_method_registered": True,
    "upstream_ci_admitted_as_pass": False,
    "actual_final_bridge_independent_a4_audit_admitted": False,
    "actual_Xi_endpoint_derivative_independently_verified": False,
    "source_coordinate_profiles_are_cartesian_velocity": False,
    "source_prepared_appendixA_pressure_materialized": False,
    "source_prepared_appendixA_full_reference_stress_pair_materialized": False,
    "source_admitted_global_kappa0_smallness": False,
    "final_bridge_cone_admissibility_independently_certified": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
    "five_moment_repair_closed": False,
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
    "stage": "candidate-side source-coordinate final-bridge plus independent-audit registration",
    "input": "A5 #937 registered Xi110 reference-stress prerequisite",
    "new_output": "A1 #940 typed F/U/E bridge plus G_i/ell_i handoff, paired with A4 #943 independent audit protocol",
    "not_output": "global Cartesian velocity/pressure/forcing candidate",
    "next_shortest_blocker": (
        "A1/A4 exact-head evidence resolution, then actual upstream five-moment discrepancy and "
        "inner-to-outer/global Cartesian leading join with matched pressure"
    ),
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
        "agent1_actual_final_bridge_xi110": copy.deepcopy(AGENT1_ACTUAL_FINAL_BRIDGE_XI110),
        "agent4_actual_final_bridge_audit": copy.deepcopy(AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT),
        "actual_final_bridge_ingest_seam": {
            "registration_only": True,
            "actual_candidate_source_coordinate_bridge_registered": True,
            "G_i_ell_i_handoff_registered": True,
            "matching_independent_a4_audit_protocol_registered": True,
            "a4_exact_head_pass_required_before_audit_admission": True,
            "queued_or_unresolved_ci_is_pass": False,
            "profile_success_substitutes_for_cartesian_delivery": False,
            "profile_or_visual_success_substitutes_for_pde_validation": False,
            "residual_defined_forcing_forbidden": True,
        },
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
        "frozen_science": copy.deepcopy(FROZEN_SCIENCE),
        "final_gate": copy.deepcopy(FINAL_GATE),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
    }
    payload["contract_sha256"] = _sha256(payload)
    return payload


def validate_contract(contract: Mapping[str, Any]) -> None:
    if not isinstance(contract, Mapping):
        raise TypeError("contract must be a mapping")
    exact_head = contract.get("agent5_exact_head")
    _require_hex40(exact_head, "agent5_exact_head")
    expected = build_contract(exact_head)
    if dict(contract) != expected:
        raise ValueError("contract drifted from frozen Agent-5 actual-final-bridge Xi110 registration")
    unsigned = dict(contract)
    digest = unsigned.pop("contract_sha256")
    if digest != _sha256(unsigned):
        raise ValueError("contract_sha256 mismatch")


def write_contract(path: str | Path, *, exact_head: str) -> dict[str, Any]:
    contract = build_contract(exact_head)
    validate_contract(contract)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return contract


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    contract = write_contract(args.output, exact_head=args.exact_head)
    print(json.dumps({"contract_sha256": contract["contract_sha256"], "output": args.output}, sort_keys=True))


if __name__ == "__main__":
    _main()
