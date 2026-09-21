"""Agent-5 ingest contract for the current Kokuno Xi five-moment discrepancy.

Integration glue only.  This seam pins the current Agent-1 candidate-side
Xi prefix-moment/discrepancy artifact, Agent-3's typed PA.16 row-transform
handoff, and Agent-4's implementation-distinct save/reload audit.

The discrepancy is a real current-lineage *source-moment* mismatch.  It is not
a defect derived from a complete Navier--Stokes candidate.  This module does
not choose ``T_sh``, solve/apply PA.16, construct a correction velocity, build
a global Cartesian leading field, fit pressure/forcing, or promote PDE state.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-xi-moment-discrepancy-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-XI-MOMENT-DISCREPANCY-INGEST-089"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 944,
    "head": "a732424661a615f24b7034f733bf5f281ac9a2d2",
    "branch": "codex/kokuno-a5-actual-final-bridge-xi110-ingest-088",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_actual_final_bridge_xi110_ingest_contract.py",
    "source_blob": "fc0f463a451fd26d0d7a5e90570a2c1698f725e6",
    "test_path": "tests/test_kokuno_a5_actual_final_bridge_xi110_ingest_contract.py",
    "test_blob": "8372a6a612839077d1d7e88f50dfe61edeab7c19",
    "workflow_path": ".github/workflows/kokuno-agent5-actual-final-bridge-xi110-ingest.yml",
    "workflow_blob": "5391fe52e78b280628e7e62c71f7686842b1a3ef",
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35555800805, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35555800653, "status": "queued", "conclusion": None},
    },
}

AGENT1_XI_PREFIX_MOMENTS = {
    "pr": 947,
    "head": "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    "branch": "agent-kokuno-1/actual-xi-moments-089",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_actual_xi_prefix_moments.py",
    "source_blob": "92dce65c9a8793e06497a7e7347c832861032238",
    "test_path": "tests/test_kokuno_pa10_actual_xi_prefix_moments.py",
    "test_blob": "169a6282e1fd0e20225c669ff4473958d6e8e014",
    "workflow_path": ".github/workflows/kokuno-agent1-actual-xi-prefix-moments.yml",
    "workflow_blob": "58a1b0c07cd7137d943a437d78f8db223d584652",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "public_source_tex_blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "corrected_reader_date": "2026-09-09",
    "X_i": 110.0,
    "source_normalization_C": 1000.0,
    "moment_names": ["M", "I", "J", "S", "C_p"],
    "moment_formulas": {
        "M": "int U dX",
        "I": "int H dX; H=sqrt(2X)E=2XF",
        "J": "int U H dX",
        "S": "int (U^2-E^2/2) dX",
        "C_p": "int E^2/(2X) dX = int F^2 dX",
    },
    "candidate_side_actual_minus_ideal_discrepancy_materialized": True,
    "G_i_ell_i_bound_into_pa16_shaped_input": True,
    "pa16_repair_applied": False,
    "source_T_sh_lower_bound_verified": False,
    "global_cartesian_velocity_materialized": False,
    "matched_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "observed_ci_at_registration": {
        "dedicated": {"run_id": 35557679528, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35557679567, "status": "queued", "conclusion": None},
    },
}

AGENT3_XI_DISCREPANCY_BRIDGE = {
    "pr": 949,
    "head": "87c1ef765d6bf3613c7c6c301669ffccf5260117",
    "branch": "codex/kokuno-a3-current-xi-moment-bridge-091",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_xi_moment_discrepancy_bridge.py",
    "source_blob": "223058e2132bb22298c6a642b0f2a55562181e25",
    "test_path": "tests/test_constrained_kokuno_current_xi_moment_discrepancy_bridge.py",
    "test_blob": "875bd551924a988b6adbc8f8f98180cddd78b04c",
    "workflow_path": ".github/workflows/kokuno-agent3-current-xi-moment-discrepancy-bridge.yml",
    "workflow_blob": "23493753b08408826d54464a2646f95a3f8c33aa",
    "upstream_agent1_pr": 947,
    "upstream_agent1_head": "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    "public_pa16_row_transform": [
        "M",
        "J-4*eta*I",
        "I",
        "S-8*eta*M",
        "C_p",
    ],
    "current_candidate_side_discrepancy_consumed": True,
    "pa16_inverse_or_fixed_point_solved": False,
    "pa16_repair_applied": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "current_real_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "observed_ci_at_registration": {
        "dedicated": {"run_id": 35558337953, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35558338012, "status": "queued", "conclusion": None},
    },
}

AGENT4_XI_MOMENT_AUDIT = {
    "pr": 950,
    "head": "d366a8aa3b36bb41825e5f5927548cc8bab58eaf",
    "branch": "codex/kokuno-a4-actual-xi-prefix-moments-independent-audit-091",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_actual_xi_prefix_moments_independent_audit.py",
    "source_blob": "b503283bfe7f5e36a4f1da82d4fec559fcf3e171",
    "test_path": "tests/test_constrained_kokuno_a4_actual_xi_prefix_moments_independent_audit.py",
    "test_blob": "6bdb06c8fbeec5a42b1af5091c46cc7f4837a61d",
    "workflow_path": ".github/workflows/kokuno-agent4-actual-xi-prefix-moments-independent-audit.yml",
    "workflow_blob": "3b3af059ac84e3ae6314a8aba25f9e19ca5d6fd0",
    "audited_agent1_pr": 947,
    "audited_agent1_head": "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    "saved_reloaded_public_profile_values_only": True,
    "production_moment_integrator_used_as_reference": False,
    "independent_reference": "piecewise composite Simpson quadrature from public F/U/E values",
    "protocol": {
        "seed": 9173631,
        "fresh_off_grid_eta_count": 10,
        "explicit_eta_probes": [0.0, 1.0e-8, -1.0e-8, 1.0e-6, -1.0e-6],
        "simpson_panels_per_public_segment": [48, 96, 192],
        "relative_rms_gate": 5.0e-3,
        "relative_sampled_max_gate": 2.0e-2,
        "refinement_ratio_gate": 6.0,
        "refinement_floor": 5.0e-9,
        "nontrivial_discrepancy_norm_floor": 1.0e-14,
    },
    "independent_xi_moment_audit_registered": True,
    "independent_xi_moment_audit_admitted": False,
    "complete_ns_residual_audit": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35559078319, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35559078353, "status": "queued", "conclusion": None},
    },
}

AGENT2_SIBLING_STATUS = {
    "pr": 948,
    "head": "fc89770a9ee2ed1897e53ca83f57b995e7b1250d",
    "role": "independent oscillatory Stokes/curl circulation consistency diagnostic",
    "closes_current_xi_moment_seam": False,
    "creates_correction_velocity": False,
    "complete_ns_residual": False,
    "observed_ci_at_registration": {
        "dedicated": {"run_id": 35557903403, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35557903382, "status": "queued", "conclusion": None},
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
    "current_candidate_side_xi_prefix_moment_discrepancy_materialized": True,
    "agent3_current_xi_discrepancy_handoff_registered": True,
    "agent4_independent_xi_moment_audit_registered": True,
    "agent4_independent_xi_moment_audit_admitted": False,
    "upstream_ci_admitted_as_pass": False,
    "source_prepared_appendixA_pressure_materialized": False,
    "source_prepared_appendixA_full_reference_stress_pair_materialized": False,
    "source_admitted_global_kappa0_smallness": False,
    "source_T_sh_lower_bound_verified": False,
    "pa16_repair_applied": False,
    "five_moment_repair_closed": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "current_real_ns_defect_five_moment_discrepancy_materialized": False,
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
    "stage": "candidate-side Xi five-moment discrepancy plus typed correction-lane handoff",
    "input": "A5 #944 registered actual final bridge through Xi=110",
    "new_output": (
        "A1 #947 current Xi moments/discrepancy -> A3 #949 PA.16 row-transform handoff, "
        "paired with A4 #950 independent moment/discrepancy audit protocol"
    ),
    "not_output": (
        "PA.16 repair, complete-NS defect, correction velocity, global Cartesian leading field, "
        "matched pressure/forcing, or full NS validation"
    ),
    "next_shortest_blocker": (
        "independent A4 evidence resolution, then source-governed PA.16/T_sh repair application and "
        "actual inner-to-outer/global Cartesian leading join with matched pressure"
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
        "agent1_xi_prefix_moments": copy.deepcopy(AGENT1_XI_PREFIX_MOMENTS),
        "agent2_sibling_status": copy.deepcopy(AGENT2_SIBLING_STATUS),
        "agent3_xi_discrepancy_bridge": copy.deepcopy(AGENT3_XI_DISCREPANCY_BRIDGE),
        "agent4_xi_moment_audit": copy.deepcopy(AGENT4_XI_MOMENT_AUDIT),
        "xi_moment_discrepancy_ingest_seam": {
            "registration_only": True,
            "current_candidate_side_discrepancy_registered": True,
            "agent3_row_transform_handoff_registered": True,
            "matching_independent_a4_audit_protocol_registered": True,
            "a4_exact_head_pass_required_before_audit_admission": True,
            "queued_or_unresolved_ci_is_pass": False,
            "source_moment_discrepancy_is_complete_ns_defect": False,
            "source_moment_discrepancy_authorizes_correction_stage": False,
            "row_transform_is_pa16_repair": False,
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
        raise ValueError("contract drifted from frozen Agent-5 current-Xi discrepancy registration")
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
