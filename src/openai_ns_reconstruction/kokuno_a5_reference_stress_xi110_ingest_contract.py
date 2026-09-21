"""Agent-5 registration seam for autonomous reference stress through Xi=110.

Integration glue only. This module binds Agent 1 PR #933's executable
repository-autonomous reference stress on 100 <= X <= 110 to Agent 4 PR #936's
implementation-distinct save/reload public-values-only audit.  It does not
materialize the actual candidate-side fixed-kappa F/U/E bridge on 100 < X < 110,
the five-moment repair, a corrected/global Cartesian leading velocity, matched
pressure, restricted forcing, correction velocity, or a complete Navier--Stokes
defect.  The Kokuno route remains a provenance-preserving reconstruction and is
not paper-exact.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-reference-stress-xi110-ingest-v1"
TASK_ID = "KOKUNO-A5-AUTONOMOUS-REFERENCE-STRESS-XI110-INGEST-087"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 929,
    "head": "744a4fbc642b861f61da91de713dd06adc357f09",
    "branch": "codex/kokuno-a5-fixed-kappa-x100-ingest-086",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_fixed_kappa_x100_ingest_contract.py",
    "source_blob": "502e64e97903853a16025815eabe31bf8b2111d0",
    "test_path": "tests/test_kokuno_a5_fixed_kappa_x100_ingest_contract.py",
    "test_blob": "cb6ab8142651b1879bd1a4d597fb7f9e3f6f9266",
    "workflow_path": ".github/workflows/kokuno-agent5-fixed-kappa-x100-ingest.yml",
    "workflow_blob": "57b90219fcc45be891b6a826bd53a0c908ca2a65",
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35549579344, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35549579368, "status": "queued", "conclusion": None},
    },
}

AGENT1_REFERENCE_STRESS_XI110 = {
    "pr": 933,
    "head": "15791891fbbcfdaa614930791a843983afde2b4d",
    "branch": "codex/kokuno-a1-reference-stress-xi110-085",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_autonomous_reference_stress_xi110.py",
    "source_blob": "d534da9efed25df9d281bd0c78fc18ec9b706c09",
    "test_path": "tests/test_kokuno_pa10_autonomous_reference_stress_xi110.py",
    "test_blob": "369f9201303375a446a3d1bf04d019f5c0029d29",
    "workflow_path": ".github/workflows/kokuno-agent1-reference-stress-xi110.yml",
    "workflow_blob": "1144c17afa4d918f213d5477cd51e2c8b8ba4ed5",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "corrected_reader_date": "2026-09-09",
    "reference_range": [100.0, 110.0],
    "X_i": 110.0,
    "quadrature_order": 64,
    "public_values": ["F_r", "U_r", "E_r", "Pi", "p1r", "N_s", "n_s", "p2r"],
    "public_radial_derivatives_registered": True,
    "Xi_handoff_registered": True,
    "semantic_identity_and_save_load_registered": True,
    "pressure_reference_kind": "repository-autonomous admissible Pi0-derived reference stress",
    "source_prepared_appendixA_pressure": False,
    "source_prepared_appendixA_full_reference_stress_pair": False,
    "actual_candidate_fixed_kappa_FUE_bridge_100_to_Xi": False,
    "actual_final_interpolation_100_to_Xi": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35550963448, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35550963482, "status": "queued", "conclusion": None},
    },
}

AGENT4_REFERENCE_STRESS_XI110_AUDIT = {
    "pr": 936,
    "head": "5580526b783500be3701fbddabd6d7c7f90a5a61",
    "branch": "codex/kokuno-a4-reference-stress-xi110-independent-audit-087",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_reference_stress_xi110_independent_audit.py",
    "source_blob": "43b73888c2ceeeb3c406cfc3c884eb0935855118",
    "test_path": "tests/test_constrained_kokuno_a4_reference_stress_xi110_independent_audit.py",
    "test_blob": "b1be7ada5ccfad758692052398efd63e695b8206",
    "workflow_path": ".github/workflows/kokuno-agent4-reference-stress-xi110-independent-audit.yml",
    "workflow_blob": "de2285dfe7575285629f1ff66dc7ba72412c7497",
    "saved_reloaded_public_values_only": True,
    "production_derivatives_used_as_reference": False,
    "independent_operator": "nonuniform degree-six local polynomial interpolation derivative on seven Chebyshev-Lobatto nodes",
    "protocol": {
        "seed": 9173611,
        "random_offgrid_count": 128,
        "total_sample_count": 138,
        "x_range": [100.35, 109.65],
        "physical_halfwidths": [0.24, 0.12, 0.06],
        "derivative_relative_rms_gate": 5.0e-3,
        "derivative_relative_max_gate": 2.0e-2,
        "refinement_ratio_gate": 20.0,
        "refinement_floor": 2.0e-9,
        "algebraic_closure_relative_gate": 5.0e-12,
        "derivative_closure_relative_gate": 5.0e-8,
        "plateau_abs_gate": 2.0e-12,
        "Xi_handoff_abs_gate": 2.0e-12,
    },
    "source_coordinate_prerequisite_audit_only": True,
    "actual_candidate_bridge_audit": False,
    "complete_ns_residual_audit": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35552047780, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35552047829, "status": "queued", "conclusion": None},
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
    "autonomous_reference_stress_to_Xi_registered": True,
    "reference_stress_Xi_handoff_registered": True,
    "independent_reference_stress_Xi_audit_registered": True,
    "upstream_ci_admitted_as_pass": False,
    "source_prepared_appendixA_pressure_materialized": False,
    "source_prepared_appendixA_full_reference_stress_pair_materialized": False,
    "actual_final_interpolation_100_to_Xi_materialized": False,
    "fixed_kappa_FUE_X100_to_Xi_bridge_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
    "actual_ell_i_at_Xi_materialized": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
    "five_moment_repair_closed": False,
    "corrected_global_cartesian_leading_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "real_agent3_correction_velocity_materialized": False,
    "complete_candidate_api_ready": False,
    "same_protocol_full_ns_residual_available": False,
    "current_candidate_eligible_for_full_ns_validation": False,
    "same_protocol_st006_comparison_available_now": False,
    "scientific_admission": False,
    "paper_exact": False,
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
        "agent1_reference_stress_xi110": copy.deepcopy(AGENT1_REFERENCE_STRESS_XI110),
        "agent4_reference_stress_xi110_audit": copy.deepcopy(AGENT4_REFERENCE_STRESS_XI110_AUDIT),
        "reference_stress_xi110_ingest_seam": {
            "registration_only": True,
            "reference_stress_Xi_prerequisite_registered": True,
            "this_is_actual_candidate_fixed_kappa_FUE_bridge": False,
            "semantic_identity_and_save_load_registered": True,
            "saved_reloaded_public_values_only_independent_audit_registered": True,
            "queued_or_unresolved_ci_is_pass": False,
            "profile_or_visual_success_substitutes_for_pde_validation": False,
            "residual_defined_forcing_forbidden": True,
        },
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
        raise ValueError("contract drifted from frozen Agent-5 reference-stress Xi110 ingest registration")
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
