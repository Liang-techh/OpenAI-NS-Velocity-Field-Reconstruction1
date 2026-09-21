"""Agent-5 registration seam for the candidate-side fixed-kappa continuation to X=100.

Integration glue only.  This module binds Agent 1 PR #925's executable
candidate-side fixed-kappa PA.10 continuation to Agent 4 PR #928's
implementation-distinct save/reload audit.  The Agent-1 continuation consumes
repository-autonomous pressure-derived reference stress.  It is not the
source-prepared Appendix-A continuation and it does not materialize the final
X=100->110 bridge, a global Cartesian leading velocity, matched pressure,
restricted forcing, correction velocity, or a complete Navier--Stokes defect.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-fixed-kappa-x100-ingest-v1"
TASK_ID = "KOKUNO-A5-FIXED-KAPPA-X100-INGEST-086"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 922,
    "head": "8bb1139f2b2d9824fd3bbac693ffa4bf2e71dea2",
    "branch": "codex/kokuno-a5-autonomous-reference-stress-ingest-085",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_autonomous_reference_stress_ingest_contract.py",
    "source_blob": "733f76ab40bd48d916301ea9ed2d4c8f9e1eff09",
    "test_path": "tests/test_kokuno_a5_autonomous_reference_stress_ingest_contract.py",
    "test_blob": "2f7d14f70ba86efa1cc7422ac0baa70eb3d65cd1",
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35545886767, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35545886780, "status": "queued", "conclusion": None},
    },
}

AGENT1_FIXED_KAPPA_X100 = {
    "pr": 925,
    "head": "1f3dec711f8ce2052276cd951b2f44d2849579bf",
    "branch": "codex/kokuno-a1-fixed-kappa-x100-084",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_fixed_kappa_x100.py",
    "source_blob": "818d965d82bd92505bd6d3e805e698e1e0da64d8",
    "test_path": "tests/test_kokuno_pa10_fixed_kappa_x100.py",
    "test_blob": "dcf0e4b02000364eb77abb85751a4f76fecb9da7",
    "workflow_path": ".github/workflows/kokuno-agent1-fixed-kappa-x100.yml",
    "workflow_blob": "87352d0a12e72a898f0a3b74bf98dccbd13c65b5",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "corrected_reader_date": "2026-09-09",
    "continuation_range": [8.0, 100.0],
    "selected_kappa0": 0.1,
    "continuation_quadrature_order": 48,
    "public_values": ["F", "U", "E"],
    "public_radial_derivatives": ["F_X", "U_X", "E_X"],
    "public_auxiliary": ["shear"],
    "composition": {
        "F": "activation-endpoint F continued by partial_X log(F)=-kappa_0*p1r/(2X)",
        "U": "activation-endpoint U continued by partial_X U=-kappa_0*nsr/2",
        "E": "sqrt(2X)*F",
    },
    "semantic_identity_and_save_load_registered": True,
    "pressure_reference_kind": "repository-autonomous admissible Pi0-derived reference stress",
    "source_prepared_appendixA_pressure": False,
    "source_exact_fixed_kappa_continuation": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35547826154, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35547826267, "status": "queued", "conclusion": None},
    },
}

AGENT4_FIXED_KAPPA_X100_AUDIT = {
    "pr": 928,
    "head": "602a8bb17d78421cddd493e2a6688a7e4c622e1e",
    "branch": "codex/kokuno-a4-fixed-kappa-x100-independent-audit-086",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_fixed_kappa_x100_independent_audit.py",
    "source_blob": "d6df5ca24c92fe8e2a608dca983719711ca2f85d",
    "test_path": "tests/test_constrained_kokuno_a4_fixed_kappa_x100_independent_audit.py",
    "test_blob": "95710d62c4889f4822dd80eb3b14912a3ca8e667",
    "workflow_path": ".github/workflows/kokuno-agent4-fixed-kappa-x100-independent-audit.yml",
    "workflow_blob": "9693c6347d4904b63690db83e8b51e630e5b190b",
    "saved_reloaded_public_values_only": True,
    "production_derivatives_used_as_reference": False,
    "independent_operator": "centered eighth-order X finite difference reconstructed from public F/U/E values",
    "protocol": {
        "seed": 9173601,
        "random_offgrid_count": 128,
        "probe_families": ["fresh-offgrid", "center/near-center", "endpoint"],
        "relative_x_steps": [1.6e-3, 8.0e-4, 4.0e-4],
        "step_scale": "max(1,abs(X))",
        "derivative_relative_rms_gate": 5.0e-3,
        "derivative_relative_max_gate": 2.0e-2,
        "refinement_ratio_gate": 20.0,
        "refinement_floor": 2.0e-9,
        "value_identity_relative_gate": 5.0e-12,
        "derivative_identity_relative_gate": 5.0e-5,
        "endpoint_composition_abs_gate": 2.0e-12,
    },
    "source_coordinate_profile_audit_only": True,
    "complete_ns_residual_audit": False,
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35548886917, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35548886935, "status": "queued", "conclusion": None},
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
    "candidate_side_fixed_kappa_continuation_to_X100_registered": True,
    "independent_public_values_audit_registered": True,
    "upstream_ci_admitted_as_pass": False,
    "source_prepared_appendixA_pressure_materialized": False,
    "source_exact_fixed_kappa_continuation_materialized": False,
    "final_X100_to_X110_interpolation_materialized": False,
    "Xi110_handoff_materialized": False,
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
        "agent1_fixed_kappa_x100": copy.deepcopy(AGENT1_FIXED_KAPPA_X100),
        "agent4_fixed_kappa_x100_audit": copy.deepcopy(AGENT4_FIXED_KAPPA_X100_AUDIT),
        "fixed_kappa_x100_ingest_seam": {
            "registration_only": True,
            "typed_candidate_side_profile_handoff_registered": True,
            "semantic_identity_and_save_load_registered": True,
            "saved_reloaded_public_values_only_independent_audit_registered": True,
            "queued_or_unresolved_ci_is_pass": False,
            "visual_or_profile_success_substitutes_for_pde_validation": False,
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
        raise ValueError("contract drifted from frozen Agent-5 fixed-kappa X100 ingest registration")
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
