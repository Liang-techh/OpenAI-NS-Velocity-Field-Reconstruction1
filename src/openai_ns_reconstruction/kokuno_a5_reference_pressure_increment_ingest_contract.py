"""Agent-5 registration seam for the first Kokuno pressure-bearing artifact.

Integration glue only.  This module binds Agent 1 PR #907's source-determined
reference pressure *increment* ``C_p = Pi_r - Pi_0`` together with Agent 4
PR #910's implementation-distinct derivative audit.  The absolute axis datum
``Pi_0(eta)`` is still missing, so this seam must not be promoted to matched
pressure, Cartesian ``grad p``, a complete NS residual, or PDE validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-reference-pressure-increment-ingest-v1"
TASK_ID = "KOKUNO-A5-REFERENCE-PRESSURE-INCREMENT-INGEST-084"

PARENT_A5 = {
    "pr": 904,
    "head": "94fe65a1b93b0a9b67e4736015ed93034c2a88df",
    "branch": "codex/kokuno-a5-canonical-quadrature-evidence-ingest-083",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_canonical_quadrature_evidence_ingest_contract.py",
    "source_blob": "14b4352ce4248cf0d6b197835aa7efad2edbd906",
}

AGENT1_PRESSURE_INCREMENT = {
    "pr": 907,
    "head": "3f43931cec1fc762678a5a8134adfceaf7edfa81",
    "branch": "codex/kokuno-a1-reference-pressure-increment-x100-081",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_reference_pressure_increment_x100.py",
    "source_blob": "4704fed4271da2c6070110a6935d2f33f7932f56",
    "workflow_path": ".github/workflows/kokuno-agent1-reference-pressure-increment-x100.yml",
    "workflow_blob": "b88fb38ef162097e9b471e67e05566dd3e2ad61a",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "corrected_reader_date": "2026-09-09",
    "surface": "C_p_reference(X,eta)=Pi_r(X,eta)-Pi_0(eta)",
    "source_identity": "partial_X C_p_reference = F_r^2 = E_r^2/(2X)",
    "radial_quadrature": {"kind": "Gauss-Legendre", "order": 96},
    "eta_derivative_realization": {"kind": "centered-FD4", "step": 2.0e-5},
    "observed_ci": {
        "dedicated": {"run_id": 35539131381, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35539131366, "status": "queued", "conclusion": None},
    },
}

AGENT4_PRESSURE_INCREMENT_AUDIT = {
    "pr": 910,
    "head": "68a9ee155f3e4e7f50ff625778d3fbe21d4021cb",
    "branch": "codex/kokuno-a4-reference-pressure-increment-audit-083",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_reference_pressure_increment_independent_audit.py",
    "source_blob": "027abf88fa1a916fcde3c6d818b19a9ff143c445",
    "test_path": "tests/test_constrained_kokuno_a4_reference_pressure_increment_independent_audit.py",
    "test_blob": "59595ca760e136649b958e96fecb48633e23de90",
    "dedicated_workflow_present": False,
    "scientific_input_surface": "serialized/reloaded public pressure_increment(X,eta) only",
    "independent_operator": "centered two-point + Richardson extrapolation D_R=(4D(h/2)-D(h))/3",
    "protocol": {
        "seed": 9173571,
        "random_offgrid_count": 128,
        "X_range": [0.08, 99.0],
        "eta_range": [-0.70, 0.70],
        "near_axis_probe_count": 3,
        "exact_axis_probe_count": 5,
        "X_steps": [8.0e-3, 4.0e-3, 2.0e-3],
        "eta_steps": [8.0e-4, 4.0e-4, 2.0e-4],
        "axis_increment_abs_gate": 1.0e-12,
        "radial_relative_rms_gate": 2.0e-4,
        "radial_relative_sampled_max_gate": 1.0e-3,
        "eta_relative_rms_gate": 2.0e-3,
        "eta_relative_sampled_max_gate": 1.0e-2,
        "refinement_ratio_gate": 6.0,
        "refinement_floor": 1.0e-10,
        "pressure_increment_rms_floor": 1.0e-10,
        "radial_gradient_rms_floor": 1.0e-10,
    },
    "observed_ci": {
        "repository_tests": {"run_id": 35539845126, "status": "queued", "conclusion": None},
    },
}

FINAL_GATE = {
    "momentum_sampled_max": 1.0e-3,
    "momentum_volume_l2": 1.0e-3,
    "divergence_sampled_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}

ST006_BASELINE = {
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


def _canonical(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(dict(payload), sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical(payload)).hexdigest()


def _without_sha(contract: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(contract))
    payload.pop("contract_sha256", None)
    return payload


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def build_contract(*, exact_head: str) -> dict[str, Any]:
    """Build the fail-closed Agent-5 registration receipt for this exact head."""
    contract: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "task_id": TASK_ID,
        "exact_head": str(exact_head),
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_pressure_increment": {
            **copy.deepcopy(AGENT1_PRESSURE_INCREMENT),
            "registered": True,
            "source_determined_pressure_increment": True,
            "absolute_axis_pressure_Pi0_materialized": False,
            "absolute_reference_pressure_materialized": False,
            "matched_global_pressure_materialized": False,
            "cartesian_pressure_gradient_materialized": False,
            "scientific_admission": False,
        },
        "agent4_pressure_increment_audit": {
            **copy.deepcopy(AGENT4_PRESSURE_INCREMENT_AUDIT),
            "registered": True,
            "implementation_distinct_from_agent1_production_derivatives": True,
            "saved_reloaded_public_surface_only": True,
            "reference_pressure_increment_independently_audited": False,
            "absolute_pressure_audited": False,
            "cartesian_pressure_gradient_assessed": False,
            "complete_ns_residual_assessed": False,
            "scientific_admission": False,
        },
        "pressure_ingest_seam": {
            "reference_pressure_increment_registered": True,
            "independent_audit_registered": True,
            "registration_is_scientific_admission": False,
            "queued_or_unresolved_ci_is_pass": False,
            "pressure_increment_is_absolute_pressure": False,
            "pressure_increment_authorized_as_matched_global_pressure": False,
            "pressure_increment_authorized_as_cartesian_grad_p": False,
            "Pi0_required_before_absolute_pressure": True,
            "complete_candidate_pressure_api_ready": False,
            "complete_candidate_api_ready": False,
            "current_candidate_eligible_for_full_ns_validation": False,
            "same_protocol_st006_comparison_available_now": False,
            "residual_defined_forcing_forbidden": True,
            "next_pressure_dependencies": [
                "materialize source-consistent Pi_0(eta) axis datum",
                "assemble absolute reference pressure",
                "materialize pressure-bearing n_{s,r} / complete reference stress pair",
                "complete actual fixed-kappa continuation and global leading join",
                "assemble matched Cartesian pressure/grad p on the global candidate",
            ],
        },
        "readiness": copy.deepcopy(READINESS),
        "final_gate": copy.deepcopy(FINAL_GATE),
        "baseline": {"st006": copy.deepcopy(ST006_BASELINE)},
        "truth_boundary": {
            "kokuno_reconstruction_provenance_bounded": True,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "matched_pressure_present": False,
            "restricted_forcing_present": False,
            "real_correction_velocity_present": False,
            "complete_full_ns_candidate_present": False,
            "heldout_complete_ns_residual_assessed": False,
            "canonical_quadrature_receipt_available": False,
            "same_protocol_st006_improvement_claimed": False,
            "pde_validated": False,
        },
    }
    contract["contract_sha256"] = _sha256(contract)
    validate_contract(contract)
    return contract


def validate_contract(contract: Mapping[str, Any]) -> None:
    """Reject provenance drift or any pressure/PDE promotion beyond current evidence."""
    c = dict(contract)
    _require(c.get("schema_name") == SCHEMA_NAME, "schema drift")
    _require(c.get("task_id") == TASK_ID, "task drift")
    _require(c.get("parent_a5") == PARENT_A5, "parent Agent-5 provenance drift")

    a1 = c.get("agent1_pressure_increment", {})
    for key, value in AGENT1_PRESSURE_INCREMENT.items():
        _require(a1.get(key) == value, f"Agent-1 pressure provenance drift: {key}")
    _require(a1.get("registered") is True, "Agent-1 pressure increment must remain registered")
    for key in (
        "absolute_axis_pressure_Pi0_materialized",
        "absolute_reference_pressure_materialized",
        "matched_global_pressure_materialized",
        "cartesian_pressure_gradient_materialized",
        "scientific_admission",
    ):
        _require(a1.get(key) is False, f"unsupported Agent-1 pressure promotion: {key}")

    a4 = c.get("agent4_pressure_increment_audit", {})
    for key, value in AGENT4_PRESSURE_INCREMENT_AUDIT.items():
        _require(a4.get(key) == value, f"Agent-4 audit provenance drift: {key}")
    _require(a4.get("registered") is True, "Agent-4 audit must remain registered")
    _require(a4.get("implementation_distinct_from_agent1_production_derivatives") is True, "A4 independence lost")
    _require(a4.get("saved_reloaded_public_surface_only") is True, "A4 public-surface boundary lost")
    for key in (
        "reference_pressure_increment_independently_audited",
        "absolute_pressure_audited",
        "cartesian_pressure_gradient_assessed",
        "complete_ns_residual_assessed",
        "scientific_admission",
    ):
        _require(a4.get(key) is False, f"unsupported Agent-4 audit promotion: {key}")

    seam = c.get("pressure_ingest_seam", {})
    _require(seam.get("reference_pressure_increment_registered") is True, "pressure increment registration missing")
    _require(seam.get("independent_audit_registered") is True, "independent audit registration missing")
    _require(seam.get("registration_is_scientific_admission") is False, "registration cannot imply admission")
    _require(seam.get("queued_or_unresolved_ci_is_pass") is False, "queued CI cannot be PASS")
    _require(seam.get("pressure_increment_is_absolute_pressure") is False, "C_p increment cannot be absolute pressure")
    _require(seam.get("pressure_increment_authorized_as_matched_global_pressure") is False, "C_p increment cannot be matched pressure")
    _require(seam.get("pressure_increment_authorized_as_cartesian_grad_p") is False, "C_p increment cannot be Cartesian grad p")
    _require(seam.get("Pi0_required_before_absolute_pressure") is True, "Pi_0 dependency cannot be dropped")
    _require(seam.get("complete_candidate_pressure_api_ready") is False, "pressure API is not complete")
    _require(seam.get("complete_candidate_api_ready") is False, "candidate API is not complete")
    _require(seam.get("current_candidate_eligible_for_full_ns_validation") is False, "candidate is not validator-eligible")
    _require(seam.get("same_protocol_st006_comparison_available_now") is False, "ST006 comparison is not available")
    _require(seam.get("residual_defined_forcing_forbidden") is True, "free residual-defined forcing remains forbidden")

    _require(c.get("readiness") == READINESS, "readiness promotion/drift")
    _require(c.get("final_gate") == FINAL_GATE, "final scientific gate drift")
    _require(c.get("baseline") == {"st006": ST006_BASELINE}, "ST006 baseline drift")

    truth = c.get("truth_boundary", {})
    _require(truth.get("kokuno_reconstruction_provenance_bounded") is True, "provenance boundary lost")
    for key in (
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "matched_pressure_present",
        "restricted_forcing_present",
        "real_correction_velocity_present",
        "complete_full_ns_candidate_present",
        "heldout_complete_ns_residual_assessed",
        "canonical_quadrature_receipt_available",
        "same_protocol_st006_improvement_claimed",
        "pde_validated",
    ):
        _require(truth.get(key) is False, f"unsupported truth promotion: {key}")

    expected_sha = _sha256(_without_sha(c))
    _require(c.get("contract_sha256") == expected_sha, "contract checksum mismatch")


def write_contract(path: str | Path, *, exact_head: str) -> dict[str, Any]:
    contract = build_contract(exact_head=exact_head)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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
