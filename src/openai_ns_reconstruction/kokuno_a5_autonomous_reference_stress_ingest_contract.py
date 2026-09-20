"""Agent-5 registration seam for the autonomous-pressure PA.10 reference stress.

Integration glue only.  This module binds Agent 1 PR #918's executable
``N_s / n_s / p_2`` continuation through ``X=100`` to Agent 4 PR #921's
implementation-distinct save/reload audit.  The pressure datum inherited by
#918 is the repository-autonomous admissible ``Pi_0`` realization from #913;
it is not the source-prepared Appendix-A datum and it is not matched/global
pressure.  Consequently this seam cannot promote the Kokuno route to a global
leading candidate, a complete NS residual, or PDE validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-autonomous-reference-stress-ingest-v1"
TASK_ID = "KOKUNO-A5-AUTONOMOUS-REFERENCE-STRESS-INGEST-085"

PARENT_A5 = {
    "pr": 911,
    "head": "dd613c3636391f23a8fc8958af6f17d024482f5f",
    "branch": "codex/kokuno-a5-reference-pressure-increment-ingest-084",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_reference_pressure_increment_ingest_contract.py",
    "source_blob": "c8c23490c62c981440e7e725f8345335df3e59f9",
}

AGENT1_AUTONOMOUS_REFERENCE_STRESS = {
    "pr": 918,
    "head": "c5ccdb8fd790c62cfb760139e5ec108c26d88a81",
    "branch": "codex/kokuno-a1-autonomous-reference-ns-x100-083",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa10_autonomous_reference_ns_x100.py",
    "source_blob": "44159b63c440ee9764261dad41fb3f9d3fec7b0d",
    "test_path": "tests/test_kokuno_pa10_autonomous_reference_ns_x100.py",
    "test_blob": "77eda6075ca185b28700c14a81d18a86b846c94d",
    "workflow_path": ".github/workflows/kokuno-agent1-autonomous-reference-ns-x100.yml",
    "workflow_blob": "242da4be9d8ab9af10901be94f4e7b1edbd7c992",
    "public_source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "public_source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "corrected_reader_date": "2026-09-09",
    "public_surface": "values(X,eta)->N_s_reference_autonomous,n_s_reference_autonomous,p2_reference_autonomous",
    "source_ode": "D_X N_s + N_s = S_n; D_X=X*partial_X; n_s=N_s/L; p_2=X*n_s/E",
    "stress_free_value_handoff": "n_s(X_0,eta)=-2*U_X(X_0,eta)",
    "pressure_provenance": "repository-autonomous admissible Pi0 seed plus source-determined C_p increment",
    "U_eta_realization": {"kind": "centered-FD4", "step": 2.0e-5},
    "Ns_radial_quadrature": {"kind": "Gauss-Legendre", "order": 64},
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35544742368, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35544742466, "status": "queued", "conclusion": None},
    },
}

AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT = {
    "pr": 921,
    "head": "6e6e0e21be5c3a547a2a8db5088eaaf33d63c48b",
    "branch": "codex/kokuno-a4-autonomous-reference-ns-independent-audit-085",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_autonomous_reference_ns_independent_audit.py",
    "source_blob": "4707b4de00d79f7ae9db6531ff625bece4ead2df",
    "test_path": "tests/test_constrained_kokuno_a4_autonomous_reference_ns_independent_audit.py",
    "test_blob": "f4c3b76606907f5b5fd6b82b89dcc5cbffcd9434",
    "workflow_path": ".github/workflows/kokuno-agent4-autonomous-reference-ns-independent-audit.yml",
    "workflow_blob": "ef632d50117370a84f26d672316e3546f46adb74",
    "scientific_input_surface": "serialized/reloaded public values(X,eta) only",
    "independent_operator": "centered Cartesian FD6 radial derivative with scale-aware steps",
    "protocol": {
        "seed": 9173591,
        "random_offgrid_count": 128,
        "random_X_fraction_of_interval": [0.03, 0.97],
        "edge_near_probe_count": 4,
        "exact_X0_probe_count": 5,
        "scale_aware_fd6_steps": [2.0e-3, 1.0e-3, 5.0e-4],
        "step_scale": "max(1,abs(X))",
        "Ns_relative_rms_gate": 5.0e-3,
        "Ns_relative_sampled_max_gate": 2.0e-2,
        "ns_relative_rms_gate": 5.0e-3,
        "ns_relative_sampled_max_gate": 2.0e-2,
        "refinement_ratio_gate": 12.0,
        "refinement_floor": 2.0e-9,
        "Ns_equals_L_ns_and_derivative_closure_relative_gate": 5.0e-10,
        "X0_composition_abs_gate": 2.0e-10,
        "nontriviality_required": True,
    },
    "observed_ci_at_registration": {
        "repository_tests": {"run_id": 35545564129, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35545564264, "status": "queued", "conclusion": None},
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
    """Build the checksum-bound fail-closed registration receipt."""
    contract: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "task_id": TASK_ID,
        "exact_head": str(exact_head),
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_autonomous_reference_stress": {
            **copy.deepcopy(AGENT1_AUTONOMOUS_REFERENCE_STRESS),
            "registered": True,
            "public_source_stress_ode_used": True,
            "autonomous_pressure_reference_stress_pair_callable": True,
            "source_prepared_appendixA_Pi0_materialized": False,
            "source_prepared_reference_nsr_materialized": False,
            "source_prepared_full_reference_stress_pair_materialized": False,
            "matched_global_pressure_materialized": False,
            "global_leading_velocity_materialized": False,
            "scientific_admission": False,
        },
        "agent4_autonomous_reference_stress_audit": {
            **copy.deepcopy(AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT),
            "registered": True,
            "saved_reloaded_public_values_only": True,
            "implementation_distinct_from_agent1_radial_derivative": True,
            "autonomous_reference_stress_independently_audited": False,
            "source_prepared_reference_stress_audited": False,
            "matched_cartesian_pressure_gradient_assessed": False,
            "complete_ns_residual_assessed": False,
            "scientific_admission": False,
        },
        "reference_stress_ingest_seam": {
            "autonomous_reference_stress_registered": True,
            "independent_audit_registered": True,
            "registration_is_scientific_admission": False,
            "queued_or_unresolved_ci_is_pass": False,
            "autonomous_pressure_realization_is_source_prepared_Pi0": False,
            "autonomous_reference_stress_authorized_as_source_prepared_pair": False,
            "autonomous_reference_stress_authorized_as_matched_global_pressure": False,
            "actual_fixed_kappa_continuation_materialized": False,
            "final_X100_to_X110_interpolation_materialized": False,
            "outer_global_leading_velocity_materialized": False,
            "complete_candidate_pressure_api_ready": False,
            "restricted_forcing_api_ready": False,
            "real_agent3_correction_velocity_present": False,
            "complete_candidate_api_ready": False,
            "current_candidate_eligible_for_full_ns_validation": False,
            "same_protocol_st006_comparison_available_now": False,
            "residual_defined_forcing_forbidden": True,
            "next_shortest_dependencies": [
                "close or explicitly preserve the autonomous-vs-source-prepared pressure provenance for the continuation used downstream",
                "materialize the actual fixed-kappa continuation and X=100 to X=110 bridge",
                "assemble the corrected/global Cartesian leading velocity",
                "assemble matched Cartesian pressure/grad-p for that exact global candidate",
                "bind preregistered restricted forcing under the same physical contract",
                "materialize the Agent-3 correction velocity and run the fixed-contract finite cycle",
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
            "source_prepared_pressure_identity_established": False,
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
    """Reject provenance drift and any unsupported global/PDE promotion."""
    c = dict(contract)
    _require(c.get("schema_name") == SCHEMA_NAME, "schema drift")
    _require(c.get("task_id") == TASK_ID, "task drift")
    _require(c.get("parent_a5") == PARENT_A5, "parent Agent-5 provenance drift")

    a1 = c.get("agent1_autonomous_reference_stress", {})
    for key, value in AGENT1_AUTONOMOUS_REFERENCE_STRESS.items():
        _require(a1.get(key) == value, f"Agent-1 reference-stress provenance drift: {key}")
    _require(a1.get("registered") is True, "Agent-1 reference stress must remain registered")
    _require(a1.get("public_source_stress_ode_used") is True, "source stress ODE provenance lost")
    _require(a1.get("autonomous_pressure_reference_stress_pair_callable") is True, "callable autonomous stress pair lost")
    for key in (
        "source_prepared_appendixA_Pi0_materialized",
        "source_prepared_reference_nsr_materialized",
        "source_prepared_full_reference_stress_pair_materialized",
        "matched_global_pressure_materialized",
        "global_leading_velocity_materialized",
        "scientific_admission",
    ):
        _require(a1.get(key) is False, f"unsupported Agent-1 promotion: {key}")

    a4 = c.get("agent4_autonomous_reference_stress_audit", {})
    for key, value in AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT.items():
        _require(a4.get(key) == value, f"Agent-4 reference-stress audit provenance drift: {key}")
    _require(a4.get("registered") is True, "Agent-4 reference-stress audit must remain registered")
    _require(a4.get("saved_reloaded_public_values_only") is True, "A4 public-surface boundary lost")
    _require(a4.get("implementation_distinct_from_agent1_radial_derivative") is True, "A4 independence lost")
    for key in (
        "autonomous_reference_stress_independently_audited",
        "source_prepared_reference_stress_audited",
        "matched_cartesian_pressure_gradient_assessed",
        "complete_ns_residual_assessed",
        "scientific_admission",
    ):
        _require(a4.get(key) is False, f"unsupported Agent-4 audit promotion: {key}")

    seam = c.get("reference_stress_ingest_seam", {})
    _require(seam.get("autonomous_reference_stress_registered") is True, "reference-stress registration missing")
    _require(seam.get("independent_audit_registered") is True, "independent audit registration missing")
    _require(seam.get("registration_is_scientific_admission") is False, "registration cannot imply admission")
    _require(seam.get("queued_or_unresolved_ci_is_pass") is False, "queued CI cannot be PASS")
    for key in (
        "autonomous_pressure_realization_is_source_prepared_Pi0",
        "autonomous_reference_stress_authorized_as_source_prepared_pair",
        "autonomous_reference_stress_authorized_as_matched_global_pressure",
        "actual_fixed_kappa_continuation_materialized",
        "final_X100_to_X110_interpolation_materialized",
        "outer_global_leading_velocity_materialized",
        "complete_candidate_pressure_api_ready",
        "restricted_forcing_api_ready",
        "real_agent3_correction_velocity_present",
        "complete_candidate_api_ready",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
    ):
        _require(seam.get(key) is False, f"unsupported reference-stress seam promotion: {key}")
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
        "source_prepared_pressure_identity_established",
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
