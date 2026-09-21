"""Agent-5 integration contract for the current PA.16 source-profile application.

This module is provenance/integration glue only.  It binds the current Agent-1
joined source profile, Agent-3 PA.16 repair receipt/application, Agent-4 audit of
the *underlying Agent-1 joined profile*, and the latest Agent-2 oscillatory
structural lane.  It intentionally does not promote the PA.16 source-moment
repair into a Navier--Stokes correction.

Important truth boundary: Agent-4 #962 does not audit Agent-3 #961 repair
coefficients or the repaired-profile source moments.  Therefore this contract
may record that the source-profile repair application is materialized, but it
must not claim five-moment repair closure, a Cartesian/global leading velocity,
a complete NS defect, or PDE validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-pa16-profile-application-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-PA16-PROFILE-APPLICATION-INGEST-091"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 956,
    "head": "77d29056ac960fa09fdc899ce596dce68f842cd1",
    "branch": "codex/kokuno-a5-current-tsh-certificate-ingest-090",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_tsh_certificate_ingest_contract.py",
    "observed_ci": {
        "dedicated": {"run_id": 35562951875, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35562951840, "status": "queued", "conclusion": None},
    },
}

AGENT1_JOINED_PROFILE = {
    "pr": 959,
    "head": "61b8730fba5623c5acfbd3ec42d54f7cc4f45748",
    "role": "current-lineage joined source-coordinate F/U/E profile through Xi=110",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_joined_profile.py",
    "source_blob": "e38b9b79378aaf97a17226e41da75b477f26ae23",
    "test_path": "tests/test_kokuno_pa16_current_joined_profile.py",
    "test_blob": "41ad4626ab668ec5ccafcfb3dad9c4dcf1d52082",
    "workflow_path": ".github/workflows/kokuno-agent1-current-joined-profile.yml",
    "workflow_blob": "0e67c1dae61fe63db56cf320aa4afda485e0ae61",
    "public_api": ["F", "U", "E", "radial_derivatives"],
    "source_coordinate_profile_only": True,
    "cartesian_global_leading_velocity": False,
    "matched_pressure": False,
    "restricted_forcing": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35564288565, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35564288440, "status": "queued", "conclusion": None},
    },
}

AGENT3_PA16_APPLICATION = {
    "pr": 961,
    "head": "670f7b71d72735d441dcf48f2391c271ee5c163c",
    "role": "bind current PA.16 repair receipt and apply it to the actual joined source profile",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_pa16_profile_application.py",
    "source_blob": "cb3dcd25e0c55a7935daee1b77ffb1741320f83b",
    "test_path": "tests/test_kokuno_current_pa16_profile_application.py",
    "test_blob": "353711260d04f1cb7139d06a6ebccec4a67c19e6",
    "workflow_path": ".github/workflows/kokuno-agent3-current-pa16-profile-application.yml",
    "workflow_blob": "c321a6f3c8ce9847a19de486ca1570ae4bd1949a",
    "audited_agent1_head_expected": "61b8730fba5623c5acfbd3ec42d54f7cc4f45748",
    "coefficient_inf_abs_error": 1.8126505079421038e-13,
    "coefficient_error_is_pde_residual": False,
    "source_profile_application_materialized": True,
    "independently_certifies_post_application_source_moments": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "real_ns_correction_velocity_materialized": False,
    "observed_ci": {
        "dedicated": {"run_id": 35564839899, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35564839903, "status": "queued", "conclusion": None},
    },
}

AGENT4_JOINED_PROFILE_AUDIT = {
    "pr": 962,
    "head": "ac0d099777e9269719302850bcc0c567b139ca88",
    "role": "implementation-distinct audit of A1 #959 public joined profile",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_joined_profile_independent_audit.py",
    "source_blob": "f3c02e1ee28e1b4f4b162be486953d1eca60a37b",
    "test_path": "tests/test_constrained_kokuno_a4_current_joined_profile_independent_audit.py",
    "test_blob": "f4519aac84e56ecefddd266aa85b8ef156a3e6be",
    "workflow_path": ".github/workflows/kokuno-agent4-current-joined-profile-independent-audit.yml",
    "workflow_blob": "631702ec74aebd459fad3cf4c317f912b43920a9",
    "audited_agent1_pr": 959,
    "audited_agent1_head": "61b8730fba5623c5acfbd3ec42d54f7cc4f45748",
    "audits_agent3_repair_coefficients": False,
    "audits_repaired_profile_source_moments": False,
    "protocol": {
        "seed": 9173641,
        "eta_probe_count": 24,
        "segment_fd6_relative_rms_gate": 5.0e-3,
        "segment_fd6_relative_max_gate": 2.0e-2,
        "refinement_ratio_floor": 8.0,
        "join_value_scale_normalized_gate": 5.0e-4,
        "join_derivative_scale_normalized_gate": 5.0e-2,
    },
    "audit_registered": True,
    "audit_admitted": False,
    "complete_ns_residual_audit": False,
    "observed_ci": {
        "dedicated": {"run_id": 35565018042, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35565018024, "status": "queued", "conclusion": None},
    },
}

AGENT2_SIBLING = {
    "pr": 960,
    "head": "6d2fb1f701a34f783dca15a267ae2ce0734ba741",
    "role": "axis-safe batch spatial differentials for frozen public oscillatory velocity",
    "source_path": "src/openai_ns_reconstruction/kokuno_oscillatory_batch_differentials.py",
    "source_blob": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "test_path": "tests/test_constrained_kokuno_oscillatory_batch_differentials.py",
    "test_blob": "df787c86e526d97591dac9be2bc6137f28dfc3f8",
    "workflow_path": ".github/workflows/kokuno-agent2-oscillatory-batch-differentials.yml",
    "workflow_blob": "14fa43d25f076486a97aecfbcf37b05cb8eea8f9",
    "coupled_into_pa16_profile_application": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35564521407, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35564521387, "status": "queued", "conclusion": None},
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
    "current_joined_source_profile_materialized": True,
    "current_pa16_repair_receipt_bound": True,
    "current_pa16_source_profile_application_materialized": True,
    "current_joined_profile_independent_a4_audit_registered": True,
    "current_joined_profile_independent_a4_audit_admitted": False,
    "current_repaired_profile_independent_a4_audit_available": False,
    "post_application_source_moments_independently_verified": False,
    "five_moment_repair_closed": False,
    "source_profile_repair_is_complete_ns_correction": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "actual_inner_to_outer_global_join_materialized": False,
    "corrected_global_cartesian_leading_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "same_protocol_full_ns_residual_available": False,
    "current_candidate_eligible_for_full_ns_validation": False,
    "same_protocol_st006_comparison_available_now": False,
    "upstream_ci_admitted_as_pass": False,
    "scientific_admission": False,
    "paper_exact": False,
}

PIPELINE_POSITION = {
    "stage": "current PA.16 source-profile repair application registration",
    "input": "A1 #959 joined source profile plus A3 #961 bound PA.16 repair receipt",
    "new_output": "checksum-bound registration that the current PA.16 repair is applied to the actual joined source-profile realization",
    "independent_audit_scope": "A4 #962 audits A1 #959 joined F/U/E only; it excludes A3 #961 repair coefficients and repaired-profile moments",
    "not_output": "independently closed five source moments, Cartesian/global leading velocity, matched pressure/forcing, complete NS correction, candidate export, or PDE validation",
    "next_shortest_blocker": "independently audit the repaired profile/source moments, then materialize the actual inner-to-outer/global Cartesian leading join",
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
        "agent1_joined_profile": copy.deepcopy(AGENT1_JOINED_PROFILE),
        "agent3_pa16_application": copy.deepcopy(AGENT3_PA16_APPLICATION),
        "agent4_joined_profile_audit": copy.deepcopy(AGENT4_JOINED_PROFILE_AUDIT),
        "agent2_sibling": copy.deepcopy(AGENT2_SIBLING),
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

    expected = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent1_joined_profile": AGENT1_JOINED_PROFILE,
        "agent3_pa16_application": AGENT3_PA16_APPLICATION,
        "agent4_joined_profile_audit": AGENT4_JOINED_PROFILE_AUDIT,
        "agent2_sibling": AGENT2_SIBLING,
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
    except (TypeError, ValueError):
        errors.append("agent5_exact_head_invalid")

    truth = observed.get("truth_boundary", {})
    readiness = observed.get("readiness", {})
    science = observed.get("frozen_science", {})
    gate = observed.get("final_gate", {})
    a4 = observed.get("agent4_joined_profile_audit", {})
    a3 = observed.get("agent3_pa16_application", {})

    if truth.get("current_repaired_profile_independent_a4_audit_available") is not False:
        errors.append("a4_scope_laundering")
    if a4.get("audits_agent3_repair_coefficients") is not False or a4.get("audits_repaired_profile_source_moments") is not False:
        errors.append("a4_scope_laundering")
    if truth.get("five_moment_repair_closed") is not False or truth.get("post_application_source_moments_independently_verified") is not False:
        errors.append("five_moment_closure_laundering")
    if truth.get("source_profile_repair_is_complete_ns_correction") is not False:
        errors.append("ns_correction_laundering")
    if a3.get("discrepancy_from_complete_ns_defect") is not False or a3.get("authorized_for_gain_gated_ns_stage") is not False:
        errors.append("complete_ns_defect_laundering")
    if science.get("residual_defined_free_forcing_forbidden") is not True:
        errors.append("free_forcing_firewall_weakened")
    if gate != FINAL_GATE:
        errors.append("final_gate_changed")
    if observed.get("st006_baseline") != ST006_BASELINE:
        errors.append("st006_baseline_changed")

    forbidden_truth_promotions = (
        "actual_inner_to_outer_global_join_materialized",
        "corrected_global_cartesian_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "same_protocol_full_ns_residual_available",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
        "upstream_ci_admitted_as_pass",
        "scientific_admission",
        "paper_exact",
    )
    for key in forbidden_truth_promotions:
        if truth.get(key) is not False:
            errors.append(f"forbidden_truth_promotion:{key}")

    if readiness != READINESS:
        errors.append("readiness_promotion_forbidden")
    if truth.get("current_joined_profile_independent_a4_audit_admitted") is not False:
        errors.append("queued_a4_evidence_cannot_be_admitted")
    return sorted(set(errors))


def materialize_receipt(exact_head: str, output: str | Path) -> dict[str, Any]:
    payload = build_contract(exact_head)
    errors = validate_contract(payload)
    if errors:
        raise ValueError("invalid contract: " + ", ".join(errors))
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = materialize_receipt(args.exact_head, args.output)
    print(receipt["contract_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
