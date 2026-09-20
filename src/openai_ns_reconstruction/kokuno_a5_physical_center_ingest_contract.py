"""Fail-closed A5 registration of Agent-1's PA.10 physical-center map.

This contract records the newly executable source formulas

    Y = Lambda X,
    F = (phi_*/C) Phi,
    U = U_* + Lambda^-1 u,

without promoting the PA.10 contraction center to the final corrected leading
profile.  Exact-head CI and Agent-4 independent audit outcomes are deliberately
kept separate from interface registration.  In particular, the source complex
normalization condition on ``C`` remains a blocking admission condition.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_candidate_artifact_contract import (
    FINAL_PROJECT_GATES,
    PIPELINE_ORDER,
    ST006_BASELINE,
)
from .kokuno_a5_phi0_profile_ingest_contract import (
    deterministic_phi0_ingest_contract,
    validate_phi0_ingest_contract,
)

SCHEMA = "kokuno-a5-physical-center-ingest-contract-v1"
TASK = "KOKUNO-A5-PHYSICAL-CENTER-INGEST-CONTRACT-069"
PARENT_A5_PR = 782
PARENT_A5_HEAD = "1d7055bbbca528f8c542a65fbf71aeff84b6af43"

AGENT1_PHYSICAL_PR = 787
AGENT1_PHYSICAL_HEAD = "dbde743d6f71146a3c591968cf60a0abe449d310"
AGENT1_PHYSICAL_SOURCE_BLOB = "f24246ce84ccb8c1c6a7b1eba31865c416f2da80"
AGENT1_PHYSICAL_DEDICATED_RUN = 35494388163
AGENT1_PHYSICAL_SCHEMA = "kokuno-pa10-physical-center-profile-contract-v1"
AGENT1_PHYSICAL_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_physical_center_profile_contract"
)
AGENT1_PHYSICAL_CLASS = "KokunoPA10PhysicalCenterProfileContract"

AGENT4_PHYSICAL_AUDIT_PR = 790
AGENT4_PHYSICAL_AUDIT_HEAD = "cfd4ca4bb49940a9754f419404730b74fc1e7623"
AGENT4_PHYSICAL_AUDIT_SOURCE_BLOB = "1d1851d07db42834b594c55c8185f66d9de010fc"
AGENT4_PHYSICAL_AUDIT_DEDICATED_RUN = 35495064305
AGENT4_PHYSICAL_AUDIT_SCHEMA = (
    "kokuno-agent4-pa10-physical-center-independent-audit-v1"
)

# Exact-head sibling workflows are queued at this freeze.  None of these
# evidence values are caller-controlled.
AGENT1_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_MAPPING_FORMULA_INDEPENDENTLY_CONSISTENT: bool | None = None
AGENT4_CONFIGURED_C_REAL_AXIS_OBSTRUCTION_DETECTED: bool | None = None

# These are known truth-boundary facts of the pinned #787/#790 source heads.
AGENT4_PHYSICAL_CENTER_PROFILE_INDEPENDENTLY_ADMITTED = False
SOURCE_COMPLEX_C_NORMALIZATION_CERTIFIED = False

REQUIRED_METHODS = (
    "zeta_primitive",
    "phi_star",
    "values",
    "derivatives",
    "incompressibility_defect",
    "configuration",
    "from_configuration",
    "save_configuration",
    "load_configuration",
    "report",
    "save_report",
)
REQUIRED_VALUE_FIELDS = (
    "X",
    "Y",
    "eta",
    "Lambda",
    "C",
    "phi_star",
    "g",
    "Phi_0",
    "u_0",
    "F_0",
    "E_0",
    "U_0",
    "M_0",
    "M_0_over_X",
    "v_0",
    "V_0",
)
REQUIRED_DERIVATIVE_FIELDS = (
    "phi_star_eta",
    "g_eta",
    "Phi_0_X",
    "Phi_0_eta",
    "F_0_X",
    "F_0_eta",
    "U_0_X",
    "U_0_eta",
    "v_0_X",
)
SOURCE_FORMULAS = {
    "source_radial_mapping": "Y=Lambda X",
    "physical_swirl_profile": "F=(phi_*/C)Phi",
    "physical_axial_profile": "U=U_*+Lambda^-1 u",
    "azimuthal_profile": "E=sqrt(2X)F",
    "radial_profile": "V_0=X v_0",
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_physical_center_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_phi0_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_phi0_ingest_contract(parent)
    if PIPELINE_ORDER[0] != "source_profile_ingest":
        raise RuntimeError("pipeline drift")

    parent_status = parent["ingest_status"]
    independently_admitted = bool(
        parent_status["phi0_profile_ingest_admitted"]
        and AGENT1_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_MAPPING_FORMULA_INDEPENDENTLY_CONSISTENT is True
        and AGENT4_CONFIGURED_C_REAL_AXIS_OBSTRUCTION_DETECTED is False
        and SOURCE_COMPLEX_C_NORMALIZATION_CERTIFIED
        and AGENT4_PHYSICAL_CENTER_PROFILE_INDEPENDENTLY_ADMITTED
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "source_physical_center_profile_ingest",
        "pipeline_order": list(PIPELINE_ORDER),
        "parent_phi0_contract_sha256": parent["contract_sha256"],
        "agent1_physical_center_binding": {
            "pr": AGENT1_PHYSICAL_PR,
            "head": AGENT1_PHYSICAL_HEAD,
            "source_blob_sha": AGENT1_PHYSICAL_SOURCE_BLOB,
            "dedicated_run": AGENT1_PHYSICAL_DEDICATED_RUN,
            "schema": AGENT1_PHYSICAL_SCHEMA,
            "module": AGENT1_PHYSICAL_MODULE,
            "class": AGENT1_PHYSICAL_CLASS,
            "required_methods": list(REQUIRED_METHODS),
            "required_value_fields": list(REQUIRED_VALUE_FIELDS),
            "required_derivative_fields": list(REQUIRED_DERIVATIVE_FIELDS),
            "source_formulas": dict(SOURCE_FORMULAS),
            "scope": {
                "contraction_center_only": True,
                "final_corrected_fixed_point": False,
                "cartesian_spacetime_velocity": False,
                "matched_global_pressure": False,
            },
        },
        "agent4_independent_audit_binding": {
            "pr": AGENT4_PHYSICAL_AUDIT_PR,
            "head": AGENT4_PHYSICAL_AUDIT_HEAD,
            "source_blob_sha": AGENT4_PHYSICAL_AUDIT_SOURCE_BLOB,
            "dedicated_run": AGENT4_PHYSICAL_AUDIT_DEDICATED_RUN,
            "schema": AGENT4_PHYSICAL_AUDIT_SCHEMA,
            "independent_paths": [
                "fresh off-grid source-profile mapping replay",
                "composite Simpson real-axis phi_* primitive",
                "FD4 public-value derivative replay",
                "FD4 center incompressibility replay",
            ],
        },
        "evidence": {
            "agent1_exact_head_ci_conclusion": AGENT1_EXACT_HEAD_CI_CONCLUSION,
            "agent4_audit_present": True,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_mapping_formula_independently_consistent": (
                AGENT4_MAPPING_FORMULA_INDEPENDENTLY_CONSISTENT
            ),
            "agent4_configured_C_real_axis_obstruction_detected": (
                AGENT4_CONFIGURED_C_REAL_AXIS_OBSTRUCTION_DETECTED
            ),
            "agent4_physical_center_profile_independently_admitted": (
                AGENT4_PHYSICAL_CENTER_PROFILE_INDEPENDENTLY_ADMITTED
            ),
            "source_complex_C_normalization_certified": (
                SOURCE_COMPLEX_C_NORMALIZATION_CERTIFIED
            ),
        },
        "mapping_registration": {
            "source_Y_equals_Lambda_X_formula_registered": True,
            "physical_F_from_Phi_formula_registered": True,
            "physical_U_from_u_formula_registered": True,
            "repository_X_equals_source_Y_asserted": False,
            "mapping_independently_admitted": independently_admitted,
        },
        "ingest_status": {
            "typed_axis_profile_interface_registered": parent_status[
                "typed_axis_profile_interface_registered"
            ],
            "axis_profile_ingest_admitted": parent_status[
                "axis_profile_ingest_admitted"
            ],
            "typed_phi0_profile_interface_registered": parent_status[
                "typed_phi0_profile_interface_registered"
            ],
            "phi0_profile_ingest_admitted": parent_status[
                "phi0_profile_ingest_admitted"
            ],
            "typed_physical_center_profile_interface_registered": True,
            "physical_center_profile_ingest_admitted": independently_admitted,
            "source_Y_equals_Lambda_X_formula_registered": True,
            "physical_F_from_Phi_formula_registered": True,
            "physical_U_from_u_formula_registered": True,
            "source_complex_C_normalization_certified": False,
            "fixed_point_correction_materialized": False,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
            "leading_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": {
            "retained_repository_baseline": dict(ST006_BASELINE),
            "kokuno_same_protocol_full_candidate_residual_available": False,
            "comparison_performed": False,
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "interface_registration_not_scientific_admission": True,
            "queued_sibling_ci_promoted": False,
            "agent4_unknown_result_invented": False,
            "repository_X_equals_source_Y_asserted": False,
            "source_complex_C_normalization_claimed": False,
            "contraction_center_promoted_to_final_fixed_point": False,
            "source_native_center_promoted_to_cartesian_velocity": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "kokuno_replay_used_as_final_independent_validation": False,
            "leading_ready": False,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
    }
    payload["contract_sha256"] = _digest(payload)
    return payload


def validate_physical_center_ingest_contract(payload: Mapping[str, Any]) -> None:
    parent = deterministic_phi0_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_phi0_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_phi0_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    binding = payload.get("agent1_physical_center_binding")
    if not isinstance(binding, Mapping):
        raise ValueError("missing Agent-1 physical-center binding")
    expected_agent1 = {
        "head": AGENT1_PHYSICAL_HEAD,
        "source_blob_sha": AGENT1_PHYSICAL_SOURCE_BLOB,
        "schema": AGENT1_PHYSICAL_SCHEMA,
        "module": AGENT1_PHYSICAL_MODULE,
        "class": AGENT1_PHYSICAL_CLASS,
    }
    for key, value in expected_agent1.items():
        if binding.get(key) != value:
            raise ValueError(f"Agent-1 physical-center {key} drift")
    if binding.get("source_formulas") != SOURCE_FORMULAS:
        raise ValueError("source formula drift")

    audit = payload.get("agent4_independent_audit_binding")
    if not isinstance(audit, Mapping):
        raise ValueError("missing Agent-4 audit binding")
    for key, value in {
        "head": AGENT4_PHYSICAL_AUDIT_HEAD,
        "source_blob_sha": AGENT4_PHYSICAL_AUDIT_SOURCE_BLOB,
        "schema": AGENT4_PHYSICAL_AUDIT_SCHEMA,
    }.items():
        if audit.get(key) != value:
            raise ValueError(f"Agent-4 physical-center {key} drift")

    expected_evidence = {
        "agent1_exact_head_ci_conclusion": None,
        "agent4_audit_present": True,
        "agent4_exact_head_ci_conclusion": None,
        "agent4_mapping_formula_independently_consistent": None,
        "agent4_configured_C_real_axis_obstruction_detected": None,
        "agent4_physical_center_profile_independently_admitted": False,
        "source_complex_C_normalization_certified": False,
    }
    if payload.get("evidence") != expected_evidence:
        raise ValueError("evidence laundering")

    mapping = payload.get("mapping_registration")
    if mapping != {
        "source_Y_equals_Lambda_X_formula_registered": True,
        "physical_F_from_Phi_formula_registered": True,
        "physical_U_from_u_formula_registered": True,
        "repository_X_equals_source_Y_asserted": False,
        "mapping_independently_admitted": False,
    }:
        raise ValueError("mapping state laundering")

    parent_status = parent["ingest_status"]
    expected_status = {
        "typed_axis_profile_interface_registered": parent_status[
            "typed_axis_profile_interface_registered"
        ],
        "axis_profile_ingest_admitted": parent_status["axis_profile_ingest_admitted"],
        "typed_phi0_profile_interface_registered": parent_status[
            "typed_phi0_profile_interface_registered"
        ],
        "phi0_profile_ingest_admitted": parent_status[
            "phi0_profile_ingest_admitted"
        ],
        "typed_physical_center_profile_interface_registered": True,
        "physical_center_profile_ingest_admitted": False,
        "source_Y_equals_Lambda_X_formula_registered": True,
        "physical_F_from_Phi_formula_registered": True,
        "physical_U_from_u_formula_registered": True,
        "source_complex_C_normalization_certified": False,
        "fixed_point_correction_materialized": False,
        "global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_composite_validated": False,
        "leading_ready": False,
    }
    if payload.get("ingest_status") != expected_status:
        raise ValueError("ingest state laundering")
    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state promotion")
    if payload.get("final_project_gates_unchanged") != FINAL_PROJECT_GATES:
        raise ValueError("gate drift")

    baseline = payload.get("baseline_vs_kokuno")
    if not isinstance(baseline, Mapping):
        raise ValueError("missing baseline comparison state")
    if baseline.get("retained_repository_baseline") != ST006_BASELINE:
        raise ValueError("ST006 baseline drift")
    if baseline.get("kokuno_same_protocol_full_candidate_residual_available") is not False:
        raise ValueError("fabricated full-candidate residual")
    if baseline.get("comparison_performed") is not False:
        raise ValueError("premature ST006 comparison")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    if truth.get("interface_registration_not_scientific_admission") is not True:
        raise ValueError("registration/admission boundary lost")
    for key in (
        "queued_sibling_ci_promoted",
        "agent4_unknown_result_invented",
        "repository_X_equals_source_Y_asserted",
        "source_complex_C_normalization_claimed",
        "contraction_center_promoted_to_final_fixed_point",
        "source_native_center_promoted_to_cartesian_velocity",
        "matched_global_pressure_invented",
        "restricted_forcing_composite_invented",
        "free_residual_defined_forcing_allowed",
        "threshold_relaxed",
        "kokuno_replay_used_as_final_independent_validation",
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError("truth-boundary promotion")

    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def save_contract(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_physical_center_ingest_contract(exact_head=exact_head)
    validate_physical_center_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = save_contract(args.output, exact_head=args.exact_head)
    print(
        json.dumps(
            {
                "exact_head": payload["exact_head"],
                "physical_center_profile_ingest_admitted": payload["ingest_status"][
                    "physical_center_profile_ingest_admitted"
                ],
                "leading_ready": payload["stage_state"]["leading_ready"],
                "contract_sha256": payload["contract_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
