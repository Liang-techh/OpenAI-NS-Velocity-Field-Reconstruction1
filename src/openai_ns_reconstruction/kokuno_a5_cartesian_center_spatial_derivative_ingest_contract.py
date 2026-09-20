"""Fail-closed A5 registration of PA.10 inner Cartesian spatial derivatives.

This seam harvests Agent-1 PR #819's analytic ``velocity_jacobian``,
``divergence`` and ``vorticity`` for the already-registered PA.10 inner
contraction-center ``velocity``/``velocity_dt`` stack.  It is deliberately an
interface/provenance registration only: at this freeze there is no dedicated
Agent-4 independent audit of #819, the parent inner velocity/velocity_dt seams
are not independently admitted on this A5 ancestry, and the field remains
inner-only with no final fixed point, outer/global join, matched pressure,
restricted forcing, full candidate residual or PDE validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_cartesian_center_velocity_dt_ingest_contract import (
    deterministic_cartesian_center_velocity_dt_ingest_contract,
    validate_cartesian_center_velocity_dt_ingest_contract,
)

SCHEMA = "kokuno-a5-cartesian-center-spatial-derivative-ingest-contract-v1"
TASK = "KOKUNO-A5-CARTESIAN-CENTER-SPATIAL-DERIVATIVE-INGEST-CONTRACT-073"
PARENT_A5_PR = 815
PARENT_A5_HEAD = "5fb23be8c30b36736d44274e53304247bc81985e"

AGENT1_SPATIAL_DERIVATIVE_PR = 819
AGENT1_SPATIAL_DERIVATIVE_HEAD = "cb54e6e1a9061cacf78e446dc66bc992f64a0f8d"
AGENT1_SPATIAL_DERIVATIVE_SOURCE_BLOB = "3dbde84d431ba72e847c631543137a5fdfb69459"
AGENT1_SPATIAL_DERIVATIVE_DEDICATED_RUN = 35505168683
AGENT1_SPATIAL_DERIVATIVE_SCHEMA = "kokuno-pa10-cartesian-center-spatial-derivatives-v1"
AGENT1_SPATIAL_DERIVATIVE_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_cartesian_center_spatial_derivatives"
)
AGENT1_SPATIAL_DERIVATIVE_CLASS = "KokunoPA10CartesianCenterSpatialDerivatives"

# #819 exact-head CI was queued at this freeze.  No dedicated Agent-4 audit of
# this new spatial-derivative surface existed in the refreshed open-PR set.
# These are immutable evidence facts for this increment, never caller switches.
AGENT1_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_DEDICATED_SPATIAL_DERIVATIVE_AUDIT_PRESENT = False
AGENT4_INNER_CARTESIAN_CENTER_SPATIAL_DERIVATIVES_INDEPENDENTLY_AUDITED = False

REQUIRED_AGENT1_MEMBERS = (
    "velocity",
    "velocity_dt",
    "coordinate_spatial_derivatives",
    "values",
    "velocity_jacobian",
    "divergence",
    "vorticity",
    "configuration",
    "field_sha256",
    "temporal_derivative_sha256",
    "spatial_derivative_sha256",
    "save_configuration",
    "from_configuration",
    "load_configuration",
    "truth_boundary",
    "report",
    "save_report",
)

SOURCE_SPATIAL_DERIVATIVE_FORMULAS = {
    "coordinates": (
        "tau=1-t; q-z^2 q^(2h)-tau=0; X=(x^2+y^2)/(2q); "
        "eta=z/q^D; L=1-2h eta^2"
    ),
    "coordinate_spatial_derivatives": (
        "q_x=q_y=0; q_z=2 z q^(2h)/L; X_x=x/q; X_y=y/q; "
        "X_z=-X q_z/q; eta_x=eta_y=0; eta_z=(1-eta^2)q^(-D)/L"
    ),
    "jacobian_convention": "J[component,axis]=partial_axis velocity_component",
    "divergence": "trace(velocity_jacobian)",
    "vorticity": "curl(velocity) from the analytic velocity_jacobian",
}

AGENT1_ENGINEERING_PROTOCOL = {
    "production_realization": "analytic-source-chain-rule-v1",
    "velocity_jacobian_vs_fd4_fine_relative_max_gate": 5.0e-5,
    "fd4_fine_vs_coarse_relative_max_gate": 5.0e-5,
    "vorticity_vs_fd4_curl_relative_max_gate": 5.0e-5,
    "analytic_divergence_relative_to_jacobian_scale_gate": 1.0e-9,
    "q_z_vs_centered_fd_relative_max_gate": 2.0e-8,
    "X_z_vs_centered_fd_relative_max_gate": 2.0e-8,
    "eta_z_vs_centered_fd_relative_max_gate": 2.0e-8,
    "rotation_covariance_relative_max_gate": 5.0e-10,
    "axis_jacobian_all_finite_required": True,
    "vorticity_nontrivial_on_probe_required": True,
    "all_probe_values_finite_required": True,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_cartesian_center_spatial_derivative_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_cartesian_center_velocity_dt_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_velocity_dt_ingest_contract(parent)
    parent_status = parent["ingest_status"]

    independently_admitted = bool(
        parent_status["inner_cartesian_center_velocity_dt_ingest_admitted"]
        and AGENT1_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_DEDICATED_SPATIAL_DERIVATIVE_AUDIT_PRESENT
        and AGENT4_INNER_CARTESIAN_CENTER_SPATIAL_DERIVATIVES_INDEPENDENTLY_AUDITED
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "inner_cartesian_center_spatial_derivative_ingest",
        "parent_velocity_dt_contract_sha256": parent["contract_sha256"],
        "agent1_spatial_derivative_binding": {
            "pr": AGENT1_SPATIAL_DERIVATIVE_PR,
            "head": AGENT1_SPATIAL_DERIVATIVE_HEAD,
            "source_blob_sha": AGENT1_SPATIAL_DERIVATIVE_SOURCE_BLOB,
            "dedicated_run": AGENT1_SPATIAL_DERIVATIVE_DEDICATED_RUN,
            "schema": AGENT1_SPATIAL_DERIVATIVE_SCHEMA,
            "module": AGENT1_SPATIAL_DERIVATIVE_MODULE,
            "class": AGENT1_SPATIAL_DERIVATIVE_CLASS,
            "required_members": list(REQUIRED_AGENT1_MEMBERS),
            "source_spatial_derivative_formulas": dict(
                SOURCE_SPATIAL_DERIVATIVE_FORMULAS
            ),
            "engineering_protocol": dict(AGENT1_ENGINEERING_PROTOCOL),
            "public_signatures": {
                "velocity_jacobian": "velocity_jacobian(x,y,z,t)->[...,3,3]",
                "divergence": "divergence(x,y,z,t)->[...]",
                "vorticity": "vorticity(x,y,z,t)->[...,3]",
            },
            "registered_time_interval": [0.25, 0.75],
            "domain_boundary": {
                "inner_source_X_only": True,
                "outside_inner_X_raises": True,
                "outer_global_continuation_materialized": False,
            },
            "identity_boundary": {
                "underlying_field_sha256_preserved": True,
                "temporal_derivative_sha256_preserved": True,
                "spatial_derivative_sha256_separate": True,
                "production_finite_difference_used": False,
            },
        },
        "agent4_spatial_derivative_audit_binding": None,
        "evidence": {
            "parent_inner_cartesian_center_velocity_dt_ingest_admitted": bool(
                parent_status["inner_cartesian_center_velocity_dt_ingest_admitted"]
            ),
            "agent1_exact_head_ci_conclusion": AGENT1_EXACT_HEAD_CI_CONCLUSION,
            "agent4_dedicated_spatial_derivative_audit_present": (
                AGENT4_DEDICATED_SPATIAL_DERIVATIVE_AUDIT_PRESENT
            ),
            "agent4_inner_cartesian_center_spatial_derivatives_independently_audited": (
                AGENT4_INNER_CARTESIAN_CENTER_SPATIAL_DERIVATIVES_INDEPENDENTLY_AUDITED
            ),
        },
        "candidate_api_handoff": dict(parent["candidate_api_handoff"]),
        "differential_operator_handoff": {
            "velocity_jacobian": (
                "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.velocity_jacobian"
            ),
            "divergence": (
                "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.divergence"
            ),
            "vorticity": (
                "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.vorticity"
            ),
            "independent_audit_available": False,
            "usable_for_final_independent_pde_validation": False,
            "reason": (
                "analytic inner-center differential operators are typed and pinned, "
                "but no dedicated Agent-4 audit of #819 exists and the candidate is "
                "still inner-only without pressure/forcing/global assembly"
            ),
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_inner_cartesian_center_spatial_derivative_interface_registered": True,
            "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered": False,
            "inner_cartesian_center_spatial_derivative_ingest_admitted": independently_admitted,
            "inner_cartesian_center_spatial_derivative_identity_registered": True,
            "outer_global_leading_join_materialized": False,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
            "leading_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "interface_registration_not_scientific_admission": True,
            "parent_velocity_dt_admission_invented": False,
            "agent1_spatial_derivative_ci_promoted": False,
            "agent4_spatial_derivative_audit_invented": False,
            "inner_spatial_derivatives_independently_admitted": False,
            "source_center_promoted_to_final_fixed_point": False,
            "outer_global_join_invented": False,
            "global_leading_velocity_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "complete_candidate_api_promoted": False,
            "analytic_inner_jacobian_used_as_final_independent_validator": False,
            "leading_only_ns_residual_invented": False,
            "full_kokuno_residual_invented": False,
            "premature_ST006_comparison_performed": False,
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


def validate_cartesian_center_spatial_derivative_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    parent = deterministic_cartesian_center_velocity_dt_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_velocity_dt_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_velocity_dt_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    a1 = payload.get("agent1_spatial_derivative_binding")
    if not isinstance(a1, Mapping):
        raise ValueError("missing Agent-1 spatial-derivative binding")
    expected_a1 = {
        "pr": AGENT1_SPATIAL_DERIVATIVE_PR,
        "head": AGENT1_SPATIAL_DERIVATIVE_HEAD,
        "source_blob_sha": AGENT1_SPATIAL_DERIVATIVE_SOURCE_BLOB,
        "dedicated_run": AGENT1_SPATIAL_DERIVATIVE_DEDICATED_RUN,
        "schema": AGENT1_SPATIAL_DERIVATIVE_SCHEMA,
        "module": AGENT1_SPATIAL_DERIVATIVE_MODULE,
        "class": AGENT1_SPATIAL_DERIVATIVE_CLASS,
        "registered_time_interval": [0.25, 0.75],
    }
    for key, value in expected_a1.items():
        if a1.get(key) != value:
            raise ValueError(f"Agent-1 spatial-derivative {key} drift")
    if tuple(a1.get("required_members", ())) != REQUIRED_AGENT1_MEMBERS:
        raise ValueError("Agent-1 spatial-derivative required-member drift")
    if a1.get("source_spatial_derivative_formulas") != SOURCE_SPATIAL_DERIVATIVE_FORMULAS:
        raise ValueError("source spatial-derivative formula drift")
    if a1.get("engineering_protocol") != AGENT1_ENGINEERING_PROTOCOL:
        raise ValueError("Agent-1 engineering protocol drift")
    if a1.get("public_signatures") != {
        "velocity_jacobian": "velocity_jacobian(x,y,z,t)->[...,3,3]",
        "divergence": "divergence(x,y,z,t)->[...]",
        "vorticity": "vorticity(x,y,z,t)->[...,3]",
    }:
        raise ValueError("Agent-1 spatial-derivative signature drift")
    if a1.get("domain_boundary") != {
        "inner_source_X_only": True,
        "outside_inner_X_raises": True,
        "outer_global_continuation_materialized": False,
    }:
        raise ValueError("Agent-1 spatial-derivative domain boundary drift")
    if a1.get("identity_boundary") != {
        "underlying_field_sha256_preserved": True,
        "temporal_derivative_sha256_preserved": True,
        "spatial_derivative_sha256_separate": True,
        "production_finite_difference_used": False,
    }:
        raise ValueError("Agent-1 spatial-derivative identity drift")

    if payload.get("agent4_spatial_derivative_audit_binding") is not None:
        raise ValueError("invented Agent-4 spatial-derivative audit binding")
    if payload.get("evidence") != {
        "parent_inner_cartesian_center_velocity_dt_ingest_admitted": False,
        "agent1_exact_head_ci_conclusion": None,
        "agent4_dedicated_spatial_derivative_audit_present": False,
        "agent4_inner_cartesian_center_spatial_derivatives_independently_audited": False,
    }:
        raise ValueError("evidence laundering")

    if payload.get("candidate_api_handoff") != parent["candidate_api_handoff"]:
        raise ValueError("candidate API drift")
    operators = payload.get("differential_operator_handoff")
    if not isinstance(operators, Mapping):
        raise ValueError("missing differential operator handoff")
    if operators.get("velocity_jacobian") != (
        "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.velocity_jacobian"
    ):
        raise ValueError("velocity_jacobian handoff drift")
    if operators.get("divergence") != (
        "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.divergence"
    ):
        raise ValueError("divergence handoff drift")
    if operators.get("vorticity") != (
        "Agent1#819.KokunoPA10CartesianCenterSpatialDerivatives.vorticity"
    ):
        raise ValueError("vorticity handoff drift")
    if operators.get("independent_audit_available") is not False:
        raise ValueError("invented independent audit availability")
    if operators.get("usable_for_final_independent_pde_validation") is not False:
        raise ValueError("premature validator promotion")

    expected_status = {
        **dict(parent["ingest_status"]),
        "typed_inner_cartesian_center_spatial_derivative_interface_registered": True,
        "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered": False,
        "inner_cartesian_center_spatial_derivative_ingest_admitted": False,
        "inner_cartesian_center_spatial_derivative_identity_registered": True,
        "outer_global_leading_join_materialized": False,
        "global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_composite_validated": False,
        "leading_ready": False,
    }
    if payload.get("ingest_status") != expected_status:
        raise ValueError("ingest state laundering")
    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state promotion")
    if payload.get("baseline_vs_kokuno") != parent["baseline_vs_kokuno"]:
        raise ValueError("baseline/comparison drift")
    if payload.get("final_project_gates_unchanged") != parent["final_project_gates_unchanged"]:
        raise ValueError("gate drift")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    if truth.get("interface_registration_not_scientific_admission") is not True:
        raise ValueError("registration/admission boundary lost")
    for key, value in truth.items():
        if key == "interface_registration_not_scientific_admission":
            continue
        if value is not False:
            raise ValueError(f"truth-boundary promotion: {key}")

    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def save_cartesian_center_spatial_derivative_ingest_contract(
    path: str | Path, *, exact_head: str | None = None
) -> dict[str, Any]:
    payload = deterministic_cartesian_center_spatial_derivative_ingest_contract(
        exact_head=exact_head
    )
    validate_cartesian_center_spatial_derivative_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_cartesian_center_spatial_derivative_ingest_contract(
        args.output, exact_head=args.exact_head
    )
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "exact_head": payload["exact_head"],
                "contract_sha256": payload["contract_sha256"],
                "inner_cartesian_center_spatial_derivative_ingest_admitted": payload[
                    "ingest_status"
                ]["inner_cartesian_center_spatial_derivative_ingest_admitted"],
                "independent_spatial_derivative_audit_available": payload[
                    "differential_operator_handoff"
                ]["independent_audit_available"],
                "complete_candidate_api_ready": payload["candidate_api_handoff"][
                    "complete_candidate_api_ready"
                ],
                "leading_ready": payload["stage_state"]["leading_ready"],
                "pde_validated": payload["stage_state"]["pde_validated"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
