"""Fail-closed A5 registration of the PA.10 inner Cartesian center velocity_dt.

This seam harvests Agent-1 PR #811's analytic fixed-Cartesian ``velocity_dt``
for the already-registered Agent-1 #803 inner PA.10 contraction-center velocity
and pins Agent-4 PR #814 as the corresponding implementation-distinct audit
surface.  Agent 4 reconstructs the time derivative from the public ``velocity``
callable with fresh fixed-Cartesian FD4 samples; it does not use Agent-1's
analytic coordinate-time helpers as a numerical oracle.

Registration remains narrower than scientific admission.  The field is still
inner-only, the parent Cartesian-center velocity itself is not independently
admitted on this A5 ancestry, and no final fixed point, outer/global join,
matched pressure, restricted forcing, complete candidate, or held-out momentum
validation is materialized here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_cartesian_center_ingest_contract import (
    deterministic_cartesian_center_ingest_contract,
    validate_cartesian_center_ingest_contract,
)

SCHEMA = "kokuno-a5-cartesian-center-velocity-dt-ingest-contract-v1"
TASK = "KOKUNO-A5-CARTESIAN-CENTER-VELOCITY-DT-INGEST-CONTRACT-072"
PARENT_A5_PR = 807
PARENT_A5_HEAD = "4c1b8b14151affd807fab9d9e3f8cac7d986b224"

AGENT1_VELOCITY_DT_PR = 811
AGENT1_VELOCITY_DT_HEAD = "6f16e837c3a95a6e54af5c2cb7725d094f667026"
AGENT1_VELOCITY_DT_SOURCE_BLOB = "c1d576820416945443bbc24cc931cb64e9e032d9"
AGENT1_VELOCITY_DT_DEDICATED_RUN = 35502463718
AGENT1_VELOCITY_DT_SCHEMA = "kokuno-pa10-cartesian-center-velocity-dt-v1"
AGENT1_VELOCITY_DT_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_cartesian_center_velocity_dt"
)
AGENT1_VELOCITY_DT_CLASS = "KokunoPA10CartesianCenterVelocityTimeDerivative"

AGENT4_VELOCITY_DT_PR = 814
AGENT4_VELOCITY_DT_HEAD = "6f2c9789921885a2c21bb69c642d7f820fdeed41"
AGENT4_VELOCITY_DT_SOURCE_BLOB = "210f445053f027034f3fe20dc7765cd32be1e74d"
AGENT4_VELOCITY_DT_DEDICATED_RUN = 35503374920
AGENT4_VELOCITY_DT_SCHEMA = (
    "kokuno-a4-pa10-cartesian-center-velocity-dt-independent-audit-v1"
)

# Exact-head workflows were queued at this freeze.  They are immutable evidence
# constants in this increment, not caller-controlled switches.
AGENT1_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INNER_CARTESIAN_CENTER_VELOCITY_DT_INDEPENDENTLY_AUDITED = False

REQUIRED_AGENT1_MEMBERS = (
    "velocity",
    "coordinate_time_derivatives",
    "v0_eta",
    "values",
    "velocity_dt",
    "configuration",
    "field_sha256",
    "derivative_sha256",
    "save_configuration",
    "from_configuration",
    "load_configuration",
    "truth_boundary",
    "report",
    "save_report",
)

SOURCE_TIME_DERIVATIVE_FORMULAS = {
    "coordinate_time_derivatives": (
        "L=1-2h eta^2; q_t=-1/L; X_t=X/(qL); eta_t=D eta/(qL)"
    ),
    "fixed_cartesian_derivative": (
        "differentiate u1,u2,u3 at fixed x,y,z using q_t,X_t,eta_t and "
        "analytic F_X,F_eta,U_X,U_eta,(v0)_X,(v0)_eta"
    ),
    "underlying_velocity_identity": (
        "Agent1#811 delegates velocity to Agent1#803 and preserves field_sha256; "
        "derivative_sha256 separately binds the analytic derivative realization"
    ),
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_cartesian_center_velocity_dt_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_cartesian_center_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_cartesian_center_ingest_contract(parent)
    parent_status = parent["ingest_status"]

    independently_admitted = bool(
        parent_status["inner_cartesian_center_velocity_ingest_admitted"]
        and AGENT1_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_INNER_CARTESIAN_CENTER_VELOCITY_DT_INDEPENDENTLY_AUDITED
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "inner_cartesian_center_velocity_dt_ingest",
        "parent_cartesian_center_contract_sha256": parent["contract_sha256"],
        "agent1_velocity_dt_binding": {
            "pr": AGENT1_VELOCITY_DT_PR,
            "head": AGENT1_VELOCITY_DT_HEAD,
            "source_blob_sha": AGENT1_VELOCITY_DT_SOURCE_BLOB,
            "dedicated_run": AGENT1_VELOCITY_DT_DEDICATED_RUN,
            "schema": AGENT1_VELOCITY_DT_SCHEMA,
            "module": AGENT1_VELOCITY_DT_MODULE,
            "class": AGENT1_VELOCITY_DT_CLASS,
            "required_members": list(REQUIRED_AGENT1_MEMBERS),
            "source_time_derivative_formulas": dict(SOURCE_TIME_DERIVATIVE_FORMULAS),
            "public_velocity_dt_signature": "velocity_dt(x,y,z,t)->[...,3]",
            "registered_time_interval": [0.25, 0.75],
            "domain_boundary": {
                "inner_source_X_only": True,
                "outside_inner_X_raises": True,
                "outer_global_continuation_materialized": False,
            },
            "identity_boundary": {
                "underlying_field_sha256_preserved": True,
                "derivative_sha256_separate": True,
                "derivative_realization": "analytic-source-chain-rule-v1",
            },
        },
        "agent4_velocity_dt_binding": {
            "pr": AGENT4_VELOCITY_DT_PR,
            "head": AGENT4_VELOCITY_DT_HEAD,
            "source_blob_sha": AGENT4_VELOCITY_DT_SOURCE_BLOB,
            "dedicated_run": AGENT4_VELOCITY_DT_DEDICATED_RUN,
            "schema": AGENT4_VELOCITY_DT_SCHEMA,
            "independent_numerical_oracle": "public velocity only; fixed-Cartesian centered FD4",
            "seed": 9173471,
            "fresh_random_sample_count": 1024,
            "fd4_time_steps": [4.0e-4, 2.0e-4, 1.0e-4],
            "fine_relative_rms_gate": 5.0e-6,
            "fine_relative_max_gate": 3.0e-5,
            "refinement_ratio_gate": 6.0,
            "refinement_floor": 2.0e-10,
            "velocity_dt_rms_nontriviality_floor": 1.0e-8,
            "exact_axis_transverse_gate": 1.0e-13,
            "axis_near_radius": 1.0e-12,
            "axis_near_required_only_finite": True,
            "mutations": ["velocity_dt_x0.99", "flip_axial_velocity_dt_sign"],
        },
        "evidence": {
            "parent_inner_cartesian_center_velocity_ingest_admitted": bool(
                parent_status["inner_cartesian_center_velocity_ingest_admitted"]
            ),
            "agent1_exact_head_ci_conclusion": AGENT1_EXACT_HEAD_CI_CONCLUSION,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_inner_cartesian_center_velocity_dt_independently_audited": (
                AGENT4_INNER_CARTESIAN_CENTER_VELOCITY_DT_INDEPENDENTLY_AUDITED
            ),
        },
        "candidate_api_handoff": {
            "velocity": parent["candidate_api_handoff"]["velocity"],
            "velocity_dt": (
                "Agent1#811.KokunoPA10CartesianCenterVelocityTimeDerivative.velocity_dt"
            ),
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
            "reason": (
                "inner velocity and analytic velocity_dt are registered, but the parent "
                "inner field is not independently admitted and no final fixed point, "
                "outer/global join, matched pressure, or restricted forcing exists"
            ),
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_inner_cartesian_center_velocity_dt_interface_registered": True,
            "typed_inner_cartesian_center_velocity_dt_independent_audit_registered": True,
            "inner_cartesian_center_velocity_dt_ingest_admitted": independently_admitted,
            "inner_cartesian_center_derivative_identity_registered": True,
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
            "parent_inner_cartesian_center_velocity_admission_invented": False,
            "agent1_velocity_dt_ci_promoted": False,
            "agent4_velocity_dt_audit_promoted": False,
            "inner_cartesian_center_velocity_dt_independently_admitted": False,
            "source_center_promoted_to_final_fixed_point": False,
            "outer_global_join_invented": False,
            "global_leading_velocity_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "complete_candidate_api_promoted": False,
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


def validate_cartesian_center_velocity_dt_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    parent = deterministic_cartesian_center_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_cartesian_center_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_cartesian_center_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    a1 = payload.get("agent1_velocity_dt_binding")
    if not isinstance(a1, Mapping):
        raise ValueError("missing Agent-1 velocity_dt binding")
    expected_a1 = {
        "pr": AGENT1_VELOCITY_DT_PR,
        "head": AGENT1_VELOCITY_DT_HEAD,
        "source_blob_sha": AGENT1_VELOCITY_DT_SOURCE_BLOB,
        "dedicated_run": AGENT1_VELOCITY_DT_DEDICATED_RUN,
        "schema": AGENT1_VELOCITY_DT_SCHEMA,
        "module": AGENT1_VELOCITY_DT_MODULE,
        "class": AGENT1_VELOCITY_DT_CLASS,
        "public_velocity_dt_signature": "velocity_dt(x,y,z,t)->[...,3]",
        "registered_time_interval": [0.25, 0.75],
    }
    for key, value in expected_a1.items():
        if a1.get(key) != value:
            raise ValueError(f"Agent-1 velocity_dt {key} drift")
    if tuple(a1.get("required_members", ())) != REQUIRED_AGENT1_MEMBERS:
        raise ValueError("Agent-1 velocity_dt required-member drift")
    if a1.get("source_time_derivative_formulas") != SOURCE_TIME_DERIVATIVE_FORMULAS:
        raise ValueError("source time-derivative formula drift")
    if a1.get("domain_boundary") != {
        "inner_source_X_only": True,
        "outside_inner_X_raises": True,
        "outer_global_continuation_materialized": False,
    }:
        raise ValueError("Agent-1 velocity_dt domain boundary drift")
    if a1.get("identity_boundary") != {
        "underlying_field_sha256_preserved": True,
        "derivative_sha256_separate": True,
        "derivative_realization": "analytic-source-chain-rule-v1",
    }:
        raise ValueError("Agent-1 derivative identity drift")

    a4 = payload.get("agent4_velocity_dt_binding")
    if not isinstance(a4, Mapping):
        raise ValueError("missing Agent-4 velocity_dt audit binding")
    expected_a4 = {
        "pr": AGENT4_VELOCITY_DT_PR,
        "head": AGENT4_VELOCITY_DT_HEAD,
        "source_blob_sha": AGENT4_VELOCITY_DT_SOURCE_BLOB,
        "dedicated_run": AGENT4_VELOCITY_DT_DEDICATED_RUN,
        "schema": AGENT4_VELOCITY_DT_SCHEMA,
        "independent_numerical_oracle": "public velocity only; fixed-Cartesian centered FD4",
        "seed": 9173471,
        "fresh_random_sample_count": 1024,
        "fd4_time_steps": [4.0e-4, 2.0e-4, 1.0e-4],
        "fine_relative_rms_gate": 5.0e-6,
        "fine_relative_max_gate": 3.0e-5,
        "refinement_ratio_gate": 6.0,
        "refinement_floor": 2.0e-10,
        "velocity_dt_rms_nontriviality_floor": 1.0e-8,
        "exact_axis_transverse_gate": 1.0e-13,
        "axis_near_radius": 1.0e-12,
        "axis_near_required_only_finite": True,
        "mutations": ["velocity_dt_x0.99", "flip_axial_velocity_dt_sign"],
    }
    for key, value in expected_a4.items():
        if a4.get(key) != value:
            raise ValueError(f"Agent-4 velocity_dt {key} drift")

    if payload.get("evidence") != {
        "parent_inner_cartesian_center_velocity_ingest_admitted": False,
        "agent1_exact_head_ci_conclusion": None,
        "agent4_exact_head_ci_conclusion": None,
        "agent4_inner_cartesian_center_velocity_dt_independently_audited": False,
    }:
        raise ValueError("evidence laundering")

    api = payload.get("candidate_api_handoff")
    if not isinstance(api, Mapping):
        raise ValueError("missing candidate API handoff")
    if api.get("velocity") != parent["candidate_api_handoff"]["velocity"]:
        raise ValueError("velocity handoff drift")
    if api.get("velocity_dt") != (
        "Agent1#811.KokunoPA10CartesianCenterVelocityTimeDerivative.velocity_dt"
    ):
        raise ValueError("velocity_dt handoff drift")
    for key in ("pressure", "forcing"):
        if api.get(key) is not None:
            raise ValueError(f"invented candidate API component: {key}")
    if api.get("complete_candidate_api_ready") is not False:
        raise ValueError("premature complete candidate API")

    parent_status = parent["ingest_status"]
    expected_status = {
        **dict(parent_status),
        "typed_inner_cartesian_center_velocity_dt_interface_registered": True,
        "typed_inner_cartesian_center_velocity_dt_independent_audit_registered": True,
        "inner_cartesian_center_velocity_dt_ingest_admitted": False,
        "inner_cartesian_center_derivative_identity_registered": True,
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


def save_cartesian_center_velocity_dt_ingest_contract(
    path: str | Path, *, exact_head: str | None = None
) -> dict[str, Any]:
    payload = deterministic_cartesian_center_velocity_dt_ingest_contract(
        exact_head=exact_head
    )
    validate_cartesian_center_velocity_dt_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_cartesian_center_velocity_dt_ingest_contract(
        args.output, exact_head=args.exact_head
    )
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "exact_head": payload["exact_head"],
                "contract_sha256": payload["contract_sha256"],
                "inner_cartesian_center_velocity_dt_ingest_admitted": payload[
                    "ingest_status"
                ]["inner_cartesian_center_velocity_dt_ingest_admitted"],
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
