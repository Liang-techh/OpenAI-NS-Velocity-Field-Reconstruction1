"""Fail-closed A5 registration of the strict-inner leading+oscillatory transport seam.

Agent 2 #840 composes the currently executable PA.10 *inner contraction-center*
velocity with the frozen oscillatory field into

    T = d_t u + (u . grad)u - nu Delta u,  nu = 0.01,

for ``u = u_inner + u_osc``.  This A5 layer binds that typed handoff into the
integration pipeline without copying Agent-1/2 mathematics and without calling
it a Navier--Stokes residual: pressure, restricted forcing, the corrected/global
leading field, and the real correction velocity are all still absent.

The seam is useful because it closes the pressure/forcing-free transport block
needed by a future actual-defect operator.  Registration is not scientific
admission.  Exact-head CI is unresolved at this integration cut and no dedicated
Agent-4 audit of #840 exists yet.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_cartesian_center_spatial_derivative_audit_ingest_contract import (
    deterministic_cartesian_center_spatial_derivative_audit_ingest_contract,
    validate_cartesian_center_spatial_derivative_audit_ingest_contract,
)

SCHEMA = "kokuno-a5-inner-leading-oscillatory-transport-ingest-contract-v1"
TASK = "KOKUNO-A5-INNER-LEADING-OSCILLATORY-TRANSPORT-INGEST-075"
PARENT_A5_PR = 833
PARENT_A5_HEAD = "c4d29ea5987e297977676bb786b3425634f58fd2"

AGENT2_TRANSPORT_PR = 840
AGENT2_TRANSPORT_HEAD = "016e3c152d7d11e4fedcd79b819df337f8942983"
AGENT2_TRANSPORT_SOURCE_BLOB = "d6082ebd346ba8150321bf3d62f104c107a70300"
AGENT2_TRANSPORT_DEDICATED_RUN = 35511160478
AGENT2_TRANSPORT_TESTS_RUN = 35511160452
AGENT2_TRANSPORT_SCHEMA = "kokuno-a2-inner-leading-oscillatory-transport-v1"
AGENT2_TRANSPORT_MODULE = (
    "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport"
)
AGENT2_TRANSPORT_EVALUATOR = "evaluate_inner_leading_oscillatory_transport"

AGENT1_LAPLACIAN_PR = 839
AGENT1_LAPLACIAN_HEAD = "e96ee90144976a992b62d76d4361a94eb16bd91e"
AGENT1_LAPLACIAN_SOURCE_BLOB = "74b0e185e8e0bbb08695d493e0a091c305a03006"
AGENT1_LAPLACIAN_DEDICATED_RUN = 35510543422
AGENT1_LAPLACIAN_TESTS_RUN = 35510543420

# Frozen facts at this integration cut.  Queued/not-yet-run is not PASS.
AGENT1_LAPLACIAN_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT2_TRANSPORT_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_TRANSPORT_AUDIT_CONCLUSION: str | None = None

TRANSPORT_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center plus frozen oscillatory field",
    "formula": "dt(u_inner+u_osc)+(u_inner+u_osc).grad(u_inner+u_osc)-nu*Delta(u_inner+u_osc)",
    "viscosity": 0.01,
    "inner_advection_spatial_step": 1.0e-3,
    "oscillatory_advection_spatial_step": 1.0e-3,
    "inner_laplacian_spatial_step": 1.0e-3,
    "oscillatory_laplacian_spatial_step": 1.0e-3,
    "independent_reference": "public-summed-velocity-only-centered-FD4",
    "independent_fd4_steps": [4.0e-3, 2.0e-3, 1.0e-3],
    "fresh_probe_count": 12,
    "source_X_range": [0.172, 0.291],
    "time_range": [0.462, 0.542],
    "successive_stabilization_ratio_gate": 4.0,
    "finest_relative_rms_gate": 2.0e-2,
    "finest_relative_sampled_max_gate": 5.0e-2,
    "transport_rms_nontriviality_floor": 1.0e-10,
    "composite_velocity_rms_nontriviality_floor": 1.0e-8,
    "additive_velocity_closure_gate": 1.0e-14,
    "time_laplacian_viscous_transport_closure_gate": 1.0e-13,
    "post_observation_retuning_allowed": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_inner_leading_oscillatory_transport_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_cartesian_center_spatial_derivative_audit_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_spatial_derivative_audit_ingest_contract(parent)
    parent_status = parent["ingest_status"]

    admitted = bool(
        parent_status["inner_cartesian_center_spatial_derivative_ingest_admitted"]
        and AGENT1_LAPLACIAN_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT2_TRANSPORT_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_TRANSPORT_AUDIT_CONCLUSION == "pass"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "oscillatory_augmentation",
        "pipeline_substage": "strict_inner_pressure_forcing_free_transport_ingest",
        "parent_contract_sha256": parent["contract_sha256"],
        "agent2_transport_binding": {
            "pr": AGENT2_TRANSPORT_PR,
            "head": AGENT2_TRANSPORT_HEAD,
            "source_blob_sha": AGENT2_TRANSPORT_SOURCE_BLOB,
            "dedicated_run": AGENT2_TRANSPORT_DEDICATED_RUN,
            "tests_run": AGENT2_TRANSPORT_TESTS_RUN,
            "schema": AGENT2_TRANSPORT_SCHEMA,
            "module": AGENT2_TRANSPORT_MODULE,
            "evaluator": AGENT2_TRANSPORT_EVALUATOR,
            "protocol": dict(TRANSPORT_PROTOCOL),
        },
        "agent1_laplacian_binding": {
            "pr": AGENT1_LAPLACIAN_PR,
            "head": AGENT1_LAPLACIAN_HEAD,
            "source_blob_sha": AGENT1_LAPLACIAN_SOURCE_BLOB,
            "dedicated_run": AGENT1_LAPLACIAN_DEDICATED_RUN,
            "tests_run": AGENT1_LAPLACIAN_TESTS_RUN,
            "production_realization": "fixed centered Cartesian FD6",
            "production_spatial_step": 1.0e-3,
        },
        "evidence": {
            "parent_inner_spatial_derivative_ingest_admitted": bool(
                parent_status["inner_cartesian_center_spatial_derivative_ingest_admitted"]
            ),
            "agent1_laplacian_exact_head_ci_conclusion": (
                AGENT1_LAPLACIAN_EXACT_HEAD_CI_CONCLUSION
            ),
            "agent2_transport_exact_head_ci_conclusion": (
                AGENT2_TRANSPORT_EXACT_HEAD_CI_CONCLUSION
            ),
            "agent4_dedicated_transport_audit_present": False,
            "agent4_transport_audit_conclusion": AGENT4_TRANSPORT_AUDIT_CONCLUSION,
            "strict_inner_transport_scientifically_admitted": False,
        },
        "candidate_api_handoff": dict(parent["candidate_api_handoff"]),
        "differential_operator_handoff": dict(parent["differential_operator_handoff"]),
        "transport_operator_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "authority": "Agent2#840",
            "evaluator": (
                "Agent2#840.evaluate_inner_leading_oscillatory_transport"
            ),
            "velocity_scope": "u_inner+u_osc only",
            "includes_velocity_dt": True,
            "includes_full_strict_inner_advection": True,
            "includes_base_viscosity": True,
            "includes_pressure_gradient": False,
            "includes_restricted_forcing": False,
            "includes_correction_velocity": False,
            "includes_outer_global_leading_join": False,
            "complete_ns_momentum_residual": False,
            "usable_as_actual_defect_input_now": False,
            "usable_for_final_independent_pde_validation": False,
            "reason": (
                "the executable transport block is inner-only and pressure/forcing-free; "
                "its exact-head CI is unresolved and no dedicated Agent-4 transport audit exists"
            ),
        },
        "staged_residual_availability": {
            "leading_only_full_ns_residual_available": False,
            "leading_plus_oscillatory_full_ns_residual_available": False,
            "after_correction_full_ns_residual_available": False,
            "same_protocol_ST006_comparison_legal": False,
            "blocking_terms": [
                "final corrected/global leading velocity",
                "outer/global join",
                "matched pressure",
                "preregistered restricted forcing",
                "real full-candidate correction velocity",
                "Agent-4 full-candidate held-out validator receipt",
            ],
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_strict_inner_leading_oscillatory_transport_registered": True,
            "strict_inner_leading_oscillatory_transport_ingest_admitted": admitted,
            "strict_inner_transport_ready": False,
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
            "strict_inner_transport_surface_registered": True,
            "registration_not_scientific_admission": True,
            "agent1_laplacian_ci_promoted": False,
            "agent2_transport_ci_promoted": False,
            "agent4_transport_audit_invented": False,
            "strict_inner_transport_independently_admitted": False,
            "transport_relabelled_as_complete_ns_residual": False,
            "source_center_promoted_to_final_fixed_point": False,
            "outer_global_join_invented": False,
            "global_leading_velocity_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "correction_velocity_invented": False,
            "complete_candidate_api_promoted": False,
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


def validate_inner_leading_oscillatory_transport_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    parent = deterministic_cartesian_center_spatial_derivative_audit_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_spatial_derivative_audit_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    a2 = payload.get("agent2_transport_binding")
    if not isinstance(a2, Mapping):
        raise ValueError("missing Agent-2 transport binding")
    expected_a2 = {
        "pr": AGENT2_TRANSPORT_PR,
        "head": AGENT2_TRANSPORT_HEAD,
        "source_blob_sha": AGENT2_TRANSPORT_SOURCE_BLOB,
        "dedicated_run": AGENT2_TRANSPORT_DEDICATED_RUN,
        "tests_run": AGENT2_TRANSPORT_TESTS_RUN,
        "schema": AGENT2_TRANSPORT_SCHEMA,
        "module": AGENT2_TRANSPORT_MODULE,
        "evaluator": AGENT2_TRANSPORT_EVALUATOR,
    }
    for key, value in expected_a2.items():
        if a2.get(key) != value:
            raise ValueError(f"Agent-2 {key} drift")
    if a2.get("protocol") != TRANSPORT_PROTOCOL:
        raise ValueError("Agent-2 transport protocol drift")

    a1 = payload.get("agent1_laplacian_binding")
    if not isinstance(a1, Mapping):
        raise ValueError("missing Agent-1 Laplacian binding")
    expected_a1 = {
        "pr": AGENT1_LAPLACIAN_PR,
        "head": AGENT1_LAPLACIAN_HEAD,
        "source_blob_sha": AGENT1_LAPLACIAN_SOURCE_BLOB,
        "dedicated_run": AGENT1_LAPLACIAN_DEDICATED_RUN,
        "tests_run": AGENT1_LAPLACIAN_TESTS_RUN,
        "production_realization": "fixed centered Cartesian FD6",
        "production_spatial_step": 1.0e-3,
    }
    if dict(a1) != expected_a1:
        raise ValueError("Agent-1 Laplacian binding drift")

    expected_evidence = {
        "parent_inner_spatial_derivative_ingest_admitted": False,
        "agent1_laplacian_exact_head_ci_conclusion": None,
        "agent2_transport_exact_head_ci_conclusion": None,
        "agent4_dedicated_transport_audit_present": False,
        "agent4_transport_audit_conclusion": None,
        "strict_inner_transport_scientifically_admitted": False,
    }
    if payload.get("evidence") != expected_evidence:
        raise ValueError("evidence laundering")

    if payload.get("candidate_api_handoff") != parent["candidate_api_handoff"]:
        raise ValueError("candidate API drift")
    if payload.get("differential_operator_handoff") != parent[
        "differential_operator_handoff"
    ]:
        raise ValueError("differential operator handoff drift")

    transport = payload.get("transport_operator_handoff")
    if not isinstance(transport, Mapping):
        raise ValueError("missing transport handoff")
    required_transport = {
        "registered": True,
        "status": "registered_unresolved",
        "authority": "Agent2#840",
        "evaluator": "Agent2#840.evaluate_inner_leading_oscillatory_transport",
        "velocity_scope": "u_inner+u_osc only",
        "includes_velocity_dt": True,
        "includes_full_strict_inner_advection": True,
        "includes_base_viscosity": True,
        "includes_pressure_gradient": False,
        "includes_restricted_forcing": False,
        "includes_correction_velocity": False,
        "includes_outer_global_leading_join": False,
        "complete_ns_momentum_residual": False,
        "usable_as_actual_defect_input_now": False,
        "usable_for_final_independent_pde_validation": False,
    }
    for key, value in required_transport.items():
        if transport.get(key) != value:
            raise ValueError(f"transport handoff {key} drift")

    staged = payload.get("staged_residual_availability")
    if not isinstance(staged, Mapping):
        raise ValueError("missing staged residual availability")
    for key in (
        "leading_only_full_ns_residual_available",
        "leading_plus_oscillatory_full_ns_residual_available",
        "after_correction_full_ns_residual_available",
        "same_protocol_ST006_comparison_legal",
    ):
        if staged.get(key) is not False:
            raise ValueError(f"premature staged residual promotion: {key}")
    if not staged.get("blocking_terms"):
        raise ValueError("missing staged residual blockers")

    status = payload.get("ingest_status")
    if not isinstance(status, Mapping):
        raise ValueError("missing ingest status")
    expected_status = dict(parent["ingest_status"])
    expected_status.update(
        {
            "typed_strict_inner_leading_oscillatory_transport_registered": True,
            "strict_inner_leading_oscillatory_transport_ingest_admitted": False,
            "strict_inner_transport_ready": False,
            "outer_global_leading_join_materialized": False,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
            "leading_ready": False,
        }
    )
    if dict(status) != expected_status:
        raise ValueError("ingest status promotion/drift")

    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state drift")
    if payload.get("baseline_vs_kokuno") != parent["baseline_vs_kokuno"]:
        raise ValueError("baseline comparison drift")
    if payload.get("final_project_gates_unchanged") != parent[
        "final_project_gates_unchanged"
    ]:
        raise ValueError("project gate drift")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    true_keys = {
        "strict_inner_transport_surface_registered",
        "registration_not_scientific_admission",
    }
    for key in true_keys:
        if truth.get(key) is not True:
            raise ValueError(f"required truth marker lost: {key}")
    for key, value in truth.items():
        if key in true_keys:
            continue
        if value is not False:
            raise ValueError(f"truth promotion: {key}")

    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = deterministic_inner_leading_oscillatory_transport_ingest_contract(
        exact_head=args.exact_head
    )
    validate_inner_leading_oscillatory_transport_ingest_contract(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"schema": SCHEMA, "contract_sha256": payload["contract_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
