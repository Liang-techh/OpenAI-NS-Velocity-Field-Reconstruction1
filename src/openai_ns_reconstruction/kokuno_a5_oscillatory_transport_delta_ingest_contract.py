"""Fail-closed A5 registration of the strict-inner oscillatory transport delta.

Agent 2 #849 exposes the before/after pressure/forcing-free transport increment

    delta_T_osc = T(u_inner + u_osc) - T(u_inner),
    T(u) = d_t u + (u . grad)u - 0.01 Delta u,

for the current PA.10 *inner contraction-center* background. Agent 4 #851
independently reconstructs the two transports from public velocity values only
with a centered FD8 realization and audits their difference.

This A5 seam records typed provenance and the preregistered staged diagnostic
without copying either sibling implementation. It does not promote the scoped
transport delta into a Navier--Stokes residual: global joining, matched pressure,
restricted forcing, and a real Agent-3 correction velocity are still absent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_inner_leading_oscillatory_transport_ingest_contract import (
    deterministic_inner_leading_oscillatory_transport_ingest_contract,
    validate_inner_leading_oscillatory_transport_ingest_contract,
)

SCHEMA = "kokuno-a5-oscillatory-transport-delta-ingest-contract-v1"
TASK = "KOKUNO-A5-OSCILLATORY-TRANSPORT-DELTA-INGEST-076"
PARENT_A5_PR = 843
PARENT_A5_HEAD = "7ba3fc1dd2aab2f09b0c4708a6b45bf08589f932"

AGENT2_DELTA_PR = 849
AGENT2_DELTA_HEAD = "616e61abd2420a7c7aff3590a08393360039cd1a"
AGENT2_DELTA_SOURCE_BLOB = "0f4d1614945688a6f5d0d6f4192979cabbbe9da8"
AGENT2_DELTA_DEDICATED_RUN = 35513939185
AGENT2_DELTA_TESTS_RUN = 35513939315
AGENT2_DELTA_SCHEMA = "kokuno-a2-inner-leading-oscillatory-transport-delta-v1"
AGENT2_DELTA_MODULE = (
    "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport_delta"
)
AGENT2_DELTA_EVALUATOR = "evaluate_inner_leading_oscillatory_transport_delta"

AGENT4_DELTA_AUDIT_PR = 851
AGENT4_DELTA_AUDIT_HEAD = "579bc7555ae3b2e74adb10669b95d75aeef3dd53"
AGENT4_DELTA_AUDIT_SOURCE_BLOB = "a1ce30cc1067a854f1d4c252dd24b5b1e9fb9a64"
AGENT4_DELTA_AUDIT_DEDICATED_RUN = 35514693095
AGENT4_DELTA_AUDIT_TESTS_RUN = 35514693112
AGENT4_DELTA_AUDIT_SCHEMA = (
    "kokuno-a4-oscillatory-transport-delta-independent-audit-v1"
)
AGENT4_DELTA_AUDIT_MODULE = (
    "openai_ns_reconstruction.kokuno_a4_oscillatory_transport_delta_independent_audit"
)

# Exact-head Actions are unresolved at this integration cut. Queued is not PASS.
AGENT2_DELTA_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_DELTA_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_DELTA_AUDIT_CONCLUSION: str | None = None

AGENT2_DELTA_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center before/after frozen oscillatory field",
    "formula": "delta_T_osc=T(u_inner+u_osc)-T(u_inner)",
    "transport_formula": "dt(u)+u.grad(u)-nu*Delta(u)",
    "viscosity": 0.01,
    "fresh_probe_count": 12,
    "source_X_range": [0.166, 0.286],
    "time_range": [0.464, 0.536],
    "independent_reference": "public-velocity-only-centered-FD4-before-after",
    "fd4_steps": [4.0e-3, 2.0e-3, 1.0e-3],
    "delta_refinement_ratio_gate": 3.0,
    "delta_finest_relative_rms_gate": 5.0e-2,
    "delta_finest_relative_sampled_max_gate": 1.0e-1,
    "inner_finest_relative_rms_gate": 2.0e-2,
    "inner_finest_relative_sampled_max_gate": 5.0e-2,
    "inner_transport_rms_nontriviality_floor": 1.0e-10,
    "delta_transport_rms_nontriviality_floor": 1.0e-10,
    "before_after_decomposition_closure_gate": 5.0e-13,
    "velocity_additive_closure_gate": 1.0e-14,
    "ratio_R_is_observation_only": True,
    "post_observation_retuning_allowed": False,
}

AGENT4_DELTA_PROTOCOL = {
    "independent_reference": "inner-and-total-public-velocity-only-centered-FD8",
    "seed": 9173511,
    "random_sample_count": 96,
    "axis_and_axis_near_probe_count": 6,
    "source_X_range": [0.155, 0.265],
    "source_eta_range": [-0.26, 0.26],
    "time_range": [0.445, 0.555],
    "fd8_steps": [2.0e-3, 1.0e-3, 5.0e-4],
    "inner_transport_consistency_relative_rms_gate": 1.5e-2,
    "inner_transport_consistency_relative_sampled_max_gate": 4.0e-2,
    "delta_transport_consistency_relative_rms_gate": 5.0e-2,
    "delta_transport_consistency_relative_sampled_max_gate": 1.0e-1,
    "inner_total_refinement_ratio_gate": 8.0,
    "delta_refinement_ratio_gate": 6.0,
    "refinement_floor": 2.0e-8,
    "divergence_sampled_max_gate": 1.0e-5,
    "divergence_sampled_rms_gate": 1.0e-5,
    "negative_controls": [
        "production-delta_T_osc-times-0.90",
        "production-delta_T_osc-sign-flip",
        "+1e-3-independent-total-divergence",
    ],
    "ratio_R_osc_is_observation_only": True,
    "viscosity_sensitivity_only": "nu +/-0.1%",
    "post_observation_retuning_allowed": False,
}

NORM_SCOPE_FIREWALL = {
    "transport_delta_is_complete_ns_momentum_residual": False,
    "transport_ratio_is_cr001_acceptance_metric": False,
    "transport_ratio_directly_comparable_to_ST006_full_residual": False,
    "agent4_divergence_rms_is_sampled_rms": True,
    "agent4_divergence_rms_is_canonical_volume_weighted_L2": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_oscillatory_transport_delta_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_inner_leading_oscillatory_transport_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_inner_leading_oscillatory_transport_ingest_contract(parent)
    parent_status = parent["ingest_status"]

    admitted = bool(
        parent_status["strict_inner_leading_oscillatory_transport_ingest_admitted"]
        and AGENT2_DELTA_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_DELTA_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_DELTA_AUDIT_CONCLUSION == "pass"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "oscillatory_augmentation",
        "pipeline_substage": "before_after_strict_inner_transport_delta",
        "parent_contract_sha256": parent["contract_sha256"],
        "agent2_delta_binding": {
            "pr": AGENT2_DELTA_PR,
            "head": AGENT2_DELTA_HEAD,
            "source_blob_sha": AGENT2_DELTA_SOURCE_BLOB,
            "dedicated_run": AGENT2_DELTA_DEDICATED_RUN,
            "tests_run": AGENT2_DELTA_TESTS_RUN,
            "schema": AGENT2_DELTA_SCHEMA,
            "module": AGENT2_DELTA_MODULE,
            "evaluator": AGENT2_DELTA_EVALUATOR,
            "protocol": dict(AGENT2_DELTA_PROTOCOL),
        },
        "agent4_delta_audit_binding": {
            "pr": AGENT4_DELTA_AUDIT_PR,
            "head": AGENT4_DELTA_AUDIT_HEAD,
            "source_blob_sha": AGENT4_DELTA_AUDIT_SOURCE_BLOB,
            "dedicated_run": AGENT4_DELTA_AUDIT_DEDICATED_RUN,
            "tests_run": AGENT4_DELTA_AUDIT_TESTS_RUN,
            "schema": AGENT4_DELTA_AUDIT_SCHEMA,
            "module": AGENT4_DELTA_AUDIT_MODULE,
            "audited_agent2_pr": AGENT2_DELTA_PR,
            "audited_agent2_head": AGENT2_DELTA_HEAD,
            "protocol": dict(AGENT4_DELTA_PROTOCOL),
            "norm_scope_firewall": dict(NORM_SCOPE_FIREWALL),
        },
        "evidence": {
            "parent_strict_inner_transport_ingest_admitted": bool(
                parent_status["strict_inner_leading_oscillatory_transport_ingest_admitted"]
            ),
            "agent2_delta_exact_head_ci_conclusion": AGENT2_DELTA_EXACT_HEAD_CI_CONCLUSION,
            "agent4_dedicated_delta_audit_present": True,
            "agent4_delta_exact_head_ci_conclusion": AGENT4_DELTA_EXACT_HEAD_CI_CONCLUSION,
            "agent4_delta_audit_conclusion": AGENT4_DELTA_AUDIT_CONCLUSION,
            "agent4_delta_scientific_receipt_admitted": False,
            "oscillatory_transport_delta_scientifically_admitted": False,
            "staged_transport_ratio_numeric_value": None,
        },
        "candidate_api_handoff": dict(parent["candidate_api_handoff"]),
        "differential_operator_handoff": dict(parent["differential_operator_handoff"]),
        "transport_operator_handoff": dict(parent["transport_operator_handoff"]),
        "oscillatory_transport_delta_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "construction_authority": "Agent2#849",
            "independent_audit_available": True,
            "independent_audit_status": "registered_unresolved",
            "independent_audit_authority": "Agent4#851",
            "evaluator": "Agent2#849.evaluate_inner_leading_oscillatory_transport_delta",
            "quantity": "T(u_inner+u_osc)-T(u_inner)",
            "pressure_included": False,
            "restricted_forcing_included": False,
            "correction_velocity_included": False,
            "outer_global_join_included": False,
            "complete_ns_momentum_residual": False,
            "ratio_R_definition": "RMS(T_inner+osc)/RMS(T_inner)",
            "ratio_R_observation_only": True,
            "ratio_R_available_now": False,
            "usable_for_staged_transport_attribution_if_passes": True,
            "usable_as_full_candidate_defect_now": False,
            "usable_for_final_independent_pde_validation": False,
        },
        "staged_comparison": {
            "leading_only_transport_surface_available": True,
            "leading_plus_oscillatory_transport_surface_available": True,
            "oscillatory_transport_delta_surface_available": True,
            "independent_before_after_delta_audit_registered": True,
            "numeric_before_after_ratio_available": False,
            "leading_only_full_ns_residual_available": False,
            "leading_plus_oscillatory_full_ns_residual_available": False,
            "after_correction_full_ns_residual_available": False,
            "same_protocol_ST006_comparison_legal": False,
            "reason": (
                "the staged quantity omits matched pressure, restricted forcing, "
                "global joining, and correction transport"
            ),
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_oscillatory_transport_delta_registered": True,
            "typed_oscillatory_transport_delta_independent_audit_registered": True,
            "oscillatory_transport_delta_ingest_admitted": admitted,
            "staged_transport_attribution_ready": False,
            "leading_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "oscillatory_transport_delta_surface_registered": True,
            "agent4_delta_audit_surface_registered": True,
            "registration_not_scientific_admission": True,
            "transport_delta_is_not_ns_residual": True,
            "transport_ratio_is_observation_only": True,
            "agent2_delta_ci_promoted": False,
            "agent4_delta_ci_promoted": False,
            "agent4_delta_receipt_invented": False,
            "oscillatory_delta_independently_admitted": False,
            "transport_ratio_laundered_as_cr001_acceptance": False,
            "sampled_divergence_rms_laundered_as_cr001_volume_l2": False,
            "full_ns_residual_invented": False,
            "premature_ST006_comparison_performed": False,
            "matched_pressure_invented": False,
            "restricted_forcing_invented": False,
            "correction_velocity_invented": False,
            "outer_global_join_invented": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "leading_ready": False,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
    }
    payload["contract_sha256"] = _digest(payload)
    return payload


def validate_oscillatory_transport_delta_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    parent = deterministic_inner_leading_oscillatory_transport_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_inner_leading_oscillatory_transport_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    a2 = payload.get("agent2_delta_binding")
    if not isinstance(a2, Mapping):
        raise ValueError("missing Agent-2 delta binding")
    for key, value in {
        "pr": AGENT2_DELTA_PR,
        "head": AGENT2_DELTA_HEAD,
        "source_blob_sha": AGENT2_DELTA_SOURCE_BLOB,
        "dedicated_run": AGENT2_DELTA_DEDICATED_RUN,
        "tests_run": AGENT2_DELTA_TESTS_RUN,
        "schema": AGENT2_DELTA_SCHEMA,
        "module": AGENT2_DELTA_MODULE,
        "evaluator": AGENT2_DELTA_EVALUATOR,
    }.items():
        if a2.get(key) != value:
            raise ValueError(f"Agent-2 delta {key} drift")
    if a2.get("protocol") != AGENT2_DELTA_PROTOCOL:
        raise ValueError("Agent-2 delta protocol drift")

    a4 = payload.get("agent4_delta_audit_binding")
    if not isinstance(a4, Mapping):
        raise ValueError("missing Agent-4 delta audit binding")
    for key, value in {
        "pr": AGENT4_DELTA_AUDIT_PR,
        "head": AGENT4_DELTA_AUDIT_HEAD,
        "source_blob_sha": AGENT4_DELTA_AUDIT_SOURCE_BLOB,
        "dedicated_run": AGENT4_DELTA_AUDIT_DEDICATED_RUN,
        "tests_run": AGENT4_DELTA_AUDIT_TESTS_RUN,
        "schema": AGENT4_DELTA_AUDIT_SCHEMA,
        "module": AGENT4_DELTA_AUDIT_MODULE,
        "audited_agent2_pr": AGENT2_DELTA_PR,
        "audited_agent2_head": AGENT2_DELTA_HEAD,
    }.items():
        if a4.get(key) != value:
            raise ValueError(f"Agent-4 delta {key} drift")
    if a4.get("protocol") != AGENT4_DELTA_PROTOCOL:
        raise ValueError("Agent-4 delta protocol drift")
    if a4.get("norm_scope_firewall") != NORM_SCOPE_FIREWALL:
        raise ValueError("Agent-4 delta norm scope drift")

    expected_evidence = {
        "parent_strict_inner_transport_ingest_admitted": False,
        "agent2_delta_exact_head_ci_conclusion": None,
        "agent4_dedicated_delta_audit_present": True,
        "agent4_delta_exact_head_ci_conclusion": None,
        "agent4_delta_audit_conclusion": None,
        "agent4_delta_scientific_receipt_admitted": False,
        "oscillatory_transport_delta_scientifically_admitted": False,
        "staged_transport_ratio_numeric_value": None,
    }
    if payload.get("evidence") != expected_evidence:
        raise ValueError("evidence laundering")

    for key in (
        "candidate_api_handoff",
        "differential_operator_handoff",
        "transport_operator_handoff",
        "stage_state",
        "baseline_vs_kokuno",
        "final_project_gates_unchanged",
    ):
        if payload.get(key) != parent[key]:
            raise ValueError(f"parent {key} drift")

    handoff = payload.get("oscillatory_transport_delta_handoff")
    if not isinstance(handoff, Mapping):
        raise ValueError("missing oscillatory delta handoff")
    required_handoff = {
        "registered": True,
        "status": "registered_unresolved",
        "construction_authority": "Agent2#849",
        "independent_audit_available": True,
        "independent_audit_status": "registered_unresolved",
        "independent_audit_authority": "Agent4#851",
        "evaluator": "Agent2#849.evaluate_inner_leading_oscillatory_transport_delta",
        "complete_ns_momentum_residual": False,
        "ratio_R_observation_only": True,
        "ratio_R_available_now": False,
        "usable_for_staged_transport_attribution_if_passes": True,
        "usable_as_full_candidate_defect_now": False,
        "usable_for_final_independent_pde_validation": False,
    }
    for key, value in required_handoff.items():
        if handoff.get(key) != value:
            raise ValueError(f"oscillatory delta handoff {key} drift")
    for key in (
        "pressure_included",
        "restricted_forcing_included",
        "correction_velocity_included",
        "outer_global_join_included",
    ):
        if handoff.get(key) is not False:
            raise ValueError(f"oscillatory delta scope promotion: {key}")

    staged = payload.get("staged_comparison")
    if not isinstance(staged, Mapping):
        raise ValueError("missing staged comparison")
    for key in (
        "leading_only_transport_surface_available",
        "leading_plus_oscillatory_transport_surface_available",
        "oscillatory_transport_delta_surface_available",
        "independent_before_after_delta_audit_registered",
    ):
        if staged.get(key) is not True:
            raise ValueError(f"lost staged transport surface: {key}")
    for key in (
        "numeric_before_after_ratio_available",
        "leading_only_full_ns_residual_available",
        "leading_plus_oscillatory_full_ns_residual_available",
        "after_correction_full_ns_residual_available",
        "same_protocol_ST006_comparison_legal",
    ):
        if staged.get(key) is not False:
            raise ValueError(f"premature staged promotion: {key}")

    expected_status = dict(parent["ingest_status"])
    expected_status.update(
        {
            "typed_oscillatory_transport_delta_registered": True,
            "typed_oscillatory_transport_delta_independent_audit_registered": True,
            "oscillatory_transport_delta_ingest_admitted": False,
            "staged_transport_attribution_ready": False,
            "leading_ready": False,
        }
    )
    if payload.get("ingest_status") != expected_status:
        raise ValueError("ingest status promotion/drift")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    true_keys = {
        "oscillatory_transport_delta_surface_registered",
        "agent4_delta_audit_surface_registered",
        "registration_not_scientific_admission",
        "transport_delta_is_not_ns_residual",
        "transport_ratio_is_observation_only",
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

    payload = deterministic_oscillatory_transport_delta_ingest_contract(
        exact_head=args.exact_head
    )
    validate_oscillatory_transport_delta_ingest_contract(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"schema": SCHEMA, "contract_sha256": payload["contract_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
