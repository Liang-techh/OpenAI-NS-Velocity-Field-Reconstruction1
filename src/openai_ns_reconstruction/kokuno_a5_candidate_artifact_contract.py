"""Fail-closed Agent-5 contract for the future Kokuno candidate artifact.

This module does not materialize a velocity candidate.  It freezes the typed
serialization/API contract that a future real Kokuno candidate must satisfy
before it may be called save/load ready or sent to Agent 4 for the formal
held-out PDE gate.

The contract is intentionally independent of any surrogate regression field.
It imports the current Agent-5 admission state, preserves the registered final
scientific gates, and keeps candidate/export/PDE readiness false until a later
pipeline supplies real velocity/pressure/restricted-forcing bindings and an
independent Agent-4 validation receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_source_leading_admission_firewall import FINAL_PROJECT_GATES
from .kokuno_a5_source_leading_admission_firewall_v4 import ST006_BASELINE
from .kokuno_a5_source_leading_admission_firewall_v5 import derive_state

SCHEMA = "kokuno-a5-candidate-artifact-contract-v1"
TASK = "KOKUNO-A5-CANDIDATE-ARTIFACT-CONTRACT-066"
PARENT_A5_PR = 756
PARENT_A5_HEAD = "44b68dc387bb4482f4a836a06b14698e93a70a71"

PINNED_KOKUNO_SOURCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "source_file": "navier-stokes/navier_stokes_workbench.tex",
    "reader_date": "2026-09-09",
    "paper_exact_claim": False,
}

PIPELINE_ORDER = [
    "source_profile_ingest",
    "leading_candidate",
    "oscillatory_augmentation",
    "mean_radial_corrections",
    "finite_correction_cycle",
    "candidate_artifact",
    "independent_validation",
    "report",
]

PUBLIC_API_CONTRACT = {
    "velocity": {
        "signature": "velocity(points, t) -> (..., 3)",
        "required_for_candidate_artifact": True,
    },
    "velocity_dt": {
        "signature": "velocity_dt(points, t) -> (..., 3)",
        "required_for_candidate_artifact": True,
    },
    "pressure": {
        "signature": "pressure(points, t) -> (...) or (..., 1)",
        "required_for_candidate_artifact": True,
    },
    "forcing": {
        "signature": "forcing(points, t) -> (..., 3)",
        "required_for_candidate_artifact": True,
        "must_be_preregistered_restricted_family": True,
        "may_not_be_defined_from_residual": True,
    },
}

PARAMETER_SERIALIZATION = {
    "format": "canonical-json-utf8",
    "sort_keys": True,
    "separators": [",", ":"],
    "allow_nan": False,
    "parameter_sha256_required": True,
    "source_and_version_provenance_required": True,
}


def _canonical(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def _digest_payload(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    return hashlib.sha256(_canonical(body)).hexdigest()


def _stage_state() -> dict[str, bool]:
    current = derive_state()
    return {
        "leading_ready": bool(current["leading_ready"]),
        "oscillatory_ready": bool(current["oscillatory_ready"]),
        "correction_ready": bool(current["correction_ready"]),
        "velocity_export_ready": bool(current["velocity_export_ready"]),
        "pde_validated": bool(current["pde_validated"]),
    }


def deterministic_contract(*, exact_head: str | None = None) -> dict[str, Any]:
    """Return the current fail-closed candidate-artifact contract.

    This is a contract artifact only.  It deliberately contains no candidate
    coefficients and no evaluator import binding, because the global leading
    field and real correction cycle are not yet available on the A5 stack.
    """

    state = _stage_state()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "source_provenance": dict(PINNED_KOKUNO_SOURCE),
        "pipeline_order": list(PIPELINE_ORDER),
        "public_api_contract": json.loads(json.dumps(PUBLIC_API_CONTRACT)),
        "parameter_serialization": dict(PARAMETER_SERIALIZATION),
        "stage_state": state,
        "candidate_payload": {
            "materialized": False,
            "candidate_id": None,
            "parameter_payload": None,
            "parameter_sha256": None,
            "evaluator_bindings": None,
            "validation_receipt": None,
        },
        "promotion_rules": {
            "leading_ready_requires": [
                "global leading velocity",
                "matched pressure",
                "independently validated restricted-forcing semantics",
            ],
            "correction_ready_requires": [
                "real full-candidate defect consumed",
                "source-certified correction inputs",
                "real correction velocity",
                "guarded finite correction cycle",
            ],
            "velocity_export_ready_requires": [
                "materialized candidate artifact",
                "all four public evaluator bindings",
                "canonical parameter digest",
                "save/load round-trip smoke",
            ],
            "pde_validated_requires": [
                "Agent-4 independent held-out validation",
                "normalized momentum max <= 1e-3",
                "normalized momentum L2 <= 1e-3",
                "divergence max <= 1e-5",
                "divergence L2 <= 1e-5",
            ],
        },
        "baseline_vs_kokuno": {
            "retained_repository_baseline": dict(ST006_BASELINE),
            "kokuno_same_protocol_full_candidate_residual_available": False,
            "comparison_performed": False,
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "contract_only_not_candidate": True,
            "surrogate_promoted_to_candidate": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "kokuno_replay_used_as_final_independent_validation": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "leading_ready": state["leading_ready"],
            "correction_ready": state["correction_ready"],
            "velocity_export_ready": state["velocity_export_ready"],
            "pde_validated": state["pde_validated"],
        },
    }
    payload["contract_sha256"] = _digest_payload(payload)
    return payload


def validate_contract(payload: Mapping[str, Any]) -> None:
    """Reject scientific-state laundering or weakening of the artifact contract."""

    if payload.get("schema") != SCHEMA:
        raise ValueError("unexpected Kokuno candidate-artifact schema")
    if payload.get("pipeline_order") != PIPELINE_ORDER:
        raise RuntimeError("pipeline order changed")
    if payload.get("public_api_contract") != PUBLIC_API_CONTRACT:
        raise RuntimeError("public velocity/pressure/forcing API contract changed")
    if payload.get("parameter_serialization") != PARAMETER_SERIALIZATION:
        raise RuntimeError("parameter serialization contract changed")
    if payload.get("source_provenance") != PINNED_KOKUNO_SOURCE:
        raise RuntimeError("Kokuno source/version provenance changed")
    if payload.get("stage_state") != _stage_state():
        raise RuntimeError("artifact contract does not match current A5 stage state")

    candidate = payload["candidate_payload"]
    if candidate != {
        "materialized": False,
        "candidate_id": None,
        "parameter_payload": None,
        "parameter_sha256": None,
        "evaluator_bindings": None,
        "validation_receipt": None,
    }:
        raise RuntimeError("nonexistent Kokuno candidate was materialized in the contract")

    gates = payload["final_project_gates_unchanged"]
    if gates != FINAL_PROJECT_GATES:
        raise RuntimeError("final project gates changed")
    baseline = payload["baseline_vs_kokuno"]
    if baseline["retained_repository_baseline"] != ST006_BASELINE:
        raise RuntimeError("ST006 baseline changed")
    if baseline["kokuno_same_protocol_full_candidate_residual_available"] is not False:
        raise RuntimeError("nonexistent full Kokuno residual was promoted")
    if baseline["comparison_performed"] is not False:
        raise RuntimeError("baseline comparison cannot run before a full Kokuno candidate")

    truth = payload["truth_boundary"]
    required_false = (
        "surrogate_promoted_to_candidate",
        "free_residual_defined_forcing_allowed",
        "threshold_relaxed",
        "kokuno_replay_used_as_final_independent_validation",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    )
    bad = [key for key in required_false if truth[key] is not False]
    if bad:
        raise RuntimeError(f"truth-boundary state was promoted: {bad}")
    if truth["contract_only_not_candidate"] is not True:
        raise RuntimeError("contract artifact was relabeled as a candidate")

    forcing = payload["public_api_contract"]["forcing"]
    if forcing["must_be_preregistered_restricted_family"] is not True:
        raise RuntimeError("forcing restriction was weakened")
    if forcing["may_not_be_defined_from_residual"] is not True:
        raise RuntimeError("residual-defined forcing escape hatch introduced")

    expected_digest = _digest_payload(payload)
    if payload.get("contract_sha256") != expected_digest:
        raise RuntimeError("candidate-artifact contract digest mismatch")


def write_contract(path: Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_contract(exact_head=exact_head)
    validate_contract(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def load_contract(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    validate_contract(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--exact-head")
    args = parser.parse_args()
    payload = deterministic_contract(exact_head=args.exact_head)
    validate_contract(payload)
    if args.output is None:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        write_contract(args.output, exact_head=args.exact_head)


if __name__ == "__main__":
    main()
