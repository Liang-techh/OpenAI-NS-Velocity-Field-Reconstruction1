"""CR002 governance for binding an ST052-M replay artifact to its validation receipt.

The exact PR #619 CI artifact is internally consistent, but its manifest builder
currently accepts candidate and validation paths independently and does not check
that ``validation.json[candidate_sha256]`` equals the candidate file checksum.
This module records that narrow seam and provides a reusable fail-closed pair
check.  It does not materialize or promote any velocity candidate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "configs"
    / "st052m_parent_capsule_linkage_governance_contract.json"
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_contract(path: str | Path | None = None) -> dict[str, Any]:
    target = CONTRACT_PATH if path is None else Path(path)
    return json.loads(target.read_text(encoding="utf-8"))


def verify_candidate_validation_binding(
    candidate_path: str | Path, validation_path: str | Path
) -> dict[str, Any]:
    """Fail closed unless the validation report names the exact candidate bytes."""
    candidate = Path(candidate_path)
    validation = Path(validation_path)
    candidate_sha = _sha256_file(candidate)
    report = json.loads(validation.read_text(encoding="utf-8"))
    report_sha = report.get("candidate_sha256")
    if report_sha != candidate_sha:
        raise ValueError(
            "candidate/validation SHA mismatch: "
            f"candidate={candidate_sha}, validation.candidate_sha256={report_sha!r}"
        )
    if report.get("pde_validated") is not False:
        raise ValueError("validation.pde_validated must remain false for this governed receipt")
    if report.get("all_numeric_gates_pass") is not False:
        raise ValueError("validation.all_numeric_gates_pass must remain false")
    gates = report.get("gates")
    if not isinstance(gates, dict):
        raise ValueError("validation.gates missing")
    if gates.get("momentum_max") is not False or gates.get("momentum_L2") is not False:
        raise ValueError("expected ST052-M momentum rejection is absent")
    return {
        "candidate_sha256": candidate_sha,
        "validation_sha256": _sha256_file(validation),
        "validation_candidate_sha256": report_sha,
        "candidate_validation_sha_pair_matches": True,
        "pde_validated": False,
    }


def _require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def audit_contract(contract: dict[str, Any] | None = None) -> dict[str, Any]:
    """Validate the immutable governance assertions for the exact upstream receipt."""
    data = load_contract() if contract is None else contract
    _require(data.get("schema_version") == 1, "schema_version drift")
    _require(
        data.get("contract_id") == "st052m_parent_capsule_linkage_governance_v1",
        "contract_id drift",
    )
    base = data["integration_base"]
    _require(base["branch"] == "codex/cr001-constraints", "integration branch drift")
    _require(
        base["head"] == "d3286b585d9e6dba6f3846ab94e9f7388d17176c",
        "integration head drift",
    )

    upstream = data["upstream_parent_capsule"]
    _require(upstream["pull_request"] == 619, "upstream PR drift")
    _require(
        upstream["head"] == "088d7263c403b3aa3a3818ea527ea0c5845e3089",
        "upstream head drift",
    )
    _require(upstream["base_branch"] == "main", "upstream base classification drift")
    _require(upstream["workflow_conclusion"] == "success", "upstream CI receipt drift")

    pair = data["exact_artifact_pair"]
    _require(
        pair["candidate_file_sha256"] == pair["validation_report_candidate_sha256"],
        "recorded candidate/validation SHA pair diverged",
    )
    _require(pair["candidate_validation_sha_pair_matches"] is True, "pair match must stay true")
    _require(pair["validation_experiment"] == "additional_diagnostic", "validation scope drift")
    _require(pair["validation_seed"] == 9175291, "diagnostic seed drift")
    for key in (
        "pde_validated",
        "all_numeric_gates_pass",
        "momentum_max_gate",
        "momentum_l2_gate",
    ):
        _require(pair[key] is False, f"exact_artifact_pair.{key} must stay false")
    _require(
        pair["replay_identity_scope"]
        == "source_commit_recipe_parent_modifier_reference_lineage_not_global_function_hash",
        "replay identity scope laundering",
    )

    finding = data["governance_finding"]
    _require(
        finding["upstream_exact_artifact_pair_is_internally_consistent"] is True,
        "exact artifact consistency drift",
    )
    _require(
        finding["upstream_builder_checks_validation_candidate_sha_against_candidate_file"]
        is False,
        "builder gap classification drift",
    )
    _require(
        finding["cross_file_substitution_gap_present_at_builder_api"] is True,
        "builder substitution gap must remain explicit",
    )
    _require(
        finding[
            "candidate_validation_sha_binding_required_before_manifest_claims_saved_file_validation"
        ]
        is True,
        "candidate-validation binding prerequisite removed",
    )
    _require(
        finding["co_generated_files_in_ci_are_execution_evidence_not_a_general_builder_invariant"]
        is True,
        "CI co-generation overreach",
    )
    _require(
        finding["replay_identity_may_not_be_laundered_into_candidate_file_identity"] is True,
        "replay identity laundering",
    )

    delivery = data["delivery_state"]
    _require(
        delivery["live_integration_contains_parent_behavioral_binding"] is True,
        "live behavioral binding state drift",
    )
    _require(
        delivery["live_integration_contains_parent_replay_capsule"] is False,
        "open PR may not be counted as integrated",
    )
    _require(
        delivery["parent_replay_artifact_evidence_exists_on_open_pr"] is True,
        "upstream artifact evidence lost",
    )
    for key in (
        "complete_temporal_child_materialized",
        "complete_candidate_save_load_ready",
        "experimental_st052_temporal_velocity_export_ready",
    ):
        _require(delivery[key] is False, f"delivery_state.{key} premature promotion")
    _require(
        delivery["canonical_eq45_velocity_export_ready"] is True,
        "canonical Eq45 readiness drift",
    )

    for key, value in data["truth_states"].items():
        _require(value is False, f"truth_states.{key} premature promotion")

    classes = data["source_classification"]
    _require(
        set(classes)
        == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "source classification categories drift",
    )

    cr = data["cr001_snapshot"]
    _require(cr["nu"] == 0.01, "CR001 nu drift")
    _require(cr["physical_domain"] == "R^3", "CR001 domain drift")
    _require(cr["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "CR001 box drift")
    _require(cr["time_interval"] == [0.25, 0.75], "CR001 time drift")
    _require(cr["formal_validation_seed"] == 914027, "CR001 formal seed drift")
    _require(
        pair["validation_seed"] != cr["formal_validation_seed"],
        "diagnostic seed laundered into formal acceptance",
    )
    _require(cr["held_out_points"] == 4096, "CR001 held-out count drift")
    _require(
        cr["derivative_steps"] == [0.02, 0.01, 0.005],
        "CR001 derivative ladder drift",
    )
    _require(
        cr["quadrature_orders_per_axis"] == [24, 48, 96],
        "CR001 quadrature ladder drift",
    )
    _require(
        cr["divergence_max"] == 1e-5 and cr["divergence_l2"] == 1e-5,
        "CR001 divergence threshold drift",
    )
    _require(
        cr["momentum_max"] == 1e-3 and cr["momentum_l2"] == 1e-3,
        "CR001 momentum threshold drift",
    )
    _require(
        cr["free_residual_defined_forcing_allowed"] is False,
        "free residual forcing shortcut enabled",
    )
    _require(
        cr["amplitude_collapse_success_allowed"] is False,
        "amplitude-collapse shortcut enabled",
    )

    return {
        "status": "pass",
        "contract_id": data["contract_id"],
        "upstream_pr": upstream["pull_request"],
        "candidate_sha256": pair["candidate_file_sha256"],
        "validation_sha256": pair["validation_file_sha256"],
        "candidate_validation_sha_bound_in_exact_artifact": True,
        "builder_cross_file_binding_gap_recorded": True,
        "complete_candidate_save_load_ready": False,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_contract(), indent=2, sort_keys=True))
