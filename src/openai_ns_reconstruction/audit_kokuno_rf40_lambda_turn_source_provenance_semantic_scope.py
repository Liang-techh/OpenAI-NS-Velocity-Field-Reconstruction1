"""Fail-closed CR002 audit for RF40 lambda-turn source-provenance identity scope.

This module is governance only. It does not alter or evaluate a velocity field.
It machine-locks the distinction exposed by A1 #998: the candidate has a
callable/save-load runtime semantic identity, while the external corrected-reader
provenance tuple is reported/pinned outside that semantic hash payload.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-rf40-lambda-turn-source-provenance-semantic-scope-v1"
TASK_ID = "CR002-KOKUNO-RF40-LAMBDA-TURN-SOURCE-PROVENANCE-IDENTITY-103"
BASE_HEAD = "bb3d943461820fe43c0bdb2063cb11b428c4ea1c"
A1_HEAD = "43b295444b1e9558222d757cb551e04385eddcb5"
A1_SOURCE_BLOB = "07743c30978360e305a8863e05e7ed322d818b32"
CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

SOURCE_PROVENANCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "path": "navier-stokes/navier_stokes_workbench.tex",
    "git_blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "corrected_release_date": "2026-09-09",
    "zenodo_record": "22678406",
    "classification": "public_corrected_reader_provenance",
    "original_openai_paper_exact_correspondence_independently_verified_here": False,
}

_REQUIRED_TRUE = (
    "a1_998_lambda_turn_leading_callable_materialized",
    "a1_998_configuration_save_load_available",
    "a1_998_runtime_semantic_sha256_available",
    "a1_998_report_carries_external_source_provenance",
    "a1_998_dedicated_ci_pins_external_source_blob",
    "a1_998_semantic_sha_binds_schema_source_formulas_numerical_realization_truth_boundary_configuration",
)

_REQUIRED_FALSE = (
    "a1_998_configuration_binds_external_source_provenance",
    "a1_998_semantic_sha_binds_external_source_provenance",
    "ci_source_blob_pin_implies_serialized_semantic_source_binding",
    "report_source_block_implies_serialized_semantic_source_binding",
    "source_provenance_semantic_identity_closed",
    "source_provenance_semantic_identity_closed_implies_paper_exact",
    "source_provenance_semantic_identity_closed_implies_openai_field_identified",
    "lambda_turn_callable_save_load_implies_lambda_turn_composite_materialized",
    "lambda_turn_callable_save_load_implies_global_project_domain_totality",
    "lambda_turn_callable_save_load_implies_pde_validation",
)

_FORBIDDEN_FALSE = (
    "kokuno_lambda_turn_velocity_export_ready",
    "kokuno_pde_validated",
    "kokuno_paper_exact",
    "kokuno_openai_field_identified",
    "kokuno_blowup_proved",
    "residual_defined_free_force_allowed",
    "candidate_collapse_allowed",
    "post_hoc_threshold_relaxation_allowed",
)

_EXPECTED_SEMANTIC_FIELDS = [
    "schema",
    "source_formulas",
    "numerical_realization",
    "truth_boundary",
    "configuration",
]


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def load_contract(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def semantic_scope_witness() -> dict[str, Any]:
    """Return a mechanics-only witness for the identity-boundary distinction.

    The artifact payload mirrors only the categories observed in A1 #998's
    semantic hash. The toy values are not Kokuno/OpenAI parameters or candidate
    data. Mutating source metadata outside the payload leaves the unbound digest
    unchanged, whereas binding the source block makes the mutation visible.
    """

    artifact_payload = {
        "schema": "toy-artifact",
        "source_formulas": {"formula": "toy"},
        "numerical_realization": {"choice": "toy"},
        "truth_boundary": {"paper_exact": False},
        "configuration": {"parameter": 1},
    }
    source_a = {
        "repository": "public/source",
        "commit": "a" * 40,
        "path": "reader.tex",
        "git_blob": "b" * 40,
        "record": "1",
    }
    source_b = copy.deepcopy(source_a)
    source_b["commit"] = "c" * 40
    source_b["git_blob"] = "d" * 40

    unbound_a = _sha256(artifact_payload)
    unbound_b = _sha256(artifact_payload)
    bound_a = _sha256({"artifact": artifact_payload, "source_provenance": source_a})
    bound_b = _sha256({"artifact": artifact_payload, "source_provenance": source_b})
    return {
        "classification": "autonomous_mechanics_only",
        "unbound_digest_a": unbound_a,
        "unbound_digest_b": unbound_b,
        "bound_digest_a": bound_a,
        "bound_digest_b": bound_b,
        "unbound_digest_detects_source_mutation": unbound_a != unbound_b,
        "bound_digest_detects_source_mutation": bound_a != bound_b,
        "not_a_public_source_fact": True,
        "not_candidate_numerical_evidence": True,
    }


def audit_scope(contract: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []

    if contract.get("schema") != SCHEMA:
        errors.append("schema drift")
    if contract.get("task_id") != TASK_ID:
        errors.append("task_id drift")

    base = contract.get("exact_base", {})
    if base.get("pr") != 1002 or base.get("head") != BASE_HEAD:
        errors.append("exact A5 #1002 base identity drift")

    active = contract.get("active_constrained_integration_observed", {})
    if active.get("branch") != "codex/cr001-constraints":
        errors.append("active constrained integration branch drift")

    upstream = contract.get("upstream_scope", {})
    a1 = upstream.get("agent1_998", {})
    if a1.get("pr") != 998 or a1.get("head") != A1_HEAD:
        errors.append("A1 #998 identity drift")
    if a1.get("source_blob") != A1_SOURCE_BLOB:
        errors.append("A1 #998 source blob drift")
    if a1.get("semantic_payload_fields_observed") != _EXPECTED_SEMANTIC_FIELDS:
        errors.append("A1 #998 semantic payload field boundary drift")
    if a1.get("report_carries_external_source_provenance") is not True:
        errors.append("A1 #998 report provenance fact lost")
    if a1.get("configuration_carries_external_source_provenance") is not False:
        errors.append("A1 #998 configuration provenance boundary promoted")
    if a1.get("semantic_sha256_carries_external_source_provenance") is not False:
        errors.append("A1 #998 semantic provenance boundary promoted")

    a5 = upstream.get("agent5_1002", {})
    if a5.get("pr") != 1002 or a5.get("head") != BASE_HEAD:
        errors.append("A5 #1002 registration identity drift")
    if a5.get("agent1_998_recorded_as_non_consumed_newer_sibling") is not True:
        errors.append("A1 #998 non-consumed sibling boundary lost")
    if a5.get("lambda_turn_composite_materialized") is not False:
        errors.append("lambda-turn composite prematurely promoted")
    if a5.get("lambda_turn_independent_audit_registered") is not False:
        errors.append("lambda-turn independent audit prematurely promoted")

    if contract.get("corrected_public_source_provenance") != SOURCE_PROVENANCE:
        errors.append("corrected-reader provenance tuple drift")

    locked = contract.get("machine_locked_distinctions", {})
    for key in _REQUIRED_TRUE:
        if locked.get(key) is not True:
            errors.append(f"{key} must remain true")
    for key in _REQUIRED_FALSE:
        if locked.get(key) is not False:
            errors.append(f"{key} must remain false")

    witness = contract.get("semantic_scope_witness", {})
    if witness.get("classification") != "autonomous_mechanics_only":
        errors.append("semantic witness classification drift")
    if witness.get("not_a_public_source_fact") is not True:
        errors.append("semantic witness may not become a public-source fact")
    if witness.get("not_candidate_numerical_evidence") is not True:
        errors.append("semantic witness may not become candidate numerical evidence")

    reqs = contract.get("future_promotion_requirements", {})
    source_reqs = reqs.get("source_provenance_semantic_identity_closed", [])
    joined = " ".join(str(item) for item in source_reqs)
    for required_term in ("repository", "commit", "source path", "source blob"):
        if required_term not in joined:
            errors.append(f"future provenance closure requirement missing {required_term!r}")
    if "recompute" not in joined or "fail closed" not in joined:
        errors.append("future provenance closure must require replay and fail-closed mutation detection")
    if "paper exactness" not in joined or "OpenAI-field identity" not in joined:
        errors.append("future provenance closure must preserve paper/OpenAI claim separation")

    cr001 = contract.get("cr001_snapshot", {})
    expected_cr001 = {
        "canonical_constraints_path": "configs/constraints.json",
        "canonical_constraints_blob": CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max_threshold": 0.001,
        "momentum_l2_threshold": 0.001,
        "divergence_max_threshold": 1e-5,
        "divergence_l2_threshold": 1e-5,
        "residual_defined_free_force_allowed": False,
        "candidate_collapse_allowed": False,
        "post_hoc_threshold_relaxation_allowed": False,
    }
    if cr001 != expected_cr001:
        errors.append("canonical CR001 snapshot drift")

    delivery = contract.get("canonical_delivery_independent_state", {})
    expected_delivery = {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_api": "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "kokuno_scope_may_downgrade_eq45_velocity_export_ready": False,
    }
    if delivery != expected_delivery:
        errors.append("canonical Eq45 delivery/state separation drift")

    forbidden = contract.get("forbidden_promotions", {})
    for key in _FORBIDDEN_FALSE:
        if forbidden.get(key) is not False:
            errors.append(f"forbidden promotion enabled: {key}")

    witness_result = semantic_scope_witness()
    if witness_result["unbound_digest_detects_source_mutation"]:
        errors.append("mechanics witness unexpectedly detects source mutation in unbound digest")
    if not witness_result["bound_digest_detects_source_mutation"]:
        errors.append("mechanics witness failed to detect source mutation after provenance binding")

    return errors


def assert_scope(contract: Mapping[str, Any]) -> None:
    errors = audit_scope(contract)
    if errors:
        raise ValueError(
            "CR002 RF40 lambda-turn source-provenance scope violation: "
            + "; ".join(errors)
        )


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "contract",
        nargs="?",
        default="configs/kokuno_rf40_lambda_turn_source_provenance_semantic_scope.json",
    )
    args = parser.parse_args()
    contract = load_contract(args.contract)
    errors = audit_scope(contract)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    witness = semantic_scope_witness()
    print("PASS: CR002 RF40 lambda-turn source-provenance semantic scope remains fail-closed")
    print(json.dumps(witness, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
