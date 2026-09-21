"""Fail-closed CR002 audit for RF40 power-law composite/audit scope.

Governance only: this module does not validate Navier--Stokes, visual
correspondence, paper exactness, OpenAI-field identity, or blow-up.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

SCHEMA = "cr002-kokuno-rf40-power-law-composite-audit-scope-v1"
TASK_ID = "CR002-KOKUNO-RF40-POWER-LAW-COMPOSITE-AUDIT-SCOPE-104"
CONTRACT_REL = Path("configs/kokuno_rf40_power_law_composite_audit_scope.json")
CANONICAL_CONSTRAINTS_REL = Path("configs/constraints.json")
CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1 = "6c559e42895a606e2ef025ade4cb448966d75814"

EXPECTED_BASE = {
    "pr": 1008,
    "head": "28860ee1f4655670c62b24880f08aa3d50459623",
    "branch": "codex/kokuno-a5-rf40-lambda-turn-composite-ingest-099",
    "source_blob": "c71499e4e0ae38c9f6eb1032c2e63be7dde58894",
}
EXPECTED_ACTIVE_INTEGRATION = {
    "branch": "codex/cr001-constraints",
    "head": "7f6658cac1a8ffea178bec5accf2fa5a784b8249",
    "role": "live constrained truth/routing reference; this narrow Kokuno governance increment is stacked on exact A5 #1008 to avoid retargeting sibling evidence",
}

EXPECTED_IDENTITIES = {
    ("upstream_scope", "agent1_1005", "pr"): 1005,
    ("upstream_scope", "agent1_1005", "head"): "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
    ("upstream_scope", "agent1_1005", "source_blob"): "36a0183352a2005bd0b5f477594e05a2e4fe3868",
    ("upstream_scope", "agent2_1004", "pr"): 1004,
    ("upstream_scope", "agent2_1004", "head"): "e9697cdfe2249e1b8bc2f9241bd16d456c60418f",
    ("upstream_scope", "agent2_1004", "source_blob"): "def2998256b2339115c2400c75b6db6ab58dfe95",
    ("upstream_scope", "agent4_1007", "pr"): 1007,
    ("upstream_scope", "agent4_1007", "head"): "8b5b337940edd568bf30ddba7afe68815bce1fd9",
    ("upstream_scope", "agent4_1007", "source_blob"): "2ae59ab9bc8480979f0b0660b1db8b801fe847c1",
    ("upstream_scope", "agent5_1008", "pr"): 1008,
    ("upstream_scope", "agent5_1008", "head"): "28860ee1f4655670c62b24880f08aa3d50459623",
    ("upstream_scope", "agent5_1008", "source_blob"): "c71499e4e0ae38c9f6eb1032c2e63be7dde58894",
}

EXPECTED_BOOLS = {
    ("upstream_scope", "agent1_1005", "deterministic_configuration_save_load"): True,
    ("upstream_scope", "agent1_1005", "fails_closed_after_X_4"): True,
    ("upstream_scope", "agent1_1005", "power_law_leading_materialized"): True,
    ("upstream_scope", "agent1_1005", "power_law_composite_materialized"): False,
    ("upstream_scope", "agent1_1005", "implementation_distinct_power_law_composite_audit_available"): False,
    ("upstream_scope", "agent1_1005", "global_complete_velocity_materialized"): False,
    ("upstream_scope", "agent1_1005", "pde_validated"): False,
    ("upstream_scope", "agent2_1004", "identity_preserving_save_load"): True,
    ("upstream_scope", "agent2_1004", "full_concrete_oscillatory_runtime_digest_bound"): True,
    ("upstream_scope", "agent2_1004", "fails_closed_after_X_3"): True,
    ("upstream_scope", "agent2_1004", "power_law_composite_materialized"): False,
    ("upstream_scope", "agent4_1007", "strict_power_law_probes_X3_to_X4"): False,
    ("upstream_scope", "agent4_1007", "uses_only_save_reloaded_public_cartesian_velocity_for_numerical_differentiation"): True,
    ("upstream_scope", "agent4_1007", "uses_agent1_source_derivative_helpers"): False,
    ("upstream_scope", "agent4_1007", "scoped_divergence_evidence_only"): True,
    ("upstream_scope", "agent4_1007", "canonical_whole_domain_admission"): False,
    ("upstream_scope", "agent4_1007", "momentum_or_full_ns_evidence"): False,
    ("upstream_scope", "agent5_1008", "lambda_turn_composite_registered"): True,
    ("upstream_scope", "agent5_1008", "lambda_turn_independent_audit_registered"): True,
    ("upstream_scope", "agent5_1008", "agent1_1005_power_law_recorded_as_non_consumed_sibling"): True,
    ("upstream_scope", "agent5_1008", "leading_plus_oscillatory_velocity_through_power_law_materialized"): False,
    ("upstream_scope", "agent5_1008", "independent_power_law_composite_audit_available"): False,
    ("upstream_scope", "agent5_1008", "complete_post_xr_rf40_current_lineage_materialized"): False,
    ("upstream_scope", "agent5_1008", "velocity_export_ready"): False,
    ("upstream_scope", "agent5_1008", "pde_validated"): False,
    ("machine_locked_distinctions", "agent1_1005_power_law_leading_only_callable_materialized"): True,
    ("machine_locked_distinctions", "agent1_1005_has_deterministic_configuration_save_load"): True,
    ("machine_locked_distinctions", "agent2_1004_lambda_turn_composite_exists"): True,
    ("machine_locked_distinctions", "agent4_1007_lambda_turn_public_cartesian_audit_exists"): True,
    ("machine_locked_distinctions", "agent4_1007_covers_power_law_X3_to_X4"): False,
    ("machine_locked_distinctions", "power_law_composite_materialized"): False,
    ("machine_locked_distinctions", "power_law_composite_independent_cartesian_public_velocity_audit_available"): False,
    ("machine_locked_distinctions", "lambda_turn_composite_or_audit_may_be_relabelled_as_power_law_composite_or_audit"): False,
    ("machine_locked_distinctions", "agent1_1005_production_or_source_coordinate_checks_may_be_relabelled_as_independent_power_law_composite_validation"): False,
    ("machine_locked_distinctions", "power_law_leading_materialization_implies_power_law_leading_plus_oscillatory_composite"): False,
    ("machine_locked_distinctions", "power_law_leading_materialization_implies_global_project_domain_totality"): False,
    ("machine_locked_distinctions", "power_law_leading_materialization_implies_complete_ns_defect"): False,
    ("machine_locked_distinctions", "power_law_leading_materialization_implies_pde_validation"): False,
    ("machine_locked_distinctions", "lambda_turn_composite_save_load_identity_implies_power_law_composite_save_load_identity"): False,
    ("scope_logic_witness", "not_a_public_source_fact"): True,
    ("scope_logic_witness", "not_candidate_numerical_evidence"): True,
    ("canonical_delivery_independence", "velocity_export_ready"): True,
    ("canonical_delivery_independence", "visualization_ready"): False,
    ("canonical_delivery_independence", "visual_correspondence_verified"): False,
    ("canonical_delivery_independence", "pde_validated"): False,
    ("canonical_delivery_independence", "paper_exact"): False,
    ("canonical_delivery_independence", "openai_field_identified"): False,
    ("kokuno_route_state", "velocity_export_ready"): False,
    ("kokuno_route_state", "visual_correspondence_verified"): False,
    ("kokuno_route_state", "pde_validated"): False,
    ("kokuno_route_state", "paper_exact"): False,
    ("kokuno_route_state", "openai_field_identified"): False,
    ("kokuno_route_state", "complete_candidate_api_ready"): False,
    ("scientific_changes", "velocity_or_profile_bytes_changed"): False,
    ("scientific_changes", "pressure_changed"): False,
    ("scientific_changes", "forcing_changed"): False,
    ("scientific_changes", "viscosity_changed"): False,
    ("scientific_changes", "candidate_coefficients_changed"): False,
    ("scientific_changes", "validation_samples_changed"): False,
    ("scientific_changes", "derivative_or_quadrature_protocol_changed"): False,
    ("scientific_changes", "scientific_thresholds_changed"): False,
}

EXPECTED_CR001 = {
    "constraints_blob_sha1": CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1,
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
    "momentum_max": 0.001,
    "momentum_l2": 0.001,
    "divergence_max": 0.00001,
    "divergence_l2": 0.00001,
    "residual_defined_free_forcing_forbidden": True,
    "candidate_collapse_forbidden": True,
    "post_hoc_threshold_relaxation_forbidden": True,
}

REQUIRED_FORBIDDEN = {
    "A2 #1004 lambda-turn composite -> power-law composite",
    "A4 #1007 lambda-turn audit -> power-law composite audit",
    "A1 #1005 leading-only power-law callable -> leading+oscillatory power-law composite",
    "A1 #1005 production/source-coordinate checks -> implementation-distinct power-law composite audit",
    "lambda-turn composite save/load -> power-law composite save/load",
    "scoped divergence evidence -> canonical whole-domain admission",
    "scoped divergence evidence -> momentum/full-NS evidence",
    "Kokuno route incompleteness -> downgrade of canonical Eq45 velocity_export_ready",
    "residual -> unrestricted pointwise forcing",
    "failed fixed gate -> post-hoc threshold relaxation",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _get(mapping: Mapping[str, Any], path: Sequence[str]) -> Any:
    value: Any = mapping
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            raise KeyError(".".join(path))
        value = value[key]
    return value


def _expect(errors: list[str], mapping: Mapping[str, Any], path: Sequence[str], expected: Any) -> None:
    label = ".".join(path)
    try:
        actual = _get(mapping, path)
    except KeyError:
        errors.append(f"missing:{label}")
        return
    if actual != expected:
        errors.append(f"mismatch:{label}: expected={expected!r} actual={actual!r}")


def _git_blob_sha1(raw: bytes) -> str:
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def load_contract(path: Path | None = None) -> dict[str, Any]:
    target = path or (_repo_root() / CONTRACT_REL)
    return json.loads(target.read_text(encoding="utf-8"))


def load_canonical_constraints(path: Path | None = None) -> tuple[dict[str, Any], str]:
    target = path or (_repo_root() / CANONICAL_CONSTRAINTS_REL)
    raw = target.read_bytes()
    return json.loads(raw.decode("utf-8")), _git_blob_sha1(raw)


def _canonical_snapshot(canonical: Mapping[str, Any], blob_sha1: str) -> dict[str, Any]:
    return {
        "constraints_blob_sha1": blob_sha1,
        "nu": canonical["nu"],
        "physical_domain": canonical["domain"]["physical"],
        "evaluation_box": canonical["domain"]["evaluation_box"],
        "support": canonical["domain"]["support"],
        "time_interval": canonical["domain"]["time_interval"],
        "forcing_mode": canonical["forcing"]["mode"],
        "reference_energy": canonical["nontriviality"]["reference_energy"],
        "reference_energy_abs_tolerance": canonical["nontriviality"]["reference_energy_abs_tolerance"],
        "validation_seed": canonical["validation"]["seed"],
        "held_out_points": canonical["validation"]["held_out_points"],
        "derivative_steps": canonical["validation"]["derivative_steps"],
        "quadrature_orders_per_axis": canonical["validation"]["quadrature_orders_per_axis"],
        "momentum_max": canonical["validation"]["thresholds"]["pde_residual_max"],
        "momentum_l2": canonical["validation"]["thresholds"]["pde_residual_L2"],
        "divergence_max": canonical["validation"]["thresholds"]["divergence_max"],
        "divergence_l2": canonical["validation"]["thresholds"]["divergence_L2"],
        "residual_defined_free_forcing_forbidden": (
            "No residual-dependent basis or pointwise free force" in canonical["forcing"]["restriction"]
        ),
        "candidate_collapse_forbidden": (
            "reject collapsed candidates" in canonical["nontriviality"]["enforcement"]
        ),
        "post_hoc_threshold_relaxation_forbidden": (
            "changing thresholds requires a new experiment version" in canonical["validation"]["failure_policy"]
        ),
    }


def audit_scope(
    contract: Mapping[str, Any],
    canonical_constraints: Mapping[str, Any] | None = None,
    canonical_blob_sha1: str | None = None,
) -> list[str]:
    errors: list[str] = []
    _expect(errors, contract, ("schema",), SCHEMA)
    _expect(errors, contract, ("task_id",), TASK_ID)
    _expect(errors, contract, ("exact_base",), EXPECTED_BASE)
    _expect(errors, contract, ("active_integration_reference",), EXPECTED_ACTIVE_INTEGRATION)

    provenance = contract.get("four_way_provenance")
    expected_buckets = {
        "user_requirements",
        "public_source_facts",
        "autonomous_repository_facts",
        "pending_or_unknown",
    }
    if not isinstance(provenance, Mapping):
        errors.append("missing:four_way_provenance")
    else:
        if set(provenance) != expected_buckets:
            errors.append(f"mismatch:four_way_provenance.keys:{sorted(provenance)}")
        for bucket in sorted(expected_buckets):
            values = provenance.get(bucket)
            if not isinstance(values, list) or not values or not all(
                isinstance(v, str) and v.strip() for v in values
            ):
                errors.append(f"invalid:four_way_provenance.{bucket}")

    for path, expected in EXPECTED_IDENTITIES.items():
        _expect(errors, contract, path, expected)
    for path, expected in EXPECTED_BOOLS.items():
        _expect(errors, contract, path, expected)

    _expect(
        errors,
        contract,
        ("upstream_scope", "agent1_1005", "materialized_domain"),
        "leading-only through RF40 constant-lambda power-law endpoint X_4",
    )
    _expect(
        errors,
        contract,
        ("upstream_scope", "agent2_1004", "materialized_domain"),
        "leading + frozen complete-curl oscillation through RF40 lambda turn X_3=X_w",
    )
    _expect(
        errors,
        contract,
        ("upstream_scope", "agent4_1007", "audited_domain"),
        "lambda-turn composite through X_3=X_w",
    )
    _expect(errors, contract, ("scope_logic_witness", "classification"), "autonomous_mechanics_only")

    for key, minimum in (
        ("power_law_composite_ready", 3),
        ("power_law_composite_independent_cartesian_public_velocity_audit_available", 4),
        ("correction_or_pde_promotion", 3),
    ):
        values = contract.get("future_promotion_requirements", {}).get(key)
        if not isinstance(values, list) or len(values) < minimum or not all(
            isinstance(v, str) and v.strip() for v in values
        ):
            errors.append(f"invalid:future_promotion_requirements.{key}")

    if contract.get("cr001_snapshot") != EXPECTED_CR001:
        errors.append("mismatch:cr001_snapshot")

    if canonical_constraints is not None:
        if canonical_blob_sha1 is None:
            errors.append("missing:canonical_blob_sha1")
        else:
            observed = _canonical_snapshot(canonical_constraints, canonical_blob_sha1)
            if observed != EXPECTED_CR001:
                errors.append("mismatch:live_canonical_constraints")
            if contract.get("cr001_snapshot") != observed:
                errors.append("mismatch:contract_vs_live_canonical_constraints")

    forbidden = contract.get("forbidden_promotions")
    if not isinstance(forbidden, list):
        errors.append("missing:forbidden_promotions")
    else:
        missing = sorted(REQUIRED_FORBIDDEN - set(forbidden))
        if missing:
            errors.append(f"missing:forbidden_promotions:{missing}")

    changes = contract.get("scientific_changes")
    if not isinstance(changes, Mapping) or not changes:
        errors.append("missing:scientific_changes")
    elif any(value is not False for value in changes.values()):
        errors.append("promotion:scientific_changes")

    return errors


def audit_repository() -> list[str]:
    contract = load_contract()
    canonical, blob = load_canonical_constraints()
    return audit_scope(contract, canonical, blob)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--canonical-constraints", type=Path)
    args = parser.parse_args(argv)

    contract = load_contract(args.contract)
    if args.canonical_constraints is None:
        canonical, blob = load_canonical_constraints()
    else:
        raw = args.canonical_constraints.read_bytes()
        canonical = json.loads(raw.decode("utf-8"))
        blob = _git_blob_sha1(raw)

    errors = audit_scope(contract, canonical, blob)
    if errors:
        for error in errors:
            print(error)
        return 1
    print("PASS: CR002 RF40 power-law composite/audit scope remains fail-closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
