"""Fail-closed CR002 audit for RF40 axial-shutdown evidence scope.

This module governs representation/evidence semantics only.  It does not validate
Navier--Stokes, visual correspondence, paper exactness, or OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

SCHEMA = "cr002-kokuno-rf40-axial-shutdown-independent-audit-scope-v1"
TASK_ID = "CR002-KOKUNO-RF40-AXIAL-SHUTDOWN-INDEPENDENT-AUDIT-SCOPE-102"
CONTRACT_REL = Path("configs/kokuno_rf40_axial_shutdown_independent_audit_scope.json")
CANONICAL_CONSTRAINTS_REL = Path("configs/constraints.json")
CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1 = "6c559e42895a606e2ef025ade4cb448966d75814"

EXPECTED_BASE = {
    "pr": 996,
    "head": "88fdafc03ab71e3f8ffe7f39add900e548627daf",
    "branch": "codex/kokuno-a5-rf40-first-turn-composite-ingest-097",
    "source_blob": "07cd42a7f393bf34b918090630a2a263c2610289",
}

EXPECTED_HEADS = {
    ("upstream_scope", "agent1_993", "pr"): 993,
    ("upstream_scope", "agent1_993", "head"): "2ac6460b483efb1c07f2fa65e7fed781a32f2718",
    ("upstream_scope", "agent1_993", "source_blob"): "057a0514c8a941c3e59922b8158f480434b4441e",
    ("upstream_scope", "agent2_992", "pr"): 992,
    ("upstream_scope", "agent2_992", "head"): "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377",
    ("upstream_scope", "agent4_995", "pr"): 995,
    ("upstream_scope", "agent4_995", "head"): "70cd44ed0b5cabbcf5103df8959ff008377d17fa",
    ("upstream_scope", "agent4_995", "source_blob"): "3881a0bb98e6369cf8fde85556c9d7ac0c2db7fc",
    ("upstream_scope", "agent5_996", "pr"): 996,
    ("upstream_scope", "agent5_996", "head"): "88fdafc03ab71e3f8ffe7f39add900e548627daf",
    ("upstream_scope", "agent5_996", "source_blob"): "07cd42a7f393bf34b918090630a2a263c2610289",
}

EXPECTED_BOOLS = {
    ("upstream_scope", "agent1_993", "deterministic_save_load_identity"): True,
    ("upstream_scope", "agent1_993", "focused_production_and_source_coordinate_checks"): True,
    ("upstream_scope", "agent1_993", "implementation_distinct_axial_shutdown_cartesian_public_velocity_audit"): False,
    ("upstream_scope", "agent1_993", "consumed_by_agent2_992"): False,
    ("upstream_scope", "agent1_993", "consumed_by_agent3_994"): False,
    ("upstream_scope", "agent1_993", "covered_by_agent4_995"): False,
    ("upstream_scope", "agent2_992", "fails_closed_after_X_1"): True,
    ("upstream_scope", "agent2_992", "axial_shutdown_composite_materialized"): False,
    ("upstream_scope", "agent4_995", "strictly_post_xr_first_turn_probes"): True,
    ("upstream_scope", "agent4_995", "axial_shutdown_probes_X1_to_X2"): False,
    ("upstream_scope", "agent4_995", "uses_only_save_reloaded_public_cartesian_velocity_for_numerical_differentiation"): True,
    ("upstream_scope", "agent4_995", "uses_agent1_source_derivative_helpers"): False,
    ("upstream_scope", "agent4_995", "scoped_divergence_evidence_only"): True,
    ("upstream_scope", "agent4_995", "canonical_whole_domain_admission"): False,
    ("upstream_scope", "agent4_995", "momentum_or_full_ns_evidence"): False,
    ("upstream_scope", "agent5_996", "agent1_993_recorded_as_non_consumed_sibling"): True,
    ("upstream_scope", "agent5_996", "first_turn_composite_registered"): True,
    ("upstream_scope", "agent5_996", "first_turn_independent_audit_registered"): True,
    ("upstream_scope", "agent5_996", "leading_plus_oscillatory_velocity_through_axial_shutdown_materialized"): False,
    ("upstream_scope", "agent5_996", "independent_axial_shutdown_composite_audit_available"): False,
    ("upstream_scope", "agent5_996", "velocity_export_ready"): False,
    ("upstream_scope", "agent5_996", "pde_validated"): False,
    ("machine_locked_distinctions", "agent1_993_axial_shutdown_leading_only_callable_materialized"): True,
    ("machine_locked_distinctions", "agent1_993_has_identity_preserving_save_load"): True,
    ("machine_locked_distinctions", "agent1_993_has_focused_production_and_source_coordinate_checks"): True,
    ("machine_locked_distinctions", "agent4_995_first_turn_cartesian_public_velocity_audit_exists"): True,
    ("machine_locked_distinctions", "agent4_995_covers_axial_shutdown_X1_to_X2"): False,
    ("machine_locked_distinctions", "axial_shutdown_leading_independent_cartesian_public_velocity_audit_available"): False,
    ("machine_locked_distinctions", "first_turn_a4_995_audit_may_be_relabelled_as_axial_shutdown_audit"): False,
    ("machine_locked_distinctions", "source_coordinate_handoff_or_derivative_replay_may_be_relabelled_as_independent_axial_shutdown_cartesian_validation"): False,
    ("machine_locked_distinctions", "axial_shutdown_leading_materialization_implies_axial_shutdown_composite_materialization"): False,
    ("machine_locked_distinctions", "axial_shutdown_leading_materialization_implies_global_project_domain_totality"): False,
    ("machine_locked_distinctions", "axial_shutdown_leading_materialization_implies_pde_validation"): False,
    ("machine_locked_distinctions", "first_turn_composite_save_load_identity_implies_axial_shutdown_composite_save_load_identity"): False,
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
    "A4 #995 first-turn audit -> axial-shutdown audit",
    "A1 #993 production/source-coordinate checks -> implementation-distinct axial-shutdown Cartesian audit",
    "A1 #993 leading-only axial-shutdown callable -> leading+oscillatory axial-shutdown composite",
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
        "residual_defined_free_forcing_forbidden": "No residual-dependent basis or pointwise free force" in canonical["forcing"]["restriction"],
        "candidate_collapse_forbidden": "reject collapsed candidates" in canonical["nontriviality"]["enforcement"],
        "post_hoc_threshold_relaxation_forbidden": "changing thresholds requires a new experiment version" in canonical["validation"]["failure_policy"],
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

    provenance = contract.get("four_way_provenance")
    if not isinstance(provenance, Mapping):
        errors.append("missing:four_way_provenance")
    else:
        expected_buckets = {"user_requirements", "public_source_facts", "autonomous_repository_facts", "pending_or_unknown"}
        if set(provenance) != expected_buckets:
            errors.append(f"mismatch:four_way_provenance.keys:{sorted(provenance)}")
        for key in sorted(expected_buckets):
            value = provenance.get(key)
            if not isinstance(value, list) or not value or not all(isinstance(item, str) and item.strip() for item in value):
                errors.append(f"invalid:four_way_provenance.{key}")

    for path, expected in EXPECTED_HEADS.items():
        _expect(errors, contract, path, expected)
    for path, expected in EXPECTED_BOOLS.items():
        _expect(errors, contract, path, expected)

    _expect(errors, contract, ("upstream_scope", "agent1_993", "materialized_domain"), "leading-only through RF40 axial shutdown X_2")
    _expect(errors, contract, ("upstream_scope", "agent2_992", "materialized_domain"), "leading + frozen complete-curl oscillation through X_1=e X_R")
    _expect(errors, contract, ("upstream_scope", "agent4_995", "audited_domain"), "through RF40 first turn X_1=e X_R")
    _expect(errors, contract, ("scope_logic_witness", "classification"), "autonomous_mechanics_only")

    for promotion_key in (
        "axial_shutdown_leading_independent_cartesian_public_velocity_audit_available",
        "axial_shutdown_composite_ready",
    ):
        items = contract.get("future_promotion_requirements", {}).get(promotion_key)
        if not isinstance(items, list) or len(items) < 3 or not all(isinstance(x, str) and x.strip() for x in items):
            errors.append(f"invalid:future_promotion_requirements.{promotion_key}")

    actual_cr001 = contract.get("cr001_snapshot")
    if actual_cr001 != EXPECTED_CR001:
        errors.append("mismatch:cr001_snapshot")

    if canonical_constraints is not None:
        if canonical_blob_sha1 is None:
            errors.append("missing:canonical_blob_sha1")
        else:
            replay = _canonical_snapshot(canonical_constraints, canonical_blob_sha1)
            if replay != EXPECTED_CR001:
                errors.append("mismatch:canonical_constraints_replay")
            if actual_cr001 != replay:
                errors.append("mismatch:contract_vs_canonical_constraints")

    _expect(
        errors,
        contract,
        ("canonical_delivery_independence", "candidate_family"),
        "eq45_supported_velocity_candidate_v1",
    )
    _expect(
        errors,
        contract,
        ("canonical_delivery_independence", "velocity_api"),
        "openai_ns_reconstruction.eq45_supported_delivery:velocity",
    )

    forbidden = contract.get("forbidden_promotions")
    if not isinstance(forbidden, list):
        errors.append("invalid:forbidden_promotions")
    else:
        missing = REQUIRED_FORBIDDEN.difference(forbidden)
        if missing:
            errors.append("missing:forbidden_promotions:" + "|".join(sorted(missing)))

    statement = contract.get("scope_logic_witness", {}).get("statement", "")
    if "X_1 < X <= X_2" not in statement or "X <= X_1" not in statement:
        errors.append("invalid:scope_logic_witness.statement")

    return errors


def audit_default() -> list[str]:
    contract = load_contract()
    canonical, blob_sha1 = load_canonical_constraints()
    return audit_scope(contract, canonical, blob_sha1)


def receipt(contract: Mapping[str, Any], errors: Sequence[str]) -> dict[str, Any]:
    return {
        "schema": "cr002-kokuno-rf40-axial-shutdown-independent-audit-scope-receipt-v1",
        "task_id": TASK_ID,
        "ok": not errors,
        "errors": list(errors),
        "exact_base_head": contract.get("exact_base", {}).get("head"),
        "axial_shutdown_leading_materialized": contract.get("machine_locked_distinctions", {}).get(
            "agent1_993_axial_shutdown_leading_only_callable_materialized"
        ),
        "axial_shutdown_independent_public_cartesian_audit_available": contract.get(
            "machine_locked_distinctions", {}
        ).get("axial_shutdown_leading_independent_cartesian_public_velocity_audit_available"),
        "canonical_eq45_velocity_export_ready": contract.get("canonical_delivery_independence", {}).get(
            "velocity_export_ready"
        ),
        "kokuno_velocity_export_ready": contract.get("kokuno_route_state", {}).get("velocity_export_ready"),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=None)
    parser.add_argument("--constraints", type=Path, default=None)
    args = parser.parse_args(argv)

    contract = load_contract(args.contract)
    canonical, blob_sha1 = load_canonical_constraints(args.constraints)
    errors = audit_scope(contract, canonical, blob_sha1)
    print(json.dumps(receipt(contract, errors), sort_keys=True))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
