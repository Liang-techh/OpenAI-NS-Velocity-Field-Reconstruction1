"""Fail-closed CR002 audit for exact A1 #1107 log-X representation promotion.

This audits representation/evidence scope only. Construction replay, save/load,
CI, or source-form agreement never imply PDE validation, visual correspondence,
paper exactness, or exact OpenAI-field identity.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT_PATH = Path("configs/kokuno_main_pulse_logx_representation_promotion_scope.json")
CHILD_PATH = Path("src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_main_pulse_logx.py")
PARENT_PATH = Path("src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_main_pulse.py")
CONSTRAINTS_PATH = Path("configs/constraints.json")
STATUS_PATH = Path("project_status.json")

EXPECTED_BASE_HEAD = "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
EXPECTED_PARENT_HEAD = "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156"
EXPECTED_CHILD_SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-logx-v1"
EXPECTED_PARENT_SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-v1"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

_CHILD_TRUE = {
    "overflow_safe_logX_similarity_materialized",
    "overflow_safe_logF_cartesian_swirl_materialized",
    "full_source_xi_11_current_cartesian_materialized",
}
_CHILD_FALSE = {
    "source_exact_amplitude_root_materialized",
    "source_pulse_end_MJ_corrections_materialized",
    "source_hidden_parameters_recovered",
    "source_terminal_tail_schedule_bound_into_current_velocity",
    "source_exterior_heat_replacement_materialized",
    "outer_global_leading_velocity_materialized",
    "unified_global_cartesian_velocity_export_ready",
    "matched_global_pressure_materialized",
    "restricted_forcing_materialized",
    "heldout_ns_residual_assessed",
    "same_protocol_comparable_to_st006",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}
_REQUIRED_METHODS = {
    "velocity", "configuration", "from_configuration", "save_configuration",
    "load_configuration", "semantic_sha256", "truth_boundary", "receipt",
}


def _read_json(path: Path) -> dict[str, Any]:
    out = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(out, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return out


def _git_blob_sha1(data: bytes) -> str:
    prefix = b"blob " + str(len(data)).encode() + b"\0"
    return hashlib.sha1(prefix + data).hexdigest()


def _literals_and_methods(path: Path, class_name: str) -> tuple[dict[str, Any], set[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    literals: dict[str, Any] = {}
    methods: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                literals[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            methods.update(
                item.name for item in node.body
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
            )
    return literals, methods


def representation_witness() -> dict[str, Any]:
    """Synthetic mechanics witness only; never source/OpenAI/PDE evidence."""
    float_max = float.fromhex("0x1.fffffffffffffp+1023")
    log_max = math.log(float_max)
    log_x = log_max + 100.0
    try:
        explicit_x_ok = math.isfinite(math.exp(log_x))
    except OverflowError:
        explicit_x_ok = False
    # X=r^2/(2q), q=1: unrepresentable X can coexist with finite radius.
    log_r = 0.5 * (log_x + math.log(2.0))
    return {
        "classification": "autonomous_mechanics_only",
        "log_X": log_x,
        "log_float64_max": log_max,
        "explicit_X_float64_representable": explicit_x_ok,
        "physical_radius_float64_finite": math.isfinite(math.exp(log_r)),
    }


def audit_contract(repo_root: str | Path, contract: Mapping[str, Any] | None = None) -> list[str]:
    root = Path(repo_root)
    errors: list[str] = []
    try:
        c = dict(contract) if contract is not None else _read_json(root / CONTRACT_PATH)
    except Exception as exc:
        return [f"contract unreadable: {exc}"]

    if c.get("schema") != "cr002-kokuno-logx-representation-promotion-v1":
        errors.append("unexpected contract schema")
    stack = c.get("exact_stack", {})
    if stack.get("base_head") != EXPECTED_BASE_HEAD:
        errors.append("exact #1107 base head drifted")
    if stack.get("parent_finite_x_head") != EXPECTED_PARENT_HEAD:
        errors.append("exact #1100 parent head drifted")

    provenance = c.get("four_way_provenance", {})
    classes = {"user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"}
    if set(provenance) != classes or any(not provenance.get(k) for k in classes):
        errors.append("four-way provenance incomplete or reclassified")

    transition = c.get("representation_transition", {})
    expected_transition = {
        "parent_schema": EXPECTED_PARENT_SCHEMA,
        "child_schema": EXPECTED_CHILD_SCHEMA,
        "new_representation_identity_required": True,
        "parent_and_child_same_semantic_candidate_identity": False,
        "parent_finite_x_gap_closed_by_child_logx_through_public_principal_xi11": True,
        "child_public_cartesian_velocity_finite_on_declared_xi11_domain": True,
        "child_save_load_semantic_identity_present": True,
        "source_exact_main_pulse_recovered": False,
    }
    for key, value in expected_transition.items():
        if transition.get(key) != value:
            errors.append(f"representation transition drift: {key}")

    transfer = c.get("evidence_transfer", {})
    for key in (
        "parent_finite_x_receipts_prove_new_logx_only_interval",
        "a2_pr1108_composite_consumes_child_pr1107",
        "a4_pr1110_divergence_audits_child_pr1107",
        "a3_pr1109_correction_evidence_belongs_to_child_pr1107",
    ):
        if transfer.get(key) is not False:
            errors.append(f"forbidden evidence transfer enabled: {key}")
    for key in (
        "fresh_matching_a2_composition_required_for_child",
        "fresh_implementation_distinct_a4_audit_required_for_child",
    ):
        if transfer.get(key) is not True:
            errors.append(f"fresh downstream evidence requirement weakened: {key}")

    states = c.get("truth_states", {})
    if states.get("logx_leading_only_public_principal_xi11_materialized") is not True:
        errors.append("log-X xi=11 leading representation not registered")
    for key in (
        "source_exact_amplitude_root_materialized", "source_pulse_end_MJ_corrections_materialized",
        "outer_global_leading_velocity_materialized", "kokuno_velocity_export_ready",
        "visual_correspondence_verified", "pde_validated", "paper_exact",
        "openai_field_identified", "blowup_proved",
    ):
        if states.get(key) is not False:
            errors.append(f"truth over-promoted: {key}")

    try:
        child, methods = _literals_and_methods(root / CHILD_PATH, "KokunoPA16CurrentCartesianMainPulseLogX")
        parent, _ = _literals_and_methods(root / PARENT_PATH, "KokunoPA16CurrentCartesianMainPulse")
    except Exception as exc:
        errors.append(f"representation implementation unreadable: {exc}")
    else:
        if child.get("SCHEMA") != EXPECTED_CHILD_SCHEMA:
            errors.append("live child schema drifted")
        if child.get("PARENT_EXACT_HEAD") != EXPECTED_PARENT_HEAD:
            errors.append("live child parent binding drifted")
        if parent.get("SCHEMA") != EXPECTED_PARENT_SCHEMA:
            errors.append("live parent schema drifted")
        if child.get("SCHEMA") == parent.get("SCHEMA"):
            errors.append("new representation did not receive a distinct schema")
        truth = child.get("_TRUTH_UPDATES")
        if not isinstance(truth, dict):
            errors.append("child truth boundary is not literal/auditable")
        else:
            for key in _CHILD_TRUE:
                if truth.get(key) is not True:
                    errors.append(f"child representation fact missing: {key}")
            for key in _CHILD_FALSE:
                if truth.get(key) is not False:
                    errors.append(f"child scientific boundary weakened: {key}")
        missing = _REQUIRED_METHODS - methods
        if missing:
            errors.append(f"child callable/save-load/identity API missing: {sorted(missing)}")
        if "np.isfinite(out)" not in (root / CHILD_PATH).read_text(encoding="utf-8"):
            errors.append("finite Cartesian output fail-close check missing")

    try:
        raw = (root / CONSTRAINTS_PATH).read_bytes()
        constraints = json.loads(raw)
    except Exception as exc:
        errors.append(f"CR001 constraints unreadable: {exc}")
    else:
        if _git_blob_sha1(raw) != EXPECTED_CONSTRAINTS_BLOB:
            errors.append("canonical CR001 constraints blob changed")
        cr = c.get("canonical_cr001", {})
        observed = {
            "nu": constraints.get("nu"),
            "physical_domain": constraints.get("domain", {}).get("physical"),
            "evaluation_box": constraints.get("domain", {}).get("evaluation_box"),
            "support": constraints.get("domain", {}).get("support"),
            "time_interval": constraints.get("domain", {}).get("time_interval"),
            "forcing_mode": constraints.get("forcing", {}).get("mode"),
            "reference_energy": constraints.get("nontriviality", {}).get("reference_energy"),
            "reference_energy_abs_tolerance": constraints.get("nontriviality", {}).get("reference_energy_abs_tolerance"),
            "validation_seed": constraints.get("validation", {}).get("seed"),
            "held_out_points": constraints.get("validation", {}).get("held_out_points"),
            "derivative_steps": constraints.get("validation", {}).get("derivative_steps"),
            "quadrature_orders_per_axis": constraints.get("validation", {}).get("quadrature_orders_per_axis"),
            "momentum_max": constraints.get("validation", {}).get("thresholds", {}).get("pde_residual_max"),
            "momentum_L2": constraints.get("validation", {}).get("thresholds", {}).get("pde_residual_L2"),
            "divergence_max": constraints.get("validation", {}).get("thresholds", {}).get("divergence_max"),
            "divergence_L2": constraints.get("validation", {}).get("thresholds", {}).get("divergence_L2"),
        }
        if cr.get("constraints_git_blob_sha1") != EXPECTED_CONSTRAINTS_BLOB:
            errors.append("CR001 blob pin weakened")
        for key, value in observed.items():
            if cr.get(key) != value:
                errors.append(f"CR001 contract drift: {key}")
        for key in ("residual_defined_free_forcing_forbidden", "candidate_collapse_forbidden", "post_hoc_threshold_relaxation_forbidden"):
            if cr.get(key) is not True:
                errors.append(f"CR001 prohibition weakened: {key}")

    try:
        status = _read_json(root / STATUS_PATH)
    except Exception as exc:
        errors.append(f"project status unreadable: {exc}")
    else:
        expected = c.get("canonical_delivery", {})
        live = status.get("states", {})
        if status.get("candidate_family") != expected.get("candidate_family"):
            errors.append("canonical delivery candidate family drifted")
        for key in ("velocity_export_ready", "visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
            if live.get(key) != expected.get(key):
                errors.append(f"canonical Eq45 state drift: {key}")

    witness = representation_witness()
    if witness["classification"] != "autonomous_mechanics_only":
        errors.append("mechanics witness misclassified")
    if witness["explicit_X_float64_representable"] is not False:
        errors.append("mechanics witness did not exceed explicit-X float64 range")
    if witness["physical_radius_float64_finite"] is not True:
        errors.append("mechanics witness did not preserve finite physical radius")
    return errors


def assert_contract(repo_root: str | Path, contract: Mapping[str, Any] | None = None) -> None:
    errors = audit_contract(repo_root, contract)
    if errors:
        raise AssertionError("\n".join(errors))


__all__ = ["audit_contract", "assert_contract", "representation_witness"]
