#!/usr/bin/env python3
"""Fail-closed audit for CR002 canonical-quadrature execution provenance.

This audit does not validate a Kokuno candidate.  It verifies a narrower truth
boundary in PR #903: a checksum-valid receipt built from caller-supplied metric
rows is an integrity object, not by itself evidence that those rows came from an
independent registered 24/48/96 quadrature execution.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCOPE_PATH = ROOT / "configs/kokuno_canonical_quadrature_provenance_scope.json"


def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _find_function(tree: ast.AST, name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"missing function: {name}")


def _call_name(node: ast.Call) -> str | None:
    if isinstance(node.func, ast.Name):
        return node.func.id
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    return None


def _builder_keyword_values(builder: ast.FunctionDef) -> dict[str, ast.AST]:
    for node in ast.walk(builder):
        if isinstance(node, ast.Call) and _call_name(node) == "CanonicalQuadratureEvidence":
            return {kw.arg: kw.value for kw in node.keywords if kw.arg is not None}
    raise AssertionError("builder does not construct CanonicalQuadratureEvidence")


def _is_literal_true(node: ast.AST | None) -> bool:
    return isinstance(node, ast.Constant) and node.value is True


def _all_arg_names(fn: ast.FunctionDef) -> set[str]:
    args = fn.args
    nodes = list(args.posonlyargs) + list(args.args) + list(args.kwonlyargs)
    return {node.arg for node in nodes}


def _assert_canonical_constraints(scope: dict[str, Any], constraints: dict[str, Any]) -> None:
    expected = scope["canonical_cr001_invariants"]
    assert constraints["nu"] == expected["nu"]
    assert constraints["domain"]["physical"] == expected["physical_domain"]
    assert constraints["domain"]["evaluation_box"] == expected["evaluation_box"]
    assert constraints["domain"]["support"] == expected["support"]
    assert constraints["domain"]["time_interval"] == expected["time_interval"]
    assert constraints["forcing"]["mode"] == expected["forcing_mode"]
    assert constraints["nontriviality"]["reference_energy"] == expected["reference_energy"]
    assert (
        constraints["nontriviality"]["reference_energy_abs_tolerance"]
        == expected["reference_energy_abs_tolerance"]
    )
    validation = constraints["validation"]
    assert validation["seed"] == expected["validation_seed"]
    assert validation["held_out_points"] == expected["held_out_points"]
    assert validation["derivative_steps"] == expected["derivative_steps"]
    assert validation["quadrature_orders_per_axis"] == expected["quadrature_orders_per_axis"]
    thresholds = validation["thresholds"]
    assert thresholds["divergence_max"] == expected["divergence_max"]
    assert thresholds["divergence_L2"] == expected["divergence_L2"]
    assert thresholds["pde_residual_max"] == expected["pde_residual_max"]
    assert thresholds["pde_residual_L2"] == expected["pde_residual_L2"]
    restriction = constraints["forcing"]["restriction"]
    assert "No residual-dependent basis or pointwise free force" in restriction
    assert "reject collapsed candidates" in constraints["nontriviality"]["enforcement"]


def audit() -> dict[str, Any]:
    scope = _load_json(SCOPE_PATH)
    base = scope["base"]
    source_path = ROOT / base["audited_source_path"]
    constraints_path = ROOT / base["canonical_constraints_path"]

    source_blob = _git_blob_sha(source_path)
    constraints_blob = _git_blob_sha(constraints_path)
    assert source_blob == base["audited_source_blob"], (
        "audited PR #903 source drifted; refresh this audit instead of silently "
        "reusing its conclusions"
    )
    assert constraints_blob == base["canonical_constraints_blob"], (
        "canonical CR001 constraints drifted; refresh invariants before reuse"
    )

    _assert_canonical_constraints(scope, _load_json(constraints_path))

    tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    builder = _find_function(tree, "build_canonical_quadrature_evidence")
    builder_args = _all_arg_names(builder)
    constructor_keywords = _builder_keyword_values(builder)

    # This is the precise seam governed by this audit.  The builder receives
    # metric rows and identities, but no external execution/producer receipt,
    # while two provenance-bearing booleans are asserted internally.
    prohibited_missing_today = {
        "execution_receipt_sha256",
        "producer_identity_sha256",
        "quadrature_output_artifact_sha256",
        "independence_evidence_sha256",
    }
    no_external_execution_provenance_input = builder_args.isdisjoint(prohibited_missing_today)
    builder_self_asserts_convergence = _is_literal_true(
        constructor_keywords.get("convergence_assessed")
    )
    builder_self_asserts_independence = _is_literal_true(
        constructor_keywords.get("independent_of_training_and_held_in_selection")
    )
    builder_blocks_residual_defined_forcing = (
        isinstance(constructor_keywords.get("residual_defined_forcing"), ast.Constant)
        and constructor_keywords["residual_defined_forcing"].value is False
    )

    assert no_external_execution_provenance_input
    assert builder_self_asserts_convergence
    assert builder_self_asserts_independence
    assert builder_blocks_residual_defined_forcing

    truth_states = scope["independent_truth_states"]
    assert truth_states["canonical_eq45_velocity_export_ready"] is True
    assert truth_states["kokuno_real_canonical_quadrature_execution_verified"] is False
    assert truth_states["kokuno_pde_validated"] is False
    assert truth_states["visual_correspondence_verified"] is False
    assert truth_states["paper_exact"] is False
    assert truth_states["openai_field_identified"] is False

    return {
        "audit_id": scope["audit_id"],
        "audited_head": base["audited_head"],
        "audited_source_blob": source_blob,
        "canonical_constraints_blob": constraints_blob,
        "builder_parameters": sorted(builder_args),
        "external_execution_provenance_input_present": not no_external_execution_provenance_input,
        "builder_self_asserts_convergence": builder_self_asserts_convergence,
        "builder_self_asserts_independence": builder_self_asserts_independence,
        "builder_blocks_residual_defined_forcing": builder_blocks_residual_defined_forcing,
        "receipt_integrity_is_execution_provenance": False,
        "checksum_valid_self_built_receipt_proves_independent_run": False,
        "scientific_state_promoted": False,
        "pde_validated": False,
        "audit_passed": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = audit()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
