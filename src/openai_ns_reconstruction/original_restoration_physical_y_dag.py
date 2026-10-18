"""Focused held-phase physical-y checks for restoration and postrestore.

The accepted physical-y DAG implementation supplies the pure hash-consed
``_DAG``/``_Oracle`` rules and root-record serialization used here.  This
module only reads the pinned C6, C5, and Q4 archives as data and checks the
two new charts ``restoration`` and ``postrestore`` (C6 calls the latter
``post_restore``).  The old post-power/long/reference 36-root receipt is
consumed as a saved lemma; its roots are deliberately not replayed here.

The returned receipt is conditional.  In particular, the native field seed
identities, smooth lazy-flat branch, and positive-loop/IFT hypotheses remain
external premises.  This packet does not claim full Q_y, a seam, or global
five-moment closure.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

try:  # Package import is used by the repository runners.
    from .original_physical_y_dag_derivative import (
        BudgetExceeded,
        GraphProofError,
        OpenRuleError,
        _DAG,
        _Oracle,
        _canonical_json,
        _first_dag_mismatch,
        _mapping_digest,
        _operation_counts,
        _positive_j_check,
        _radius_relation,
        _read_gzip_json,
        _reachable,
        _root_check,
        _sha256,
    )
except ImportError:  # pragma: no cover - direct source execution fallback
    from openai_ns_reconstruction.original_physical_y_dag_derivative import (
        BudgetExceeded,
        GraphProofError,
        OpenRuleError,
        _DAG,
        _Oracle,
        _canonical_json,
        _first_dag_mismatch,
        _mapping_digest,
        _operation_counts,
        _positive_j_check,
        _radius_relation,
        _read_gzip_json,
        _reachable,
        _root_check,
        _sha256,
    )


Q4_STEM = "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions"
C5_STEM = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions"
C6_STEM = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions"

PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
OLD_PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
OLD_PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
GENERIC_MODULE_RELATIVE = "src/openai_ns_reconstruction/original_physical_y_dag_derivative.py"
GENERIC_MODULE_SHA256 = "1fd508af8b3261006649de58a433bb05b20ae0d077852e847d8026f9f153b49a"

SCOPE_CHARTS = ("restoration", "postrestore")
# The C6 archive uses the source-function spelling ``post_restore``.
C6_CHART_FOR_Q4 = {"restoration": "restoration", "postrestore": "post_restore"}
FIELD_NAMES = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
ORDERS = tuple(range(6))
AXIAL_SCOPE = {
    "Z_lower": "37/100",
    "Z_upper": "93/250",
    "Z_interval": ["37/100", "93/250"],
    "ordinary_Z_orders": [0, 1],
    "source": "accepted inherited source domain; no extension to all axial Z",
    "inherited": True,
    "global_five_moment_closure": False,
}

# These are the exact source/cache/check bytes consumed by the new worker.
# The old generic worker's source guard and parent metadata are read from the
# accepted parent receipt below; no old worker build is invoked.
DEPENDENCY_HASHES = {
    f"experiments/root_st073/{C5_STEM}.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    f"experiments/root_st073/{C5_STEM}.json.gz": "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    f"experiments/root_st073/{C5_STEM}_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    f"experiments/root_st073/{C6_STEM}.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    f"experiments/root_st073/{C6_STEM}.json.gz": "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    f"experiments/root_st073/{C6_STEM}_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    f"experiments/root_st073/{Q4_STEM}.py": "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    f"experiments/root_st073/{Q4_STEM}.json.gz": "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    f"experiments/root_st073/{Q4_STEM}_check.json": "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
    f"experiments/root_st073/{OLD_PARENT_NAME}": OLD_PARENT_SHA256,
    GENERIC_MODULE_RELATIVE: GENERIC_MODULE_SHA256,
}
# Keep the accepted worker's public constant spellings available to callers
# that compose this focused packet with the earlier DAG receipt.
EXPECTED_PARENT_SHA256 = PARENT_SHA256
EXPECTED_INPUT_SHA256 = DEPENDENCY_HASHES


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_parent(root: str | Path, supplied: Any) -> tuple[Path, Mapping[str, Any], dict[str, Any]]:
    """Load the exact accepted modified-seams parent and verify its mapping."""

    path = Path(root) / "experiments" / "root_st073" / PARENT_NAME
    actual_hash = _sha256(path)
    disk = _read_json(path)
    if not isinstance(disk, Mapping):
        raise GraphProofError("accepted parent report is not a mapping")
    supplied_path: str | None = None
    selected = supplied
    if isinstance(supplied, (str, Path)):
        supplied_path = str(supplied)
        selected = _read_json(Path(supplied))
    mapping_exact = isinstance(selected, Mapping) and dict(selected) == dict(disk)
    receipt = {
        "path": str(path),
        "expected_sha256": PARENT_SHA256,
        "actual_sha256": actual_hash,
        "mapping_sha256": _mapping_digest(disk),
        "supplied_path": supplied_path,
        "supplied_mapping_sha256": _mapping_digest(selected) if isinstance(selected, Mapping) else None,
        "mapping_exact": mapping_exact,
        "passed": actual_hash == PARENT_SHA256 and mapping_exact,
        "source_executed": False,
    }
    return path, disk, receipt


def _pin_one(root: str | Path, relative: str, expected: str, *, kind: str) -> dict[str, Any]:
    try:
        path = Path(root) / Path(relative)
        actual = _sha256(path)
        passed = actual == expected
        reason = None if passed else "sha256 mismatch"
    except Exception as exc:  # pragma: no cover - malformed checkout receipt
        path = None
        actual = None
        passed = False
        reason = str(exc)
    return {
        "path": relative,
        "kind": kind,
        "expected_sha256": expected,
        "actual_sha256": actual,
        "passed": passed,
        "source_executed": False,
        "reason": reason,
    }


def _source_pins(root: str | Path, module_path: Path) -> tuple[dict[str, Any], bool]:
    pins: dict[str, Any] = {}
    passed = True
    for relative, expected in DEPENDENCY_HASHES.items():
        kind = "implementation" if relative.endswith(".py") else "cache" if relative.endswith((".json", ".json.gz")) else "input"
        pin = _pin_one(root, relative, expected, kind=kind)
        pins[relative] = pin
        passed = passed and bool(pin["passed"])
    parent_pin = _pin_one(root, f"experiments/root_st073/{PARENT_NAME}", PARENT_SHA256, kind="parent")
    pins[f"experiments/root_st073/{PARENT_NAME}"] = parent_pin
    passed = passed and bool(parent_pin["passed"])
    self_relative = "src/openai_ns_reconstruction/original_restoration_physical_y_dag.py"
    self_hash = _sha256(module_path)
    pins[self_relative] = {
        "path": self_relative,
        "kind": "implementation",
        "expected_sha256": self_hash,
        "actual_sha256": self_hash,
        "passed": True,
        "source_executed": False,
        "self_binding": True,
    }
    return pins, passed


def _saved_generic_lemma(parent: Mapping[str, Any], source_family: Mapping[str, Any]) -> dict[str, Any]:
    """Consume the old 36-root receipt without executing its build/checker."""

    saved = parent.get("actual_physical_y_dag_derivative")
    if not isinstance(saved, Mapping):
        return {
            "passed": False,
            "consumed_as_saved_lemma": True,
            "replayed": False,
            "reason": "parent has no saved actual_physical_y_dag_derivative receipt",
        }
    guard = saved.get("source_guard")
    snippets = saved.get("source_guard_snippets")
    gates = saved.get("gates")
    old_parent_validation = saved.get("parent_validation")
    expected_gate_keys = (
        "all_36_actual_roots",
        "all_root_structural_checks",
        "all_radius_relations",
        "all_inverse_J_bindings_and_independence",
        "all_J_independence",
        "all_held_phase_cache_zero",
        "no_unsupported_reachable_operations",
        "no_active_open_rules",
        "all_required_rule_checks",
        "lazy_metadata_and_active_body_lemmas",
        "binder_scope_and_capture",
    )
    guard_ok = bool(
        isinstance(guard, Mapping)
        and guard.get("all_rule_bindings_present") is True
        and isinstance(snippets, list)
        and all(isinstance(row, (list, tuple)) and len(row) == 2 for row in snippets)
    )
    gates_ok = isinstance(gates, Mapping) and all(gates.get(key) is True for key in expected_gate_keys)
    old_parent_ok = isinstance(old_parent_validation, Mapping) and old_parent_validation.get("passed") is True
    source_family_ok = saved.get("source_family") == dict(source_family)
    passed = bool(
        saved.get("passed") is True
        and saved.get("conditional_passed") is True
        and saved.get("source_executed") is False
        and saved.get("accepted_ancestor_execution") is False
        and saved.get("physical_derivative_identity_proved") is False
        and saved.get("seed_identity_discharged") is False
        and saved.get("module_sha256") == GENERIC_MODULE_SHA256
        and saved.get("parent_report") == OLD_PARENT_NAME
        and saved.get("parent_report_sha256") == OLD_PARENT_SHA256
        and guard_ok
        and gates_ok
        and old_parent_ok
        and source_family_ok
    )
    return {
        "passed": passed,
        "consumed_as_saved_lemma": True,
        "replayed": False,
        "old_root_count": 36,
        "old_scope": saved.get("scope"),
        "old_module_sha256": saved.get("module_sha256"),
        "old_parent_report": saved.get("parent_report"),
        "old_parent_report_sha256": saved.get("parent_report_sha256"),
        "old_parent_validation": copy.deepcopy(dict(old_parent_validation)) if isinstance(old_parent_validation, Mapping) else None,
        "source_guard": copy.deepcopy(dict(guard)) if isinstance(guard, Mapping) else None,
        "source_guard_snippets": [tuple(row) for row in snippets] if isinstance(snippets, list) else [],
        "saved_gates": copy.deepcopy(dict(gates)) if isinstance(gates, Mapping) else {},
        "source_family_match": source_family_ok,
        "reason": None if passed else "saved generic source/parent lemma is incomplete or mismatched",
    }


def _seed_bindings_new(
    oracle: _Oracle,
    dag: _DAG,
    chart: str,
    c6_chart_name: str,
    held_phase: int,
    c6_bindings: Mapping[str, Any],
    q4_bindings: Mapping[str, Any],
) -> tuple[dict[int, int], dict[str, Any], list[dict[str, Any]]]:
    """Build native C6 field seeds only; finite-N modified rows are excluded."""

    c6_chart = c6_bindings["charts"][c6_chart_name]
    q4_chart = q4_bindings["charts"][chart]
    J = int(c6_chart["physical_measure"])
    # C6 owns the native coordinate/phase handle.  Q4 stores only the held
    # actual phase on its radial port, so use the C6 phase for the native
    # variable identity (it equals the C5 coordinate on these charts).
    coordinate = int(c6_chart["phase"])
    native_alias = oracle.name(coordinate)
    radial = q4_chart["radial_background"]
    native_fields = c6_chart["fields"]
    seed_map: dict[int, int] = {}
    field_rows: dict[str, list[dict[str, Any]]] = {}
    table: list[dict[str, Any]] = []
    for field in FIELD_NAMES:
        source_rows = native_fields.get(field)
        radial_rows = radial.get(field)
        if not isinstance(source_rows, list) or len(source_rows) != 7:
            raise GraphProofError(f"native C6 field row length mismatch {c6_chart_name}.{field}")
        if not isinstance(radial_rows, list) or len(radial_rows) != 7:
            raise GraphProofError(f"Q4 radial background row length mismatch {chart}.{field}")
        receipts: list[dict[str, Any]] = []
        for order, (source, radial_row) in enumerate(zip(source_rows, radial_rows)):
            source = int(source)
            radial_row = int(radial_row)
            expression = dag.mul(oracle.canonical(J), oracle.canonical(radial_row))
            old = seed_map.get(source)
            if old is not None and old != expression:
                raise GraphProofError(f"conflicting native seed for source node {source}")
            seed_map[source] = expression
            row = {
                "chart": chart,
                "c6_native_chart": c6_chart_name,
                "native_alias": native_alias,
                "J_handle": J,
                "field_name": field,
                "ordinary_order": order,
                "source_field_node": source,
                "background_field_handle": source,
                "radial_background_handle": radial_row,
                "seed_source": "C6 native fields",
                "finite_N_modified_fields_used": False,
                "seed_Dq_formula": "J*background_radial_row",
                "seed_dag_digest": dag.digest(expression),
                "seed_dag_node_count": len(dag.nodes),
            }
            table.append(row)
            receipts.append({
                "ordinary_order": order,
                "source_field_node": source,
                "physical_radial_row_node": radial_row,
                "seed_expression_sha256": dag.digest(expression),
                "seed_source": "C6 native fields",
            })
        field_rows[field] = receipts
    receipt = {
        "variable_phase_handle": int(c6_chart["phase"]),
        "native_alias": native_alias,
        "held_phase_handle": int(held_phase),
        "c6_chart_name": c6_chart_name,
        "physical_J_node": J,
        "native_field_row_count": 7,
        "native_fields": list(FIELD_NAMES),
        "finite_N_modified_fields_used": False,
        "field_rows": field_rows,
        "seed_identity_obligation": "D_native(field)=J*background_radial_row for every C6 native field row",
        "seed_identity_discharged": False,
    }
    return seed_map, receipt, table


def _chart_rows(c5_chart: Mapping[str, Any], q4_chart: Mapping[str, Any], c6_chart: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "native_fields": {
            "fields": list(FIELD_NAMES),
            "ordinary_rows_per_field": 7,
            "row_count": len(FIELD_NAMES) * 7,
            "source": "C6 native fields",
        },
        "C5_loop": {
            "A_rows": len(c5_chart["loop"]["A"]),
            "B_over_Pstar_rows": len(c5_chart["loop"]["B_over_Pstar"]),
            "ordinary_orders": list(ORDERS),
        },
        "Q4_slow_physical_y": {
            "slow_A_y_rows": len(q4_chart["radial"]["slow_A_y"]),
            "slow_B_y_rows": len(q4_chart["radial"]["slow_B_y"]),
            "ordinary_orders": list(ORDERS),
        },
        "J": int(c6_chart["physical_measure"]),
        "new_root_rows": 2 * len(ORDERS),
        "finite_N_modified_fields_used": False,
    }


def _rule_coverage(
    nodes: list[dict[str, Any]],
    source_reachable: set[int],
    actual_reachable: set[int],
    oracle: _Oracle,
) -> dict[str, Any]:
    source_ops = _operation_counts(nodes, source_reachable)
    actual_ops = _operation_counts(nodes, actual_reachable)
    supported = {
        "bound_variable", "sum", "product", "negative", "positive_quotient",
        "analytic_unary", "exact_original_source_flat_sigma", "function_substitution",
        "definite_integral", "substitute_original_inverse_angle", "original_monotone_phase_inverse",
        "original_lazy_flat_branch", "current_original_unique_analytic_core_component",
        "current_original_leading_function_recipe", "current_original_axis_pressure_function_derivative",
        "current_original_source_parameter", "exact_real_comparison", "logical_or",
        "exact_rational", "exact_positive_integer_expression", "mathematical_pi",
        "mathematical_positive_infinity",
    }
    unsupported = sorted(operation for operation in set(source_ops) | set(actual_ops) if operation not in supported)
    unsupported_clear = not unsupported
    open_rules_clear = not oracle.open_rules
    rule_checks = {
        "sum": source_ops.get("sum", 0) > 0 and actual_ops.get("sum", 0) > 0,
        "product": source_ops.get("product", 0) > 0 and actual_ops.get("product", 0) > 0,
        "quotient": source_ops.get("positive_quotient", 0) > 0 and actual_ops.get("positive_quotient", 0) > 0,
        "unary": source_ops.get("analytic_unary", 0) > 0 and actual_ops.get("analytic_unary", 0) > 0,
        "squared_flat_sigma": source_ops.get("original_lazy_flat_branch", 0) > 0 and bool(oracle.lazy_lemma_receipts),
        "function_substitution_chain": source_ops.get("function_substitution", 0) > 0 and actual_ops.get("function_substitution", 0) > 0,
        "definite_integral_partial_and_endpoints": source_ops.get("definite_integral", 0) > 0 and actual_ops.get("definite_integral", 0) > 0,
        "inverse_angle_chain": source_ops.get("substitute_original_inverse_angle", 0) > 0 and actual_ops.get("substitute_original_inverse_angle", 0) > 0,
        "implicit_inverse": source_ops.get("original_monotone_phase_inverse", 0) > 0 and actual_ops.get("original_monotone_phase_inverse", 0) > 0,
        "lazy_flat_branch": source_ops.get("original_lazy_flat_branch", 0) > 0 and actual_ops.get("original_lazy_flat_branch", 0) > 0,
    }
    return {
        "source_operation_counts": source_ops,
        "actual_derivative_operation_counts": actual_ops,
        "required_rules": list(rule_checks),
        "rule_checks": rule_checks,
        "rule_checks_passed": all(rule_checks.values()),
        "unsupported_reachable_operations": unsupported,
        "unsupported_reachable_operations_clear": unsupported_clear,
        "active_open_rules": sorted(oracle.open_rules),
        "active_open_rules_clear": open_rules_clear,
        "lazy_metadata_issues": list(oracle.lazy_metadata_issues),
        "binder_scope_issues": list(oracle.binder_scope_issues),
        "lazy_lemma_receipts": list(oracle.lazy_lemma_receipts),
        "nested_tuple_expansion": False,
        "reachable_graph_only": True,
        "smooth_control_remains_external": True,
    }


def build(root: str | Path, parent: Any = None) -> dict[str, Any]:
    """Build the conditional restoration/postrestore 24-root receipt."""

    started = time.monotonic()
    root_path = Path(root)
    module_path = root_path / "src" / "openai_ns_reconstruction" / "original_restoration_physical_y_dag.py"
    parent_path, parent_mapping, parent_validation = _load_parent(root_path, parent)

    experiment_root = root_path / "experiments" / "root_st073"
    archive_paths = {
        "q4": experiment_root / f"{Q4_STEM}.json.gz",
        "c5": experiment_root / f"{C5_STEM}.json.gz",
        "c6": experiment_root / f"{C6_STEM}.json.gz",
    }
    pins, pins_ok = _source_pins(root_path, module_path)

    # Archives are intentionally read as JSON data.  No producer, owner,
    # checker, recurrence, or quadrature module is imported or called.
    q4 = _read_gzip_json(archive_paths["q4"])
    c5 = _read_gzip_json(archive_paths["c5"])
    c6 = _read_gzip_json(archive_paths["c6"])
    nodes = q4.get("exact_graph_nodes")
    c5_nodes = c5.get("exact_graph_nodes")
    c6_nodes = c6.get("exact_graph_nodes")
    if not isinstance(nodes, list) or not isinstance(c5_nodes, list) or not isinstance(c6_nodes, list):
        raise GraphProofError("missing exact graph node arrays")
    q4_prefix = nodes[: len(c5_nodes)] == c5_nodes
    c5_prefix = c5_nodes[: len(c6_nodes)] == c6_nodes

    q4_bindings = q4.get("modified_Q4_pressure_C3_bindings")
    c5_bindings = c5.get("prescribed_loop_C5_bindings")
    c6_bindings = c6.get("prescribed_R100_Rm_C6_bindings")
    if not all(isinstance(value, Mapping) for value in (q4_bindings, c5_bindings, c6_bindings)):
        raise GraphProofError("missing saved Q4/C5/C6 source bindings")
    q4_charts = q4_bindings.get("charts")
    c5_charts = c5_bindings.get("charts")
    c6_charts = c6_bindings.get("charts")
    if not all(isinstance(value, Mapping) for value in (q4_charts, c5_charts, c6_charts)):
        raise GraphProofError("missing saved chart binding maps")

    source_family = q4.get("source_family")
    family_match = bool(
        isinstance(source_family, Mapping)
        and source_family == c5.get("source_family")
        and source_family == c6.get("source_family")
        and parent_mapping.get("source_family") == dict(source_family)
    )
    generic_lemma = _saved_generic_lemma(parent_mapping, source_family if isinstance(source_family, Mapping) else {})

    dag = _DAG()
    oracle = _Oracle(nodes, dag)
    chart_receipts: dict[str, Any] = {}
    seed_table: list[dict[str, Any]] = []
    seed_tables: dict[str, Any] = {}
    all_roots: list[dict[str, Any]] = []
    all_radius = True
    all_j = True
    all_j_independence = True
    all_phase_zero = True
    all_unsupported_clear = True
    all_open_rules_clear = True
    all_rule_checks = True
    all_roots_passed = True
    all_chart_rows = True

    for chart in SCOPE_CHARTS:
        c5_chart = c5_charts[chart]
        c6_chart_name = C6_CHART_FOR_Q4[chart]
        c6_chart = c6_charts[c6_chart_name]
        q4_chart = q4_charts[chart]
        if int(c5_chart["physical_measure"]) != int(c6_chart["physical_measure"]):
            raise GraphProofError(f"C5/C6 physical J mismatch {chart}")
        if int(c6_chart["physical_measure"]) != 1:
            raise GraphProofError(f"new chart {chart} does not have J=1")

        chart_seed_map, seed_receipt, rows = _seed_bindings_new(
            oracle, dag, chart, c6_chart_name, int(c5_chart["actual_phase"]), c6_bindings, q4_bindings
        )
        seed_table.extend(rows)
        seed_tables[chart] = rows
        coordinate = int(c5_chart["coordinate"])
        native_name = oracle.name(coordinate)
        phase = int(c5_chart["actual_phase"])
        source_rows = list(c5_chart["loop"]["A"]) + list(c5_chart["loop"]["B_over_Pstar"])
        actual_rows = list(q4_chart["radial"]["slow_A_y"]) + list(q4_chart["radial"]["slow_B_y"])
        source_reachable = _reachable(nodes, source_rows)
        actual_reachable = _reachable(nodes, actual_rows)
        row_checks: list[dict[str, Any]] = []
        for field, source_key, actual_key in (
            ("A", "A", "slow_A_y"),
            ("B_over_Pstar", "B_over_Pstar", "slow_B_y"),
        ):
            source_jet = list(c5_chart["loop"][source_key])
            actual_jet = list(q4_chart["radial"][actual_key])
            if len(source_jet) != 6 or len(actual_jet) != 6:
                raise GraphProofError(f"{chart}.{field} must expose six roots")
            for order, (source, actual) in enumerate(zip(source_jet, actual_jet)):
                receipt = _root_check(
                    oracle,
                    dag,
                    int(source),
                    int(actual),
                    native_name,
                    coordinate,
                    phase,
                    chart_seed_map,
                )
                receipt.update({
                    "chart": chart,
                    "field": field,
                    "ordinary_order": order,
                    "source_kind": "C5 loop A/B",
                    "actual_kind": "Q4 slow A/B physical-y",
                })
                row_checks.append(receipt)
                all_roots.append(receipt)

        radius = _radius_relation(oracle, c5_chart, dag)
        positive_j = _positive_j_check(oracle, chart, c5_chart, c6_chart, q4_chart)
        phase_zero = oracle.derivative(
            phase,
            native_name,
            ("node", coordinate),
            phase,
            chart_seed_map,
            ("node", coordinate),
            oracle.seed_profile(chart_seed_map),
        ) == dag.zero
        coverage = _rule_coverage(nodes, source_reachable, actual_reachable, oracle)
        rows_receipt = _chart_rows(c5_chart, q4_chart, c6_chart)
        rows_ok = bool(
            rows_receipt["native_fields"]["ordinary_rows_per_field"] == 7
            and rows_receipt["C5_loop"]["A_rows"] == 6
            and rows_receipt["C5_loop"]["B_over_Pstar_rows"] == 6
            and rows_receipt["Q4_slow_physical_y"]["slow_A_y_rows"] == 6
            and rows_receipt["Q4_slow_physical_y"]["slow_B_y_rows"] == 6
            and rows_receipt["J"] == 1
        )
        all_chart_rows = all_chart_rows and rows_ok
        all_radius = all_radius and bool(radius["R_q_equals_JR"])
        all_j = all_j and bool(positive_j["passed"])
        all_j_independence = all_j_independence and bool(positive_j["J_independence_gate_passed"])
        all_phase_zero = all_phase_zero and phase_zero
        all_roots_passed = all_roots_passed and all(bool(row["passed"]) for row in row_checks)
        all_unsupported_clear = all_unsupported_clear and bool(coverage["unsupported_reachable_operations_clear"])
        all_open_rules_clear = all_open_rules_clear and bool(coverage["active_open_rules_clear"])
        all_rule_checks = all_rule_checks and bool(coverage["rule_checks_passed"])
        chart_receipts[chart] = {
            "native_coordinate_handle": coordinate,
            "native_alias": native_name,
            "c6_native_chart": c6_chart_name,
            "held_actual_phase_handle": phase,
            "physical_J_handle": int(c6_chart["physical_measure"]),
            "J": 1,
            "seed_profile": oracle.seed_profile(chart_seed_map),
            "source_loop_roots": {
                "A": [int(value) for value in c5_chart["loop"]["A"]],
                "B_over_Pstar": [int(value) for value in c5_chart["loop"]["B_over_Pstar"]],
            },
            "actual_slow_roots": {
                "A": [int(value) for value in q4_chart["radial"]["slow_A_y"]],
                "B_over_Pstar": [int(value) for value in q4_chart["radial"]["slow_B_y"]],
                "slow_A_y": [int(value) for value in q4_chart["radial"]["slow_A_y"]],
                "slow_B_y": [int(value) for value in q4_chart["radial"]["slow_B_y"]],
            },
            "source_reachable_node_count": len(source_reachable),
            "actual_reachable_node_count": len(actual_reachable),
            "seed_bindings": seed_receipt,
            "root_receipts": row_checks,
            "radius_relation": radius,
            "exact_q_Rq_JR": radius,
            "positive_J": positive_j,
            "held_phase_cache_zero": phase_zero,
            "J_independence_gate_passed": bool(positive_j["J_independence_gate_passed"]),
            "new_chart_rows": rows_receipt,
            "rule_coverage": coverage,
        }

    # Native seed equations, smoothness, and IFT/positive-loop controls are
    # intentionally external.  Structural 24-root success is therefore a
    # conditional result, not a full physical-field or seam theorem.
    conditional = bool(
        len(all_roots) == 24
        and all_roots_passed
        and all_radius
        and all_j
        and all_j_independence
        and all_phase_zero
        and all_chart_rows
        and all_unsupported_clear
        and all_open_rules_clear
        and all_rule_checks
        and not oracle.lazy_metadata_issues
        and not oracle.binder_scope_issues
        and q4_prefix
        and c5_prefix
        and pins_ok
        and parent_validation["passed"]
        and family_match
        and generic_lemma["passed"]
    )
    module_hash = _sha256(module_path)
    result: dict[str, Any] = {
        "passed": conditional,
        "conditional_passed": conditional,
        "all_passed": False,
        "parent_report_sha256": PARENT_SHA256,
        "parent_report": PARENT_NAME,
        "module_sha256": module_hash,
        "generic_dependency_sha256": GENERIC_MODULE_SHA256,
        "physical_derivative_identity_proved": False,
        "seed_identity_discharged": False,
        "source_executed": False,
        "pure_generic_rule_engine_new24_only": True,
        "generic_rule_engine_audit": {
            "applied_root_count": 24,
            "accepted_old_36_root_arithmetic_replayed": False,
            "accepted_old_builders_or_checkers_replayed": False,
            "old_generic_receipt_consumed_as_saved_lemma": True,
        },
        "accepted_ancestor_execution": False,
        "accepted_builders_executed": False,
        "accepted_checkers_executed": False,
        "accepted_producers_executed": False,
        "owners_executed": False,
        "norms_executed": False,
        "recurrences_executed": False,
        "quadrature_executed": False,
        "source_family": copy.deepcopy(source_family),
        "source_family_match": family_match,
        "scope": "actual Q4/C5 held-fast-phase slow A/B physical-y roots for restoration and postrestore",
        "axial_scope": copy.deepcopy(AXIAL_SCOPE),
        "parent_validation": parent_validation,
        "saved_generic_lemma": generic_lemma,
        "source_pins": pins,
        "source_pins_passed": pins_ok,
        "source_guard": copy.deepcopy(generic_lemma.get("source_guard")),
        "source_guard_snippets": [tuple(row) for row in generic_lemma.get("source_guard_snippets", [])],
        "sourceguards_tuple": all(
            isinstance(row, tuple) and len(row) == 2
            for row in generic_lemma.get("source_guard_snippets", [])
        ),
        "gates": {
            "all_24_actual_roots": len(all_roots) == 24,
            "all_root_structural_checks": all_roots_passed,
            "all_radius_relations": all_radius,
            "all_inverse_J_bindings_and_independence": all_j,
            "all_J_independence": all_j_independence,
            "all_held_phase_cache_zero": all_phase_zero,
            "all_new_chart_rows": all_chart_rows,
            "no_nested_tuple_expansion": True,
            "no_unsupported_reachable_operations": all_unsupported_clear,
            "no_active_open_rules": all_open_rules_clear,
            "all_required_rule_checks": all_rule_checks,
            "lazy_metadata_and_active_body_lemmas": not oracle.lazy_metadata_issues,
            "binder_scope_and_capture": not oracle.binder_scope_issues,
            "parent_bytes_hash_mapping": bool(parent_validation["passed"]),
            "source_family_strict": family_match,
            "saved_old_generic_lemma": bool(generic_lemma["passed"]),
            "seed_and_smooth_external": True,
            "positive_loop_IFT_external": True,
        },
        "lazy_lemma_receipts": list(oracle.lazy_lemma_receipts),
        "lazy_metadata_issues": list(oracle.lazy_metadata_issues),
        "binder_scope_issues": list(oracle.binder_scope_issues),
        "binder_scope_receipts": list(oracle.binder_scope_receipts),
        "graph_prefix": {
            "Q4_over_C5": q4_prefix,
            "C5_over_C6": c5_prefix,
            "Q4_nodes": len(nodes),
            "C5_nodes": len(c5_nodes),
            "C6_nodes": len(c6_nodes),
        },
        "charts": chart_receipts,
        "new_chart_rows": {chart: chart_receipts[chart]["new_chart_rows"] for chart in SCOPE_CHARTS},
        "seed_tables": seed_tables,
        "seed_binding_table": seed_table,
        "counts": {
            "Q4_graph_nodes": len(nodes),
            "C5_graph_nodes": len(c5_nodes),
            "C6_graph_nodes": len(c6_nodes),
            "scoped_charts": len(SCOPE_CHARTS),
            "new_chart_count": len(SCOPE_CHARTS),
            "native_field_count_per_chart": len(FIELD_NAMES),
            "native_field_rows_per_chart": len(FIELD_NAMES) * 7,
            "C5_A_rows_per_chart": 6,
            "C5_B_rows_per_chart": 6,
            "Q4_slow_A_y_rows_per_chart": 6,
            "Q4_slow_B_y_rows_per_chart": 6,
            "slow_derivative_roots_total": len(all_roots),
            "slow_derivative_roots_passed": sum(1 for row in all_roots if row["passed"]),
            "seed_rows_total": len(seed_table),
            "hash_cons_dag_nodes": len(dag.nodes),
            "hash_cons_dag_unique_keys": len(dag.keys),
            "hash_cons_dag_source_leaf_nodes": sum(1 for node in dag.nodes if node[0] == "source_leaf"),
            "derivative_cache_entries": len(oracle._derivative_cache),
            "elapsed_seconds": time.monotonic() - started,
        },
        "limits": {
            "full_Qy": False,
            "full_modified_Qy": False,
            "seam_or_global_closure": False,
            "global_five_moment_closure": False,
            "no_numeric_values": True,
            "accepted_execution": False,
            "actual_native_seed_identity_external": True,
            "smooth_lazy_branch_control_external": True,
            "positive_loop_controls_and_IFT_external": True,
            "finite_N_modified_fields_used": False,
            "unscoped_Q4_charts": sorted(set(q4_charts) - set(SCOPE_CHARTS)),
        },
        "external_premises": [
            "actual C6 native field identities D_native(field)=J*background_radial_row",
            "smooth actual sigma/lazy-flat branch and binder semantics",
            "positive loop controls and implicit-function-theorem smoothness",
            "the accepted parent source/transport lemmas on Z in [37/100,93/250]",
        ],
        "open_limits": [
            "native physical seed equations remain an upstream obligation",
            "smooth lazy-flat branch and positive-loop IFT controls remain external",
            "full Q_y and any restoration/postrestore seam are outside this DAG packet",
            "global five-moment closure, heat, stress, recursion, and corrected Cartesian NS are not claimed",
            "no old post_power/long_reshape/reference roots were replayed",
        ],
    }
    # A stable receipt hash is useful to downstream composition, but excludes
    # itself from the hashed payload.
    result["receipt_sha256"] = hashlib.sha256(
        _canonical_json(result).encode("utf-8")
    ).hexdigest()
    return result


__all__ = [
    "BudgetExceeded",
    "GraphProofError",
    "OpenRuleError",
    "build",
]
