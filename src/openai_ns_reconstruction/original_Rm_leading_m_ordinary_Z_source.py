"""Admit the source provenance of the leading physical ``m`` Z jet.
This module is a source bridge.  It reads a frozen graph and caller-supplied
copies of the already accepted producer sources as data, then checks their
hashes, AST declarations, receipt flags, and the six saved graph rows.  It
does not execute a producer, evaluate a field, or construct a new derivative
graph.  A function-jet array or a C5/C6 label without these checks is never a
derivative premise.
"""
from __future__ import annotations
import ast
from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
from typing import Any
SCHEMA = "original_Rm_leading_m_ordinary_Z_source_v1"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300
Z_NODE = 20
PSTAR_NODE = 33572
# These are the parent-calculus input pins, not caller assertions.  The
# runner must provide exactly these UTF-8 source/receipt records.
RECORD_HASHES = {
    PREFIX + "prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "true_Rm_patch_C6_functions.py": "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    PREFIX + "true_Rm_patch_C6_functions_check.json": "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7",
    PREFIX + "core_inlet_axial6_generic_source.py": "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
    PREFIX + "core_inlet_axial6_generic_source_check.json": "9d8b746ef78502956cbb2b53f96d452be636eecbf645d021aaca3f9c2e21885d",
}
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
CONTEXTS = {
    "postrestore": {"coordinate": 84528, "endpoint": 1},
    "patch": {"coordinate": 185676, "endpoint": 1},
}
EDGES = {
    "postrestore": (96436, 96437, 96438),
    "patch": (195000, 195001, 195002),
}
def canonical_digest(value: Any) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()
def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)
def _row(nodes: Sequence[Any], index: int) -> Mapping[str, Any]:
    _require(type(index) is int and 0 <= index < len(nodes), f"graph handle out of range: {index}")
    row = nodes[index]
    _require(isinstance(row, Mapping), f"graph handle is not a mapping: {index}")
    return row
def _expect_row(nodes: Sequence[Any], index: int, operation: str, **fields: Any) -> dict[str, Any]:
    row = _row(nodes, index)
    _require(row.get("operation") == operation, f"operation changed at graph handle {index}")
    for key, value in fields.items():
        _require(row.get(key) == value, f"source field {key} changed at graph handle {index}")
    return {"handle": index, "operation": operation, "fields": deepcopy(fields)}
def _record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> str | None:
    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        gates.append({"gate": "source_record_missing", "record": name})
        return None
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        gates.append({"gate": "source_record_path_changed", "record": name, "expected": expected_path})
        return None
    text = record.get("source_text")
    actual = hashlib.sha256(text.encode("utf-8")).hexdigest() if isinstance(text, str) else None
    if record.get("sha256") != expected or actual != expected:
        gates.append({"gate": "source_record_hash_changed", "record": name, "expected_sha256": expected})
        return None
    return text
def _definitions(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as exc:
        gates.append({"gate": "source_ast_invalid", "record": name, "error": str(exc)})
        return {}
    found: dict[str, dict[str, int]] = {}
    for item in ast.walk(tree):
        if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
            found[item.name] = {"line": item.lineno, "end_line": item.end_lineno or item.lineno}
    for item in wanted:
        if item not in found:
            gates.append({"gate": "source_declaration_missing", "record": name, "declaration": item})
    return found
def _text_tokens(text: str, name: str, tokens: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in tokens if token not in text]
    if missing:
        gates.append({"gate": "source_constructor_token_missing", "record": name, "tokens": missing})
        return False
    return True
def _receipt(
    text: str,
    name: str,
    required: tuple[tuple[tuple[str, ...], bool], ...],
    gates: list[dict[str, Any]],
) -> bool:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        gates.append({"gate": "accepted_source_receipt_invalid", "record": name, "error": str(exc)})
        return False
    ok = value.get("all_passed") is True and value.get("source_family") == SOURCE_FAMILY
    if not ok:
        gates.append({"gate": "accepted_source_receipt_header_missing", "record": name})
    for path, expected in required:
        current: Any = value
        for key in path:
            current = current.get(key) if isinstance(current, Mapping) else None
        if current is not expected:
            gates.append({
                "gate": "accepted_source_receipt_flag_missing",
                "record": name,
                "path": ".".join(path),
                "expected": expected,
            })
            ok = False
    return ok
def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    texts = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(text is None for text in texts.values()):
        return False, {"records": {name: RECORD_HASHES[name] for name in RECORD_HASHES}}
    bridge_name = PREFIX + "prescribed_bridge_C6_functions.py"
    r100_name = PREFIX + "prescribed_R100_Rm_C6_functions.py"
    r100_check = PREFIX + "prescribed_R100_Rm_C6_functions_check.json"
    patch_name = PREFIX + "true_Rm_patch_C6_functions.py"
    patch_check = PREFIX + "true_Rm_patch_C6_functions_check.json"
    core_name = PREFIX + "core_inlet_axial6_generic_source.py"
    bridge_check = PREFIX + "prescribed_bridge_C6_functions_check.json"
    core_check = PREFIX + "core_inlet_axial6_generic_source_check.json"
    expected_declarations = {
        bridge_name: {"BridgeAlgebra", "histories", "generic_fields"},
        r100_name: {"CurrentOriginalPrescribedR100RmC6Functions", "add_chart", "end"},
        patch_name: {"implicit_rows", "patch_fields", "CurrentOriginalTrueRmPatchC6Functions"},
        core_name: {"FunctionJet", "JetAlgebra"},
    }
    declarations = {
        name: _definitions(texts[name], name, tuple(sorted(expected)))
        for name, expected in expected_declarations.items()
    }
    declaration_ok = True
    for name, expected in expected_declarations.items():
        actual = set(declarations[name])
        if actual != expected:
            gates.append({"gate": "source_declaration_set_changed", "record": name,
                          "expected": sorted(expected), "actual": sorted(actual)})
            declaration_ok = False
    token_ok = all((
        _text_tokens(texts[bridge_name], bridge_name, (
            "original_incoming_memory_retained=True", "a.integral(", "m=a.scale(J['M'],inverse)",
        ), gates),
        _text_tokens(texts[r100_name], r100_name, ("Z_independent_substitution=True", "previous.histories", "previous.generic_fields"), gates),
        _text_tokens(texts[patch_name], patch_name, (
            "same_full_Jacobian_used=True", "mass=a.add(a.scale(z,4)", "fields=dict(E=a.mul(am,H),V=",
            "m=a.scale(mass,inv)", "Z_independent_substitution=True",
        ), gates),
        _text_tokens(texts[core_name], core_name, (
            "ordinary_slow_Z_derivative_order=k", "endpoint_independent_of_Z=True", "def mul(", "def integral(",
        ), gates),
    ))
    r100_ok = _receipt(texts[r100_check], r100_check, (
        (("checks", "actual_source", "passed"), True),
        (("checks", "actual_source", "all_genuine_seven_source_and_six_generic_rows"), True),
        (("checks", "actual_source", "every_nonzero_incoming_history_function_retained_by_identity"), True),
        (("checks", "actual_source", "endpoint_phi_U_and_physical_y_derivatives_exposed"), True),
        (("checks", "actual_source", "finite_N_modified_source_and_quantitative_C6_not_claimed"), True),
        (("checks", "geometry", "complete_radius_order_proved_by_actual_source_positive_guards"), True),
    ), gates)
    patch_ok = _receipt(texts[patch_check], patch_check, (
        (("checks", "actual_source", "passed"), True),
        (("checks", "actual_source", "ordinary_implicit_rows_not_padded"), True),
        (("checks", "actual_source", "unique_branch_defining_residual_bound_by_theorem_not_numeric_zero"), True),
        (("checks", "actual_source", "original_log_Jacobian_and_finite_patch_endpoints"), True),
        (("checks", "actual_source", "no_finite_N_history_reset_or_numeric_containment_claim"), True),
        (("checks", "branch_theorem", "same_exact_original_unique_branch_theorem_bound"), True),
        (("checks", "branch_theorem", "original_common_box_strict_Neumann_bound_applies_to_same_exact_Jacobian"), True),
        (("checks", "branch_theorem", "sixth_function_rows_follow_implicit_function_theorem_not_sixth_numeric_norm"), True),
        (("checks", "implicit_calculus", "exact_full_Jacobian_inverse_and_ordinary_recurrence_checked"), True),
        (("checks", "normalized_conjugacy", "normalized_Neumann_theorem_proves_raw_Jacobian_invertible"), True),
    ), gates)
    bridge_receipt_ok = _receipt(texts[bridge_check], bridge_check, (
        (("checks", "actual_source", "passed"), True),
        (("checks", "actual_source", "preserved_core_C6_before_actual_C7_comparison_source"), True),
        (("checks", "actual_source", "F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho"), True),
        (("checks", "own_six_histories", "passed"), True),
        (("checks", "own_six_histories", "all_nonzero_incoming_histories_and_sixth_derivatives_retained"), True),
        (("checks", "own_six_histories", "nonunit_physical_jacobian_applied_once"), True),
    ), gates)
    core_receipt_ok = _receipt(texts[core_check], core_check, (
        (("checks", "source_rows", "passed"), True),
        (("checks", "source_rows", "original_P0_C3_preserved_and_higher_rows_are_explicit_function_derivatives"), True),
        (("checks", "physical_units", "generic_V_mean_k_invPstar_applied_once"), True),
        (("checks", "independent_calculus", "nonunit_Pstar_without_extra_half_and_single_radius_factor_checked"), True),
    ), gates)
    ok = token_ok and declaration_ok and r100_ok and patch_ok and bridge_receipt_ok and core_receipt_ok
    return ok, {
        "records": {name: RECORD_HASHES[name] for name in RECORD_HASHES},
        "declarations": declarations,
        "accepted_receipts": {
            r100_check: {"sha256": RECORD_HASHES[r100_check], "passed": r100_ok},
            patch_check: {"sha256": RECORD_HASHES[patch_check], "passed": patch_ok},
            bridge_check: {"sha256": RECORD_HASHES[bridge_check], "passed": bridge_receipt_ok},
            core_check: {"sha256": RECORD_HASHES[core_check], "passed": core_receipt_ok},
        },
        "constructor_semantics": {
            "bridge": "BridgeAlgebra.histories -> JetAlgebra.integral; generic_fields scales M by inverse Pstar",
            "core": "JetAlgebra ordinary rows use Leibniz/quotient recurrences and Z-independent endpoint substitution",
            "patch": "implicit_rows uses the same full Jacobian ordinary recurrence; patch_fields scales physical mass by inverse Pstar",
        },
    }
def _pstar_constant_closure(nodes: Sequence[Any]) -> dict[str, Any]:
    seen: set[int] = set()
    def visit(index: int) -> None:
        if index in seen:
            return
        row = _row(nodes, index)
        operation = row.get("operation")
        if operation == "exact_rational":
            _require(type(row.get("numerator")) is int and type(row.get("denominator")) is int
                     and row["denominator"] != 0, "malformed constant rational")
            children = []
        elif operation == "exact_positive_integer_expression":
            _require(row.get("Z_independent") is True and row.get("exact_positive_integer") is True,
                     "source integer lacks explicit constant premise")
            children = []
        elif operation in ("sum", "product"):
            children = row.get("arguments")
        elif operation in ("negative", "analytic_unary"):
            children = [row.get("argument")]
        elif operation == "positive_quotient":
            _require(bool(row.get("source_positive_certificate")), "constant quotient lacks guard")
            children = [row.get("numerator"), row.get("denominator")]
        else:
            raise ValueError("Pstar constant closure rejects operation " + str(operation))
        _require(isinstance(children, list), "constant source children are malformed")
        seen.add(index)
        for child in children:
            _require(type(child) is int and 0 <= child < index, "constant source child is not backward")
            visit(child)
    visit(PSTAR_NODE)
    return {"root": PSTAR_NODE, "argument_handle": 5, "Z_independent": True,
            "used_nodes": sorted(seen), "source_rows_digest": canonical_digest([[i, nodes[i]] for i in sorted(seen)]),
            "inverse_rows": [69245, 183811]}


def _graph_witnesses(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    try:
        _expect_row(nodes, Z_NODE, "bound_variable", name="Z")
        _expect_row(nodes, PSTAR_NODE, "analytic_unary", name="exp", argument=5)
        pstar_closure = _pstar_constant_closure(nodes)
        _expect_row(nodes, 69245, "positive_quotient", numerator=1, denominator=PSTAR_NODE,
                    source_positive_certificate="same strictly positive Pstar")
        _expect_row(nodes, 183811, "positive_quotient", numerator=1, denominator=PSTAR_NODE,
                    source_positive_certificate="same positive Pstar")
        _expect_row(nodes, 13, "exact_rational", numerator=4, denominator=1)
        _expect_row(nodes, 35566, "product", arguments=[Z_NODE, 13])
        _expect_row(nodes, 186855, "positive_quotient", numerator=1, denominator=185676,
                    source_positive_certificate="positive original x=R/Rm")
        _expect_row(nodes, 185676, "bound_variable", name="true_Rm_patch_x")
        post = []
        for i, (primitive, moment, integral, substitution, order) in enumerate((
            (96366, 96373, 84756, 96226, 0),
            (96367, 96374, 84757, 96227, 1),
            (96368, 96375, 84758, 96228, 2),
        )):
            _expect_row(nodes, primitive, "product", arguments=[substitution, 84740])
            _expect_row(nodes, substitution, "function_substitution", expression=95309 + i,
                        variable=83010, value=1, Z_independent_substitution=True)
            _expect_row(nodes, integral, "definite_integral", integrand=84746 + i, variable="true_post_restore_history_t",
                        lower=0, upper=84528, actual_history_name="M", actual_six_history_rate=1,
                        exact_function_integral=True, numerical_cap_not_defining_integrand=True,
                        numerical_value_not_installed=True, endpoint_independent_of_Z=True,
                        ordinary_slow_Z_derivative_order=order, original_incoming_memory_retained=True,
                        prescribed_bridge_chart="true_post_restore", source_measure_applied_once=True,
                        true_core_inlet_complete_history=True)
            _expect_row(nodes, moment, "sum", arguments=[primitive, integral])
            post.append({"primitive": primitive, "history": integral, "moment": moment, "ordinary_order": order})
        for target, source in zip(EDGES["postrestore"], (96373, 96374, 96375)):
            _expect_row(nodes, target, "product", arguments=[source, 69245])
        _expect_row(nodes, 194630, "sum", arguments=[184044, 194616, 194623])
        _expect_row(nodes, 194631, "sum", arguments=[184045, 194617, 194624])
        _expect_row(nodes, 194632, "sum", arguments=[183815, 194618, 194625])
        _expect_row(nodes, 194840, "product", arguments=[194630, 186855])
        _expect_row(nodes, 194841, "product", arguments=[194631, 186855])
        _expect_row(nodes, 194842, "product", arguments=[194632, 186855])
        _expect_row(nodes, 194847, "sum", arguments=[35566, 194840])
        _expect_row(nodes, 194848, "sum", arguments=[13, 194841])
        _expect_row(nodes, 195000, "product", arguments=[194847, 183811])
        _expect_row(nodes, 195001, "product", arguments=[194848, 183811])
        _expect_row(nodes, 195002, "product", arguments=[194842, 183811])
        return True, {
            "pstar_constant_closure": pstar_closure,
            "postrestore": post,
            "patch": {"physical_mass": [194847, 194848, 194842], "scaled_m": list(EDGES["patch"])},
        }
    except ValueError as exc:
        gates.append({"gate": "saved_source_graph_constructor_changed", "error": str(exc)})
        return False, {}
def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed source receipt for the two leading m charts."""
    _require(expected_nodes_digest == EXPECTED_NODES_DIGEST, "unexpected source graph digest assertion")
    _require(len(nodes) == EXPECTED_NODE_COUNT and canonical_digest(nodes) == EXPECTED_NODES_DIGEST,
             "saved source graph changed")
    gates: list[dict[str, Any]] = []
    source_ok, source_witnesses = _source_witnesses(source_records, gates)
    graph_ok, graph_witnesses = _graph_witnesses(nodes, gates)
    admitted = source_ok and graph_ok
    rows = []
    for side in ("postrestore", "patch"):
        context = deepcopy(CONTEXTS[side])
        handles = EDGES[side]
        for order in (0, 1):
            rows.append({
                "side": side,
                "field": "m",
                "ordinary_order": order,
                "source_handle": handles[order],
                "derivative_handle": handles[order + 1],
                "source_derivative_proved": admitted,
                "regularity_proved": admitted,
                "context": context,
                "proof_origin": "accepted C6 source producer plus exact saved ordinary jet constructor",
                "source_receipt_ref": (
                    PREFIX + "prescribed_R100_Rm_C6_functions_check.json"
                    if side == "postrestore"
                    else PREFIX + "true_Rm_patch_C6_functions_check.json"
                ),
                "source_recipe_refs": (
                    [PREFIX + "prescribed_bridge_C6_functions.py", PREFIX + "prescribed_R100_Rm_C6_functions.py"]
                    if side == "postrestore"
                    else [PREFIX + "true_Rm_patch_C6_functions.py"]
                ),
            })
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "Z_variable": Z_NODE,
        "contexts": deepcopy(CONTEXTS),
        "rows": rows,
        "passed": admitted,
        "witnesses": {
            "producer_sources": source_witnesses,
            "saved_graph": graph_witnesses,
            "dependency_handles": {
                "Pstar": PSTAR_NODE,
                "postrestore_moment_rows": [96373, 96374, 96375],
                "patch_physical_mass_rows": [194847, 194848, 194842],
            },
        },
        "unresolved_gates": gates,
        "numeric_values_computed": False,
        "operators_reexecuted": False,
        "old_operators_replayed": False,
        "cross_side_value_equality_used": False,
        "full_goal_complete": False,
        "full_Qy1_or_full_goal_claimed": False,
        "source_record_hashes": dict(RECORD_HASHES),
    }
