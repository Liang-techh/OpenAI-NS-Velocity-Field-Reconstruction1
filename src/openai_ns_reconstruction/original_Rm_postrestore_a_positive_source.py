"""Bounded analytic source proof for the post-restore loop root ``a``.

The proof binds the actual post-restore C6 constructors to the saved generic
C5 rows.  It only inspects caller-supplied source text/receipts and a small
set of graph rows; it never runs a producer, an ancestor operator, a
quadrature, or a numerical field evaluator.  The saved whole-cover bound is
used to authenticate the stated Z domain, while positivity is proved by the
constructor identity ``E=Utheta/Pstar``.
"""
from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_postrestore_a_positive_source_v1"
PREFIX = "lei_ren_part1_paper_compliant_current_"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300

PARENT_REPORT = PREFIX + "original_Rm_physical_Qy_Z_calculus.json"
PARENT_REPORT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
LATEST_DELTA = PREFIX + "original_Rm_leading_m_ordinary_Z_sources.json"
LATEST_DELTA_SHA256 = "02cb39ce0272143655566f2a876bfdc1a083ea90c7aaf869506d401ed10ffcf8"

Z_NODE = 20
PSTAR_NODE = 33572
COORDINATE_NODE = 84528
ACTUAL_E = 84911
ACTUAL_EY = 85146
ACTUAL_A = 119435
OPTIONAL_A_Z1 = 119439
OPTIONAL_A_Z2 = 119445

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

# These are fixed source/receipt pins from the authoritative parent input
# ledger.  Gzip archives are intentionally not inputs to this checker.
RECORD_HASHES = {
    PREFIX + "original_prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "original_prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "original_prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "original_prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "original_prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "original_prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.py": "a590cbf346ddbf3cc6298e16b4c28fe7ed330d0acdb14e60d280b83e9ddf7adf",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json": "2817a7849a3bd943acb2e0ffaebae53d0f4d2c396b0c15f0f9b4d6effa8364bd",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds_check.json": "7350b72291cb6278afaafd407685eb006fd17b503abb7c5f5f852aab878b804c",
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
    expected = RECORD_HASHES[name]
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
    for declaration in wanted:
        if declaration not in found:
            gates.append({"gate": "source_declaration_missing", "record": name, "declaration": declaration})
    return found


def _tokens(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in wanted if token not in text]
    if missing:
        gates.append({"gate": "source_constructor_token_missing", "record": name, "tokens": missing})
        return False
    return True


def _path(value: Any, path: tuple[str, ...]) -> Any:
    for key in path:
        value = value.get(key) if isinstance(value, Mapping) else None
    return value


def _receipt(text: str, name: str, paths: tuple[tuple[str, ...], ...], gates: list[dict[str, Any]]) -> tuple[bool, Any]:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        gates.append({"gate": "accepted_source_receipt_invalid", "record": name, "error": str(exc)})
        return False, None
    ok = value.get("all_passed") is True and value.get("source_family") == SOURCE_FAMILY
    if not ok:
        gates.append({"gate": "accepted_source_receipt_header_missing", "record": name})
    for path in paths:
        if _path(value, path) is not True:
            gates.append({"gate": "accepted_source_receipt_flag_missing", "record": name, "path": ".".join(path)})
            ok = False
    return ok, value


def _generic_domain(text: str, name: str, gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        gates.append({"gate": "generic_bounds_receipt_invalid", "record": name, "error": str(exc)})
        return False, {}
    provenance = _path(value, ("current_original_source_log_bound_charts", "restore_buffer", "original_source_provenance"))
    z_box = provenance.get("Z_box") if isinstance(provenance, Mapping) else None
    theorem = _path(value, ("positive_noncore_source_denominator_theorem", "source_charts", "restore_buffer"))
    source_flag = theorem.get("source_function_positivity_not_inferred_from_saved_box") if isinstance(theorem, Mapping) else None
    ok = (
        value.get("source_family") == SOURCE_FAMILY
        and isinstance(z_box, Mapping)
        and z_box.get("lower") == "-1.0"
        and z_box.get("upper") == "1.0"
        and source_flag is True
    )
    if not ok:
        gates.append({"gate": "whole_Z_source_domain_or_independence_changed", "record": name})
    return ok, {
        "chart": "restore_buffer",
        "Z_box": {"lower": "-1", "upper": "1"},
        "source_function_positivity_not_inferred_from_saved_box": source_flag is True,
        "saved_box_used_for_E_positivity": False,
    }


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    texts = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(text is None for text in texts.values()):
        return False, {"records": dict(RECORD_HASHES)}

    bridge = PREFIX + "original_prescribed_bridge_C6_functions.py"
    bridge_check = PREFIX + "original_prescribed_bridge_C6_functions_check.json"
    r100 = PREFIX + "original_prescribed_R100_Rm_C6_functions.py"
    r100_check = PREFIX + "original_prescribed_R100_Rm_C6_functions_check.json"
    loop = PREFIX + "original_prescribed_loop_C5_functions.py"
    loop_check = PREFIX + "original_prescribed_loop_C5_functions_check.json"
    bounds_py = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.py"
    bounds = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json"
    bounds_check = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds_check.json"

    declarations = {
        bridge: _definitions(texts[bridge], bridge, ("BridgeAlgebra", "generic_fields"), gates),
        r100: _definitions(texts[r100], r100, ("CurrentOriginalPrescribedR100RmC6Functions", "phi_from_Utheta", "add_chart"), gates),
        loop: _definitions(texts[loop], loop, ("generic_roots", "CurrentOriginalPrescribedLoopC5Functions"), gates),
        bounds_py: _definitions(texts[bounds_py], bounds_py, ("positive_denominators", "CurrentGenericSourceBounds"), gates),
    }
    phi_methods = [item for item in ast.walk(ast.parse(texts[r100]))
                   if isinstance(item, ast.FunctionDef) and item.name == "phi_from_Utheta"]
    expected_phi_return = ast.parse(
        "a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R))),"
        "'same positive sqrt(2R)*F0; F0 is Z-only')", mode="eval").body
    phi_return_ok = (len(phi_methods) == 1 and isinstance(phi_methods[0].body[-1], ast.Return)
                     and ast.dump(phi_methods[0].body[-1].value, include_attributes=False)
                     == ast.dump(expected_phi_return, include_attributes=False))
    if not phi_return_ok:
        gates.append({"gate": "actual_phi_from_Utheta_return_constructor_changed", "record": r100})
    token_ok = all((
        _tokens(texts[bridge], bridge, (
            "def generic_fields", "inverse=g.quotient(g.one,Pstar", "W=a.scale(F0",
            "fields=dict(E=a.mul(W,actual['phi'])", "g.symbol('Z')",
        ), gates),
        _tokens(texts[r100], r100, (
            "p=g.symbol('true_prescribed_post_restore_fraction')", "Utheta=a.exp(logUtheta)",
            "a.scale(phi,'-2/5')", "self.add_chart('post_restore',p,R,g.one,actual",
            "self.charts[name]=dict(phase=phase,domain=[0,1]", "same positive sqrt(2R)*F0; F0 is Z-only",
        ), gates),
        _tokens(texts[loop], loop, (
            "def generic_roots", "C=a.sub(E,a.scale(a.trim(radial['E_y']),2))",
            "divide=lambda jet:a.div(jet,E,'same strictly positive actual prescribed swirl E')",
            "'postrestore':background.charts['post_restore']",
        ), gates),
    ))
    bridge_ok, _ = _receipt(texts[bridge_check], bridge_check, (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho"),
    ), gates)
    r100_ok, r100_receipt = _receipt(texts[r100_check], r100_check, (
        ("checks", "actual_source", "passed"),
        ("checks", "geometry", "passed"),
        ("checks", "actual_source", "endpoint_phi_U_and_physical_y_derivatives_exposed"),
        ("checks", "geometry", "complete_radius_order_proved_by_actual_source_positive_guards"),
    ), gates)
    if isinstance(r100_receipt, Mapping):
        charts = _path(r100_receipt, ("checks", "actual_source", "actual_prescribed_charts_checked"))
        if not isinstance(charts, list) or "post_restore" not in charts:
            gates.append({"gate": "post_restore_chart_missing_from_C6_receipt", "record": r100_check})
            r100_ok = False
    loop_ok, _ = _receipt(texts[loop_check], loop_check, (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "genuine_C6_inputs_and_C5_generic_roots"),
    ), gates)
    bounds_check_ok, _ = _receipt(texts[bounds_check], bounds_check, (
        ("current_original_noncore_generic_quotient_derivative_log_bounds_certified",),
    ), gates)
    domain_ok, domain_witness = _generic_domain(texts[bounds], bounds, gates)
    ok = token_ok and phi_return_ok and bridge_ok and r100_ok and loop_ok and bounds_check_ok and domain_ok
    return ok, {
        "records": dict(RECORD_HASHES),
        "declarations": declarations,
        "accepted_receipts": {
            bridge_check: {"sha256": RECORD_HASHES[bridge_check], "passed": bridge_ok},
            r100_check: {"sha256": RECORD_HASHES[r100_check], "passed": r100_ok},
            loop_check: {"sha256": RECORD_HASHES[loop_check], "passed": loop_ok},
            bounds_check: {"sha256": RECORD_HASHES[bounds_check], "passed": bounds_check_ok},
        },
        "constructor_tokens_checked": token_ok,
        "phi_from_Utheta_return_ast_bound": phi_return_ok,
        "phi_return_ast_digest": canonical_digest(ast.dump(expected_phi_return, include_attributes=False)),
        "generic_bounds_domain": domain_witness,
        "source_family": deepcopy(SOURCE_FAMILY),
    }


def _graph_witnesses(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    try:
        rows: dict[str, Any] = {}
        rows["Z"] = _expect_row(nodes, Z_NODE, "bound_variable", name="Z")
        rows["Pstar"] = _expect_row(nodes, PSTAR_NODE, "analytic_unary", name="exp", argument=5)
        rows["Pstar_inverse"] = _expect_row(
            nodes, 69245, "positive_quotient", numerator=1, denominator=PSTAR_NODE,
            source_positive_certificate="same strictly positive Pstar",
        )
        rows["coordinate"] = _expect_row(nodes, COORDINATE_NODE, "bound_variable", name="true_prescribed_post_restore_fraction")
        rows["radius_110"] = _expect_row(nodes, 79968, "exact_rational", numerator=110, denominator=1)
        rows["radius_minus_seven"] = _expect_row(nodes, 84529, "exact_rational", numerator=-7, denominator=1)
        rows["radius_offset"] = _expect_row(nodes, 84530, "sum", arguments=[84529, COORDINATE_NODE])
        rows["radius_log"] = _expect_row(nodes, 84531, "sum", arguments=[81514, 84530])
        rows["radius_exp"] = _expect_row(nodes, 84532, "analytic_unary", name="exp", argument=84531)
        rows["radius"] = _expect_row(nodes, 84533, "product", arguments=[79968, 84532])
        rows["sqrt_radius_argument"] = _expect_row(nodes, 84564, "product", arguments=[48, 84533])
        rows["sqrt_radius"] = _expect_row(nodes, 84565, "analytic_unary", name="sqrt", argument=84564)
        rows["W_inverse_Pstar"] = _expect_row(nodes, 84895, "product", arguments=[84565, 69245])
        rows["W"] = _expect_row(nodes, 84896, "product", arguments=[35657, 84895])
        rows["F0_real_exponential"] = _expect_row(nodes, 35657, "analytic_unary", name="exp", argument=35656)
        rows["F0_exponent"] = _expect_row(nodes, 35656, "sum", arguments=[35647, 35652])
        rows["Utheta_real_exponential"] = _expect_row(nodes, 84537, "analytic_unary", name="exp", argument=84536)
        rows["Utheta_exponent"] = _expect_row(nodes, 84536, "sum", arguments=[79981, 84535])
        rows["Phi_denominator"] = _expect_row(nodes, 84566, "product", arguments=[35657, 84565])
        rows["Phi"] = _expect_row(
            nodes, 84573, "positive_quotient", numerator=84537, denominator=84566,
            source_positive_certificate="same positive sqrt(2R)*F0; F0 is Z-only",
        )
        rows["E"] = _expect_row(nodes, ACTUAL_E, "product", arguments=[84896, 84573])
        rows["half_phi"] = _expect_row(nodes, 85132, "product", arguments=[84573, 28])
        rows["half"] = _expect_row(nodes, 28, "exact_rational", numerator=1, denominator=2)
        rows["phi_y"] = _expect_row(nodes, 84618, "product", arguments=[84573, 76958])
        rows["minus_two_fifths"] = _expect_row(nodes, 76958, "exact_rational", numerator=-2, denominator=5)
        rows["Ey_inner"] = _expect_row(nodes, 85139, "sum", arguments=[85132, 84618])
        rows["Ey"] = _expect_row(nodes, ACTUAL_EY, "product", arguments=[84896, 85139])
        rows["two_Ey"] = _expect_row(nodes, 119417, "product", arguments=[85146, 48])
        rows["two"] = _expect_row(nodes, 48, "exact_rational", numerator=2, denominator=1)
        rows["minus_two_Ey"] = _expect_row(nodes, 119423, "negative", argument=119417)
        rows["C"] = _expect_row(nodes, 119429, "sum", arguments=[84911, 119423])
        rows["a"] = _expect_row(
            nodes, ACTUAL_A, "positive_quotient", numerator=119429, denominator=ACTUAL_E,
            source_positive_certificate="same strictly positive actual prescribed swirl E",
        )
        optional = {}
        for handle in (OPTIONAL_A_Z1, OPTIONAL_A_Z2):
            row = _row(nodes, handle)
            _require(row.get("operation") == "positive_quotient", f"optional a root changed at graph handle {handle}")
            _require(row.get("denominator") == ACTUAL_E, f"optional a denominator changed at graph handle {handle}")
            _require(row.get("source_positive_certificate") == "same strictly positive actual prescribed swirl E",
                     f"optional a positivity guard changed at graph handle {handle}")
            optional[str(handle)] = {"handle": handle, "numerator": row.get("numerator"), "denominator": ACTUAL_E}
        half = Fraction(1, 2)
        phi_y_ratio = Fraction(-2, 5)
        ey_ratio = half + phi_y_ratio
        c_ratio = 1 - 2 * ey_ratio
        _require(ey_ratio == Fraction(1, 10) and c_ratio == Fraction(4, 5), "exact source rational identity changed")
        return True, {
            "rows": rows,
            "optional_same_denominator_roots": optional,
            "rational_atoms": {"one_half": "1/2", "phi_y_over_phi": "-2/5", "E_y_over_E": "1/10", "C_over_E": "4/5"},
            "positive_source_factors": {
                "F0": "exp(35656)>0 on the accepted real C6 source domain",
                "Utheta": "exp(84536)>0 on the accepted real C6 source domain",
                "Pstar": "exp(5)>0, where 5 is an argument handle",
                "sqrt_2R": "sqrt(2*110*exp(84531))>0",
                "Phi_denominator": "same actual F0*sqrt(2R), both explicitly positive",
                "E": "exact factor cancellation gives Utheta/Pstar>0",
                "positivity_guard_strings_used_as_premise": False,
            },
            "source_rows_digest": canonical_digest([[name, rows[name]] for name in sorted(rows)]),
        }
    except (TypeError, ValueError, KeyError) as exc:
        gates.append({"gate": "saved_actual_postrestore_source_rows_changed", "error": str(exc)})
        return False, {}


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed analytic source receipt for the post-restore ``a`` root."""
    _require(expected_nodes_digest == EXPECTED_NODES_DIGEST, "unexpected source graph digest assertion")
    _require(len(nodes) == EXPECTED_NODE_COUNT and canonical_digest(nodes) == EXPECTED_NODES_DIGEST,
             "saved source graph changed")
    gates: list[dict[str, Any]] = []
    source_ok, source = _source_witnesses(source_records, gates)
    graph_ok, graph = _graph_witnesses(nodes, gates)
    admitted = source_ok and graph_ok
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "passed": admitted,
        "actual_target_handle": ACTUAL_A,
        "whole_native_domain": ["0", "1"],
        "native_context": {
            "coordinate_handle": COORDINATE_NODE,
            "coordinate_name": "true_prescribed_post_restore_fraction",
            "interval": [0, 1],
            "phase_held_fixed": True,
            "source_ref": PREFIX + "original_prescribed_R100_Rm_C6_functions.py",
        },
        "whole_Z_domain": ["-1", "1"],
        "Z_context": {
            "variable_handle": Z_NODE,
            "variable_name": "Z",
            "interval": [-1, 1],
            "source_ref": "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json",
        },
        "strict_positive_proved": admitted,
        "lower_bound_exact": "4/5",
        "identity_source_witness": {
            "formula": "W=F0*sqrt(2R)/Pstar; Phi=Utheta/(sqrt(2R)*F0); E=W*Phi=Utheta/Pstar",
            "radial_formula": "E_y=W*(Phi/2+Phi_y), Phi_y=-2*Phi/5, C=E-2*E_y, a=C/E",
            "exact_reduction": "E_y=E/10; C=4*E/5; a=4/5",
            "positivity": "Utheta=exp(log_Utheta)>0, Pstar=exp(logP)>0, hence E>0 on the whole source domain",
            "handles": {"W": 84896, "Phi": 84573, "E": ACTUAL_E, "E_y": ACTUAL_EY, "C": 119429, "a": ACTUAL_A},
            "source_constructor_checks": source,
            "saved_graph_checks": graph,
            "saved_box_used_for_positivity": False,
        },
        "authoritative_parent": {"path": "experiments/root_st073/" + PARENT_REPORT, "sha256": PARENT_REPORT_SHA256},
        "latest_delta": {"path": "experiments/root_st073/" + LATEST_DELTA, "sha256": LATEST_DELTA_SHA256},
        "unresolved_gates": gates,
        "old_operators_replayed": False,
        "numeric_values_computed": False,
        "full_goal_complete": False,
        "source_record_hashes": dict(RECORD_HASHES),
    }
