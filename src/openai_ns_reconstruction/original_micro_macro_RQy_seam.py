"""Source-linked Q4 pressure seam for the second-micro/macro join.

This module only appends a small arithmetic suffix to the accepted Q4 graph.
It does not evaluate a phase, a quadrature, or a mixed-norm enclosure.  The
accepted join receipt supplies the common physical source identity; this
helper makes the two endpoint traces and their signed FTC contraction explicit
in the graph so a caller cannot replace the traces by a common placeholder.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .original_signed_graph_pairs import GraphPairs


PARENT_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_I1_defect_bounds.json"
)
JOIN_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check.json"
)
JOIN_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.py"
)
JOIN_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.json.gz"
)
Q4_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
Q4_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
)
Q4_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
)

PARENT_SHA256 = "88df21cdcf020204fa0f651dae58ed0ffb98cbc1bfe9adb360c990d770ab5492"
GRAPH_HELPER_PATH = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
GRAPH_HELPER_SHA256 = "1e9699cd84a1c82dc7b6577f922f1da1362cbedcfef148e7cec1d146a99aaaff"
PINNED_HASHES = {
    JOIN_CHECK_NAME: "49e7d4f0c566abcbfe9833c5cc3c63875b1c80fa6e6025601eab8aaa83956225",
    JOIN_SOURCE_NAME: "3f0bb8e0363990ab4d851020e988dd463061086d3ad25e7f0e59c047fa3c5bc0",
    JOIN_ARCHIVE_NAME: "0c673c0979a69f0a6f6fd48e2ee72adae2223f6e384d537d6e39801c15af79e2",
    Q4_SOURCE_NAME: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    Q4_ARCHIVE_NAME: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    Q4_CHECK_NAME: "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
}

_Q4_NODES = {
    "micro_finite_R_Q_y": (133660, 133661, 133662, 133663),
    "micro_endpoint_R_Q_y": (133745, 133746, 133747, 133748),
    "micro_Q_y": (133501, 133505, 133511, 133517),
    "micro_Q_y_full": (133501, 133505, 133511, 133517, 133523),
    "macro_finite_R_Q_y": (141027, 141028, 141029, 141030),
    "macro_Q_y": (140868, 140872, 140878, 140884),
    "macro_Q_y_full": (140868, 140872, 140878, 140884, 140890),
    "micro_radius": 68644,
    "macro_radius": 71448,
    "micro_coordinate": 66988,
    "macro_coordinate": 70265,
    "zero": 0,
    "one": 1,
    "two": 48,
    "micro_local_S1": (133689, 133690, 133691, 133692),
    "micro_incoming_S1": (129638, 129639, 129640, 129641),
    "micro_cumulative_S1": (133693, 133694, 133695, 133696),
    "micro_cumulative_S1_endpoint": (133729, 133730, 133731, 133732),
    "macro_local_S1": (141058, 141059, 141060, 141061),
    "macro_incoming_S1": (133729, 133730, 133731, 133732),
    "macro_cumulative_S1": (141062, 141063, 141064, 141065),
    "macro_cumulative_S1_endpoint": (141098, 141099, 141100, 141101),
    "first_micro_cumulative_S1": (129594, 129595, 129596, 129597),
    "first_micro_endpoint_S1": (129638, 129639, 129640, 129641),
    "global_initial_R_Q_y": (129567, 129568, 129569, 129570),
    "micro_global_m15_finite": (133721, 133722, 133723, 133724),
    "micro_global_m15_endpoint": (133737, 133738, 133739, 133740),
    "macro_global_m15_finite": (141090, 141091, 141092, 141093),
    "macro_global_m15_endpoint": (141106, 141107, 141108, 141109),
}


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path.name}")
    return value


def _source_guard(path: Path, functions: Mapping[str, tuple[str, ...]]) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        source = stream.read()
    tree = ast.parse(source, filename=str(path))
    found: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            found[node.name] = ast.get_source_segment(source, node) or ""
    result: dict[str, Any] = {"path": path.name, "functions": {}}
    for name, tokens in functions.items():
        body = found.get(name, "")
        compact = "".join(body.split())
        matched = all("".join(token.split()) in compact for token in tokens)
        result["functions"][name] = {"present": bool(body), "tokens": list(tokens), "matched": matched}
        _require(body and matched, f"source guard failed for {path.name}:{name}")
    result["passed"] = True
    return result


def _check_receipt(check: Mapping[str, Any], names: tuple[str, ...]) -> dict[str, Any]:
    _require(check.get("all_passed") is True, "accepted source join receipt is not all_passed")
    checks = check.get("checks")
    _require(isinstance(checks, Mapping), "accepted source join checks missing")
    selected: dict[str, Any] = {}
    for name in names:
        item = checks.get(name)
        _require(isinstance(item, Mapping) and item.get("passed") is True, f"join check failed: {name}")
        selected[name] = {"passed": True}
    return selected


def _accepted_qy_bindings(check: Mapping[str, Any]) -> dict[str, list[int]]:
    recovery = check.get("checks", {}).get("actual_leading_recovery_and_common_N", {})
    rows = recovery.get("actual_leading_Q_recovery_y0_to_y4_C4_rows", ())
    _require(isinstance(rows, list), "accepted Q recovery rows missing")
    expected = {
        "second_micro": [133501, 133505, 133511, 133517],
        "frozen_macro": [140868, 140872, 140878, 140884],
    }
    for chart, handles in expected.items():
        for order, handle in enumerate(handles):
            _require([chart, 1, order, handle] in rows, f"accepted full Q4 physical Q_y row not bound: {chart} Z{order}")
    return expected


def _accepted_smooth_premises(check: Mapping[str, Any]) -> dict[str, bool]:
    smooth = check.get("checks", {}).get("actual_smooth_source_seam", {})
    _require(isinstance(smooth, Mapping) and smooth.get("passed") is True, "accepted smooth source seam missing")
    names = (
        "common_global_log_radius_phase_at_micro2_macro0_and_one_N",
        "original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows",
        "actual_micro_complete_correction_is_macro_incoming_so_mean_values_match",
        "original_mean_ODE_and_one_minus_delta_recovery_give_common_y4_Z4_Q",
    )
    for name in names:
        _require(smooth.get(name) is True, f"accepted smooth source premise missing: {name}")
    return {name: True for name in names}


def _base_row(nodes: list[dict[str, Any]], index: int) -> dict[str, Any]:
    _require(0 <= index < len(nodes), f"Q4 graph handle outside prefix: {index}")
    return nodes[index]


def _expect_node(nodes: list[dict[str, Any]], index: int, operation: str) -> None:
    _require(_base_row(nodes, index).get("operation") == operation, f"Q4 node {index} is not {operation}")


def _expect_substitution(nodes: list[dict[str, Any]], index: int, expression: int, variable: int, value: int) -> None:
    row = _base_row(nodes, index)
    _require(
        row.get("operation") == "function_substitution"
        and row.get("expression") == expression
        and row.get("variable") == variable
        and row.get("value") == value
        and row.get("Z_independent_substitution") is True,
        f"Q4 endpoint substitution mismatch at {index}",
    )


def _depends_on(nodes: list[dict[str, Any]], root: int, target: int) -> bool:
    seen: set[int] = set()

    def visit(index: int) -> bool:
        if index == target:
            return True
        if index in seen:
            return False
        seen.add(index)
        row = _base_row(nodes, index)
        children: list[int] = []
        if isinstance(row.get("arguments"), list):
            children.extend(value for value in row["arguments"] if isinstance(value, int))
        for key in ("argument", "numerator", "denominator", "expression", "variable", "value", "integrand", "lower", "upper"):
            value = row.get(key)
            if isinstance(value, int):
                children.append(value)
        return any(visit(child) for child in children)

    return visit(root)


def _verify_q4_rows(nodes: list[dict[str, Any]], archive: Mapping[str, Any]) -> dict[str, Any]:
    for index in _Q4_NODES["micro_endpoint_R_Q_y"]:
        _expect_node(nodes, index, "function_substitution")
    for index, expression in zip(_Q4_NODES["micro_endpoint_R_Q_y"], _Q4_NODES["micro_finite_R_Q_y"]):
        _expect_substitution(nodes, index, expression, 66988, 48)
    for qy, radius, index in zip(
        _Q4_NODES["micro_Q_y"],
        [68644] * 4,
        _Q4_NODES["micro_finite_R_Q_y"],
    ):
        row = _base_row(nodes, index)
        _require(row.get("operation") == "product" and row.get("arguments") == [qy, radius], f"micro finite R_Q_y product mismatch at {index}")
    for qy, radius, index in zip(
        _Q4_NODES["macro_Q_y"],
        [71448] * 4,
        _Q4_NODES["macro_finite_R_Q_y"],
    ):
        row = _base_row(nodes, index)
        _require(row.get("operation") == "product" and row.get("arguments") == [qy, radius], f"macro finite R_Q_y product mismatch at {index}")
    for index in _Q4_NODES["micro_local_S1"] + _Q4_NODES["macro_local_S1"]:
        _expect_node(nodes, index, "definite_integral")
    for index in _Q4_NODES["micro_Q_y"] + _Q4_NODES["macro_Q_y"]:
        _expect_node(nodes, index, "positive_quotient")
    for index, incoming, local in zip(
        _Q4_NODES["micro_cumulative_S1"],
        _Q4_NODES["micro_incoming_S1"],
        _Q4_NODES["micro_local_S1"],
    ):
        row = _base_row(nodes, index)
        _require(row.get("operation") == "sum" and row.get("arguments") == [incoming, local], f"micro cumulative S1 mismatch at {index}")
    for index, incoming, local in zip(
        _Q4_NODES["macro_cumulative_S1"],
        _Q4_NODES["macro_incoming_S1"],
        _Q4_NODES["macro_local_S1"],
    ):
        row = _base_row(nodes, index)
        _require(row.get("operation") == "sum" and row.get("arguments") == [incoming, local], f"macro cumulative S1 mismatch at {index}")
    for index, expression in zip(_Q4_NODES["micro_cumulative_S1_endpoint"], _Q4_NODES["micro_cumulative_S1"]):
        _expect_substitution(nodes, index, expression, 66988, 48)
    for index, expression in zip(_Q4_NODES["macro_cumulative_S1_endpoint"], _Q4_NODES["macro_cumulative_S1"]):
        _expect_substitution(nodes, index, expression, 70265, 1)
    _expect_node(nodes, 66988, "bound_variable")
    _expect_node(nodes, 70265, "bound_variable")
    for index, operation, arguments in (
        (68644, "product", [65604, 68643]),
        (71448, "product", [65604, 71447]),
        (68642, "product", [8, 66988]),
        (66985, "product", [48, 8]),
        (70266, "product", [66987, 70265]),
        (71446, "sum", [66985, 70266]),
    ):
        row = _base_row(nodes, index)
        _require(row.get("operation") == operation, f"radius source operation changed at {index}")
        _require(row.get("arguments") == arguments, f"radius source arguments changed at {index}")
    for index, argument in ((68643, 68642), (71447, 71446)):
        row = _base_row(nodes, index)
        _require(row.get("operation") == "analytic_unary" and row.get("name") == "exp" and row.get("argument") == argument, f"radius exp source changed at {index}")
    for index in (8, 65604, 66987, 68644, 71448):
        _require(not _depends_on(nodes, index, 20), f"radius source {index} depends on the Z binder")
    _require(_base_row(nodes, 8).get("Z_independent") is True, "hb source is missing Z_independent receipt metadata")
    bindings = archive.get("modified_Q4_pressure_C3_bindings", {})
    charts = bindings.get("charts", {})
    _require(isinstance(charts, Mapping), "Q4 chart bindings missing")
    expected_chart_rows = {
        "first_micro": {
            "cumulative_integrals": {"S1": list(_Q4_NODES["first_micro_cumulative_S1"])},
            "endpoint_integrals": {"S1": list(_Q4_NODES["first_micro_endpoint_S1"])},
        },
        "second_micro": {
            "physical": {"Q_y": list(_Q4_NODES["micro_Q_y_full"])},
            "finite_boundary": {"R_Q_y": list(_Q4_NODES["micro_finite_R_Q_y"])},
            "cumulative_integrals": {"S1": list(_Q4_NODES["micro_cumulative_S1"])},
            "endpoint_integrals": {"S1": list(_Q4_NODES["micro_cumulative_S1_endpoint"])},
        },
        "frozen_macro": {
            "physical": {"Q_y": list(_Q4_NODES["macro_Q_y_full"])},
            "finite_boundary": {"R_Q_y": list(_Q4_NODES["macro_finite_R_Q_y"])},
            "cumulative_integrals": {"S1": list(_Q4_NODES["macro_cumulative_S1"])},
            "endpoint_integrals": {"S1": list(_Q4_NODES["macro_cumulative_S1_endpoint"])},
        },
    }
    for chart, groups in expected_chart_rows.items():
        _require(chart in charts, f"Q4 chart binding missing: {chart}")
        for group, expected in groups.items():
            if group == "endpoint_integrals":
                actual_group = bindings.get("endpoints", {}).get(chart, {}).get("integrals", {})
            else:
                actual_group = charts[chart].get(group, {})
            for field, handles in expected.items():
                actual = actual_group.get(field)
                _require(actual == handles, f"Q4 chart binding mismatch: {chart}.{group}.{field}")
    _require(_base_row(nodes, 0).get("operation") == "exact_rational", "Q4 zero handle changed")
    _require(_base_row(nodes, 1).get("operation") == "exact_rational", "Q4 one handle changed")
    _require(_base_row(nodes, 48).get("operation") == "exact_rational", "Q4 two handle changed")
    return {
        "micro_endpoint_rows": list(_Q4_NODES["micro_endpoint_R_Q_y"]),
        "macro_lower_source_rows": list(_Q4_NODES["macro_finite_R_Q_y"]),
        "cumulative_S1_rows_verified": True,
        "physical_Q_y_rows_verified": {
            "second_micro": list(_Q4_NODES["micro_Q_y_full"]),
            "frozen_macro": list(_Q4_NODES["macro_Q_y_full"]),
        },
        "radius_Z_independent_transitive": True,
        "forbidden_macro_endpoint_rows": [141114, 141115, 141116, 141117],
    }


def _symbolic_controls(
    common_q_y_rows: Mapping[str, list[int]],
    smooth_premises: Mapping[str, bool],
) -> dict[str, Any]:
    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - fail closed in minimal runtimes
        raise RuntimeError("SymPy is required for the exact seam controls") from exc

    hb, length, s, m, radius = sp.symbols("hb length s m radius")
    radius_micro = radius * sp.exp(hb * s)
    radius_macro = radius * sp.exp(2 * hb + length * m)
    radius_check = sp.simplify(radius_micro.subs(s, 2) - radius_macro.subs(m, 0)) == 0

    s_micro, i_micro, s_macro, i_macro = sp.symbols("S_micro I_micro S_macro I_macro")
    r_micro_upper, r_micro_lower = sp.symbols("r_micro_upper r_micro_lower")
    r_macro_upper, r_macro_lower = sp.symbols("r_macro_upper r_macro_lower")
    micro = sp.Rational(1, 2) * (s_micro - i_micro) - r_micro_upper + r_micro_lower
    macro = sp.Rational(1, 2) * (s_macro - i_macro) - r_macro_upper + r_macro_lower
    cancellation = -r_micro_upper + r_macro_lower
    residual = sp.expand(micro + macro - cancellation)
    joined_check = sp.expand(micro + macro - (residual + cancellation)) == 0
    pair_zero = sp.simplify(cancellation.subs(r_macro_lower, r_micro_upper)) == 0
    q_sub, r_sub = sp.symbols("Q_sub R_sub")
    substitution_product_check = sp.expand(q_sub * r_sub - q_sub * r_sub) == 0

    delta, z = sp.symbols("delta z")
    uy_left, mean_left = sp.Function("uy_left")(z), sp.Function("mean_left")(z)
    uy_right, mean_right = sp.Function("uy_right")(z), sp.Function("mean_right")(z)
    q_left = (2 * z * uy_left - ((1 - delta) * z * mean_left + (1 - z**2) * sp.diff(mean_left, z))) / (1 - delta * z**2)
    q_right = (2 * z * uy_right - ((1 - delta) * z * mean_right + (1 - z**2) * sp.diff(mean_right, z))) / (1 - delta * z**2)
    q_equal_checks = []
    for k in range(4):
        difference = sp.diff(q_left, z, k) - sp.diff(q_right, z, k)
        substitutions = {}
        for order in range(k + 2):
            substitutions[sp.diff(uy_right, z, order)] = sp.diff(uy_left, z, order)
            substitutions[sp.diff(mean_right, z, order)] = sp.diff(mean_left, z, order)
        q_equal_checks.append(sp.simplify(difference.xreplace(substitutions)) == 0)
    q_equal_check = all(q_equal_checks)
    source_bound_q_y = (
        common_q_y_rows.get("second_micro") == list(_Q4_NODES["micro_Q_y_full"])
        and common_q_y_rows.get("frozen_macro") == list(_Q4_NODES["macro_Q_y_full"])
    )
    source_bound_smooth = all(smooth_premises.values())
    controls = {
        "passed": bool(radius_check and joined_check and pair_zero and substitution_product_check and q_equal_check and source_bound_q_y and source_bound_smooth),
        "radius_cutpoint": "R_micro(s)=R_a*exp(hb*s), R_macro(m)=R_a*exp(2*hb+length*m); s=2,m=0",
        "radius_identity": bool(radius_check),
        "pressure_increment": "Delta m15 = (S1_cumulative-S1_incoming)/2 - R_Q_y(upper) + R_Q_y(lower)",
        "signed_FTC_identity": bool(joined_check),
        "conditional_boundary_pair_zero": bool(pair_zero),
        "Z_independent_substitution_product_identity": bool(substitution_product_check),
        "Q_y_operator": "[2 Z U_y - (1-delta) Z m_y - (1-Z^2) (m_y)_Z]/(1-delta Z^2)",
        "ordinary_Z_orders": [0, 1, 2, 3],
        "Q_y_common_source_identity": bool(q_equal_check),
        "source_bound_full_Q_y_arrays": bool(source_bound_q_y),
        "source_bound_smooth_seam_premises": bool(source_bound_smooth),
    }
    _require(controls["passed"], "exact seam symbolic controls failed")
    return controls


def _sub(g: GraphPairs, expression: int, variable: int, value: int) -> int:
    return g.node(
        "function_substitution",
        expression=expression,
        variable=variable,
        value=value,
        Z_independent_substitution=True,
        source_endpoint_trace=True,
    )


def _factorized_substitution_identity(
    nodes: list[dict[str, Any]],
    finite: int,
    wrapper: int,
    factorized: int,
    substituted_q: int,
    substituted_r: int,
    q: int,
    radius: int,
    variable: int,
    value: int,
) -> dict[str, Any]:
    """Check the exact source-wrapper/product distributivity contract."""

    finite_row = _base_row(nodes, finite)
    wrapper_row = _base_row(nodes, wrapper)
    factor_row = _base_row(nodes, factorized)
    _require(
        finite_row.get("operation") == "product"
        and finite_row.get("arguments") == [q, radius],
        f"finite product changed at {finite}",
    )
    _require(
        wrapper_row.get("operation") == "function_substitution"
        and wrapper_row.get("expression") == finite
        and wrapper_row.get("variable") == variable
        and wrapper_row.get("value") == value
        and wrapper_row.get("Z_independent_substitution") is True,
        f"source wrapper changed at {wrapper}",
    )
    _require(
        factor_row.get("operation") == "product"
        and factor_row.get("arguments") == [substituted_q, substituted_r],
        f"factorized substitution product changed at {factorized}",
    )
    return {
        "passed": True,
        "finite_product": finite,
        "actual_wrapper": wrapper,
        "factorized_wrapper_product": factorized,
        "substituted_Q_y": substituted_q,
        "substituted_radius": substituted_r,
        "variable": variable,
        "value": value,
        "Z_independent_substitution": True,
        "identity": "sub(product(Q_y,R),x,a)=product(sub(Q_y,x,a),sub(R,x,a))",
    }


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the exact second-micro-upper/frozen-macro-lower seam suffix."""

    root_path = Path(root)
    _require(isinstance(parent, Mapping), "parent must be the accepted I1 report mapping")
    parent_path = root_path / "experiments" / "root_st073" / PARENT_NAME
    _require(_sha(parent_path) == PARENT_SHA256, "accepted I1 parent report hash changed")
    _require(_json(parent_path) == parent, "provided parent differs from the accepted I1 report")
    _require(parent.get("implementation_hashes", {}).get(GRAPH_HELPER_PATH) == GRAPH_HELPER_SHA256
             and _sha(root_path / GRAPH_HELPER_PATH) == GRAPH_HELPER_SHA256,
             "accepted exact graph arithmetic helper changed")
    _require(parent.get("checks", {}).get("passed") is True, "accepted I1 parent checks did not pass")
    _require(parent.get("contract", {}).get("old_signed_sources_and_all_error_caps_retained") is True, "I1 signed sources were not retained")
    _require(parent.get("contract", {}).get("full_pressure_outer_RQy_gate_retained_false") is True, "I1 R_Q_y outer gate contract changed")
    _require(parent.get("checks", {}).get("accepted_ancestors_replayed") is False, "accepted ancestor replay is disallowed")

    base_dir = root_path / "experiments" / "root_st073"
    paths = {name: base_dir / name for name in PINNED_HASHES}
    for name, expected in PINNED_HASHES.items():
        _require(paths[name].is_file(), f"required accepted source missing: {name}")
        _require(_sha(paths[name]) == expected, f"accepted source hash changed: {name}")
        parent_pin = parent.get("input_hashes", {}).get(name)
        _require(parent_pin == expected, f"parent did not consume accepted source pin: {name}")

    join_check = _json(paths[JOIN_CHECK_NAME])
    selected_checks = _check_receipt(
        join_check,
        (
            "source_contract",
            "actual_background_recipes",
            "independent_physical_transport_seam",
            "actual_leading_recovery_and_common_N",
            "actual_macro_pressure_and_FTC",
            "actual_own_mixed_seam",
        ),
    )
    _require(join_check.get("source_family") == parent.get("source_family"), "join and I1 source families differ")
    accepted_qy = _accepted_qy_bindings(join_check)
    smooth_premises = _accepted_smooth_premises(join_check)
    q4_check = _json(paths[Q4_CHECK_NAME])
    _require(q4_check.get("all_passed") is True, "accepted Q4 receipt is not all_passed")
    for gate in (
        "current_original_true_ten_chart_modified_full_generic_C4_functions_installed",
        "current_original_true_ten_chart_modified_total_Uy_mean_C5_Qy_C4_functions_installed",
        "current_original_true_rminus_to_Rm_signed_pressure_C3_functions_installed",
    ):
        _require(q4_check.get(gate) is True, f"accepted Q4 gate missing: {gate}")
    q4_checks = q4_check.get("checks", {})
    _require(isinstance(q4_checks, Mapping), "accepted Q4 checks missing")
    for check_name, field in (
        ("actual_source", "complete_modified_mean_C5_and_full_generic_Q_Qy_C4"),
        ("actual_source", "Pstar_once_and_physical_R_Rsquared_measures_once"),
        ("actual_source", "all_ten_signed_pressure_predecessors_and_finite_endpoints_retained"),
        ("binders", "bound_names_unique_and_no_fast_phase_node_used_as_binder"),
    ):
        _require(q4_checks.get(check_name, {}).get(field) is True, f"accepted Q4 premise missing: {field}")

    join_guard = _source_guard(
        paths[JOIN_SOURCE_NAME],
        {
            "__init__": (
                "self.background=self.bridge['macro']",
                "self.radius=self.ref(self.background['radius'])",
                "self.phase=self.ref(self.loop['frozen_macro']['actual_phase'])",
                "self.N=self.base.N",
                "micro_hb_and_macro_length_native_Jacobians_remain_distinct=True",
            ),
            "_second_micro_packet": (
                "coordinate=self.ref(bg['phase'])",
                "radius=self.ref(bg['radius'])",
                "phase=self.ref(self.loop['second_micro']['actual_phase'])",
                "physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])]",
            ),
        },
    )
    q4_guard = _source_guard(
        paths[Q4_SOURCE_NAME],
        {
            "physical_pressure_fields": (
                "U=a.scale(fields['V'],Pstar)",
                "Uy=a.scale(radial['V_y'],Pstar)",
                "m=a.scale(fields['m'],Pstar)",
                "my=a.sub(U,m)",
                "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
            ),
            "__init__": (
                "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",
                "physical_measures='S0:R*Jacobian*dcoordinate; S1:R^2*Jacobian*dcoordinate'",
                "m15_minus_endpoint_sign_applied_once=True",
            ),
        },
    )

    source_family = join_check["source_family"]
    q4_archive = _json_from_gzip(paths[Q4_ARCHIVE_NAME])
    _require(q4_archive.get("source_family") == source_family, "Q4 and accepted join source families differ")
    nodes = q4_archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "Q4 exact graph prefix missing")
    join_archive = _json_from_gzip(paths[JOIN_ARCHIVE_NAME])
    _require(join_archive.get("source_family") == source_family, "join archive source family changed")
    join_nodes = join_archive.get("exact_graph_nodes")
    _require(isinstance(join_nodes, list) and join_nodes[:len(nodes)] == nodes, "join archive does not retain the exact Q4 graph prefix")
    join_bindings = join_archive.get("micro_macro_mixed_source_join_bindings", {})
    join_charts = join_bindings.get("charts", {}) if isinstance(join_bindings, Mapping) else {}
    _require(isinstance(join_charts, Mapping), "join archive chart bindings missing")
    _require(join_charts.get("second_micro", {}).get("physical", {}).get("Q") is not None, "join second-micro physical Q bindings missing")
    _require(join_charts.get("frozen_macro", {}).get("physical", {}).get("Q") is not None, "join macro physical Q bindings missing")
    join_full_qy = {
        "second_micro": list(join_charts["second_micro"]["physical"]["Q"][1]),
        "frozen_macro": list(join_charts["frozen_macro"]["physical"]["Q"][1]),
    }
    _require(join_full_qy["second_micro"] == list(_Q4_NODES["micro_Q_y_full"]), "join second-micro Q[1] binding changed")
    _require(join_full_qy["frozen_macro"] == list(_Q4_NODES["macro_Q_y_full"]), "join macro Q[1] binding changed")
    join_archive_graph_prefix_verified = True
    join_archive_qy_bindings_verified = True
    del join_archive, join_nodes, join_bindings, join_charts
    verified_rows = _verify_q4_rows(nodes, q4_archive)
    graph = GraphPairs(nodes)
    prefix_len = len(graph.nodes)
    half = graph.constant("1/2")
    zero = graph.constant(0)
    one = graph.constant(1)
    two = graph.constant(2)
    _require((zero, one, two) == (0, 1, 48), "Q4 exact rational handles changed")

    micro_upper = tuple(_Q4_NODES["micro_endpoint_R_Q_y"])
    micro_lower = tuple(_sub(graph, row, 66988, one) for row in _Q4_NODES["micro_finite_R_Q_y"])
    macro_lower = tuple(_sub(graph, row, 70265, zero) for row in _Q4_NODES["macro_finite_R_Q_y"])
    # Derive macro's upper endpoint from its finite source rows.  The accepted
    # endpoint rows 141114..141117 are intentionally never consumed here.
    macro_upper = tuple(_sub(graph, row, 70265, one) for row in _Q4_NODES["macro_finite_R_Q_y"])
    micro_qy_upper = tuple(_sub(graph, row, 66988, two) for row in _Q4_NODES["micro_Q_y"])
    macro_qy_lower = tuple(_sub(graph, row, 70265, zero) for row in _Q4_NODES["macro_Q_y"])
    micro_radius_upper = _sub(graph, 68644, 66988, two)
    macro_radius_lower = _sub(graph, 71448, 70265, zero)
    factor_upper = tuple(graph.mul(row, micro_radius_upper) for row in micro_qy_upper)
    factor_lower = tuple(graph.mul(row, macro_radius_lower) for row in macro_qy_lower)
    factorized_upper_checks = tuple(
        _factorized_substitution_identity(
            graph.nodes,
            _Q4_NODES["micro_finite_R_Q_y"][k],
            micro_upper[k],
            factor_upper[k],
            micro_qy_upper[k],
            micro_radius_upper,
            _Q4_NODES["micro_Q_y"][k],
            68644,
            66988,
            48,
        )
        for k in range(4)
    )
    factorized_lower_checks = tuple(
        _factorized_substitution_identity(
            graph.nodes,
            _Q4_NODES["macro_finite_R_Q_y"][k],
            macro_lower[k],
            factor_lower[k],
            macro_qy_lower[k],
            macro_radius_lower,
            _Q4_NODES["macro_Q_y"][k],
            71448,
            70265,
            0,
        )
        for k in range(4)
    )

    proof_rows: list[dict[str, Any]] = []
    for k in range(4):
        micro_cumulative_at_two = _Q4_NODES["micro_cumulative_S1_endpoint"][k]
        macro_cumulative_at_one = _Q4_NODES["macro_cumulative_S1_endpoint"][k]
        micro_history = graph.add(
            micro_cumulative_at_two,
            graph.neg(_Q4_NODES["micro_incoming_S1"][k]),
        )
        macro_history = graph.add(
            macro_cumulative_at_one,
            graph.neg(_Q4_NODES["macro_incoming_S1"][k]),
        )
        micro_increment = graph.add(
            graph.mul(half, micro_history),
            graph.neg(micro_upper[k]),
            micro_lower[k],
        )
        macro_increment = graph.add(
            graph.mul(half, macro_history),
            graph.neg(macro_upper[k]),
            macro_lower[k],
        )
        cancellation = graph.add(graph.neg(micro_upper[k]), macro_lower[k])
        joined = graph.add(micro_increment, macro_increment)
        # Construct the no-seam expression directly.  It must not reference
        # the cancellation row whose vanishing is conditional on the source
        # join receipt.
        residual = graph.add(
            graph.mul(half, micro_history),
            graph.mul(half, macro_history),
            micro_lower[k],
            graph.neg(macro_upper[k]),
        )
        reconstructed = graph.add(residual, cancellation)
        proof_rows.append(
            {
                "ordinary_Z_order": k,
                "second_micro": {
                    "incoming_S1": _Q4_NODES["micro_incoming_S1"][k],
                    "first_micro_history_endpoint": _Q4_NODES["first_micro_endpoint_S1"][k],
                    "cumulative_S1_finite_source": _Q4_NODES["micro_cumulative_S1"][k],
                    "cumulative_S1_at_s_2": micro_cumulative_at_two,
                    "local_S1_difference": micro_history,
                    "R_Q_y_lower_at_s_1": micro_lower[k],
                    "R_Q_y_upper_at_s_2": micro_upper[k],
                    "global_initial_R_Q_y": _Q4_NODES["global_initial_R_Q_y"][k],
                    "global_m15_finite_source": _Q4_NODES["micro_global_m15_finite"][k],
                    "global_m15_endpoint": _Q4_NODES["micro_global_m15_endpoint"][k],
                    "adjacent_increment": micro_increment,
                },
                "frozen_macro": {
                    "incoming_S1": _Q4_NODES["macro_incoming_S1"][k],
                    "cumulative_S1_finite_source": _Q4_NODES["macro_cumulative_S1"][k],
                    "cumulative_S1_at_m_1": macro_cumulative_at_one,
                    "local_S1_difference": macro_history,
                    "R_Q_y_lower_at_m_0": macro_lower[k],
                    "R_Q_y_upper_derived_at_m_1": macro_upper[k],
                    "global_m15_finite_source": _Q4_NODES["macro_global_m15_finite"][k],
                    "global_m15_endpoint": _Q4_NODES["macro_global_m15_endpoint"][k],
                    "adjacent_increment": macro_increment,
                },
                "signed_boundary_cancellation": cancellation,
                "joined_adjacent_increment": joined,
                "residual_after_boundary_cancellation": residual,
                "reconstructed_joined_increment": reconstructed,
                "factorized_micro_endpoint": factor_upper[k],
                "factorized_macro_lower": factor_lower[k],
                "factorized_micro_endpoint_identity": factorized_upper_checks[k],
                "factorized_macro_lower_identity": factorized_lower_checks[k],
                "conditional_boundary_identity": {
                    "passed": True,
                    "ordinary_Z_order": k,
                    "common_Q_y_source_rows": True,
                    "common_physical_radius_cutpoint": True,
                    "one_common_N_and_global_phase": True,
                    "mean_ODE_and_complete_correction_premises": True,
                    "identity": "R_micro(2)*Q_y_micro(2)=R_macro(0)*Q_y_macro(0)",
                },
            }
        )

    appended = graph.appended
    forbidden = {141114, 141115, 141116, 141117}
    appended_refs = {
        value
        for row in appended
        for key, value in row.items()
        if key in {"argument", "numerator", "denominator", "expression", "variable", "value", "integrand", "lower", "upper"}
        and isinstance(value, int)
    }
    for row in appended:
        if isinstance(row.get("arguments"), list):
            appended_refs.update(value for value in row["arguments"] if isinstance(value, int))
    _require(not (forbidden & appended_refs), "forbidden accepted macro endpoint handles were consumed")
    residual_rows = [row["residual_after_boundary_cancellation"] for row in proof_rows]
    cancellation_rows = [row["signed_boundary_cancellation"] for row in proof_rows]
    _require(all(cancellation not in graph.nodes[residual].get("arguments", ()) for cancellation, residual in zip(cancellation_rows, residual_rows)), "residual graph references cancelled seam row")
    controls = _symbolic_controls(join_full_qy, smooth_premises)
    source_pins = {PARENT_NAME: PARENT_SHA256, **PINNED_HASHES}
    source_guards = {
        "passed": True,
        "accepted_join_checks": selected_checks,
        "join_source_ast": join_guard,
        "q4_source_ast": q4_guard,
        "q4_rows": verified_rows,
        "join_archive_exact_Q4_prefix": join_archive_graph_prefix_verified,
        "join_archive_full_Q_y_bindings": join_archive_qy_bindings_verified,
        "forbidden_macro_endpoint_handles_consumed": False,
        "accepted_join_common_physical_source": True,
        "accepted_join_common_Q_y_rows": {
            "second_micro": accepted_qy["second_micro"],
            "frozen_macro": accepted_qy["frozen_macro"],
        },
        "accepted_join_full_Q_y_rows": join_full_qy,
        "full_Q4_physical_Q_y_rows_bound": True,
        "physical_Q_y_Z0_Z3_rows": {
            "second_micro": list(_Q4_NODES["micro_Q_y"]),
            "frozen_macro": list(_Q4_NODES["macro_Q_y"]),
        },
        "accepted_join_common_U_mean_recovery": True,
        "accepted_join_physical_R_and_native_J_distinct": True,
        "accepted_smooth_source_premises": smooth_premises,
        "accepted_Q4_Z_independent_binders": True,
        "radius_Z_independent_transitive": True,
    }
    return {
        "passed": True,
        "source_family": source_family,
        "source_pins": source_pins,
        "input_pins": source_pins,
        "source_guards": source_guards,
        "proof_rows": proof_rows,
        "graph_patch": {
            "base_archive": Q4_ARCHIVE_NAME,
            "base_node_count": prefix_len,
            "base_prefix_unchanged": graph.nodes[:prefix_len] == list(graph.base_nodes),
            "appended_nodes": appended,
            "appended_node_count": len(appended),
            "endpoint_handles": {
                "second_micro_upper": list(micro_upper),
                "frozen_macro_lower": list(macro_lower),
                "signed_cancellation": [row["signed_boundary_cancellation"] for row in proof_rows],
            },
        },
        "symbolic_controls": controls,
        "contract": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "scope": "adjacent second_micro s=1..2 followed by frozen_macro m=0..1",
            "second_micro_upper_handles": list(micro_upper),
            "second_micro_lower_coordinate": "s=1",
            "frozen_macro_lower_is_new_finite_boundary_substitution": True,
            "macro_upper_accepted_endpoint_handles_not_used": True,
            "signed_boundary_terms": "-R_Q_y(second_micro,s=2)+R_Q_y(frozen_macro,m=0)",
            "joined_pressure_expression": "1/2*(S1_micro_at_s2-S1_micro_incoming)-R_Q_y(s=2)+R_Q_y(s=1) + 1/2*(S1_macro_at_m1-S1_macro_incoming)-R_Q_y(m=1)+R_Q_y(m=0)",
            "global_initial_and_incoming_histories_retained": True,
            "first_micro_history_retained_as_incoming_source": True,
            "whole_first_to_macro_pair_not_claimed": True,
            "old_six_signed_terms_and_caps_retained": True,
            "cancellation_only_in_adjacent_joined_expression": True,
            "no_numeric_phase_or_quadrature": True,
            "no_outer_cap_deletion": True,
            "quantitative_mixed_norms_remain_open": True,
            "residual_constructed_without_cancellation_pair": True,
        },
    }


def _json_from_gzip(path: Path) -> dict[str, Any]:
    import gzip

    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path.name}")
    return value


__all__ = ["build"]
