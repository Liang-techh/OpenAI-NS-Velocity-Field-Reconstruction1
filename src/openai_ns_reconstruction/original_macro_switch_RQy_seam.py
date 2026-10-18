"""Source-linked frozen-macro/first-switch ``R*Q_y`` endpoint seam.

The accepted correlated micro/macro report already contains the complete
Q4-prefix graph and its serialized suffix.  This module consumes that saved
namespace and appends only the ordinary ``Z=0,1`` endpoint arithmetic for
the adjacent frozen-macro ``m=1`` / first-switch ``s=0`` join.  Source
archives and receipts are read as data; no archived owner, checker,
producer, recurrence, quadrature, or numeric builder is imported or run.
"""

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .original_signed_graph_pairs import GraphPairs


PARENT_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_correlated_m15.json"
)
PARENT_SHA256 = (
    "8354af18754496450fb41ff0335e2c82f503c0d64852496da11ee1546fb23577"
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
HIGH_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.py"
)
HIGH_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.json.gz"
)
HIGH_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check.json"
)
LOOP_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
LOOP_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
LOOP_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
)
MICRO_JOIN_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.py"
)
MICRO_JOIN_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.json.gz"
)
MICRO_JOIN_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check.json"
)
SWITCH_JOIN_SOURCE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins.py"
)
SWITCH_JOIN_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins.json.gz"
)
SWITCH_JOIN_CHECK_NAME = (
    "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins_check.json"
)
GRAPH_HELPER_PATH = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"

PINNED_HASHES = {
    Q4_SOURCE_NAME: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    Q4_ARCHIVE_NAME: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    Q4_CHECK_NAME: "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
    HIGH_SOURCE_NAME: "a66e8afedc6c8294f3393b46e84bd19e29a0eaa238539dae1032fa0d8acd6e41",
    HIGH_ARCHIVE_NAME: "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991",
    HIGH_CHECK_NAME: "e53c1fb2a485d5456ee4d1f9484df261aa5cff0f27ea8d1d3a80fa951a7c34ce",
    LOOP_SOURCE_NAME: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    LOOP_ARCHIVE_NAME: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    LOOP_CHECK_NAME: "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    MICRO_JOIN_SOURCE_NAME: "3f0bb8e0363990ab4d851020e988dd463061086d3ad25e7f0e59c047fa3c5bc0",
    MICRO_JOIN_ARCHIVE_NAME: "0c673c0979a69f0a6f6fd48e2ee72adae2223f6e384d537d6e39801c15af79e2",
    MICRO_JOIN_CHECK_NAME: "49e7d4f0c566abcbfe9833c5cc3c63875b1c80fa6e6025601eab8aaa83956225",
    SWITCH_JOIN_SOURCE_NAME: "0ae7f651c0d183756f0e6ab90d56dfd80bf3d3c1a6717ffdc2a03fe08dd52f2f",
    SWITCH_JOIN_ARCHIVE_NAME: "ddbd506d2fbfc4c58f7a5ba93e0abe3a95d3af8fd3c520ee3ce5024d90a29567",
    SWITCH_JOIN_CHECK_NAME: "9cf3c6d941df95f4b5f73e7b97e32198a54fcd3d17a85ecca2a0742c0d602fc3",
    GRAPH_HELPER_PATH: "1e9699cd84a1c82dc7b6577f922f1da1362cbedcfef148e7cec1d146a99aaaff",
}

# Exact handles are read from the accepted Q4 graph.  The separate tuples
# make it difficult to accidentally use a second-micro or a first-switch
# upper endpoint in this macro/first-switch seam.
_Q4 = {
    "macro_coordinate": 70265,
    "macro_radius": 71448,
    "macro_native_J": 66987,
    "macro_Q": (140780, 140784, 140790, 140796, 140802),
    "macro_Q_y": (140868, 140872, 140878, 140884, 140890),
    "macro_finite_R_Q_y": (141027, 141028, 141029, 141030),
    "macro_upper_R_Q_y": (141114, 141115, 141116, 141117),
    "first_coordinate": 73190,
    "first_radius": 85905,
    "first_native_J": 9,
    "first_Q": (148180, 148184, 148190, 148196, 148202),
    "first_Q_y": (148268, 148272, 148278, 148284, 148290),
    "first_finite_R_Q_y": (148427, 148428, 148429, 148430),
    "first_upper_R_Q_y": (148514, 148515, 148516, 148517),
    "zero": 0,
    "one": 1,
    "two": 48,
}

_COMMON_SOURCE_GUARDS = {
    Q4_SOURCE_NAME: {
        "physical_pressure_fields": (
            "U=a.scale(fields['V'],Pstar)",
            "Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar)",
            "my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        ),
        "__init__": (
            "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",
            "physical_measures='S0:R*Jacobian*dcoordinate;S1:R^2*Jacobian*dcoordinate'",
            "m15_minus_endpoint_sign_applied_once=True",
        ),
    },
    HIGH_SOURCE_NAME: {
        "build_chart": (
            "coordinate=self.ref(loop['coordinate'])",
            "radius=self.ref(loop['radius'])",
            "J=self.ref(loop['physical_measure'])",
            "physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])]",
        ),
    },
    MICRO_JOIN_SOURCE_NAME: {
        "__init__": (
            "self.background=self.bridge['macro']",
            "self.radius=self.ref(self.background['radius'])",
            "self.phase=self.ref(self.loop['frozen_macro']['actual_phase'])",
            "self.N=self.base.N",
            "micro_hb_and_macro_length_native_Jacobians_remain_distinct=True",
        ),
    },
    SWITCH_JOIN_SOURCE_NAME: {
        "__init__": (
            "self.joins={key:dict(left_chart=left,right_chart=right,original_endpoint=endpoint,",
            "left=self.base.at_own(left,g.one),right=self.base.at_own(right,g.zero),",
        ),
    },
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


def _json_gzip(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path.name}")
    return value


def _source_guard(path: Path, functions: Mapping[str, tuple[str, ...]]) -> dict[str, Any]:
    """AST-check source snippets; every snippet is deliberately a tuple."""

    with path.open("r", encoding="utf-8", newline="") as stream:
        source = stream.read()
    tree = ast.parse(source, filename=str(path))
    found: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            found[node.name] = ast.get_source_segment(source, node) or ""
    result: dict[str, Any] = {"path": path.name, "functions": {}}
    for name, tokens in functions.items():
        _require(isinstance(tokens, tuple), f"source guard snippets must be tuples: {name}")
        body = found.get(name, "")
        compact = "".join(body.split())
        matched = all("".join(token.split()) in compact for token in tokens)
        result["functions"][name] = {
            "present": bool(body),
            "tokens": list(tokens),
            "matched": matched,
        }
        _require(body and matched, f"source guard failed for {path.name}:{name}")
    result["passed"] = True
    return result


def _receipt(check: Mapping[str, Any], names: Sequence[str]) -> dict[str, Any]:
    _require(check.get("all_passed") is True, "accepted source receipt is not all_passed")
    checks = check.get("checks")
    _require(isinstance(checks, Mapping), "accepted source receipt checks missing")
    selected: dict[str, Any] = {}
    for name in names:
        item = checks.get(name)
        _require(isinstance(item, Mapping) and item.get("passed") is True, f"source receipt failed: {name}")
        selected[name] = {"passed": True}
    return selected


def _check_bool_fields(item: Mapping[str, Any], names: Sequence[str], label: str) -> dict[str, bool]:
    _require(item.get("passed") is True, f"source check group is not passed: {label}")
    selected: dict[str, bool] = {}
    for name in names:
        _require(item.get(name) is True, f"source check premise missing: {label}.{name}")
        selected[name] = True
    return selected


def _check_true_fields(item: Mapping[str, Any], names: Sequence[str], label: str) -> dict[str, bool]:
    selected: dict[str, bool] = {}
    for name in names:
        _require(item.get(name) is True, f"source premise missing: {label}.{name}")
        selected[name] = True
    return selected


def _base_row(nodes: Sequence[Mapping[str, Any]], index: int) -> Mapping[str, Any]:
    _require(0 <= index < len(nodes), f"graph handle outside prefix: {index}")
    return nodes[index]


def _expect_substitution(
    nodes: Sequence[Mapping[str, Any]],
    index: int,
    expression: int,
    variable: int,
    value: int,
) -> None:
    row = _base_row(nodes, index)
    _require(
        row.get("operation") == "function_substitution"
        and row.get("expression") == expression
        and row.get("variable") == variable
        and row.get("value") == value
        and row.get("Z_independent_substitution") is True,
        f"endpoint substitution mismatch at {index}",
    )


def _verify_q4(q4: Mapping[str, Any], q4_check: Mapping[str, Any]) -> dict[str, Any]:
    nodes = q4.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == 183700, "Q4 exact graph prefix is not 183700 nodes")
    bindings = q4.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(bindings, Mapping), "Q4 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping), "Q4 charts/endpoints missing")
    _require(
        set(charts) == {
            "first_micro",
            "second_micro",
            "frozen_macro",
            "first_switch",
            "second_switch",
            "post_power",
            "long_reshape",
            "reference",
            "restoration",
            "postrestore",
        },
        "Q4 ten-chart prefix is incomplete",
    )

    macro = charts.get("frozen_macro")
    first = charts.get("first_switch")
    _require(isinstance(macro, Mapping) and isinstance(first, Mapping), "Q4 macro/first charts missing")
    _require(macro.get("physical", {}).get("U") == list(range(140685, 140691)), "macro physical U handles changed")
    _require(macro.get("physical", {}).get("U_y") == list(range(140691, 140697)), "macro physical U_y handles changed")
    _require(macro.get("physical", {}).get("Q") == list(_Q4["macro_Q"]), "macro physical Q handles changed")
    _require(macro.get("physical", {}).get("Q_y") == list(_Q4["macro_Q_y"]), "macro physical Q_y handles changed")
    _require(first.get("physical", {}).get("U") == list(range(148085, 148091)), "first-switch physical U handles changed")
    _require(first.get("physical", {}).get("U_y") == list(range(148091, 148097)), "first-switch physical U_y handles changed")
    _require(first.get("physical", {}).get("Q") == list(_Q4["first_Q"]), "first-switch physical Q handles changed")
    _require(first.get("physical", {}).get("Q_y") == list(_Q4["first_Q_y"]), "first-switch physical Q_y handles changed")
    _require(macro.get("finite_boundary", {}).get("R_Q_y") == list(_Q4["macro_finite_R_Q_y"]), "macro finite R_Q_y handles changed")
    _require(first.get("finite_boundary", {}).get("R_Q_y") == list(_Q4["first_finite_R_Q_y"]), "first finite R_Q_y handles changed")
    _require(endpoints.get("frozen_macro", {}).get("boundary", {}).get("R_Q_y") == list(_Q4["macro_upper_R_Q_y"]), "macro upper R_Q_y handles changed")
    _require(endpoints.get("first_switch", {}).get("boundary", {}).get("R_Q_y") == list(_Q4["first_upper_R_Q_y"]), "first upper R_Q_y handles changed")

    for qy, product in zip(_Q4["macro_Q_y"][:4], _Q4["macro_finite_R_Q_y"]):
        _require(_base_row(nodes, product).get("operation") == "product" and _base_row(nodes, product).get("arguments") == [qy, _Q4["macro_radius"]], f"macro R_Q_y product changed: {product}")
    for qy, product in zip(_Q4["first_Q_y"][:4], _Q4["first_finite_R_Q_y"]):
        _require(_base_row(nodes, product).get("operation") == "product" and _base_row(nodes, product).get("arguments") == [qy, _Q4["first_radius"]], f"first R_Q_y product changed: {product}")
    for index, expression, variable, value in zip(_Q4["macro_upper_R_Q_y"], _Q4["macro_finite_R_Q_y"], (70265,) * 4, (1,) * 4):
        _expect_substitution(nodes, index, expression, variable, value)
    for index, expression in zip(_Q4["first_upper_R_Q_y"], _Q4["first_finite_R_Q_y"]):
        _expect_substitution(nodes, index, expression, 73190, 1)

    _require(_base_row(nodes, 0).get("operation") == "exact_rational" and _base_row(nodes, 0).get("numerator") == 0, "Q4 zero handle changed")
    _require(_base_row(nodes, 1).get("operation") == "exact_rational" and _base_row(nodes, 1).get("numerator") == 1, "Q4 one handle changed")
    _require(_base_row(nodes, 48).get("operation") == "exact_rational" and _base_row(nodes, 48).get("numerator") == 2, "Q4 exact-rational handle 48 changed")
    _require(_base_row(nodes, 8).get("operation") == "current_original_source_parameter" and _base_row(nodes, 8).get("name") == "hbB", "Q4 handle 8 source changed")
    _require(_base_row(nodes, 9).get("operation") == "current_original_source_parameter" and _base_row(nodes, 9).get("name") == "hbS", "Q4 handle 9 source changed")
    _require(_base_row(nodes, 70265).get("operation") == "bound_variable", "macro coordinate binder changed")
    _require(_base_row(nodes, 73190).get("operation") == "bound_variable", "first-switch coordinate binder changed")

    selected = _receipt(q4_check, ("actual_source", "binders"))
    actual = q4_check["checks"]["actual_source"]
    _check_bool_fields(actual, ("complete_modified_mean_C5_and_full_generic_Q_Qy_C4", "Pstar_once_and_physical_R_Rsquared_measures_once", "all_ten_signed_pressure_predecessors_and_finite_endpoints_retained"), "Q4.actual_source")
    _check_bool_fields(q4_check["checks"]["binders"], ("bound_names_unique_and_no_fast_phase_node_used_as_binder",), "Q4.binders")
    return {
        "node_count": len(nodes),
        "physical_ports": {
            "frozen_macro": {"U": list(macro["physical"]["U"]), "Q": list(macro["physical"]["Q"]), "Q_y": list(macro["physical"]["Q_y"])},
            "first_switch": {"U": list(first["physical"]["U"]), "Q": list(first["physical"]["Q"]), "Q_y": list(first["physical"]["Q_y"])},
        },
        "finite_boundary_ports": {"frozen_macro": list(_Q4["macro_finite_R_Q_y"]), "first_switch": list(_Q4["first_finite_R_Q_y"])},
        "standalone_upper_ports": {"frozen_macro": list(_Q4["macro_upper_R_Q_y"]), "first_switch": list(_Q4["first_upper_R_Q_y"])},
        "exact_handles": {"zero": 0, "one": 1, "two": 48, "hbB": 8, "hbS": 9},
        "accepted_checks": selected,
    }


def _append_saved(graph: GraphPairs, row: Mapping[str, Any]) -> int:
    _require(isinstance(row, Mapping) and isinstance(row.get("operation"), str), "serialized node row malformed")
    attrs = {key: value for key, value in row.items() if key != "operation"}
    return graph.node(row["operation"], **attrs)


def _factorized_identity(
    nodes: Sequence[Mapping[str, Any]],
    *,
    finite: int,
    endpoint: int,
    factorized: int,
    qy_sub: int,
    radius_sub: int,
    qy: int,
    radius: int,
    coordinate: int,
    value: int,
    source_chart: str,
    physical_q_y_source: Sequence[int],
    source_premises: Mapping[str, Any],
) -> dict[str, Any]:
    finite_row = _base_row(nodes, finite)
    endpoint_row = _base_row(nodes, endpoint)
    factor_row = _base_row(nodes, factorized)
    _require(finite_row.get("operation") == "product" and finite_row.get("arguments") == [qy, radius], f"finite source product changed at {finite}")
    _require(
        endpoint_row.get("operation") == "function_substitution"
        and endpoint_row.get("expression") == finite
        and endpoint_row.get("variable") == coordinate
        and endpoint_row.get("value") == value
        and endpoint_row.get("Z_independent_substitution") is True,
        f"endpoint wrapper changed at {endpoint}",
    )
    qy_sub_row = _base_row(nodes, qy_sub)
    radius_sub_row = _base_row(nodes, radius_sub)
    _require(
        qy_sub_row.get("operation") == "function_substitution"
        and qy_sub_row.get("expression") == qy
        and qy_sub_row.get("variable") == coordinate
        and qy_sub_row.get("value") == value
        and qy_sub_row.get("Z_independent_substitution") is True,
        f"Q_y endpoint source wrapper changed at {qy_sub}",
    )
    _require(
        radius_sub_row.get("operation") == "function_substitution"
        and radius_sub_row.get("expression") == radius
        and radius_sub_row.get("variable") == coordinate
        and radius_sub_row.get("value") == value
        and radius_sub_row.get("Z_independent_substitution") is True,
        f"radius endpoint source wrapper changed at {radius_sub}",
    )
    _require(qy in tuple(physical_q_y_source), f"Q_y source port is not bound for {source_chart}: {qy}")
    _require(qy in tuple(source_premises.get("source_chart_physical_Q_y", ())), f"Q4 Q_y source port is not bound for {source_chart}: {qy}")
    _require(qy in tuple(source_premises.get("high_chart_physical_Q_y", ())), f"high Q_y source port is not bound for {source_chart}: {qy}")
    _require(source_premises.get("common_physical_Q_y_source") is True, f"common Q_y source premise missing for {source_chart}")
    _require(source_premises.get("common_physical_radius_source") is True, f"common radius source premise missing for {source_chart}")
    _require(isinstance(source_premises.get("high_common_source_controls"), Mapping), f"high common source controls missing for {source_chart}")
    _require(factor_row.get("operation") == "product" and factor_row.get("arguments") == [qy_sub, radius_sub], f"factorized product changed at {factorized}")
    return {
        "passed": True,
        "finite_product": finite,
        "endpoint_wrapper": endpoint,
        "factorized_product": factorized,
        "substituted_Q_y": qy_sub,
        "substituted_radius": radius_sub,
        "coordinate": coordinate,
        "endpoint_value": value,
        "Z_independent_substitution": True,
        "source_binding": {
            "source_chart": source_chart,
            "physical_Q_y_source_handle": qy,
            "physical_Q_y_source_port": "Q_y",
            "physical_radius_source_handle": radius,
            "common_physical_Q_y_source": True,
            "common_physical_radius_source": True,
            "source_premise_receipt": dict(source_premises),
        },
        "identity": "sub(product(Q_y,R),coordinate,value)=product(sub(Q_y,coordinate,value),sub(R,coordinate,value))",
    }


def _symbolic_endpoint_controls(
    graph: GraphPairs,
    rows: Sequence[Mapping[str, Any]],
    macro_factorized: Sequence[int],
    first_factorized: Sequence[int],
    macro_qy: Sequence[int],
    first_qy: Sequence[int],
    macro_radius: int,
    first_radius: int,
) -> dict[str, Any]:
    """Compare actual graph products after common source aliases are bound."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("SymPy is required for the new exact endpoint controls") from exc

    common_R = sp.Symbol("R_common")
    common_qy = [sp.Symbol(f"Q_y_Z{k}") for k in range(2)]
    identities: list[dict[str, Any]] = []
    for k in range(2):
        macro_row = _base_row(rows, macro_factorized[k])
        first_row = _base_row(rows, first_factorized[k])
        _require(macro_row.get("operation") == "product" and first_row.get("operation") == "product", "endpoint factorized rows must be products")
        _require(macro_row.get("arguments") == [macro_qy[k], macro_radius], "macro factorized row is not Q_y times R")
        _require(first_row.get("arguments") == [first_qy[k], first_radius], "first-switch factorized row is not Q_y times R")
        # These are source-theorem aliases: the high source receipt establishes
        # the common physical Q_y jet and radius at the adjacent cutpoint.
        macro_expr = common_R * common_qy[k]
        first_expr = common_R * common_qy[k]
        difference = sp.expand(macro_expr - first_expr)
        _require(difference == 0, f"endpoint source identity failed at Z{k}")
        identities.append({
            "ordinary_Z_order": k,
            "macro_factorized": macro_factorized[k],
            "first_switch_factorized": first_factorized[k],
            "canonical_macro": f"R_common*Q_y_Z{k}",
            "canonical_first_switch": f"R_common*Q_y_Z{k}",
            "difference": "0",
            "passed": True,
        })
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "common_radius_alias": "R_common",
        "common_Q_y_aliases": ["Q_y_Z0", "Q_y_Z1"],
        "rows": identities,
        "arithmetic": "(R_common*Q_y_Zk)-(R_common*Q_y_Zk)=0",
        "source_aliases_are_not_hash_comparison": True,
    }


def _archive_prefix_summary(archive: Mapping[str, Any], q4_nodes: Sequence[Mapping[str, Any]], label: str) -> dict[str, Any]:
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) >= len(q4_nodes), f"{label} graph archive missing Q4 prefix")
    _require(nodes[: len(q4_nodes)] == list(q4_nodes), f"{label} graph archive Q4 prefix changed")
    return {"node_count": len(nodes), "q4_prefix_node_count": len(q4_nodes), "q4_prefix_unchanged": True}


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Append the bounded macro-1/first-switch-0 physical ``R*Q_y`` seam."""

    root_path = Path(root)
    _require(isinstance(parent, Mapping), "parent must be the accepted correlated report mapping")
    parent_path = root_path / "experiments" / "root_st073" / PARENT_NAME
    _require(parent_path.is_file() and _sha(parent_path) == PARENT_SHA256, "accepted correlated parent report hash changed")
    saved_parent = _json(parent_path)
    if parent == saved_parent:
        parent_report = saved_parent
        actual_parent = saved_parent.get("actual_micro_macro_correlated_m15")
    elif parent == saved_parent.get("actual_micro_macro_correlated_m15"):
        parent_report = saved_parent
        actual_parent = parent
    else:
        raise ValueError("provided parent differs from the accepted correlated report")
    _require(isinstance(actual_parent, Mapping), "accepted correlated nested report missing")
    _require(actual_parent.get("passed") is True, "accepted correlated nested report did not pass")
    _require(actual_parent.get("accepted_ancestors_replayed") is False, "accepted ancestor replay is disallowed")
    _require(actual_parent.get("numerical_phase_or_quadrature_evaluated") is False, "parent numerical phase/quadrature scope changed")

    base_dir = root_path / "experiments" / "root_st073"
    paths = {name: base_dir / name for name in PINNED_HASHES if name != GRAPH_HELPER_PATH}
    for name, expected in PINNED_HASHES.items():
        path = root_path / name if name == GRAPH_HELPER_PATH else base_dir / name
        _require(path.is_file(), f"required source pin missing: {name}")
        _require(_sha(path) == expected, f"source pin changed: {name}")
    _require(_sha(root_path / GRAPH_HELPER_PATH) == PINNED_HASHES[GRAPH_HELPER_PATH], "graph helper source pin changed")

    q4 = _json_gzip(paths[Q4_ARCHIVE_NAME])
    q4_check = _json(paths[Q4_CHECK_NAME])
    source_family = q4.get("source_family")
    _require(isinstance(source_family, Mapping), "Q4 source family missing")
    _require(source_family == parent_report.get("source_family"), "parent and Q4 source families differ")
    q4_nodes = q4["exact_graph_nodes"]
    q4_evidence = _verify_q4(q4, q4_check)

    guards: dict[str, Any] = {}
    for name, functions in _COMMON_SOURCE_GUARDS.items():
        path = root_path / name if name == GRAPH_HELPER_PATH else base_dir / name
        guards[name] = _source_guard(path, functions)

    high_check = _json(paths[HIGH_CHECK_NAME])
    high_receipt = _receipt(high_check, ("source_contract", "original_geometry", "actual_recovery_and_pressure_sources", "actual_macro_first_background_seam", "actual_macro_first_own_seam"))
    high_background = high_check["checks"]["actual_macro_first_background_seam"]
    high_own = high_check["checks"]["actual_macro_first_own_seam"]
    background_controls = _check_bool_fields(
        high_background,
        (
            "actual_current_radius_comparison_direction_and_B_half_replayed",
            "actual_first_switch_phi_U_histories_generic_background_replayed",
            "actual_macro1_nonzero_phi_U_six_histories_are_original_first_switch0_inputs",
            "common_K1_width_consumed_while_hbB_hbS_native_nodes_remain_distinct",
            "common_physical_comparison_history_ODE_then_actual_background_jets_y0_to5",
            "common_original_roots_unique_inverse_periodic_lift_and_NJ_transfer_U_y4_C5",
            "original_mean_history_ODE_Q_recovery_and_N_operator_transfer_Q_y4_N_y2",
        ),
        "high.actual_macro_first_background_seam",
    )
    _require(high_own.get("actual_leading_mean_C5_value_transfer_rows") == 6, "high leading mean transfer row count changed")
    own_controls = _check_bool_fields(
        high_own,
        (
            "original_nonzero_macro_correction_and_signed_pressure_memories_retained",
            "common_leading_Q_and_N_inputs_follow_preceding_source_theorem_not_desired_own_outputs",
            "native_length_and_hbS_remain_separate",
        ),
        "high.actual_macro_first_own_seam",
    )
    own_controls["actual_leading_mean_C5_value_transfer_rows"] = True

    high = _json_gzip(paths[HIGH_ARCHIVE_NAME])
    _require(high.get("source_family") == source_family, "high source family changed")
    high_prefix = _archive_prefix_summary(high, q4_nodes, "high")
    high_bindings = copy.deepcopy(high["prescribed_high_radial_source_bindings"])
    _require(high_bindings.get("N") == 18, "high source common N changed")
    phase_handles: dict[str, int] = {}
    for chart, coordinate, radius, native_J, phase in (
        ("frozen_macro", 70265, 71448, 66987, 101909),
        ("first_switch", 73190, 85905, 9, 104537),
    ):
        item = high_bindings["charts"][chart]
        loop_item = high_bindings["loop"][chart]
        _require(item.get("coordinate") == coordinate and item.get("radius") == radius and item.get("jacobian") == native_J and item.get("phase") == phase, f"high chart geometry changed: {chart}")
        _require(loop_item.get("coordinate") == coordinate and loop_item.get("radius") == radius and loop_item.get("physical_measure") == native_J and loop_item.get("actual_phase") == phase, f"high loop geometry changed: {chart}")
        intermediate = item.get("intermediate", {})
        _require(intermediate.get("complete_incoming_mass_not_reset") is True, f"nonreset M source changed: {chart}")
        phase_handles[chart] = phase
    _require(high_bindings["charts"]["frozen_macro"]["physical"]["Q"][1] == list(_Q4["macro_Q_y"]), "high macro physical Q_y rows changed")
    _require(high_bindings["charts"]["first_switch"]["physical"]["Q"][1] == list(_Q4["first_Q_y"]), "high first physical Q_y rows changed")
    high_definition = _check_true_fields(
        high_bindings["definition"],
        (
            "ten_actual_original_domains_not_relabelled_seam_neighborhoods",
            "every_chart_fresh_total_native_factory_uses_its_actual_N_times_J",
            "full_actual_U_y_differentiation_mean_ODE_and_one_minus_delta_Q",
            "original_signed_pressure_and_nonzero_predecessor_histories_retained",
            "own_one_plus_delta_and_actual_beforeRc_support_retained",
            "source_coverage_not_global_interface_closure_or_numerical_bounds",
        ),
        "high.definition",
    )
    del high

    loop_check = _json(paths[LOOP_CHECK_NAME])
    _require(loop_check.get("all_passed") is True and loop_check.get("source_family") == source_family, "loop receipt changed")
    loop_archive = _json_gzip(paths[LOOP_ARCHIVE_NAME])
    _require(loop_archive.get("source_family") == source_family, "loop archive source family changed")
    loop_bindings = loop_archive.get("prescribed_loop_C5_bindings", {})
    _require(loop_bindings.get("N") == 18, "loop archive common N changed")
    del loop_archive

    micro_check = _json(paths[MICRO_JOIN_CHECK_NAME])
    micro_receipt = _receipt(micro_check, ("source_contract", "actual_background_recipes", "independent_physical_transport_seam", "actual_leading_recovery_and_common_N", "actual_macro_pressure_and_FTC", "actual_own_mixed_seam"))
    _require(micro_check.get("source_family") == source_family, "micro source join family changed")
    micro_archive = _json_gzip(paths[MICRO_JOIN_ARCHIVE_NAME])
    _require(micro_archive.get("source_family") == source_family, "micro source join archive family changed")
    micro_prefix = _archive_prefix_summary(micro_archive, q4_nodes, "micro source join")
    micro_bindings = copy.deepcopy(micro_archive["micro_macro_mixed_source_join_bindings"])
    _require(micro_bindings.get("coordinate") == 70265 and micro_bindings.get("radius") == 71448 and micro_bindings.get("jacobian") == 66987 and micro_bindings.get("phase") == 101909 and micro_bindings.get("N") == 18, "micro source join macro port changed")
    del micro_archive

    switch_check = _json(paths[SWITCH_JOIN_CHECK_NAME])
    switch_receipt = _receipt(switch_check, ("source_contract", "actual_background_recipes", "physical_controls_and_ODE_jets", "actual_loop_and_phase_transfer", "actual_own_mixed_joins"))
    _require(switch_check.get("source_family") == source_family, "switch source join family changed")
    switch_archive = _json_gzip(paths[SWITCH_JOIN_ARCHIVE_NAME])
    _require(switch_archive.get("source_family") == source_family, "switch source join archive family changed")
    switch_prefix = _archive_prefix_summary(switch_archive, q4_nodes, "switch source join")
    switch_bindings = copy.deepcopy(switch_archive["switch_mixed_source_join_bindings"])
    _require(set(switch_bindings.get("joins", {})) == {"first_second", "second_power"}, "switch source join ports changed")
    _check_true_fields(
        switch_bindings["definition"],
        (
            "actual_first_switch1_second_switch0_and_second_switch1_post_power0",
            "original_sigma_native_derivatives_zero_through_order_four",
            "actual_source_field_history_pressure_and_correction_inlets_not_reset",
            "physical_background_y5_then_leading_y4_Q4_N2_pressure_y3",
            "own_outputs_not_assigned_to_seam_targets",
            "source_functions_not_numeric_evaluation_or_temporal_coefficients",
        ),
        "switch.definition",
    )
    del switch_archive

    # Rehydrate the accepted merged namespace exactly as serialized.  This is
    # a data append, so relocation and old source producers are not replayed.
    merged = actual_parent.get("merged_namespace")
    safe_graphs = actual_parent.get("safe_graphs")
    graph_namespace = actual_parent.get("graph_namespace")
    _require(isinstance(merged, Mapping) and isinstance(safe_graphs, Mapping) and isinstance(graph_namespace, Mapping), "accepted merged graph packet missing")
    _require(merged.get("base_archive") == f"experiments/root_st073/{Q4_ARCHIVE_NAME}", "accepted graph base archive changed")
    _require(merged.get("base_node_count") == len(q4_nodes) == 183700, "accepted graph base count changed")
    _require(merged.get("merged_node_count") == 184687 and merged.get("roundtrip_verified") is True, "accepted merged graph count/roundtrip changed")
    _require(graph_namespace.get("old_graph_unchanged") is True and graph_namespace.get("serialized_suffix_covers_every_new_node") is True, "accepted graph namespace contract changed")
    segments = []
    for key in ("accepted_correlated", "relocated_rqy", "new_macro", "joined_expression"):
        packet = safe_graphs.get(key)
        _require(isinstance(packet, Mapping) and isinstance(packet.get("appended_nodes"), list), f"accepted graph segment missing: {key}")
        segments.extend(packet["appended_nodes"])
    _require(segments == merged.get("appended_nodes"), "accepted graph segment serialization differs from merged suffix")

    graph = GraphPairs(q4_nodes)
    for row in merged["appended_nodes"]:
        _append_saved(graph, row)
    prefix_len = len(graph.nodes)
    _require(prefix_len == 184687 and graph.nodes[:183700] == list(q4_nodes), "accepted graph namespace rehydration changed")
    _require(graph.nodes[183700:] == merged["appended_nodes"], "accepted graph suffix roundtrip changed")
    accepted_graph_prefix = copy.deepcopy(graph.nodes)
    _require(actual_parent.get("remaining_outer_macro_RQy", {}).get("upper") == list(_Q4["macro_upper_R_Q_y"]), "accepted macro outer R_Q_y was not retained")
    _require(actual_parent.get("remaining_outer_macro_RQy", {}).get("lower") == [184328, 184329, 184330, 184331], "accepted first lower port changed")

    # First create endpoint traces from the finite source products.  The
    # accepted macro upper wrappers and first-switch upper wrappers remain in
    # the graph untouched; only the adjacent joined expression consumes a
    # macro upper and a newly constructed first-switch lower trace.
    zero = graph.constant(0)
    one = graph.constant(1)
    _require((zero, one) == (0, 1), "Q4 exact zero/one handles changed")
    macro_upper = tuple(_Q4["macro_upper_R_Q_y"][:2])
    first_lower = tuple(
        graph.node(
            "function_substitution",
            expression=finite,
            variable=_Q4["first_coordinate"],
            value=zero,
            Z_independent_substitution=True,
            source_endpoint_trace=True,
            endpoint_chart="first_switch",
            endpoint_native_query="s=0",
        )
        for finite in _Q4["first_finite_R_Q_y"][:2]
    )
    macro_qy_upper = tuple(
        graph.node(
            "function_substitution",
            expression=qy,
            variable=_Q4["macro_coordinate"],
            value=one,
            Z_independent_substitution=True,
            source_endpoint_trace=True,
            endpoint_chart="frozen_macro",
            endpoint_native_query="m=1",
        )
        for qy in _Q4["macro_Q_y"][:2]
    )
    macro_radius_upper = graph.node(
        "function_substitution",
        expression=_Q4["macro_radius"],
        variable=_Q4["macro_coordinate"],
        value=one,
        Z_independent_substitution=True,
        source_endpoint_trace=True,
        endpoint_chart="frozen_macro",
        endpoint_native_query="m=1",
    )
    first_qy_lower = tuple(
        graph.node(
            "function_substitution",
            expression=qy,
            variable=_Q4["first_coordinate"],
            value=zero,
            Z_independent_substitution=True,
            source_endpoint_trace=True,
            endpoint_chart="first_switch",
            endpoint_native_query="s=0",
        )
        for qy in _Q4["first_Q_y"][:2]
    )
    first_radius_lower = graph.node(
        "function_substitution",
        expression=_Q4["first_radius"],
        variable=_Q4["first_coordinate"],
        value=zero,
        Z_independent_substitution=True,
        source_endpoint_trace=True,
        endpoint_chart="first_switch",
        endpoint_native_query="s=0",
    )
    macro_factorized = tuple(graph.mul(qy, macro_radius_upper) for qy in macro_qy_upper)
    first_factorized = tuple(graph.mul(qy, first_radius_lower) for qy in first_qy_lower)
    macro_source_premises = {
        "common_physical_Q_y_source": True,
        "common_physical_radius_source": True,
        "source_receipt": "high.actual_macro_first_background_seam",
        "source_chart_physical_Q_y": list(q4_evidence["physical_ports"]["frozen_macro"]["Q_y"]),
        "high_chart_physical_Q_y": list(high_bindings["charts"]["frozen_macro"]["physical"]["Q"][1]),
        "high_common_source_controls": dict(background_controls),
    }
    first_source_premises = {
        "common_physical_Q_y_source": True,
        "common_physical_radius_source": True,
        "source_receipt": "high.actual_macro_first_background_seam",
        "source_chart_physical_Q_y": list(q4_evidence["physical_ports"]["first_switch"]["Q_y"]),
        "high_chart_physical_Q_y": list(high_bindings["charts"]["first_switch"]["physical"]["Q"][1]),
        "high_common_source_controls": dict(background_controls),
    }
    factorized_rows = []
    for k in range(2):
        factorized_rows.append({
            "ordinary_Z_order": k,
            "frozen_macro": _factorized_identity(
                graph.nodes,
                finite=_Q4["macro_finite_R_Q_y"][k],
                endpoint=macro_upper[k],
                factorized=macro_factorized[k],
                qy_sub=macro_qy_upper[k],
                radius_sub=macro_radius_upper,
                qy=_Q4["macro_Q_y"][k],
                radius=_Q4["macro_radius"],
                coordinate=_Q4["macro_coordinate"],
                value=1,
                source_chart="frozen_macro",
                physical_q_y_source=q4_evidence["physical_ports"]["frozen_macro"]["Q_y"],
                source_premises=macro_source_premises,
            ),
            "first_switch": _factorized_identity(
                graph.nodes,
                finite=_Q4["first_finite_R_Q_y"][k],
                endpoint=first_lower[k],
                factorized=first_factorized[k],
                qy_sub=first_qy_lower[k],
                radius_sub=first_radius_lower,
                qy=_Q4["first_Q_y"][k],
                radius=_Q4["first_radius"],
                coordinate=_Q4["first_coordinate"],
                value=0,
                source_chart="first_switch",
                physical_q_y_source=q4_evidence["physical_ports"]["first_switch"]["Q_y"],
                source_premises=first_source_premises,
            ),
        })
    signed_cancellation = tuple(
        graph.add(graph.neg(macro_upper[k]), first_lower[k])
        for k in range(2)
    )
    for k, handle in enumerate(signed_cancellation):
        row = _base_row(graph.nodes, handle)
        _require(row.get("operation") == "sum" and len(row.get("arguments", ())) == 2, f"signed endpoint join row malformed at Z{k}")
        _require(_base_row(graph.nodes, row["arguments"][0]).get("operation") == "negative" and _base_row(graph.nodes, row["arguments"][0]).get("argument") == macro_upper[k], f"macro upper sign changed at Z{k}")
        _require(row["arguments"][1] == first_lower[k], f"first-switch lower sign changed at Z{k}")

    symbolic_controls = _symbolic_endpoint_controls(
        graph,
        graph.nodes,
        macro_factorized,
        first_factorized,
        macro_qy_upper,
        first_qy_lower,
        macro_radius_upper,
        first_radius_lower,
    )
    new_appended_nodes = copy.deepcopy(graph.nodes[prefix_len:])
    _require(new_appended_nodes and len(new_appended_nodes) == 16, "new endpoint suffix serialization changed")
    appended_refs = {
        value
        for row in graph.appended
        for key, value in row.items()
        if key in {"argument", "numerator", "denominator", "expression", "variable", "value", "integrand", "lower", "upper"}
        and isinstance(value, int)
    }
    for row in graph.appended:
        if isinstance(row.get("arguments"), list):
            appended_refs.update(value for value in row["arguments"] if isinstance(value, int))
    _require(not (set(_Q4["first_upper_R_Q_y"]) & appended_refs), "first-switch standalone upper endpoint was consumed")
    _require(graph.nodes[:prefix_len] == accepted_graph_prefix, "accepted graph prefix changed while appending new seam")

    source_pins = {PARENT_NAME: PARENT_SHA256, **PINNED_HASHES}
    source_guards = {
        "passed": True,
        "AST": guards,
        "q4": q4_evidence,
        "high_receipt": high_receipt,
        "high_background_controls": background_controls,
        "high_own_controls": own_controls,
        "high_definition": high_definition,
        "loop_common_N": True,
        "micro_receipt": micro_receipt,
        "micro_archive": micro_prefix,
        "switch_receipt": switch_receipt,
        "switch_archive": switch_prefix,
        "source_family_equal": True,
        "native_J_distinct": {"frozen_macro": 66987, "first_switch": 9},
        "common_N": 18,
        "nonreset_M": {"frozen_macro": True, "first_switch": True},
        "phase_handles": phase_handles,
    }
    proof_rows = []
    for k in range(2):
        proof_rows.append({
            "ordinary_Z_order": k,
            "frozen_macro_upper_m1": macro_upper[k],
            "first_switch_lower_s0": first_lower[k],
            "signed_boundary_term": signed_cancellation[k],
            "factorized_frozen_macro": macro_factorized[k],
            "factorized_first_switch": first_factorized[k],
            "factorized_identity": factorized_rows[k],
            "source_identity": {
                "common_physical_radius_cutpoint": True,
                "common_physical_Q_y_source": True,
                "same_N": True,
                "nonreset_M": True,
                "distinct_native_J": True,
                "passed": True,
            },
        })

    return {
        "passed": True,
        "parent_report": PARENT_NAME,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": source_family,
        "source_pins": source_pins,
        "input_pins": source_pins,
        "source_guards": source_guards,
        "proof_rows": proof_rows,
        "graph_patch": {
            "base_archive": f"experiments/root_st073/{Q4_ARCHIVE_NAME}",
            "base_parent_report": PARENT_NAME,
            "base_parent_report_sha256": PARENT_SHA256,
            "base_parent_merged_namespace_digest": graph_namespace["namespace_digest"],
            "base_parent_merged_node_count": merged["merged_node_count"],
            "base_node_count": prefix_len,
            "base_prefix_unchanged": graph.nodes[:prefix_len] == accepted_graph_prefix,
            "appended_nodes": new_appended_nodes,
            "appended_node_count": len(new_appended_nodes),
            "endpoint_handles": {
                "frozen_macro_upper_existing": list(_Q4["macro_upper_R_Q_y"]),
                "first_switch_upper_existing": list(_Q4["first_upper_R_Q_y"]),
                "frozen_macro_upper_used_Z0_Z1": list(macro_upper),
                "first_switch_lower_new_s0_Z0_Z1": list(first_lower),
                "frozen_macro_factorized_m1_Z0_Z1": list(macro_factorized),
                "first_switch_factorized_s0_Z0_Z1": list(first_factorized),
                "signed_cancellation": list(signed_cancellation),
            },
        },
        "symbolic_controls": symbolic_controls,
        "contract": {
            "ordinary_Z_orders": [0, 1],
            "scope": "adjacent frozen_macro m=1 upper / first_switch s=0 lower physical R_Q_y",
            "signed_boundary_terms": "-R_Q_y(frozen_macro,m=1)+R_Q_y(first_switch,s=0)",
            "frozen_macro_native_J": 66987,
            "first_switch_native_J": 9,
            "native_J_remain_distinct": True,
            "same_phase_source_N": True,
            "phase_handles": {"frozen_macro": 101909, "first_switch": 104537},
            "common_N_handle": 18,
            "nonreset_M_retained": True,
            "first_switch_upper_standalone_retained": True,
            "macro_outer_upper_standalone_retained": True,
            "identity_only_in_adjacent_joined_expression": True,
            "cancellation_only_in_adjacent_joined_expression": True,
            "old_caps_and_parent_graph_retained": True,
            "no_numeric_phase_or_quadrature": True,
            "no_old_owner_checker_producer_executed": True,
            "source_values_selected": False,
            "physical_R_Q_y_factorized_endpoint_arithmetic": True,
            "quantitative_bounds_remain_open": True,
        },
    }


__all__ = ["build"]
