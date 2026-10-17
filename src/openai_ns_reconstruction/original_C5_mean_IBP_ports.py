"""Source-bound signed C5 mean-history ports for the outer O2/O3 charts.

This module reads the accepted inlet receipt and the saved true outer C6/C5
graph.  It adds only a small local expression graph: ``C_phi`` is the
accepted phase primitive of the actual signed ``B`` rows, ``C_y`` is the
phase primitive of the actual fixed-phase ``slow_B_y`` rows, and the two
rate-one Duhamel updates are written by integration by parts.

No producer, owner, checker, fixture, quadrature, cap, midpoint, physical
scale, or selected integer is executed or materialized.  The returned H rows
are exact source-bound expressions for ``H = N * delta_m``; they are not
numeric signed H values.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_signed_density_expression import ScalarExpression, SignedDensityExpressionGraph


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
)
PARENT_REPORT_SHA256 = "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e"

OUTER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
OUTER_NAME = Path(OUTER_GRAPH).name

MODIFIED_Q4_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
MODIFIED_Q4_NAME = Path(MODIFIED_Q4_SOURCE).name

TRUE_OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
TRUE_OUTER_NAME = Path(TRUE_OUTER_SOURCE).name

ZERO_MEAN_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_signed_averaging_check.json"
)
ZERO_MEAN_NAME = Path(ZERO_MEAN_RECEIPT).name
ZERO_MEAN_RECEIPT_SHA256 = "3f72fde117295a935d8873578f70617cf18f2225029ed8427d55ef9de0a9f6c2"

GATE = "current_original_C5_signed_mean_IBP_ports_installed"
ORDINARY_Z_ROWS = 4
CHARTS = ("O2_buffer", "O3_transition")

# These are source-archive anchors, not selected field values.  The binding
# tables are checked against them so a later graph cannot silently redirect a
# port to a cap, midpoint, or another chart.
ANCHORS = {
    "O2_buffer": {
        "coordinate": 232895,
        "domain": (0, 4),
        "phase": 233900,
        "radius": 233136,
        "jacobian": 1,
        "B": (235485, 235486, 235487, 235488),
        "slow_B_y": (239919, 239953, 240006, 240080),
        "incoming": (232684, 232685, 232686, 232687),
        "outgoing": (241325, 241326, 241327, 241328),
        "density": (235727, 235728, 235729, 235730),
        "memory": 232976,
        "history0": 235936,
        "history_prefix0": 235930,
        "contribution0": 235924,
        "integral_integrand0": 235918,
        "integral_variable": "true_outer_transport_O2_buffer_5",
        "phase_offset": 232900,
    },
    "O3_transition": {
        "coordinate": 241536,
        "domain": (0, 1),
        "phase": 242551,
        "radius": 241787,
        "jacobian": 1,
        "B": (244136, 244137, 244138, 244139),
        "slow_B_y": (248598, 248632, 248685, 248759),
        "incoming": (241325, 241326, 241327, 241328),
        "outgoing": (250004, 250005, 250006, 250007),
        "density": (244378, 244379, 244380, 244381),
        "memory": 241627,
        "history0": 244587,
        "history_prefix0": 244581,
        "contribution0": 244575,
        "integral_integrand0": 244569,
        "integral_variable": "true_outer_transport_O3_transition_5",
        "phase_offset": 241542,
    },
}


def sha(path: Any) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load(path: Any) -> Any:
    raw = Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo(root: Any) -> Path:
    path = Path(root)
    if path.is_file():
        path = path.parent
    for candidate in (path, path.parent, Path(__file__).resolve().parents[2]):
        if (candidate / "src").is_dir() and (candidate / "experiments").is_dir():
            return candidate.resolve()
    raise ValueError("repository root with src/ and experiments/ is required")


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _source_family(value: Any) -> dict[str, Any]:
    family = _mapping(value)
    return {
        key: family.get(key)
        for key in (
            "actual_five_defect_family_sha256",
            "implicit_source_sha256",
            "datum_enclosure_sha256",
        )
    }


def _load_parent(repo: Path, supplied: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = repo / PARENT_REPORT
    actual_sha = sha(path)
    file_parent = load(path)
    pin = {
        "path": PARENT_REPORT,
        "sha256": actual_sha,
        "expected_sha256": PARENT_REPORT_SHA256,
        "passed": actual_sha == PARENT_REPORT_SHA256,
        "source_executed": False,
    }
    _require(pin["passed"], "Accepted original inlet-source parent hash mismatch")
    if supplied is None:
        parent = file_parent
        supplied_path = None
    elif isinstance(supplied, Mapping):
        parent = supplied
        supplied_path = None
    else:
        supplied_path = str(supplied)
        candidate = Path(supplied)
        if not candidate.is_absolute():
            candidate = repo / candidate
        parent = load(candidate)
    pin.update(
        {
            "supplied_path": supplied_path,
            "accepted_mapping_sha256": digest(file_parent),
            "supplied_mapping_sha256": digest(parent),
            "mapping_match": digest(parent) == digest(file_parent),
        }
    )
    _require(pin["mapping_match"], "Supplied inlet-source parent does not match the accepted report")
    return parent, pin


def _pin(parent: Mapping[str, Any], repo: Path, relative: str, kind: str) -> dict[str, Any]:
    name = Path(relative).name
    expected = _mapping(parent.get("input_hashes")).get(name)
    path = repo / relative
    actual = sha(path)
    return {
        "name": name,
        "kind": kind,
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": isinstance(expected, str) and actual == expected,
        "source_executed": False,
    }


def _ast_guard(repo: Path, relative: str, function_name: str,
               required: tuple[str, ...]) -> dict[str, Any]:
    source = (repo / relative).read_text(encoding="utf-8")
    tree = ast.parse(source)
    if "." in function_name:
        class_name, method_name = function_name.split(".", 1)
        matches = [
            method
            for node in ast.walk(tree)
            if isinstance(node, ast.ClassDef) and node.name == class_name
            for method in node.body
            if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef))
            and method.name == method_name
        ]
    else:
        matches = [
            node for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == function_name
        ]
    _require(len(matches) == 1, "Unique source AST function required: " + function_name)
    function = matches[0]
    segment = ast.get_source_segment(source, function) or ""
    compact = "".join(segment.split())
    checks = {item: "".join(item.split()) in compact for item in required}
    return {
        "passed": all(checks.values()),
        "path": relative,
        "function": function_name,
        "required_text": checks,
        "AST_sha256": hashlib.sha256(
            ast.dump(function, include_attributes=False).encode()
        ).hexdigest(),
        "source_executed": False,
    }


def _exact_one(node: Mapping[str, Any]) -> bool:
    return (
        node.get("operation") == "exact_rational"
        and node.get("numerator") == 1
        and node.get("denominator") == 1
    )


def _exact_value(node: Mapping[str, Any], value: int) -> bool:
    return (
        node.get("operation") == "exact_rational"
        and node.get("numerator") == value
        and node.get("denominator") == 1
    )


def _phase_chain_shape(
    nodes: list[Any],
    expected: Mapping[str, Any],
    bindings: Mapping[str, Any],
) -> dict[str, Any]:
    phase_id = expected["phase"]
    phase = _mapping(nodes[phase_id])
    phase_argument_id = phase.get("argument")
    phase_argument = (
        _mapping(nodes[phase_argument_id])
        if isinstance(phase_argument_id, int)
        else {}
    )
    phase_args = phase_argument.get("arguments", ())
    n_id = bindings.get("N")
    origin_id = bindings.get("phase_origin")
    product_shape = (
        phase.get("operation") == "analytic_unary"
        and phase.get("name") == "fractional_part"
        and phase_argument.get("operation") == "product"
        and isinstance(n_id, int)
        and set(phase_args) == {n_id, next((x for x in phase_args if x != n_id), None)}
    )
    inner_id = next((x for x in phase_args if x != n_id), None)
    inner = _mapping(nodes[inner_id]) if isinstance(inner_id, int) else {}
    log_id = next(
        (
            x
            for x in inner.get("arguments", ())
            if isinstance(x, int)
            and _mapping(nodes[x]).get("operation") == "analytic_unary"
            and _mapping(nodes[x]).get("name") == "log"
        ),
        None,
    )
    negative_id = next(
        (
            x
            for x in inner.get("arguments", ())
            if isinstance(x, int) and _mapping(nodes[x]).get("operation") == "negative"
        ),
        None,
    )
    log_node = _mapping(nodes[log_id]) if isinstance(log_id, int) else {}
    ratio_id = log_node.get("argument")
    ratio = _mapping(nodes[ratio_id]) if isinstance(ratio_id, int) else {}
    negative = _mapping(nodes[negative_id]) if isinstance(negative_id, int) else {}
    checks = {
        "fractional_part": (
            phase.get("operation") == "analytic_unary"
            and phase.get("name") == "fractional_part"
        ),
        "argument_is_product_with_common_N": product_shape,
        "argument_inner_is_log_minus_phase_origin": (
            inner.get("operation") == "sum"
            and len(inner.get("arguments", ())) == 2
            and isinstance(log_id, int)
            and isinstance(negative_id, int)
        ),
        "log_is_actual_radius_over_Ra": (
            log_node.get("operation") == "analytic_unary"
            and log_node.get("name") == "log"
            and ratio.get("operation") == "positive_quotient"
            and ratio.get("numerator") == expected["radius"]
            and ratio.get("denominator") == 65604
            and ratio.get("source_positive_certificate") == "same original positive actual R/Ra"
        ),
        "phase_origin_is_literal_source_node": (
            negative.get("operation") == "negative"
            and negative.get("argument") == origin_id
            and origin_id == 24954
        ),
        "common_N_and_phase_origin_bindings": n_id == 18 and origin_id == 24954,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "phase_node": phase_id,
        "phase_argument_node": phase_argument_id,
        "inner_node": inner_id,
        "log_node": log_id,
        "ratio_node": ratio_id,
        "phase_origin_node": origin_id,
        "source_literal_chain": True,
    }


def _memory_shape(nodes: list[Any], expected: Mapping[str, Any]) -> dict[str, Any]:
    memory = _mapping(nodes[expected["memory"]])
    memory_argument_id = memory.get("argument")
    memory_argument = (
        _mapping(nodes[memory_argument_id])
        if isinstance(memory_argument_id, int)
        else {}
    )
    width_id = next(
        (
            x
            for x in memory_argument.get("arguments", ())
            if isinstance(x, int)
            and _mapping(nodes[x]).get("operation") == "sum"
        ),
        None,
    )
    width = _mapping(nodes[width_id]) if isinstance(width_id, int) else {}
    negative_id = next(
        (
            x
            for x in width.get("arguments", ())
            if isinstance(x, int) and _mapping(nodes[x]).get("operation") == "negative"
        ),
        None,
    )
    negative = _mapping(nodes[negative_id]) if isinstance(negative_id, int) else {}
    substitution_id = negative.get("argument")
    substitution = (
        _mapping(nodes[substitution_id])
        if isinstance(substitution_id, int)
        else {}
    )
    checks = {
        "negative_rate_is_exact_minus_one": _exact_value(nodes[220], -1),
        "memory_is_exp": memory.get("operation") == "analytic_unary" and memory.get("name") == "exp",
        "memory_argument_is_negative_width": (
            memory_argument.get("operation") == "product"
            and set(memory_argument.get("arguments", ()))
            == {220, width_id}
        ),
        "width_is_offset_minus_lower_offset": (
            width.get("operation") == "sum"
            and set(width.get("arguments", ()))
            == {expected["phase_offset"], negative_id}
        ),
        "lower_offset_is_source_function_substitution": (
            negative.get("operation") == "negative"
            and substitution.get("operation") == "function_substitution"
            and substitution.get("expression") == expected["phase_offset"]
            and substitution.get("variable") == expected["coordinate"]
            and substitution.get("value") == expected["domain"][0]
            and substitution.get("Z_independent_substitution") is True
        ),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "memory_node": expected["memory"],
        "width_node": width_id,
        "lower_offset_substitution_node": substitution_id,
        "source_memory_shape": "exp(-(offset(query)-offset(domain_lower)))",
    }


def _source_anchor_checks(nodes: list[Any], bindings: Mapping[str, Any]) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    charts = _mapping(bindings.get("charts"))
    backgrounds = _mapping(bindings.get("background"))
    endpoints = _mapping(bindings.get("endpoints"))
    for chart_name in CHARTS:
        expected = ANCHORS[chart_name]
        chart = _mapping(charts.get(chart_name))
        endpoint = _mapping(endpoints.get(chart_name))
        loop = _mapping(chart.get("loop"))
        radial = _mapping(chart.get("radial"))
        incoming = _mapping(chart.get("incoming"))
        finite = _mapping(chart.get("finite_N"))
        histories = _mapping(chart.get("history"))
        contributions = _mapping(chart.get("contributions"))
        memory = _mapping(chart.get("memory"))
        actual = {
            "coordinate": chart.get("coordinate"),
            "domain": tuple(chart.get("domain", ())),
            "phase": chart.get("actual_phase"),
            "radius": chart.get("radius"),
            "jacobian": chart.get("physical_measure"),
            "B": tuple(_mapping(loop.get("B_over_Pstar")).values())
            if isinstance(loop.get("B_over_Pstar"), Mapping)
            else tuple(loop.get("B_over_Pstar", ())[:4]),
            "slow_B_y": tuple(radial.get("slow_B_y", ())[:4]),
            "incoming": tuple(incoming.get("m", ())[:4]),
            "outgoing": tuple(_mapping(endpoint.get("correction")).get("m", ())[:4]),
            "density": tuple(_mapping(finite.get("five_signed_density_increments")).get("m", ())[:4]),
            "memory": memory.get("m"),
            "phase_offset": _mapping(backgrounds.get(chart_name)).get("log_offset"),
            "history0": histories.get("m", ())[0] if histories.get("m") else None,
            "contribution0": contributions.get("m", ())[0] if contributions.get("m") else None,
        }
        # The source binding stores B as a list, while the fallback above is
        # only defensive for malformed archives.  Keep the exact four rows.
        if isinstance(loop.get("B_over_Pstar"), list):
            actual["B"] = tuple(loop["B_over_Pstar"][:ORDINARY_Z_ROWS])
        checks[chart_name] = {
            "literal_anchor_ids": {
                key: list(value) if isinstance(value, tuple) else value
                for key, value in expected.items()
            },
            "binding_values": {
                key: list(value) if isinstance(value, tuple) else value
                for key, value in actual.items()
            },
            "binding_matches_literal_anchors": all(
                actual.get(key) == value for key, value in expected.items()
                if key in actual
            ),
        }
        checks[chart_name]["node_shape_checks"] = {
            "coordinate_bound": (
                isinstance(expected["coordinate"], int)
                and nodes[expected["coordinate"]].get("operation") == "bound_variable"
                and nodes[expected["coordinate"]].get("name") == "true_outer_" + chart_name
            ),
            "phase_fractional_part": (
                nodes[expected["phase"]].get("operation") == "analytic_unary"
                and nodes[expected["phase"]].get("name") == "fractional_part"
            ),
            "jacobian_exact_one": _exact_one(nodes[expected["jacobian"]]),
            "memory_exp": (
                nodes[expected["memory"]].get("operation") == "analytic_unary"
                and nodes[expected["memory"]].get("name") == "exp"
            ),
        }
        checks[chart_name]["domain_endpoint_checks"] = {
            "domain_lower_is_literal_zero": (
                expected["domain"][0] == 0 and _exact_value(nodes[0], 0)
            ),
            "domain_upper_is_exact_source_value": (
                expected["domain"][1] == (4 if chart_name == "O2_buffer" else 1)
                and _exact_value(
                    nodes[expected["domain"][1]],
                    11 if chart_name == "O2_buffer" else 1,
                )
            ),
            "endpoint_native_node_is_domain_upper": expected["domain"][1]
            in (4, 1),
        }
        checks[chart_name]["phase_chain"] = _phase_chain_shape(
            nodes, expected, bindings
        )
        checks[chart_name]["memory_shape"] = _memory_shape(nodes, expected)
        checks[chart_name]["density_rows_are_B_over_N"] = all(
            nodes[density].get("operation") == "product"
            and set(nodes[density].get("arguments", ())) == {b_row, 99411}
            for density, b_row in zip(expected["density"], expected["B"])
        )
        checks[chart_name]["history_and_integral_anchors"] = {
            "history0_is_prefix_plus_integral": (
                nodes[expected["history0"]].get("operation") == "sum"
                and set(nodes[expected["history0"]].get("arguments", ()))
                == {expected["history_prefix0"], expected["contribution0"]}
            ),
            "integral_shape": (
                nodes[expected["contribution0"]].get("operation") == "definite_integral"
                and nodes[expected["contribution0"]].get("integrand") == expected["integral_integrand0"]
                and nodes[expected["contribution0"]].get("variable")
                == expected["integral_variable"]
                and nodes[expected["contribution0"]].get("lower") == 0
                and nodes[expected["contribution0"]].get("upper") == expected["coordinate"]
                and nodes[expected["contribution0"]].get("physical_log_offset")
                == expected["phase_offset"]
                and nodes[expected["contribution0"]].get("physical_Jacobian_applied_once") is True
                and nodes[expected["contribution0"]].get("original_rate") == "1"
            ),
        }
        checks[chart_name]["passed"] = (
            checks[chart_name]["binding_matches_literal_anchors"]
            and all(checks[chart_name]["node_shape_checks"].values())
            and all(checks[chart_name]["domain_endpoint_checks"].values())
            and checks[chart_name]["phase_chain"]["passed"]
            and checks[chart_name]["memory_shape"]["passed"]
            and checks[chart_name]["density_rows_are_B_over_N"]
            and all(checks[chart_name]["history_and_integral_anchors"].values())
        )
    checks["O3_incoming_is_O2_outgoing"] = (
        tuple(_mapping(charts.get("O3_transition")).get("incoming", {}).get("m", ())[:4])
        == tuple(_mapping(endpoints.get("O2_buffer")).get("correction", {}).get("m", ())[:4])
        == ANCHORS["O3_transition"]["incoming"]
    )
    checks["passed"] = all(
        value.get("passed", value) if isinstance(value, Mapping) else value
        for key, value in checks.items()
        if key != "passed"
    )
    return checks


def _phase_primitive(
    graph: SignedDensityExpressionGraph,
    body: ScalarExpression,
    phase: ScalarExpression,
    coordinate: ScalarExpression,
    selected_N: ScalarExpression,
    chart: str,
    order: int,
    *,
    slow_y: bool = False,
) -> ScalarExpression:
    operation = (
        "current_original_C5_phase_integral_of_slow_B_y"
        if slow_y
        else "current_original_n1_periodic_phase_primitive"
    )
    metadata = {
        "operation": operation,
        "body": body.node,
        "integrand": "actual slow_B_y row %d" % order if slow_y else "actual B row %d" % order,
        "native_variable": "true_outer_" + chart,
        "native_coordinate_node": coordinate.node,
        "phase_variable": "original_global_phase",
        "source_phase_node": phase.node,
        "derived_phase_not_bound_variable": True,
        "period": "1",
        "phase_integral_lower": 0,
        "phase_endpoint_zero_mean": True,
        "primitive_endpoint_identity": "C(0)=C(1)=0 from accepted B parity",
        "selected_N": selected_N.node,
        "same_original_P0_and_global_phase": True,
        "source_zero_mean_receipt": ZERO_MEAN_NAME,
        "ordinary_Z_order": order,
        "signed_function_not_a_range_endpoint": True,
        "phase_integral_definition": (
            "C_y(y,Z,phi)=integral_0^phi slow_B_y(y,Z,s)ds"
            if slow_y
            else "C_phi(y,Z,phi)=integral_0^phi B(y,Z,s)ds"
        ),
        "original_graph_nodes_added_or_replaced": False,
    }
    if slow_y:
        metadata.update(
            {
                "body_role": "actual_slow_B_y_integrand",
                "fixed_phase": True,
                "total_y_chain_excluded": True,
                "phase_primitive_commutes_with_slow_y": True,
                "qualitative_C5_source_identity_only": True,
                "genuine_mixed_slow_y_bound_available": False,
                "not_a_derivative_of_slow_B_y": True,
            }
        )
    return graph.insert(metadata)


def _native_substitution(
    graph: SignedDensityExpressionGraph,
    expression: ScalarExpression,
    coordinate: ScalarExpression,
    value: ScalarExpression,
    phase: ScalarExpression,
    name: str,
    *,
    outer_integration_variable: bool = False,
) -> ScalarExpression:
    return graph.insert(
        {
            "operation": "function_substitution",
            "expression": expression.node,
            "variable": coordinate.node,
            "value": value.node,
            "name": name,
            "source_phase_node": phase.node,
            "derived_phase_recomputed_at_value": True,
            "derived_phase_not_bound_variable": True,
            "capture_avoiding": True,
            "preserve_distinct_outer_integration_variable": True,
            "outer_integration_variable": outer_integration_variable,
            "original_graph_nodes_added_or_replaced": False,
        }
    )


def _integral(
    graph: SignedDensityExpressionGraph,
    chart: str,
    coordinate: ScalarExpression,
    phase: ScalarExpression,
    lower: ScalarExpression,
    query: ScalarExpression,
    jacobian: ScalarExpression,
    primitive: ScalarExpression,
    slow_primitive: ScalarExpression,
    order: int,
) -> tuple[ScalarExpression, ScalarExpression, ScalarExpression]:
    outer = graph.insert(
        {
            "operation": "bound_variable",
            "name": "original_C5_mean_IBP_outer_" + chart,
            "role": "distinct physical dnative integration coordinate",
            "original_graph_nodes_added_or_replaced": False,
        }
    )
    C_outer = _native_substitution(
        graph, primitive, coordinate, outer, phase,
        "C_phi_at_outer_" + chart, outer_integration_variable=True,
    )
    Cy_outer = _native_substitution(
        graph, slow_primitive, coordinate, outer, phase,
        "C_y_slow_at_outer_" + chart, outer_integration_variable=True,
    )
    kernel = graph.unary("exp", -(query - outer))
    integrand = kernel * (C_outer + Cy_outer) * jacobian
    integral = graph.insert(
        {
            "operation": "definite_integral",
            "integrand": integrand.node,
            "variable": outer.node,
            "lower": lower.node,
            "upper": query.node,
            "variable_name": "original_C5_mean_IBP_outer_" + chart,
            "upper_variable_name": "true_outer_" + chart,
            "measure": "physical dy = Jacobian * dnative exactly once",
            "original_rate": 1,
            "rate_one_kernel": True,
            "physical_Jacobian_applied_once": True,
            "jacobian_source_node": jacobian.node,
            "ordinary_slow_Z_derivative_order": order,
            "free_query_endpoint": query.node,
            "bound_inner_native_variable": coordinate.node,
            "bound_phase_variable": "original_global_phase",
            "source_phase_node": phase.node,
            "derived_phase_not_bound_variable": True,
            "endpoint_terms_retained": True,
            "inner_and_outer_native_variables_distinct": True,
            "source_B_and_slow_B_y_are_signed_functions": True,
            "numerical_value_not_installed": True,
            "original_graph_nodes_added_or_replaced": False,
        }
    )
    return outer, integrand, integral


def scaled_history_update(
    H_in: Any,
    memory: Any,
    C_out: Any,
    C_in: Any,
    integral: Any,
    inverse_N: Any,
) -> Any:
    """Rate-one source-bound identity for ``H = N * delta_m``.

    The function intentionally contains only arithmetic.  It works for
    ordinary numeric arguments in small independent controls and for
    ``ScalarExpression`` arguments in the local source graph.
    """
    return memory * H_in + inverse_N * (C_out - memory * C_in - integral)


def scaled_history_update_control() -> dict[str, Any]:
    """Check the transformed transport equation, independently of its formula.

    C denotes the issued phase primitive along the physical path. Thus
    d_y C = C_y_slow + N*B, while J is the rate-one slow convolution.
    The control uses arbitrary functions and a nonzero symbolic inlet.
    """
    import sympy as s
    y = s.symbols("y", real=True)
    N, inlet = s.symbols("N inlet", nonzero=True)
    C_in = s.symbols("C_in")
    C, slow, B, J = (s.Function(name)(y) for name in ("C", "slow", "B", "J"))
    H = scaled_history_update(inlet, s.exp(-y), C, C_in, J, 1/N)
    derivative = (s.diff(H, y) + H).subs({
        s.diff(C, y): slow + N*B,
        s.diff(J, y): C + slow - J,
    }, simultaneous=True)
    at_inlet = H.subs(y, 0).subs({s.Function("C")(0): C_in, s.Function("J")(0): 0})
    checks = {
        "Nscaled_transport_derivative_plus_H_equals_actual_B": s.simplify(derivative-B) == 0,
        "zero_length_update_recovers_arbitrary_nonzero_inlet": s.simplify(at_inlet-inlet) == 0,
        "incoming_response_is_exact_memory": s.simplify(s.diff(H, inlet)-s.exp(-y)) == 0,
    }
    return {
        "passed": all(checks.values()), "checks": checks,
        "original_mean_equation": "delta_m_y+delta_m=B/N",
        "scaled_mean_equation": "H_y+H=B; H=N*delta_m",
        "phase_chain": "d_y C=C_y_slow+N*B",
        "slow_convolution_equation": "J_y+J=C+C_y_slow; J(0)=0",
        "nonzero_incoming_retained": True,
        "source_chart_H_values_not_evaluated": True,
        "symbolic_control_only": True,
    }


def _chart_update(
    graph: SignedDensityExpressionGraph,
    nodes: list[Any],
    bindings: Mapping[str, Any],
    chart: str,
    incoming_H: list[ScalarExpression],
    N: ScalarExpression,
    inverse_N: ScalarExpression,
    source_sha: str,
) -> tuple[list[ScalarExpression], list[ScalarExpression], dict[str, Any]]:
    archive = ANCHORS[chart]
    chart_binding = _mapping(_mapping(bindings.get("charts")).get(chart))
    coordinate_id = archive["coordinate"]
    phase_id = archive["phase"]
    query = graph.reference(OUTER_NAME, source_sha, coordinate_id, source_nodes=nodes)
    phase = graph.reference(OUTER_NAME, source_sha, phase_id, source_nodes=nodes)
    lower = graph.reference(OUTER_NAME, source_sha, 0, source_nodes=nodes)
    endpoint = graph.reference(
        OUTER_NAME,
        source_sha,
        archive["domain"][1],
        source_nodes=nodes,
    )
    jacobian = graph.reference(OUTER_NAME, source_sha, archive["jacobian"], source_nodes=nodes)
    memory = graph.reference(OUTER_NAME, source_sha, archive["memory"], source_nodes=nodes)
    B = [
        graph.reference(OUTER_NAME, source_sha, value, source_nodes=nodes)
        for value in archive["B"]
    ]
    slow_B_y = [
        graph.reference(OUTER_NAME, source_sha, value, source_nodes=nodes)
        for value in archive["slow_B_y"]
    ]
    primitive = [
        _phase_primitive(graph, body, phase, query, N, chart, order)
        for order, body in enumerate(B)
    ]
    primitive_slow_y = [
        _phase_primitive(graph, body, phase, query, N, chart, order, slow_y=True)
        for order, body in enumerate(slow_B_y)
    ]
    derivative_identities = [
        graph.insert(
            {
                "operation": "phase_primitive_derivative_identity",
                "primitive": C.node,
                "phase_derivative": body.node,
                "identity": "d_phi C_phi = B",
                "source_C5_periodic_primitive": True,
                "ordinary_Z_order": order,
                "source_zero_mean_receipt": ZERO_MEAN_NAME,
                "original_graph_nodes_added_or_replaced": False,
            }
        )
        for order, (C, body) in enumerate(zip(primitive, B))
    ]
    derivative_slow_identities = [
        graph.insert(
            {
                "operation": "phase_primitive_slow_y_phi_identity",
                "primitive": C.node,
                "phase_derivative": body.node,
                "identity": "d_phi C_y_slow = slow_B_y at fixed phase; C_y is not d_y(slow_B_y)",
                "phase_primitive_commutes_with_slow_y": True,
                "qualitative_source_identity_only": True,
                "not_a_derivative_of_slow_B_y": True,
                "ordinary_Z_order": order,
                "original_graph_nodes_added_or_replaced": False,
            }
        )
        for order, (C, body) in enumerate(zip(primitive_slow_y, slow_B_y))
    ]
    rows: list[dict[str, Any]] = []
    output: list[ScalarExpression] = []
    for order, (H, C, Cy) in enumerate(zip(incoming_H, primitive, primitive_slow_y)):
        C_in = _native_substitution(
            graph, C, query, lower, phase, "C_phi_in_" + chart + "_Z" + str(order)
        )
        C_out = _native_substitution(
            graph, C, query, query, phase, "C_phi_out_" + chart + "_Z" + str(order)
        )
        outer, integrand, integral = _integral(
            graph, chart, query, phase, lower, query, jacobian, C, Cy, order
        )
        memory_H_in = memory * H
        bracket = C_out - memory * C_in - integral
        H_out = scaled_history_update(H, memory, C_out, C_in, integral, inverse_N)
        output.append(H_out)
        rows.append(
            {
                "ordinary_Z_order": order,
                "incoming_H": H.node,
                "C_phi": C.node,
                "C_y_slow": Cy.node,
                "C_phi_derivative_identity": derivative_identities[order].node,
                "C_y_slow_derivative_identity": derivative_slow_identities[order].node,
                "C_in": C_in.node,
                "C_out": C_out.node,
                "memory": memory.node,
                "memory_H_in": memory_H_in.node,
                "IBP_integral": integral.node,
                "IBP_bracket": bracket.node,
                "H_out": H_out.node,
                "endpoint_terms_retained": True,
                "incoming_history_not_reset": True,
                "rate": 1,
                "physical_Jacobian_factor_count": 1,
                "selected_N_division_count": 1,
            }
        )
    endpoint_output = [
        _native_substitution(
            graph,
            value,
            query,
            endpoint,
            phase,
            "H_out_at_endpoint_" + chart + "_Z" + str(order),
        )
        for order, value in enumerate(output)
    ]
    return output, endpoint_output, {
        "chart": chart,
        "coordinate": query.node,
        "domain_lower_node": lower.node,
        "domain_upper_node": endpoint.node,
        "domain_upper_value": 11 if chart == "O2_buffer" else 1,
        "phase": phase.node,
        "domain": list(archive["domain"]),
        "radius": archive["radius"],
        "jacobian": jacobian.node,
        "N": N.node,
        "inverse_N": inverse_N.node,
        "source_B_rows": list(archive["B"]),
        "source_slow_B_y_rows": list(archive["slow_B_y"]),
        "source_incoming_m_rows": list(archive["incoming"]),
        "source_outgoing_m_rows": list(archive["outgoing"]),
        "source_history_rows": list(_mapping(chart_binding.get("history")).get("m", ())[:4]),
        "source_integral_anchor": archive["contribution0"],
        "rows": rows,
        "H_out_rows": [value.node for value in output],
        "H_out_moving_query_rows": [value.node for value in output],
        "H_out_endpoint_rows": [value.node for value in endpoint_output],
        "endpoint_substitution_is_capture_avoiding": True,
        "endpoint_is_separate_from_moving_query": True,
        "physical_measure": "dy = Jacobian * dnative exactly once",
        "phase_semantics": {
            "source_phase_node": phase.node,
            "derived_phase_not_bound_variable": True,
            "phase_fractional_part_is_from_original_graph": True,
            "phase_primitive_gauge": "integral_0^phi; C(0)=C(1)=0 by accepted parity",
        },
        "source_slow_y_semantics": {
            "fixed_phase": True,
            "total_y_chain_excluded": True,
            "phase_primitive_commutes_with_slow_y": True,
            "quantitative_mixed_slow_bound": False,
        },
        "source_bound_expression_not_numeric": True,
    }


def _receipt_controls(receipt: Mapping[str, Any], family: Mapping[str, Any]) -> dict[str, Any]:
    theorem = _mapping(receipt.get("exact_source_averaging_theorem"))
    parity = _mapping(theorem.get("original_primitive_parity_theorem"))
    identities = _mapping(parity.get("exact_identities"))
    checks = {
        "receipt_all_passed": receipt.get("all_passed") is True,
        "receipt_source_family": _source_family(receipt.get("source_family")) == dict(family),
        "theorem_passed": theorem.get("passed") is True,
        "B_reflection_zero_mean": identities.get("B_reflection_zero_mean") is True,
        "periodic_B_derivative": identities.get("periodic_B_derivative") is True,
        "uniform_fractional_phi_mean": theorem.get(
            "uniform_fractional_phi_mean_of_leading_and_fixed_phase_slow_jets_exact_zero"
        ) is True,
        "primitive_endpoint_identity": theorem.get("primitive")
        == "G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s)ds; G_j(0)=G_j(1)=0",
        "source_IBP_identity": theorem.get("source_IBP")
        == "integral_a^b K*f1(y,Z,frac(N*y+phi0))/N dy = ([K*G]_a^b-integral_a^b K*(G_y_slow+lambda*G)dy)/N^2",
        "slow_y_not_total_y": theorem.get("slow_y_not_total_y") is True,
        "endpoint_terms_retained": theorem.get("each_cell_endpoint_terms_retained") is True,
        "incoming_rate0_memory_retained": theorem.get("incoming_pressure_rate0_memory_retained") is True,
    }
    return {"checks": checks, "passed": all(checks.values()), "source_executed": False}


def build(root: Any, parent: Any = None) -> dict[str, Any]:
    """Build exact local C5 mean-history expressions without numeric evaluation."""
    repo = _repo(root)
    parent_map, parent_pin = _load_parent(repo, parent)
    family = _source_family(parent_map.get("source_family"))
    outer_pin = _pin(parent_map, repo, OUTER_GRAPH, "input")
    q4_pin = _pin(parent_map, repo, MODIFIED_Q4_SOURCE, "input")
    true_outer_source_pin = _pin(parent_map, repo, TRUE_OUTER_SOURCE, "source")
    _require(outer_pin["passed"], "Pinned true outer C6/C5 source changed")
    _require(q4_pin["passed"], "Pinned modified Q4 source changed")
    _require(true_outer_source_pin["passed"], "Pinned true outer C6/C5 source text changed")
    outer = load(repo / OUTER_GRAPH)
    outer_family = _source_family(outer.get("source_family"))
    _require(outer_family == family, "Outer C5 source family differs from accepted parent")
    bindings = _mapping(outer.get("true_outer_C6_C5_bindings"))
    nodes = list(outer.get("exact_graph_nodes", ()))
    _require(nodes, "Saved true outer C5 graph is empty")
    zero_receipt = load(repo / ZERO_MEAN_RECEIPT)
    zero_receipt_pin = {
        "path": ZERO_MEAN_RECEIPT,
        "sha256": sha(repo / ZERO_MEAN_RECEIPT),
        "expected_sha256": ZERO_MEAN_RECEIPT_SHA256,
        "passed": sha(repo / ZERO_MEAN_RECEIPT) == ZERO_MEAN_RECEIPT_SHA256,
        "source_executed": False,
    }
    _require(zero_receipt_pin["passed"], "Accepted native signed-averaging receipt changed")
    receipt_controls = _receipt_controls(zero_receipt, family)
    _require(receipt_controls["passed"], "Accepted signed phase-primitive receipt is not admissible")
    ast_guards = {
        "NativeDerivative.__init__": _ast_guard(
            repo,
            MODIFIED_Q4_SOURCE,
            "NativeDerivative.__init__",
            (
                "self.cache={phase.node:g.zero};self.core_extensions=[]",
            ),
        ),
        "total_radial_functions": _ast_guard(
            repo,
            MODIFIED_Q4_SOURCE,
            "total_radial_functions",
            (
                "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
                "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))forrefinrows(chart['loop'][key])))forkeyin('A','B_over_Pstar')}",
                "slow_B_y=slow['B_over_Pstar']",
                "physical_y_equals_native_partial_divided_by_Jacobian=True",
            ),
        ),
        "transport_geometry": _ast_guard(
            repo,
            TRUE_OUTER_SOURCE,
            "transport_geometry",
            (
                "Duhamel with actual physical log coordinate, once in each measure.",
                "for key,rate in RATES.items():",
                "body=a.scale(a.at(density[key],x,t),g.mul(measure,kernel))",
                "original_outer_Duhamel=True",
                "original_rate=str(rate)",
                "physical_Jacobian_applied_once=True",
                "predecessor_memory_never_reset=True",
            ),
        ),
    }
    _require(all(row["passed"] for row in ast_guards.values()), "C5 radial source AST guard failed")
    anchors = _source_anchor_checks(nodes, bindings)
    _require(anchors["passed"], "True outer C5 source anchors changed")
    selected_N = bindings.get("N")
    _require(selected_N == 18, "Expected common selected N source node 18")
    _require(
        nodes[selected_N].get("operation") == "exact_positive_integer_expression"
        and nodes[selected_N].get("one_common_N_for_all_charts_and_control_denominators") is True,
        "Common exact N source node changed",
    )
    _require(nodes[99411].get("operation") == "positive_quotient", "Exact inverse-N source node changed")
    _require(
        nodes[99411].get("numerator") == 1
        and nodes[99411].get("denominator") == selected_N
        and nodes[99411].get("source_positive_certificate")
        == "same exact positive integer N on every source chart",
        "Exact inverse-N quotient is not the common source quotient",
    )
    graph = SignedDensityExpressionGraph()
    source_sha = outer_pin["sha256"]
    ref = lambda node: graph.reference(
        OUTER_NAME, source_sha, node, source_nodes=nodes
    )
    N = ref(selected_N)
    inverse_N = ref(99411)
    initial_H = [N * ref(node) for node in ANCHORS["O2_buffer"]["incoming"]]
    O2_H, O2_endpoint_H, O2 = _chart_update(
        graph,
        nodes,
        bindings,
        "O2_buffer",
        initial_H,
        N,
        inverse_N,
        source_sha,
    )
    # The O3 inlet is the complete O2 moving-history expression evaluated at
    # the actual O2 upper native endpoint (node 4, exact value 11).  Keeping
    # this as an explicit capture-avoiding substitution prevents an O2 query
    # variable from leaking into the next chart.
    O3_H, O3_endpoint_H, O3 = _chart_update(
        graph,
        nodes,
        bindings,
        "O3_transition",
        O2_endpoint_H,
        N,
        inverse_N,
        source_sha,
    )
    scaling_control = scaled_history_update_control()
    _require(scaling_control["passed"], "Bounded scaled-history identity control failed")
    graph_packet = graph.packet()
    result = {
        "passed": True,
        GATE: True,
        "source_binding_passed": True,
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_report_pin": parent_pin,
        "source_family": family,
        "source_pins": {
            "parent_report": parent_pin,
            "outer_graph": outer_pin,
            "modified_Q4_source": q4_pin,
            "true_outer_source": true_outer_source_pin,
            "zero_mean_receipt": zero_receipt_pin,
        },
        "source_pins_ok": True,
        "AST": ast_guards,
        "accepted_zero_mean_receipt": receipt_controls,
        "source_anchor_checks": anchors,
        "selected_N": {
            "source_node": N.node,
            "inverse_source_node": inverse_N.node,
            "source_integer_not_materialized": True,
            "same_common_N_for_both_charts": True,
        },
        "source_ports": {
            "O2_buffer": O2,
            "O3_transition": O3,
            "O2_initial_H_rows": [value.node for value in initial_H],
            "O2_buffer_moving_H_rows": list(O2["H_out_moving_query_rows"]),
            "O2_buffer_endpoint_H_rows": list(O2["H_out_endpoint_rows"]),
            "O3_transition_incoming_H_rows": [value.node for value in O2_endpoint_H],
            "O3_transition_moving_H_rows": list(O3["H_out_moving_query_rows"]),
            "O3_transition_endpoint_H_rows": list(O3["H_out_endpoint_rows"]),
            "O3_incoming_H_rows_are_O2_endpoint_substitution": True,
        },
        "H": {
            "definition": "H=N*delta_m",
            "initial": "N*original incoming m; nonzero predecessor retained",
            "O2_buffer_moving_query_Z0_to_Z3": list(O2["H_out_moving_query_rows"]),
            "O2_buffer_endpoint_Z0_to_Z3": list(O2["H_out_endpoint_rows"]),
            "O2_buffer_outgoing_Z0_to_Z3": list(O2["H_out_endpoint_rows"]),
            "O3_transition_incoming_endpoint_Z0_to_Z3": [value.node for value in O2_endpoint_H],
            "O3_transition_moving_query_Z0_to_Z3": list(O3["H_out_moving_query_rows"]),
            "O3_transition_endpoint_Z0_to_Z3": list(O3["H_out_endpoint_rows"]),
            "O3_transition_incoming_Z0_to_Z3": [value.node for value in O2_endpoint_H],
            "O3_transition_outgoing_Z0_to_Z3": list(O3["H_out_endpoint_rows"]),
            "source_archive_outgoing_anchor_Z0_to_Z3": list(ANCHORS["O3_transition"]["outgoing"]),
            "O2_endpoint_native_node": 4,
            "O2_endpoint_exact_value": 11,
            "O3_endpoint_native_node": 1,
            "O3_endpoint_exact_value": 1,
            "endpoint_substitution_capture_avoiding": True,
            "source_bound_exact_expression": True,
            "numeric_signed_H_values_evaluated": False,
        },
        "scaled_history_update_control": scaling_control,
        "identities": {
            "C_phi_equals_B": True,
            "C_y_slow_equals_phase_integral_of_actual_slow_B_y": True,
            "C_phase_gauge": "C(0)=C(1)=0 from accepted B parity",
            "rate": 1,
            "physical_dy_jacobian_once": True,
            "separate_N_scaling": "H=N*delta_m; one inverse-N factor per IBP update",
            "all_endpoints_retained": True,
            "phase_is_derived_source_function_not_bound_variable": True,
            "slow_y_is_fixed_phase_not_total_y": True,
            "no_seam_cancellation_claim": True,
        },
        "numeric_value_available": False,
        "signed_numeric_H_values_installed": False,
        "cap_contraction_or_improvement": False,
        "physical_scale_materialized": False,
        "archived_owner_builder_checker_fixture_or_quadrature_executed": False,
        "full_n1_or_terminal_closure": False,
        "original_sources_numerically_executed": False,
        "exact_local_expression_nodes": graph_packet["exact_local_expression_nodes"],
        "local_nodes": graph_packet["local_nodes"],
        "graph_packet": graph_packet,
        "root": str(repo),
        "directory": str(repo / "experiments" / "root_st073"),
        "gaps": [
            "No signed numeric H values or caps are selected; root must keep this source-bound expression receipt.",
            "No mixed slow-y numerical bound is claimed.",
            "No full n1 or terminal closure is installed.",
        ],
    }
    return result


def build_original_C5_mean_IBP_ports(root: Any, parent: Any = None) -> dict[str, Any]:
    return build(root, parent)


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "OUTER_GRAPH",
    "MODIFIED_Q4_SOURCE",
    "TRUE_OUTER_SOURCE",
    "ZERO_MEAN_RECEIPT",
    "ANCHORS",
    "scaled_history_update",
    "scaled_history_update_control",
    "build",
    "build_original_C5_mean_IBP_ports",
]
