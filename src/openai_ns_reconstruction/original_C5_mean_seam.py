"""Source-only C5 seam identity for the actual O2/O3 outer charts.

The saved true outer graph contains separate graph nodes for the O2 buffer
and O3 transition loop expressions.  Their angle variables are bound with
different names, so node identity is not a proof that the functions agree at
the physical seam.  This module checks the source constructors and the
actual endpoint bindings instead.  The resulting receipt proves the
constructor identity for an arbitrary issued phase and ordinary Z orders
zero through three, and consequently the anchored phase primitive identity

    C_buffer(11, Z, phi) = integral_0^phi B_buffer(11, Z, theta) dtheta
                         = C_transition(0, Z, phi).

No source module is imported or executed.  No C, H, cap, quadrature, or
contraction value is evaluated here.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_radial_IBP.json"
)
PARENT_REPORT_SHA256 = (
    "67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0"
)

OUTER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
LOOP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
LOOP_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
GENERIC_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py"
)
GENERIC_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.json"
)
GENERIC_VIEWS = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources_views.json.gz"
)
CORE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
)
CORE_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.json.gz"
)
CONTROL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py"
)
CONTROL_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.json.gz"
)
OUTER_RC_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"
)
OUTER_RC_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.json.gz"
)
Q4_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
ZERO_MEAN_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_signed_averaging_check.json"
)

GATE = "current_original_C5_mean_seam_identity_installed"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
RATES = ("m", "h", "k", "e", "p")
CHARTS = ("O2_buffer", "O3_transition")

# These are archive node anchors.  They identify the actual source ports;
# they are never interpreted as numeric field values.
ANCHORS: dict[str, dict[str, Any]] = {
    "O2_buffer": {
        "coordinate": 232895,
        "domain": (0, 4),
        "upper_value": 11,
        "phase": 233900,
        "radius": 233136,
        "log_offset": 232900,
        "physical_measure": 1,
        "B": (235485, 235486, 235487, 235488),
        "slow_B_y": (239919, 239953, 240006, 240080),
        "incoming_m": (232684, 232685, 232686, 232687),
        "outgoing_m": (241325, 241326, 241327, 241328),
        "angle": 234382,
        "endpoint_radius": 233745,
        "background_E": (232905, 232906, 232907, 232908),
        "background_E_y": (232912, 232913, 232914, 232915),
    },
    "O3_transition": {
        "coordinate": 241536,
        "domain": (0, 1),
        "upper_value": 1,
        "phase": 242551,
        "radius": 241787,
        "log_offset": 241542,
        "physical_measure": 1,
        "B": (244136, 244137, 244138, 244139),
        "slow_B_y": (248598, 248632, 248685, 248759),
        "incoming_m": (241325, 241326, 241327, 241328),
        "outgoing_m": (250004, 250005, 250006, 250007),
        "angle": 243033,
        "endpoint_radius": 242396,
        "background_E": (241552, 241553, 241554, 241555),
        "background_E_y": (241563, 241564, 241565, 241566),
    },
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _sha(path: Any) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load(path: Any) -> Any:
    raw = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _repo(root: Any) -> Path:
    path = Path(root)
    if path.is_file():
        path = path.parent
    candidates = (path, path.parent, Path(__file__).resolve().parents[2])
    for candidate in candidates:
        if (candidate / "src").is_dir() and (candidate / "experiments").is_dir():
            return candidate.resolve()
    raise ValueError("repository root with src/ and experiments/ is required")


def _family(value: Any) -> dict[str, Any]:
    source = _mapping(value)
    return {
        key: source.get(key)
        for key in (
            "actual_five_defect_family_sha256",
            "implicit_source_sha256",
            "datum_enclosure_sha256",
        )
    }


def _exact(node: Any, numerator: int, denominator: int = 1) -> bool:
    item = _mapping(node)
    return (
        item.get("operation") == "exact_rational"
        and item.get("numerator") == numerator
        and item.get("denominator") == denominator
    )


def _node(nodes: Sequence[Any], node_id: Any) -> Mapping[str, Any]:
    if isinstance(node_id, int) and 0 <= node_id < len(nodes):
        return _mapping(nodes[node_id])
    return {}


def _load_parent(repo: Path, supplied: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = repo / PARENT_REPORT
    actual = _sha(path)
    saved = _load(path)
    pin: dict[str, Any] = {
        "path": PARENT_REPORT,
        "sha256": actual,
        "expected_sha256": PARENT_REPORT_SHA256,
        "passed": actual == PARENT_REPORT_SHA256,
        "source_executed": False,
    }
    if supplied is None:
        selected = saved
        supplied_path = None
    elif isinstance(supplied, Mapping):
        selected = supplied
        supplied_path = None
    else:
        supplied_path = str(supplied)
        candidate = Path(supplied)
        if not candidate.is_absolute():
            candidate = repo / candidate
        selected = _load(candidate)
    pin.update(
        {
            "supplied_path": supplied_path,
            "accepted_mapping_sha256": _digest(saved),
            "supplied_mapping_sha256": _digest(selected),
            "mapping_match": _digest(saved) == _digest(selected),
        }
    )
    return selected, pin


def _pin(
    parent: Mapping[str, Any],
    repo: Path,
    relative: str,
    kind: str,
    *,
    expected_override: str | None = None,
) -> dict[str, Any]:
    path = repo / relative
    expected = expected_override
    if expected is None:
        expected = _mapping(parent.get("input_hashes")).get(Path(relative).name)
    try:
        actual: str | None = _sha(path)
    except (FileNotFoundError, OSError):
        actual = None
    return {
        "path": relative,
        "name": Path(relative).name,
        "kind": kind,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": isinstance(expected, str) and actual == expected,
        "source_executed": False,
    }


def _ast_guard(
    repo: Path,
    relative: str,
    symbol: str,
    required: Sequence[str],
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": relative,
        "function": symbol,
        "required_text": {item: False for item in required},
        "AST_sha256": None,
        "passed": False,
        "source_executed": False,
    }
    try:
        source = (repo / relative).read_text(encoding="utf-8")
        tree = ast.parse(source)
    except (FileNotFoundError, OSError, SyntaxError) as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)
        return result
    matches: list[ast.AST] = []
    if "." in symbol:
        class_name, method_name = symbol.split(".", 1)
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                matches.extend(
                    child
                    for child in node.body
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and child.name == method_name
                )
    else:
        matches = [
            node
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == symbol
        ]
    if len(matches) != 1:
        result["error"] = "expected one AST function, found " + str(len(matches))
        return result
    function = matches[0]
    segment = ast.get_source_segment(source, function) or ""
    compact = "".join(segment.split())
    checks = {item: "".join(item.split()) in compact for item in required}
    result.update(
        {
            "required_text": checks,
            "AST_sha256": hashlib.sha256(
                ast.dump(function, include_attributes=False).encode()
            ).hexdigest(),
            "passed": all(checks.values()),
        }
    )
    return result


def _parent_controls(parent: Mapping[str, Any], family: Mapping[str, Any]) -> dict[str, Any]:
    accepted = _mapping(parent.get("actual_C5_signed_history_IBP"))
    source_pins = _mapping(accepted.get("source_pins"))
    graph_packet = _mapping(accepted.get("graph_packet"))
    graph_nodes = accepted.get("exact_local_expression_nodes")
    checks = {
        "accepted_parent_receipt_present": bool(accepted),
        "accepted_parent_passed": accepted.get("passed") is True,
        "accepted_parent_source_only": accepted.get("source_only") is True,
        "accepted_parent_source_not_executed": accepted.get("source_executed") is False,
        "accepted_parent_source_family": _family(accepted.get("source_family")) == dict(family),
        "accepted_parent_source_pins_ok": accepted.get("source_pins_ok") is True,
        "accepted_parent_local_graph_present": (
            isinstance(graph_nodes, list)
            and isinstance(graph_packet.get("exact_local_expression_nodes"), list)
        ),
    }
    primitive_checks: dict[str, Any] = {
        "four_ordinary_orders": False,
        "all_B_primitives": False,
        "all_slow_B_y_primitives": False,
        "all_lower_endpoints_zero": False,
        "all_source_phase_derived": False,
        "all_endpoint_zero_mean": False,
        "all_fixed_phase_slow_rows": False,
    }
    if isinstance(graph_nodes, list):
        B_rows = [
            node
            for node in graph_nodes
            if isinstance(node, Mapping)
            and node.get("operation") == "current_original_n1_periodic_phase_primitive"
        ]
        slow_rows = [
            node
            for node in graph_nodes
            if isinstance(node, Mapping)
            and node.get("operation") == "current_original_C5_phase_integral_of_slow_B_y"
        ]
        primitive_checks.update(
            {
                "four_ordinary_orders": sorted(
                    {
                        node.get("ordinary_Z_order")
                        for node in B_rows
                        if node.get("native_variable") in {
                            "true_outer_O2_buffer",
                            "true_outer_O3_transition",
                        }
                    }
                )
                == list(ORDINARY_Z_ORDERS),
                "all_B_primitives": all(
                    any(
                        node.get("native_variable") == "true_outer_" + chart
                        and node.get("ordinary_Z_order") == order
                        and node.get("phase_integral_lower") == 0
                        and node.get("source_phase_node") is not None
                        and node.get("derived_phase_not_bound_variable") is True
                        and node.get("phase_endpoint_zero_mean") is True
                        and node.get("phase_integral_definition")
                        == "C_phi(y,Z,phi)=integral_0^phi B(y,Z,s)ds"
                        for node in B_rows
                    )
                    for chart in CHARTS
                    for order in ORDINARY_Z_ORDERS
                ),
                "all_slow_B_y_primitives": all(
                    any(
                        node.get("native_variable") == "true_outer_" + chart
                        and node.get("ordinary_Z_order") == order
                        and node.get("phase_integral_lower") == 0
                        and node.get("source_phase_node") is not None
                        and node.get("derived_phase_not_bound_variable") is True
                        and node.get("body_role") == "actual_slow_B_y_integrand"
                        and node.get("fixed_phase") is True
                        and node.get("total_y_chain_excluded") is True
                        and node.get("not_a_derivative_of_slow_B_y") is True
                        for node in slow_rows
                    )
                    for chart in CHARTS
                    for order in ORDINARY_Z_ORDERS
                ),
                "all_lower_endpoints_zero": all(
                    node.get("phase_integral_lower") == 0
                    for node in B_rows + slow_rows
                    if node.get("native_variable") in {"true_outer_O2_buffer", "true_outer_O3_transition"}
                ),
                "all_source_phase_derived": all(
                    node.get("phase_variable") == "original_global_phase"
                    and node.get("derived_phase_not_bound_variable") is True
                    for node in B_rows + slow_rows
                    if node.get("native_variable") in {"true_outer_O2_buffer", "true_outer_O3_transition"}
                ),
                "all_endpoint_zero_mean": all(
                    node.get("phase_endpoint_zero_mean") is True
                    for node in B_rows + slow_rows
                    if node.get("native_variable") in {"true_outer_O2_buffer", "true_outer_O3_transition"}
                ),
                "all_fixed_phase_slow_rows": all(
                    node.get("fixed_phase") is True
                    and node.get("total_y_chain_excluded") is True
                    for node in slow_rows
                    if node.get("native_variable") in {"true_outer_O2_buffer", "true_outer_O3_transition"}
                ),
            }
        )
    primitive_checks["passed"] = all(primitive_checks.values())
    checks["accepted_parent_primitive_metadata"] = primitive_checks["passed"]
    # The path is kept explicit so the old local graph cannot be confused
    # with the new outer archive graph.  Only its existence and namespace are
    # checked here; no old build or proof runner is called.
    graph_reference = {
        "graph_path": "actual_C5_signed_history_IBP.graph_packet.exact_local_expression_nodes",
        "namespace": "accepted_parent_local_expression_graph",
        "node_count": len(graph_nodes) if isinstance(graph_nodes, list) else None,
        "source_pins_count": len(source_pins),
        "not_reused_as_outer_graph": True,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "graph_reference": graph_reference,
        "primitive_metadata": primitive_checks,
    }


def _zero_mean_controls(repo: Path, parent: Mapping[str, Any], family: Mapping[str, Any]) -> dict[str, Any]:
    accepted = _mapping(parent.get("actual_C5_signed_history_IBP"))
    old_pin = _mapping(_mapping(accepted.get("source_pins")).get("zero_mean_receipt"))
    expected = old_pin.get("expected_sha256")
    pin = _pin(parent, repo, ZERO_MEAN_RECEIPT, "accepted_source_receipt", expected_override=expected)
    checks: dict[str, Any] = {
        "receipt_pin": pin.get("passed") is True,
        "receipt_loaded": False,
        "receipt_all_passed": False,
        "source_family": False,
        "theorem_passed": False,
        "B_reflection_zero_mean": False,
        "periodic_B_derivative": False,
        "primitive_endpoint_identity": False,
        "slow_y_not_total_y": False,
        "seam_cancellation_not_claimed": False,
    }
    receipt: Mapping[str, Any] = {}
    try:
        receipt = _load(repo / ZERO_MEAN_RECEIPT)
        checks["receipt_loaded"] = True
        theorem = _mapping(receipt.get("exact_source_averaging_theorem"))
        parity = _mapping(theorem.get("original_primitive_parity_theorem"))
        identities = _mapping(parity.get("exact_identities"))
        checks.update(
            {
                "receipt_all_passed": receipt.get("all_passed") is True,
                "source_family": _family(receipt.get("source_family")) == dict(family),
                "theorem_passed": theorem.get("passed") is True,
                "B_reflection_zero_mean": identities.get("B_reflection_zero_mean") is True,
                "periodic_B_derivative": identities.get("periodic_B_derivative") is True,
                "primitive_endpoint_identity": theorem.get("primitive")
                == "G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s)ds; G_j(0)=G_j(1)=0",
                "slow_y_not_total_y": theorem.get("slow_y_not_total_y") is True,
                "seam_cancellation_not_claimed": theorem.get("seam_cancellation_not_claimed") is True,
            }
        )
    except (FileNotFoundError, OSError, ValueError, json.JSONDecodeError):
        pass
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "pin": pin,
        "source_executed": False,
        "source_receipt_not_replayed": True,
    }


def _phase_chain(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    chart: str,
    expected: Mapping[str, Any],
) -> dict[str, Any]:
    phase_id = expected["phase"]
    phase = _node(nodes, phase_id)
    phase_arg_id = phase.get("argument")
    phase_arg = _node(nodes, phase_arg_id)
    args = list(phase_arg.get("arguments", ()))
    n_id = bindings.get("N")
    origin_id = bindings.get("phase_origin")
    inner_id = next((value for value in args if value != n_id), None)
    inner = _node(nodes, inner_id)
    inner_args = list(inner.get("arguments", ()))
    log_id = next(
        (
            value
            for value in inner_args
            if _node(nodes, value).get("operation") == "analytic_unary"
            and _node(nodes, value).get("name") == "log"
        ),
        None,
    )
    negative_id = next(
        (
            value
            for value in inner_args
            if _node(nodes, value).get("operation") == "negative"
        ),
        None,
    )
    log = _node(nodes, log_id)
    ratio = _node(nodes, log.get("argument"))
    negative = _node(nodes, negative_id)
    checks = {
        "phase_is_fractional_part": phase.get("operation") == "analytic_unary"
        and phase.get("name") == "fractional_part",
        "phase_argument_is_product": phase_arg.get("operation") == "product"
        and n_id in args
        and len(args) == 2,
        "phase_uses_common_N": n_id == 18,
        "phase_inner_is_log_minus_origin": inner.get("operation") == "sum"
        and len(inner_args) == 2
        and isinstance(log_id, int)
        and isinstance(negative_id, int),
        "log_is_actual_R_over_Ra": log.get("operation") == "analytic_unary"
        and log.get("name") == "log"
        and ratio.get("operation") == "positive_quotient"
        and ratio.get("numerator") == expected["radius"]
        and ratio.get("denominator") == 65604
        and ratio.get("source_positive_certificate") == "same original positive actual R/Ra",
        "origin_is_common_literal": negative.get("operation") == "negative"
        and negative.get("argument") == origin_id
        and origin_id == 24954,
    }
    return {
        "chart": chart,
        "phase_node": phase_id,
        "phase_argument_node": phase_arg_id,
        "inner_node": inner_id,
        "log_node": log_id,
        "ratio_node": log.get("argument"),
        "checks": checks,
        "passed": all(checks.values()),
        "unwrapped_phase_formula": "N*(log(R/Ra)-phase_origin)",
        "fractional_phase_is_derived": True,
    }


def _radius_chain(nodes: Sequence[Any], bindings: Mapping[str, Any], chart: str, expected: Mapping[str, Any]) -> dict[str, Any]:
    radius = _node(nodes, expected["radius"])
    args = list(radius.get("arguments", ()))
    exp_id = next(
        (
            value
            for value in args
            if _node(nodes, value).get("operation") == "analytic_unary"
            and _node(nodes, value).get("name") == "exp"
        ),
        None,
    )
    exp_node = _node(nodes, exp_id)
    checks = {
        "radius_is_common_Rm_times_exp_offset": radius.get("operation") == "product"
        and 85743 in args
        and isinstance(exp_id, int),
        "radius_uses_chart_offset": exp_node.get("argument") == expected["log_offset"],
        "binding_offset_matches_saved_background": _mapping(bindings.get("background", {})).get(chart, {}).get("log_offset")
        == expected["log_offset"],
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "radius_node": expected["radius"],
        "common_Rm_node": 85743,
        "offset_node": expected["log_offset"],
        "exp_node": exp_id,
    }


def _anchor_checks(nodes: Sequence[Any], bindings: Mapping[str, Any]) -> dict[str, Any]:
    charts = _mapping(bindings.get("charts"))
    endpoints = _mapping(bindings.get("endpoints"))
    background = _mapping(bindings.get("background"))
    result: dict[str, Any] = {}
    for chart, expected in ANCHORS.items():
        actual = _mapping(charts.get(chart))
        loop = _mapping(actual.get("loop"))
        radial = _mapping(actual.get("radial"))
        incoming = _mapping(actual.get("incoming"))
        endpoint = _mapping(endpoints.get(chart))
        endpoint_corr = _mapping(endpoint.get("correction"))
        checks = {
            "coordinate": actual.get("coordinate") == expected["coordinate"],
            "domain": tuple(actual.get("domain", ())) == expected["domain"],
            "domain_lower_node_zero": expected["domain"][0] == 0 and _exact(_node(nodes, 0), 0),
            "domain_upper_exact_value": _exact(_node(nodes, expected["domain"][1]), expected["upper_value"]),
            "phase": actual.get("actual_phase") == expected["phase"],
            "radius": actual.get("radius") == expected["radius"],
            "physical_measure": actual.get("physical_measure") == expected["physical_measure"]
            and _exact(_node(nodes, expected["physical_measure"]), 1),
            "B_rows": tuple(loop.get("B_over_Pstar", ())[:4]) == expected["B"],
            "slow_B_y_rows": tuple(radial.get("slow_B_y", ())[:4]) == expected["slow_B_y"],
            "incoming_m_rows": tuple(incoming.get("m", ())[:4]) == expected["incoming_m"],
            "outgoing_m_rows": tuple(endpoint_corr.get("m", ())[:4]) == expected["outgoing_m"],
            "angle_anchor": loop.get("angle") == expected["angle"],
            "background_log_offset": _mapping(background.get(chart)).get("log_offset") == expected["log_offset"],
            "ordinary_rows_exist": all(
                isinstance(loop.get("B_over_Pstar", ())[order], int)
                and isinstance(radial.get("slow_B_y", ())[order], int)
                for order in ORDINARY_Z_ORDERS
            ),
        }
        phase = _phase_chain(nodes, bindings, chart, expected)
        radius = _radius_chain(nodes, bindings, chart, expected)
        result[chart] = {
            "checks": checks,
            "phase_chain": phase,
            "radius_chain": radius,
            "passed": all(checks.values()) and phase["passed"] and radius["passed"],
            "source_ports": {
                "B": list(expected["B"]),
                "slow_B_y": list(expected["slow_B_y"]),
                "angle": expected["angle"],
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            },
        }
    result["passed"] = all(_mapping(result[chart]).get("passed") is True for chart in CHARTS)
    return result


def _source_phase_and_offset_identity(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    anchors: Mapping[str, Any],
    ast_guards: Mapping[str, Any],
) -> dict[str, Any]:
    o2 = ANCHORS["O2_buffer"]
    o3 = ANCHORS["O3_transition"]
    o2_offset = _node(nodes, o2["log_offset"])
    o3_offset = _node(nodes, o3["log_offset"])
    o2_offset_args = list(o2_offset.get("arguments", ()))
    o3_offset_args = list(o3_offset.get("arguments", ()))
    o2_exp_terms = [
        value
        for value in o2_offset_args
        if _node(nodes, value).get("operation") == "analytic_unary"
        and _node(nodes, value).get("name") == "exp"
    ]
    o3_exp_terms = [
        value
        for value in o3_offset_args
        if _node(nodes, value).get("operation") == "analytic_unary"
        and _node(nodes, value).get("name") == "exp"
    ]
    o2_const = next(
        (value for value in o2_offset_args if _exact(_node(nodes, value), 6)), None
    )
    o3_const = next(
        (value for value in o3_offset_args if _exact(_node(nodes, value), 17)), None
    )
    o2_coordinate_count = o2_offset_args.count(o2["coordinate"])
    o3_coordinate_count = o3_offset_args.count(o3["coordinate"])
    affine_checks = {
        "O2_offset_is_sum": o2_offset.get("operation") == "sum",
        "O3_offset_is_sum": o3_offset.get("operation") == "sum",
        "O2_constant_is_six": o2_const is not None,
        "O3_constant_is_seventeen": o3_const is not None,
        "common_exp40_term": len(o2_exp_terms) == 1
        and len(o3_exp_terms) == 1
        and o2_exp_terms[0] == o3_exp_terms[0]
        and _node(nodes, o2_exp_terms[0]).get("argument") == 2
        and _exact(_node(nodes, 2), 40),
        "O2_native_coefficient_one": o2_coordinate_count == 1,
        "O3_native_coefficient_one": o3_coordinate_count == 1,
        "endpoint_constant_arithmetic": (
            6 + (11 if _exact(_node(nodes, 4), 11) else 0)
            == 17 + (0 if _exact(_node(nodes, 0), 0) else 1)
        ),
    }
    checks = {
        "geometry_ast": _mapping(ast_guards.get("outer_geometry")).get("passed") is True,
        "outer_rc_phase_ast": _mapping(ast_guards.get("outer_Rc_init")).get("passed") is True,
        "O2_upper_exact_node_is_11": _exact(_node(nodes, 4), 11),
        "O3_lower_exact_node_is_0": _exact(_node(nodes, 0), 0),
        "common_Rm_factor": _node(nodes, o2["radius"]).get("arguments", [None])[0]
        == _node(nodes, o3["radius"]).get("arguments", [None])[0]
        == 85743,
        "O2_radius_uses_offset": _radius_chain(nodes, bindings, "O2_buffer", o2)["passed"],
        "O3_radius_uses_offset": _radius_chain(nodes, bindings, "O3_transition", o3)["passed"],
        "O2_phase_chain": _mapping(anchors.get("O2_buffer")).get("phase_chain", {}).get("passed") is True,
        "O3_phase_chain": _mapping(anchors.get("O3_transition")).get("phase_chain", {}).get("passed") is True,
        "common_N_node": bindings.get("N") == 18 and _node(nodes, 18).get("exact_positive_integer") is True,
        "common_phase_origin_node": bindings.get("phase_origin") == 24954,
        "both_physical_Jacobians_one": all(
            _mapping(bindings.get("charts", {}).get(chart)).get("physical_measure") == 1
            for chart in CHARTS
        ),
        "affine_graph_offset_shape": all(affine_checks.values()),
    }
    # This is the source algebra for the two offsets, evaluated only as a
    # symbolic affine identity.  The constants are source constructors, not
    # field values.
    offset_identity = {
        "O2_buffer_at_11": "6 + exp(40) + 11",
        "O3_transition_at_0": "17 + exp(40) + 0",
        "canonical_common_offset": "17 + exp(40)",
        "equal": affine_checks["endpoint_constant_arithmetic"]
        and affine_checks["common_exp40_term"],
        "graph_offset_checks": affine_checks,
    }
    checks["source_offset_affine_seam"] = offset_identity["equal"]
    checks["d_offset_dnative_is_one"] = (
        affine_checks["O2_native_coefficient_one"]
        and affine_checks["O3_native_coefficient_one"]
    )
    checks["unwrapped_phase_derivative_is_common_N"] = (
        checks["d_offset_dnative_is_one"]
        and checks["both_physical_Jacobians_one"]
        and checks["common_N_node"]
    )
    checks["fractional_part_applied_after_common_unwrapped_phase"] = all(
        _mapping(anchors.get(chart, {}).get("phase_chain")).get("checks", {}).get("phase_is_fractional_part")
        is True
        for chart in CHARTS
    )
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "offset_identity": offset_identity,
        "unwrapped_phase_identity": (
            "N*(log(Rm*exp(17+exp(40))/Ra)-phase_origin)"
        ),
        "derivative_identity": "d_native unwrapped_phase = N; d_y unwrapped_phase = N/Jacobian = N",
        "global_phase_identity": "fractional_part(common unwrapped phase)",
        "actual_global_phase_not_numeric": True,
    }


def _branch_continuity(nodes: Sequence[Any], bindings: Mapping[str, Any], ast_guards: Mapping[str, Any]) -> dict[str, Any]:
    background = _mapping(bindings.get("background"))
    endpoints = _mapping(bindings.get("background_endpoints"))
    o2_bg = _mapping(background.get("O2_buffer"))
    o2_ep = _mapping(endpoints.get("O2_buffer"))
    o3_bg = _mapping(background.get("O3_transition"))
    o2_fields = _mapping(o2_ep.get("fields"))
    o2_radial = _mapping(o2_ep.get("radial"))
    o3_fields = _mapping(o3_bg.get("fields"))
    o3_radial = _mapping(o3_bg.get("radial"))
    sigma_t = _node(nodes, 241545)
    sigma_x = _node(nodes, 241559)
    transition_integral = _node(nodes, 241546)
    exponent = _node(nodes, 241550)
    exponent_linear = _node(nodes, 241547)
    exponent_integral = _node(nodes, 241548)
    exponent_negative = _node(nodes, 241549)
    ey_sum = _node(nodes, 241561)
    ey_sigma_product = _node(nodes, 241560)
    o3_v_rows = list(o3_fields.get("V", ()))[:4]
    o3_vy_rows = list(o3_radial.get("V_y", ()))[:4]
    checks: dict[str, Any] = {
        "outer_constructor_ast": _mapping(ast_guards.get("outer_init")).get("passed") is True,
        "O2_endpoint_fields_present": bool(o2_fields),
        "O3_lower_source_fields_present": bool(o3_fields),
        "O3_lower_E_rows_product_shape": len(o3_fields.get("E", ())) >= 4
        and len(o2_fields.get("E", ())) >= 4,
        "O3_lower_E_y_rows_product_shape": len(o3_radial.get("E_y", ())) >= 4
        and len(o2_fields.get("E", ())) >= 4,
        "O2_endpoint_E_y_rows_source_shape": len(o2_radial.get("E_y", ())) >= 4
        and len(o2_bg.get("fields", {}).get("E", ())) >= 4,
        "sigma_zero_branch": sigma_x.get("operation") == "exact_original_source_flat_sigma"
        and sigma_x.get("argument") == 241536
        and "0 for x<=0" in str(sigma_x.get("definition")),
        "sigma_integrand_source_shape": sigma_t.get("operation") == "exact_original_source_flat_sigma"
        and sigma_t.get("argument") == 241544
        and _node(nodes, 241544).get("operation") == "bound_variable"
        and _node(nodes, 241544).get("name") == "true_outer_transition_sigma_t"
        and "0 for x<=0" in str(sigma_t.get("definition")),
        "zero_width_sigma_integral": transition_integral.get("operation") == "definite_integral"
        and transition_integral.get("integrand") == 241545
        and transition_integral.get("variable") == "true_outer_transition_sigma_t"
        and transition_integral.get("lower") == 0
        and transition_integral.get("upper") == 241536,
        "O3_exponent_source_shape": exponent.get("operation") == "sum"
        and exponent.get("arguments") == [241547, 241549]
        and exponent_linear.get("operation") == "product"
        and exponent_linear.get("arguments") == [24, 241536]
        and _exact(_node(nodes, 24), -1, 2)
        and exponent_negative.get("operation") == "negative"
        and exponent_negative.get("argument") == 241548
        and exponent_integral.get("operation") == "product"
        and exponent_integral.get("arguments") == [17, 241546],
        "O3_factor_exp_shape": _node(nodes, 241551).get("operation") == "analytic_unary"
        and _node(nodes, 241551).get("name") == "exp"
        and _node(nodes, 241551).get("argument") == 241550,
        "O3_Ey_factor_source_shape": _node(nodes, 241562).get("operation") == "negative"
        and _node(nodes, 241562).get("argument") == 241561
        and ey_sum.get("operation") == "sum"
        and ey_sum.get("arguments") == [28, 241560]
        and _exact(_node(nodes, 28), 1, 2)
        and ey_sigma_product.get("operation") == "product"
        and ey_sigma_product.get("arguments") == [17, 241559],
        "O3_V_fixed_zero": bool(o3_v_rows)
        and all(_exact(_node(nodes, value), 0) for value in o3_v_rows),
        "O3_Vy_fixed_zero": bool(o3_vy_rows)
        and all(_exact(_node(nodes, value), 0) for value in o3_vy_rows),
        "O2_endpoint_V_zero_wrappers": bool(o2_fields.get("V", ()))
        and all(
            _node(nodes, value).get("operation") == "function_substitution"
            and _node(nodes, value).get("expression") == 0
            for value in list(o2_fields.get("V", ()))[:4]
        ),
        "O2_endpoint_Vy_zero_wrappers": bool(o2_radial.get("V_y", ()))
        and all(
            _node(nodes, value).get("operation") == "function_substitution"
            and _node(nodes, value).get("expression") == 0
            for value in list(o2_radial.get("V_y", ()))[:4]
        ),
    }
    for order in ORDINARY_Z_ORDERS:
        try:
            o2_e = o2_fields["E"][order]
            o2_ey = o2_radial["E_y"][order]
            o3_e = o3_fields["E"][order]
            o3_ey = o3_radial["E_y"][order]
            e_node = _node(nodes, o3_e)
            ey_node = _node(nodes, o3_ey)
            checks["O3_lower_E_rows_product_shape"] &= (
                e_node.get("operation") == "product"
                and list(e_node.get("arguments", ())) == [o2_e, 241551]
            )
            checks["O3_lower_E_y_rows_product_shape"] &= (
                ey_node.get("operation") == "product"
                and list(ey_node.get("arguments", ())) == [o3_e, 241562]
            )
            endpoint_ey = _node(nodes, o2_ey)
            expression = endpoint_ey.get("expression")
            expression_node = _node(nodes, expression)
            checks["O2_endpoint_E_y_rows_source_shape"] &= (
                endpoint_ey.get("operation") == "function_substitution"
                and expression_node.get("operation") == "product"
                and list(expression_node.get("arguments", ())) == [o2_bg.get("fields", {}).get("E", [None])[order], 24]
            )
        except (IndexError, KeyError, TypeError):
            checks["O3_lower_E_rows_product_shape"] = False
            checks["O3_lower_E_y_rows_product_shape"] = False
            checks["O2_endpoint_E_y_rows_source_shape"] = False
    # The source branch has sigma(0)=0 and J(0)=0; therefore exp(-1/2*0-mu*0)=1
    # and -(1/2+mu*0)=-1/2.  These are exact source substitutions, not values
    # obtained by evaluating a field or quadrature.
    checks["E_seam_factor_is_one_by_zero_width_source_algebra"] = (
        checks["sigma_zero_branch"]
        and checks["sigma_integrand_source_shape"]
        and checks["zero_width_sigma_integral"]
        and checks["O3_exponent_source_shape"]
        and checks["O3_factor_exp_shape"]
    )
    checks["Ey_seam_factor_is_minus_half_by_zero_width_source_algebra"] = (
        checks["sigma_zero_branch"] and checks["O3_Ey_factor_source_shape"]
    )
    checks["V_and_Vy_seam_identity"] = (
        checks["O3_V_fixed_zero"]
        and checks["O3_Vy_fixed_zero"]
        and checks["O2_endpoint_V_zero_wrappers"]
        and checks["O2_endpoint_Vy_zero_wrappers"]
    )
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "source_identity": (
            "E_O3(0,Z)=E_O2(11,Z), E_y_O3(0,Z)=E_y_O2(11,Z), V=V_y=0"
        ),
        "no_numeric_endpoint_evaluation": True,
    }


def _zero_width_transport(bindings: Mapping[str, Any], ast_guards: Mapping[str, Any]) -> dict[str, Any]:
    background = _mapping(bindings.get("background"))
    background_endpoints = _mapping(bindings.get("background_endpoints"))
    charts = _mapping(bindings.get("charts"))
    o2_ep = _mapping(background_endpoints.get("O2_buffer"))
    o3_bg = _mapping(background.get("O3_transition"))
    o2_correction = _mapping(_mapping(bindings.get("endpoints")).get("O2_buffer")).get("correction", {})
    o3_chart = _mapping(charts.get("O3_transition"))
    checks = {
        "transport_ast": _mapping(ast_guards.get("transport_geometry")).get("passed") is True,
        "O3_lower_domain_is_zero": list(o3_chart.get("domain", ()))[:1] == [0],
        "background_history_incoming_equal": all(
            _mapping(o3_bg.get("incoming")).get(rate) == _mapping(o2_ep.get("history")).get(rate)
            for rate in RATES
        ),
        "correction_history_incoming_equal": all(
            _mapping(o3_chart.get("incoming")).get(rate) == _mapping(o2_correction).get(rate)
            for rate in RATES
        ),
        "all_six_background_rows_retained": all(
            len(_mapping(o3_bg.get("incoming")).get(rate, ())) >= 6
            and len(_mapping(o2_ep.get("history")).get(rate, ())) >= 6
            for rate in RATES
        ),
        "all_six_correction_rows_retained": all(
            len(_mapping(o3_chart.get("incoming")).get(rate, ())) >= 6
            and len(_mapping(o2_correction).get(rate, ())) >= 6
            for rate in RATES
        ),
        "physical_measure_applied_once_source": _mapping(ast_guards.get("transport_geometry")).get("passed") is True,
        "predecessor_memory_not_reset_source": _mapping(ast_guards.get("transport_geometry")).get("passed") is True,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "rates": list(RATES),
        "zero_width_semantics": {
            "memory": "exp(-rate*(0-0))=1",
            "integral": "integral_0^0 = 0, with physical dy Jacobian present once in the source body",
            "output": "outgoing=incoming",
            "background_rows": "six ordinary source rows are retained",
            "correction_rows": "five rate families, six source rows each, are retained",
        },
        "source_history_aliases": {
            "background": "O3_transition.incoming <- O2_buffer.background_endpoint.history",
            "correction": "O3_transition.incoming <- O2_buffer.endpoint.correction",
        },
        "numeric_transport_not_evaluated": True,
    }


_ZERO_EXPR = ("rat", 0, 1)
_ONE_EXPR = ("rat", 1, 1)


def _reduced_fraction(numerator: Any, denominator: Any) -> tuple[str, int, int]:
    """Return a small exact rational token for the source normalizer."""
    try:
        num = int(numerator)
        den = int(denominator)
    except (TypeError, ValueError):
        return ("rat", 0, 1)
    if den == 0:
        return ("rat", num, den)
    import math

    divisor = math.gcd(abs(num), abs(den)) or 1
    num //= divisor
    den //= divisor
    if den < 0:
        num = -num
        den = -den
    return ("rat", num, den)


class _EndpointNormalizer:
    """Bounded symbolic normalizer for the two actual chart endpoints.

    It expands only graph operations needed for endpoint congruence.  A
    definite integral is discharged before its integrand is visited whenever
    its normalized bounds agree; consequently this never evaluates a saved
    source integrand or quadrature.  O2 field ports are named atoms so the
    O3 lower endpoint must reduce its actual transport sum to those ports.
    """

    def __init__(self, nodes: Sequence[Any], bindings: Mapping[str, Any]) -> None:
        self.nodes = nodes
        self.aliases: dict[int, tuple[Any, ...]] = {}
        self.cache: dict[tuple[Any, ...], tuple[Any, ...]] = {}
        charts = _mapping(bindings.get("background"))
        endpoints = _mapping(bindings.get("background_endpoints"))
        o2_background = _mapping(charts.get("O2_buffer"))
        o2_endpoint = _mapping(endpoints.get("O2_buffer"))
        self.buffer_coordinate = ANCHORS["O2_buffer"]["coordinate"]
        self.moving_aliases: set[int] = set()
        global_pair_ids: set[int] = set()
        # A moving buffer field may receive its endpoint atom only after
        # the actual native11 substitution. Validate each endpoint wrapper
        # before granting that alias; unrelated ports cannot be equated by
        # assigning them the same descriptive token.
        for name, moving_rows in _mapping(o2_background.get("fields")).items():
            endpoint_rows = _mapping(o2_endpoint.get("fields")).get(name, ())
            if name not in ("E", "V", "m", "h", "k", "e", "p", "P0"):
                continue
            if len(moving_rows) != len(endpoint_rows):
                raise ValueError("Incomplete actual buffer endpoint field: " + name)
            for moving, endpoint in zip(moving_rows, endpoint_rows):
                wrapper = _node(nodes, endpoint)
                if not (wrapper.get("operation") == "function_substitution"
                        and wrapper.get("expression") == moving
                        and wrapper.get("variable") == self.buffer_coordinate
                        and wrapper.get("value") == 4
                        and wrapper.get("Z_independent_substitution") is True):
                    raise ValueError("Unproved buffer field endpoint alias: " + name)
                # Earlier graph nodes cannot depend on the later-created
                # native coordinate; global P0 ports remain global atoms.
                if moving >= self.buffer_coordinate:
                    self.moving_aliases.add(moving)
                else:
                    # Global data retain their actual source-node identity;
                    # even a validated endpoint wrapper gets no field label.
                    global_pair_ids.update((moving, endpoint))
        # These are actual O2 source ports.  O3 moving rows deliberately stay
        # expanded so their memory and zero-width integral are checked.
        candidates: dict[int, tuple[Any, ...]] = {}
        duplicate_ids: set[int] = set()
        for source in (o2_background, o2_endpoint):
            fields = _mapping(source.get("fields"))
            for name in ("E", "V", "m", "h", "k", "e", "p", "P0"):
                for order, node_id in enumerate(fields.get(name, ())):
                    if isinstance(node_id, int):
                        token = ("field", name, order)
                        if node_id in candidates:
                            duplicate_ids.add(node_id)
                        else:
                            candidates[node_id] = token
        # Shared zero nodes and repeated source wrappers (notably V) must be
        # expanded, otherwise one graph node could be given several order
        # labels and the endpoint check would become tautological.
        self.aliases = {
            node_id: token
            for node_id, token in candidates.items()
            if node_id not in duplicate_ids and node_id not in global_pair_ids
        }

    @staticmethod
    def _bound_name(name: Any) -> Any:
        if not isinstance(name, str):
            return name
        # The chart suffix is an alpha-renamed source binder, not a physical
        # input.  Canonicalizing it lets the two actual constructors compare.
        return (
            name.replace("O2_buffer", "<outer_chart>")
            .replace("O3_transition", "<outer_chart>")
        )

    def _key(self, node_id: Any, environment: Mapping[Any, Any]) -> tuple[Any, ...]:
        return (
            node_id,
            tuple(sorted((repr(key), repr(value)) for key, value in environment.items())),
        )

    def _payload(self, value: Any, environment: Mapping[Any, Any]) -> Any:
        if isinstance(value, int):
            return self.normalize(value, environment)
        if isinstance(value, list):
            return tuple(self._payload(item, environment) for item in value)
        if isinstance(value, Mapping):
            return tuple(
                (key, self._payload(item, environment))
                for key, item in sorted(value.items(), key=lambda item: str(item[0]))
            )
        return value

    @staticmethod
    def _freeze_metadata(value: Any) -> Any:
        if isinstance(value, Mapping):
            return tuple(
                (key, _EndpointNormalizer._freeze_metadata(item))
                for key, item in sorted(value.items(), key=lambda item: str(item[0]))
            )
        if isinstance(value, (list, tuple)):
            return tuple(_EndpointNormalizer._freeze_metadata(item) for item in value)
        return value

    def normalize(self, node_id: Any, environment: Mapping[Any, Any] | None = None) -> tuple[Any, ...]:
        env = dict(environment or {})
        key = self._key(node_id, env)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        if isinstance(node_id, int) and node_id in env:
            result = env[node_id]
        elif (isinstance(node_id, int) and 0 <= node_id < self.buffer_coordinate
              and _node(self.nodes, node_id).get("operation") != "exact_rational"):
            # The graph is acyclic. Earlier nodes cannot depend on either
            # later native variable being substituted. Preserve their exact
            # identity, including global P0, rather than a descriptive label.
            result = ("shared_global_source_node", node_id)
        elif node_id in self.aliases and (
            node_id not in self.moving_aliases
            or env.get(self.buffer_coordinate) == _reduced_fraction(11, 1)
        ):
            result = self.aliases[node_id]
        elif not isinstance(node_id, int) or not 0 <= node_id < len(self.nodes):
            result = ("raw", repr(node_id))
        else:
            node = _mapping(self.nodes[node_id])
            operation = node.get("operation")
            if operation == "exact_rational":
                result = _reduced_fraction(node.get("numerator"), node.get("denominator", 1))
            elif operation == "sum":
                terms: list[tuple[Any, ...]] = []
                rational_total = (0, 1)
                for argument in node.get("arguments", ()):
                    term = self.normalize(argument, env)
                    if term == _ZERO_EXPR:
                        continue
                    if term and term[0] == "rat" and term[2] != 0:
                        import fractions

                        total = fractions.Fraction(*rational_total) + fractions.Fraction(term[1], term[2])
                        rational_total = (total.numerator, total.denominator)
                        continue
                    if term and term[0] == "negative" and term[1] and term[1][0] == "sum":
                        terms.extend(
                            (
                                _reduced_fraction(-item[1], item[2])
                                if item[0] == "rat"
                                else (item[1] if item[0] == "negative" else ("negative", item))
                            )
                            for item in term[1][1:]
                        )
                        continue
                    if term and term[0] == "sum":
                        terms.extend(term[1:])
                    else:
                        terms.append(term)
                # Combine rational terms introduced while distributing a
                # negative sum, then cancel structurally opposite source
                # terms such as -(offset(0)) + offset(0).
                import fractions

                combined: list[tuple[Any, ...]] = []
                signed: dict[tuple[Any, ...], int] = {}
                for term in terms:
                    if term and term[0] == "rat" and term[2] != 0:
                        total = fractions.Fraction(*rational_total) + fractions.Fraction(term[1], term[2])
                        rational_total = (total.numerator, total.denominator)
                    elif term and term[0] == "negative":
                        signed[term[1]] = signed.get(term[1], 0) - 1
                    else:
                        signed[term] = signed.get(term, 0) + 1
                if rational_total != (0, 1):
                    combined.append(("rat", rational_total[0], rational_total[1]))
                for term, count in signed.items():
                    if count == 1:
                        combined.append(term)
                    elif count == -1:
                        combined.append(("negative", term))
                    elif count:
                        combined.append(("product", ("rat", count, 1), term))
                combined.sort(key=repr)
                result = combined[0] if len(combined) == 1 else ("sum", *combined) if combined else _ZERO_EXPR
            elif operation == "product":
                factors: list[tuple[Any, ...]] = []
                zero = False
                rational_total = (1, 1)
                for argument in node.get("arguments", ()):
                    factor = self.normalize(argument, env)
                    if factor == _ZERO_EXPR:
                        zero = True
                        break
                    if factor == _ONE_EXPR:
                        continue
                    if factor and factor[0] == "rat" and factor[2] != 0:
                        import fractions

                        total = fractions.Fraction(*rational_total) * fractions.Fraction(factor[1], factor[2])
                        rational_total = (total.numerator, total.denominator)
                        continue
                    if factor and factor[0] == "product":
                        factors.extend(factor[1:])
                    else:
                        factors.append(factor)
                if zero:
                    result = _ZERO_EXPR
                else:
                    if rational_total != (1, 1):
                        factors.append(("rat", rational_total[0], rational_total[1]))
                    factors.sort(key=repr)
                    result = factors[0] if len(factors) == 1 else ("product", *factors) if factors else _ONE_EXPR
            elif operation == "negative":
                value = self.normalize(node.get("argument"), env)
                if value and value[0] == "rat":
                    result = _reduced_fraction(-value[1], value[2])
                elif value and value[0] == "negative":
                    result = value[1]
                else:
                    result = ("negative", value)
            elif operation == "analytic_unary":
                value = self.normalize(node.get("argument"), env)
                if node.get("name") == "exp" and value == _ZERO_EXPR:
                    result = _ONE_EXPR
                elif node.get("name") == "fractional_part" and value == _ZERO_EXPR:
                    result = _ZERO_EXPR
                else:
                    result = ("unary", node.get("name"), value)
            elif operation in {"positive_quotient", "quotient"}:
                numerator = self.normalize(node.get("numerator"), env)
                denominator = self.normalize(node.get("denominator"), env)
                if numerator == _ZERO_EXPR:
                    result = _ZERO_EXPR
                elif denominator == _ONE_EXPR:
                    result = numerator
                else:
                    result = ("quotient", numerator, denominator)
            elif operation == "function_substitution":
                variable = node.get("variable")
                value = self.normalize(node.get("value"), env)
                substituted = dict(env)
                substituted[variable] = value
                result = self.normalize(node.get("expression"), substituted)
            elif operation == "definite_integral":
                lower = self.normalize(node.get("lower"), env)
                upper = self.normalize(node.get("upper"), env)
                # Deliberately return before touching the integrand.
                if lower == upper:
                    result = _ZERO_EXPR
                else:
                    result = (
                        "integral",
                        self._bound_name(node.get("variable")),
                        lower,
                        upper,
                        self.normalize(node.get("integrand"), env),
                    )
            elif operation == "bound_variable":
                result = ("bound", self._bound_name(node.get("name")))
            elif operation == "exact_original_source_flat_sigma":
                argument = self.normalize(node.get("argument"), env)
                if argument == _ZERO_EXPR and "0 for x<=0" in str(node.get("definition")):
                    result = _ZERO_EXPR
                else:
                    result = (operation, argument, node.get("definition"))
            else:
                # Preserve source metadata while normalizing known graph
                # references.  This is structural source comparison, not a
                # numerical evaluation of the operation.
                payload: list[tuple[Any, Any]] = []
                for name, value in sorted(node.items(), key=lambda item: str(item[0])):
                    if name == "operation":
                        continue
                    if name in {
                        "argument",
                        "arguments",
                        "expression",
                        "integrand",
                        "lower",
                        "numerator",
                        "denominator",
                        "upper",
                        "value",
                    }:
                        payload.append((name, self._payload(value, env)))
                    else:
                        payload.append((name, self._freeze_metadata(value)))
                result = (operation, tuple(payload))
        self.cache[key] = result
        return result


def _endpoint_environment(chart: str, bindings: Mapping[str, Any]) -> dict[int, tuple[Any, ...]]:
    expected = ANCHORS[chart]
    return {
        expected["coordinate"]: _reduced_fraction(
            expected["upper_value"] if chart == "O2_buffer" else 0,
            1,
        )
    }


def _endpoint_rows(
    normalizer: _EndpointNormalizer,
    bindings: Mapping[str, Any],
) -> dict[str, Any]:
    """Normalize actual O3 lower rows against actual O2 endpoint ports."""
    charts = _mapping(bindings.get("charts"))
    background = _mapping(bindings.get("background"))
    endpoints = _mapping(bindings.get("background_endpoints"))
    o2_background = _mapping(background.get("O2_buffer"))
    o2_endpoint = _mapping(endpoints.get("O2_buffer"))
    o3_background = _mapping(background.get("O3_transition"))
    o3_incoming = _mapping(o3_background.get("incoming"))
    o2_env = _endpoint_environment("O2_buffer", bindings)
    o3_env = _endpoint_environment("O3_transition", bindings)
    field_rows: dict[str, Any] = {}
    for rate in RATES:
        rows: list[dict[str, Any]] = []
        for order in ORDINARY_Z_ORDERS:
            o2_id = _mapping(o2_endpoint.get("fields")).get(rate, ())[order]
            o3_id = _mapping(o3_background.get("fields")).get(rate, ())[order]
            incoming_id = o3_incoming.get(rate, ())[order]
            o2_value = normalizer.normalize(o2_id, o2_env)
            o3_value = normalizer.normalize(o3_id, o3_env)
            incoming_value = normalizer.normalize(incoming_id, o3_env)
            rows.append(
                {
                    "order": order,
                    "O2_buffer_endpoint_node": o2_id,
                    "O3_transition_lower_node": o3_id,
                    "O3_transition_incoming_node": incoming_id,
                    "O2_normal_form": o2_value,
                    "O3_normal_form": o3_value,
                    "O3_incoming_normal_form": incoming_value,
                    "equal": o2_value == o3_value,
                    "equal_to_actual_incoming": o3_value == incoming_value,
                }
            )
        field_rows[rate] = rows
    radial_rows: dict[str, Any] = {}
    for name in ("E_y", "V_y"):
        rows = []
        for order in ORDINARY_Z_ORDERS:
            o2_id = _mapping(o2_endpoint.get("radial")).get(name, ())[order]
            o3_id = _mapping(o3_background.get("radial")).get(name, ())[order]
            o2_value = normalizer.normalize(o2_id, o2_env)
            o3_value = normalizer.normalize(o3_id, o3_env)
            rows.append(
                {
                    "order": order,
                    "O2_buffer_endpoint_node": o2_id,
                    "O3_transition_lower_node": o3_id,
                    "O2_normal_form": o2_value,
                    "O3_normal_form": o3_value,
                    "equal": o2_value == o3_value,
                }
            )
        radial_rows[name] = rows
    o2_radius = _mapping(o2_endpoint).get("radius")
    o3_radius = _mapping(o3_background).get("radius")
    radius_o2 = normalizer.normalize(o2_radius, o2_env)
    radius_o3 = normalizer.normalize(o3_radius, o3_env)
    checks = {
        "actual_chart_nodes_bound": charts.get("O2_buffer", {}).get("coordinate") == ANCHORS["O2_buffer"]["coordinate"]
        and charts.get("O3_transition", {}).get("coordinate") == ANCHORS["O3_transition"]["coordinate"],
        "all_actual_background_fields_equal": all(
            row["equal"] for rows in field_rows.values() for row in rows
        ),
        "all_actual_background_fields_equal_incoming": all(
            row["equal_to_actual_incoming"] for rows in field_rows.values() for row in rows
        ),
        "all_actual_radial_rows_equal": all(
            row["equal"] for rows in radial_rows.values() for row in rows
        ),
        "actual_radius_equal": radius_o2 == radius_o3,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "rates": list(RATES),
        "field_rows": field_rows,
        "radial_rows": radial_rows,
        "radius": {
            "O2_buffer_endpoint_node": o2_radius,
            "O3_transition_lower_node": o3_radius,
            "O2_normal_form": radius_o2,
            "O3_normal_form": radius_o3,
            "equal": radius_o2 == radius_o3,
        },
        "normalization": {
            "O2_buffer_coordinate": o2_env,
            "O3_transition_coordinate": o3_env,
            "definite_integral_equal_bounds": "integral_0^0 is discharged before integrand traversal",
            "numeric_values_evaluated": False,
        },
    }


def _generic_and_root_rows(
    normalizer: _EndpointNormalizer,
    bindings: Mapping[str, Any],
    endpoint_rows: Mapping[str, Any],
) -> dict[str, Any]:
    charts = _mapping(bindings.get("charts"))
    background = _mapping(bindings.get("background"))
    endpoints = _mapping(bindings.get("background_endpoints"))
    o2_chart = _mapping(charts.get("O2_buffer"))
    o3_chart = _mapping(charts.get("O3_transition"))
    o2_background = _mapping(background.get("O2_buffer"))
    o3_background = _mapping(background.get("O3_transition"))
    o2_endpoint = _mapping(endpoints.get("O2_buffer"))
    o2_env = _endpoint_environment("O2_buffer", bindings)
    o3_env = _endpoint_environment("O3_transition", bindings)
    generic_rows: dict[str, Any] = {}
    for name in ("theta_linear", "theta_quadratic", "axial_linear", "axial_quadratic", "p2"):
        rows = []
        for order in ORDINARY_Z_ORDERS:
            o2_id = _mapping(o2_endpoint.get("generic")).get(name, ())[order]
            o3_id = _mapping(o3_background.get("generic")).get(name, ())[order]
            o2_value = normalizer.normalize(o2_id, o2_env)
            o3_value = normalizer.normalize(o3_id, o3_env)
            rows.append(
                {
                    "order": order,
                    "O2_buffer_endpoint_node": o2_id,
                    "O3_transition_lower_node": o3_id,
                    "O2_normal_form_digest": _digest(o2_value),
                    "O3_normal_form_digest": _digest(o3_value),
                    "equal": o2_value == o3_value,
                }
            )
        generic_rows[name] = rows
    root_rows: dict[str, Any] = {}
    for name in ("E", "V", "a", "b", "p2", "t0", "Delta", "kappa"):
        rows = []
        for order in ORDINARY_Z_ORDERS:
            o2_id = _mapping(o2_chart.get("roots")).get(name, ())[order]
            o3_id = _mapping(o3_chart.get("roots")).get(name, ())[order]
            o2_value = normalizer.normalize(o2_id, o2_env)
            o3_value = normalizer.normalize(o3_id, o3_env)
            rows.append(
                {
                    "order": order,
                    "O2_buffer_endpoint_node": o2_id,
                    "O3_transition_lower_node": o3_id,
                    "O2_normal_form_digest": _digest(o2_value),
                    "O3_normal_form_digest": _digest(o3_value),
                    "equal": o2_value == o3_value,
                }
            )
        root_rows[name] = rows
    p2_root_alias_checks = {
        chart: list(_mapping(charts.get(chart)).get("roots", {}).get("p2", ()))[:4]
        == list(_mapping(background.get(chart)).get("generic", {}).get("p2", ()))[:4]
        for chart in CHARTS
    }
    checks = {
        "generic_constructor_inputs_are_seamed": endpoint_rows.get("passed") is True,
        "actual_generic_theta_axial_p2_equal": all(
            row["equal"] for rows in generic_rows.values() for row in rows
        ),
        "actual_root_rows_equal": all(row["equal"] for rows in root_rows.values() for row in rows),
        "actual_root_p2_points_to_background_generic_p2": all(p2_root_alias_checks.values()),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "generic_rows": generic_rows,
        "root_rows": root_rows,
        "actual_root_p2_points_to_background_generic_p2": p2_root_alias_checks,
        "constructor_inputs": {
            "fields": "actual O3 lower field rows reduced to O2 endpoint field ports",
            "radial": "actual O3 lower radial rows reduced to O2 endpoint radial ports",
            "radius": "actual O3 lower radius reduced to O2 endpoint radius",
            "delta": "shared source constructor parameter",
            "Pstar": "shared source constructor argument",
        },
        "normalization": "structural source normalization after actual endpoint substitutions; no source execution",
    }


def _root_input_congruence(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    ast_guards: Mapping[str, Any],
    branch: Mapping[str, Any],
    phase: Mapping[str, Any],
    transport: Mapping[str, Any],
) -> dict[str, Any]:
    charts = _mapping(bindings.get("charts"))
    background = _mapping(bindings.get("background"))
    o2 = _mapping(charts.get("O2_buffer"))
    o3 = _mapping(charts.get("O3_transition"))
    o2_background = _mapping(background.get("O2_buffer"))
    o3_background = _mapping(background.get("O3_transition"))
    o2_p0 = list(_mapping(o2_background.get("fields")).get("P0", ()))
    o3_p0 = list(_mapping(o3_background.get("fields")).get("P0", ()))
    endpoint_normalizer = _EndpointNormalizer(nodes, bindings)
    endpoint_rows = _endpoint_rows(endpoint_normalizer, bindings)
    generic_and_root_rows = _generic_and_root_rows(
        endpoint_normalizer,
        bindings,
        endpoint_rows,
    )
    common = {
        "N": bindings.get("N"),
        "eta": bindings.get("eta"),
        "dstar": bindings.get("dstar"),
        "sc": bindings.get("sc"),
        "phase_origin": bindings.get("phase_origin"),
    }
    checks = {
        "generic_roots_ast": _mapping(ast_guards.get("generic_roots")).get("passed") is True,
        "core_generic_source_ast": _mapping(ast_guards.get("core_generic_functions")).get("passed") is True,
        "branch_continuity": branch.get("passed") is True,
        "common_radius_phase": phase.get("passed") is True,
        "zero_width_transport": transport.get("passed") is True,
        "common_loop_parameters": common == {
            "N": 18,
            "eta": 21,
            "dstar": 22,
            "sc": 10,
            "phase_origin": 24954,
        },
        "roots_have_ordinary_rows": all(
            len(_mapping(_mapping(charts.get(chart)).get("roots")).get(key, ())) >= 4
            for chart in CHARTS
            for key in ("E", "V", "a", "b", "p2", "t0", "Delta")
        ),
        "p2_source_rows_present": all(
            len(_mapping(background.get(chart)).get("generic", {}).get("p2", ())) >= 4
            for chart in CHARTS
        ),
        "same_P0_actual_input_rows": len(o2_p0) >= 4 and o2_p0 == o3_p0,
        "same_Pstar_constructor_argument": _mapping(ast_guards.get("outer_init")).get("passed") is True,
        "actual_endpoint_field_radial_radius_equalities": endpoint_rows.get("passed") is True,
        "actual_generic_theta_axial_p2_equalities": generic_and_root_rows.get("passed") is True,
        "generic_source_depends_on_actual_seamed_inputs": (
            _mapping(ast_guards.get("core_generic_functions")).get("passed") is True
            and _mapping(ast_guards.get("generic_roots")).get("passed") is True
            and endpoint_rows.get("passed") is True
            and generic_and_root_rows.get("passed") is True
        ),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "common_source_parameters": common,
        "constructor": "generic_roots(fields, radial, generic, radius, Pstar)",
        "root_identity": "equal actual root inputs imply equal generic roots by the same deterministic source constructor",
        "actual_root_input_fields": {
            "buffer_background_path": "background.O2_buffer.fields/radial/generic/radius",
            "transition_background_path": "background.O3_transition.fields/radial/generic/radius",
            "P0_rows": {"O2_buffer": o2_p0[:4], "O3_transition": o3_p0[:4]},
        },
        "actual_endpoint_equalities": endpoint_rows,
        "actual_generic_and_root_equalities": generic_and_root_rows,
        "actual_generic_input_congruence": {
            "fields": "actual O3 lower fields reduce to actual O2 endpoint fields for m,h,k,e,p and E,V,P0",
            "radial": "actual O3 lower E_y,V_y reduce to actual O2 endpoint radial rows",
            "radius": "actual O3 lower radius reduces to actual O2 endpoint radius",
            "delta_and_Pstar": "shared actual source constructor arguments are retained by the AST guard",
        },
        "literal_graph_root_node_equality_claimed": False,
        "root_inputs_not_numeric": True,
    }


def _B_constructor_identity(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    ast_guards: Mapping[str, Any],
    root_inputs: Mapping[str, Any],
) -> dict[str, Any]:
    charts = _mapping(bindings.get("charts"))
    o2 = _mapping(charts.get("O2_buffer"))
    o3 = _mapping(charts.get("O3_transition"))
    o2_loop = _mapping(o2.get("loop"))
    o3_loop = _mapping(o3.get("loop"))
    o2_rows = list(o2_loop.get("B_over_Pstar", ()))[:4]
    o3_rows = list(o3_loop.get("B_over_Pstar", ()))[:4]
    o2_angle = _node(nodes, o2_loop.get("angle"))
    o3_angle = _node(nodes, o3_loop.get("angle"))
    checks = {
        "root_inputs_congruent": root_inputs.get("passed") is True,
        "prescribed_loop_ast": _mapping(ast_guards.get("prescribed_loop_functions")).get("passed") is True,
        "generic_arbitrary_phase_ast": _mapping(ast_guards.get("generic_build")).get("passed") is True,
        "whole_Z_arbitrary_phase_ast": _mapping(ast_guards.get("original_loop")).get("passed") is True,
        "same_non_first_micro_branch": (
            "first_micro" not in CHARTS
            and _mapping(ast_guards.get("prescribed_loop_functions")).get("passed") is True
            and o2_loop.get("original_collar_upper") is None
            and o3_loop.get("original_collar_upper") is None
        ),
        "ordinary_B_rows_present": len(o2_rows) == 4 and len(o3_rows) == 4,
        "angle_nodes_are_alpha_renameable": (
            o2_loop.get("angle") != o3_loop.get("angle")
            and o2_angle.get("operation") == "bound_variable"
            and o3_angle.get("operation") == "bound_variable"
            and str(o2_angle.get("name", "")).startswith("true_loop_C5_angle_")
            and str(o3_angle.get("name", "")).startswith("true_loop_C5_angle_")
            and o2_angle.get("name") != o3_angle.get("name")
        ),
    }
    rows = []
    for order, (left, right) in enumerate(zip(o2_rows, o3_rows)):
        rows.append(
            {
                "ordinary_Z_order": order,
                "buffer_B_source_node": left,
                "transition_B_source_node": right,
                "literal_graph_nodes_equal": left == right,
                "canonical_constructor": (
                    "lazy(1/2*E*M), M=-(a*T1hat)/(2*pi)-b*phi_issued, "
                    "angle bound variable alpha-renamed"
                ),
                "same_roots_and_parameters": root_inputs.get("passed") is True,
                "arbitrary_issued_phase": (
                    _mapping(ast_guards.get("generic_build")).get("passed") is True
                    and _mapping(ast_guards.get("original_loop")).get("passed") is True
                ),
            }
        )
    literal_equal = o2_rows == o3_rows
    checks["literal_graph_B_node_equality_is_not_used"] = literal_equal is False
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "rows": rows,
        "formula": "M=-(a*T1hat)/(2*pi)-b*phi; B=lazy(E*M/2)",
        "arbitrary_issued_phase": (
            checks["generic_arbitrary_phase_ast"]
            and checks["whole_Z_arbitrary_phase_ast"]
        ),
        "phase_argument": "phi_issued independent of the derived fractional global phase",
        "literal_graph_B_node_equality": literal_equal,
        "literal_graph_B_node_equality_claimed": False,
        "alpha_renaming_required": checks["angle_nodes_are_alpha_renameable"],
        "source_function_identity_not_numeric": (
            checks["root_inputs_congruent"] and checks["ordinary_B_rows_present"]
        ),
    }


def _primitive_identity(
    zero_mean: Mapping[str, Any],
    b_identity: Mapping[str, Any],
    anchors: Mapping[str, Any],
    parent_controls: Mapping[str, Any],
) -> dict[str, Any]:
    primitive_metadata = _mapping(parent_controls.get("primitive_metadata"))
    checks = {
        "B_function_identity": b_identity.get("passed") is True,
        "accepted_zero_mean_parity": zero_mean.get("passed") is True,
        "accepted_parent_primitive_metadata": primitive_metadata.get("passed") is True,
        "same_lower_gauge": primitive_metadata.get("all_lower_endpoints_zero") is True,
        "ordinary_orders": primitive_metadata.get("four_ordinary_orders") is True,
        "source_phase_is_derived": primitive_metadata.get("all_source_phase_derived") is True,
        "endpoint_zero_mean": primitive_metadata.get("all_endpoint_zero_mean") is True,
        "fixed_phase_slow_rows": primitive_metadata.get("all_fixed_phase_slow_rows") is True,
        "actual_slow_B_y_ports_present": all(
            _mapping(anchors.get(chart)).get("passed") is True for chart in CHARTS
        ),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "definition": "C_chart(Z,phi)=integral_0^phi B_chart(Z,theta)dtheta",
        "seam_identity": "C_buffer(11,Z,phi)=C_transition(0,Z,phi)",
        "gauge": "same lower endpoint 0; accepted parity gives C(0)=C(1)=0",
        "C_y_definition": "C_y_slow=integral_0^phi slow_B_y(y,Z,theta)dtheta at fixed phase",
        "C_y_not_formal_dy_of_slow_B_y": primitive_metadata.get("all_fixed_phase_slow_rows") is True,
        "qualitative_C5_commutation_only": primitive_metadata.get("all_fixed_phase_slow_rows") is True,
        "numeric_C_not_evaluated": True,
    }


def _ast_guards(repo: Path) -> dict[str, Any]:
    return {
        "outer_geometry": _ast_guard(
            repo,
            OUTER_SOURCE,
            "geometry",
            (
                "O2_buffer':g.add(c(6),Y,x)",
                "O3_transition':g.add(c(17),Y,x)",
                "domains={'Rh_reference':(c(-5),g.zero),'O2_buffer':(g.zero,c(11))",
            ),
        ),
        "transport_geometry": _ast_guard(
            repo,
            OUTER_SOURCE,
            "transport_geometry",
            (
                "Duhamel with actual physical log coordinate, once in each measure.",
                "for key,rate in RATES.items():",
                "body=a.scale(a.at(density[key],x,t),g.mul(measure,kernel))",
                "original_outer_Duhamel=True",
                "physical_Jacobian_applied_once=True",
                "predecessor_memory_never_reset=True",
            ),
        ),
        "outer_init": _ast_guard(
            repo,
            OUTER_SOURCE,
            "CurrentOriginalTrueOuterC6C5Functions.__init__",
            (
                "elifname=='O2_buffer':",
                "elifname=='O3_transition':",
                "self.background_endpoints['O2_buffer']['fields']['E']",
                "roots=engine.generic_roots(g,bg,self.Pstar)",
                "phase=g.unary('fractional_part',g.mul(self.N,g.sub(g.unary('log',g.quotient(R,self.Ra,'same original positive actual R/Ra')),self.phase_origin)))",
                "loop=engine.loop_functions(g,roots,phase,self.eta,self.dstar,name,x,self.sc)",
            ),
        ),
        "prescribed_loop_functions": _ast_guard(
            repo,
            LOOP_SOURCE,
            "loop_functions",
            (
                "angle=g.symbol('true_loop_C5_angle_'+chart)",
                "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
                "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
                "phase_endpoints=[0,1]",
                "branch_inverse_is_original_not_support_cap=True",
            ),
        ),
        "generic_roots": _ast_guard(
            repo,
            LOOP_SOURCE,
            "generic_roots",
            (
                "field=chart['fields'];radial=chart['radial'];generic=chart['generic']",
                "p2=generic['p2']",
                "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
                "return dict(E=E,V=V,a=aa,b=bb,p1=p1,p2=p2,t0=t0,Delta=a.sub(kap,a.fixed(2)),kappa=kap)",
            ),
        ),
        "generic_build": _ast_guard(
            repo,
            GENERIC_SOURCE,
            "build",
            (
                "phi=g.node('function_variable',name='phi',domain='periodic_modulo1')",
                "inverse=g.node('monotone_phase_inverse',phase_function=Phi,angle_variable='psi',target_phase=phi,",
                "M=g.sub(g.neg(g.mul(a,T1hat,inv2pi)),g.mul(b,phi))",
                "BB=g.flat(g.mul(half,E,M),Delta,eta)",
            ),
        ),
        "original_loop": _ast_guard(
            repo,
            CONTROL_SOURCE,
            "original_loop",
            (
                "psi_name='loop_angle_'+chart",
                "M=sub(neg(div(mul(a.value,T1hat),twopi)),mul(b.value,phase))",
                "phase_Z_exact_zero=True",
                "branch_inverse_is_original_not_support_cap=True",
            ),
        ),
        "core_generic_functions": _ast_guard(
            repo,
            CORE_SOURCE,
            "generic_functions",
            (
                "E,V,m,h,k,e,p,P0=(a.trim(fields[name])fornamein('E','V','m','h','k','e','p','P0'))",
                "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,'same actual positive core inlet swirl E'),radius)",
            ),
        ),
        "NativeDerivative_init": _ast_guard(
            repo,
            Q4_SOURCE,
            "NativeDerivative.__init__",
            ("self.cache={phase.node:g.zero};self.core_extensions=[]",),
        ),
        "total_radial_functions": _ast_guard(
            repo,
            Q4_SOURCE,
            "total_radial_functions",
            (
                "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
                "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))forrefinrows(chart['loop'][key])))forkeyin('A','B_over_Pstar')}",
                "slow_B_y=slow['B_over_Pstar']",
                "physical_y_equals_native_partial_divided_by_Jacobian=True",
            ),
        ),
        "outer_Rc_init": _ast_guard(
            repo,
            OUTER_RC_SOURCE,
            "WholeZAllNOuterRcFunctions.__init__",
            (
                "pure_parent_chain='axial/buffer <- slope(1); O3 transition <- buffer(11); power <- transition(1)'",
                "original_global_phase='frac(N*(logRm+original_outer_offset-logRa-hb*s_c/2))'",
                "outer_offsets=dict(Rh_reference='6+t',O2_slope='6+t',O2_axial='6+exp(40*t)',O2_buffer='6+exp40+t',O3_transition='17+exp40+t',O3_power='18+exp40+t')",
            ),
        ),
    }


def _source_pins(repo: Path, parent: Mapping[str, Any], old: Mapping[str, Any]) -> dict[str, Any]:
    paths = {
        "outer_graph": (OUTER_GRAPH, "input"),
        "outer_source": (OUTER_SOURCE, "source"),
        "prescribed_loop_source": (LOOP_SOURCE, "source"),
        "generic_loop_source": (GENERIC_SOURCE, "source"),
        "core_source": (CORE_SOURCE, "source"),
        "control_source": (CONTROL_SOURCE, "source"),
        "outer_Rc_source": (OUTER_RC_SOURCE, "source"),
        "modified_Q4_source": (Q4_SOURCE, "source"),
    }
    result = {
        key: _pin(parent, repo, relative, kind)
        for key, (relative, kind) in paths.items()
    }
    old_pin = _mapping(_mapping(old.get("source_pins")).get("zero_mean_receipt"))
    result["zero_mean_receipt"] = _pin(
        parent,
        repo,
        ZERO_MEAN_RECEIPT,
        "accepted_source_receipt",
        expected_override=old_pin.get("expected_sha256"),
    )
    result["passed"] = all(_mapping(value).get("passed") is True for value in result.values())
    return result


def _failure(root: Any, message: str) -> dict[str, Any]:
    return {
        "passed": False,
        GATE: False,
        "source_only": True,
        "source_executed": False,
        "accepted_parent_sha256": PARENT_REPORT_SHA256,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "arbitrary_issued_phase_B_identity": False,
        "anchored_primitive_identity": False,
        "common_unwrapped_phase_identity": False,
        "numeric_C_H_evaluated": False,
        "numeric_C_values": False,
        "numeric_H_values": False,
        "errors": [message],
        "root": str(root),
        "limitations": ["Source seam proof failed closed before any source execution."],
    }


def _build(root: Any, parent: Any = None) -> dict[str, Any]:
    repo = _repo(root)
    parent_map, parent_pin = _load_parent(repo, parent)
    family = _family(parent_map.get("source_family"))
    old = _mapping(parent_map.get("actual_C5_signed_history_IBP"))
    source_pins = _source_pins(repo, parent_map, old)
    ast_guards = _ast_guards(repo)
    outer = _load(repo / OUTER_GRAPH)
    nodes = list(outer.get("exact_graph_nodes", ()))
    bindings = _mapping(outer.get("true_outer_C6_C5_bindings"))
    outer_family = _family(outer.get("source_family"))
    parent_controls = _parent_controls(parent_map, family)
    zero_mean = _zero_mean_controls(repo, parent_map, family)
    anchors = _anchor_checks(nodes, bindings)
    source_phase = _source_phase_and_offset_identity(nodes, bindings, anchors, ast_guards)
    branch = _branch_continuity(nodes, bindings, ast_guards)
    transport = _zero_width_transport(bindings, ast_guards)
    root_inputs = _root_input_congruence(nodes, bindings, ast_guards, branch, source_phase, transport)
    b_identity = _B_constructor_identity(nodes, bindings, ast_guards, root_inputs)
    primitive = _primitive_identity(zero_mean, b_identity, anchors, parent_controls)

    checks = {
        "parent_hash": parent_pin.get("passed") is True,
        "parent_mapping": parent_pin.get("mapping_match") is True,
        "parent_controls": parent_controls.get("passed") is True,
        "source_family_matches_outer": outer_family == family,
        "source_pins": source_pins.get("passed") is True,
        "AST": all(_mapping(value).get("passed") is True for value in ast_guards.values()),
        "zero_mean": zero_mean.get("passed") is True,
        "anchors": anchors.get("passed") is True,
        "phase": source_phase.get("passed") is True,
        "branch": branch.get("passed") is True,
        "transport": transport.get("passed") is True,
        "root_inputs": root_inputs.get("passed") is True,
        "B_identity": b_identity.get("passed") is True,
        "primitive": primitive.get("passed") is True,
    }
    graph_reference = {
        "outer_graph_path": "true_outer_C6_C5_bindings.exact_graph_nodes",
        "outer_graph_namespace": "actual_true_outer_C6_C5_source_graph",
        "accepted_parent_graph_path": parent_controls["graph_reference"]["graph_path"],
        "accepted_parent_graph_namespace": parent_controls["graph_reference"]["namespace"],
        "local_nodes": [],
        "exact_local_expression_nodes": [],
        "new_local_graph_nodes_added": False,
    }
    passed = all(checks.values())
    chart_ports = {
        chart: {
            "coordinate_node": ANCHORS[chart]["coordinate"],
            "domain": list(ANCHORS[chart]["domain"]),
            "domain_lower_node": 0,
            "domain_upper_node": ANCHORS[chart]["domain"][1],
            "domain_upper_exact_value": ANCHORS[chart]["upper_value"],
            "phase_node": ANCHORS[chart]["phase"],
            "radius_node": ANCHORS[chart]["radius"],
            "log_offset_node": ANCHORS[chart]["log_offset"],
            "physical_measure_node": ANCHORS[chart]["physical_measure"],
            "B_rows": list(ANCHORS[chart]["B"]),
            "slow_B_y_rows": list(ANCHORS[chart]["slow_B_y"]),
            "angle_node": ANCHORS[chart]["angle"],
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        }
        for chart in CHARTS
    }
    return {
        "passed": passed,
        GATE: passed,
        "source_binding_passed": passed,
        "source_only": True,
        "source_executed": False,
        "owners_replayed": False,
        "archived_owner_builder_checker_fixture_or_quadrature_executed": False,
        "accepted_parent_sha256": PARENT_REPORT_SHA256,
        "parent_report": PARENT_REPORT,
        "parent_report_pin": parent_pin,
        "source_family": family,
        "source_pins": source_pins,
        "source_pins_ok": source_pins.get("passed") is True,
        "AST": ast_guards,
        "checks": checks,
        "accepted_parent_controls": parent_controls,
        "accepted_zero_mean_controls": zero_mean,
        "source_anchor_checks": anchors,
        "phase_and_offset_identity": source_phase,
        "background_branch_continuity": branch,
        "zero_width_transport_identity": transport,
        "root_input_congruence": root_inputs,
        "arbitrary_issued_phase_B_identity": b_identity.get("passed") is True,
        "B_constructor_identity": b_identity,
        "anchored_primitive_identity": primitive.get("passed") is True,
        "primitive_identity": primitive,
        "common_unwrapped_phase_identity": source_phase.get("passed") is True,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "actual_ports": chart_ports,
        "graph_reference": graph_reference,
        "numeric_C_H_evaluated": False,
        "numeric_C_values": False,
        "numeric_H_values": False,
        "cap_contraction_or_improvement": False,
        "full_n1_or_terminal_closure": False,
        "physical_scale_materialized": False,
        "limitations": [
            "This receipt proves only the source-function seam for O2_buffer at 11 and O3_transition at 0.",
            "The proof covers ordinary Z orders 0 through 3 and an arbitrary issued phase.",
            "No numeric C, H, cap, contraction, full n1, or terminal closure is installed.",
            "The literal B graph nodes remain distinct because the bound angle names are alpha-renamed.",
            "C_y is the fixed-phase phase primitive of actual slow_B_y; no mixed slow-y numeric bound is claimed.",
        ],
        "root": str(repo),
        "directory": str(repo / "experiments" / "root_st073"),
    }


def build(root: Any, parent: Any = None) -> dict[str, Any]:
    """Return a fail-closed source seam receipt without executing source owners."""
    try:
        return _build(root, parent)
    except Exception as exc:  # fail closed for changed archives or incomplete source guards
        return _failure(root, type(exc).__name__ + ": " + str(exc))


def build_original_C5_mean_seam(root: Any, parent: Any = None) -> dict[str, Any]:
    return build(root, parent)


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "OUTER_GRAPH",
    "OUTER_SOURCE",
    "LOOP_SOURCE",
    "GENERIC_SOURCE",
    "CONTROL_SOURCE",
    "CORE_SOURCE",
    "OUTER_RC_SOURCE",
    "Q4_SOURCE",
    "ZERO_MEAN_RECEIPT",
    "ORDINARY_Z_ORDERS",
    "ANCHORS",
    "build",
    "build_original_C5_mean_seam",
]
