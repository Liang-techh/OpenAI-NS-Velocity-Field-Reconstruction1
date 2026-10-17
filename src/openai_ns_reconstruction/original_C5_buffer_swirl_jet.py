"""Source bound O2 buffer swirl jets and a small directed interval adapter.

The saved C5 source graph gives the buffer background in the ordinary
``Z`` jet as a positive axial amplitude times ``(1 + Z**2)**-1``.  This
module binds that graph to a closed form only after checking the actual graph
nodes, the source AST, and the saved reflected-sigma receipt.  The new
evaluator is deliberately local: it computes directed enclosures of the
source E jet and the source roots ``a`` and ``t0``.  It does not evaluate a
C5 phase inverse, a pressure p2 jet, C, H, or a quadrature.

All public interval inputs must be exact integers, Fractions, decimal strings,
or directed intervals.  Ordinary binary floats are rejected so a caller
cannot silently turn this source adapter into a point sampler.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

from mpmath.ctx_iv import MPIntervalContext
import mpmath as mp
from .original_n1_complex_source_bounds import _serialization_roundtrip


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json"
)
PARENT_SHA256 = (
    "588d221fd082d5f84d94f1812eb6342bb6647847b114dac458c0165324a063b1"
)
OUTER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
OUTER_GRAPH_SHA256 = (
    "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
)
OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
OUTER_SOURCE_SHA256 = (
    "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3"
)
LOOP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
LOOP_SOURCE_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)
SIGMA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_interval_outer_slope_field.py"
)
SIGMA_SOURCE_SHA256 = (
    "906e5efa703f53f2d61762d3b7dca38d521ee3649e00962a1f0a81ddfd3254d5"
)
SIGMA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.json"
)
SIGMA_RECEIPT_SHA256 = (
    "c96ece0110427b141f870771d07e1a2db4a7b85ae0a4bd0b3ad7d8ff38bf8286"
)
RH_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join.py"
)
RH_SOURCE_SHA256 = (
    "2a4097aef5ab72d931c91405eb8b9e672abf3a2baa8392f46a6dde2a6a748b5e"
)
RH_NAME = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join.json"
)
RH_NAME_SHA256 = (
    "5775e328db5147e7b28d0b5eda46ae7a64a840cf72c7c71413d08480993709c4"
)
RH_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join_check.json"
)
RH_RECEIPT_SHA256 = (
    "6cd643534012b8d265a74b61be7a3ef90ac76b69f94a773c61686dc49a383dbf"
)
RH_BINDINGS_DIGEST = (
    "59f362c396fd08e447ab17b1d7ecb595e7cacc472c2585a41fe92b9bda991a56"
)
RH_SOURCE_BINDINGS_AST_SHA256 = (
    "2ab5431d4e144c953724eab7460c1ba6463a8dcf7c8f24c4aa92f5388b80aa23"
)
RH_JOIN_INIT_AST_SHA256 = (
    "c369aed654fee6cc9ab24d45454bcb3838111e0dea2c81aa3cb39e88b5e9145a"
)
RH_REFERENCE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N.py"
)
RH_REFERENCE_SOURCE_SHA256 = (
    "8c5114beb4fe896d0622d4b1069d0507240a5a97174dfd653aa1124c13652834"
)
RH_REFERENCE_AST_SHA256 = (
    "089024c68c2cba935cdb9a0da48d33c55aaf092344cfa63d77a1a60944f0547c"
)

ORDINARY_Z_ORDERS = (0, 1, 2, 3)
MIXED_Y_ORDERS = (0, 1)
O2_COORDINATE = 232895
O2_DOMAIN_NODE = (0, 4)
O2_DOMAIN_VALUE = (0, 11)
O2_PHASE_NODE = 233900
O2_EXP_NODE = 232904
O2_EXP_ARGUMENT_NODE = 232903
O2_E_ROWS = (232905, 232906, 232907, 232908, 232909, 232910)
O2_EY_ROWS = (232912, 232913, 232914, 232915, 232916, 232917)
O2_ENDPOINT_AXIAL_E_ROWS = (224038, 224039, 224040, 224041, 224042, 224043)
O2_A_ROWS = (233825, 233829, 233835, 233842, 233850, 233859)
O2_T0_ZERO = 0
EXACT_MINUS_HALF_NODE = 24
EXACT_TWO_NODE = 48


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _load(path: Path) -> Any:
    raw = path.read_bytes()
    if path.name.endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _repo(root: Any) -> Path:
    path = Path(root)
    if path.is_file():
        path = path.parent
    candidates = (path, path.parent, Path(__file__).resolve().parents[2])
    for candidate in candidates:
        if (candidate / "src").is_dir() and (candidate / "experiments").is_dir():
            return candidate.resolve()
    raise ValueError("repository root with src/ and experiments/ is required")


def _node(nodes: Sequence[Any], node_id: Any) -> Mapping[str, Any]:
    if isinstance(node_id, int) and 0 <= node_id < len(nodes):
        item = nodes[node_id]
        return item if isinstance(item, Mapping) else {}
    return {}


def _exact(node: Mapping[str, Any], numerator: int, denominator: int = 1) -> bool:
    return (
        node.get("operation") == "exact_rational"
        and node.get("numerator") == numerator
        and node.get("denominator") == denominator
    )


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    try:
        actual = _sha(repo / relative)
    except (FileNotFoundError, OSError):
        actual = None
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
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
        for item in ast.walk(tree):
            if isinstance(item, ast.ClassDef) and item.name == class_name:
                matches.extend(
                    child
                    for child in item.body
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and child.name == method_name
                )
    else:
        matches = [
            item
            for item in ast.walk(tree)
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
            and item.name == symbol
        ]
    if len(matches) != 1:
        result["error"] = "expected one AST function, found " + str(len(matches))
        return result
    function = matches[0]
    segment = ast.get_source_segment(source, function) or ""
    compact = "".join(segment.split())
    result["required_text"] = {
        item: "".join(item.split()) in compact for item in required
    }
    result["AST_sha256"] = hashlib.sha256(
        ast.dump(function, include_attributes=False).encode()
    ).hexdigest()
    result["passed"] = all(result["required_text"].values())
    return result


def _parent_pin(repo: Path, parent: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    saved = _load(repo / PARENT)
    if parent is None:
        selected = saved
    elif isinstance(parent, (str, Path)):
        selected = _load(repo / str(parent))
    else:
        selected = parent
    mapping_match = _digest(saved) == _digest(selected)
    actual = _sha(repo / PARENT)
    pin = {
        "path": PARENT,
        "sha256": actual,
        "expected_sha256": PARENT_SHA256,
        "passed": actual == PARENT_SHA256,
        "mapping_match": mapping_match,
        "source_executed": False,
    }
    return selected, pin


def _endpoints(value: Any) -> tuple[Any, Any]:
    lo, hi = value._mpi_
    return mp.make_mpf(lo), mp.make_mpf(hi)


def _record(value: Any) -> dict[str, Any]:
    lo, hi = _endpoints(value)
    return {
        "lower_exact_mpf_tuple": list(map(int, lo._mpf_)),
        "upper_exact_mpf_tuple": list(map(int, hi._mpf_)),
        "lower_decimal_display": mp.nstr(lo, 28),
        "upper_decimal_display": mp.nstr(hi, 28),
        "directed_interval": True,
    }


def _scalar(ctx: MPIntervalContext, value: Any) -> Any:
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, bool):
        return ctx.mpf(int(value))
    if isinstance(value, int):
        return ctx.mpf(value)
    if isinstance(value, Fraction):
        return ctx.mpf(value.numerator) / value.denominator
    if isinstance(value, float):
        raise TypeError("binary float inputs are not accepted")
    if isinstance(value, str):
        try:
            rational = Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise TypeError("string input must be a finite decimal rational") from exc
        return ctx.mpf(rational.numerator) / rational.denominator
    raise TypeError("interval endpoints must be exact integers, Fractions, or strings")


def _interval(ctx: MPIntervalContext, value: Any) -> Any:
    """Convert an exact endpoint or directed interval record to an ivmpf."""
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        lower_tuple = value.get("lower_exact_mpf_tuple")
        upper_tuple = value.get("upper_exact_mpf_tuple")
        if lower_tuple is not None and upper_tuple is not None:
            left = mp.make_mpf(tuple(int(part) for part in lower_tuple))
            right = mp.make_mpf(tuple(int(part) for part in upper_tuple))
            return ctx.mpf([left, right])
        if "lower" in value and "upper" in value:
            return ctx.mpf([_scalar(ctx, value["lower"]), _scalar(ctx, value["upper"])])
        if "upper" in value:
            return _scalar(ctx, value["upper"])
        raise TypeError("interval mapping needs lower/upper endpoints")
    if isinstance(value, (tuple, list)):
        if len(value) != 2:
            raise TypeError("interval tuple/list must contain exactly two endpoints")
        return ctx.mpf([_scalar(ctx, value[0]), _scalar(ctx, value[1])])
    return _scalar(ctx, value)


def _jet_key(z_order: int, y_order: int) -> str:
    return f"Z{z_order}_y{y_order}"


def _derivative_factors(ctx: MPIntervalContext, z: Any) -> tuple[Any, ...]:
    denominator = 1 + z * z
    if _endpoints(denominator)[0] <= 0:
        raise ValueError("Z interval is outside the positive denominator tube")
    d2 = denominator * denominator
    d3 = d2 * denominator
    d4 = d3 * denominator
    return (
        1 / denominator,
        -2 * z / d2,
        (6 * z * z - 2) / d3,
        (24 * z - 24 * z * z * z) / d4,
    )


def evaluate_buffer_swirl(
    ctx: MPIntervalContext,
    xbox: Any,
    zbox: Any,
) -> dict[str, Any]:
    """Return directed O2-buffer E, a, t0, and b jet boxes.

    ``xbox`` is the physical O2 buffer coordinate and ``zbox`` is the
    ordinary axial variable.  The source Jacobian is one, so no extra
    physical dy factor is introduced here.  The result contains raw
    ``ivmpf`` values for programmatic composition; :func:`build` converts
    them to exact mpf tuple records for JSON.
    """
    if not isinstance(ctx, MPIntervalContext):
        raise TypeError("ctx must be an mpmath MPIntervalContext")
    x = _interval(ctx, xbox)
    z = _interval(ctx, zbox)
    xlo, xhi = _endpoints(x)
    zlo, zhi = _endpoints(z)
    if xlo < 0 or xhi > 11:
        raise ValueError("O2 buffer coordinate interval must lie in [0,11]")
    if zlo <= -1 or zhi >= 1:
        raise ValueError("Z interval must lie strictly inside (-1,1)")

    exp40 = ctx.exp(ctx.mpf(40))
    log_scale = ctx.mpf(3) / 10 - exp40 / 2 - x / 2
    scale = ctx.exp(log_scale)
    factors = _derivative_factors(ctx, z)
    rows: dict[str, Any] = {}
    a_rows: dict[str, Any] = {}
    b_rows: dict[str, Any] = {}
    t0_rows: dict[str, Any] = {}
    for order in ORDINARY_Z_ORDERS:
        value = scale * factors[order]
        for y_order in MIXED_Y_ORDERS:
            key = _jet_key(order, y_order)
            rows[key] = value if y_order == 0 else -value / 2
            a_rows[key] = ctx.mpf(2) if order == 0 and y_order == 0 else ctx.mpf(0)
            b_rows[key] = ctx.mpf(0)
            t0_rows[key] = ctx.mpf(0)
    return {
        "xbox": x,
        "zbox": z,
        "log_scale": log_scale,
        "scale": scale,
        "log_E": log_scale - ctx.log(1 + z * z),
        "E_rows": rows,
        "a_rows": a_rows,
        "b_rows": b_rows,
        "t0_rows": t0_rows,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "mixed_y_orders": list(MIXED_Y_ORDERS),
        "physical_jacobian": ctx.mpf(1),
        "source_formula": (
            "E=exp(3/10-exp(40)/2-x/2)/(1+Z^2); "
            "a=2; b=t0=0"
        ),
        "source_executed": False,
    }


def _record_rows(values: Mapping[str, Any]) -> dict[str, Any]:
    return {key: _record(value) for key, value in values.items()}


def _graph_checks(graph: Mapping[str, Any]) -> dict[str, Any]:
    nodes = graph.get("exact_graph_nodes", ())
    bindings = _mapping(graph.get("true_outer_C6_C5_bindings"))
    charts = _mapping(bindings.get("charts"))
    chart = _mapping(charts.get("O2_buffer"))
    background = _mapping(_mapping(bindings.get("background")).get("O2_buffer"))
    endpoint_axial = _mapping(
        _mapping(bindings.get("background_endpoints")).get("O2_axial")
    )
    fields = _mapping(background.get("fields"))
    radial = _mapping(background.get("radial"))
    roots = _mapping(chart.get("roots"))
    checks: dict[str, Any] = {
        "graph_has_exact_nodes": isinstance(nodes, list) and len(nodes) > O2_E_ROWS[-1],
        "actual_chart_coordinate": chart.get("coordinate") == O2_COORDINATE,
        "actual_domain_node_ids": chart.get("domain") == list(O2_DOMAIN_NODE),
        "actual_domain_values": _exact(_node(nodes, 0), 0) and _exact(_node(nodes, 4), 11),
        "actual_physical_jacobian_one": chart.get("physical_measure") == 1,
        "actual_phase_node": chart.get("actual_phase") == O2_PHASE_NODE,
        "actual_E_rows": tuple(fields.get("E", ())[:6]) == O2_E_ROWS,
        "actual_Ey_rows": tuple(radial.get("E_y", ())[:6]) == O2_EY_ROWS,
        "actual_V_rows_exact_zero": tuple(fields.get("V", ())[:6]) == (0,) * 6,
        "actual_Vy_rows_exact_zero": tuple(radial.get("V_y", ())[:6]) == (0,) * 6,
        "actual_root_E_rows": tuple(roots.get("E", ())[:6]) == O2_E_ROWS,
        "actual_root_a_rows": tuple(roots.get("a", ())[:6]) == O2_A_ROWS,
        "actual_root_b_rows_exact_zero": tuple(roots.get("b", ())[:6]) == (0,) * 6,
        "actual_root_t0_rows_exact_zero": tuple(roots.get("t0", ())[:6]) == (0,) * 6,
        "actual_source_eta_node": bindings.get("eta") == 21,
        "actual_source_dstar_node": bindings.get("dstar") == 22,
        "actual_common_N_node": bindings.get("N") == 18,
        "seed_Rh_identity": _mapping(bindings.get("seed_identity")).get(
            "canonical_normalized_Rh"
        )
        == "E=exp(-1/2)/(1+Z^2); V=m=4Z/Pstar; h=5E/8; k=h*V; e=V^2-5E^2/12; p=5E^2/2",
    }
    e_rows_ok = True
    ey_rows_ok = True
    for index, row_id in enumerate(O2_E_ROWS):
        row = _node(nodes, row_id)
        e_rows_ok &= row.get("operation") == "product"
        e_rows_ok &= row.get("arguments") == [O2_ENDPOINT_AXIAL_E_ROWS[index], O2_EXP_NODE]
    for index, row_id in enumerate(O2_EY_ROWS):
        row = _node(nodes, row_id)
        ey_rows_ok &= row.get("operation") == "product"
        ey_rows_ok &= row.get("arguments") == [O2_E_ROWS[index], EXACT_MINUS_HALF_NODE]
    checks["E_rows_are_actual_axial_endpoint_times_buffer_exp"] = bool(e_rows_ok)
    checks["Ey_rows_are_actual_E_times_minus_half"] = bool(ey_rows_ok)
    exp_argument = _node(nodes, O2_EXP_ARGUMENT_NODE)
    exp_node = _node(nodes, O2_EXP_NODE)
    checks["buffer_exponent_source_chain"] = (
        exp_argument.get("operation") == "product"
        and exp_argument.get("arguments") == [EXACT_MINUS_HALF_NODE, O2_COORDINATE]
        and exp_node.get("operation") == "analytic_unary"
        and exp_node.get("name") == "exp"
        and exp_node.get("argument") == O2_EXP_ARGUMENT_NODE
    )
    slope_endpoint_e = _node(nodes, 215023)
    slope_factor = _node(nodes, 214131)
    slope_argument = _node(nodes, 214130)
    slope_integral_factor = _node(nodes, 214128)
    slope_integral = _node(nodes, 214126)
    axial_endpoint_e = _node(nodes, 224038)
    axial_expression = _node(nodes, 223070)
    axial_exp = _node(nodes, 223069)
    axial_exponent = _node(nodes, 223068)
    axial_distance = _node(nodes, 223067)
    axial_exp40 = _node(nodes, 223061)
    checks["slope_endpoint_source_chain"] = (
        tuple(_mapping(_mapping(bindings.get("background_endpoints")).get("O2_slope"))
              .get("fields", {}).get("E", ())[:1])
        == (215023,)
        and slope_endpoint_e.get("operation") == "function_substitution"
        and slope_endpoint_e.get("expression") == 214132
        and slope_endpoint_e.get("value") == 1
        and _node(nodes, 214132).get("operation") == "product"
        and _node(nodes, 214132).get("arguments") == [206109, 214131]
        and slope_factor.get("operation") == "analytic_unary"
        and slope_factor.get("name") == "exp"
        and slope_factor.get("argument") == 214130
        and slope_argument.get("operation") == "sum"
        and slope_argument.get("arguments") == [214127, 214129]
        and slope_integral_factor.get("operation") == "product"
        and _node(nodes, 214129).get("operation") == "negative"
        and _node(nodes, 214129).get("argument") == 214128
        and _node(nodes, 214127).get("arguments") == [31355, 214116]
        and slope_integral_factor.get("arguments") == [31357, 214126]
        and _exact(_node(nodes, 31355), 1, 10)
        and _exact(_node(nodes, 31357), 3, 5)
        and slope_integral.get("operation") == "definite_integral"
        and slope_integral.get("lower") == 0
        and slope_integral.get("upper") == 214116
        and slope_integral.get("exact_function_integral") is True
        and slope_integral.get("exact_original_sigma_integral") is True
        and slope_integral.get("numerical_value_not_installed") is True
    )
    checks["axial_endpoint_source_chain"] = (
        tuple(endpoint_axial.get("fields", {}).get("E", ())[:1]) == (224038,)
        and axial_endpoint_e.get("operation") == "function_substitution"
        and axial_endpoint_e.get("expression") == 223070
        and axial_endpoint_e.get("value") == 1
        and axial_expression.get("operation") == "product"
        and axial_expression.get("arguments") == [215023, 223069]
        and axial_exp.get("operation") == "analytic_unary"
        and axial_exp.get("name") == "exp"
        and axial_exp.get("argument") == 223068
        and axial_exponent.get("operation") == "product"
        and axial_exponent.get("arguments") == [EXACT_MINUS_HALF_NODE, 223067]
        and axial_distance.get("operation") == "sum"
        and axial_distance.get("arguments") == [223061, 25085]
        and axial_exp40.get("operation") == "analytic_unary"
        and axial_exp40.get("name") == "exp"
        and _node(nodes, 223060).get("arguments") == [2, 223058]
        and _exact(_node(nodes, 2), 40)
        and _node(nodes, 25085).get("operation") == "negative"
        and _node(nodes, 25085).get("argument") == 1
        and _exact(_node(nodes, 1), 1)
    )
    checks["canonical_Rh_endpoint_normalization_anchor"] = (
        _node(nodes, 206109).get("operation") == "function_substitution"
        and _node(nodes, 206109).get("value") == 0
        and _node(nodes, 206109).get("Z_independent_substitution") is True
    )
    checks["a_rows_are_source_positive_quotients"] = all(
        _node(nodes, row_id).get("operation") == "positive_quotient"
        and _node(nodes, row_id).get("denominator") == O2_E_ROWS[0]
        and _node(nodes, row_id).get("source_positive_certificate")
        == "same strictly positive actual prescribed swirl E"
        for row_id in O2_A_ROWS
    )
    checks["a_order_zero_uses_E_minus_two_Ey"] = (
        _node(nodes, 233825).get("numerator") == 233819
        and _node(nodes, 233819).get("operation") == "sum"
        and _node(nodes, 233819).get("arguments") == [O2_E_ROWS[0], 233813]
        and _node(nodes, 233813).get("operation") == "negative"
        and _node(nodes, 233813).get("argument") == 233807
        and _node(nodes, 233807).get("operation") == "product"
        and _node(nodes, 233807).get("arguments") == [O2_EY_ROWS[0], EXACT_TWO_NODE]
    )
    checks["all_source_root_checks"] = all(
        checks[key]
        for key in (
            "actual_root_a_rows",
            "actual_root_b_rows_exact_zero",
            "actual_root_t0_rows_exact_zero",
            "a_rows_are_source_positive_quotients",
            "a_order_zero_uses_E_minus_two_Ey",
            "slope_endpoint_source_chain",
            "axial_endpoint_source_chain",
            "canonical_Rh_endpoint_normalization_anchor",
        )
    )
    checks["passed"] = all(checks.values())
    return checks


def build(root: Any, parent: Any = None) -> dict[str, Any]:
    """Build a source-only receipt and a fresh closed-form E jet enclosure."""
    repo = _repo(root)
    errors: list[str] = []
    try:
        accepted, parent_pin = _parent_pin(repo, parent)
    except Exception as exc:  # fail closed while retaining a reviewable receipt
        accepted = {}
        parent_pin = {
            "path": PARENT,
            "expected_sha256": PARENT_SHA256,
            "passed": False,
            "mapping_match": False,
            "source_executed": False,
            "error": type(exc).__name__ + ": " + str(exc),
        }
        errors.append(parent_pin["error"])
    source_pins = {
        "parent": parent_pin,
        "outer_graph": _pin(repo, OUTER_GRAPH, OUTER_GRAPH_SHA256),
        "outer_source": _pin(repo, OUTER_SOURCE, OUTER_SOURCE_SHA256),
        "prescribed_loop_source": _pin(repo, LOOP_SOURCE, LOOP_SOURCE_SHA256),
        "sigma_source": _pin(repo, SIGMA_SOURCE, SIGMA_SOURCE_SHA256),
        "sigma_receipt": _pin(repo, SIGMA_RECEIPT, SIGMA_RECEIPT_SHA256),
        "Rh_source": _pin(repo, RH_SOURCE, RH_SOURCE_SHA256),
        "Rh_name": _pin(repo, RH_NAME, RH_NAME_SHA256),
        "Rh_receipt": _pin(repo, RH_RECEIPT, RH_RECEIPT_SHA256),
        "Rh_reference_source": _pin(
            repo, RH_REFERENCE_SOURCE, RH_REFERENCE_SOURCE_SHA256
        ),
    }
    try:
        graph = _load(repo / OUTER_GRAPH)
    except Exception as exc:
        graph = {}
        errors.append("graph load: " + type(exc).__name__ + ": " + str(exc))
    graph_checks = _graph_checks(graph)
    graph_family = graph.get("source_family")
    accepted_parent_sections = {
        key: _mapping(accepted.get(key)).get("passed") is True
        for key in (
            "checks",
            "actual_C5_signed_history_IBP",
            "actual_C5_O2_O3_seam",
            "joined_C5_signed_history",
            "actual_C5_terminal_flat",
            "directed_C5_angle_primitive",
        )
    }
    sigma_ast = _ast_guard(
        repo,
        SIGMA_SOURCE,
        "transition_integrals",
        ("if y==1:j=c.mpf('.5')",),
    )
    outer_ast = _ast_guard(
        repo,
        OUTER_SOURCE,
        "CurrentOriginalTrueOuterC6C5Functions.__init__",
        (
            "for module in (original_leading,original_outer,original_geometry,original_Rh):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))",
            "original_Rh_leading_identity=original_Rh.source_bindings()",
            "Rh_receipt=parent.linear.read(original_Rh.RECEIPT)",
            "Rh_source=parent.linear.read(original_Rh.NAME)",
            "original_whole_Z_leading_identity_receipt=original_Rh.RECEIPT",
            "ifname=='O2_slope':",
            "exact_original_sigma_integral=True",
            "ifname=='O2_axial':",
            "distance=g.sub(g.unary('exp',g.mul(g.constant(40),x)),g.one)",
            "ifname=='O2_buffer':",
            "g.unary('exp',g.mul(g.constant('-1/2'),x))",
            "Ey=a.scale(E,'-1/2')",
        ),
    )
    rh_bindings_ast = _ast_guard(
        repo,
        RH_SOURCE,
        "source_bindings",
        (
            "reference_binding=reference.source_bindings()",
            "returndict(partial_primitives=",
            "current_all_N_outer_transport=require_statements(",
        ),
    )
    rh_join_ast = _ast_guard(
        repo,
        RH_SOURCE,
        "WholeZRhFunctionalJoin.__init__",
        ("self.bindings=source_bindings()", "self.proof=exact_join_proof()"),
    )
    rh_reference_ast = _ast_guard(
        repo,
        RH_REFERENCE_SOURCE,
        "background_cell",
        ("logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10",),
    )
    rh_reference_ast["expected_AST_sha256"] = RH_REFERENCE_AST_SHA256
    rh_reference_ast["AST_hash_matches"] = (
        rh_reference_ast.get("AST_sha256") == RH_REFERENCE_AST_SHA256
    )
    rh_reference_ast["passed"] = (
        rh_reference_ast.get("passed") is True
        and rh_reference_ast["AST_hash_matches"]
    )
    for guard, expected in (
        (rh_bindings_ast, RH_SOURCE_BINDINGS_AST_SHA256),
        (rh_join_ast, RH_JOIN_INIT_AST_SHA256),
    ):
        guard["expected_AST_sha256"] = expected
        guard["AST_hash_matches"] = guard.get("AST_sha256") == expected
        guard["passed"] = guard.get("passed") is True and guard["AST_hash_matches"]
    loop_ast = _ast_guard(
        repo,
        LOOP_SOURCE,
        "generic_roots",
        (
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2))",
            "B=a.scale(a.trim(radial['V_y']),2)",
            "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
            "kap=a.add(aa,a.div(a.mul(bb,bb),aa,'same original positive a'))",
        ),
    )
    try:
        sigma_receipt = _load(repo / SIGMA_RECEIPT)
    except Exception as exc:
        sigma_receipt = {}
        errors.append("sigma receipt load: " + type(exc).__name__ + ": " + str(exc))
    sigma_identity = _mapping(sigma_receipt.get("accepted_exact_sigma_reflection_and_J1"))
    sigma_input_pin = _mapping(sigma_receipt.get("input_hashes")).get(
        Path(SIGMA_SOURCE).name
    ) == SIGMA_SOURCE_SHA256
    sigma_checks = {
        "reflection_receipt": sigma_identity.get("original_sigma_reflection_identity") is True,
        "J1_half_receipt": sigma_identity.get("original_J1_half_by_reflected_exact_integral")
        is True,
        "sigma_source_hash_in_receipt": sigma_input_pin,
        "sigma_source_ast": sigma_ast.get("passed") is True,
        "passed": all(
            (
                sigma_identity.get("original_sigma_reflection_identity") is True,
                sigma_identity.get("original_J1_half_by_reflected_exact_integral") is True,
                sigma_input_pin,
                sigma_ast.get("passed") is True,
            )
        ),
    }
    try:
        rh_name_report = _load(repo / RH_NAME)
        rh_receipt = _load(repo / RH_RECEIPT)
    except Exception as exc:
        rh_name_report = {}
        rh_receipt = {}
        errors.append("Rh identity receipt load: " + type(exc).__name__ + ": " + str(exc))
    graph_bindings = _mapping(graph.get("true_outer_C6_C5_bindings"))
    graph_rh_bindings = _mapping(
        _mapping(graph_bindings.get("bindings")).get("original_Rh_leading_identity")
    )
    rh_actual_bindings = _mapping(rh_name_report.get("actual_source_bindings"))
    rh_reference_binding = _mapping(rh_actual_bindings.get("reference_background"))
    rh_reference_statements = rh_reference_binding.get("required_defining_statements", ())
    rh_receipt_inputs = _mapping(rh_receipt.get("input_hashes"))
    rh_report_inputs = _mapping(rh_name_report.get("input_hashes"))
    rh_identity = {
        "source_report_all_passed": rh_name_report.get("all_passed") is True,
        "source_receipt_all_passed": rh_receipt.get("all_passed") is True,
        "source_receipt_identity_gate": rh_receipt.get(
            "current_whole_Z_repaired_leading_Rh_reference_function_identity_checked"
        )
        is True,
        "source_receipt_c1_gate": rh_receipt.get(
            "current_whole_Z_Rh_C1_correction_preserving_reference_continuation_installed"
        )
        is True,
        "source_report_family_matches_graph": rh_name_report.get("source_family")
        == graph_family,
        "source_receipt_family_matches_graph": rh_receipt.get("source_family")
        == graph_family,
        "source_report_bindings_digest": _digest(rh_actual_bindings)
        == RH_BINDINGS_DIGEST,
        "graph_report_bindings_equal": graph_rh_bindings == rh_actual_bindings,
        "source_report_input_hash": rh_report_inputs.get(Path(RH_SOURCE).name)
        == RH_SOURCE_SHA256,
        "source_report_name_hash": source_pins["Rh_name"]["passed"],
        "source_receipt_input_hash": rh_receipt_inputs.get(Path(RH_SOURCE).name)
        == RH_SOURCE_SHA256,
        "source_receipt_name_hash": rh_receipt_inputs.get(Path(RH_NAME).name)
        == RH_NAME_SHA256,
        "reference_background_function_identity": rh_reference_binding.get("function")
        == "background_cell",
        "reference_background_module_hash": rh_reference_binding.get("module_sha256")
        == RH_REFERENCE_SOURCE_SHA256,
        "reference_background_AST_hash": rh_reference_binding.get("AST_sha256")
        == RH_REFERENCE_AST_SHA256,
        "reference_background_E_formula_binding": (
            "logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10"
            in rh_reference_statements
        ),
        "reference_background_AST_guard": rh_reference_ast.get("passed") is True,
        "fixed_leading_input_N0_binding": _node(
            graph.get("exact_graph_nodes", ()), 21
        ).get("fixed_leading_input_N0")
        == rh_name_report.get("fixed_leading_input_N0"),
    }
    rh_identity["passed"] = all(rh_identity.values())
    family_match = graph_family == accepted.get("source_family")
    checks = {
        "accepted_parent_sha": parent_pin.get("passed") is True,
        "accepted_parent_mapping": parent_pin.get("mapping_match") is True,
        "accepted_parent_sections": all(accepted_parent_sections.values()),
        "source_family_matches_parent": family_match,
        "outer_graph_hash": source_pins["outer_graph"]["passed"],
        "outer_source_hash": source_pins["outer_source"]["passed"],
        "loop_source_hash": source_pins["prescribed_loop_source"]["passed"],
        "sigma_source_hash": source_pins["sigma_source"]["passed"],
        "sigma_receipt_hash": source_pins["sigma_receipt"]["passed"],
        "outer_ast": outer_ast.get("passed") is True,
        "loop_ast": loop_ast.get("passed") is True,
        "Rh_identity": rh_identity["passed"],
        "Rh_source_bindings_ast": rh_bindings_ast.get("passed") is True,
        "Rh_join_init_ast": rh_join_ast.get("passed") is True,
        "Rh_reference_ast": rh_reference_ast.get("passed") is True,
        "sigma_identity": sigma_checks["passed"],
        "graph": graph_checks.get("passed") is True,
    }
    evaluation: dict[str, Any]
    try:
        ctx = MPIntervalContext()
        ctx.dps = 80
        raw = evaluate_buffer_swirl(ctx, (Fraction(0), Fraction(11)), (Fraction(-1, 2), Fraction(1, 2)))
        evaluation = {
            "passed": True,
            "x_domain_test": [0, 11],
            "Z_domain_test": ["-1/2", "1/2"],
            "log_scale": _record(raw["log_scale"]),
            "scale": _record(raw["scale"]),
            "log_E": _record(raw["log_E"]),
            "E_rows": _record_rows(raw["E_rows"]),
            "a_rows": _record_rows(raw["a_rows"]),
            "b_rows": _record_rows(raw["b_rows"]),
            "t0_rows": _record_rows(raw["t0_rows"]),
            "ordinary_Z_orders": raw["ordinary_Z_orders"],
            "mixed_y_orders": raw["mixed_y_orders"],
            "physical_jacobian": _record(raw["physical_jacobian"]),
            "source_formula": raw["source_formula"],
            "source_executed": False,
        }
    except Exception as exc:
        evaluation = {
            "passed": False,
            "source_executed": False,
            "error": type(exc).__name__ + ": " + str(exc),
        }
        errors.append("fresh interval evaluator: " + evaluation["error"])
    checks["fresh_interval_evaluator"] = evaluation.get("passed") is True
    checks["passed"] = all(checks.values())
    # An error is a failed source certificate, even when it is retained in
    # the receipt for diagnosis.  Keep this fail closed for missing or stale
    # saved inputs.
    checks["errors_recorded"] = not errors
    result = {
        "passed": checks["passed"],
        "checks": checks,
        "source_only": True,
        "source_executed": False,
        "accepted_parent_sha256": PARENT_SHA256,
        "parent_report": PARENT,
        "parent_report_pin": parent_pin,
        "source_family": graph_family,
        "source_pins": source_pins,
        "AST": {
            "outer": outer_ast,
            "prescribed_loop": loop_ast,
            "sigma": sigma_ast,
            "Rh_source_bindings": rh_bindings_ast,
            "Rh_join_init": rh_join_ast,
            "Rh_reference": rh_reference_ast,
        },
        "graph_checks": graph_checks,
        "sigma_identity": sigma_checks,
        "Rh_identity": rh_identity,
        "Rh_reference_formula": {
            "source_module": RH_REFERENCE_SOURCE,
            "formula": "logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10",
            "source_executed": False,
        },
        "accepted_parent_sections": accepted_parent_sections,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "mixed_y_orders": list(MIXED_Y_ORDERS),
        "actual_chart": {
            "name": "O2_buffer",
            "coordinate_node": O2_COORDINATE,
            "domain_node_ids": list(O2_DOMAIN_NODE),
            "domain_values": list(O2_DOMAIN_VALUE),
            "physical_jacobian": 1,
            "phase_node": O2_PHASE_NODE,
            "E_rows": list(O2_E_ROWS),
            "E_y_rows": list(O2_EY_ROWS),
            "a_rows": list(O2_A_ROWS),
            "b_rows": [0] * 6,
            "t0_rows": [0] * 6,
        },
        "source_formula": (
            "E(x,Z)=exp(3/10-exp(40)/2-x/2)/(1+Z^2); "
            "E_y=-E/2; a=2; b=t0=0"
        ),
        "q_dstar": {
            "q_formula": "q=sqrt(eta/2)",
            "eta_node": 21,
            "dstar_node": 22,
            "N_node": 18,
            "eta_and_dstar_numeric_values_evaluated": False,
        },
        "evaluation": evaluation,
        "numeric_C_H_evaluated": False,
        "numeric_p2_mixed_evaluated": False,
        "numeric_q_dstar_evaluated": False,
        "leading_background_swirl_not_modified_E_N": True,
        "cap_contraction_or_improvement": False,
        "archived_owner_builder_checker_fixture_or_quadrature_executed": False,
        "limitations": [
            "This receipt supplies source E/a/t0 interval jets only.",
            "The p2 mixed jet, phase inverse, C, H, and all signed history values remain unevaluated.",
            "The huge common N is retained as source node 18 and is not materialized.",
        ],
        "errors": errors,
        "root": str(repo),
    }
    if evaluation.get("passed") is True:
        result["exact_interval_records_roundtripped"] = _serialization_roundtrip(ctx, evaluation)
        checks["interval_serialization_roundtrip"] = result["exact_interval_records_roundtripped"] > 0
        checks["passed"] = all(checks.values())
        result["passed"] = checks["passed"]
    return result


__all__ = ["build", "evaluate_buffer_swirl"]
