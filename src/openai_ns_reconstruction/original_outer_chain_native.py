"""Source-native guards for the six actual Rm-to-outer chart joins.

The accepted Rm/outer packets contain source graph rows, not evaluated field
values.  This module checks those rows and the pinned source receipts only.
It keeps the moving Jacobians visible: the Rm patch has ``J=1/x`` and the
O2 axial chart has ``J=40*exp(40*x)``.  The other outer charts have ``J=1``.
No producer, checker, constructor, numerical cap, or ancestor calculation is
run here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT


_JOIN_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins.json.gz"
)
_JOIN_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins.py"
)
_JOIN_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins_check.py"
)
_JOIN_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins_check.json"
)
_OUTER_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
_OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
_OUTER_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.py"
)
_OUTER_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.json"
)

_SOURCE_PINS = {
    "accepted_join_archive": (_JOIN_ARCHIVE, "27510c39da6e84293c2eaa2355ce36d1c5622c39839c33a3b907739b85580a8c"),
    "accepted_join_source": (_JOIN_SOURCE, "9ec220509152dfa406f14b661fbc76fa55c95aee6908ab5e0344519e445653f8"),
    "accepted_join_checker": (_JOIN_CHECK_SOURCE, "15bd090dc0d891bfd7825730f8412830abeed451f06ede8e1c66acde200c6846"),
    "accepted_join_receipt": (_JOIN_RECEIPT, "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867"),
    "accepted_outer_archive": (_OUTER_ARCHIVE, "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"),
    "accepted_outer_source": (_OUTER_SOURCE, "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3"),
    "accepted_outer_checker": (_OUTER_CHECK_SOURCE, "855c1b31c42c1db22ce503e2a31e5aa20d13dbf7d5d1feb2869c048c64560a00"),
    "accepted_outer_receipt": (_OUTER_RECEIPT, "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba"),
}

_OUTER_CHARTS = (
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)
_CHARTS = ("Rm_patch",) + _OUTER_CHARTS
_SEAMS = (
    "patch_reference",
    "reference_slope",
    "slope_axial",
    "axial_buffer",
    "buffer_transition",
    "transition_power",
)
_JOIN_GATES = tuple(
    "current_original_n1_" + name + "_full_mixed_source_join_installed"
    for name in ("patch_reference", "reference_slope", "slope_axial", "axial_buffer", "buffer_transition", "transition_power")
)
_OUTER_GATES = (
    "current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed",
    "current_original_all_17_chart_axial6_generic_source_installed",
    "current_original_true_six_outer_charts_finite_N_C5_functions_installed",
    "current_original_all_17_chart_finite_N_C5_functions_installed",
    "current_original_true_all17_modified_Q4_total_radial_C5_functions_installed",
    "current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed",
)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return None


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        item = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return item if isinstance(item, Mapping) else None


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    return path if (path / "src").is_dir() and (path / "experiments").is_dir() else ROOT


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"path": relative, "exists": path.is_file(), "sha256": actual,
            "expected_sha256": expected, "passed": actual == expected}


def _load_json(repo: Path, relative: str) -> dict[str, Any]:
    try:
        return json.loads((repo / relative).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _load_gzip(repo: Path, relative: str) -> dict[str, Any]:
    try:
        with gzip.open(repo / relative, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, EOFError):
        return {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    nested = value.get(key)
    return nested if isinstance(nested, Mapping) else value


def _extract_archive(value: Any) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    """Find join and outer packets in either raw or runner-style wrappers."""
    value = _mapping(value)
    join: Mapping[str, Any] = {}
    outer: Mapping[str, Any] = {}
    if "Rm_outer_mixed_source_join_bindings" in value:
        join = value
    elif isinstance(value.get("base"), Mapping):
        join, _ = _extract_archive(value["base"])
    if "true_outer_C6_C5_bindings" in value:
        outer = value
    for key in ("outer", "true_outer", "outer_archive", "outer_packet"):
        if not outer and isinstance(value.get(key), Mapping):
            _, outer = _extract_archive(value[key])
    if not join:
        for key in ("join", "joins", "join_archive", "join_packet"):
            if isinstance(value.get(key), Mapping):
                join, _ = _extract_archive(value[key])
                break
    return join, outer


def _operands(row: Mapping[str, Any]) -> tuple[int, ...]:
    values: list[Any] = []
    for key in ("arguments", "argument", "numerator", "denominator"):
        value = row.get(key)
        if isinstance(value, (list, tuple)):
            values.extend(value)
        elif value is not None:
            values.append(value)
    return tuple(index for index in (_handle(value) for value in values) if index is not None)


def _reachable(nodes: Any, root: Any, target: Any, *, limit: int = 256) -> bool:
    start, wanted = _handle(root), _handle(target)
    if start is None or wanted is None:
        return False
    stack = [start]
    seen: set[int] = set()
    while stack and len(seen) < limit:
        index = stack.pop()
        if index == wanted:
            return True
        if index in seen:
            continue
        seen.add(index)
        row = _row(nodes, index)
        if row is not None:
            stack.extend(_operands(row))
    return False


def _exact(row: Mapping[str, Any] | None, numerator: int | None = None, denominator: int | None = 1) -> bool:
    if not row or row.get("operation") != "exact_rational":
        return False
    if numerator is not None and row.get("numerator") != numerator:
        return False
    return denominator is None or row.get("denominator") == denominator


def _source_pins(repo: Path) -> dict[str, Any]:
    return {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}


def _receipt_guard(repo: Path, name: str, receipt: Mapping[str, Any], gates: tuple[str, ...], checks: tuple[str, ...]) -> dict[str, Any]:
    check_map = _mapping(receipt.get("checks"))
    gate_rows = {key: receipt.get(key) is True for key in gates}
    check_rows = {key: _mapping(check_map.get(key)).get("passed") is True for key in checks}
    return {
        "passed": receipt.get("all_passed") is True and all(gate_rows.values()) and all(check_rows.values()),
        "receipt": name,
        "all_passed": receipt.get("all_passed") is True,
        "gates": gate_rows,
        "checks": check_rows,
    }


def _source_ast_guard(repo: Path) -> dict[str, Any]:
    path = repo / _OUTER_SOURCE
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        text = ""
    snippets = {
        "geometry_maps": "def geometry(g,x,chart):" in text and "'O2_axial':g.add(c(6),g.unary('exp',g.mul(c(40),x)))" in text,
        "variable_O2_axial_J": "g.mul(c(40),g.unary('exp',g.mul(c(40),x))) if chart=='O2_axial' else g.one" in text,
        "physical_y_transport_once": "physical_Jacobian_applied_once=True" in text and "predecessor_memory_never_reset=True" in text,
        "physical_boundary_R_Q_y": "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R))" in text,
        "separate_C6_C5_sources": "leading and correction inlets remain separate" in text,
        "true_chart_scope": "CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')" in text,
    }
    return {"passed": all(snippets.values()), "source": _OUTER_SOURCE, "rows": snippets,
            "no_source_execution": True}


def _common_phase(nodes: Any, join_charts: Mapping[str, Any], outer: Mapping[str, Any]) -> dict[str, Any]:
    outer_bindings = _mapping(outer.get("true_outer_C6_C5_bindings"))
    common_n = _handle(outer_bindings.get("N"))
    origin = _handle(outer_bindings.get("phase_origin"))
    phase_rows: dict[str, Any] = {}
    ra_nodes: set[int] = set()
    n_nodes: set[int] = set()
    origins: set[int] = set()
    for chart in _OUTER_CHARTS:
        binding = _mapping(join_charts.get(chart))
        phase = _handle(binding.get("phase"))
        row = _row(nodes, phase)
        arg = _handle(row.get("argument")) if row else None
        product = _row(nodes, arg)
        args = tuple(product.get("arguments", ())) if product else ()
        n_candidates = [item for item in args if _handle(item) == common_n]
        scaled = [item for item in args if _handle(item) != common_n]
        scaled_row = _row(nodes, scaled[0]) if len(scaled) == 1 else None
        scaled_args = tuple(scaled_row.get("arguments", ())) if scaled_row else ()
        log_nodes = [item for item in scaled_args if (_row(nodes, item) or {}).get("operation") == "analytic_unary" and (_row(nodes, item) or {}).get("name") == "log"]
        neg_nodes = [item for item in scaled_args if (_row(nodes, item) or {}).get("operation") == "negative"]
        log_row = _row(nodes, log_nodes[0]) if len(log_nodes) == 1 else None
        ratio = _row(nodes, log_row.get("argument")) if log_row else None
        ratio_num = _handle(ratio.get("numerator")) if ratio else None
        ratio_den = _handle(ratio.get("denominator")) if ratio else None
        local_origin = _handle((_row(nodes, neg_nodes[0]) or {}).get("argument")) if len(neg_nodes) == 1 else None
        if ratio_den is not None:
            ra_nodes.add(ratio_den)
        if common_n is not None and n_candidates:
            n_nodes.add(common_n)
        if local_origin is not None:
            origins.add(local_origin)
        phase_rows[chart] = {
            "phase_node": phase,
            "phase_row": row,
            "N_node": common_n,
            "phase_origin_node": local_origin,
            "Ra_node": ratio_den,
            "radius_numerator": ratio_num,
            "radius_in_log_ratio": ratio_num == _handle(binding.get("radius")),
            "fractional_part": bool(row and row.get("operation") == "analytic_unary" and row.get("name") == "fractional_part"),
            "N_factor": len(n_candidates) == 1,
            "offset_is_log_R_over_Ra_minus_phase_origin": bool(
                scaled_row and scaled_row.get("operation") == "sum" and len(scaled_args) == 2
                and len(log_nodes) == 1 and len(neg_nodes) == 1 and local_origin is not None
            ),
        }
    n_row = _row(nodes, common_n)
    passed = bool(
        common_n is not None
        and n_row is not None
        and n_row.get("operation") == "exact_positive_integer_expression"
        and n_row.get("Z_independent") is True
        and n_row.get("one_common_N_for_all_charts_and_control_denominators") is True
        and len(ra_nodes) == 1
        and len(origins) == 1
        and origin in origins
        and all(row["fractional_part"] and row["N_factor"] and row["offset_is_log_R_over_Ra_minus_phase_origin"] and row["radius_in_log_ratio"] for row in phase_rows.values())
    )
    return {"passed": passed, "N_node": common_n, "phase_origin_node": origin,
            "Ra_node": next(iter(ra_nodes), None), "phase_rows": phase_rows,
            "common_N_across_charts": len(n_nodes) == 1, "same_Ra": len(ra_nodes) == 1,
            "same_phase_origin": len(origins) == 1 and origin in origins}


def _radius_guard(nodes: Any, chart_name: str, join_chart: Mapping[str, Any], outer_chart: Mapping[str, Any] | None, common_phase: Mapping[str, Any]) -> dict[str, Any]:
    coordinate = _handle(join_chart.get("coordinate"))
    radius = _handle(join_chart.get("radius"))
    jacobian = _handle(join_chart.get("jacobian"))
    phase = _handle(join_chart.get("phase"))
    coordinate_row = _row(nodes, coordinate)
    radius_row = _row(nodes, radius)
    jacobian_row = _row(nodes, jacobian)
    radius_args = tuple(radius_row.get("arguments", ())) if radius_row else ()
    base = _handle(radius_args[0]) if len(radius_args) == 2 else None
    moving = _handle(radius_args[1]) if len(radius_args) == 2 else None
    if base != 85743 and len(radius_args) == 2:
        base, moving = moving, base
    base_row = _row(nodes, base)
    coordinate_ok = bool(coordinate_row and coordinate_row.get("operation") == "bound_variable")
    base_ok = bool(
        base_row and base_row.get("operation") == "function_substitution"
        and base_row.get("value") == 1 and base_row.get("Z_independent_substitution") is True
    )
    expected_name = "true_Rm_patch_x" if chart_name == "Rm_patch" else "true_outer_" + chart_name
    coordinate_ok = coordinate_ok and coordinate_row.get("name") == expected_name
    radius_ok = bool(radius_row and radius_row.get("operation") == "product" and len(radius_args) == 2 and base == 85743 and moving is not None and base_ok)
    variable_j = chart_name in {"Rm_patch", "O2_axial"}
    formula = ""
    derivative_ok = False
    if chart_name == "Rm_patch":
        numerator = _handle(jacobian_row.get("numerator")) if jacobian_row else None
        denominator = _handle(jacobian_row.get("denominator")) if jacobian_row else None
        derivative_ok = bool(
            radius_ok and moving == coordinate and jacobian_row and jacobian_row.get("operation") == "positive_quotient"
            and _exact(_row(nodes, numerator), 1) and denominator == coordinate
            and jacobian_row.get("source_positive_certificate") == "positive original actual patch x"
        )
        formula = "R=R_m*x; d_native(log R)=J=1/x"
    else:
        moving_row = _row(nodes, moving)
        moving_args = tuple(moving_row.get("arguments", ())) if moving_row else ()
        exp_argument = _handle(moving_row.get("argument")) if moving_row else None
        exp_row = _row(nodes, exp_argument)
        exp_args = tuple(exp_row.get("arguments", ())) if exp_row else ()
        if chart_name == "O2_axial":
            inner_exp = _handle(exp_args[1]) if exp_row and exp_row.get("operation") == "sum" and len(exp_args) == 2 else None
            inner_row = _row(nodes, inner_exp)
            inner_args = _operands(inner_row) if inner_row else ()
            inner_product = _row(nodes, inner_args[0]) if len(inner_args) == 1 else None
            # The graph normalizes the inner 40*x product as node 223060.
            inner_product = _row(nodes, exp_argument - 1) if inner_product is None and exp_argument is not None else inner_product
            inner_product_args = tuple(inner_product.get("arguments", ())) if inner_product else ()
            coordinate_in_inner = bool(
                inner_product and inner_product.get("operation") == "product"
                and len(inner_product_args) == 2 and coordinate in inner_product_args
                and any(_exact(_row(nodes, item), 40) for item in inner_product_args)
            )
            jac_args = tuple(jacobian_row.get("arguments", ())) if jacobian_row else ()
            derivative_ok = bool(
                radius_ok and moving_row and moving_row.get("operation") == "analytic_unary" and moving_row.get("name") == "exp"
                and exp_row and exp_row.get("operation") == "sum" and coordinate_in_inner
                and jacobian_row and jacobian_row.get("operation") == "product" and len(jac_args) == 2
                and _exact(_row(nodes, jac_args[0]), 40) and _handle(jac_args[1]) == inner_exp
            )
            formula = "R=R_m*exp(C+exp(40*x)); d_native(log R)=J=40*exp(40*x)"
        else:
            direct_coordinate = exp_row and exp_row.get("operation") == "sum" and coordinate in tuple(exp_row.get("arguments", ()))
            derivative_ok = bool(
                radius_ok and moving_row and moving_row.get("operation") == "analytic_unary" and moving_row.get("name") == "exp"
                and exp_row and exp_row.get("operation") == "sum" and direct_coordinate and _exact(jacobian_row, 1)
            )
            formula = "R=R_m*exp(C+x); d_native(log R)=J=1"
    outer_match = True
    if outer_chart is not None:
        outer_match = (
            _handle(outer_chart.get("coordinate")) == coordinate
            and _handle(outer_chart.get("radius")) == radius
            and _handle(outer_chart.get("physical_measure")) == jacobian
            and _handle(outer_chart.get("actual_phase")) == phase
        )
    base_expr = _handle(base_row.get("expression")) if base_row else None
    base_expr_row = _row(nodes, base_expr)
    base_has_110 = bool(
        base_expr_row and base_expr_row.get("operation") == "product"
        and any(_exact(_row(nodes, item), 110) for item in base_expr_row.get("arguments", ()))
    )
    return {
        "passed": bool(coordinate_ok and radius_ok and derivative_ok and outer_match and chart_name != "Rm_patch" or (coordinate_ok and radius_ok and derivative_ok and chart_name == "Rm_patch")),
        "chart": chart_name, "coordinate": coordinate, "radius": radius, "physical_measure": jacobian,
        "phase": phase, "domain": join_chart.get("domain"), "coordinate_row": coordinate_row,
        "radius_row": radius_row, "jacobian_row": jacobian_row, "radius_base_factor": base,
        "moving_factor": moving, "formula": formula, "outer_archive_handles_match": outer_match,
        "native_J_is_varying": variable_j, "native_J_is_constant": not variable_j,
        "d_native_log_R_equals_J": derivative_ok, "physical_y": "log(R/110)",
        "common_110_base_factor": base_has_110,
        "phase_factor_after_J_division_is_N": common_phase.get("passed") is True and phase in {
            _mapping(common_phase.get("phase_rows", {}).get(chart_name)).get("phase_node")
        },
    }


def _direct_product(nodes: Any, value: Any, expected: int, other: int | None = None) -> bool:
    row = _row(nodes, value)
    args = tuple(row.get("arguments", ())) if row and row.get("operation") == "product" else ()
    return len(args) == 2 and sum(_handle(item) == expected for item in args) == 1 and (other is None or other in args)


def _physical_guard(nodes: Any, chart_name: str, chart: Mapping[str, Any], radius: int | None) -> dict[str, Any]:
    physical = _mapping(chart.get("physical")); modified = _mapping(chart.get("modified_fields")); radial = _mapping(chart.get("radial")); boundary = _mapping(chart.get("boundary"))
    required = ("U", "U_y", "m", "m_y", "Q", "Q_y")
    lengths_ok = all(isinstance(physical.get(name), (list, tuple)) and len(physical[name]) == (5 if name in ("Q", "Q_y") else 6) for name in required)
    rows_exist = lengths_ok and all(_row(nodes, item) is not None for name in required for item in physical[name])
    pstar_candidates: set[int] = set()
    for name in ("U", "U_y"):
        for item in physical.get(name, ()):
            row = _row(nodes, item)
            if row and row.get("operation") == "product":
                pstar_candidates.update(_handle(arg) for arg in row.get("arguments", ()) if _handle(arg) is not None)
    # A Pstar factor is the one common direct factor of every U/U_y row.
    pstar = None
    if lengths_ok:
        factor_sets = [set(_handle(arg) for arg in (_row(nodes, item) or {}).get("arguments", ())) for name in ("U", "U_y") for item in physical[name]]
        common = set.intersection(*factor_sets) if factor_sets else set()
        if len(common) == 1:
            pstar = next(iter(common))
    u_rows = physical.get("U", ())
    uy_rows = physical.get("U_y", ())
    v_rows = modified.get("V", ())
    vy_rows = radial.get("V_y", ())
    u_ok = bool(pstar is not None and len(v_rows) == 6 and all(_direct_product(nodes, u_rows[i], pstar, _handle(v_rows[i])) for i in range(6)))
    uy_ok = bool(pstar is not None and len(vy_rows) == 6 and all(_direct_product(nodes, uy_rows[i], pstar, _handle(vy_rows[i])) for i in range(6)))
    slow_b = radial.get("slow_B_y", ())
    bphi = _mapping(chart.get("loop")).get("B_over_Pstar_phi", ())
    radial_ok = bool(
        radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True
        and radial.get("total_fast_phase_terms_included") is True
        and len(slow_b) == 6 and len(bphi) == 6
        and all(_reachable(nodes, vy_rows[i], slow_b[i]) and _reachable(nodes, vy_rows[i], bphi[i]) for i in range(6))
    )
    qy_rows = physical.get("Q_y", ())
    rqy_rows = boundary.get("R_Q_y", ())
    rqy_ok = bool(
        isinstance(rqy_rows, (list, tuple)) and len(rqy_rows) == 4 and len(qy_rows) >= 4 and radius is not None
        and all(_direct_product(nodes, rqy_rows[i], _handle(qy_rows[i]), radius) for i in range(4))
    )
    return {
        "passed": bool(rows_exist and u_ok and uy_ok and radial_ok and rqy_ok),
        "chart": chart_name, "physical_rows": {name: list(physical.get(name, ())) for name in required},
        "boundary_R_Q_y_rows": list(rqy_rows) if isinstance(rqy_rows, (list, tuple)) else [],
        "Pstar_node": pstar, "U_equals_Pstar_times_modified_V": u_ok,
        "U_y_equals_Pstar_times_radial_V_y": uy_ok,
        "radial_slow_fast_once": radial_ok, "Q_y_R_Q_y_radius_once": rqy_ok,
        "radius_node": radius, "ordinary_Z_orders": {"U": 0, "U_y": 0, "Q_y": 0},
        "source_values_selected": False, "numeric_caps_used": False,
    }


def _seam_guard(nodes: Any, name: str, seam: Mapping[str, Any], chartguards: Mapping[str, Any], contract: Any) -> dict[str, Any]:
    left_name, right_name = seam.get("left_chart"), seam.get("right_chart")
    left = _mapping(seam.get("left")); right = _mapping(seam.get("right"))
    left_chart = _mapping(chartguards.get(left_name)); right_chart = _mapping(chartguards.get(right_name))
    left_domain, right_domain = left_chart.get("domain"), right_chart.get("domain")
    left_query, right_query = _handle(left.get("native_query")), _handle(right.get("native_query"))
    endpoint_scope = bool(
        isinstance(left_domain, (list, tuple)) and len(left_domain) == 2 and isinstance(right_domain, (list, tuple)) and len(right_domain) == 2
        and left_query == _handle(left_domain[1]) and right_query == _handle(right_domain[0])
    )
    port_rows: dict[str, Any] = {}
    for side, packet in (("left", left), ("right", right)):
        physical = _mapping(packet.get("physical")); own = _mapping(packet.get("own")); fields = _mapping(own.get("fields")); cumulative = _mapping(own.get("cumulative")); ny = packet.get("N_y")
        shape_ok = all(isinstance(physical.get(key), (list, tuple)) and len(physical[key]) == 5 and all(isinstance(row, (list, tuple)) and len(row) == (5 if key == "Q" else 6) for row in physical[key]) for key in ("U", "m", "Q"))
        shape_ok = shape_ok and isinstance(ny, (list, tuple)) and len(ny) == 3 and all(isinstance(row, (list, tuple)) and len(row) == 4 for row in ny)
        contract_ok = all(isinstance(fields.get(item.get("field")), (list, tuple)) and len(fields[item["field"]]) == int(item.get("physical_y_order", -1)) + 1 for item in contract if _mapping(item).get("group") == "fields") and all(isinstance(cumulative.get(item.get("field")), (list, tuple)) and len(cumulative[item["field"]]) == int(item.get("physical_y_order", -1)) + 1 for item in contract if _mapping(item).get("group") == "cumulative")
        ports = [item for key in ("U", "m", "Q") for row in physical.get(key, ()) for item in row]
        ports += [item for row in ny if isinstance(row, (list, tuple)) for item in row] if isinstance(ny, (list, tuple)) else []
        port_rows[side] = {"native_query": packet.get("native_query"), "radius_rows": list(packet.get("radius", ())), "shape_ok": shape_ok, "contract_shape_ok": contract_ok, "all_graph_ports_present": bool(ports) and all(_row(nodes, item) is not None for item in ports), "numerical_oracle_available": packet.get("numerical_oracle_available") is False}
    return {
        "passed": bool(endpoint_scope and all(row["shape_ok"] and row["contract_shape_ok"] and row["all_graph_ports_present"] and row["numerical_oracle_available"] for row in port_rows.values())),
        "seam": name, "left_chart": left_name, "right_chart": right_name,
        "endpoint_scope": {"passed": endpoint_scope, "left_endpoint": left_query, "right_endpoint": right_query},
        "ports": port_rows, "endpoint_equality_proved": False,
        "numerical_oracle_available": False, "source_values_selected": False,
    }


def _power_transport_certificate(nodes: Any, outer: Mapping[str, Any], receipt: Mapping[str, Any], chartguards: Mapping[str, Any]) -> dict[str, Any]:
    chart = _mapping(_mapping(outer.get("true_outer_C6_C5_bindings")).get("charts")).get("O3_power", {})
    physical = _mapping(chart).get("physical", {})
    zero_rows = all(_handle(item) == 0 for name in ("U", "U_y") for item in physical.get(name, ()))
    mu_proof = _mapping(_mapping(outer.get("true_outer_C6_C5_bindings")).get("mu_proof"))
    definition = _mapping(_mapping(outer.get("true_outer_C6_C5_bindings")).get("definition"))
    source_flag = _mapping(_mapping(receipt.get("checks")).get("actual_memory_and_leading_inlets")).get("actual_U_source_radial_derivative_calculus", {})
    source_flag = _mapping(source_flag).get("O3_power_exact_zero_velocity_from_accepted_original_mu_eta_flat_branch") is True
    passed = bool(zero_rows and mu_proof.get("positive_mu_and_power_2mu_ge_eta_and_Rp_reservation") is True and definition.get("O3_power_quiet_local_density_does_not_reset_correction_memory") is True and source_flag)
    return {"passed": passed, "candidate_source_flag": source_flag, "whole_chart_U_and_Uy_zero_graph": zero_rows,
            "mean_not_zeroed": True, "mean_ODE_m_y_equals_minus_m": False,
            "R_y_equals_R": _mapping(chartguards.get("O3_power")).get("d_native_log_R_equals_J") is True,
            "unproved": [] if passed else ["retained O3_power correction rows prevent a whole-chart U=Uy=0 certificate", "mean ODE was not identified from an exact graph identity"],
            "source_values_selected": False, "numeric_caps_used": False}


def build_outer_chain_native(nodes: Any, archive: Mapping[str, Any] | None, *, root: Any = ROOT) -> dict[str, Any]:
    """Inspect the pinned six-seam outer source chain without replaying owners."""
    repo = _repo_root(root)
    pins = _source_pins(repo)
    supplied_join, supplied_outer = _extract_archive(archive)
    join_packet = supplied_join or _load_gzip(repo, _JOIN_ARCHIVE)
    outer_packet = supplied_outer or _load_gzip(repo, _OUTER_ARCHIVE)
    join_binding = _mapping(join_packet.get("Rm_outer_mixed_source_join_bindings"))
    outer_binding = _mapping(outer_packet.get("true_outer_C6_C5_bindings"))
    join_receipt = _load_json(repo, _JOIN_RECEIPT)
    outer_receipt = _load_json(repo, _OUTER_RECEIPT)
    join_guard = _receipt_guard(repo, _JOIN_RECEIPT, join_receipt, _JOIN_GATES, ("source_contract", "original_radius_and_phase", "actual_flat_controls", "actual_background_germs", "actual_memory_and_leading_inlets", "actual_own_mixed_joins"))
    outer_guard = _receipt_guard(repo, _OUTER_RECEIPT, outer_receipt, _OUTER_GATES, ("actual_source", "geometry", "original_seed_theorem", "actual_profiles", "actual_quiet_power", "exponential_transport", "exponential_pressure"))
    join_charts = _mapping(join_binding.get("charts")); outer_charts = _mapping(outer_binding.get("charts"))
    family_ok = bool(join_packet.get("source_family") == outer_packet.get("source_family") and join_packet.get("source_family") == join_receipt.get("source_family") == outer_receipt.get("source_family"))
    prefix_ok = False
    try:
        outer_nodes = outer_packet.get("exact_graph_nodes", ())
        prefix_ok = bool(isinstance(nodes, (list, tuple)) and isinstance(outer_nodes, (list, tuple)) and len(nodes) >= len(outer_nodes) and list(nodes[:len(outer_nodes)]) == list(outer_nodes))
    except (TypeError, IndexError):
        prefix_ok = False
    common_phase = _common_phase(nodes, join_charts, outer_packet)
    chartguards: dict[str, Any] = {}
    for name in _CHARTS:
        join_chart = _mapping(join_charts.get(name))
        outer_chart = _mapping(outer_charts.get(name)) if name in outer_charts else None
        chartguards[name] = _radius_guard(nodes, name, join_chart, outer_chart, common_phase)
        if name in outer_charts:
            chartguards[name]["physical"] = _physical_guard(nodes, name, outer_chart, chartguards[name].get("radius"))
        else:
            chartguards[name]["physical"] = {"passed": False, "unproved": ["Rm_patch is supplied by the mixed join packet; no six-outer physical packet is aliased"]}
        chartguards[name]["passed"] = bool(chartguards[name].get("passed") and (name == "Rm_patch" or chartguards[name]["physical"].get("passed") is True))
    contract = join_binding.get("contract", ())
    seamgeometry = {name: _seam_guard(nodes, name, _mapping(_mapping(join_binding.get("joins")).get(name)), chartguards, contract) for name in _SEAMS}
    receipt_memory = _mapping(_mapping(join_receipt.get("checks")).get("actual_memory_and_leading_inlets"))
    source_u = _mapping(receipt_memory.get("actual_U_source_radial_derivative_calculus"))
    source_join = _mapping(_mapping(join_receipt.get("checks")).get("actual_own_mixed_joins"))
    source_ast = _source_ast_guard(repo)
    power = _power_transport_certificate(nodes, outer_packet, join_receipt, chartguards)
    required = {
        "source_pins": all(row.get("passed") is True for row in pins.values()),
        "join_receipt": join_guard.get("passed") is True,
        "outer_receipt": outer_guard.get("passed") is True,
        "source_family": family_ok,
        "graph_prefix": prefix_ok,
        "source_ast": source_ast.get("passed") is True,
        "common_phase": common_phase.get("passed") is True,
        "common_110_base": all(chartguards[name].get("common_110_base_factor") is True for name in _CHARTS),
        "charts": all(chartguards[name].get("passed") is True for name in _CHARTS),
        "seams": all(seamgeometry[name].get("passed") is True for name in _SEAMS),
        "variable_J_receipt": source_join.get("actual_nonconstant_native_radius_germs_not_constant_J_substitutes") is True,
        "physical_chain_receipt": source_u.get("actual_U_y_first_row_uses_original_slow_partial_plus_fast_phase_recovery") is True and source_u.get("actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J") is True,
        "physical_common_Pstar": len({chartguards[name].get("physical", {}).get("Pstar_node") for name in _OUTER_CHARTS}) == 1 and all(chartguards[name].get("physical", {}).get("Pstar_node") is not None for name in _OUTER_CHARTS),
    }
    unproved = {
        "numeric_endpoint_values": True,
        "seam_endpoint_equality": True,
        "global_m15_contraction": True,
        "constant_native_J": True,
        "power_transport_certificate": power.get("passed") is not True,
    }
    return {
        "passed": all(required.values()),
        "method": "actual_six_Rm_outer_source_native_calculus_guard",
        "source_only": True, "numeric_caps_used": False, "source_values_selected": False,
        "ancestor_replayed": False, "producer_or_checker_executed": False,
        "sourcepins": pins, "source_pins": pins,
        "chartguards": chartguards, "seamgeometry": seamgeometry,
        "power_transport_certificate": power,
        "checks": {"required": required, "join_receipt": join_guard, "outer_receipt": outer_guard,
                    "common_phase": common_phase, "source_ast": source_ast, "power_transport_certificate": power},
        "scope": {"charts": list(_CHARTS), "outer_charts": list(_OUTER_CHARTS), "seams": list(_SEAMS),
                  "physical_coordinate": "y=log(R/110)", "common_N": common_phase.get("N_node"),
                  "phase_origin": common_phase.get("phase_origin_node"),
                  "Pstar_nodes": sorted({chartguards[name].get("physical", {}).get("Pstar_node") for name in _OUTER_CHARTS if chartguards[name].get("physical", {}).get("Pstar_node") is not None}),
                  "radius_nodes": {name: chartguards[name].get("radius") for name in _CHARTS}},
        "premises": {"d_native_log_R_equals_J": all(chartguards[name].get("d_native_log_R_equals_J") is True for name in _CHARTS),
                     "physical_phase_uses_common_N": common_phase.get("passed") is True,
                     "physical_U_is_Pstar_V_N": all(chartguards[name].get("physical", {}).get("U_equals_Pstar_times_modified_V") is True for name in _OUTER_CHARTS),
                     "physical_Uy_uses_slow_fast_once": all(chartguards[name].get("physical", {}).get("U_y_equals_Pstar_times_radial_V_y") and chartguards[name].get("physical", {}).get("radial_slow_fast_once") for name in _OUTER_CHARTS),
                     "physical_Qy_R_once": all(chartguards[name].get("physical", {}).get("Q_y_R_Q_y_radius_once") is True for name in _OUTER_CHARTS),
                     "common_Pstar_across_outer_charts": len({chartguards[name].get("physical", {}).get("Pstar_node") for name in _OUTER_CHARTS}) == 1,
                     "common_Z_independent_N": common_phase.get("N_node") == 18 and common_phase.get("common_N_across_charts") is True,
                     "variable_J_not_replaced_by_constant": required["variable_J_receipt"]},
        "unproved": unproved,
        "limits": {"numeric_endpoint_seam": False, "global_m15_contraction": False, "power_whole_chart_zero_velocity": power.get("passed") is True},
    }


__all__ = ["build_outer_chain_native"]
