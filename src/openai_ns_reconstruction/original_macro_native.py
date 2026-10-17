"""Archive-only native derivative guards for the macro endpoint joins.

This module binds the actual serialized C5 native derivative graph for the
three charts needed by the macro endpoint source proof.  It deliberately
does not instantiate an owner, run a producer/checker, evaluate a field, or
interpret a numerical enclosure as a source value.

The high radial archive contains the complete graph prefix used by the
micro--macro source packet and the macro--first-switch packet.  The ordinary
loop archive is consumed only for the shared N and phase-origin metadata;
the selected chart bindings must be byte-for-byte equal to the high packet's
loop and modified bindings.  This keeps the physical conversion

    d_y = d_native / J

explicit for each chart, including the distinct macro J = Y - 2 hb source
expression.  Native J nodes are never identified across charts.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_native_derivative_guards import (
    _ast_bundle,
    _contains,
    _handle,
    _int_row,
    _mapping,
    _native_expr,
    _phase_binding,
    _poly_serial,
    _slow_b_source_binding,
)
from .original_reference_germ_guards import _function_evidence
from .original_upstream_native import _common_N, _pressure_port


_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_PREFIX = "experiments/root_st073/"
_CHARTS = ("second_micro", "frozen_macro", "first_switch")
_ORDINARY_Z = (0, 1)
_FIELD_KEYS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")

_EXPECTED = {
    "second_micro": {
        "coordinate": 66988,
        "domain": [1, 48],
        "radius": 68644,
        "native_J": 8,
        "global_offset": 97299,
        "actual_phase": 97301,
        "J_source_name": "hbB",
    },
    "frozen_macro": {
        "coordinate": 70265,
        "domain": [0, 1],
        "radius": 71448,
        "native_J": 66987,
        "global_offset": 101907,
        "actual_phase": 101909,
        "J_source_name": None,
    },
    "first_switch": {
        "coordinate": 73190,
        "domain": [0, 1],
        "radius": 85905,
        "native_J": 9,
        "global_offset": 104535,
        "actual_phase": 104537,
        "J_source_name": "hbS",
    },
}

_ARCHIVE_KEYS = {
    "loop": "prescribed_loop_C5_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
    "high": "prescribed_high_radial_source_bindings",
}

_DEFAULTS = {
    "loop": _PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz",
    "modified": _PREFIX + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz",
    "high": _PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.json.gz",
    "micro": _PREFIX + "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.json.gz",
    "switch": _PREFIX + "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins.json.gz",
    "parent": _PREFIX + "lei_ren_part1_paper_compliant_current_original_switch_m15_numeric.json",
}

_PARENT_SHA256 = "b4575dd765e43f9f02aac18b5fa11e3825380bd42f6a30aefd5b420bc86ba0e5"

# These pins are source provenance, never numerical value selections.  The
# parent report is checked separately because its bytes are updated by the
# current owner as the source stage advances.
_PINS = {
    "loop_archive": (_DEFAULTS["loop"], "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "loop_source": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py", "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"),
    "loop_receipt": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json", "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"),
    "modified_archive": (_DEFAULTS["modified"], "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "modified_source": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py", "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
    "modified_receipt": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json", "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"),
    "high_archive": (_DEFAULTS["high"], "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991"),
    "high_source": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.py", "a66e8afedc6c8294f3393b46e84bd19e29a0eaa238539dae1032fa0d8acd6e41"),
    "high_receipt": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check.json", "e53c1fb2a485d5456ee4d1f9484df261aa5cff0f27ea8d1d3a80fa951a7c34ce"),
    "micro_archive": (_DEFAULTS["micro"], "0c673c0979a69f0a6f6fd48e2ee72adae2223f6e384d537d6e39801c15af79e2"),
    "micro_source": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.py", "3f0bb8e0363990ab4d851020e988dd463061086d3ad25e7f0e59c047fa3c5bc0"),
    "micro_receipt": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check.json", "49e7d4f0c566abcbfe9833c5cc3c63875b1c80fa6e6025601eab8aaa83956225"),
    "switch_archive": (_DEFAULTS["switch"], "ddbd506d2fbfc4c58f7a5ba93e0abe3a95d3af8fd3c520ee3ce5024d90a29567"),
    "switch_source": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins.py", "0ae7f651c0d183756f0e6ab90d56dfd80bf3d3c1a6717ffdc2a03fe08dd52f2f"),
    "switch_receipt": (_PREFIX + "lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins_check.json", "9cf3c6d941df95f4b5f73e7b97e32198a54fcd3d17a85ecca2a0742c0d602fc3"),
}

_RECEIPT_CHECKS = {
    "micro": (
        "source_contract",
        "actual_background_recipes",
        "independent_physical_transport_seam",
        "actual_smooth_source_seam",
        "actual_leading_recovery_and_common_N",
        "actual_macro_pressure_and_FTC",
        "actual_own_mixed_seam",
    ),
    "high": (
        "source_contract",
        "original_geometry",
        "actual_recovery_and_pressure_sources",
        "actual_signed_transport_sources",
        "actual_supported_own_flow_and_queries",
        "actual_macro_first_background_seam",
        "actual_macro_first_own_seam",
    ),
}


def _mapping_or_empty(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "experiments").is_dir() and (path / "src").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _nodes(value: Any) -> Any:
    candidate = getattr(value, "nodes", None)
    return candidate if isinstance(candidate, (list, tuple, Mapping)) else value


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def _load_packet(value: Any, repo: Path, relative: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"path": "<supplied mapping>", "passed": True, "supplied_mapping": True}
    path = Path(value) if value is not None else repo / relative
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {"path": str(path), "passed": False, "error": str(exc), "supplied_mapping": False}
    return (data if isinstance(data, Mapping) else {}), {
        "path": str(path), "passed": isinstance(data, Mapping), "supplied_mapping": False,
    }


def _extract(value: Any, key: str, depth: int = 0) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or depth > 4:
        return {}
    child = value.get(key)
    if isinstance(child, Mapping):
        return child
    for item in value.values():
        if isinstance(item, Mapping):
            found = _extract(item, key, depth + 1)
            if found:
                return found
    return {}


def _binding(packet: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    found = _extract(packet, key)
    if found:
        return found
    # Direct bindings are accepted only when their distinguishing structure is
    # present; a random mapping must not become a source alias.
    if key == _ARCHIVE_KEYS["loop"] and "charts" in packet and "N" in packet:
        return packet
    if key == _ARCHIVE_KEYS["modified"] and "charts" in packet and "endpoints" in packet:
        return packet
    if key == _ARCHIVE_KEYS["high"] and "loop" in packet and "modified" in packet and "backgrounds" in packet:
        return packet
    return {}


def _pin(repo: Path, key: str) -> dict[str, Any]:
    relative, expected = _PINS[key]
    path = repo / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        return {"key": key, "path": relative, "expected_sha256": expected, "actual_sha256": None,
                "passed": False, "error": str(exc)}
    return {"key": key, "path": relative, "expected_sha256": expected,
            "actual_sha256": actual, "passed": actual == expected}


def _load_receipt(repo: Path, key: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    relative, _ = _PINS[f"{key}_receipt"]
    path = repo / relative
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, {"path": relative, "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {"path": relative, "passed": isinstance(data, Mapping)}


def _receipt_check(value: Any) -> bool:
    if value is True:
        return True
    if isinstance(value, Mapping):
        return value.get("passed") is True
    return False


def _receipt_guard(repo: Path, key: str, receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _mapping_or_empty(receipt.get("checks"))
    rows = {name: _receipt_check(checks.get(name)) for name in _RECEIPT_CHECKS[key]}
    family = receipt.get("source_family")
    passed = bool(receipt.get("all_passed") is True and family == _FAMILY and all(rows.values()))
    return {"passed": passed, "all_passed": receipt.get("all_passed") is True,
            "source_family": family, "checks": rows, "path": _PINS[f"{key}_receipt"][0]}


def _prefix(nodes: Any, packet: Mapping[str, Any]) -> dict[str, Any]:
    current = _nodes(nodes)
    archived = packet.get("exact_graph_nodes")
    if not isinstance(current, (list, tuple)) or not isinstance(archived, list):
        return {"passed": False, "current_node_count": len(current) if isinstance(current, (list, tuple)) else 0,
                "archived_node_count": len(archived) if isinstance(archived, list) else 0}
    passed = len(current) >= len(archived) and list(current[:len(archived)]) == archived
    return {"passed": bool(passed), "current_node_count": len(current), "archived_node_count": len(archived),
            "exact_prefix": bool(passed)}


def _family(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    value = packet.get("source_family")
    return value if isinstance(value, Mapping) else {}


def _source_family_ok(*packets: Mapping[str, Any]) -> bool:
    return bool(packets) and all(_family(packet) == _FAMILY for packet in packets)


def _same_selected_charts(left: Mapping[str, Any], right: Mapping[str, Any], key: str) -> dict[str, Any]:
    rows = {}
    for name in _CHARTS:
        a = _mapping_or_empty(left.get(name))
        b = _mapping_or_empty(right.get(name))
        rows[name] = {
            "passed": bool(a and b and a == b),
            "left_present": bool(a), "right_present": bool(b),
            "key": key,
        }
    return {"passed": bool(rows) and all(row["passed"] for row in rows.values()), "charts": rows}


def _parent_source_row(parent: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    native = _mapping_or_empty(parent.get("native_source_guards"))
    charts = _mapping_or_empty(native.get("chart_guards"))
    return _mapping_or_empty(charts.get(name))


def _canonical_parent_guard(repo: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Require the pinned parent report and its selected native rows.

    A caller-supplied mapping is only an input view; it cannot establish the
    parent native-J premise by itself.  The canonical saved report must still
    have the accepted bytes and the two rows used by this three-chart proof
    must match it exactly.
    """
    path = repo / _DEFAULTS["parent"]
    try:
        actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        canonical = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {"passed": False, "path": _DEFAULTS["parent"], "error": str(exc)}
    canonical_rows = {
        name: _parent_source_row(canonical, name)
        for name in ("frozen_macro", "first_switch")
    }
    supplied_rows = {
        name: _parent_source_row(parent, name)
        for name in ("frozen_macro", "first_switch")
    }
    checks = {
        "sha256": actual_sha == _PARENT_SHA256,
        "report_checks": _mapping_or_empty(canonical.get("checks")).get("passed") is True,
        "native_guards": _mapping_or_empty(canonical.get("native_source_guards")).get("passed") is True,
        "selected_rows": canonical_rows == supplied_rows,
    }
    return {"passed": bool(all(checks.values())), "path": _DEFAULTS["parent"],
            "expected_sha256": _PARENT_SHA256, "actual_sha256": actual_sha,
            "checks": checks, "selected_rows": canonical_rows,
            "source_values_selected": False}


def _positive_jacobian(
    nodes: Any,
    name: str,
    chart: Mapping[str, Any],
    parent: Mapping[str, Any],
    micro_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    handle = _handle(chart.get("physical_measure"))
    row = _int_row(nodes, handle)
    expected = _EXPECTED[name]
    parent_row = _parent_source_row(parent, name)
    parent_bound = bool(
        parent_row
        and parent_row.get("passed") is True
        and parent_row.get("native_J_positive") is True
        and parent_row.get("positive_native_J_source_premise") is True
        and parent_row.get("native_J") == handle
        and parent_row.get("native_J_row") == row
    )
    if name == "second_micro":
        source_bound = bool(
            row
            and row.get("operation") == "current_original_source_parameter"
            and row.get("name") == expected["J_source_name"]
            and row.get("Z_independent") is True
            and "positive" in str(row.get("definition", "")).lower()
            and _receipt_check(_mapping_or_empty(micro_receipt.get("checks")).get("source_contract"))
        )
    else:
        source_bound = parent_bound
    passed = bool(handle == expected["native_J"] and row and (source_bound or parent_bound))
    return {
        "passed": passed,
        "handle": handle,
        "node": row,
        "native_J_is_constant_in_chart": True,
        "source_positive_premise": "micro positive source ledger" if name == "second_micro" else "canonical parent native source guard",
        "parent_binding": parent_bound,
        "source_binding": source_bound,
        "not_cover_width_alias": name == "frozen_macro" and not (row and row.get("operation") == "current_original_source_parameter"),
        "source_values_selected": False,
    }


def _radius_phase_guard(
    nodes: Any,
    chart: Mapping[str, Any],
    common_n: Any,
    phase_origin: Any,
) -> dict[str, Any]:
    coordinate = _handle(chart.get("coordinate"))
    radius = _handle(chart.get("radius"))
    jacobian = _handle(chart.get("physical_measure"))
    radius_row = _int_row(nodes, radius)
    coordinate_row = _int_row(nodes, coordinate)
    jacobian_row = _int_row(nodes, jacobian)
    if not (coordinate is not None and radius is not None and jacobian is not None and radius_row and coordinate_row and jacobian_row):
        return {"passed": False, "reason": "missing coordinate, radius, or physical measure"}
    radius_args = tuple(_handle(item) for item in radius_row.get("arguments", ()))
    exp_nodes = [item for item in radius_args if (_int_row(nodes, item) or {}).get("operation") == "analytic_unary"
                 and (_int_row(nodes, item) or {}).get("name") == "exp"]
    if radius_row.get("operation") != "product" or len(radius_args) != 2 or len(exp_nodes) != 1:
        return {"passed": False, "reason": "radius must have one exponential factor", "radius": radius}
    exp_node = exp_nodes[0]
    exp_row = _int_row(nodes, exp_node)
    exponent = _handle(exp_row.get("argument")) if exp_row else None
    memo: dict[int, Any] = {}
    jac_expr = _native_expr(nodes, jacobian, coordinate, memo)
    exponent_expr = _native_expr(nodes, exponent, coordinate, memo)
    factor_rows = []
    for item in radius_args:
        row = _int_row(nodes, item)
        expr = _native_expr(nodes, item, coordinate, memo)
        factor_rows.append({"handle": item, "operation": row.get("operation") if row else None,
                            "native_derivative_zero": bool(expr is not None and not expr[1])})
    other_factors_independent = all(row["native_derivative_zero"] for row in factor_rows if row["handle"] != exp_node)
    log_derivative_equals_J = bool(jac_expr is not None and exponent_expr is not None and exponent_expr[1] == jac_expr[0])
    native_J_constant = bool(jac_expr is not None and not jac_expr[1])

    offset = _handle(chart.get("global_offset"))
    offset_row = _int_row(nodes, offset)
    offset_terms = tuple(_handle(item) for item in (offset_row or {}).get("arguments", ()))
    log_term = None
    origin_term = None
    radius_ratio = None
    for item in offset_terms:
        item_row = _int_row(nodes, item)
        if item_row and item_row.get("operation") == "analytic_unary" and item_row.get("name") == "log":
            log_term = item
            ratio = _int_row(nodes, item_row.get("argument"))
            if ratio and ratio.get("operation") == "positive_quotient" and _handle(ratio.get("numerator")) == radius:
                radius_ratio = _handle(ratio.get("denominator"))
        elif item_row and item_row.get("operation") == "negative" and _handle(item_row.get("argument")) == _handle(phase_origin):
            origin_term = item
    ratio_expr = _native_expr(nodes, radius_ratio, coordinate, memo) if radius_ratio is not None else None
    origin_expr = _native_expr(nodes, phase_origin, coordinate, memo) if phase_origin is not None else None
    offset_ok = bool(
        offset_row and offset_row.get("operation") == "sum" and len(offset_terms) == 2
        and log_term is not None and radius_ratio is not None and ratio_expr is not None and not ratio_expr[1]
        and origin_term is not None and origin_expr is not None and not origin_expr[1]
    )
    phase = _phase_binding(nodes, chart, common_n)
    phase_ok = phase.get("passed") is True
    passed = bool(other_factors_independent and native_J_constant and log_derivative_equals_J and offset_ok and phase_ok)
    return {
        "passed": passed,
        "coordinate": coordinate,
        "radius": radius,
        "native_J": jacobian,
        "coordinate_row": coordinate_row,
        "radius_row": radius_row,
        "native_J_row": jacobian_row,
        "radius_exponential_node": exp_node,
        "radius_exponent_node": exponent,
        "radius_factor_rows": factor_rows,
        "other_radius_factors_native_independent": other_factors_independent,
        "canonical_J_expression": _poly_serial(jac_expr[0]) if jac_expr is not None else None,
        "canonical_log_R_native_derivative": _poly_serial(exponent_expr[1]) if exponent_expr is not None else None,
        "d_native_log_R_equals_J": log_derivative_equals_J,
        "native_J_is_constant": native_J_constant,
        "global_offset": offset,
        "global_offset_log_ratio": log_term,
        "global_offset_radius_denominator": radius_ratio,
        "phase_origin": _handle(phase_origin),
        "global_offset_exact_log_ratio_minus_phase_origin": offset_ok,
        "phase": phase,
        "phase_factor_after_J_division_is_N": bool(phase_ok and log_derivative_equals_J),
        "physical_y_equals_native_divided_by_J": bool(log_derivative_equals_J),
        "source_values_selected": False,
    }


def _inverse_N_row(nodes: Any, values: Any, common_n: int) -> tuple[int | None, dict[str, Any]]:
    rows = values if isinstance(values, (list, tuple)) else ()
    if not rows:
        return None, {"passed": False, "reason": "missing A/N source rows"}
    first = _int_row(nodes, rows[0])
    args = tuple(_handle(item) for item in (first or {}).get("arguments", ()))
    candidates = []
    for item in args:
        row = _int_row(nodes, item)
        if row and row.get("operation") == "positive_quotient" and row.get("numerator") == 1 and _handle(row.get("denominator")) == common_n:
            candidates.append(item)
    certificate = (_int_row(nodes, candidates[0]) or {}).get("source_positive_certificate") if candidates else None
    passed = bool(first and first.get("operation") == "product" and len(candidates) == 1
                  and certificate in {
                      "same exact positive integer N on every source chart",
                      "same original exact positive integer N on every source chart",
                  })
    return (candidates[0] if candidates else None), {"passed": passed, "A_over_N_row": _handle(rows[0]),
                                                        "inverse_N_node": candidates[0] if candidates else None,
                                                        "common_N": common_n}


def _radial_Bphi_guard(
    nodes: Any,
    name: str,
    modified: Mapping[str, Any],
    loop_chart: Mapping[str, Any],
    background: Mapping[str, Any],
    common_n: int,
) -> dict[str, Any]:
    radial = _mapping_or_empty(modified.get("radial"))
    loop = _mapping_or_empty(loop_chart.get("loop"))
    inverse_n, n_guard = _inverse_N_row(nodes, _mapping_or_empty(loop_chart.get("finite_N")).get("A_over_N"), common_n)
    vy = radial.get("V_y")
    slow = radial.get("slow_B_y")
    bphi = loop.get("B_over_Pstar_phi")
    bg_vy = _mapping_or_empty(background.get("radial")).get("V_y")
    rows = []
    for order in range(6):
        vy_i = _handle(vy[order]) if isinstance(vy, (list, tuple)) and len(vy) > order else None
        slow_i = _handle(slow[order]) if isinstance(slow, (list, tuple)) and len(slow) > order else None
        bphi_i = _handle(bphi[order]) if isinstance(bphi, (list, tuple)) and len(bphi) > order else None
        bg_i = _handle(bg_vy[order]) if isinstance(bg_vy, (list, tuple)) and len(bg_vy) > order else None
        vy_row = _int_row(nodes, vy_i)
        args = tuple(_handle(item) for item in (vy_row or {}).get("arguments", ()))
        slow_terms = []
        for item in args:
            item_row = _int_row(nodes, item)
            item_args = tuple(_handle(x) for x in (item_row or {}).get("arguments", ()))
            if item_row and item_row.get("operation") == "product" and slow_i in item_args:
                # The source graph has a chart-local positive 1/N node in
                # this radial branch.  It is equivalent to the A/N node but
                # is intentionally not aliased by handle.  Bind its own
                # quotient and certificate instead of requiring graph-node
                # reuse across independently built source expressions.
                quotients = []
                for factor in item_args:
                    factor_row = _int_row(nodes, factor)
                    cert = str((factor_row or {}).get("source_positive_certificate", ""))
                    if (factor_row and factor_row.get("operation") == "positive_quotient"
                            and factor_row.get("numerator") == 1
                            and _handle(factor_row.get("denominator")) == common_n
                            and cert in {
                                "same exact positive common N",
                                "same original exact positive common N",
                                "same exact positive integer N on every source chart",
                                "same original exact positive integer N on every source chart",
                            }):
                        quotients.append(factor)
                if len(quotients) == 1 and len(item_args) == 2:
                    slow_terms.append({"product": item, "inverse_N": quotients[0]})
        rows.append({
            "ordinary_Z_order": order,
            "V_y_node": vy_i,
            "background_V_y_node": bg_i,
            "slow_B_y_node": slow_i,
            "B_over_Pstar_phi_node": bphi_i,
            "operation": vy_row.get("operation") if vy_row else None,
            "background_term_once": bg_i is not None and args.count(bg_i) == 1,
            "slow_B_over_N_term_once": len(slow_terms) == 1,
            "slow_B_over_N_node": slow_terms[0]["inverse_N"] if len(slow_terms) == 1 else None,
            "Bphi_term_once": bphi_i is not None and args.count(bphi_i) == 1,
            "Bphi_has_no_extra_Pstar": bool(bphi_i is not None and _int_row(nodes, bphi_i)
                                              and _int_row(nodes, bphi_i).get("operation") == "original_lazy_flat_branch"),
        })
    ast_rows = {key: bool(value) for key, value in {
        "total_fast_phase_terms_included": radial.get("total_fast_phase_terms_included") is True,
        "physical_y_equals_native_partial_divided_by_Jacobian": radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
    }.items()}
    passed = bool(n_guard.get("passed") and rows and all(
        row["operation"] == "sum" and row["background_term_once"] and row["slow_B_over_N_term_once"]
        and row["Bphi_term_once"] and row["Bphi_has_no_extra_Pstar"] for row in rows
    ) and all(ast_rows.values()))
    return {"passed": passed, "chart": name, "inverse_N": n_guard, "rows": rows,
            "ast_radial_formula": ast_rows, "one_fast_phase_addition": passed,
            "source_values_selected": False}


def _request_guard(nodes: Any, name: str, chart: Mapping[str, Any], modified: Mapping[str, Any], requests: Any) -> dict[str, Any]:
    reqs = requests if isinstance(requests, list) else []
    physical = _mapping_or_empty(modified.get("physical"))
    uy = physical.get("U_y") if isinstance(physical.get("U_y"), (list, tuple)) else []
    expected_orders = [2] * 6 + [3] * 6 + [4] * 6
    J = _handle(chart.get("physical_measure"))
    rows = []
    chain_ok = True
    for index, request in enumerate(reqs):
        if not isinstance(request, Mapping):
            rows.append({"index": index, "passed": False, "reason": "request is not a mapping"})
            chain_ok = False
            continue
        input_i = _handle(request.get("input")); native_i = _handle(request.get("native_output")); physical_i = _handle(request.get("physical_output"))
        physical_row = _int_row(nodes, physical_i)
        args = tuple(_handle(item) for item in (physical_row or {}).get("arguments", ()))
        inverse_rows = []
        for item in args:
            item_row = _int_row(nodes, item)
            if item_row and item_row.get("operation") == "positive_quotient" and item_row.get("numerator") == 1 and _handle(item_row.get("denominator")) == J:
                inverse_rows.append(item)
        if index < 6:
            expected_input = _handle(uy[index]) if index < len(uy) else None
        elif index < 12:
            expected_input = _handle(reqs[index - 6].get("physical_output")) if isinstance(reqs[index - 6], Mapping) else None
        else:
            expected_input = _handle(reqs[index - 6].get("physical_output")) if isinstance(reqs[index - 6], Mapping) else None
        order = request.get("physical_y_order")
        passed = bool(
            order == (expected_orders[index] if index < len(expected_orders) else None)
            and input_i == expected_input
            and _int_row(nodes, input_i) is not None
            and _int_row(nodes, native_i) is not None
            and physical_row and physical_row.get("operation") == "product"
            and native_i in args and len(args) == 2 and len(inverse_rows) == 1
            and (_int_row(nodes, inverse_rows[0]) or {}).get("source_positive_certificate", "").startswith("same positive original")
        )
        rows.append({"index": index, "physical_y_order": order, "input": input_i,
                     "native_output": native_i, "physical_output": physical_i,
                     "inverse_J_node": inverse_rows[0] if inverse_rows else None, "passed": passed})
        chain_ok = chain_ok and passed
    passed = bool(len(reqs) == 18 and chain_ok and all(row.get("passed") is True for row in rows))
    return {"passed": passed, "chart": name, "request_count": len(reqs), "rows": rows,
            "orders": [row.get("physical_y_order") for row in rows],
            "one_native_derivative_then_one_inverse_J": passed, "source_values_selected": False}


def _seed_product_index(nodes: Any) -> dict[tuple[int, int], tuple[int, ...]]:
    """Index the serialized two-factor products used by high-chart seeding."""
    index: dict[tuple[int, int], list[int]] = {}
    iterable = enumerate(nodes) if isinstance(nodes, (list, tuple)) else (
        ((int(key), value) for key, value in nodes.items()) if isinstance(nodes, Mapping) else ()
    )
    for handle, row in iterable:
        if not isinstance(row, Mapping) or row.get("operation") != "product":
            continue
        args = row.get("arguments")
        if not isinstance(args, (list, tuple)) or len(args) != 2:
            continue
        left, right = _handle(args[0]), _handle(args[1])
        if left is None or right is None:
            continue
        index.setdefault(tuple(sorted((left, right))), []).append(handle)
    return {key: tuple(value) for key, value in index.items()}


def _high_seed_recipe_guard(repo: Path) -> dict[str, Any]:
    """Bind the serialized high producer's actual build_chart seed recipe."""
    relative, expected_sha = _PINS["high_source"]
    path = repo / relative
    statements = (
        "coordinate=self.ref(loop['coordinate'])",
        "radius=self.ref(loop['radius'])",
        "J=self.ref(loop['physical_measure'])",
        "phase=self.ref(loop['actual_phase'])",
        "self.background_y[name]={key:self.jet(ids) for key,ids in old['radial_background'].items()}",
        "D=previous.previous.PeriodicTotalNativeDerivative(g,coordinate,phase,coordinate,J,self.N,self.active.guard_core)",
        "self.factories[name]=D",
        "for node,dy in zip(ids,outer.rows(self.background_y[name][key])):D.seed(self.ref(node),g.mul(J,dy))",
        "physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}",
        "refs.append(g.mul(invJ,native))",
    )
    try:
        actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        return {"passed": False, "source": relative, "expected_sha256": expected_sha,
                "actual_sha256": None, "error": str(exc), "statements": []}
    evidence = _function_evidence(path, "CurrentOriginalPrescribedHighRadialSourceFunctions.build_chart", statements)
    return {
        "passed": bool(actual_sha == expected_sha and evidence.get("passed") is True),
        "source": relative,
        "expected_sha256": expected_sha,
        "actual_sha256": actual_sha,
        "scope": evidence.get("scope"),
        "function_ast_sha256": evidence.get("function_ast_sha256", []),
        "statements": evidence.get("statements", []),
        "source_values_selected": False,
    }


def _background_seed_guard(
    nodes: Any,
    modified: Mapping[str, Any],
    background: Mapping[str, Any],
    background_y: Mapping[str, Any],
    high_chart: Mapping[str, Any],
    loop_chart: Mapping[str, Any],
    common_n: int,
    seed_products: Mapping[tuple[int, int], tuple[int, ...]],
    high_seed_recipe: Mapping[str, Any],
) -> dict[str, Any]:
    radial_background = _mapping_or_empty(modified.get("radial_background"))
    fields = _mapping_or_empty(background.get("fields"))
    radial = _mapping_or_empty(background.get("radial"))
    rows = {}
    seed_rows = {}
    J = _handle(loop_chart.get("physical_measure"))
    for key in _FIELD_KEYS:
        left = radial_background.get(key)
        right = background_y.get(key)
        field_values = fields.get(key)
        field_rows = []
        list_shape = bool(
            isinstance(left, (list, tuple)) and isinstance(right, (list, tuple))
            and isinstance(field_values, (list, tuple))
            and len(left) == len(right) == len(field_values)
        )
        if list_shape:
            for order, (field_node, dy_node) in enumerate(zip(field_values, right)):
                field_i, dy_i = _handle(field_node), _handle(dy_node)
                field_row = _int_row(nodes, field_i)
                dy_row = _int_row(nodes, dy_i)
                candidates = list(seed_products.get(tuple(sorted((J, dy_i))), ())) if J is not None and dy_i is not None else []
                dy_zero = bool(
                    dy_row and dy_row.get("operation") == "exact_rational"
                    and dy_row.get("numerator") == 0 and dy_row.get("denominator") == 1
                )
                product_i = candidates[0] if len(candidates) == 1 else (0 if dy_zero else None)
                product_row = _int_row(nodes, product_i)
                product_args = tuple(_handle(item) for item in (product_row or {}).get("arguments", ()))
                product_ok = bool(
                    dy_zero and product_i == 0
                    or (len(candidates) == 1 and product_row
                        and product_row.get("operation") == "product"
                        and product_args == tuple(sorted((J, dy_i))))
                )
                field_rows.append({
                    "ordinary_Z_order": order,
                    "field_node": field_i,
                    "field_row_present": field_row is not None,
                    "background_y_node": dy_i,
                    "background_y_row_present": dy_row is not None,
                    "seed_product_node": product_i,
                    "seed_product_arguments": list(product_args),
                    "seed_product_candidates": candidates,
                    "seed_zero_simplification": dy_zero,
                    "passed": bool(field_row is not None and dy_row is not None and product_ok),
                })
        rows[key] = {"passed": bool(list_shape and left == right and field_rows and all(row["passed"] for row in field_rows)),
                      "modified_radial_background": list(left) if isinstance(left, (list, tuple)) else [],
                      "source_background_y": list(right) if isinstance(right, (list, tuple)) else [],
                      "source_field_nodes": list(field_values) if isinstance(field_values, (list, tuple)) else [],
                      "seed_rows": field_rows}
        seed_rows[key] = field_rows
    vy = radial.get("V_y")
    ey = radial.get("E_y")
    by_e = background_y.get("E")
    by_v = background_y.get("V")
    radial_rows = {
        "radial_E_y_source": {"passed": isinstance(ey, (list, tuple)) and list(ey) == list(by_e) if isinstance(by_e, (list, tuple)) else False,
                              "rows": list(ey) if isinstance(ey, (list, tuple)) else []},
        "radial_V_y_source": {"passed": isinstance(vy, (list, tuple)) and len(vy) >= 6 and list(vy) == list(by_v) if isinstance(by_v, (list, tuple)) else False,
                              "rows": list(vy) if isinstance(vy, (list, tuple)) else []},
    }
    rows.update(radial_rows)
    high_physical = _mapping_or_empty(high_chart.get("physical"))
    modified_physical = _mapping_or_empty(modified.get("physical"))
    physical_rows = {}
    for key in ("U", "m", "Q"):
        high_rows = high_physical.get(key)
        high_y_rows = high_physical.get(key)[1] if isinstance(high_physical.get(key), (list, tuple)) and len(high_physical.get(key)) > 1 else None
        modified_rows = modified_physical.get(key)
        modified_y_rows = modified_physical.get(key + "_y")
        physical_rows[key] = {
            "passed": isinstance(high_rows, (list, tuple)) and len(high_rows) > 1
            and list(high_rows[0]) == list(modified_rows) if isinstance(modified_rows, (list, tuple)) else False,
            "high_chart_order0": list(high_rows[0]) if isinstance(high_rows, (list, tuple)) and high_rows else [],
            "modified_order0": list(modified_rows) if isinstance(modified_rows, (list, tuple)) else [],
        }
        physical_rows[key + "_y"] = {
            "passed": isinstance(high_y_rows, (list, tuple))
            and list(high_y_rows) == list(modified_y_rows) if isinstance(modified_y_rows, (list, tuple)) else False,
            "high_chart_order1": list(high_y_rows) if isinstance(high_y_rows, (list, tuple)) else [],
            "modified_order1": list(modified_y_rows) if isinstance(modified_y_rows, (list, tuple)) else [],
        }
    geometry_rows = {
        "coordinate": _handle(high_chart.get("coordinate")) == _handle(loop_chart.get("coordinate")),
        "phase": _handle(high_chart.get("phase")) == _handle(loop_chart.get("actual_phase")),
        "J": _handle(high_chart.get("jacobian")) == J and J is not None,
        "N": common_n == 18,
    }
    rows["physical_port_links"] = {"passed": bool(physical_rows) and all(row["passed"] for row in physical_rows.values()),
                                    "rows": physical_rows}
    rows["geometry_seed_inputs"] = {"passed": all(geometry_rows.values()), "rows": geometry_rows,
                                     "coordinate": _handle(loop_chart.get("coordinate")),
                                     "phase": _handle(loop_chart.get("actual_phase")),
                                     "J": J, "N": common_n}
    rows["high_build_chart_seed_recipe"] = dict(high_seed_recipe)
    passed = bool(
        high_seed_recipe.get("passed") is True
        and rows["physical_port_links"]["passed"]
        and rows["geometry_seed_inputs"]["passed"]
        and rows.get("radial_E_y_source", {}).get("passed") is True
        and rows.get("radial_V_y_source", {}).get("passed") is True
        and all(rows[key].get("passed") is True for key in _FIELD_KEYS)
    )
    return {"passed": passed, "rows": rows, "seed_rows": seed_rows,
            "high_build_chart_seed_recipe": dict(high_seed_recipe),
            "source_values_selected": False}


def _chart_guard(
    nodes: Any,
    name: str,
    loop: Mapping[str, Any],
    high_loop: Mapping[str, Any],
    modified: Mapping[str, Any],
    high_modified: Mapping[str, Any],
    high_chart: Mapping[str, Any],
    background: Mapping[str, Any],
    background_y: Mapping[str, Any],
    requests: Any,
    parent: Mapping[str, Any],
    micro_receipt: Mapping[str, Any],
    common_n: Mapping[str, Any],
    phase_origin: Any,
    pstar: Any,
    ast: Mapping[str, Any],
    seed_products: Mapping[tuple[int, int], tuple[int, ...]],
    high_seed_recipe: Mapping[str, Any],
) -> dict[str, Any]:
    expected = _EXPECTED[name]
    structure = {
        "coordinate": _handle(loop.get("coordinate")) == expected["coordinate"],
        "domain": list(loop.get("domain", ())) == expected["domain"],
        "radius": _handle(loop.get("radius")) == expected["radius"],
        "native_J": _handle(loop.get("physical_measure")) == expected["native_J"],
        "global_offset": _handle(loop.get("global_offset")) == expected["global_offset"],
        "actual_phase": _handle(loop.get("actual_phase")) == expected["actual_phase"],
        "high_loop_exact_binding": bool(high_loop and high_loop == loop),
        "high_modified_exact_binding": bool(high_modified and high_modified == modified),
        "high_chart_physical_present": bool(high_chart and _mapping_or_empty(high_chart.get("physical"))),
    }
    radius_phase = _radius_phase_guard(nodes, loop, common_n.get("node"), phase_origin)
    positive = _positive_jacobian(nodes, name, loop, parent, micro_receipt)
    modified_shape = {
        "U": modified.get("physical", {}).get("U") if isinstance(modified.get("physical"), Mapping) else None,
        "U_y": modified.get("physical", {}).get("U_y") if isinstance(modified.get("physical"), Mapping) else None,
        "m": modified.get("physical", {}).get("m") if isinstance(modified.get("physical"), Mapping) else None,
        "m_y": modified.get("physical", {}).get("m_y") if isinstance(modified.get("physical"), Mapping) else None,
        "Q": modified.get("physical", {}).get("Q") if isinstance(modified.get("physical"), Mapping) else None,
        "Q_y": modified.get("physical", {}).get("Q_y") if isinstance(modified.get("physical"), Mapping) else None,
    }
    shape_ok = all(isinstance(value, (list, tuple)) and len(value) == (5 if key in ("Q", "Q_y") else 6)
                   and all(_int_row(nodes, item) is not None for item in value) for key, value in modified_shape.items())
    pressure = _pressure_port(nodes, name, modified, loop, pstar)
    slow = _slow_b_source_binding(nodes, _mapping_or_empty(modified.get("radial")).get("slow_B_y"),
                                   _mapping_or_empty(loop.get("loop")).get("B_over_Pstar"), loop.get("physical_measure"))
    radial = _radial_Bphi_guard(nodes, name, modified, loop, background, common_n.get("node"))
    requests_guard = _request_guard(nodes, name, loop, modified, requests)
    seed = _background_seed_guard(
        nodes, modified, background, background_y, high_chart, loop,
        common_n.get("node"), seed_products, high_seed_recipe,
    )
    source_ast = {key: _mapping_or_empty(ast.get(key)) for key in (
        "native_init", "native_seed", "native_substitution", "native_derivative", "native_constructor",
        "radial_formula", "roots", "loop_phase", "finite_n", "pressure_formula",
    )}
    ast_ok = bool(source_ast) and all(row.get("passed") is True for row in source_ast.values())
    chart_rows = {
        "structure": structure,
        "radius_phase_J": radius_phase,
        "positive_native_J": positive,
        "modified_shape": {key: list(value) if isinstance(value, (list, tuple)) else [] for key, value in modified_shape.items()},
        "slow_B_native_derivative": slow,
        "radial_Bphi_and_fast_phase": radial,
        "native_requests": requests_guard,
        "background_seed": seed,
        "physical_pressure_ports": pressure,
        "source_ast": {key: {"passed": row.get("passed") is True, "source": row.get("source"),
                              "statements": row.get("statements", [])} for key, row in source_ast.items()},
    }
    passed = bool(
        all(structure.values()) and radius_phase.get("passed") is True and positive.get("passed") is True
        and shape_ok and slow.get("passed") is True and radial.get("passed") is True
        and requests_guard.get("passed") is True and seed.get("passed") is True
        and pressure.get("passed") is True and ast_ok
    )
    chart_rows.update({
        "passed": passed,
        "chart": name,
        "coordinate": expected["coordinate"],
        "radius": expected["radius"],
        "native_J": expected["native_J"],
        "physical_measure": expected["native_J"],
        "actual_phase": expected["actual_phase"],
        "phase": radius_phase.get("phase"),
        "native_J_is_constant": radius_phase.get("native_J_is_constant") is True,
        "native_J_is_varying": False,
        "native_J_positive": positive.get("passed") is True,
        "positive_native_J_source_premise": {
            "passed": positive.get("passed") is True,
            "native_J_positive": positive.get("passed") is True,
            "numeric_value_selected": False,
            "source": positive.get("source_positive_premise"),
        },
        "d_native_equals_J_d_y": radius_phase.get("d_native_log_R_equals_J") is True,
        "d_y_equals_d_native_over_J": radius_phase.get("physical_y_equals_native_divided_by_J") is True,
        "physical_y_equals_native_over_J": radius_phase.get("physical_y_equals_native_divided_by_J") is True,
        "phase_N_once": radius_phase.get("phase_factor_after_J_division_is_N") is True and radial.get("one_fast_phase_addition") is True,
        "physical_Pstar_once": pressure.get("Pstar_once") is True,
        "physical_R_once": pressure.get("R_Q_y_radius_once") is True,
        "ordinary_Z1": pressure.get("ordinary_Z1") is True and radial.get("passed") is True,
        "source_values_selected": False,
    })
    return chart_rows


def build_macro_native(
    nodes: Any,
    loop_archive: Any,
    modified_archive: Any,
    high_archive: Any,
    parent_report: Any,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Return source-only native derivative evidence for macro endpoint joins."""
    repo = _repo_root(root)
    graph = _nodes(nodes)
    loop_packet, loop_load = _load_packet(loop_archive, repo, _DEFAULTS["loop"])
    modified_packet, modified_load = _load_packet(modified_archive, repo, _DEFAULTS["modified"])
    high_packet, high_load = _load_packet(high_archive, repo, _DEFAULTS["high"])
    micro_packet, micro_load = _load_packet(None, repo, _DEFAULTS["micro"])
    switch_packet, switch_load = _load_packet(None, repo, _DEFAULTS["switch"])
    loop = _binding(loop_packet, _ARCHIVE_KEYS["loop"])
    modified = _binding(modified_packet, _ARCHIVE_KEYS["modified"])
    high = _binding(high_packet, _ARCHIVE_KEYS["high"])
    high_loop = _mapping_or_empty(high.get("loop"))
    high_modified = _mapping_or_empty(high.get("modified"))
    high_charts = _mapping_or_empty(high.get("charts"))
    high_backgrounds = _mapping_or_empty(high.get("backgrounds"))
    high_background_y = _mapping_or_empty(high.get("background_y"))
    high_requests = _mapping_or_empty(high.get("requests"))

    parent_info: dict[str, Any]
    if isinstance(parent_report, Mapping):
        parent = parent_report
        parent_info = {"supplied_mapping": True, "passed": True, "sha256": _digest(parent)}
    else:
        parent, parent_load = _load_packet(parent_report, repo, _DEFAULTS["parent"])
        parent_info = dict(parent_load)
        parent_info["sha256"] = hashlib.sha256((repo / _DEFAULTS["parent"]).read_bytes()).hexdigest() if (repo / _DEFAULTS["parent"]).is_file() else None

    parent_guard = _canonical_parent_guard(repo, parent)
    parent_info["canonical_guard"] = parent_guard

    pins = {key: _pin(repo, key) for key in _PINS}
    pin_passed = all(row.get("passed") is True for row in pins.values())
    source_family_passed = _source_family_ok(loop_packet, modified_packet, high_packet, micro_packet, switch_packet)
    micro_binding = _extract(micro_packet, "micro_macro_mixed_source_join_bindings")
    switch_binding = _extract(switch_packet, "switch_mixed_source_join_bindings")
    graph_prefixes = {
        "full_switch": _prefix(graph, switch_packet),
        "high": _prefix(graph, high_packet),
        "micro_in_high": _prefix(high_packet.get("exact_graph_nodes"), micro_packet),
    }
    # `_prefix` accepts either a graph sequence or an object exposing nodes;
    # high_packet itself contains the serialized high graph.
    if isinstance(high_packet.get("exact_graph_nodes"), list):
        graph_prefixes["micro_in_high"] = {
            "passed": bool(len(high_packet["exact_graph_nodes"]) >= len(micro_packet.get("exact_graph_nodes", ()))
                         and high_packet["exact_graph_nodes"][:len(micro_packet.get("exact_graph_nodes", ()))] == micro_packet.get("exact_graph_nodes")),
            "current_node_count": len(high_packet["exact_graph_nodes"]),
            "archived_node_count": len(micro_packet.get("exact_graph_nodes", ())) if isinstance(micro_packet.get("exact_graph_nodes"), list) else 0,
        }
    loop_match = _same_selected_charts(_mapping_or_empty(loop.get("charts")), high_loop, "loop")
    modified_match = _same_selected_charts(_mapping_or_empty(modified.get("charts")), high_modified, "modified")

    micro_receipt, _ = _load_receipt(repo, "micro")
    high_receipt, _ = _load_receipt(repo, "high")
    micro_receipt_guard = _receipt_guard(repo, "micro", micro_receipt)
    high_receipt_guard = _receipt_guard(repo, "high", high_receipt)
    receipt_passed = bool(micro_receipt_guard.get("passed") and high_receipt_guard.get("passed"))

    common = _common_N(graph, loop)
    common_node = common.get("node")
    common_passed = bool(common.get("passed") and _handle(high.get("N")) == common_node == 18)
    phase_origin = loop.get("phase_origin")
    pstar = _handle(_mapping_or_empty(high.get("downstream")).get("Pstar"))
    ast = _ast_bundle(repo)
    seed_products = _seed_product_index(graph)
    high_seed_recipe = _high_seed_recipe_guard(repo)
    chart_guards: dict[str, Any] = {}
    charts: dict[str, Any] = {}
    for name in _CHARTS:
        chart_guards[name] = _chart_guard(
            graph,
            name,
            _mapping_or_empty(_mapping_or_empty(loop.get("charts")).get(name)),
            _mapping_or_empty(high_loop.get(name)),
            _mapping_or_empty(_mapping_or_empty(modified.get("charts")).get(name)),
            _mapping_or_empty(high_modified.get(name)),
            _mapping_or_empty(high_charts.get(name)),
            _mapping_or_empty(high_backgrounds.get(name)),
            _mapping_or_empty(high_background_y.get(name)),
            high_requests.get(name),
            parent,
            micro_receipt,
            {"passed": common_passed, "node": common_node},
            phase_origin,
            pstar,
            ast,
            seed_products,
            high_seed_recipe,
        )
        # Keep the public chart packet compact.  The full graph stays in the
        # immutable archive and every selected handle remains explicit here.
        charts[name] = {
            "passed": chart_guards[name].get("passed") is True,
            "coordinate": chart_guards[name].get("coordinate"),
            "radius": chart_guards[name].get("radius"),
            "native_J": chart_guards[name].get("native_J"),
            "actual_phase": chart_guards[name].get("actual_phase"),
            "modified_physical": chart_guards[name].get("physical_pressure_ports", {}).get("physical_rows", {}),
            "modified_boundary": chart_guards[name].get("physical_pressure_ports", {}).get("finite_boundary_rows", {}),
            "requests": chart_guards[name].get("native_requests", {}).get("rows", []),
            "slow_B": chart_guards[name].get("slow_B_native_derivative", {}).get("rows", []),
            "radial_Bphi": chart_guards[name].get("radial_Bphi_and_fast_phase", {}).get("rows", []),
            "background_seed_sha256": _digest(chart_guards[name].get("background_seed", {})),
            "high_seed_recipe_sha256": _digest(chart_guards[name].get("background_seed", {}).get("high_build_chart_seed_recipe", {})),
            "radius_phase_J_sha256": _digest(chart_guards[name].get("radius_phase_J", {})),
            "positive_native_J_sha256": _digest(chart_guards[name].get("positive_native_J", {})),
            "source_ast_sha256": _digest(chart_guards[name].get("source_ast", {})),
            "source_values_selected": False,
        }

    ordinary_z1 = bool(all(row.get("ordinary_Z1") is True for row in chart_guards.values()))
    physical_r_once = bool(all(row.get("physical_R_once") is True for row in chart_guards.values()))
    phase_N_once = bool(all(row.get("phase_N_once") is True for row in chart_guards.values()))
    chart_passed = bool(all(row.get("passed") is True for row in chart_guards.values()))
    source_contract = _mapping_or_empty(high.get("definition"))
    definition_required = (
        "ten_actual_original_domains_not_relabelled_seam_neighborhoods",
        "every_chart_fresh_total_native_factory_uses_its_actual_N_times_J",
        "full_actual_U_y_differentiation_mean_ODE_and_one_minus_delta_Q",
        "original_signed_pressure_and_nonzero_predecessor_histories_retained",
        "source_coverage_not_global_interface_closure_or_numerical_bounds",
    )
    definition_passed = all(source_contract.get(key) is True for key in definition_required)
    required = {
        "source_pins": pin_passed,
        "source_family": source_family_passed,
        "graph_prefixes": all(row.get("passed") is True for row in graph_prefixes.values()),
        "source_receipts": receipt_passed,
        "loop_modified_bindings": loop_match.get("passed") is True and modified_match.get("passed") is True,
        "common_N": common_passed,
        "parent_native_guards": parent_guard.get("passed") is True,
        "high_definition": definition_passed,
        "charts": chart_passed,
        "ordinary_Z1": ordinary_z1,
        "physical_R_once": physical_r_once,
        "phase_N_once": phase_N_once,
    }
    passed = bool(all(required.values()))
    source_pins = {key: {"path": row.get("path"), "sha256": row.get("actual_sha256"),
                         "expected_sha256": row.get("expected_sha256"),
                         "actual_sha256": row.get("actual_sha256"), "passed": row.get("passed") is True}
                   for key, row in pins.items()}
    output = {
        "passed": passed,
        "method": "archive_only_macro_second_micro_frozen_macro_first_switch_native_derivative_guard",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_family": dict(_FAMILY),
        "source_pins": source_pins,
        "archive_loads": {"loop": loop_load, "modified": modified_load, "high": high_load, "micro": micro_load, "switch": switch_load, "parent": parent_info},
        "archive_prefixes": graph_prefixes,
        "archive_binding_matches": {"loop": loop_match, "modified": modified_match},
        "source_receipts": {"micro": micro_receipt_guard, "high": high_receipt_guard},
        "common_N": {"node": common_node, "row": common.get("row"), "passed": common_passed,
                     "Z_independent": common.get("Z_independent") is True, "one_common_N": common.get("one_common_N") is True},
        "charts": charts,
        "chart_guards": chart_guards,
        "pressure_ports": {name: chart_guards[name].get("physical_pressure_ports", {}) for name in _CHARTS},
        "phase_and_J": {name: chart_guards[name].get("radius_phase_J", {}) for name in _CHARTS},
        "native_derivative_congruence": {name: chart_guards[name].get("native_requests", {}) for name in _CHARTS},
        "source_ast": {key: {"passed": _mapping_or_empty(value).get("passed") is True,
                              "source": _mapping_or_empty(value).get("source"),
                              "statements": _mapping_or_empty(value).get("statements", [])}
                     for key, value in ast.items()},
        "high_seed_recipe": dict(high_seed_recipe),
        "ordinary_Z1": ordinary_z1,
        "physical_R_once": physical_r_once,
        "phase_N_once": phase_N_once,
        "physical_Pstar_once": bool(all(chart_guards[name].get("physical_Pstar_once") is True for name in _CHARTS)),
        "charts_sha256": _digest(charts),
        "checks": {"required": required, "high_definition": {"passed": definition_passed, "rows": {key: source_contract.get(key) is True for key in definition_required}},
                    "loop_modified_bindings": {"loop": loop_match, "modified": modified_match},
                    "receipts": {"micro": micro_receipt_guard, "high": high_receipt_guard},
                    "parent_native_guards": parent_guard},
        "limits": {
            "source_values_selected": False,
            "numeric_caps_used": False,
            "caps_as_values": False,
            "ancestor_constructor_or_producer_replay": False,
            "endpoint_function_identity": False,
            "signed_budget_contraction": False,
        },
    }
    # Digest the complete native certificate before attaching its own receipt,
    # matching the source helper's non-self-referential certificate convention.
    output["certificate_sha256"] = _digest(output)
    return output


__all__ = ["build_macro_native"]
