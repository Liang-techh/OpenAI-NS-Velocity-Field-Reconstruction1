"""Saved-packet numeric contraction for the two macro endpoint joins.

The adapter consumes the accepted switch packet and two source-only endpoint
certificates.  It deliberately does not construct an ancestor, call a
producer, replay a checker, or evaluate a source function.  The raw packet is
copied verbatim and the signed endpoint ledger removes only the three sectors
whose common source identity has been proved:

``-D(second_micro, 2) + D(frozen_macro, 0) - D(frozen_macro, 1)
   + D(first_switch, 0) = 0``.

Here ``D = R * (Qy_modified - Qy_original_leading)``.  The endpoint terms are
kept as provenance after the packet and source checks; no positive endpoint
budget is manufactured for the exact identity.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import PREFIX, ROOT


SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PARENT_REDUCED_COUNT = 186
NEW_REDUCED_COUNT = 183
UNCHANGED_COUNT = 183

PARENT = PREFIX + "compliant_current_original_switch_m15_numeric.json"
PARENT_SHA = "b4575dd765e43f9f02aac18b5fa11e3825380bd42f6a30aefd5b420bc86ba0e5"
MACRO_REPORT = PREFIX + "compliant_current_original_macro_m15_numeric.json"
MACRO_REPORT_SHA = "076900a16aaea965e229d652b311dd1bcbd2fc73645061fe56aa58575acaf54c"
SECOND_REPORT = PREFIX + "compliant_current_original_second_micro_m15_numeric.json"
SECOND_REPORT_SHA = "57155e42c649dcdc3b30bfaee20b13914038aea62c13300b8ea678d80ef63dc3"
FIRST_REPORT = PREFIX + "compliant_current_original_first_switch_m15_numeric.json"
FIRST_REPORT_SHA = "fd427c4ba14a541f3ec042b644c2e9c65b66077c8350f6807a96c957a36bf83c"

# These are the immutable producer sources used only for AST/report
# decomposition.  They are never imported or executed here.
MACRO_SOURCE = "src/openai_ns_reconstruction/original_macro_m15.py"
MACRO_SOURCE_SHA = "5d2c162e6c61bee3d05955004114c9be73b44471555c0638a06b0bfb3d1febf5"
MACRO_ROWS_SOURCE = "src/openai_ns_reconstruction/original_macro_pressure_rows.py"
MACRO_ROWS_SOURCE_SHA = "6a96ec977ae8f4235bdb874ce9d2c7b5acc15b91f55710fef1ec040e025d6753"
MACRO_RUNNER = "scripts/run_original_macro_m15.py"
MACRO_RUNNER_SHA = "e5bc884f72b3f829c0ed5eb4dbd4fdd79159283890c859e17ceeee1999983319"
SIGNED_TRACE_SOURCE = "src/openai_ns_reconstruction/original_micro_viscous_m15.py"
SIGNED_TRACE_SOURCE_SHA = "fe514784962f2d2ab968b66c537d6e907dd5f8c1d30183123e6498ba030a33b7"
# The generic signed trace receipt fixes an orientation convention.  These
# two accepted source files provide the physical m15 sign and the call path
# that transports it to each prescribed chart.  They are read and parsed
# only; this adapter never imports or executes them.
PRESSURE_DENSITY_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py"
)
PRESSURE_DENSITY_SOURCE_SHA = "58f0ffb862704e8d566e1105fd98128bfa8a669854e978f3fb11b5278a286498"
HIGH_RADIAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.py"
)
HIGH_RADIAL_SOURCE_SHA = "a66e8afedc6c8294f3393b46e84bd19e29a0eaa238539dae1032fa0d8acd6e41"

MICRO_TERM = "micro_shared_finite_N_viscous_outer_endpoint_error"
MACRO_TERM = "macro_finite_N_viscous_error"
FIRST_TERM = "first_switch_second_finite_N_viscous_external_endpoints"

SECOND_MICRO = "second_micro"
FROZEN_MACRO = "frozen_macro"
FIRST_SWITCH = "first_switch"

# Graph handles are intentionally kept separate from the endpoint values.
SCOPES = {
    "micro_macro": {
        "left_chart": SECOND_MICRO, "right_chart": FROZEN_MACRO,
        "left_value": 2, "right_value": 0,
        "left_query": 48, "right_query": 0,
        "left_coordinate": 66988, "right_coordinate": 70265,
        "left_radius": 68644, "right_radius": 71448,
        "left_jacobian": 8, "right_jacobian": 66987,
        "left_phase": 97301, "right_phase": 101909,
    },
    "macro_first": {
        "left_chart": FROZEN_MACRO, "right_chart": FIRST_SWITCH,
        "left_value": 1, "right_value": 0,
        "left_query": 1, "right_query": 0,
        "left_coordinate": 70265, "right_coordinate": 73190,
        "left_radius": 71448, "right_radius": 85905,
        "left_jacobian": 66987, "right_jacobian": 9,
        "left_phase": 101909, "right_phase": 104537,
    },
}


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _copy_mapping(value: Any, label: str) -> dict[str, Any]:
    return copy.deepcopy(dict(_mapping(value, label)))


def _find(rows: Any, *, label: str, **criteria: Any) -> Mapping[str, Any]:
    if not isinstance(rows, (list, tuple)):
        raise ValueError(f"Saved {label} rows are missing")
    for row in rows:
        if isinstance(row, Mapping) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise ValueError(f"Saved {label} not found: {criteria}")


def _path(owner: Any, relative: str) -> Path:
    """Find a pinned repository file without relying on the owner's cwd."""
    rel = Path(relative)
    candidates: list[Path] = []
    directory = getattr(owner, "directory", None)
    if directory is not None:
        base = Path(directory)
        candidates.extend((base / relative, base / rel.name, base.parent / relative))
        candidates.append(base.parent.parent / relative)
    candidates.extend((ROOT / relative, ROOT / "experiments" / "root_st073" / rel.name))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0] if candidates else ROOT / relative


def _pin_file(owner: Any, relative: str, expected: str) -> dict[str, Any]:
    path = _path(owner, relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    hashes = _mapping(getattr(owner, "hashes", {}), "owner.hashes")
    code_hashes = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes")
    basename = Path(relative).name
    known = (hashes.get(relative), hashes.get(basename), code_hashes.get(relative), code_hashes.get(basename))
    for value in known:
        if value is not None and value != actual:
            raise ValueError(f"Changed accepted source: {relative}")
    if actual != expected:
        raise ValueError(f"Unexpected accepted source: {relative}")
    return {"path": relative, "sha256": actual, "expected_sha256": expected, "passed": True}


def _load_json(owner: Any, relative: str, expected: str) -> Mapping[str, Any]:
    _pin_file(owner, relative, expected)
    try:
        value = json.loads(_path(owner, relative).read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load pinned report {relative}") from exc
    return _mapping(value, relative)


def _candidate_report(owner: Any, names: tuple[str, ...], pinned: Mapping[str, Any]) -> Mapping[str, Any]:
    for name in names:
        candidate = getattr(owner, name, None)
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list):
            if _canonical(candidate) != _canonical(pinned):
                raise ValueError(f"In-memory {name} report differs from pinned report")
            return candidate
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list):
            if _canonical(report) != _canonical(pinned):
                raise ValueError(f"In-memory {name} report differs from pinned report")
            return report
    return pinned


def _parent_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, PARENT, PARENT_SHA)
    return _candidate_report(owner, ("parent", "parent_report", "switch_report"), pinned)


def _macro_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, MACRO_REPORT, MACRO_REPORT_SHA)
    return _candidate_report(owner, ("macro_report", "macro_m15_report", "endpoint_report"), pinned)


def _second_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, SECOND_REPORT, SECOND_REPORT_SHA)
    return _candidate_report(owner, ("second_micro_report", "second_report", "micro2_report"), pinned)


def _first_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, FIRST_REPORT, FIRST_REPORT_SHA)
    return _candidate_report(owner, ("first_switch_report", "first_report", "first_m15_report"), pinned)


def _term_names(rows: Any, label: str) -> list[str]:
    if not isinstance(rows, list):
        raise ValueError(f"{label} must be a list")
    names = [row.get("name") if isinstance(row, Mapping) else None for row in rows]
    if any(name is None for name in names):
        raise ValueError(f"{label} has an unnamed term")
    return names


def _exact(label: str, left: Any, right: Any) -> None:
    if left != right:
        raise ValueError(f"{label} full dictionary equality changed")


def _family(owner: Any, parent: Mapping[str, Any]) -> Any:
    family = getattr(owner, "family", parent.get("source_family"))
    if family != parent.get("source_family"):
        raise ValueError("Macro source family differs from accepted switch parent")
    return family


def _unwrap(value: Any, *keys: str) -> Mapping[str, Any]:
    result = _mapping(value, "source binding")
    for key in keys:
        nested = result.get(key)
        if isinstance(nested, Mapping):
            return nested
    return result


def _node(nodes: Any, index: int, label: str) -> Mapping[str, Any]:
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes[index]
    except (IndexError, KeyError, TypeError):
        raise ValueError(f"{label} graph node is missing") from None
    return _mapping(row, label)


def _chart_bindings(owner: Any) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]:
    loop = _unwrap(getattr(owner, "l", None), "prescribed_loop_C5_bindings")
    modified = _unwrap(getattr(owner, "m", None), "modified_Q4_pressure_C3_bindings")
    high_value = getattr(owner, "high", None)
    if high_value is None:
        high_value = getattr(owner, "highbindings", None)
    high = _unwrap(high_value, "prescribed_high_radial_source_bindings")
    for label, value in (("loop", loop), ("modified", modified), ("high", high)):
        if not isinstance(value.get("charts"), Mapping):
            raise ValueError(f"Actual {label} chart bindings are required")
    return loop, modified, high


def _graph_binding(owner: Any) -> dict[str, Any]:
    loop, modified, high = _chart_bindings(owner)
    nodes = getattr(owner, "nodes", None)
    if not isinstance(nodes, (list, tuple, Mapping)):
        raise ValueError("Actual macro graph nodes are required")
    ports: dict[str, dict[str, Any]] = {}
    for chart in (SECOND_MICRO, FROZEN_MACRO, FIRST_SWITCH):
        spec = SCOPES["micro_macro"] if chart == SECOND_MICRO else SCOPES["macro_first"] if chart == FIRST_SWITCH else None
        if chart == FROZEN_MACRO:
            expected = {
                "coordinate": 70265, "radius": 71448, "jacobian": 66987,
                "phase": 101909,
            }
        elif chart == SECOND_MICRO:
            expected = {"coordinate": 66988, "radius": 68644, "jacobian": 8, "phase": 97301}
        else:
            expected = {"coordinate": 73190, "radius": 85905, "jacobian": 9, "phase": 104537}
        lchart = _mapping(loop["charts"].get(chart), f"loop {chart} chart")
        # The prescribed-loop archive owns the coordinate/radius/Jacobian;
        # its phase handle is carried by the high-radial chart registry.
        for key, value in expected.items():
            actual_key = "physical_measure" if key == "jacobian" else key
            if key == "phase" and actual_key not in lchart:
                continue
            if lchart.get(actual_key) != value:
                raise ValueError(f"Actual {chart} {actual_key} changed")
        mchart = _mapping(modified["charts"].get(chart), f"modified {chart} chart")
        physical = _mapping(mchart.get("physical"), f"modified {chart} physical")
        qy = physical.get("Q_y")
        finite = _mapping(mchart.get("finite_boundary"), f"modified {chart} finite boundary")
        rqy = finite.get("R_Q_y")
        if not isinstance(qy, list) or len(qy) < 2 or not isinstance(rqy, list) or len(rqy) < 2:
            raise ValueError(f"Actual {chart} Qy/RQy rows are required")
        qy = list(qy[:2]); rqy = list(rqy[:2])
        for order, qnode in enumerate(qy):
            qrow = _node(nodes, qnode, f"{chart} Qy[{order}]")
            if qrow.get("operation") != "positive_quotient":
                raise ValueError(f"{chart} Qy[{order}] is not an archived pressure function")
        for order, rnode in enumerate(rqy):
            rrow = _node(nodes, rnode, f"{chart} RQy[{order}]")
            if rrow.get("operation") != "product" or list(rrow.get("arguments", [])).count(expected["radius"]) != 1:
                raise ValueError(f"{chart} RQy[{order}] must contain one physical radius")
        prow = _node(nodes, expected["phase"], f"{chart} phase")
        if prow.get("operation") != "analytic_unary":
            raise ValueError(f"{chart} phase source changed")
        ports[chart] = {
            "chart": chart, "coordinate_node": expected["coordinate"],
            "radius_node": expected["radius"], "native_Jacobian_node": expected["jacobian"],
            "phase_node": expected["phase"], "Qy_function_nodes": qy,
            "RQy_function_nodes": rqy, "ordinary_Z_orders": [0, 1],
            "physical_R_once": True, "actual_graph_scope": True,
        }
        hchart = high.get("charts", {}).get(chart)
        if isinstance(hchart, Mapping):
            for key, value in expected.items():
                actual_key = "physical_measure" if key == "jacobian" else key
                if hchart.get(actual_key) not in (None, value):
                    raise ValueError(f"High source {chart} {actual_key} changed")
    # The query nodes are exact graph constants, not coordinate node numbers.
    for value in (0, 1, 2):
        qrow = _node(nodes, {0: 0, 1: 1, 2: 48}[value], f"native query {value}")
        if qrow.get("operation") != "exact_rational" or qrow.get("denominator") != 1 or qrow.get("numerator") != value:
            raise ValueError("Native endpoint query graph value changed")
    return {"passed": True, "charts": ports, "query_value_nodes": {"0": 0, "1": 1, "2": 48},
            "source_function_nodes_not_endpoint_values": True, "physical_R_once": True,
            "ordinary_Z_orders": [0, 1], "nodes_are_actual_archive_handles": True}


def _scope(proof: Mapping[str, Any], side: str, expected: Mapping[str, Any], chart: str,
           endpoint: int, query_node: int) -> Mapping[str, Any]:
    value = proof.get(f"{side}_scope")
    if not isinstance(value, Mapping):
        raise ValueError(f"{proof.get('left_chart')}->{proof.get('right_chart')} {side} scope missing")
    exact = {"chart": chart, "coordinate_node": expected[f"{side}_coordinate"],
             "native_query_node": query_node, "radius_node": expected[f"{side}_radius"],
             "physical_measure_node": expected[f"{side}_jacobian"],
             "ordinary_parameter": "Z", "ordinary_Z_orders": [0, 1, 2, 3]}
    for key, wanted in exact.items():
        if key in value and value.get(key) != wanted:
            raise ValueError(f"{chart} {side} scope has wrong {key}")
    if value.get("native_query_value") not in (None, str(endpoint), endpoint):
        raise ValueError(f"{chart} {side} scope has wrong endpoint value")
    # A coordinate graph handle and a native query value are distinct objects.
    if value.get("endpoint_node") not in (None, endpoint):
        raise ValueError(f"{chart} {side} endpoint scope changed")
    return value


def _proofs(seam_proofs: Mapping[str, Any], graph: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    if not isinstance(seam_proofs, Mapping) or seam_proofs.get("passed") is not True:
        raise ValueError("Macro source seam proofs must pass")
    joins = seam_proofs.get("joins")
    if not isinstance(joins, Mapping):
        joins = seam_proofs.get("seams")
    if not isinstance(joins, Mapping):
        raise ValueError("Macro source seam joins are missing")
    out: dict[str, Mapping[str, Any]] = {}
    for name, expected in SCOPES.items():
        proof = joins.get(name)
        if not isinstance(proof, Mapping) or proof.get("passed") is not True:
            raise ValueError(f"{name} source proof must pass")
        if proof.get("left_chart") != expected["left_chart"] or proof.get("right_chart") != expected["right_chart"]:
            raise ValueError(f"{name} chart identities changed")
        if proof.get("left_endpoint") != expected["left_value"] or proof.get("right_endpoint") != expected["right_value"]:
            raise ValueError(f"{name} endpoint values changed")
        if proof.get("left_endpoint_node") != expected["left_query"] or proof.get("right_endpoint_node") != expected["right_query"]:
            raise ValueError(f"{name} endpoint graph handles changed")
        if proof.get("left_coordinate") != expected["left_coordinate"] or proof.get("right_coordinate") != expected["right_coordinate"]:
            raise ValueError(f"{name} coordinate handles changed")
        _scope(proof, "left", expected, expected["left_chart"], expected["left_value"], expected["left_query"])
        _scope(proof, "right", expected, expected["right_chart"], expected["right_value"], expected["right_query"])
        traces = proof.get("trace_rows")
        if not isinstance(traces, list) or len(traces) != 15 or not all(isinstance(row, Mapping) and row.get("passed") is True for row in traces):
            raise ValueError(f"{name} complete 15-row source trace required")
        counts: dict[str, int] = {}
        for row in traces:
            counts[str(row.get("quantity"))] = counts.get(str(row.get("quantity")), 0) + 1
        if counts != {"U": 3, "Uy": 2, "m": 3, "Q": 3, "Qy": 2, "RQy": 2}:
            raise ValueError(f"{name} source trace quantity rows are incomplete")
        modified = proof.get("modified_source_rows")
        leading = proof.get("leading_source_rows")
        if not isinstance(modified, list) or len(modified) != 9 or not all(isinstance(row, Mapping) and row.get("passed") is True for row in modified):
            raise ValueError(f"{name} complete modified source rows required")
        if not isinstance(leading, list) or len(leading) != 8 or not all(isinstance(row, Mapping) and row.get("passed") is True for row in leading):
            raise ValueError(f"{name} independent leading source rows required")
        modified_counts: dict[str, int] = {}
        for row in modified:
            modified_counts[str(row.get("quantity"))] = modified_counts.get(str(row.get("quantity")), 0) + 1
        leading_counts: dict[str, int] = {}
        for row in leading:
            leading_counts[str(row.get("quantity"))] = leading_counts.get(str(row.get("quantity")), 0) + 1
        if modified_counts != {"U": 3, "Uy": 2, "m": 4} or leading_counts != {"V": 3, "m": 3, "Vy": 2}:
            raise ValueError(f"{name} source row decomposition changed")
        algebra = _mapping(proof.get("pressure_algebra_witness"), f"{name} pressure algebra")
        if algebra.get("passed") is not True:
            raise ValueError(f"{name} pressure algebra witness must pass")
        algebra_rows = algebra.get("rows")
        if not isinstance(algebra_rows, list) or len(algebra_rows) != 5 or not all(isinstance(row, Mapping) and row.get("passed") is True for row in algebra_rows):
            raise ValueError(f"{name} pressure Q/Qy rows are incomplete")
        for row in algebra_rows:
            if row.get("reduced_difference_monomials") not in (None, 0):
                raise ValueError(f"{name} pressure algebra has an unreduced difference")
            if row.get("ordinary_derivative_product_rule") not in (None, True):
                raise ValueError(f"{name} ordinary Z product rule is not certified")
        leading_common = _mapping(proof.get("leading_Qy_common"), f"{name} leading Qy source")
        if leading_common.get("passed") is not True or leading_common.get("source_values_selected") is True or leading_common.get("leading_output_alias_assigned") is True:
            raise ValueError(f"{name} independent leading Qy source is not bound")
        if not isinstance(leading_common.get("input_rows"), list) or len(leading_common["input_rows"]) != 8:
            raise ValueError(f"{name} leading Qy input rows are incomplete")
        bindings = proof.get("endpoint_binding_rows")
        if not isinstance(bindings, list) or len(bindings) != 18 or not all(isinstance(row, Mapping) and row.get("passed") is True for row in bindings):
            raise ValueError(f"{name} endpoint function bindings are incomplete")
        if proof.get("endpoint_wrapper_binding_passed") is not True:
            raise ValueError(f"{name} endpoint wrapper binding is missing")
        drows = proof.get("D_equality_rows")
        if not isinstance(drows, list) or len(drows) != 2 or not all(isinstance(row, Mapping) and row.get("passed") is True and row.get("physical_R_once") is True for row in drows):
            raise ValueError(f"{name} physical D rows are incomplete")
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if proof.get(key) is True:
                raise ValueError(f"{name} source proof promoted a value, cap, or ancestor")
        # Bind the actual Qy and RQy graph functions in the rows, when the
        # proof exposes them.  This rejects a proof that merely relabels a
        # foreign endpoint node.
        for row in traces:
            for side, chart in (("left", expected["left_chart"]), ("right", expected["right_chart"])):
                side_row = row.get(side)
                if not isinstance(side_row, Mapping):
                    continue
                node_value = side_row.get("full_chart_function_node")
                if row.get("quantity") == "Qy" and node_value is not None and node_value not in graph["charts"][chart]["Qy_function_nodes"]:
                    raise ValueError(f"{name} Qy row is not the actual chart function")
                if row.get("quantity") == "RQy" and node_value is not None and node_value not in graph["charts"][chart]["RQy_function_nodes"]:
                    raise ValueError(f"{name} RQy row is not the actual chart function")
        out[name] = copy.deepcopy(dict(proof))
    return out


def _pin_ast_composition(owner: Any) -> dict[str, Any]:
    pins = {
        "macro_m15_source": _pin_file(owner, MACRO_SOURCE, MACRO_SOURCE_SHA),
        "macro_pressure_rows_source": _pin_file(owner, MACRO_ROWS_SOURCE, MACRO_ROWS_SOURCE_SHA),
        "macro_m15_runner_source": _pin_file(owner, MACRO_RUNNER, MACRO_RUNNER_SHA),
    }
    path = _path(owner, MACRO_SOURCE)
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    cls = next((node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "OriginalMacroM15"), None)
    if cls is None:
        raise ValueError("OriginalMacroM15 producer class missing from pinned source")
    methods = {node.name: node for node in cls.body if isinstance(node, ast.FunctionDef)}
    for name in ("correction_densities", "finite_errors", "evaluate"):
        if name not in methods:
            raise ValueError(f"OriginalMacroM15.{name} source method missing")
    correction_text = ast.get_source_segment(text, methods["correction_densities"]) or ""
    finite_text = ast.get_source_segment(text, methods["finite_errors"]) or ""
    evaluate_text = ast.get_source_segment(text, methods["evaluate"]) or ""
    # The Q_y error rows are selected in correction_densities(), then finite_errors()
    # combines only the upper/lower endpoint rows for the viscous branch.  Keep
    # that check scoped to the method that actually contains each source fact;
    # requiring errors['Q_y'] in finite_errors() would falsely reject the pinned
    # producer even though the AST composition is exactly the required one.
    required_correction = ("viscous=errors['Q_y']",)
    required_finite = ("for name, power in (('endpoint', 2), ('viscous', 1))", "result['volume']")
    required_eval = ("for name in ('endpoint', 'volume', 'viscous')", "append('macro_finite_N_'+name+'_error'", "errors[name]")
    if (any(fragment not in correction_text for fragment in required_correction)
            or any(fragment not in finite_text for fragment in required_finite)
            or any(fragment not in evaluate_text for fragment in required_eval)):
        raise ValueError("Macro finite-N producer composition does not isolate the viscous endpoint row")
    methods_out = []
    for name in ("correction_densities", "finite_errors", "evaluate"):
        methods_out.append({"method": f"OriginalMacroM15.{name}",
                            "method_ast_sha256": _digest(ast.dump(methods[name], include_attributes=False)),
                            "source_not_executed": True})
    return {"passed": True, "source_pins": pins, "methods": methods_out,
            "finite_viscous_component": {
                "source_method": "OriginalMacroM15.finite_errors -> correction_densities()['Q_y']",
                "term_factory": "OriginalMacroM15.evaluate append('macro_finite_N_'+name+'_error')",
                "component_names": ["macro_finite_N_endpoint_error", "macro_finite_N_volume_error", MACRO_TERM],
                "viscous_is_endpoint_Qy_only": True, "leading_viscous_rows_are_separate": True,
                "source_not_executed": True}}


def _pin_signed_trace_source(owner: Any) -> dict[str, Any]:
    """Pin the immutable lower-minus-upper FTC orientation and D definition.

    The numeric adapter never calls this source.  Its AST is the accepted
    source receipt for the sign of each chart contribution and for the single
    physical ``R`` factor in ``delta(RQy)``.  The positive Qy error packet is
    therefore not being used to infer a signed function.
    """
    pin = _pin_file(owner, SIGNED_TRACE_SOURCE, SIGNED_TRACE_SOURCE_SHA)
    path = _path(owner, SIGNED_TRACE_SOURCE)
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    trace = next((node for node in tree.body
                  if isinstance(node, ast.FunctionDef) and node.name == "contract_trace_chain"), None)
    cls = next((node for node in tree.body
                if isinstance(node, ast.ClassDef) and node.name == "OriginalMicroViscousM15"), None)
    if trace is None or cls is None:
        raise ValueError("Pinned signed trace source methods are missing")
    interface = next((node for node in cls.body
                      if isinstance(node, ast.FunctionDef) and node.name == "prove_interface"), None)
    if interface is None:
        raise ValueError("Pinned OriginalMicroViscousM15.prove_interface is missing")
    trace_text = ast.get_source_segment(text, trace) or ""
    interface_text = ast.get_source_segment(text, interface) or ""
    if ("for side,sign in (('lower',1),('upper',-1))" not in trace_text
            or "definition='delta(RQy)=R*(Qy_modified-Qy_original_leading)" not in interface_text):
        raise ValueError("Pinned signed FTC/RQy source orientation changed")
    return {
        "passed": True,
        "source_pin": pin,
        "methods": [
            {"method": "contract_trace_chain",
             "method_ast_sha256": _digest(ast.dump(trace, include_attributes=False)),
             "source_not_executed": True},
            {"method": "OriginalMicroViscousM15.prove_interface",
             "method_ast_sha256": _digest(ast.dump(interface, include_attributes=False)),
             "source_not_executed": True},
        ],
        "signed_orientation": {"lower": 1, "upper": -1,
                               "source_statement": "contract_trace_chain lower-minus-upper source IDs"},
        "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
        "physical_radius_power": 1,
        "physical_R_once": True,
        "ordinary_Z_orders": [0, 1],
    }


def _pin_physical_density_source(owner: Any, graph: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the actual m15 density sign to the macro chart ports.

    ``contract_trace_chain`` records the lower-minus-upper convention, but it
    does not establish that convention from the physical m15 source.  This
    receipt therefore parses the accepted pulse density/FTC source and the
    accepted high-radial chart builder.  No source function is imported or
    called.  The high builder is also checked against the already scoped
    frozen-macro graph so a foreign Qy function cannot be relabelled as the
    macro endpoint.
    """
    density_pin = _pin_file(owner, PRESSURE_DENSITY_SOURCE, PRESSURE_DENSITY_SOURCE_SHA)
    high_pin = _pin_file(owner, HIGH_RADIAL_SOURCE, HIGH_RADIAL_SOURCE_SHA)
    density_path = _path(owner, PRESSURE_DENSITY_SOURCE)
    high_path = _path(owner, HIGH_RADIAL_SOURCE)
    density_text = density_path.read_text(encoding="utf-8")
    high_text = high_path.read_text(encoding="utf-8")
    density_tree = ast.parse(density_text, filename=str(density_path))
    high_tree = ast.parse(high_text, filename=str(high_path))
    density_methods = {
        node.name: node for node in density_tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    high_class = next((node for node in high_tree.body
                       if isinstance(node, ast.ClassDef)
                       and node.name == "CurrentOriginalPrescribedHighRadialSourceFunctions"), None)
    if high_class is None:
        raise ValueError("Pinned high-radial source class is missing")
    high_methods = {
        node.name: node for node in high_class.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    density = density_methods.get("pressure_density")
    ftc = density_methods.get("pressure_FTC")
    build_chart = high_methods.get("build_chart")
    if density is None or ftc is None or build_chart is None:
        raise ValueError("Pinned physical density, FTC, and high build_chart methods are required")
    density_source = ast.get_source_segment(density_text, density) or ""
    ftc_source = ast.get_source_segment(density_text, ftc) or ""
    build_source = ast.get_source_segment(high_text, build_chart) or ""
    compact = lambda value: "".join(value.split())
    density_compact = compact(density_source)
    ftc_compact = compact(ftc_source)
    build_compact = compact(build_source)
    required_density = (
        "Py=a.add(",
        "a.trim(fields['Q_y']),a.trim(fields['Q_yy']))",
        "N1p=a.scale(Py,g.quotient(g.one,R,'positiveactualpulseradius'))",
        "m15_y=a.scale(Py,g.neg(R))",
    )
    required_ftc = ("R_q_y=a.scale(a.add(Qy,Qyy),R)",)
    required_build = (
        "radius=self.ref(loop['radius'])",
        "density=rc_source.radial.pressure_density(g,f,self.parameters['delta'],radius)",
        "Ny=[density['N1p'],higher['N1p_y'],higher['N1p_yy']];n0,n1,n2=Ny",
        "negative_R2=g.neg(g.mul(radius,radius))",
        "m15=[b.scale(n0,negative_R2)",
        "raw_m15_y=[b.scale(density['P1_y'],g.neg(radius)),higher['m15_yy'],higher['m15_yyy']]",
    )
    if any(fragment not in density_compact for fragment in required_density):
        raise ValueError("Pinned pressure_density physical m15 sign changed")
    if any(fragment not in ftc_compact for fragment in required_ftc):
        raise ValueError("Pinned pressure_FTC RQy identity changed")
    if any(fragment not in build_compact for fragment in required_build):
        raise ValueError("Pinned high-radial m15 density call path changed")

    macro = _mapping(graph.get("charts", {}).get(FROZEN_MACRO), "frozen_macro graph binding")
    expected_geometry = {
        "coordinate_node": 70265, "radius_node": 71448,
        "native_Jacobian_node": 66987, "phase_node": 101909,
    }
    for key, expected in expected_geometry.items():
        if macro.get(key) != expected:
            raise ValueError(f"Frozen macro physical {key} changed")
    qy_nodes = list(macro.get("Qy_function_nodes", []))
    rqy_nodes = list(macro.get("RQy_function_nodes", []))
    if qy_nodes != [140868, 140872] or rqy_nodes != [141027, 141028]:
        raise ValueError("Frozen macro physical Qy/RQy graph ports changed")
    _, _, high = _chart_bindings(owner)
    high_charts = _mapping(high.get("charts"), "high-radial chart bindings")
    high_macro = _mapping(high_charts.get(FROZEN_MACRO), "high frozen_macro chart")
    for key, expected in (("coordinate", 70265), ("radius", 71448),
                          ("jacobian", 66987), ("phase", 101909)):
        # The high archive carries the actual derivative phase handle; the
        # chart coordinate remains the separate loop coordinate node.
        if key in high_macro and high_macro.get(key) != expected:
            raise ValueError(f"High frozen_macro {key} geometry changed")
    high_physical = _mapping(high_macro.get("physical"), "high frozen_macro physical ports")
    high_q = high_physical.get("Q")
    if (not isinstance(high_q, list) or len(high_q) < 2
            or list(high_q[1][:2]) != qy_nodes):
        raise ValueError("High-radial frozen_macro Qy ports are not the actual macro functions")

    methods = [
        {"method": "pressure_density",
         "method_ast_sha256": _digest(ast.dump(density, include_attributes=False)),
         "source_not_executed": True},
        {"method": "pressure_FTC",
         "method_ast_sha256": _digest(ast.dump(ftc, include_attributes=False)),
         "source_not_executed": True},
        {"method": "CurrentOriginalPrescribedHighRadialSourceFunctions.build_chart",
         "method_ast_sha256": _digest(ast.dump(build_chart, include_attributes=False)),
         "source_not_executed": True},
    ]
    return {
        "passed": True,
        "source_pins": {
            "pulse_radial_pressure_C3_source": density_pin,
            "prescribed_high_radial_source": high_pin,
        },
        "methods": methods,
        "source_not_executed": True,
        "physical_density_sign": {
            "pressure_density_terms": "Py=...+Q_y+Q_yy",
            "pressure_density_N1p": "N1p=Py/R",
            "pressure_density_m15_y": "m15_y=-R*Py",
            "pressure_FTC_R_q_y": "R_q_y=R*(Qy+Qyy)",
            "viscous_identity": "m15_y(viscous)=-R*(Q_y+Q_yy)=-d_y(R*Q_y)",
            "lower_minus_upper": True,
            "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
            "physical_R_once": True,
            "ordinary_Z_orders": [0, 1],
            "source_function_not_cap_value": True,
        },
        "high_radial_call_path": {
            "builder": "CurrentOriginalPrescribedHighRadialSourceFunctions.build_chart",
            "pressure_density_call": "rc_source.radial.pressure_density(g,f,self.parameters['delta'],radius)",
            "density_P1_y_to_raw_m15_y": "raw_m15_y[0]=b.scale(density['P1_y'],g.neg(radius))",
            "N1p_to_m15": "n0=density['N1p']; m15[0]=b.scale(n0,-radius**2)",
            "radius_source_expression": "radius=self.ref(loop['radius'])",
            "source_R_y_equals_R": True,
            "source_not_executed": True,
        },
        "macro_chart_binding": {
            "chart": FROZEN_MACRO,
            "coordinate_node": expected_geometry["coordinate_node"],
            "radius_node": expected_geometry["radius_node"],
            "native_Jacobian_node": expected_geometry["native_Jacobian_node"],
            "phase_node": expected_geometry["phase_node"],
            "Qy_function_nodes": qy_nodes,
            "RQy_function_nodes": rqy_nodes,
            "high_radial_Qy_function_nodes": list(high_q[1][:2]),
            "physical_R_once": True,
            "R_y_equals_R": True,
            "actual_source_ports": True,
        },
    }


def _macro_decomposition(macro: Mapping[str, Any], ast_binding: Mapping[str, Any], Z: str) -> dict[str, Any]:
    row0 = _find(macro.get("samples"), label="macro decomposition sample", Z=Z,
                 coordinate={"kind": "exact_decimal_macro_fraction", "value": ".25"})
    terms = row0.get("m15_C0_C1_terms")
    if not isinstance(terms, list) or len(terms) != 28 or terms[27].get("name") != MACRO_TERM:
        raise ValueError("Macro raw term 27 is not the finite-N viscous sector")
    term = _mapping(terms[27], "macro raw term 27")
    binding = _mapping(term.get("actual_function_binding"), "macro raw term 27 binding")
    for key in ("both_nonflat_pressure_endpoints_retained", "complete_first_second_macro_history",
                "quadratic_cross_errors_retained", "ordinary_C0_C1_bounds_not_derivatives_of_caps"):
        if binding.get(key) is not True:
            raise ValueError("Macro raw term 27 lacks its endpoint-only source binding")
    logs = term.get("log_scale_terms")
    if not isinstance(logs, list) or len(logs) != 1:
        raise ValueError("Macro raw term 27 must have one combined directed endpoint atom")
    atom = _mapping(logs[0], "macro raw term 27 atom")
    if atom.get("positive_error_bound_atom_only") is not True or atom.get("source_value_selected") is not False:
        raise ValueError("Macro raw term 27 atom is not an error-only enclosure")
    if row0.get("finite_N_error_logs_in_complete_physical_units", {}).get("viscous") is None:
        raise ValueError("Macro raw term 27 viscous physical rows are missing")
    return {
        "passed": True, "raw_index": 27, "name": MACRO_TERM,
        "producer_ast": copy.deepcopy(dict(ast_binding)),
        "report_path": MACRO_REPORT, "report_sha256": MACRO_REPORT_SHA,
        "report_sample": {"Z": Z, "coordinate": copy.deepcopy(row0.get("coordinate"))},
        "component_source": "finite_N_error_logs_in_complete_physical_units['viscous']",
        "component_terms": [MACRO_TERM],
        "producer_branch_siblings": ["macro_finite_N_endpoint_error", "macro_finite_N_volume_error"],
        "leading_viscous_indices_retained": [23, 24],
        "raw27_contains_only_viscous_endpoint_components": True,
        "source_function_not_cap_value": True, "directed_atom_only": True,
    }


def _packet_witness(parent_sample: Mapping[str, Any], macro: Mapping[str, Any], second: Mapping[str, Any],
                    first: Mapping[str, Any], family: Any) -> dict[str, Any]:
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Parent raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PARENT_REDUCED_COUNT:
        raise ValueError("Parent switch 186-term view required")
    _term_names(raw, "parent raw terms"); _term_names(previous, "parent reduced terms")
    if raw[27].get("name") != MACRO_TERM:
        raise ValueError("Parent raw macro viscous index 27 changed")
    for name in (MICRO_TERM, MACRO_TERM, FIRST_TERM):
        if sum(term.get("name") == name for term in previous) != 1:
            raise ValueError(f"Parent reduced sector {name} must occur exactly once")
    parent_terms = {name: next(term for term in previous if term.get("name") == name)
                    for name in (MICRO_TERM, MACRO_TERM, FIRST_TERM)}
    micro_binding = _mapping(parent_terms[MICRO_TERM].get("actual_function_binding"),
                             "parent micro signed endpoint binding")
    if (micro_binding.get("function") != "-delta(RQy)(second_micro s=2)"
            or micro_binding.get("ordinary_Z1") is not True
            or micro_binding.get("error_enclosure_only") is not True):
        raise ValueError("Parent micro endpoint sign or physical scale changed")
    first_binding = _mapping(parent_terms[FIRST_TERM].get("actual_function_binding"),
                             "parent first signed endpoint binding")
    if (first_binding.get("function") != "D(first_switch native0)"
            or first_binding.get("delta_definition") != "R*(Qy_modified-Qy_leading)"
            or first_binding.get("ordinary_Z1") is not True
            or first_binding.get("physical_R_once") is not True
            or first_binding.get("error_enclosure_only") is not True):
        raise ValueError("Parent first endpoint sign or physical scale changed")
    selected = _mapping(parent_sample.get("selected_N_evidence"), "parent selected N evidence")
    Z = str(parent_sample.get("Z"))
    macro0 = _find(macro.get("samples"), label="macro native0", Z=Z,
                   coordinate={"kind": "exact_decimal_macro_fraction", "value": "0"})
    macro1 = _find(macro.get("samples"), label="macro native1", Z=Z,
                   coordinate={"kind": "exact_decimal_macro_fraction", "value": "1"})
    second2 = _find(second.get("samples"), label="second micro native2", Z=Z,
                    coordinate={"kind": "exact_decimal_phase", "value": "2"})
    first0 = _find(first.get("samples"), label="first switch native0", Z=Z,
                   coordinate={"kind": "exact_decimal_first_switch_phase", "value": "0"})
    for label, row in (("macro0", macro0), ("macro1", macro1), ("second2", second2), ("first0", first0)):
        if row.get("source_family") != family or row.get("actual_selected_integer") != selected:
            raise ValueError(f"{label} source family or selected N changed")
    macro0_terms = macro0.get("m15_C0_C1_terms"); macro0_in = macro0.get("incoming_second_micro_m15_terms")
    macro1_terms = macro1.get("m15_C0_C1_terms"); macro1_in = macro1.get("incoming_second_micro_m15_terms")
    second_terms = second2.get("m15_C0_C1_terms"); first_in = first0.get("incoming_macro_m15_terms")
    if not isinstance(macro0_terms, list) or len(macro0_terms) != 28 or not isinstance(macro0_in, list) or len(macro0_in) != 17:
        raise ValueError("Macro native0 packet counts changed")
    if not isinstance(macro1_terms, list) or len(macro1_terms) != 28 or not isinstance(macro1_in, list) or len(macro1_in) != 17:
        raise ValueError("Macro native1 packet counts changed")
    if not isinstance(second_terms, list) or len(second_terms) != 17 or not isinstance(first_in, list) or len(first_in) != 28:
        raise ValueError("Adjacent endpoint packet counts changed")
    _exact("macro0 incoming versus second2 terms", macro0_in, second_terms)
    _exact("macro0 incoming versus parent raw prefix17", macro0_in, raw[:17])
    _exact("macro0 term prefix versus incoming", macro0_terms[:17], macro0_in)
    _exact("macro1 incoming versus second2 terms", macro1_in, second_terms)
    _exact("macro1 incoming versus parent raw prefix17", macro1_in, raw[:17])
    _exact("macro1 term prefix versus incoming", macro1_terms[:17], macro1_in)
    _exact("macro1 complete versus first0 incoming", macro1_terms, first_in)
    _exact("first0 incoming versus parent raw prefix28", first_in, raw[:28])
    if macro0.get("exact_empty_increment") is not True or macro1.get("exact_empty_increment") is not False:
        raise ValueError("Macro endpoint emptiness changed")
    if second2.get("exact_empty_increment") is not False or first0.get("exact_empty_increment") is not True:
        raise ValueError("Adjacent endpoint emptiness changed")
    return {
        "passed": True, "Z": Z, "selected_integer": copy.deepcopy(dict(selected)),
        "macro0": macro0, "macro1": macro1, "second2": second2, "first0": first0,
        "full_dictionary_equality": True,
        "raw_prefix_equals_second2_incoming": True, "raw_prefix_equals_macro0_incoming": True,
        "raw_prefix_equals_macro1_incoming": True, "raw_prefix_equals_first0_incoming": True,
        "macro0_terms_prefix_equals_incoming": True, "macro1_terms_prefix_equals_incoming": True,
        "macro1_complete_equals_first0_incoming": True,
        "parent_signed_endpoint_bindings": {
            "micro": copy.deepcopy(dict(micro_binding)),
            "macro": copy.deepcopy(dict(parent_terms[MACRO_TERM].get("actual_function_binding", {}))),
            "first": copy.deepcopy(dict(first_binding)),
        },
        "second2_packet_prefix_count": 17, "macro0_packet_prefix_count": 17,
        "macro1_packet_prefix_count": 28, "first0_packet_prefix_count": 28,
        "raw17_sha256": _digest(raw[:17]), "raw28_sha256": _digest(raw[:28]),
        "macro1_complete_sha256": _digest(macro1_terms),
    }


def _trace_chain(proofs: Mapping[str, Mapping[str, Any]], signed_source: Mapping[str, Any],
                 physical_source: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "pieces": [
            {"name": "micro_macro", "lower": "second_micro_native2", "upper": "frozen_macro_native0",
             "signed_piece": "-D(second_micro native2)+D(frozen_macro native0)"},
            {"name": "macro_first", "lower": "frozen_macro_native1", "upper": "first_switch_native0",
             "signed_piece": "-D(frozen_macro native1)+D(first_switch native0)"},
        ],
        "chart_contributions": [
            {"name": "frozen_macro", "lower": "frozen_macro native0", "upper": "frozen_macro native1",
             "signed_piece": "+D(frozen_macro native0)-D(frozen_macro native1)",
             "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
             "physical_R_once": True, "ordinary_Z_orders": [0, 1]},
        ],
        "function": "-D(second_micro native2)+D(frozen_macro native0)-D(frozen_macro native1)+D(first_switch native0)=0",
        "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
        "ordinary_Z_orders": [0, 1], "physical_R_once": True,
        "leading_subtraction_retained": True, "exact_function_identity_zero": True,
        "endpoint_values_are_not_graph_nodes": True,
        "source_functions_not_cap_values": True,
        "macro_chart_signed_orientation": {
            "lower": "frozen_macro native0", "upper": "frozen_macro native1",
            "signed_piece": "D(frozen_macro native0)-D(frozen_macro native1)",
            "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
            "physical_R_once": True, "ordinary_Z_orders": [0, 1],
        },
        "signed_endpoint_source": copy.deepcopy(dict(signed_source)),
        "physical_density_source": copy.deepcopy(dict(physical_source)),
        "source_proof_sha256": {name: _digest(value) for name, value in proofs.items()},
    }


def _source_binding(owner: Any, parent: Mapping[str, Any], family: Any,
                    graph: Mapping[str, Any], proofs: Mapping[str, Mapping[str, Any]],
                    ast_binding: Mapping[str, Any], signed_source: Mapping[str, Any],
                    physical_source: Mapping[str, Any]) -> dict[str, Any]:
    parent_source = _mapping(parent.get("numerical_source_binding") or parent.get("source_binding"), "parent numerical source binding")
    inherited = {key: copy.deepcopy(parent_source[key]) for key in (
        "delta_definition", "physical_radius_power", "ordinary_Z_orders", "physical_R_once",
        "leading_subtraction_retained", "cap_is_not_signed_function_value", "all_other_errors_retained",
        "source_functions_values_or_caps_not_selected", "accepted_ancestor_constructors_not_run")
        if key in parent_source}
    source = {
        **inherited, "passed": True, "source_family": copy.deepcopy(family),
        "function": "-D(second_micro native2)+D(frozen_macro native0)-D(frozen_macro native1)+D(first_switch native0)",
        "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
        "physical_radius_power": 1, "ordinary_Z_orders": [0, 1], "physical_R_once": True,
        "leading_subtraction_retained": True, "cap_is_not_signed_function_value": True,
        "all_other_errors_retained": True, "source_functions_values_or_caps_not_selected": True,
        "accepted_ancestor_constructors_not_run": True,
        "actual_macro_endpoint_graph_binding": copy.deepcopy(dict(graph)),
        "micro_macro_source_join": copy.deepcopy(dict(proofs["micro_macro"])),
        "macro_first_source_join": copy.deepcopy(dict(proofs["macro_first"])),
        "micro_macro_source_join_proved": True, "macro_first_source_join_proved": True,
        "source_proof_sha256": {name: _digest(value) for name, value in proofs.items()},
        "parent_switch_report": {"path": PARENT, "sha256": PARENT_SHA},
        "parent_source_binding_sha256": _digest(parent_source),
        "raw_macro_viscous_index": 27,
        "previous_reduced_removed_terms": [MICRO_TERM, MACRO_TERM, FIRST_TERM],
        "deleted_endpoint_atoms_provenance_only": True,
        "no_active_endpoint_budget": True,
        "no_cap_subtraction_or_midpoint_closure": True,
        "macro_chart_signed_orientation": {
            "lower": "frozen_macro native0", "upper": "frozen_macro native1",
            "signed_piece": "D(frozen_macro native0)-D(frozen_macro native1)",
            "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
            "physical_R_once": True, "ordinary_Z_orders": [0, 1],
        },
        "generic_nodes_reused_only_under_explicit_chart_scope": True,
        "generic_history_bounds_inherited_not_recertified": True,
        "source_ast_composition": copy.deepcopy(dict(ast_binding)),
        "signed_endpoint_source": copy.deepcopy(dict(signed_source)),
        "physical_density_source": copy.deepcopy(dict(physical_source)),
        "macro_endpoint_chart_scopes": {name: copy.deepcopy(dict(value)) for name, value in graph["charts"].items()},
    }
    return source


def _evaluate_sample(owner: Any, parent: Mapping[str, Any], macro: Mapping[str, Any], second: Mapping[str, Any],
                     first: Mapping[str, Any], family: Any, proofs: Mapping[str, Mapping[str, Any]],
                     source: Mapping[str, Any], graph: Mapping[str, Any], ast_binding: Mapping[str, Any],
                     signed_source: Mapping[str, Any],
                     physical_source: Mapping[str, Any],
                     parent_sample: Mapping[str, Any]) -> dict[str, Any]:
    packet = _packet_witness(parent_sample, macro, second, first, family)
    raw = parent_sample["raw_m15_C0_C1_terms"]
    previous = parent_sample["correlated_reduced_m15_C0_C1_terms"]
    removed_indices = {name: next(i for i, term in enumerate(previous) if term.get("name") == name)
                      for name in (MICRO_TERM, MACRO_TERM, FIRST_TERM)}
    reduced = [copy.deepcopy(term) for term in previous if term.get("name") not in removed_indices]
    if len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 183-term view required")
    old_unchanged = [term for term in previous if term.get("name") not in removed_indices]
    if len(old_unchanged) != UNCHANGED_COUNT or reduced != old_unchanged:
        raise ValueError("Other 183 correlated sectors must remain verbatim")
    target = copy.deepcopy(parent_sample.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")
    decomposition = _macro_decomposition(macro, ast_binding, str(parent_sample.get("Z")))
    provenance = []
    source_packets = {
        MICRO_TERM: {"left": "second_micro native2", "right": "frozen_macro native0", "packet": "second2 -> macro0"},
        MACRO_TERM: {"left": "frozen_macro native0", "right": "frozen_macro native1", "packet": "macro finite-N viscous raw27"},
        FIRST_TERM: {"left": "frozen_macro native1", "right": "first_switch native0", "packet": "macro1 -> first0"},
    }
    for name in (MICRO_TERM, MACRO_TERM, FIRST_TERM):
        row = copy.deepcopy(previous[removed_indices[name]])
        provenance.append({"name": name, "previous_reduced_index": removed_indices[name],
                           "term": row, "active": False, "provenance_only": True,
                           "source_packet": source_packets[name], "same_source_proof_passed": True,
                           "exact_packet_prefix_and_selected_N_passed": True,
                           "no_cap_value_or_midpoint": True,
                           "macro_raw27_decomposition": copy.deepcopy(decomposition) if name == MACRO_TERM else None})
    selected = packet["selected_integer"]
    try:
        caps = owner.largest_caps(reduced, target)
    except Exception as exc:
        raise ValueError("Owner largest_caps pure utility is required") from exc
    return {
        "Z": str(parent_sample.get("Z")), "source_family": copy.deepcopy(family),
        "raw_m15_C0_C1_terms": copy.deepcopy(raw), "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT, "previous_reduced_terms": PARENT_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT, "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": _trace_chain(proofs, signed_source, physical_source),
        "signed_endpoint_source": copy.deepcopy(dict(signed_source)),
        "physical_density_source": copy.deepcopy(dict(physical_source)),
        "source_seam_proofs": {name: copy.deepcopy(dict(value)) for name, value in proofs.items()},
        "actual_signed_function_ports": [copy.deepcopy(graph["charts"][SECOND_MICRO]),
                                          copy.deepcopy(graph["charts"][FROZEN_MACRO]),
                                          copy.deepcopy(graph["charts"][FIRST_SWITCH])],
        "removed_endpoint_provenance": provenance,
        "old_removed_endpoint_provenance": copy.deepcopy(provenance),
        "macro_raw27_decomposition": decomposition,
        "packet_prefix_witness": packet,
        "selected_N_evidence": copy.deepcopy(selected), "independent_target": target,
        "largest_remaining_individual_log_caps": caps,
        "contract": {
            "raw_213_verbatim": True, "other_183_verbatim": True, "independent_target_unchanged": True,
            "reserved_retained_once": sum(term.get("name") == "reserved_n1_own_swirl_once" for term in reduced) == 1,
            "six_signed_inner_lower_retained": True, "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True, "micro_macro_source_join_proved": True,
            "macro_first_source_join_proved": True, "macro_raw27_endpoint_only_proved": True,
            "three_endpoint_sectors_removed_once": True, "deleted_endpoint_atoms_provenance_only": True,
            "no_active_endpoint_budget": True, "source_functions_not_cap_values": True,
            "no_error_cap_subtraction_or_midpoint": True, "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True, "source_bound_signed_trace_contraction": True,
            "previous_186_to_new_183": True,
        },
        "limits": {
            "micro_macro_internal_trace_resolved": True, "macro_first_internal_trace_resolved": True,
            "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False, "numeric_Rp_closure": False,
            "global_five_moment_common_N_closure": False, "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False, "temporal_n_recursion_and_summation": False,
            "mean_oscillatory_correction": False, "independent_corrected_Cartesian_NS": False,
        },
    }


def build_macro_numeric(owner: Any, seam_proofs: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved raw 213 -> correlated 183 macro endpoint ledger."""
    parent = _parent_report(owner)
    family = _family(owner, parent)
    graph = _graph_binding(owner)
    proofs = _proofs(seam_proofs, graph)
    macro = _macro_report(owner); second = _second_report(owner); first = _first_report(owner)
    for label, report in (("macro", macro), ("second micro", second), ("first switch", first)):
        if report.get("source_family") != family or _mapping(report.get("checks"), f"{label} checks").get("passed") is not True:
            raise ValueError(f"Accepted {label} report is not passed")
    ast_binding = _pin_ast_composition(owner)
    signed_source = _pin_signed_trace_source(owner)
    physical_source = _pin_physical_density_source(owner, graph)
    samples: list[dict[str, Any]] = []
    for Z in SAVED_Z:
        parent_sample = _find(parent.get("samples"), label="accepted switch parent sample", Z=Z)
        source = _source_binding(owner, parent, family, graph, proofs, ast_binding, signed_source, physical_source)
        samples.append(_evaluate_sample(owner, parent, macro, second, first, family, proofs, source, graph,
                                        ast_binding, signed_source, physical_source, parent_sample))
    selected = [sample["selected_N_evidence"] for sample in samples]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All macro samples must use one selected N")
    source = _source_binding(owner, parent, family, graph, proofs, ast_binding, signed_source, physical_source)
    source.update({"selected_N_evidence": copy.deepcopy(selected[0]),
                   "micro_macro_source_proof_sha256": _digest(proofs["micro_macro"]),
                   "macro_first_source_proof_sha256": _digest(proofs["macro_first"]),
                   "macro_raw27_source_proof_sha256": _digest(ast_binding),
                   "signed_endpoint_source_sha256": _digest(signed_source),
                   "physical_density_source_sha256": _digest(physical_source),
                   "parent_signed_endpoint_bindings": copy.deepcopy(
                       samples[0]["packet_prefix_witness"]["parent_signed_endpoint_bindings"])})
    return {
        "source_binding": source, "samples": samples,
        "source_seam_proofs": {"passed": True, "joins": {name: copy.deepcopy(dict(value)) for name, value in proofs.items()}},
        "contract": {
            "raw_213_retained": True, "previous_186_to_new_183": True, "unchanged_other_183": True,
            "three_signed_endpoint_sectors_removed": True, "micro_macro_source_join_required": True,
            "macro_first_source_join_required": True, "macro_raw27_endpoint_only_required": True,
            "deleted_endpoint_provenance_only": True, "no_active_endpoint_budget": True,
            "source_functions_not_cap_values": True, "ordinary_Z1_and_physical_R_once": True,
            "raw_target_reserved_and_signed_inner_lower_retained": True, "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True,
        },
        "limits": {
            "numeric_Rp_closure": False, "global_all_five_common_N_closure": False,
            "full_goal_complete": False, "endpoint_signed_values_resolved": False,
        },
    }


__all__ = ["build_macro_numeric"]
