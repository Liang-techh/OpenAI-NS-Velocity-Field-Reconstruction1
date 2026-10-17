"""Contract the saved second-switch to long upstream viscous chain.

This module is a saved-packet adapter.  It consumes the accepted tail report,
the pinned R110 packet, and two caller-supplied source seams.  It never loads
an historical owner, constructor, producer, or report checker.  The two
interfaces are required before the three endpoint sectors are contracted:

``second_switch native 1 = post_power native 0`` and
``post_power native 1 = long_reshape native 0``.

The returned packet keeps the raw 213-term state and all unrelated terms
verbatim.  It removes the second-switch and post-power finite-N viscous terms
and the prior long-to-compact external term, replacing them with one directed
positive error enclosure for the surviving second-switch native 0 endpoint.
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
PREVIOUS_REDUCED_COUNT = 189
NEW_REDUCED_COUNT = 187
UNCHANGED_COUNT = 186

PARENT = PREFIX + "compliant_current_original_tail_m15_numeric.json"
PARENT_SHA = "a394d3b538e420936e382b36b63dc7405b7a098b57a557c668deba7a2cc60c9a"
R110_REPORT = PREFIX + "compliant_current_original_r110_m15_numeric.json"
R110_REPORT_SHA = "7226a0546b338c6834d035f6daa86707091c4c1e2901f631ea8a2351231bada1"
R110_PRESSURE_SOURCE = "src/openai_ns_reconstruction/original_r110_pressure_rows.py"
R110_PRESSURE_SOURCE_SHA = "9f60fc44291c16f0906ce73b19a414645fa1b2a2af1bec3ca2d3fc76e77761ef"
R110_M15_SOURCE = "src/openai_ns_reconstruction/original_r110_m15.py"
R110_M15_SOURCE_SHA = "ca6e74ff8ca4eaaf48e2693328c1b922ba8b683b28b87712941541f497cec20d"

SECOND_VISCOUS_INDEX = 68
POST_POWER_VISCOUS_INDEX = 87
SECOND_INCOMING_COUNT = 51
SECOND_TERM_COUNT = 69
POST_POWER_INCOMING_COUNT = 69
POST_POWER_TERM_COUNT = 88

OLD_SECOND_TERM = "second_switch_finite_N_viscous_error"
OLD_POST_POWER_TERM = "post_power_finite_N_viscous_error"
OLD_EXTERNAL_TERM = "long_compact2Rc_finite_N_viscous_external_endpoints"
NEW_EXTERNAL_TERM = "second_switch_long_finite_N_viscous_external_endpoints"

SECOND_Q0 = {"kind": "exact_decimal_second_switch_coordinate", "value": "0"}
SECOND_Q1 = {"kind": "exact_decimal_second_switch_coordinate", "value": "1"}
POST_POWER_Q0 = {"kind": "exact_decimal_post_power_coordinate", "value": "0"}
POST_POWER_Q1 = {"kind": "exact_decimal_post_power_coordinate", "value": "1"}

# These are the archived graph handles.  Endpoint values are deliberately
# kept separate from graph handles in source-seam metadata.
SECOND_COORDINATE_NODE = 75949
POST_POWER_COORDINATE_NODE = 78472
LONG_COORDINATE_NODE = 79966
SECOND_ENDPOINT_NODE = 0
SECOND_UPPER_ENDPOINT_NODE = 1
SECOND_RADIUS_NODE = 88600
SECOND_JACOBIAN_NODE = 9
SECOND_QY_NODES = (154430, 154434)
SECOND_RQY_NODES = (154589, 154590)


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


def _coerce_coordinate(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return value


def _path(owner: Any, relative: str) -> Path:
    candidates: list[Path] = []
    directory = getattr(owner, "directory", None)
    if directory is not None:
        base = Path(directory)
        candidates.extend((base / relative, base / "experiments" / "root_st073" / relative))
    candidates.extend((ROOT / relative, ROOT / "experiments" / "root_st073" / relative))
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
    known_values = (
        hashes.get(relative),
        hashes.get(basename),
        code_hashes.get(relative),
        code_hashes.get(basename),
    )
    for known in known_values:
        if known is not None and known != actual:
            raise ValueError(f"Changed accepted source: {relative}")
    if actual != expected:
        raise ValueError(f"Unexpected accepted source: {relative}")
    return {"path": relative, "sha256": actual, "expected_sha256": expected, "passed": True}


def _load_json(owner: Any, relative: str, expected: str) -> Mapping[str, Any]:
    _pin_file(owner, relative, expected)
    value = json.loads(_path(owner, relative).read_text(encoding="utf-8"))
    return _mapping(value, relative)


def _parent_report(owner: Any) -> Mapping[str, Any]:
    for candidate in (getattr(owner, "parent", None), getattr(owner, "parent_report", None)):
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list):
            return candidate
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list):
            return report
    return _load_json(owner, PARENT, PARENT_SHA)


def _r110_report(owner: Any) -> Mapping[str, Any]:
    # A caller may provide the already hydrated accepted report.  Pin and
    # compare it to the accepted file before allowing it to feed the adapter;
    # an in-memory candidate is not an independent source of evidence.
    pinned = _load_json(owner, R110_REPORT, R110_REPORT_SHA)
    for name in ("r110", "r110_report", "second_report", "second"):
        candidate = getattr(owner, name, None)
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list) and isinstance(candidate.get("pressure_row_samples"), list):
            if _canonical(candidate) != _canonical(pinned):
                raise ValueError("In-memory R110 report differs from pinned accepted report")
            return candidate
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list) and isinstance(report.get("pressure_row_samples"), list):
            if _canonical(report) != _canonical(pinned):
                raise ValueError("In-memory R110 report differs from pinned accepted report")
            return report
    return pinned


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


def _term_is_zero(term: Mapping[str, Any]) -> bool:
    coefficients = term.get("coefficient_C0_C1")
    if not isinstance(coefficients, list) or len(coefficients) != 2:
        return False
    for row in coefficients:
        if not isinstance(row, Mapping):
            return False
        if row.get("approximation") not in ("0.0", "0") or row.get("absolute_error_upper") not in ("0.0", "0"):
            return False
        if row.get("raw_directed_error_exact_mpf_tuple") != [0, 0, 0, 0]:
            return False
        if row.get("source_enclosure_exact_mpf_tuples") != [[0, 0, 0, 0], [0, 0, 0, 0]]:
            return False
    return True


def _positive_atom(value: Any, label: str) -> Mapping[str, Any]:
    atom = _mapping(value, label)
    if atom.get("positive_error_bound_atom_only") is not True or atom.get("source_value_selected") is not False:
        raise ValueError(f"{label} must remain an error-only positive atom")
    if not isinstance(atom.get("source_enclosure_exact_mpf_tuples"), list):
        raise ValueError(f"{label} source enclosure is missing")
    return atom


def _unwrapped(value: Any, key: str) -> Mapping[str, Any]:
    mapping = _mapping(value, key)
    if key in mapping and isinstance(mapping.get(key), Mapping):
        return _mapping(mapping[key], key)
    return mapping


def _pin_r110_ast(owner: Any) -> dict[str, Any]:
    pressure_pin = _pin_file(owner, R110_PRESSURE_SOURCE, R110_PRESSURE_SOURCE_SHA)
    m15_pin = _pin_file(owner, R110_M15_SOURCE, R110_M15_SOURCE_SHA)
    path = _path(owner, R110_M15_SOURCE)
    tree = ast.parse(path.read_bytes(), filename=str(path))
    cls = next((node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "OriginalR110M15"), None)
    if cls is None:
        raise ValueError("OriginalR110M15 class missing from pinned source")
    method = next((node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == "finite_errors"), None)
    if method is None:
        raise ValueError("OriginalR110M15.finite_errors method missing")
    calls = {
        node.func.attr
        for node in ast.walk(method)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    required = {"background", "correction_densities", "sum_logs"}
    if not required.issubset(calls):
        raise ValueError("Pinned R110 finite_errors AST lacks complete endpoint/viscous rows")
    return {
        "passed": True,
        "pressure_rows_source": pressure_pin,
        "finite_errors_source": m15_pin,
        "class": "OriginalR110M15",
        "method": "OriginalR110M15.finite_errors",
        "method_ast_sha256": hashlib.sha256(ast.dump(method, include_attributes=False).encode("utf-8")).hexdigest(),
        "required_calls": sorted(required),
        "complete_prefix_error_source": True,
        "source_not_executed": True,
    }


def _proofs(seam_proofs: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    if not isinstance(seam_proofs, Mapping) or seam_proofs.get("passed") is not True:
        raise ValueError("Second-switch upstream source seam proofs must pass")
    candidate = seam_proofs.get("joins")
    if not isinstance(candidate, Mapping):
        candidate = seam_proofs.get("seams")
    if not isinstance(candidate, Mapping):
        candidate = seam_proofs
    expected = {
        "second_power": ("second_switch", "post_power", 1, 0, SECOND_COORDINATE_NODE, POST_POWER_COORDINATE_NODE),
        "power_long": ("post_power", "long_reshape", 1, 0, POST_POWER_COORDINATE_NODE, LONG_COORDINATE_NODE),
    }
    result: dict[str, Mapping[str, Any]] = {}
    for name, (left_chart, right_chart, left_value, right_value, left_coord, right_coord) in expected.items():
        proof = candidate.get(name)
        if not isinstance(proof, Mapping) or proof.get("passed") is not True:
            raise ValueError(f"Accepted {name} source proof must pass")
        for key, value in (("left_chart", left_chart), ("right_chart", right_chart),
                           ("left_endpoint", left_value), ("right_endpoint", right_value),
                           ("left_native_query", left_value), ("right_native_query", right_value)):
            if key in proof and proof.get(key) != value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        for key, value in (("left_coordinate_node", left_coord), ("right_coordinate_node", right_coord),
                           ("left_native_coordinate_node", left_coord), ("right_native_coordinate_node", right_coord)):
            if key in proof and proof.get(key) != value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        # Endpoint values are 1 and 0; source helpers may additionally expose
        # the actual query graph handle for the left endpoint.  Keep those
        # namespaces separate and accept the archived query handle 48 when it
        # is present, while still rejecting a value/node mix-up.
        node_values = {
            "left_endpoint_node": (left_value, 48),
            "right_endpoint_node": (right_value,),
            "left_native_query_node": (left_value, 48),
            "right_native_query_node": (right_value,),
        }
        for key, allowed in node_values.items():
            if key in proof and proof.get(key) not in allowed:
                raise ValueError(f"{name} source proof has the wrong {key}")
        for side, coord, endpoint, endpoint_value in (("left", left_coord, left_value, left_value),
                                                        ("right", right_coord, right_value, right_value)):
            scopes = [proof[key] for key in (f"{side}_scope", f"{side}_endpoint_scope")
                      if isinstance(proof.get(key), Mapping)]
            for scope in scopes:
                if "coordinate_node" in scope and scope.get("coordinate_node") != coord:
                    raise ValueError(f"{name} source proof has wrong {side} coordinate scope")
                if "native_coordinate_node" in scope and scope.get("native_coordinate_node") != coord:
                    raise ValueError(f"{name} source proof has wrong {side} native coordinate scope")
                if "endpoint_node" in scope and scope.get("endpoint_node") != endpoint:
                    raise ValueError(f"{name} source proof has wrong {side} endpoint scope")
                for key in ("native_query_node", "query_node"):
                    if key in scope and scope.get(key) not in (endpoint, 48):
                        raise ValueError(f"{name} source proof has wrong {side} query node")
                if "native_query_value" in scope and str(scope.get("native_query_value")) != str(endpoint_value):
                    raise ValueError(f"{name} source proof has wrong {side} query value")
        if (proof.get("source_values_selected") is True or proof.get("numeric_caps_used") is True
                or proof.get("source_value_selected") is True or proof.get("ancestor_replayed") is True):
            raise ValueError(f"{name} source proof selected a value or cap")
        result[name] = copy.deepcopy(dict(proof))
    return result


def _chart_bindings(owner: Any) -> tuple[Mapping[str, Any] | None, Mapping[str, Any] | None]:
    loop = getattr(owner, "l", None)
    modified = getattr(owner, "m", None)
    if loop is None:
        loop = getattr(owner, "loop", None)
    if modified is None:
        modified = getattr(owner, "modified", None)
    try:
        loop_map = _unwrapped(loop, "prescribed_loop_C5_bindings") if loop is not None else None
    except ValueError:
        loop_map = None
    try:
        modified_map = _unwrapped(modified, "modified_Q4_pressure_C3_bindings") if modified is not None else None
    except ValueError:
        modified_map = None
    return loop_map, modified_map


def _second_native_guard(owner: Any, proofs: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    guards = getattr(owner, "native_guards", None)
    if not isinstance(guards, Mapping):
        guards = getattr(owner, "native", None)
    if isinstance(guards, Mapping):
        chart_rows = guards.get("charts")
        if isinstance(chart_rows, Mapping) and isinstance(chart_rows.get("second_switch"), Mapping):
            row = chart_rows["second_switch"]
            if guards.get("passed") is True and row.get("passed") is True:
                return {
                    "passed": True,
                    "chart": "second_switch",
                    "guard": copy.deepcopy(dict(row)),
                    "source_guard_packet": copy.deepcopy(dict(guards)),
                    "ordinary_Z_orders": [0, 1],
                    "physical_RQy_factor_once": True,
                    "physical_Uy_Pstar_factor_once": True,
                }
    for proof in proofs.values():
        for key in ("second_native_guard", "native_guard", "actual_native_guard"):
            value = proof.get(key)
            if isinstance(value, Mapping) and value.get("passed") is True:
                return {"passed": True, "chart": "second_switch", "guard": copy.deepcopy(dict(value)),
                        "ordinary_Z_orders": [0, 1], "physical_RQy_factor_once": True,
                        "physical_Uy_Pstar_factor_once": True}
    raise ValueError("Passed second_switch physical native guard required")


def _graph_binding(owner: Any, parent: Mapping[str, Any], proofs: Mapping[str, Mapping[str, Any]], guard: Mapping[str, Any]) -> dict[str, Any]:
    prior: Mapping[str, Any] | None = None
    for sample in parent.get("samples", []):
        if not isinstance(sample, Mapping):
            continue
        budget = sample.get("external_endpoint_budget_witness")
        if not isinstance(budget, Mapping):
            continue
        own = budget.get("own_longreshape0")
        if isinstance(own, Mapping) and isinstance(own.get("actual_long0_graph_envelope_binding"), Mapping):
            prior = own["actual_long0_graph_envelope_binding"]
            break
    if prior is None:
        raise ValueError("Parent long native0 graph envelope certificate required")
    prior = _mapping(prior, "parent long graph envelope")
    if prior.get("passed") is not True:
        raise ValueError("Parent long graph envelope certificate must pass")
    for key in ("exact_numeric_Qy_coefficients", "ordinary_Z_product_rule",
                "actual_archived_pressure_formula_certificate", "actual_archived_radial_formula_certificate"):
        if key not in prior:
            raise ValueError(f"Parent universal Qy certificate missing {key}")

    loop, modified = _chart_bindings(owner)
    loop_chart: Mapping[str, Any] | None = None
    modified_chart: Mapping[str, Any] | None = None
    if loop is not None:
        charts = loop.get("charts")
        if isinstance(charts, Mapping) and isinstance(charts.get("second_switch"), Mapping):
            loop_chart = charts["second_switch"]
    if modified is not None:
        charts = modified.get("charts")
        if isinstance(charts, Mapping) and isinstance(charts.get("second_switch"), Mapping):
            modified_chart = charts["second_switch"]
    coordinate = loop_chart.get("coordinate") if loop_chart is not None else SECOND_COORDINATE_NODE
    domain = loop_chart.get("domain") if loop_chart is not None else [SECOND_ENDPOINT_NODE, SECOND_UPPER_ENDPOINT_NODE]
    radius = loop_chart.get("radius") if loop_chart is not None else SECOND_RADIUS_NODE
    physical_measure = loop_chart.get("physical_measure") if loop_chart is not None else SECOND_JACOBIAN_NODE
    if coordinate != SECOND_COORDINATE_NODE or domain != [SECOND_ENDPOINT_NODE, SECOND_UPPER_ENDPOINT_NODE]:
        raise ValueError("Second-switch native coordinate/domain changed")
    if radius != SECOND_RADIUS_NODE or physical_measure != SECOND_JACOBIAN_NODE:
        raise ValueError("Second-switch physical radius/Jacobian changed")

    qy_nodes = list(SECOND_QY_NODES)
    rqy_nodes = list(SECOND_RQY_NODES)
    if modified_chart is not None:
        physical = _mapping(modified_chart.get("physical"), "second-switch modified physical")
        finite = _mapping(modified_chart.get("finite_boundary"), "second-switch finite boundary")
        if list(physical.get("Q_y", [])[:2]) != qy_nodes or list(finite.get("R_Q_y", [])[:2]) != rqy_nodes:
            raise ValueError("Second-switch Qy/RQy graph rows changed")

    nodes = getattr(owner, "nodes", None)
    graph_rows: list[dict[str, Any]] = []
    fingerprints: dict[int, str] = {}
    if isinstance(nodes, (list, tuple, Mapping)):
        for node in rqy_nodes:
            try:
                row = nodes[node] if not isinstance(nodes, Mapping) else nodes[node]
            except (IndexError, KeyError, TypeError):
                raise ValueError("Second-switch RQy graph node missing") from None
            if not isinstance(row, Mapping) or row.get("operation") != "product":
                raise ValueError("Second-switch RQy must be an actual product")
            arguments = row.get("arguments", [])
            if not isinstance(arguments, list) or arguments.count(radius) != 1:
                raise ValueError("Second-switch RQy must contain one physical R factor")
        from .original_active_viscous_m15 import ActiveEndpointFingerprint
        fingerprint = ActiveEndpointFingerprint(nodes, dict(coordinate=coordinate, loop=loop.get("charts", {}) if loop else {}))
        scope = {coordinate: fingerprint.at(domain[0], {})}
        for node in (*qy_nodes, *rqy_nodes):
            fingerprints[node] = fingerprint.at(node, scope)[0]

    supplied_rows = None
    for proof in proofs.values():
        for key in ("second0_graph_function_rows", "graph_function_rows", "actual_graph_function_rows"):
            if isinstance(proof.get(key), list) and len(proof[key]) >= 2:
                supplied_rows = proof[key][:2]
                break
        if supplied_rows is not None:
            break
    for order in (0, 1):
        row = {
            "ordinary_Z_order": order,
            "Qy_function_node": qy_nodes[order],
            "RQy_function_node": rqy_nodes[order],
            "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": SECOND_ENDPOINT_NODE},
            "physical_radius_node": radius,
            "physical_R_explicit_once": True,
        }
        if qy_nodes[order] in fingerprints:
            row["Qy_at_second0_sha256"] = fingerprints[qy_nodes[order]]
        if rqy_nodes[order] in fingerprints:
            row["RQy_at_second0_sha256"] = fingerprints[rqy_nodes[order]]
        if supplied_rows is not None:
            supplied = _mapping(supplied_rows[order], "second graph function row")
            for key in ("Qy_at_second0_sha256", "RQy_at_second0_sha256"):
                if key in supplied and key in row and supplied[key] != row[key]:
                    raise ValueError("Second-switch graph fingerprint changed")
            row["source_proof_row"] = copy.deepcopy(dict(supplied))
        graph_rows.append(row)
    inherited_keys = ("exact_numeric_Qy_coefficients", "exact_coefficient_geometry_statements",
                      "ordinary_Z_product_rule", "actual_archived_pressure_formula_certificate",
                      "actual_archived_radial_formula_certificate")
    return {
        "passed": True,
        "method": "source_bound_linear_Qy_envelope_at_actual_second_switch_native0",
        "inherited_long0_envelope_sha256": _digest(prior),
        **{key: copy.deepcopy(prior[key]) for key in inherited_keys},
        "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": SECOND_ENDPOINT_NODE},
        "graph_function_rows": graph_rows,
        "chart": "second_switch",
        "derivative_chart": "switch_second",
        "physical_radius_node": radius,
        "native_Jacobian_node": physical_measure,
        "ordinary_Z_orders": [0, 1],
        "physical_RQy_factor_once": True,
        "physical_Uy_Pstar_factor_once": True,
        "second_native_guard": copy.deepcopy(dict(guard)),
        "generic_nodes_reused_only_under_explicit_scope": True,
        "generic_history_bounds_inherited_not_recertified": True,
        "source_function_values_or_caps_not_selected": True,
        "source_function_not_identified_with_post_power1": True,
        "actual_second_switch0_scope_only": True,
    }


def _validate_packets(owner: Any, parent: Mapping[str, Any], report: Mapping[str, Any], parent_sample: Mapping[str, Any]) -> dict[str, Any]:
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Parent raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PREVIOUS_REDUCED_COUNT:
        raise ValueError("Parent tail 189-term view required")
    _term_names(raw, "parent raw terms")
    _term_names(previous, "parent reduced terms")
    if raw[SECOND_VISCOUS_INDEX].get("name") != OLD_SECOND_TERM or raw[POST_POWER_VISCOUS_INDEX].get("name") != OLD_POST_POWER_TERM:
        raise ValueError("Parent raw second/post-power viscous indices changed")
    family = getattr(owner, "family", parent.get("source_family"))
    if report.get("source_family") != family:
        raise ValueError("R110 source family changed")
    selected_parent = _mapping(parent_sample.get("selected_N_evidence"), "parent selected-N evidence")
    Z = str(parent_sample.get("Z"))
    rows: dict[str, Mapping[str, Any]] = {}
    specs = {
        "second0": ("second_switch", SECOND_Q0, SECOND_INCOMING_COUNT, SECOND_TERM_COUNT),
        "second1": ("second_switch", SECOND_Q1, SECOND_INCOMING_COUNT, SECOND_TERM_COUNT),
        "post0": ("post_power", POST_POWER_Q0, POST_POWER_INCOMING_COUNT, POST_POWER_TERM_COUNT),
        "post1": ("post_power", POST_POWER_Q1, POST_POWER_INCOMING_COUNT, POST_POWER_TERM_COUNT),
    }
    for name, (chart, coordinate, prefix_count, term_count) in specs.items():
        row = _find(report.get("samples"), label=f"R110 {name}", Z=Z, chart=chart, coordinate=coordinate)
        if row.get("source_family") != family:
            raise ValueError(f"R110 {name} source family changed")
        terms = row.get("m15_C0_C1_terms")
        incoming = row.get("incoming_m15_terms")
        if not isinstance(terms, list) or len(terms) != term_count:
            raise ValueError(f"R110 {name} complete packet count changed")
        if not isinstance(incoming, list) or len(incoming) != prefix_count:
            raise ValueError(f"R110 {name} incoming prefix count changed")
        _term_names(terms, f"R110 {name} terms")
        _term_names(incoming, f"R110 {name} incoming")
        if row.get("actual_selected_integer") != selected_parent:
            raise ValueError(f"R110 {name} selected N changed")
        if row.get("exact_empty_increment") != name.endswith("0"):
            raise ValueError(f"R110 {name} exact-empty endpoint flag changed")
        if name.endswith("0") and not all(_term_is_zero(term) for term in terms[prefix_count:]):
            raise ValueError(f"R110 {name} tails are not exact zero packets")
        rows[name] = row

    raw_second_prefix = raw[:SECOND_INCOMING_COUNT]
    _exact("second0 incoming versus parent raw prefix", rows["second0"]["incoming_m15_terms"], raw_second_prefix)
    _exact("second1 incoming versus parent raw prefix", rows["second1"]["incoming_m15_terms"], raw_second_prefix)
    _exact("second0 terms prefix versus incoming", rows["second0"]["m15_C0_C1_terms"][:SECOND_INCOMING_COUNT], rows["second0"]["incoming_m15_terms"])
    _exact("second1 terms prefix versus incoming", rows["second1"]["m15_C0_C1_terms"][:SECOND_INCOMING_COUNT], rows["second1"]["incoming_m15_terms"])
    _exact("second1 complete packet versus post0 incoming", rows["second1"]["m15_C0_C1_terms"], rows["post0"]["incoming_m15_terms"])
    _exact("second1 complete packet versus post1 incoming", rows["second1"]["m15_C0_C1_terms"], rows["post1"]["incoming_m15_terms"])
    _exact("post0 terms prefix versus incoming", rows["post0"]["m15_C0_C1_terms"][:POST_POWER_INCOMING_COUNT], rows["post0"]["incoming_m15_terms"])
    _exact("post1 terms prefix versus incoming", rows["post1"]["m15_C0_C1_terms"][:POST_POWER_INCOMING_COUNT], rows["post1"]["incoming_m15_terms"])

    return {
        "Z": Z,
        "selected_integer": copy.deepcopy(dict(selected_parent)),
        "rows": rows,
        "exact_prefix_witness": {
            "passed": True,
            "second_prefix_count": SECOND_INCOMING_COUNT,
            "second_term_count": SECOND_TERM_COUNT,
            "post_power_prefix_count": POST_POWER_INCOMING_COUNT,
            "post_power_term_count": POST_POWER_TERM_COUNT,
            "raw_prefix_equals_second0_incoming": True,
            "raw_prefix_equals_second1_incoming": True,
            "second0_terms_prefix_equals_incoming": True,
            "second1_terms_prefix_equals_incoming": True,
            "second1_complete_equals_post0_incoming": True,
            "second1_complete_equals_post1_incoming": True,
            "post0_terms_prefix_equals_incoming": True,
            "post1_terms_prefix_equals_incoming": True,
            "full_dictionary_equality": True,
            "raw_prefix_sha256": _digest(raw_second_prefix),
            "second1_complete_sha256": _digest(rows["second1"]["m15_C0_C1_terms"]),
            "post0_incoming_sha256": _digest(rows["post0"]["incoming_m15_terms"]),
        },
    }


def _second0_pressure_witness(owner: Any, report: Mapping[str, Any], Z: str, selected: Mapping[str, Any], graph: Mapping[str, Any]) -> tuple[list[Any], dict[str, Any]]:
    row = _find(report.get("pressure_row_samples"), label="second-switch native0 pressure row", Z=Z, coordinate=SECOND_Q0)
    if row.get("source_family") != getattr(owner, "family", report.get("source_family")):
        raise ValueError("Second-switch native0 pressure source family changed")
    if row.get("selected_integer") != selected:
        raise ValueError("Second-switch native0 selected N changed")
    cover = _mapping(row.get("source_cover"), "second-switch source cover")
    if cover.get("generic_derivative_chart") != "switch_second":
        raise ValueError("Second-switch native0 source derivative chart changed")
    histories = cover.get("cumulative_history_charts")
    if not isinstance(histories, list) or not histories or histories[-1] != "second_switch":
        raise ValueError("Second-switch cumulative source history changed")
    contract = _mapping(row.get("contract"), "second-switch pressure contract")
    for key in ("original_Pstar_once", "finite_R_Qy_and_ordinary_Z1_retained", "actual_local_chart_measure"):
        if contract.get(key) is not True:
            raise ValueError(f"Second-switch pressure contract lacks {key}")
    if cover.get("all_axial_frames_retained") is not True:
        raise ValueError("Second-switch source cover must retain all axial frames")
    physical_rows = _mapping(row.get("physical_pressure_rows"), "second-switch physical pressure rows")
    qy = _mapping(physical_rows.get("Q_y"), "second-switch Qy rows")
    endpoint = _mapping(row.get("R_Q_y_endpoint"), "second-switch RQy endpoint")
    errors = qy.get("finite_N_absolute_error_log_bounds")
    endpoint_errors = endpoint.get("finite_N_coefficient_error_log_bounds")
    if errors != endpoint_errors or not isinstance(errors, list) or len(errors) != 2:
        raise ValueError("Second-switch Qy C0/Z1 finite-N rows required")
    radius = endpoint.get("finite_N_log_physical_factor")
    geometry = _mapping(row.get("geometry"), "second-switch pressure geometry")
    if radius != geometry.get("log_physical_radius"):
        raise ValueError("Second-switch physical R factor must be sourced once")
    weighted: list[Any] = []
    for error in errors:
        error_map = _mapping(error, "second-switch Qy error row")
        if error_map.get("exact_zero") is True:
            weighted.append(None)
        else:
            upper = _mapping(error_map.get("log_absolute_upper"), "second-switch Qy log error")
            weighted.append(owner.record(owner.decode(upper) + owner.decode(radius)))
    own_rows = {
        name: copy.deepcopy(_mapping(physical_rows.get(name), f"second-switch {name} rows").get("finite_N_absolute_error_log_bounds"))
        for name in ("U", "m", "U_y")
    }
    if any(not isinstance(value, list) for value in own_rows.values()):
        raise ValueError("Second-switch U/m/Uy source error rows missing")
    slow_fast = copy.deepcopy(row.get("finite_N_axial_radial_error_components"))
    if not isinstance(slow_fast, Mapping):
        raise ValueError("Second-switch slow/fast source error components missing")
    witness = {
        "source_report": R110_REPORT,
        "source_packet_sha256": _digest(row),
        "coordinate": copy.deepcopy(row.get("coordinate")),
        "selected_integer": copy.deepcopy(dict(selected)),
        "source_cover": copy.deepcopy(dict(cover)),
        "actual_own_Qy_error_rows": copy.deepcopy(errors),
        "actual_own_log_radius": copy.deepcopy(radius),
        "actual_linear_error_input_rows": own_rows,
        "own_U_m_Uy_error_rows": copy.deepcopy(own_rows),
        "actual_slow_fast_error_components": slow_fast,
        "own_slow_fast_components": copy.deepcopy(slow_fast),
        "radius_weighted_C0_C1_log_bounds": copy.deepcopy(weighted),
        "actual_second0_graph_envelope_binding": copy.deepcopy(dict(graph)),
        "generic_derivative_chart": "switch_second",
        "second_switch_native0_source_packet": True,
        "physical_R_once": True,
        "physical_Pstar_already_in_error": True,
        "source_function_values_evaluated": False,
        "post_power1_packet_not_reinterpreted": True,
        "no_caps_used": True,
    }
    return weighted, witness


def _parent_external_term(previous: list[Mapping[str, Any]]) -> Mapping[str, Any]:
    rows = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    if len(rows) != 1:
        raise ValueError("Parent long-to-compact external replacement required")
    logs = rows[0].get("log_scale_terms")
    if not isinstance(logs, list) or len(logs) != 1:
        raise ValueError("Parent long-to-compact external replacement needs one positive atom")
    _positive_atom(logs[0], "parent long-to-compact external atom")
    return rows[0]


def _trace_chain() -> dict[str, Any]:
    return {
        "pieces": [
            {"name": "second_switch", "lower": "second_switch_native0", "upper": "second_switch_native1"},
            {"name": "post_power", "lower": "post_power_native0", "upper": "post_power_native1"},
            {"name": "long_reshape", "lower": "long_reshape_native0", "upper": "long_reshape_native1"},
        ],
        "aliases": {
            "second_switch_native1": "post_power_native0",
            "post_power_native1": "long_reshape_native0",
        },
        "cancelled_interfaces_explicit": [
            "second_switch native1 = post_power native0",
            "post_power native1 = long_reshape native0",
        ],
        "external_endpoints": {"lower": "second_switch native0", "upper": "long native1 / parent compact2Rc chain"},
        "function": "D(second_switch native0)",
        "terminal_upper_zero_by_source_identity": True,
        "delta_definition": "D=R*(Qy_modified-Qy_leading)",
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": True,
    }


def _source_binding(owner: Any, parent: Mapping[str, Any], proofs: Mapping[str, Mapping[str, Any]], graph: Mapping[str, Any], guard: Mapping[str, Any], ast_binding: Mapping[str, Any]) -> dict[str, Any]:
    parent_source = parent.get("numerical_source_binding") or parent.get("source_binding")
    parent_source = _mapping(parent_source, "parent numerical source binding")
    inherited_keys = (
        "delta_definition", "physical_radius_power", "ordinary_Z_orders", "physical_R_once",
        "leading_subtraction_retained", "cap_is_not_signed_function_value", "all_other_errors_retained",
        "source_functions_values_or_caps_not_selected", "accepted_ancestor_constructors_not_run",
    )
    source = {key: copy.deepcopy(parent_source[key]) for key in inherited_keys if key in parent_source}
    source.update(
        {
            "passed": True,
            "function": "D=R*(Qy_modified-Qy_leading); second_switch native0 after compact2Rc source-identity zero",
            "fallback_source_function": "D(second_switch native0)",
            "delta_definition": "R*(Qy_modified-Qy_leading)",
            "physical_radius_power": 1,
            "ordinary_Z_orders": [0, 1],
            "physical_R_once": True,
            "leading_subtraction_retained": True,
            "cap_is_not_signed_function_value": True,
            "source_functions_values_or_caps_not_selected": True,
            "accepted_ancestor_constructors_not_run": True,
            "actual_second0_graph_envelope_binding": copy.deepcopy(dict(graph)),
            "second_switch_native_guard": copy.deepcopy(dict(guard)),
            "second_switch_pressure_rows_source_binding": copy.deepcopy(dict(ast_binding)),
            "second_power_source_join": copy.deepcopy(dict(proofs["second_power"])),
            "power_long_source_join": copy.deepcopy(dict(proofs["power_long"])),
            "second_power_source_join_proved": True,
            "power_long_source_join_proved": True,
            "second_switch_native0_scope": {
                "chart": "second_switch",
                "coordinate_node": SECOND_COORDINATE_NODE,
                "endpoint_node": SECOND_ENDPOINT_NODE,
                "native_coordinate_value": "0",
                "Qy_function_nodes": list(SECOND_QY_NODES),
                "RQy_function_nodes": list(SECOND_RQY_NODES),
                "ordinary_Z_orders": [0, 1],
                "physical_R_once": True,
            },
            "source_family": copy.deepcopy(parent.get("source_family")),
            "parent_tail_report": {"path": PARENT, "sha256": PARENT_SHA},
            "parent_source_binding_sha256": _digest(parent_source),
            "raw_viscous_indices": {
                OLD_SECOND_TERM: SECOND_VISCOUS_INDEX,
                OLD_POST_POWER_TERM: POST_POWER_VISCOUS_INDEX,
            },
            "no_assumed_second_power_or_power_long_seam": False,
            "raw_prefix_source_kind_preserved": True,
            "generic_nodes_reused_only_under_explicit_scope": True,
            "generic_history_bounds_inherited_not_recertified": True,
        }
    )
    terminal = parent.get("terminal_D_zero_certificate")
    if not isinstance(terminal, Mapping):
        terminal = parent_source.get("terminal_D_zero_certificate")
    if not isinstance(terminal, Mapping) or terminal.get("passed") is not True:
        raise ValueError("Parent terminal D-zero source certificate required")
    source["terminal_D_zero_certificate"] = copy.deepcopy(dict(terminal))
    source["terminal_D_zero_certificate_sha256"] = _digest(terminal)
    source["terminal_D_zero_used"] = True
    source["parent_compact2Rc_D_zero_by_source_identity"] = True
    source["old_compact_external_is_provenance_only"] = True
    source["external_endpoints_not_zeroed"] = False
    source["surviving_second_switch_native0_not_zeroed"] = True
    source["endpoint_zero_only_by_source_identity"] = True
    source["positive_union_uses_own_second_switch0_only"] = True
    source["source_functions_not_cap_values"] = True
    source["old_compact_external_removed_from_active_view"] = True
    source["new_second_switch_external_once"] = True
    return source


def _evaluate_sample(owner: Any, parent: Mapping[str, Any], report: Mapping[str, Any], source: Mapping[str, Any], proofs: Mapping[str, Mapping[str, Any]], graph: Mapping[str, Any], parent_sample: Mapping[str, Any]) -> dict[str, Any]:
    Z = str(parent_sample.get("Z"))
    raw = parent_sample["raw_m15_C0_C1_terms"]
    previous = parent_sample["correlated_reduced_m15_C0_C1_terms"]
    packet = _validate_packets(owner, parent, report, parent_sample)
    old_external = _parent_external_term(previous)
    selected = packet["selected_integer"]
    inlet_atoms, inlet_witness = _second0_pressure_witness(owner, report, Z, selected, graph)
    # The parent tail source certificate proves the compact 2 Rc endpoint is
    # exactly zero.  Its old external atom is retained as provenance only;
    # activating it here would double count the already cancelled long leg.
    positive_inputs = [owner.decode(atom) for atom in inlet_atoms if atom is not None]
    bound = owner.positive_union(positive_inputs)
    replacement = copy.deepcopy(old_external)
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True,
        "function": "D(second_switch native0)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": False,
        "surviving_second_switch_native0_not_zeroed": True,
        "endpoint_zero_only_by_source_identity": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_second_switch_native0_only": True,
        "independent_own_second_switch0_budget": True,
        "parent_long_compact2Rc_external_budget_retained": False,
        "parent_long_compact2Rc_external_provenance_retained": True,
        "parent_compact2Rc_D_zero_by_source_identity": True,
        "parent_terminal_D_zero_certificate_used": True,
        "second_power_source_join_proved": True,
        "power_long_source_join_proved": True,
        "source_functions_values_or_caps_not_selected": True,
        "source_seam_proofs": {name: copy.deepcopy(dict(value)) for name, value in proofs.items()},
    }

    reduced: list[dict[str, Any]] = []
    removed: list[str] = []
    for term in previous:
        name = term.get("name")
        if name in (OLD_SECOND_TERM, OLD_POST_POWER_TERM):
            removed.append(name)
            continue
        if name == OLD_EXTERNAL_TERM:
            removed.append(name)
            reduced.append(copy.deepcopy(replacement))
        else:
            reduced.append(copy.deepcopy(term))
    if sorted(removed) != sorted((OLD_SECOND_TERM, OLD_POST_POWER_TERM, OLD_EXTERNAL_TERM)) or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 187-term view required")
    old_unchanged = [term for term in previous if term.get("name") not in {OLD_SECOND_TERM, OLD_POST_POWER_TERM, OLD_EXTERNAL_TERM}]
    new_unchanged = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(old_unchanged) != UNCHANGED_COUNT or new_unchanged != old_unchanged:
        raise ValueError("Other 186 correlated sectors must remain verbatim")
    if sum(term.get("name") == NEW_EXTERNAL_TERM for term in reduced) != 1:
        raise ValueError("Exactly one upstream external replacement required")

    target = copy.deepcopy(parent_sample.get("independent_target"))
    if target is None:
        target = copy.deepcopy(parent.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")
    return {
        "Z": Z,
        "source_family": getattr(owner, "family", parent.get("source_family")),
        "raw_m15_C0_C1_terms": copy.deepcopy(raw),
        "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT,
        "previous_reduced_terms": PREVIOUS_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT,
        "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": _trace_chain(),
        "source_seam_proofs": {name: copy.deepcopy(dict(value)) for name, value in proofs.items()},
        "actual_signed_function_ports": [
            {
                "chart": "second_switch",
                "side": "lower",
                "coordinate_node": SECOND_COORDINATE_NODE,
                "endpoint_node": SECOND_ENDPOINT_NODE,
                "endpoint_scope": {"coordinate_node": SECOND_COORDINATE_NODE, "endpoint_node": SECOND_ENDPOINT_NODE},
                "modified_Qy_nodes": list(SECOND_QY_NODES),
                "modified_RQy_nodes": list(SECOND_RQY_NODES),
                "definition": "D=R*(Qy_modified-Qy_leading)",
                "ordinary_Z_orders": [0, 1],
                "physical_R_once": True,
                "graph_scope_is_explicit": True,
            }
        ],
        "external_endpoint_budget_witness": {
            "own_second_switch_native0": inlet_witness,
            "parent_long_compact2Rc_external_term": copy.deepcopy(old_external),
            "parent_long_compact2Rc_external_term_active": False,
            "positive_union": owner.record(bound),
            "positive_union_inputs": copy.deepcopy(inlet_atoms),
            "source_binding": copy.deepcopy(dict(source)),
            "full_second_switch_raw68_aggregate_used": False,
            "raw68_and_raw87_retained_as_provenance_only": True,
            "no_cap_subtraction_value_midpoint_or_overlap": True,
        },
        "old_removed_upstream_viscous_terms": [
            copy.deepcopy(raw[SECOND_VISCOUS_INDEX]),
            copy.deepcopy(raw[POST_POWER_VISCOUS_INDEX]),
        ],
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_external),
        "replacement_outer_endpoints": copy.deepcopy(replacement),
        "packet_prefix_witness": copy.deepcopy(packet["exact_prefix_witness"]),
        "selected_N_evidence": copy.deepcopy(selected),
        "independent_target": target,
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True,
            "other_186_verbatim": True,
            "independent_target_unchanged": True,
            "reserved_retained_once": sum(term.get("name") == "reserved_n1_own_swirl_once" for term in reduced) == 1,
            "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True,
            "second_power_and_power_long_source_joins_proved": True,
            "second_and_post_power_viscous_sectors_removed_once": True,
            "parent_external_term_replaced_once": True,
            "second_switch_native0_budget_owned": True,
            "source_functions_not_cap_values": True,
            "second_power_source_join_proved": True,
            "power_long_source_join_proved": True,
            "old_compact_external_removed": True,
            "new_second_switch_external_once": True,
            "previous_189_to_new_187": True,
            "raw68_full_aggregate_not_used_as_new_bound": True,
            "raw87_full_aggregate_not_used_as_new_bound": True,
            "no_error_cap_as_value": True,
            "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True,
            "source_bound_signed_trace_contraction": True,
        },
        "limits": {
            "upstream_internal_trace_resolved": True,
            "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "mean_oscillatory_correction": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


def build_upstream_numeric(owner: Any, seam_proofs: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved second-switch to long 187-term contraction packet."""
    proofs = _proofs(seam_proofs)
    parent = _parent_report(owner)
    report = _r110_report(owner)
    family = getattr(owner, "family", parent.get("source_family"))
    if parent.get("source_family") != family:
        raise ValueError("Parent source family changed")
    ast_binding = _pin_r110_ast(owner)
    guard = _second_native_guard(owner, proofs)
    graph = _graph_binding(owner, parent, proofs, guard)
    source = _source_binding(owner, parent, proofs, graph, guard, ast_binding)
    samples: list[dict[str, Any]] = []
    for Z in SAVED_Z:
        parent_sample = _find(parent.get("samples"), label="latest tail parent sample", Z=Z)
        samples.append(_evaluate_sample(owner, parent, report, source, proofs, graph, parent_sample))
    selected = [sample["selected_N_evidence"] for sample in samples]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All upstream samples must use one selected N")
    source.update(
        {
            "source_proof_sha256": {name: _digest(proof) for name, proof in proofs.items()},
            "r110_report": {"path": R110_REPORT, "sha256": R110_REPORT_SHA},
            "r110_packet_prefix_counts": {
                "second_switch": SECOND_INCOMING_COUNT,
                "post_power": POST_POWER_INCOMING_COUNT,
            },
            "actual_second0_Qy_nodes": list(SECOND_QY_NODES),
            "actual_second0_RQy_nodes": list(SECOND_RQY_NODES),
            "actual_second0_endpoint_scope": {"coordinate_node": SECOND_COORDINATE_NODE, "endpoint_node": SECOND_ENDPOINT_NODE},
            "selected_N_evidence": copy.deepcopy(selected[0]),
        }
    )
    return {
        "source_binding": source,
        "samples": samples,
        "source_seam_proofs": {"passed": True, "joins": copy.deepcopy(proofs)},
        "contract": {
            "raw_213_retained": True,
            "previous_189_to_new_187": True,
            "unchanged_other_186": True,
            "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True,
            "cap_is_not_signed_function_value": True,
            "second_power_source_join_required": True,
            "power_long_source_join_required": True,
            "second_switch_native0_budget_required": True,
            "source_functions_not_cap_values": True,
            "second_power_source_join_proved": True,
            "power_long_source_join_proved": True,
            "old_compact_external_removed": True,
            "new_second_switch_external_once": True,
            "raw68_raw87_provenance_only": True,
            "compact2Rc_terminal_zero_by_source_identity": True,
        },
        "limits": {
            "numeric_Rp_closure": False,
            "global_all_five_closure": False,
            "full_goal_complete": False,
        },
    }


__all__ = ["build_upstream_numeric"]
