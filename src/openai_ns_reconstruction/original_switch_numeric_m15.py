"""Contract the saved first-switch to second-switch viscous chain.

This is a saved-packet adapter.  It consumes the accepted upstream packet,
the pinned first-switch and R110 packets, and the caller's actual
``first_switch native 1 -> second_switch native 0`` source seam.  It does
not construct an ancestor or evaluate an original producer.  The raw 213
term state and independent target remain immutable; the correlated view
removes the first-switch finite-N sector and the previous external
replacement, retaining a positive enclosure for first-switch native 0.
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
PARENT_REDUCED_COUNT = 187
NEW_REDUCED_COUNT = 186
UNCHANGED_COUNT = 185

PARENT = PREFIX + "compliant_current_original_upstream_m15_numeric.json"
PARENT_SHA = "6028d49ac9ec1b94ae967b15d698366061b4a15b9e0ddef9b45be2d3971e355e"
FIRST_REPORT = PREFIX + "compliant_current_original_first_switch_m15_numeric.json"
FIRST_REPORT_SHA = "fd427c4ba14a541f3ec042b644c2e9c65b66077c8350f6807a96c957a36bf83c"
FIRST_PRESSURE_REPORT = PREFIX + "compliant_current_original_first_switch_pressure_rows.json"
FIRST_PRESSURE_REPORT_SHA = "6bcbd35931ddc35a3d743ae0a46c879f9479d4c6d77aba43ea3412d16881ae33"
R110_REPORT = PREFIX + "compliant_current_original_r110_m15_numeric.json"
R110_REPORT_SHA = "7226a0546b338c6834d035f6daa86707091c4c1e2901f631ea8a2351231bada1"

FIRST_M15_SOURCE = "src/openai_ns_reconstruction/original_first_switch_m15.py"
FIRST_M15_SOURCE_SHA = "0dffed08283abdfbc6b358f7bb59159dd5daba7b6c2be4b9e602bda510c27baf"
FIRST_PRESSURE_SOURCE = "src/openai_ns_reconstruction/original_first_switch_pressure_rows.py"
FIRST_PRESSURE_SOURCE_SHA = "fa175264d587321b81b58002ce1afe35b6de4bfc440e57698e587964190c7b44"
R110_M15_SOURCE = "src/openai_ns_reconstruction/original_r110_m15.py"
R110_M15_SOURCE_SHA = "ca6e74ff8ca4eaaf48e2693328c1b922ba8b683b28b87712941541f497cec20d"
R110_PRESSURE_SOURCE = "src/openai_ns_reconstruction/original_r110_pressure_rows.py"
R110_PRESSURE_SOURCE_SHA = "9f60fc44291c16f0906ce73b19a414645fa1b2a2af1bec3ca2d3fc76e77761ef"

FIRST_SWITCH_VISCOUS_INDEX = 50
FIRST_INCOMING_COUNT = 28
FIRST_TERM_COUNT = 51
SECOND_INCOMING_COUNT = 51
SECOND_TERM_COUNT = 69
OLD_FIRST_TERM = "first_switch_finite_N_viscous_error"
OLD_EXTERNAL_TERM = "second_switch_long_finite_N_viscous_external_endpoints"
NEW_EXTERNAL_TERM = "first_switch_second_finite_N_viscous_external_endpoints"

FIRST_Q0 = {"kind": "exact_decimal_first_switch_phase", "value": "0"}
FIRST_Q1 = {"kind": "exact_decimal_first_switch_phase", "value": "1"}
SECOND_Q0 = {"kind": "exact_decimal_second_switch_coordinate", "value": "0"}
SECOND_Q1 = {"kind": "exact_decimal_second_switch_coordinate", "value": "1"}

# These are graph handles.  Endpoint values remain separate in seam scopes.
FIRST_COORDINATE_NODE = 73190
SECOND_COORDINATE_NODE = 75949
FIRST_ENDPOINT_NODE = 0
FIRST_UPPER_ENDPOINT_NODE = 1
SECOND_ENDPOINT_NODE = 0
SECOND_UPPER_ENDPOINT_NODE = 1
FIRST_RADIUS_NODE = 85905
FIRST_JACOBIAN_NODE = 9
FIRST_QY_NODES = (148268, 148272)
FIRST_RQY_NODES = (148427, 148428)


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
    known_values = (hashes.get(relative), hashes.get(basename), code_hashes.get(relative), code_hashes.get(basename))
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
    return _candidate_report(owner, ("parent", "parent_report"), pinned)


def _first_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, FIRST_REPORT, FIRST_REPORT_SHA)
    return _candidate_report(owner, ("first", "first_switch", "first_m15", "first_switch_report"), pinned)


def _first_pressure_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, FIRST_PRESSURE_REPORT, FIRST_PRESSURE_REPORT_SHA)
    return _candidate_report(owner, ("first_pressure", "first_switch_pressure", "pressure_rows"), pinned)


def _r110_report(owner: Any) -> Mapping[str, Any]:
    pinned = _load_json(owner, R110_REPORT, R110_REPORT_SHA)
    return _candidate_report(owner, ("r110", "r110_report", "second_report", "second"), pinned)


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
    if atom.get("cap_as_value", False) or atom.get("signed_function_value", False):
        raise ValueError(f"{label} is a cap/value rather than an atom")
    return atom


def _unwrapped(value: Any, key: str) -> Mapping[str, Any]:
    mapping = _mapping(value, key)
    if key in mapping and isinstance(mapping.get(key), Mapping):
        return _mapping(mapping[key], key)
    return mapping


def _proofs(seam_proofs: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    if not isinstance(seam_proofs, Mapping) or seam_proofs.get("passed") is not True:
        raise ValueError("First-switch source seam proofs must pass")
    candidate = seam_proofs.get("joins")
    if not isinstance(candidate, Mapping):
        candidate = seam_proofs.get("seams")
    if not isinstance(candidate, Mapping):
        candidate = seam_proofs
    proof = candidate.get("first_second")
    if not isinstance(proof, Mapping) or proof.get("passed") is not True:
        raise ValueError("Accepted first_second source proof must pass")
    expected = {
        "left_chart": "first_switch", "right_chart": "second_switch",
        "left_endpoint": 1, "right_endpoint": 0,
        "left_endpoint_node": 1, "right_endpoint_node": 0,
        "left_coordinate": FIRST_COORDINATE_NODE, "right_coordinate": SECOND_COORDINATE_NODE,
        "left_coordinate_node": FIRST_COORDINATE_NODE, "right_coordinate_node": SECOND_COORDINATE_NODE,
    }
    for key, value in expected.items():
        if key in proof and proof.get(key) != value:
            raise ValueError(f"first_second source proof has the wrong {key}")
    for side, coordinate, endpoint in (("left", FIRST_COORDINATE_NODE, 1), ("right", SECOND_COORDINATE_NODE, 0)):
        scopes = [proof[key] for key in (f"{side}_scope", f"{side}_endpoint_scope") if isinstance(proof.get(key), Mapping)]
        for scope in scopes:
            if "coordinate_node" in scope and scope.get("coordinate_node") != coordinate:
                raise ValueError(f"first_second {side} coordinate scope changed")
            if "native_coordinate_node" in scope and scope.get("native_coordinate_node") != coordinate:
                raise ValueError(f"first_second {side} native coordinate scope changed")
            if "endpoint_node" in scope and scope.get("endpoint_node") != endpoint:
                raise ValueError(f"first_second {side} endpoint scope changed")
            for key in ("native_query_node", "query_node"):
                if key in scope and scope.get(key) != endpoint:
                    raise ValueError(f"first_second {side} query node changed")
            for key in ("native_query_value", "query_value"):
                if key in scope and str(scope.get(key)) != str(endpoint):
                    raise ValueError(f"first_second {side} query value changed")
    for key in ("source_values_selected", "numeric_caps_used", "source_value_selected", "ancestor_replayed"):
        if proof.get(key) is True:
            raise ValueError("first_second source proof selected a value or cap")
    return {"first_second": copy.deepcopy(dict(proof))}


def _chart_bindings(owner: Any) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    loop = getattr(owner, "l", None) or getattr(owner, "loop", None)
    modified = getattr(owner, "m", None) or getattr(owner, "modified", None)
    loop_map = _unwrapped(loop, "prescribed_loop_C5_bindings")
    modified_map = _unwrapped(modified, "modified_Q4_pressure_C3_bindings")
    return loop_map, modified_map


def _first_native_guard(owner: Any, proof: Mapping[str, Any]) -> dict[str, Any]:
    guards = getattr(owner, "native_guards", None)
    if not isinstance(guards, Mapping):
        guards = getattr(owner, "native_guard", None)
    if not isinstance(guards, Mapping):
        guards = getattr(owner, "native", None)
    charts = guards.get("charts") if isinstance(guards, Mapping) else None
    if not isinstance(charts, Mapping):
        charts = guards.get("chart_guards") if isinstance(guards, Mapping) else None
    row = charts.get("first_switch") if isinstance(charts, Mapping) else None
    if not isinstance(guards, Mapping) or guards.get("passed") is not True or not isinstance(row, Mapping) or row.get("passed") is not True:
        candidate = proof.get("first_native_guard", proof.get("native_guard"))
        if isinstance(candidate, Mapping) and candidate.get("passed") is True:
            row = candidate
            guards = candidate
        else:
            raise ValueError("Passed first_switch physical native guard required")
    loop, _ = _chart_bindings(owner)
    chart = _mapping(loop.get("charts", {}).get("first_switch"), "first_switch loop chart")
    if row.get("coordinate") not in (None, chart.get("coordinate")) or row.get("radius") not in (None, chart.get("radius")):
        raise ValueError("First-switch native guard is bound to a foreign chart")
    if row.get("native_J") not in (None, chart.get("physical_measure")):
        raise ValueError("First-switch native Jacobian guard changed")
    return {
        "passed": True, "chart": "first_switch", "guard": copy.deepcopy(dict(row)),
        "source_guard_packet": copy.deepcopy(dict(guards)), "ordinary_Z_orders": [0, 1],
        "physical_RQy_factor_once": True, "physical_Uy_Pstar_factor_once": True,
    }


def _parent_second_graph(parent: Mapping[str, Any]) -> Mapping[str, Any]:
    source = parent.get("numerical_source_binding") or parent.get("source_binding")
    if isinstance(source, Mapping) and isinstance(source.get("actual_second0_graph_envelope_binding"), Mapping):
        return _mapping(source["actual_second0_graph_envelope_binding"], "parent second graph envelope")
    for sample in parent.get("samples", []):
        if not isinstance(sample, Mapping):
            continue
        budget = sample.get("external_endpoint_budget_witness")
        if isinstance(budget, Mapping):
            own = budget.get("own_second_switch_native0")
            if isinstance(own, Mapping) and isinstance(own.get("actual_second0_graph_envelope_binding"), Mapping):
                return _mapping(own["actual_second0_graph_envelope_binding"], "parent second graph envelope")
    raise ValueError("Parent second-switch graph envelope certificate required")


def _graph_binding(owner: Any, parent: Mapping[str, Any], guard: Mapping[str, Any], proof: Mapping[str, Any]) -> dict[str, Any]:
    prior = _parent_second_graph(parent)
    if prior.get("passed") is not True:
        raise ValueError("Parent universal graph envelope certificate must pass")
    inherited_keys = (
        "exact_numeric_Qy_coefficients", "exact_coefficient_geometry_statements",
        "ordinary_Z_product_rule", "actual_archived_pressure_formula_certificate",
        "actual_archived_radial_formula_certificate",
    )
    for key in inherited_keys:
        if key not in prior:
            raise ValueError(f"Parent universal Qy certificate missing {key}")
    loop, modified = _chart_bindings(owner)
    charts = _mapping(loop.get("charts"), "loop charts")
    loop_chart = _mapping(charts.get("first_switch"), "first-switch loop chart")
    mod_charts = _mapping(modified.get("charts"), "modified charts")
    modified_chart = _mapping(mod_charts.get("first_switch"), "first-switch modified chart")
    coordinate = loop_chart.get("coordinate")
    domain = loop_chart.get("domain")
    radius = loop_chart.get("radius")
    physical_measure = loop_chart.get("physical_measure")
    if (coordinate, domain, radius, physical_measure) != (FIRST_COORDINATE_NODE, [0, 1], FIRST_RADIUS_NODE, FIRST_JACOBIAN_NODE):
        raise ValueError("First-switch native coordinate/domain/radius/Jacobian changed")
    physical = _mapping(modified_chart.get("physical"), "first-switch modified physical")
    finite = _mapping(modified_chart.get("finite_boundary"), "first-switch finite boundary")
    qy_nodes = list(physical.get("Q_y", [])[:2])
    rqy_nodes = list(finite.get("R_Q_y", [])[:2])
    if qy_nodes != list(FIRST_QY_NODES) or rqy_nodes != list(FIRST_RQY_NODES):
        raise ValueError("First-switch Qy/RQy graph rows changed")
    nodes = getattr(owner, "nodes", None)
    if not isinstance(nodes, (list, tuple, Mapping)):
        raise ValueError("Actual switch graph nodes required")
    graph_rows: list[dict[str, Any]] = []
    fingerprints: dict[int, str] = {}
    for node in rqy_nodes:
        try:
            row = nodes[node] if not isinstance(nodes, Mapping) else nodes[node]
        except (IndexError, KeyError, TypeError):
            raise ValueError("First-switch RQy graph node missing") from None
        if not isinstance(row, Mapping) or row.get("operation") != "product":
            raise ValueError("First-switch RQy must be an actual product")
        arguments = row.get("arguments", [])
        if not isinstance(arguments, list) or arguments.count(radius) != 1:
            raise ValueError("First-switch RQy must contain one physical R factor")
    try:
        from .original_active_viscous_m15 import ActiveEndpointFingerprint
        fingerprint = ActiveEndpointFingerprint(nodes, dict(coordinate=coordinate, loop=charts))
        scope = {coordinate: fingerprint.at(domain[0], {})}
        for node in (*qy_nodes, *rqy_nodes):
            fingerprints[node] = fingerprint.at(node, scope)[0]
    except Exception as exc:
        raise ValueError("First-switch graph fingerprints could not be bound") from exc
    supplied_rows = None
    for key in ("first_switch0_graph_function_rows", "first0_graph_function_rows", "graph_function_rows"):
        value = proof.get(key)
        if isinstance(value, list) and len(value) >= 2:
            supplied_rows = value[:2]
            break
    graph_rows = []
    for order in (0, 1):
        row: dict[str, Any] = {
            "ordinary_Z_order": order,
            "Qy_function_node": qy_nodes[order], "RQy_function_node": rqy_nodes[order],
            "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": FIRST_ENDPOINT_NODE},
            "physical_radius_node": radius, "physical_R_explicit_once": True,
        }
        if qy_nodes[order] in fingerprints:
            row["Qy_at_first0_sha256"] = fingerprints[qy_nodes[order]]
        if rqy_nodes[order] in fingerprints:
            row["RQy_at_first0_sha256"] = fingerprints[rqy_nodes[order]]
        if supplied_rows is not None:
            supplied = _mapping(supplied_rows[order], "first graph function row")
            for key in ("Qy_at_first0_sha256", "RQy_at_first0_sha256"):
                if key in supplied and key in row and supplied[key] != row[key]:
                    raise ValueError("First-switch graph fingerprint changed")
            row["source_proof_row"] = copy.deepcopy(dict(supplied))
        graph_rows.append(row)
    return {
        "passed": True,
        "method": "source_bound_linear_Qy_envelope_at_actual_first_switch_native0",
        "inherited_second0_envelope_sha256": _digest(prior),
        **{key: copy.deepcopy(prior[key]) for key in inherited_keys},
        "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": FIRST_ENDPOINT_NODE},
        "graph_function_rows": graph_rows,
        "chart": "first_switch", "derivative_chart": "switch_first",
        "physical_radius_node": radius, "native_Jacobian_node": physical_measure,
        "ordinary_Z_orders": [0, 1], "physical_RQy_factor_once": True,
        "physical_Uy_Pstar_factor_once": True, "first_native_guard": copy.deepcopy(dict(guard)),
        "first_second_source_scope": copy.deepcopy(dict(proof)),
        "generic_nodes_reused_only_under_explicit_scope": True,
        "generic_history_bounds_inherited_not_recertified": True,
        "source_function_values_or_caps_not_selected": True,
        "actual_first_switch0_scope_only": True,
        "no_unproved_whole_chart_alias": True,
    }


def _pin_ast_sources(owner: Any) -> dict[str, Any]:
    pins = {
        "first_pressure_rows_source": _pin_file(owner, FIRST_PRESSURE_SOURCE, FIRST_PRESSURE_SOURCE_SHA),
        "first_m15_source": _pin_file(owner, FIRST_M15_SOURCE, FIRST_M15_SOURCE_SHA),
        "r110_pressure_rows_source": _pin_file(owner, R110_PRESSURE_SOURCE, R110_PRESSURE_SOURCE_SHA),
        "r110_m15_source": _pin_file(owner, R110_M15_SOURCE, R110_M15_SOURCE_SHA),
    }
    checks = []
    specs = ((FIRST_M15_SOURCE, "OriginalFirstSwitchM15", {"background", "correction_densities", "sum_logs"}),
             (R110_M15_SOURCE, "OriginalR110M15", {"background", "correction_densities", "sum_logs"}))
    for relative, class_name, required in specs:
        path = _path(owner, relative)
        tree = ast.parse(path.read_bytes(), filename=str(path))
        cls = next((node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == class_name), None)
        if cls is None:
            raise ValueError(f"{class_name} class missing from pinned source")
        method = next((node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == "finite_errors"), None)
        if method is None:
            raise ValueError(f"{class_name}.finite_errors method missing")
        calls = {node.func.attr for node in ast.walk(method) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)}
        if not required.issubset(calls):
            raise ValueError(f"{class_name}.finite_errors lacks complete finite-N source rows")
        checks.append({"class": class_name, "method": f"{class_name}.finite_errors",
                       "method_ast_sha256": hashlib.sha256(ast.dump(method, include_attributes=False).encode("utf-8")).hexdigest(),
                       "required_calls": sorted(required), "source_not_executed": True})
    return {"passed": True, "source_pins": pins, "finite_errors_AST": checks,
            "source_not_executed": True, "Qy_coefficient_geometry_source_bound": True}


def _validate_packets(parent_sample: Mapping[str, Any], first: Mapping[str, Any], r110: Mapping[str, Any], family: Any) -> dict[str, Any]:
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Parent raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PARENT_REDUCED_COUNT:
        raise ValueError("Parent upstream 187-term view required")
    _term_names(raw, "parent raw terms")
    _term_names(previous, "parent reduced terms")
    if raw[FIRST_SWITCH_VISCOUS_INDEX].get("name") != OLD_FIRST_TERM:
        raise ValueError("Parent raw first-switch viscous index changed")
    if sum(term.get("name") == OLD_EXTERNAL_TERM for term in previous) != 1:
        raise ValueError("Parent second-switch external replacement required")
    selected_parent = _mapping(parent_sample.get("selected_N_evidence"), "parent selected-N evidence")
    Z = str(parent_sample.get("Z"))
    first_rows = {
        "first0": _find(first.get("samples"), label="first-switch native0", Z=Z, coordinate=FIRST_Q0),
        "first1": _find(first.get("samples"), label="first-switch native1", Z=Z, coordinate=FIRST_Q1),
    }
    for name, row in first_rows.items():
        if row.get("source_family") != family:
            raise ValueError(f"First-switch {name} source family changed")
        if row.get("actual_selected_integer") != selected_parent:
            raise ValueError(f"First-switch {name} selected N changed")
        terms, incoming = row.get("m15_C0_C1_terms"), row.get("incoming_macro_m15_terms")
        if not isinstance(terms, list) or len(terms) != FIRST_TERM_COUNT or not isinstance(incoming, list) or len(incoming) != FIRST_INCOMING_COUNT:
            raise ValueError(f"First-switch {name} packet count changed")
        _term_names(terms, f"first-switch {name} terms")
        _term_names(incoming, f"first-switch {name} incoming")
        _exact(f"{name} incoming versus parent raw prefix", incoming, raw[:FIRST_INCOMING_COUNT])
        _exact(f"{name} terms prefix versus incoming", terms[:FIRST_INCOMING_COUNT], incoming)
    if first_rows["first0"].get("exact_empty_increment") is not True or first_rows["first1"].get("exact_empty_increment") is not False:
        raise ValueError("First-switch endpoint emptiness changed")
    second_rows = {
        "second0": _find(r110.get("samples"), label="R110 second-switch native0", Z=Z, chart="second_switch", coordinate=SECOND_Q0),
        "second1": _find(r110.get("samples"), label="R110 second-switch native1", Z=Z, chart="second_switch", coordinate=SECOND_Q1),
    }
    for name, row in second_rows.items():
        if row.get("source_family") != family or row.get("actual_selected_integer") != selected_parent:
            raise ValueError(f"R110 {name} family or selected N changed")
        terms, incoming = row.get("m15_C0_C1_terms"), row.get("incoming_m15_terms")
        if not isinstance(terms, list) or len(terms) != SECOND_TERM_COUNT or not isinstance(incoming, list) or len(incoming) != SECOND_INCOMING_COUNT:
            raise ValueError(f"R110 {name} packet count changed")
        _term_names(terms, f"R110 {name} terms")
        _term_names(incoming, f"R110 {name} incoming")
        _exact(f"{name} incoming versus parent raw prefix", incoming, raw[:SECOND_INCOMING_COUNT])
        _exact(f"{name} terms prefix versus incoming", terms[:SECOND_INCOMING_COUNT], incoming)
    if second_rows["second0"].get("exact_empty_increment") is not True or second_rows["second1"].get("exact_empty_increment") is not False:
        raise ValueError("R110 endpoint emptiness changed")
    _exact("first-switch native1 complete packet versus second-switch native0 incoming", first_rows["first1"]["m15_C0_C1_terms"], second_rows["second0"]["incoming_m15_terms"])
    _exact("first-switch native1 complete packet versus second-switch native1 incoming", first_rows["first1"]["m15_C0_C1_terms"], second_rows["second1"]["incoming_m15_terms"])
    return {
        "Z": Z, "selected_integer": copy.deepcopy(dict(selected_parent)), "first_rows": first_rows, "second_rows": second_rows,
        "exact_prefix_witness": {
            "passed": True, "first_prefix_count": FIRST_INCOMING_COUNT, "first_term_count": FIRST_TERM_COUNT,
            "second_prefix_count": SECOND_INCOMING_COUNT, "second_term_count": SECOND_TERM_COUNT,
            "raw_prefix_equals_first0_incoming": True, "raw_prefix_equals_first1_incoming": True,
            "raw_prefix_equals_second0_incoming": True, "raw_prefix_equals_second1_incoming": True,
            "first0_terms_prefix_equals_incoming": True, "first1_terms_prefix_equals_incoming": True,
            "second0_terms_prefix_equals_incoming": True, "second1_terms_prefix_equals_incoming": True,
            "first1_complete_equals_second0_incoming": True, "first1_complete_equals_second1_incoming": True,
            "full_dictionary_equality": True, "raw_prefix_sha256": _digest(raw[:FIRST_INCOMING_COUNT]),
            "first1_complete_sha256": _digest(first_rows["first1"]["m15_C0_C1_terms"]),
        },
    }


def _first0_pressure_witness(owner: Any, report: Mapping[str, Any], Z: str, selected: Mapping[str, Any], graph: Mapping[str, Any]) -> tuple[list[Any], dict[str, Any]]:
    row = _find(report.get("samples"), label="first-switch native0 pressure row", Z=Z, coordinate=FIRST_Q0)
    family = getattr(owner, "family", report.get("source_family"))
    if row.get("source_family") != family or row.get("selected_integer") != selected:
        raise ValueError("First-switch native0 pressure family or selected N changed")
    cover = _mapping(row.get("source_cover"), "first-switch source cover")
    if cover.get("generic_derivative_chart") != "switch_first":
        raise ValueError("First-switch source derivative chart changed")
    histories = cover.get("cumulative_history_charts")
    if not isinstance(histories, list) or not histories or histories[-1] != "first_switch":
        raise ValueError("First-switch cumulative source history changed")
    if cover.get("all_axial_frames_retained") is not True or cover.get("native_source_domain") != ["0", "1"]:
        raise ValueError("First-switch source cover must retain all native axial frames")
    contract = _mapping(row.get("contract"), "first-switch pressure contract")
    for key in ("original_Pstar_once", "finite_R_Qy_and_ordinary_Z1_retained"):
        if contract.get(key) is not True:
            raise ValueError(f"First-switch pressure contract lacks {key}")
    physical_rows = _mapping(row.get("physical_pressure_rows"), "first-switch physical pressure rows")
    qy = _mapping(physical_rows.get("Q_y"), "first-switch Qy rows")
    endpoint = _mapping(row.get("R_Q_y_endpoint"), "first-switch RQy endpoint")
    errors = qy.get("finite_N_absolute_error_log_bounds")
    endpoint_errors = endpoint.get("finite_N_coefficient_error_log_bounds")
    if errors != endpoint_errors or not isinstance(errors, list) or len(errors) != 2:
        raise ValueError("First-switch Qy C0/Z1 finite-N rows required")
    radius = endpoint.get("finite_N_log_physical_factor")
    geometry = _mapping(row.get("geometry"), "first-switch pressure geometry")
    if radius != geometry.get("log_physical_radius"):
        raise ValueError("First-switch physical R factor must be sourced once")
    atoms: list[Any] = []
    for error in errors:
        error_map = _mapping(error, "first-switch Qy error row")
        if error_map.get("exact_zero") is True:
            raise ValueError("First-switch own C0/Z1 envelope unexpectedly empty")
        upper = _mapping(error_map.get("log_absolute_upper"), "first-switch Qy log error")
        atoms.append(owner.record(owner.decode(upper) + owner.decode(radius)))
    own_rows: dict[str, Any] = {}
    for name in ("U", "m", "U_y"):
        value = _mapping(physical_rows.get(name), f"first-switch {name} rows").get("finite_N_absolute_error_log_bounds")
        if not isinstance(value, list):
            raise ValueError("First-switch U/m/Uy source error rows missing")
        own_rows[name] = copy.deepcopy(value)
    slow_fast = copy.deepcopy(row.get("finite_N_axial_radial_error_components"))
    if not isinstance(slow_fast, Mapping):
        raise ValueError("First-switch slow/fast source error components missing")
    witness = {
        "source_report": FIRST_PRESSURE_REPORT, "source_packet_sha256": _digest(row),
        "coordinate": copy.deepcopy(row.get("coordinate")), "selected_integer": copy.deepcopy(dict(selected)),
        "source_cover": copy.deepcopy(dict(cover)), "actual_own_Qy_error_rows": copy.deepcopy(errors),
        "actual_own_log_radius": copy.deepcopy(radius), "actual_linear_error_input_rows": own_rows,
        "own_U_m_Uy_error_rows": copy.deepcopy(own_rows), "actual_slow_fast_error_components": slow_fast,
        "own_slow_fast_components": copy.deepcopy(slow_fast), "radius_weighted_C0_C1_log_bounds": copy.deepcopy(atoms),
        "actual_first0_graph_envelope_binding": copy.deepcopy(dict(graph)),
        "generic_derivative_chart": "switch_first", "first_switch_native0_source_packet": True,
        "physical_R_once": True, "physical_Pstar_already_in_error": True,
        "source_function_values_evaluated": False, "no_caps_used": True,
        "all_axial_frames_retained": True, "ordinary_Z_orders": [0, 1],
    }
    return atoms, witness


def _parent_external_term(previous: list[Mapping[str, Any]]) -> Mapping[str, Any]:
    rows = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    if len(rows) != 1:
        raise ValueError("Parent second-switch external replacement required")
    logs = rows[0].get("log_scale_terms")
    if not isinstance(logs, list) or len(logs) != 1:
        raise ValueError("Parent external replacement needs one positive atom")
    _positive_atom(logs[0], "parent second-switch external atom")
    return rows[0]


def _trace_chain() -> dict[str, Any]:
    return {
        "pieces": [{"name": "first_switch", "lower": "first_switch_native0", "upper": "first_switch_native1"},
                   {"name": "second_switch", "lower": "second_switch_native0", "upper": "second_switch_native1"}],
        "aliases": {"first_switch_native1": "second_switch_native0"},
        "cancelled_interfaces_explicit": ["first_switch native1 = second_switch native0"],
        "external_endpoints": {"lower": "first_switch native0", "upper": "second_switch native0 / parent upstream chain"},
        "function": "D(first_switch native0)", "delta_definition": "D=R*(Qy_modified-Qy_leading)",
        "ordinary_Z_orders": [0, 1], "physical_R_once": True,
        "terminal_upper_zero_by_source_identity": False,
    }


def _source_binding(owner: Any, parent: Mapping[str, Any], proofs: Mapping[str, Mapping[str, Any]], graph: Mapping[str, Any], guard: Mapping[str, Any], ast_binding: Mapping[str, Any]) -> dict[str, Any]:
    parent_source = _mapping(parent.get("numerical_source_binding") or parent.get("source_binding"), "parent numerical source binding")
    inherited = {key: copy.deepcopy(parent_source[key]) for key in (
        "delta_definition", "physical_radius_power", "ordinary_Z_orders", "physical_R_once",
        "leading_subtraction_retained", "cap_is_not_signed_function_value", "all_other_errors_retained",
        "source_functions_values_or_caps_not_selected", "accepted_ancestor_constructors_not_run") if key in parent_source}
    source = {
        **inherited, "passed": True, "function": "D=R*(Qy_modified-Qy_leading); first_switch native0",
        "fallback_source_function": "D(first_switch native0)", "delta_definition": "R*(Qy_modified-Qy_leading)",
        "physical_radius_power": 1, "ordinary_Z_orders": [0, 1], "physical_R_once": True,
        "leading_subtraction_retained": True, "cap_is_not_signed_function_value": True,
        "source_functions_values_or_caps_not_selected": True, "accepted_ancestor_constructors_not_run": True,
        "actual_first0_graph_envelope_binding": copy.deepcopy(dict(graph)),
        "first_switch_native_guard": copy.deepcopy(dict(guard)),
        "first_switch_pressure_rows_source_binding": copy.deepcopy(dict(ast_binding)),
        "first_second_source_join": copy.deepcopy(dict(proofs["first_second"])),
        "first_second_source_join_proved": True,
        "first_switch_native0_scope": {"chart": "first_switch", "coordinate_node": FIRST_COORDINATE_NODE,
            "endpoint_node": FIRST_ENDPOINT_NODE, "native_coordinate_value": "0", "Qy_function_nodes": list(FIRST_QY_NODES),
            "RQy_function_nodes": list(FIRST_RQY_NODES), "ordinary_Z_orders": [0, 1], "physical_R_once": True},
        "source_family": copy.deepcopy(parent.get("source_family")),
        "parent_upstream_report": {"path": PARENT, "sha256": PARENT_SHA},
        "parent_source_binding_sha256": _digest(parent_source),
        "raw_viscous_indices": {OLD_FIRST_TERM: FIRST_SWITCH_VISCOUS_INDEX},
        "first_switch_native1_aliases_second_switch_native0": True,
        "raw_prefix_source_kind_preserved": True,
        "generic_nodes_reused_only_under_explicit_scope": True,
        "generic_history_bounds_inherited_not_recertified": True,
        "source_proof_sha256": {name: _digest(proof) for name, proof in proofs.items()},
        "first_report": {"path": FIRST_REPORT, "sha256": FIRST_REPORT_SHA},
        "first_pressure_report": {"path": FIRST_PRESSURE_REPORT, "sha256": FIRST_PRESSURE_REPORT_SHA},
        "r110_report": {"path": R110_REPORT, "sha256": R110_REPORT_SHA},
        "first_packet_prefix_count": FIRST_INCOMING_COUNT, "second_packet_prefix_count": SECOND_INCOMING_COUNT,
        "actual_first0_Qy_nodes": list(FIRST_QY_NODES), "actual_first0_RQy_nodes": list(FIRST_RQY_NODES),
        "actual_first0_endpoint_scope": {"coordinate_node": FIRST_COORDINATE_NODE, "endpoint_node": FIRST_ENDPOINT_NODE},
    }
    source["old_second_switch_external_provenance_only"] = True
    source["old_second_switch_external_removed_from_active_view"] = True
    source["external_endpoints_not_zeroed"] = True
    source["surviving_first_switch_native0_not_zeroed"] = True
    source["endpoint_zero_only_by_source_identity"] = True
    source["positive_union_uses_own_first_switch0_only"] = True
    source["new_first_switch_external_once"] = True
    return source


def _evaluate_sample(owner: Any, parent: Mapping[str, Any], first: Mapping[str, Any], first_pressure: Mapping[str, Any], r110: Mapping[str, Any], source: Mapping[str, Any], proof: Mapping[str, Mapping[str, Any]], graph: Mapping[str, Any], parent_sample: Mapping[str, Any]) -> dict[str, Any]:
    Z = str(parent_sample.get("Z"))
    raw = parent_sample["raw_m15_C0_C1_terms"]
    previous = parent_sample["correlated_reduced_m15_C0_C1_terms"]
    packet = _validate_packets(parent_sample, first, r110, getattr(owner, "family", parent.get("source_family")))
    old_external = _parent_external_term(previous)
    selected = packet["selected_integer"]
    inlet_atoms, inlet_witness = _first0_pressure_witness(owner, first_pressure, Z, selected, graph)
    positive_inputs = [owner.decode(atom) for atom in inlet_atoms]
    if len(positive_inputs) != 2:
        raise ValueError("Exactly two own first-switch native0 atoms are required")
    bound = owner.positive_union(positive_inputs)
    replacement = copy.deepcopy(old_external)
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True, "function": "D(first_switch native0)",
        "delta_definition": "R*(Qy_modified-Qy_leading)", "ordinary_Z1": True, "physical_R_once": True,
        "external_endpoints_not_zeroed": True, "surviving_first_switch_native0_not_zeroed": True,
        "endpoint_zero_only_by_source_identity": True, "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True, "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_first_switch_native0_only": True, "independent_own_first_switch0_budget": True,
        "parent_second_switch_external_budget_retained": False,
        "parent_second_switch_external_provenance_retained": True,
        "first_second_source_join_proved": True, "source_functions_not_cap_values": True,
        "source_seam_proofs": {name: copy.deepcopy(dict(value)) for name, value in proof.items()},
    }
    reduced: list[dict[str, Any]] = []
    removed: list[str] = []
    for term in previous:
        name = term.get("name")
        if name == OLD_FIRST_TERM:
            removed.append(name)
            continue
        if name == OLD_EXTERNAL_TERM:
            removed.append(name)
            reduced.append(copy.deepcopy(replacement))
        else:
            reduced.append(copy.deepcopy(term))
    if sorted(removed) != sorted((OLD_FIRST_TERM, OLD_EXTERNAL_TERM)) or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 186-term view required")
    old_unchanged = [term for term in previous if term.get("name") not in {OLD_FIRST_TERM, OLD_EXTERNAL_TERM}]
    new_unchanged = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(old_unchanged) != UNCHANGED_COUNT or new_unchanged != old_unchanged:
        raise ValueError("Other 185 correlated sectors must remain verbatim")
    if sum(term.get("name") == NEW_EXTERNAL_TERM for term in reduced) != 1:
        raise ValueError("Exactly one first-switch external replacement required")
    target = copy.deepcopy(parent_sample.get("independent_target"))
    if target is None:
        target = copy.deepcopy(parent.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")
    return {
        "Z": Z, "source_family": getattr(owner, "family", parent.get("source_family")),
        "raw_m15_C0_C1_terms": copy.deepcopy(raw), "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT, "previous_reduced_terms": PARENT_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT, "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": _trace_chain(), "source_seam_proofs": {"first_second": copy.deepcopy(dict(proof["first_second"]))},
        "actual_signed_function_ports": [{"chart": "first_switch", "side": "lower", "coordinate_node": FIRST_COORDINATE_NODE,
            "endpoint_node": FIRST_ENDPOINT_NODE, "endpoint_scope": {"coordinate_node": FIRST_COORDINATE_NODE, "endpoint_node": FIRST_ENDPOINT_NODE},
            "modified_Qy_nodes": list(FIRST_QY_NODES), "modified_RQy_nodes": list(FIRST_RQY_NODES),
            "definition": "D=R*(Qy_modified-Qy_leading)", "ordinary_Z_orders": [0, 1], "physical_R_once": True,
            "graph_scope_is_explicit": True}],
        "external_endpoint_budget_witness": {
            "own_first_switch_native0": inlet_witness,
            "parent_second_switch_external_term": copy.deepcopy(old_external),
            "parent_second_switch_external_term_active": False,
            "positive_union": owner.record(bound), "positive_union_inputs": copy.deepcopy(inlet_atoms),
            "source_binding": copy.deepcopy(dict(source)), "full_first_switch_raw50_aggregate_used": False,
            "raw50_and_previous_external_provenance_only": True,
            "no_cap_subtraction_value_midpoint_or_overlap": True,
        },
        "old_removed_upstream_viscous_terms": [copy.deepcopy(raw[FIRST_SWITCH_VISCOUS_INDEX])],
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_external),
        "replacement_outer_endpoints": copy.deepcopy(replacement),
        "packet_prefix_witness": copy.deepcopy(packet["exact_prefix_witness"]),
        "selected_N_evidence": copy.deepcopy(selected), "independent_target": target,
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True, "other_185_verbatim": True, "independent_target_unchanged": True,
            "reserved_retained_once": sum(term.get("name") == "reserved_n1_own_swirl_once" for term in reduced) == 1,
            "ordinary_Z1_and_physical_R_once": True, "leading_subtraction_retained": True,
            "first_second_source_join_proved": True, "first_switch_viscous_sector_removed_once": True,
            "parent_external_term_replaced_once": True, "first_switch_native0_budget_owned": True,
            "source_functions_not_cap_values": True, "old_second_switch_external_removed": True,
            "new_first_switch_external_once": True, "previous_187_to_new_186": True,
            "positive_union_own_first_switch0_only": True, "second_switch_external_not_active_after_seam": True,
            "no_error_cap_as_value": True, "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True, "source_bound_signed_trace_contraction": True,
        },
        "limits": {"first_second_internal_trace_resolved": True, "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False, "numeric_Rp_closure": False, "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False, "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False, "mean_oscillatory_correction": False,
            "independent_corrected_Cartesian_NS": False},
    }


def build_switch_numeric(owner: Any, seam_proofs: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved first-switch to second-switch 186-term contraction."""
    proofs = _proofs(seam_proofs)
    parent = _parent_report(owner)
    family = getattr(owner, "family", parent.get("source_family"))
    if parent.get("source_family") != family:
        raise ValueError("Parent source family changed")
    first = _first_report(owner)
    first_pressure = _first_pressure_report(owner)
    r110 = _r110_report(owner)
    for report_name, report in (("first", first), ("first pressure", first_pressure), ("R110", r110)):
        if report.get("source_family") != family or _mapping(report.get("checks"), f"{report_name} checks").get("passed") is not True:
            raise ValueError(f"Accepted {report_name} report is not passed")
    ast_binding = _pin_ast_sources(owner)
    guard = _first_native_guard(owner, proofs["first_second"])
    graph = _graph_binding(owner, parent, guard, proofs["first_second"])
    source = _source_binding(owner, parent, proofs, graph, guard, ast_binding)
    samples: list[dict[str, Any]] = []
    for Z in SAVED_Z:
        parent_sample = _find(parent.get("samples"), label="latest upstream parent sample", Z=Z)
        samples.append(_evaluate_sample(owner, parent, first, first_pressure, r110, source, proofs, graph, parent_sample))
    selected = [sample["selected_N_evidence"] for sample in samples]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All switch samples must use one selected N")
    source.update({"selected_N_evidence": copy.deepcopy(selected[0]), "first_second_source_proof_sha256": _digest(proofs["first_second"])})
    return {
        "source_binding": source, "samples": samples,
        "source_seam_proofs": {"passed": True, "joins": {"first_second": copy.deepcopy(proofs["first_second"])}},
        "contract": {"raw_213_retained": True, "previous_187_to_new_186": True,
            "unchanged_other_185": True, "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True, "cap_is_not_signed_function_value": True,
            "first_second_source_join_required": True, "first_switch_native0_budget_required": True,
            "source_functions_not_cap_values": True, "old_second_switch_external_removed": True,
            "new_first_switch_external_once": True, "raw50_and_old_external_provenance_only": True},
        "limits": {"numeric_Rp_closure": False, "global_all_five_closure": False, "full_goal_complete": False},
    }


__all__ = ["build_switch_numeric"]
