"""Contract the saved long--Rm/Rh finite-N viscous endpoint packet.

This adapter consumes the already saved long/reference 197-term packet and
the checked Rm m15 packet.  It does not construct an Rm producer, replay an
ancestor, or evaluate a signed endpoint.  The only changed correlated view
is the replacement of the saved long/reference external viscous term and the
removal of the now-used Rm patch viscous aggregate; the other 195 entries are
copied bitwise.
"""

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import PREFIX, ROOT
from .original_micro_viscous_m15 import contract_trace_chain


SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PREVIOUS_REDUCED_COUNT = 197
NEW_REDUCED_COUNT = 196
UNCHANGED_COUNT = 195
RM_VISCOUS_INDEX = 167
INCOMING_COUNT = 156

OLD_RM_TERM = "Rm_patch_finite_N_viscous_error"
OLD_EXTERNAL_TERM = "long_reference_restore_finite_N_viscous_external_endpoints"
NEW_EXTERNAL_TERM = "long_Rm_Rh_finite_N_viscous_external_endpoints"

RM_REPORT = PREFIX + "compliant_current_original_Rm_m15_numeric.json"
RM_PRESSURE_REPORT = PREFIX + "compliant_current_original_Rm_pressure_rows_numeric.json"
RM_ARCHIVE = PREFIX + "compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
RM_SOURCE = PREFIX + "compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
RM_C6_SOURCE = PREFIX + "compliant_current_original_true_Rm_patch_C6_functions.py"
RM_M15_SOURCE = "src/openai_ns_reconstruction/original_Rm_m15.py"

RM_REPORT_SHA256 = "cfcb87fd7af7451dc936409beb7d6515e1eabd52f8f85373d28ff667c8dc3452"
RM_PRESSURE_REPORT_SHA256 = "7d4ba34214321f03521cb04858542e6afa70436db458f93685dad1a0d42111dd"
RM_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
RM_SOURCE_SHA256 = "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15"
RM_C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
RM_M15_SOURCE_SHA256 = "b7f2e5dc8f372fa3d09db3498f4a9bb477777e6b095037d046d5698032e8c5cb"


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


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value, key)
    nested = value.get(key)
    return _mapping(nested, key) if isinstance(nested, Mapping) else value


def _path(root: Any, relative: str) -> Path:
    candidates = []
    if root is not None:
        candidates.append(Path(root) / relative)
    candidates.extend((ROOT / relative, ROOT / "experiments" / "root_st073" / relative))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0] if candidates else ROOT / relative


def _pin_file(owner: Any, relative: str, expected: str) -> dict[str, Any]:
    path = _path(getattr(owner, "directory", None), relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    code_hashes = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes")
    known_code = code_hashes.get(relative)
    if known_code is not None and known_code != actual:
        raise ValueError(f"Changed accepted source: {relative}")
    # Serialized reports and archives are registered in the owner's input
    # hash map under their basename.  Require that registration whenever it
    # exists, so an in-memory or alternate-path object cannot bypass the
    # accepted input receipt.
    input_hashes = _mapping(getattr(owner, "hashes", {}), "owner.hashes")
    basename = Path(relative).name
    known_input = input_hashes.get(basename)
    if known_input is not None and known_input != actual:
        raise ValueError(f"Changed accepted input: {relative}")
    if expected is not None and actual != expected:
        raise ValueError(f"Unexpected accepted source: {relative}")
    return {
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _load_rm_archive(owner: Any) -> Mapping[str, Any]:
    candidates = (
        getattr(owner, "rm_bindings", None),
        getattr(owner, "rm_c5", None),
        getattr(owner, "finite_N_bindings", None),
        getattr(owner, "bindings", None),
    )
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            if isinstance(candidate.get("Rm_patch_finite_N_C5_bindings"), Mapping):
                return candidate
            if isinstance(candidate.get("chart"), Mapping) and isinstance(candidate.get("endpoints"), Mapping):
                return {"Rm_patch_finite_N_C5_bindings": candidate}
    path = _path(getattr(owner, "directory", None), RM_ARCHIVE)
    if not path.is_file():
        raise ValueError("Saved Rm finite-N C5 archive required")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, Mapping) or not isinstance(value.get("Rm_patch_finite_N_C5_bindings"), Mapping):
        raise ValueError("Saved Rm finite-N C5 bindings required")
    return value


def _rm_report(owner: Any) -> Mapping[str, Any]:
    # Always load and hash the accepted serialized report first.  An owner
    # candidate is usable only when it is byte-for-byte equivalent to that
    # pinned report; a convenient in-memory ``samples`` mapping must never
    # silently replace the accepted input.
    path = _path(getattr(owner, "directory", None), RM_REPORT)
    if not path.is_file():
        raise ValueError("Saved original Rm m15 report required")
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != RM_REPORT_SHA256:
        raise ValueError("Changed saved original Rm m15 report")
    value = json.loads(raw)
    if not isinstance(value, Mapping) or not isinstance(value.get("samples"), list):
        raise ValueError("Saved original Rm m15 samples required")

    candidates = (
        getattr(owner, "rm", None),
        getattr(owner, "rm_m15", None),
        getattr(owner, "Rm", None),
        getattr(owner, "rm_report", None),
    )
    for candidate in candidates:
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list):
            if isinstance(candidate.get("checks"), Mapping) and candidate["checks"].get("passed") is not True:
                raise ValueError("Checked saved original Rm m15 report required")
            if _canonical(candidate) != _canonical(value):
                raise ValueError("In-memory Rm m15 report does not match the pinned report")
            return value
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list):
            if isinstance(report.get("checks"), Mapping) and report["checks"].get("passed") is not True:
                raise ValueError("Checked saved original Rm m15 report required")
            if _canonical(report) != _canonical(value):
                raise ValueError("In-memory Rm m15 report does not match the pinned report")
            return value
    return value


def _parent_report(owner: Any) -> Mapping[str, Any]:
    parent = getattr(owner, "parent", None)
    if isinstance(parent, Mapping) and isinstance(parent.get("samples"), list):
        return parent
    raise ValueError("Latest long/reference parent report required")


def _parent_source(parent: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in ("source_binding", "numerical_source_binding"):
        if isinstance(parent.get(key), Mapping):
            return parent[key]
    raise ValueError("Latest long/reference source binding required")


def _rm_finite_errors_ast(owner: Any) -> dict[str, Any]:
    """Pin the old Rm ``finite_errors`` implementation without executing it."""
    pin = _pin_file(owner, RM_M15_SOURCE, RM_M15_SOURCE_SHA256)
    tree = ast.parse(_path(getattr(owner, "directory", None), RM_M15_SOURCE).read_bytes())
    cls = next((node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "OriginalRmM15"), None)
    if cls is None:
        raise ValueError("OriginalRmM15 class missing")
    method = next((node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == "finite_errors"), None)
    if method is None:
        raise ValueError("OriginalRmM15.finite_errors method missing")
    calls = {
        node.func.attr
        for node in ast.walk(method)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    required_calls = {"background", "evaluate", "correction_densities", "prefix_error_packet", "sum_logs"}
    if not required_calls.issubset(calls):
        raise ValueError("Rm finite_errors AST no longer carries complete prefix/endpoint bounds")
    method_dump = ast.dump(method, include_attributes=False)
    return {
        "passed": pin["passed"],
        "source": RM_M15_SOURCE,
        "source_sha256": pin["sha256"],
        "class": "OriginalRmM15",
        "method": "OriginalRmM15.finite_errors",
        "method_ast_sha256": hashlib.sha256(method_dump.encode("utf-8")).hexdigest(),
        "required_calls": sorted(required_calls),
        "complete_prefix_error_source": True,
        "endpoint_and_viscous_power_rows": [2, 1],
        "source_not_executed": True,
    }


def _rm_graph_binding(owner: Any, archive: Mapping[str, Any]) -> dict[str, Any]:
    rm = _unwrap(archive, "Rm_patch_finite_N_C5_bindings")
    chart = _mapping(rm.get("chart"), "Rm C5 chart")
    # The signed external port is evaluated at Rh, but its Q_y/R_Q_y
    # functions must remain the original chart functions.  The endpoint
    # entries under ``endpoints`` are substitutions at x=1 and are reserved
    # for the inlet seam proof below; they are not the external Rh function.
    physical = _mapping(rm.get("physical"), "Rm C5 chart physical binding")
    qy = physical.get("Q_y")
    if not isinstance(qy, list) or len(qy) < 2:
        raise ValueError("Rm chart Q_y C0/Z1 graph rows required")
    boundary = _mapping(rm.get("boundary"), "Rm C5 boundary")
    rqy = boundary.get("R_Q_y")
    if not isinstance(rqy, list) or len(rqy) < 2:
        raise ValueError("Rm chart R_Q_y C0/Z1 graph rows required")
    nodes = getattr(owner, "nodes", None)
    if not isinstance(nodes, (list, tuple, Mapping)):
        raise ValueError("Full Rm C5 graph nodes required")
    radius = chart.get("radius")
    if not isinstance(radius, int):
        raise ValueError("Rm physical radius node required")
    domain = chart.get("domain")
    if domain != [1, 2948]:
        raise ValueError("Rm chart native domain must expose x=1 and x=exp(1) node 2948")
    endpoint_scopes = {
        "inlet_x1": {
            "chart": "actual_patch",
            "coordinate_node": chart.get("coordinate"),
            "endpoint_node": 1,
            "native_argument": {"kind": "exact_decimal_Rm_patch", "value": "1"},
            "native_x": "1",
            "endpoint": "Rm",
            "source_seam_role": "postrestore native1 to Rm patch native x=1",
        },
        "outlet_Rh": {
            "chart": "actual_patch",
            "coordinate_node": chart.get("coordinate"),
            "endpoint_node": 2948,
            "native_argument": {"kind": "exact_exp1", "value": "exp(1)"},
            "native_x": "Rh",
            "endpoint": "Rh",
            "external_port_role": "upper signed endpoint",
        },
    }

    graph_rows = []
    for order, node in enumerate(rqy[:2]):
        try:
            row = nodes[node] if not isinstance(nodes, Mapping) else nodes[node]
        except (IndexError, KeyError, TypeError):
            raise ValueError("Rm R_Q_y graph node missing") from None
        if not isinstance(row, Mapping) or row.get("operation") != "product":
            raise ValueError("Rm R_Q_y must be an actual product")
        arguments = row.get("arguments", [])
        if not isinstance(arguments, list) or arguments.count(radius) != 1:
            raise ValueError("Rm R_Q_y must contain one physical R factor")
        try:
            qy_row = nodes[qy[order]] if not isinstance(nodes, Mapping) else nodes[qy[order]]
        except (IndexError, KeyError, TypeError):
            raise ValueError("Rm chart Q_y graph node missing") from None
        if not isinstance(qy_row, Mapping):
            raise ValueError("Rm chart Q_y must be an actual graph function")
        graph_rows.append(
            {
                "ordinary_Z_order": order,
                "Qy_function_node": qy[order],
                "RQy_function_node": node,
                "physical_radius_node": radius,
                "physical_R_explicit_once": True,
                "endpoint_scopes": copy.deepcopy(endpoint_scopes),
            }
        )
    return {
        "passed": True,
        "chart": {
            key: copy.deepcopy(chart.get(key))
            for key in ("coordinate", "domain", "radius", "physical_measure", "global_offset", "actual_phase")
        },
        "physical_measure": "dy=dx/x once",
        "radius_formula": "Rm*x",
        "native_domain": ["1", "e"],
        "boundary_Q_plus_Q_y": copy.deepcopy(boundary.get("Q_plus_Q_y")),
        "boundary_R_Q_y": copy.deepcopy(rqy),
        "initial_Q_plus_Q_y": copy.deepcopy(_mapping(rm.get("initial"), "Rm C5 initial").get("Q_plus_Q_y")),
        "initial_R_Q_y": copy.deepcopy(_mapping(rm.get("initial"), "Rm C5 initial").get("R_Q_y")),
        "chart_physical": copy.deepcopy(physical),
        "graph_function_rows": graph_rows,
        "endpoint_scopes": endpoint_scopes,
        "inlet_source_seam_scope": copy.deepcopy(endpoint_scopes["inlet_x1"]),
        "external_outlet_scope": copy.deepcopy(endpoint_scopes["outlet_Rh"]),
        "original_phase": _mapping(rm.get("bindings"), "Rm C5 bindings").get("original_all_N_source", {}).get(
            "original_global_phase"
        ),
        "whole_source_archive_scope_only": True,
    }


def _rm_source_binding(owner: Any, parent: Mapping[str, Any], seam: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    if seam.get("passed") is not True:
        raise ValueError("Postrestore native1/Rm patch x=1 source seam must pass")
    archive = _load_rm_archive(owner)
    graph = _rm_graph_binding(owner, archive)
    ast_binding = _rm_finite_errors_ast(owner)
    pins = {
        "rm_m15_report": _pin_file(owner, RM_REPORT, RM_REPORT_SHA256),
        "rm_pressure_rows_report": _pin_file(owner, RM_PRESSURE_REPORT, RM_PRESSURE_REPORT_SHA256),
        "rm_finite_N_archive": _pin_file(owner, RM_ARCHIVE, RM_ARCHIVE_SHA256),
        "rm_finite_N_source": _pin_file(owner, RM_SOURCE, RM_SOURCE_SHA256),
        "rm_patch_C6_source": _pin_file(owner, RM_C6_SOURCE, RM_C6_SOURCE_SHA256),
        "rm_m15_finite_errors": ast_binding,
    }
    parent_source = _parent_source(parent)
    inherited_keys = (
        "producer_rows",
        "function",
        "physical_radius_power",
        "ordinary_Z_orders",
        "positive_sum_domination",
        "ordinary_Z1_product_rule_retained",
        "physical_Pstar_already_in_Qy_error",
        "slow_B_y_over_N_and_fast_B_phi_both_retained",
        "postrestore_entire_original_155_cap_retained",
        "inherited_pressure_bound_validity_not_independently_recertified",
    )
    source = {
        key: copy.deepcopy(parent_source[key])
        for key in inherited_keys
        if key in parent_source
    }
    source.update(
        {
            "passed": True,
            "function": "D=R*(Qy_modified-Qy_leading); long native0 minus Rm/Rh upper source",
            "physical_radius_power": 1,
            "ordinary_Z_orders": [0, 1],
            "ordinary_Z1_product_rule_retained": True,
            "physical_Pstar_already_in_Qy_error": True,
            "rm_patch_physical_measure_once": True,
            "rm_patch_radius_formula": "Rm*x",
            "actual_Rm_patch_finite_N_source_binding": graph,
            "actual_Rm_finite_errors_AST_binding": ast_binding,
            "source_pins": pins,
            "inherited_parent_source_binding_sha256": _digest(parent_source),
            "inherited_parent_source_binding_provenance": copy.deepcopy(dict(parent_source)),
            "inherited_parent_report_sha256": _digest(parent),
            "postrestore_native1_Rm_patch_x1_seam_required": True,
            "actual_source_chain": "long native0 -> reference0 -> restoration0 -> postrestore1=Rm -> Rh",
            "new_external_endpoint_function": "D(long native0)-D(RmRh)",
            "delta_definition": "R*(Qy_modified-Qy_leading)",
            "full_original_Rm_viscous_aggregate_used_for_Rh": True,
            "source_functions_values_or_caps_not_selected": True,
            "no_error_cap_as_function_value": True,
            "no_ancestor_constructors_or_producers_run": True,
            "generic_history_bounds_inherited_not_recertified": True,
        }
    )
    rows = list(source.get("producer_rows", []))
    rows.append(
        {
            "source": RM_M15_SOURCE,
            "source_sha256": ast_binding["source_sha256"],
            "method": ast_binding["method"],
            "method_ast_sha256": ast_binding["method_ast_sha256"],
            "complete_prefix_error_source": True,
        }
    )
    source["producer_rows"] = rows
    seams = {
        "passed": True,
        "postrestore_Rm_source_seam": copy.deepcopy(dict(seam)),
        "inlet_source_seam_scope": copy.deepcopy(graph["endpoint_scopes"]["inlet_x1"]),
        "external_outlet_scope": copy.deepcopy(graph["endpoint_scopes"]["outlet_Rh"]),
        "function": "D(long native0)-D(RmRh)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": True,
        "leading_subtraction_retained": True,
        "source_functions_not_error_caps": True,
        "native_Jacobians_not_equated": True,
        "no_phase_reset": True,
    }
    return source, seams


def _rm_report_sample(report: Mapping[str, Any], Z: str, native_x: str) -> Mapping[str, Any]:
    return _find(report.get("samples"), label="Rm m15 sample", Z=Z, chart="actual_patch", native_x=native_x)


def _term_is_zero(term: Mapping[str, Any]) -> bool:
    coefficients = term.get("coefficient_C0_C1")
    if not isinstance(coefficients, list) or len(coefficients) != 2:
        return False
    for row in coefficients:
        if row.get("approximation") not in ("0.0", "0") or row.get("absolute_error_upper") not in ("0.0", "0"):
            return False
        if row.get("raw_directed_error_exact_mpf_tuple") != [0, 0, 0, 0]:
            return False
        if row.get("source_enclosure_exact_mpf_tuples") != [[0, 0, 0, 0], [0, 0, 0, 0]]:
            return False
    return True


def _validate_rm_packet(owner: Any, parent_sample: Mapping[str, Any], rm_report: Mapping[str, Any], Z: str) -> dict[str, Any]:
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Raw 213 m15 state must remain verbatim")
    expected_family = getattr(owner, "family", None)
    if expected_family is not None and rm_report.get("source_family") != expected_family:
        raise ValueError("Saved Rm m15 source family changed")
    if raw[RM_VISCOUS_INDEX].get("name") != OLD_RM_TERM:
        raise ValueError("Raw index 167 must be Rm patch finite-N viscous error")
    rm_x1 = _rm_report_sample(rm_report, Z, "1")
    rm_rh = _rm_report_sample(rm_report, Z, "Rh")
    x1_terms = rm_x1.get("m15_C0_C1_terms")
    incoming = rm_x1.get("incoming_m15_terms")
    if not isinstance(x1_terms, list) or len(x1_terms) != 168:
        raise ValueError("Saved Rm x=1 168-term packet required")
    if not isinstance(incoming, list) or len(incoming) != INCOMING_COUNT or x1_terms[:INCOMING_COUNT] != incoming:
        raise ValueError("Rm x=1 must transfer the raw 156-term incoming prefix exactly")
    if incoming != raw[:INCOMING_COUNT]:
        raise ValueError("Rm x=1 incoming prefix must equal the long/reference raw 156-term prefix")
    if rm_x1.get("exact_empty_increment") is not True:
        raise ValueError("Rm x=1 must carry an exact empty increment")
    tail = x1_terms[INCOMING_COUNT:]
    expected_tail = [
        "Rm_patch_constant_background_integral",
        "Rm_patch_R2F_upper_t1",
        "Rm_patch_R2F_lower_t1",
        "Rm_patch_G_integral_t1",
        "Rm_patch_R2F_upper_t2",
        "Rm_patch_R2F_lower_t2",
        "Rm_patch_G_integral_t2",
        "Rm_patch_RQy_upper_t1",
        "Rm_patch_RQy_lower_t1",
        "Rm_patch_finite_N_endpoint_error",
        "Rm_patch_finite_N_volume_error",
        OLD_RM_TERM,
    ]
    if [term.get("name") for term in tail] != expected_tail or not all(_term_is_zero(term) for term in tail):
        raise ValueError("Rm x=1 twelve empty tails must be exact zero packets")
    rh_terms = rm_rh.get("m15_C0_C1_terms")
    if (not isinstance(rh_terms, list) or len(rh_terms) != 168 or
            rh_terms[:INCOMING_COUNT] != incoming or
            rh_terms[:INCOMING_COUNT] != raw[:INCOMING_COUNT] or
            rh_terms[RM_VISCOUS_INDEX] != raw[RM_VISCOUS_INDEX]):
        raise ValueError("Saved Rm Rh packet must retain the raw 156-term prefix and viscous aggregate")
    rm_term = raw[RM_VISCOUS_INDEX]
    binding = _mapping(rm_term.get("actual_function_binding"), "Rm viscous actual function binding")
    required = {
        "chart": "actual_patch",
        "incoming156_terms_retained_once": True,
        "scale_degree_is_original_fixed_axial_t_not_epsilon": True,
        "error_enclosure_only": True,
        "both_nonflat_endpoints": True,
        "volume_uses_complete_prefix_not_final_decayed_error": True,
    }
    if any(binding.get(key) != value for key, value in required.items()):
        raise ValueError("Rm finite-N viscous producer binding changed")
    logs = rm_term.get("log_scale_terms")
    if not isinstance(logs, list) or len(logs) != 1:
        raise ValueError("Entire Rm viscous aggregate must have one positive log bound")
    atom = logs[0]
    if atom.get("positive_error_bound_atom_only") is not True or atom.get("source_value_selected") is not False:
        raise ValueError("Rm viscous aggregate must remain an error-only positive atom")
    if rm_x1.get("actual_selected_integer") != rm_rh.get("actual_selected_integer"):
        raise ValueError("Rm x=1 and Rh must use the same selected N")
    parent_selected = _mapping(parent_sample.get("selected_N_evidence"), "parent selected N evidence")
    parent_selected_postrestore = parent_selected.get("postrestore_native1")
    if not isinstance(parent_selected_postrestore, Mapping):
        raise ValueError("Parent selected N evidence must expose postrestore native1")
    if any(parent_selected.get(key) != parent_selected_postrestore for key in
           ("long_reshape_native1", "reference_native0", "postrestore_native1")):
        raise ValueError("Parent selected N evidence is not common across the joined charts")
    if rm_x1.get("actual_selected_integer") != parent_selected_postrestore:
        raise ValueError("Rm selected N must equal the parent postrestore selected N evidence")
    return {
        "x1_sample_sha256": _digest(rm_x1),
        "Rh_sample_sha256": _digest(rm_rh),
        "x1_native_coordinate": "1",
        "Rh_native_coordinate": "Rh",
        "incoming156_transfer": True,
        "empty_tail_count": len(tail),
        "empty_tail_names": expected_tail,
        "aggregate_raw_index": RM_VISCOUS_INDEX,
        "aggregate_name": OLD_RM_TERM,
        "aggregate_source_binding": copy.deepcopy(dict(binding)),
        "aggregate_log_atom": copy.deepcopy(dict(atom)),
        "aggregate_matches_saved_Rm_Rh": True,
        "aggregate_dominates_Rh": True,
        "error_enclosure_only": True,
        "source_value_selected": False,
        "selected_integer": copy.deepcopy(rm_x1.get("actual_selected_integer")),
        "parent_selected_N_evidence": copy.deepcopy(parent_selected_postrestore),
        "source_report": RM_REPORT,
    }


def _parent_long_witness(parent_sample: Mapping[str, Any]) -> Mapping[str, Any]:
    budget = _mapping(parent_sample.get("external_endpoint_budget_witness"), "parent external endpoint budget")
    witness = _mapping(budget.get("own_longreshape0"), "parent own long native0 witness")
    rows = witness.get("radius_weighted_C0_C1_log_bounds")
    if not isinstance(rows, list) or len(rows) != 2:
        raise ValueError("Parent own long native0 C0/Z1 budget rows required")
    return witness


def _trace_chain() -> dict[str, Any]:
    pieces = [
        {"name": "long_reshape", "lower": "long_native0", "upper": "long_native1"},
        {"name": "reference", "lower": "reference_native0", "upper": "reference_native1"},
        {"name": "restoration", "lower": "restoration_native0", "upper": "restoration_native1"},
        {"name": "postrestore", "lower": "postrestore_native0", "upper": "Rm"},
        {"name": "Rm_patch", "lower": "Rm", "upper": "Rh"},
    ]
    aliases = {
        "long_native1": "reference_native0",
        "reference_native1": "restoration_native0",
        "restoration_native1": "postrestore_native0",
        "postrestore_native1": "Rm",
    }
    result = contract_trace_chain(pieces, aliases)
    result.update(
        {
            "external_endpoints": {"lower": "long native0", "upper": "Rm/Rh"},
            "cancelled_interfaces_explicit": [
                "long native1 = reference native0",
                "reference native1 = restoration native0",
                "restoration native1 = postrestore native0",
                "postrestore native1 = Rm patch x=1",
            ],
            "function": "D(long native0)-D(RmRh)",
        }
    )
    return result


def _rm_ports(source: Mapping[str, Any]) -> list[dict[str, Any]]:
    graph = _mapping(source.get("actual_Rm_patch_finite_N_source_binding"), "Rm source binding")
    chart = _mapping(graph.get("chart"), "Rm chart source binding")
    scopes = _mapping(graph.get("endpoint_scopes"), "Rm endpoint scopes")
    outlet = _mapping(scopes.get("outlet_Rh"), "Rm Rh outlet scope")
    rows = graph.get("graph_function_rows")
    if not isinstance(rows, list) or len(rows) != 2:
        raise ValueError("Rm chart graph function rows required for external port")
    return [
        {
            "chart": "actual_patch",
            "side": "upper",
            "coordinate_node": chart.get("coordinate"),
            "endpoint_node": outlet.get("endpoint_node"),
            "native_argument": copy.deepcopy(outlet.get("native_argument")),
            "endpoint_scope": copy.deepcopy(outlet),
            "modified_Qy_nodes": [row.get("Qy_function_node") for row in rows],
            "modified_RQy_nodes": [row.get("RQy_function_node") for row in rows],
            "physical_radius_node": chart.get("radius"),
            "physical_measure_node": chart.get("physical_measure"),
            "definition": "D=R*(Qy_modified-Qy_leading)",
            "ordinary_Z_orders": [0, 1],
            "physical_R_explicit_once_in_both_Z_rows": True,
            "source_upper": "Rh=exp(1) complete original patch aggregate",
        }
    ]


def _evaluate_sample(
    owner: Any,
    source: Mapping[str, Any],
    seams: Mapping[str, Any],
    rm_report: Mapping[str, Any],
    Z: str,
) -> dict[str, Any]:
    parent = _parent_report(owner)
    parent_sample = _find(parent.get("samples"), label="latest long/reference parent sample", Z=Z)
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PREVIOUS_REDUCED_COUNT:
        raise ValueError("Latest long/reference 197-term view required")
    rm_packet = _validate_rm_packet(owner, parent_sample, rm_report, Z)
    witness = _parent_long_witness(parent_sample)
    long_inputs = [value for value in witness["radius_weighted_C0_C1_log_bounds"] if value is not None]
    if len(long_inputs) != 2:
        raise ValueError("Both parent long native0 C0/Z1 positive bounds required")
    rm_bound = owner.decode(rm_packet["aggregate_log_atom"])
    bound = owner.positive_union([owner.decode(value) for value in long_inputs] + [rm_bound])

    old_rm = [term for term in previous if term.get("name") == OLD_RM_TERM]
    old_external = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    if len(old_rm) != 1 or len(old_external) != 1:
        raise ValueError("Latest 197-term view must contain one Rm viscous term and one long/reference external term")
    replacement = copy.deepcopy(old_external[0])
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True,
        "function": "D(long native0)-D(RmRh)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_long_native0_then_entire_Rm_viscous_167": True,
        "parent_long_native0_budget_retained": True,
        "entire_original_Rm_viscous_aggregate_used": True,
        "aggregate_dominates_Rh": True,
        "incoming156_transfer_checked": True,
        "x1_empty_twelve_tails_checked": True,
        "postrestore_Rm_source_seam_proved": True,
        "source_functions_values_or_caps_not_selected": True,
        "new_Rm_Rh_source_seam": copy.deepcopy(seams["postrestore_Rm_source_seam"]),
    }

    reduced: list[dict[str, Any]] = []
    removed: list[str] = []
    for term in previous:
        name = term.get("name")
        if name == OLD_RM_TERM:
            removed.append(name)
            continue
        if name == OLD_EXTERNAL_TERM:
            removed.append(name)
            reduced.append(replacement)
        else:
            reduced.append(copy.deepcopy(term))
    if sorted(removed) != sorted((OLD_RM_TERM, OLD_EXTERNAL_TERM)) or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 196-term view required")
    unchanged_old = [
        term for term in previous if term.get("name") not in {OLD_RM_TERM, OLD_EXTERNAL_TERM}
    ]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(unchanged_old) != UNCHANGED_COUNT or unchanged_new != unchanged_old:
        raise ValueError("Other 195 correlated sectors must remain verbatim")

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
        "source_seam_proofs": copy.deepcopy(dict(seams)),
        "actual_signed_function_ports": _rm_ports(source),
        "external_endpoint_budget_witness": {
            "own_longreshape0": copy.deepcopy(dict(witness)),
            "original_Rm_viscous_aggregate_167": copy.deepcopy(rm_packet),
            "positive_union": owner.record(bound),
            "positive_union_inputs": copy.deepcopy(long_inputs) + [copy.deepcopy(rm_packet["aggregate_log_atom"])],
            "source_binding": copy.deepcopy(dict(source)),
            "postrestore_upper_separate_budget_claimed": False,
            "entire_Rm_viscous_aggregate_dominates_Rh": True,
            "no_cap_subtraction_value_midpoint_or_overlap": True,
        },
        "old_budget_terms_retained_as_provenance": [copy.deepcopy(raw[RM_VISCOUS_INDEX]), copy.deepcopy(old_rm[0])],
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_external[0]),
        "replacement_outer_endpoints": copy.deepcopy(replacement),
        "independent_target": target,
        "selected_N_evidence": copy.deepcopy(rm_packet["selected_integer"]),
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True,
            "other_195_verbatim": True,
            "independent_target_unchanged": True,
            "reserved_added_this_stage": 0,
            "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True,
            "long_native0_and_RmRh_external_source_functions_retained": True,
            "parent_long_native0_uses_own_saved_budget": True,
            "entire_original_Rm_viscous_aggregate_167_used": True,
            "aggregate_dominates_Rh": True,
            "raw156_incoming_transfer_exact": True,
            "actual_x1_empty_twelve_tails": True,
            "no_error_cap_as_value": True,
            "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True,
            "source_bound_signed_trace_contraction": True,
        },
        "limits": {
            "long_Rm_Rh_internal_trace_resolved": True,
            "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False,
            "fresh_upper_pressure_packet_proof": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "mean_oscillatory_correction": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


def build_Rm_inlet_numeric(owner: Any, seam_proof: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved long--Rm/Rh 196-term contraction packet."""
    if not isinstance(seam_proof, Mapping):
        raise ValueError("A postrestore native1/Rm patch x=1 seam mapping is required")
    parent = _parent_report(owner)
    source, seams = _rm_source_binding(owner, parent, seam_proof)
    report = _rm_report(owner)
    samples = [_evaluate_sample(owner, source, seams, report, Z) for Z in SAVED_Z]
    selected = [sample["selected_N_evidence"] for sample in samples]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All saved Rm samples must use one selected N")
    return {
        "source_binding": source,
        "samples": samples,
        "source_seam_proofs": seams,
        "contract": {
            "raw_213_retained": True,
            "previous_197_to_new_196": True,
            "unchanged_other_195": True,
            "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True,
            "cap_is_not_signed_function_value": True,
            "postrestore_Rm_source_scope_explicit": True,
            "incoming156_and_x1_empty_tails_verified": True,
        },
        "limits": {
            "numeric_Rp_closure": False,
            "global_all_five_closure": False,
            "full_goal_complete": False,
        },
    }


__all__ = ["build_Rm_inlet_numeric"]
