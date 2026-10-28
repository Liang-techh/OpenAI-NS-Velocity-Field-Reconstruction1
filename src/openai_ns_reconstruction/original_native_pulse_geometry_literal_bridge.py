"""Strict, source-only geometry bridge draft for the current physical pulse.

This module is deliberately independent of the existing source helper and
consumer chain.  The low-level :func:`build_geometry_bridge` accepts literal
rows and authenticated metadata, while the one-shot
:func:`load_geometry_bridge_once` reads only the frozen source reports named
below.  The loader hashes every raw input before it decodes either gzip
archive; it never calls an earlier helper, graph consumer, evaluator, or
checker.

The bridge has a deliberately narrow theorem surface:

* accept the saved radius-stage identity 31162 = 263349 only in its recorded
  ``absolute-radius log recipe only`` scope;
* prove the old and current Rp-squared products are congruent by rewriting the
  two literal product rows with that scoped source identity;
* prove the three literal chart-radius rewrites
  ``exp(log(Rp) + y) = Rp * exp(y)`` for the same y port and positive-radius
  source law;
* record common Pstar, finite N, ordinary-Z and Z-independent assumptions;
* preserve the saved physical pressure correction ports without reducing them.

Hashes pin provenance only.  They never serve as an identity proof.  This
draft intentionally does not assert complete pressure m15 equality, endpoint
or measure equality, or Z 0--Z 3 correction cancellation.  Missing correction
ports and missing JSON pointers fail closed; no saved correction constants are
used as defaults.
"""


import gzip
import hashlib
import json
import os


SCHEMA = "current_physical_pulse_geometry_bridge_kernel_v2"

BASE_GRAPH = "5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e"
PHYSICAL_GRAPH = "29d2257340a9bb9399448dfe0ab6bd2cd31e1815d3646741a242dac506a0a84c"

MASTER_RAW = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
MASTER_NAME = "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
OWN_RAW = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.json.gz"
)
PHYSICAL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json"
)
PHYSICAL_SOURCE = "src/openai_ns_reconstruction/original_native_physical_pulse_source.py"
SCALE_SOURCE = "src/openai_ns_reconstruction/original_native_physical_scale_source_identity.py"
SCALE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json"
)
SCALE_EXECUTION = (
    "experiments/root_st073/"
    "original_native_physical_scale_source_identity_execution.json"
)
TERMINAL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_terminal_log_repair.json"
)

# These are the authenticated bytes used to pin the literal inputs.  The
# caller must pass the same maps to build_geometry_bridge.  A hash mismatch is
# a hard failure, but a hash match is never used as a source equality proof.
EXPECTED_RAW_HASHES = {
    MASTER_RAW: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    OWN_RAW: "2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b",
    PHYSICAL_REPORT: "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    SCALE_REPORT: "b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760",
    SCALE_EXECUTION: "3edc1abde1256177688a4d8e7e0a78ee21057607b01d56d690f753288a249504",
    PHYSICAL_SOURCE: "e8b33ad22d2e50f9f0509e210425f8f37156e88a1a651615e3dbd21387000463",
    SCALE_SOURCE: "bceca092aca0062214cd14933fc97c60397f600c3740f18edff92118359f8266",
    TERMINAL_REPORT: "9566d63baaf92856b8210e1ce07824a5b12c97b6af7a641bfa4c5cab36f5f040",
}

EXPECTED_DECODED_HASHES = {
    MASTER_RAW: "56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1",
    OWN_RAW: "8cfee4eb2ad845c98df7cb8e7dafac6a32df473c56f02eef56751902293715a3",
}

RADIUS_LOG_DIFFERENCE = [
    {"monomial": ["native_core_logC"], "coefficient": "10"},
    {"monomial": ["old_selected_logCstar"], "coefficient": "-10"},
]

OLD_RP = 31162
OLD_LOG_RP = 31161
OLD_RP2 = 34349
CURRENT_RP = 263349
CURRENT_LOG_RP = 280182
CURRENT_RP2 = 286794
PSTAR = 33572
FINITE_N_SOURCE_HANDLE = 18
ORDINARY_Z_ORDERS = [0, 1, 2, 3]

GEOMETRY_SUBSTITUTION = {
    "31161": 280182,
    "31162": 263349,
    "33805": 280141,
    "33982": 280155,
    "34272": 280169,
}

# Literal rows copied from the authenticated source graphs.  They are kept as
# expected rows rather than generated expressions so that a caller must supply
# the actual source rows and pass exact structural comparison.
EXPECTED_OLD_ROWS = {
    31157: {"operation": "sum", "arguments": [24956, 25304]},
    31160: {"operation": "product", "arguments": [31158, 31159]},
    31161: {"operation": "sum", "arguments": [31157, 31160]},
    31162: {"operation": "analytic_unary", "name": "exp", "argument": 31161},
    33804: {"operation": "sum", "arguments": [31161, 33579]},
    33805: {"operation": "analytic_unary", "name": "exp", "argument": 33804},
    33981: {"operation": "sum", "arguments": [31161, 33830]},
    33982: {"operation": "analytic_unary", "name": "exp", "argument": 33981},
    34271: {"operation": "sum", "arguments": [31161, 34006]},
    34272: {"operation": "analytic_unary", "name": "exp", "argument": 34271},
    34349: {"operation": "product", "arguments": [31162, 31162]},
}

EXPECTED_CURRENT_ROWS = {
    33572: {"operation": "analytic_unary", "name": "exp", "argument": 5},
    18: {
        "operation": "exact_positive_integer_expression",
        "type": "DoubleDyadicInteger",
        "J": 1358356628656378313,
        "definition": "N=2^(2^J)",
        "exact_positive_integer": True,
        "Z_independent": True,
        "outer_integer_and_inner_power_not_materialized": True,
        "one_common_N_for_all_charts_and_control_denominators": True,
    },
    263344: {"operation": "bound_variable", "name": "true_post_2Rc_quiet_s"},
    263346: {"operation": "sum", "arguments": [263345, 12628]},
    263348: {"operation": "product", "arguments": [262918, 263347]},
    263349: {
        "operation": "function_substitution",
        "expression": 263348,
        "variable": 263344,
        "value": 263346,
        "Z_independent_substitution": True,
    },
    280140: {"operation": "analytic_unary", "name": "exp", "argument": 33579},
    280141: {"operation": "product", "arguments": [263349, 280140]},
    280154: {"operation": "analytic_unary", "name": "exp", "argument": 33830},
    280155: {"operation": "product", "arguments": [263349, 280154]},
    280168: {"operation": "analytic_unary", "name": "exp", "argument": 34006},
    280169: {"operation": "product", "arguments": [263349, 280168]},
}

EXPECTED_PHYSICAL_ROWS = {
    280182: {
        "operation": "analytic_unary",
        "name": "log",
        "argument": 263349,
        "source_positive_certificate": "actual native Rp>0",
        "ordinary_Z1_exact_zero": True,
    },
    286794: {"operation": "product", "arguments": [263349, 263349]},
}

CHARTS = {
    "main": {
        "query": 33578,
        "y": 33579,
        "old_log_radius": 33804,
        "old_radius": 33805,
        "current_exp_y": 280140,
        "current_radius": 280141,
        "coordinate": "xi=mu*log(R/Rp)",
    },
    "gap": {
        "query": 33829,
        "y": 33830,
        "old_log_radius": 33981,
        "old_radius": 33982,
        "current_exp_y": 280154,
        "current_radius": 280155,
        "coordinate": "xi=mu*log(R/Rp)",
    },
    "end": {
        "query": 34005,
        "y": 34006,
        "old_log_radius": 34271,
        "old_radius": 34272,
        "current_exp_y": 280168,
        "current_radius": 280169,
        "coordinate": "s=log(R/Rv)",
    },
}


def _row_pointers(path, handles):
    """Return exact JSON pointers for a source graph's integer handles."""
    return {
        str(handle): {"path": path, "pointer": "/exact_graph_nodes/" + str(handle)}
        for handle in handles
    }


EXPECTED_SOURCE_POINTERS = {
    "physical": {
        "schema": "/component_result/schema",
        "geometry": "/component_result/physical_pulse_bindings/geometry",
        "coordinates": "/component_result/physical_pulse_bindings/coordinates",
        "geometry_substitution": "/component_result/physical_pulse_bindings/geometry_substitution",
        "current_logRp": "/component_result/physical_pulse_bindings/current_logRp",
        "pressure_m15": "/component_result/physical_pulse_bindings/pressure/m15",
        "negative_complete_m15_target": "/component_result/physical_pulse_bindings/negative_complete_m15_target",
        "pulse_geometry_correction_m15": "/component_result/physical_pulse_bindings/pulse_geometry_correction_m15",
        "definition_graph": "/component_result/definition_graph",
        "current_logRp_row": "/component_result/definition_graph/appended_nodes/0",
        "current_Rp2_row": "/component_result/definition_graph/appended_nodes/6612",
        "physical_scale_receipt": "/component_result/physical_scale_receipt",
        "source_prerequisite_audit": "/component_result/source_prerequisite_audit",
        "source_inspection_audit": "/component_result/source_inspection_audit",
        "source_scope": "/source_scope",
    },
    "scale": {
        "schema": "/component_result/schema",
        "radius_alias": "/component_result/source_interface/source_aliases/Rp",
        "pstar_alias": "/component_result/source_interface/source_aliases/Pstar",
        "source_identities": "/component_result/source_interface/source_identities",
        "physical_pulse_source_ports": "/component_result/source_interface/physical_pulse_source_ports",
        "remaining_chart_ports": "/component_result/source_interface/remaining_physical_C1_ports",
        "source_prerequisite_audit": "/component_result/source_prerequisite_audit",
        "source_inspection_audit": "/component_result/source_inspection_audit",
        "accepted_radius_stage": "/component_result/native_original_absolute_Rc_Rp_identity_proved",
    },
    "terminal": {
        "schema": "/component_result/schema",
        "radius_stage": "/component_result/absolute_radius_source_audit",
        "source_interface": "/component_result/source_interface",
        "physical_pulse_source": "/component_result/source_interface/physical_pulse_source",
        "native_base_source": "/component_result/source_interface/native_base_source",
        "comparison_identities": "/component_result/comparison_identities",
        "source_prerequisite_audit": "/component_result/source_prerequisite_audit",
        "source_inspection_audit": "/component_result/source_inspection_audit",
    },
    "master": {
        "N_row": "/exact_graph_nodes/18",
        "coordinates": "/Rp_native_C5_pulse_source_bindings/pulse_coordinates",
        "geometry": "/Rp_native_C5_pulse_source_bindings/geometry",
        "pressure_m15": "/Rp_native_C5_pulse_source_bindings/pressure/m15",
    },
    "own": {
        "N_row": "/exact_graph_nodes/18",
        "coordinates": "/own_main_gap_end_fields",
    },
}

EXPECTED_SOURCE_POINTERS["rows"] = {
    "old_rows": {
        **_row_pointers(MASTER_RAW, (31157, 31160, 31161, 31162, 34349)),
        **_row_pointers(OWN_RAW, (33804, 33805, 33981, 33982, 34271, 34272)),
    },
    "current_rows": _row_pointers(MASTER_RAW, tuple(EXPECTED_CURRENT_ROWS)),
    "physical_rows": {
        "280182": {
            "path": PHYSICAL_REPORT,
            "pointer": "/component_result/definition_graph/appended_nodes/0",
        },
        "286794": {
            "path": PHYSICAL_REPORT,
            "pointer": "/component_result/definition_graph/appended_nodes/6612",
        },
    },
}

EXPECTED_CORRECTION_SOURCE_POINTERS = {
    "current_physical_m15": {
        "path": PHYSICAL_REPORT,
        "pointer": "/component_result/physical_pulse_bindings/pressure/m15",
    },
    "frame_only_m15": {
        "path": TERMINAL_REPORT,
        "pointer": "/component_result/source_interface/native_base_source/old_frame_only",
    },
    "negative_complete_m15_target": {
        "path": PHYSICAL_REPORT,
        "pointer": "/component_result/physical_pulse_bindings/negative_complete_m15_target",
    },
    "geometry_correction_m15": {
        "path": PHYSICAL_REPORT,
        "pointer": "/component_result/physical_pulse_bindings/pulse_geometry_correction_m15",
    },
    "source_identity": {
        "path": PHYSICAL_REPORT,
        "pointer": "/component_result/terminal_comparison/0/source_identity",
    },
}

EXPECTED_POSITIVE_RADIUS_POINTERS = {
    "old_Rp": {
        "path": TERMINAL_REPORT,
        "pointer": "/component_result/absolute_radius_source_audit",
    },
    "current_Rp": {
        "path": PHYSICAL_REPORT,
        "pointer": "/component_result/definition_graph/appended_nodes/0",
    },
}

# The correction is intentionally retained, not compared or cancelled here.
SAVED_CURRENT_M15 = [287417, 287491, 287663, 287956]
SAVED_FRAME_ONLY_M15 = [279426, 279542, 279770, 280139]
SAVED_GEOMETRY_CORRECTION_M15 = [290917, 290919, 290921, 290923]


def _row(rows, handle):
    """Return one exact integer-keyed literal row, or None on malformed input."""
    if not isinstance(rows, dict) or not isinstance(handle, int):
        return None
    value = rows.get(handle)
    return value if isinstance(value, dict) else None


def _require(checks, failures, name, condition, detail):
    """Record one closed-world check and retain a short audit detail."""
    passed = condition is True
    checks[name] = {"passed": passed, "detail": detail}
    if not passed:
        failures.append(name)
    return passed


def _match_required(mapping, expected):
    """Require every expected field without accepting missing fields."""
    if not isinstance(mapping, dict):
        return False
    return all(mapping.get(key) == value for key, value in expected.items())


def _preserve_correction(correction):
    """Return the saved correction record without simplification or deletion."""
    if not isinstance(correction, dict):
        return {
            "status": "unavailable",
            "reason": "correction record was not supplied",
            "zero_proved": False,
            "full_pressure_identity_proved": False,
        }
    return {
        "status": "retained",
        "current_physical_m15": correction.get("current_physical_m15"),
        "frame_only_m15": correction.get("frame_only_m15"),
        "negative_complete_m15_target": correction.get("negative_complete_m15_target"),
        "geometry_correction_m15": correction.get("geometry_correction_m15"),
        "source_identity": correction.get("source_identity"),
        "source_pointers": correction.get("source_pointers"),
        "zero_proved": False,
        "full_pressure_identity_proved": False,
    }


def _expected_radius_stage():
    """Return the exact saved radius-stage facts that must be re-presented."""
    return {
        "source_archive": MASTER_NAME,
        "source_graph_digest": BASE_GRAPH,
        "native_Rc": 251066,
        "old_Rc": 30828,
        "native_Rp": CURRENT_RP,
        "old_Rp": OLD_RP,
        "logRc_difference": list(RADIUS_LOG_DIFFERENCE),
        "logRp_difference": list(RADIUS_LOG_DIFFERENCE),
        "quiet_length_difference": [],
        "substituted_endpoint_handle": 48,
        "substituted_endpoint_exact_value": "2",
        "exact_positive_radius_log_laws_only": True,
        "source_values_or_bounds_evaluated": False,
        "old_math_consumers_called": 0,
    }


def _chart_bridge(chart, rows_old, rows_current, rows_physical):
    """Build one chart proof only after all four literal row shapes match."""
    spec = CHARTS[chart]
    old_log = _row(rows_old, spec["old_log_radius"])
    old_radius = _row(rows_old, spec["old_radius"])
    current_exp = _row(rows_current, spec["current_exp_y"])
    current_radius = _row(rows_current, spec["current_radius"])
    current_log = _row(rows_physical, CURRENT_LOG_RP)
    return {
        "chart": chart,
        "same_query_port": spec["query"],
        "same_y_port": spec["y"],
        "coordinate": spec["coordinate"],
        "old_log_radius": {"node": spec["old_log_radius"], "row": old_log},
        "old_radius": {"node": spec["old_radius"], "row": old_radius},
        "current_exp_y": {"node": spec["current_exp_y"], "row": current_exp},
        "current_radius": {"node": spec["current_radius"], "row": current_radius},
        "current_log_Rp": {"node": CURRENT_LOG_RP, "row": current_log},
        "identity": "exp(log(Rp)+y) = Rp*exp(y)",
        "source_level_only": True,
        "positive_radius_assumption": True,
        "chart_measure_proved": False,
        "chart_endpoint_proved": False,
    }


_LOAD_STATE = {"called": False}
_LOADER_ATTESTATION = object()


def _pointer_value(document, pointer):
    """Read one RFC 6901 JSON pointer without interpreting its value."""
    if not isinstance(pointer, str) or (pointer and not pointer.startswith("/")):
        raise KeyError("invalid JSON pointer")
    value = document
    if not pointer:
        return value
    for escaped in pointer[1:].split("/"):
        segment = escaped.replace("~1", "/").replace("~0", "~")
        if isinstance(value, dict):
            if segment not in value:
                raise KeyError(pointer)
            value = value[segment]
        elif isinstance(value, list):
            if segment == "-" or not segment.isdigit():
                raise KeyError(pointer)
            index = int(segment)
            if index >= len(value):
                raise KeyError(pointer)
            value = value[index]
        else:
            raise KeyError(pointer)
    return value


def _workspace_path(workspace_root, relative_path):
    """Resolve a repository-relative frozen input path."""
    root = os.fspath(workspace_root)
    if not isinstance(root, str) or not isinstance(relative_path, str):
        raise TypeError("workspace_root and relative_path must be path-like strings")
    return os.path.join(root, *relative_path.split("/"))


def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _authenticated_documents(workspace_root):
    """Hash every input first, then decode only the two authenticated gzip files."""
    raw_hashes = {}
    raw_bytes = {}
    for relative_path, expected_hash in EXPECTED_RAW_HASHES.items():
        path = _workspace_path(workspace_root, relative_path)
        with open(path, "rb") as stream:
            data = stream.read()
        actual_hash = _sha256_bytes(data)
        raw_hashes[relative_path] = actual_hash
        if actual_hash != expected_hash:
            raise ValueError("raw SHA-256 mismatch: " + relative_path)
        raw_bytes[relative_path] = data

    decoded_hashes = {}
    documents = {}
    for relative_path in (MASTER_RAW, OWN_RAW):
        decoded = gzip.decompress(raw_bytes[relative_path])
        decoded_hash = _sha256_bytes(decoded)
        decoded_hashes[relative_path] = decoded_hash
        if decoded_hash != EXPECTED_DECODED_HASHES[relative_path]:
            raise ValueError("decoded SHA-256 mismatch: " + relative_path)
        documents[relative_path] = json.loads(decoded.decode("utf-8"))

    for relative_path in (PHYSICAL_REPORT, SCALE_REPORT, TERMINAL_REPORT):
        documents[relative_path] = json.loads(
            raw_bytes[relative_path].decode("utf-8")
        )
    return documents, raw_hashes, decoded_hashes


def _extract_pointer_records(documents):
    """Extract graph rows using the pinned source pointer for every row."""
    row_receipts = EXPECTED_SOURCE_POINTERS["rows"]
    extracted = {"old_rows": {}, "current_rows": {}, "physical_rows": {}}
    for group, receipts in row_receipts.items():
        for handle_text, receipt in receipts.items():
            source = documents[receipt["path"]]
            value = _pointer_value(source, receipt["pointer"])
            if not isinstance(value, dict):
                raise TypeError("graph row is not a literal mapping: " + handle_text)
            extracted[group][int(handle_text)] = value
    return extracted


def _chart_receipt(entry, path, pointer, own_fields=False):
    """Normalize one actual chart port while retaining its source pointer."""
    if not isinstance(entry, dict):
        raise TypeError("chart port is not a mapping")
    if own_fields:
        domain = entry["domain"]
        coordinate = domain["coordinate"]
        radius = entry["radius"]
        old_log_radius = entry["log_radius"]
    else:
        coordinate = entry["coordinate"]
        radius = entry["radius"]
        old_log_radius = None
    return {
        "query": entry["query"],
        "y": entry["y"],
        "coordinate": coordinate,
        "radius": radius,
        "old_log_radius": old_log_radius,
        "original_mu_and_actual_current_Rp_used": entry.get(
            "original_mu_and_actual_current_Rp_used"
        ),
        "source_path": path,
        "source_pointer": pointer,
    }


def _extract_chart_metadata(documents):
    """Read main/gap/end ports from physical, master and own namespaces."""
    physical = documents[PHYSICAL_REPORT]
    master = documents[MASTER_RAW]
    own = documents[OWN_RAW]
    physical_base = EXPECTED_SOURCE_POINTERS["physical"]["coordinates"]
    master_base = EXPECTED_SOURCE_POINTERS["master"]["coordinates"]
    own_base = EXPECTED_SOURCE_POINTERS["own"]["coordinates"]
    physical_coordinates = _pointer_value(physical, physical_base)
    master_coordinates = _pointer_value(master, master_base)
    own_coordinates = _pointer_value(own, own_base)
    output = {"physical": {}, "master": {}, "own": {}}
    for chart in CHARTS:
        output["physical"][chart] = _chart_receipt(
            physical_coordinates[chart],
            PHYSICAL_REPORT,
            physical_base + "/" + chart,
        )
        output["master"][chart] = _chart_receipt(
            master_coordinates[chart],
            MASTER_RAW,
            master_base + "/" + chart,
        )
        output["own"][chart] = _chart_receipt(
            own_coordinates[chart],
            OWN_RAW,
            own_base + "/" + chart,
            own_fields=True,
        )
    return output


def _loader_blocked(reason, extra=None):
    """Return a blocked loader result without fabricating correction ports."""
    result = {
        "schema": SCHEMA,
        "status": "blocked",
        "passed": False,
        "checks": {},
        "failures": [reason],
        "geometry_bridge": None,
        "preserved_correction": _preserve_correction(None),
        "scope_limits": {
            "complete_pressure_m15": False,
            "chart_measures": False,
            "chart_endpoints": False,
            "ordinary_Z_0_to_Z_3_cancellation": False,
        },
    }
    if isinstance(extra, dict):
        result["loader_receipt"] = extra
    return result


def load_geometry_bridge_once(workspace_root):
    """Load frozen source reports exactly once and build the geometry bridge.

    All raw files are read and SHA-256 checked before either gzip stream is
    decoded.  The loader then extracts only literal rows and named source
    ports, retaining a JSON pointer for every extracted row and chart port.
    A second call, or any missing/mismatched source, is blocked.
    """
    if _LOAD_STATE["called"]:
        return _loader_blocked("one-shot source loader was already called")
    _LOAD_STATE["called"] = True
    try:
        documents, raw_hashes, decoded_hashes = _authenticated_documents(
            workspace_root
        )
        rows = _extract_pointer_records(documents)
        chart_metadata = _extract_chart_metadata(documents)

        physical = documents[PHYSICAL_REPORT]
        scale = documents[SCALE_REPORT]
        terminal = documents[TERMINAL_REPORT]
        master = documents[MASTER_RAW]
        own = documents[OWN_RAW]
        physical_component = physical["component_result"]
        scale_component = scale["component_result"]
        physical_coordinates = _pointer_value(
            physical, EXPECTED_SOURCE_POINTERS["physical"]["coordinates"]
        )
        physical_definition_graph = _pointer_value(
            physical, EXPECTED_SOURCE_POINTERS["physical"]["definition_graph"]
        )
        physical_scale_receipt = _pointer_value(
            physical, EXPECTED_SOURCE_POINTERS["physical"]["physical_scale_receipt"]
        )
        physical_prerequisites = _pointer_value(
            physical, EXPECTED_SOURCE_POINTERS["physical"]["source_prerequisite_audit"]
        )
        physical_inspection = _pointer_value(
            physical, EXPECTED_SOURCE_POINTERS["physical"]["source_inspection_audit"]
        )
        scale_radius_alias = _pointer_value(
            scale, EXPECTED_SOURCE_POINTERS["scale"]["radius_alias"]
        )
        scale_pstar_alias = _pointer_value(
            scale, EXPECTED_SOURCE_POINTERS["scale"]["pstar_alias"]
        )
        scale_pulse_ports = _pointer_value(
            scale, EXPECTED_SOURCE_POINTERS["scale"]["physical_pulse_source_ports"]
        )
        terminal_radius = _pointer_value(
            terminal, EXPECTED_SOURCE_POINTERS["terminal"]["radius_stage"]
        )
        terminal_interface = _pointer_value(
            terminal, EXPECTED_SOURCE_POINTERS["terminal"]["source_interface"]
        )
        terminal_physical_ports = _pointer_value(
            terminal, EXPECTED_SOURCE_POINTERS["terminal"]["physical_pulse_source"]
        )
        terminal_native_ports = _pointer_value(
            terminal, EXPECTED_SOURCE_POINTERS["terminal"]["native_base_source"]
        )
        master_pressure = _pointer_value(
            master, EXPECTED_SOURCE_POINTERS["master"]["pressure_m15"]
        )
        master_N_definition = _pointer_value(
            master, EXPECTED_SOURCE_POINTERS["master"]["N_row"]
        )
        own_N_definition = _pointer_value(
            own, EXPECTED_SOURCE_POINTERS["own"]["N_row"]
        )
        if (
            master_N_definition != own_N_definition
            or master_N_definition != rows["current_rows"][FINITE_N_SOURCE_HANDLE]
        ):
            raise ValueError("N source handle differs across master and own reports")

        radius_stage = dict(terminal_radius)
        radius_stage["identity_scope"] = scale_radius_alias["identity_scope"]
        radius_stage["adapter_alias"] = scale_radius_alias["adapter_alias"]

        source_identity = {
            "current_Pstar": physical_scale_receipt["Pstar"],
            "old_Pstar": physical_scale_receipt["old_pulse_Pstar"],
            "underlying_formula_node": scale_pstar_alias["underlying_formula_node"],
            "provider_node": scale_pstar_alias["provider_node"],
            "raw_quantity_path": scale_pstar_alias["raw_quantity_path"],
            "materialized_value_alias": scale_pstar_alias["materialized_value_alias"],
            "extra_Pstar_factor_applied": physical_scale_receipt[
                "extra_Pstar_factor_applied"
            ],
            "radius_stage_accepted": (
                scale_component["native_original_absolute_Rc_Rp_identity_proved"]
                is True
                and terminal_radius.get("identity")
                == "log(Rc_native/Rc_old)=log(Rp_native/Rp_old)=10*(native_core_logC-old_selected_logCstar)"
            ),
            "finite_N": True,
            "finite_N_source_handle": FINITE_N_SOURCE_HANDLE,
            "finite_N_definition": rows["current_rows"][FINITE_N_SOURCE_HANDLE],
            "ordinary_Z_orders": physical_component["ordinary_Z_orders"],
            "source_only": (
                physical.get("source_only") is True
                and scale.get("source_only") is True
                and terminal.get("source_only") is True
            ),
            "geometry_substitution_Z_independent": physical_prerequisites[
                "geometry_binding_Z_independent"
            ],
            "source_values_or_bounds_evaluated": terminal_radius[
                "source_values_or_bounds_evaluated"
            ],
            "numeric_values_computed": (
                physical.get("numeric_values_computed") is True
                or scale.get("numeric_values_computed") is True
                or terminal.get("numeric_values_computed") is True
                or physical_inspection.get("numeric_values_computed") is True
            ),
            "old_math_consumers_called": terminal_radius["old_math_consumers_called"],
            "pulse_identity_consumed": False,
            "global_node_alias": False,
        }

        current_log_row = rows["physical_rows"][CURRENT_LOG_RP]
        positive_radius = {
            "old_Rp": {
                "node": OLD_RP,
                "positive": terminal_radius["exact_positive_radius_log_laws_only"],
                "source_values_or_bounds_evaluated": terminal_radius[
                    "source_values_or_bounds_evaluated"
                ],
                "source_pointer": EXPECTED_POSITIVE_RADIUS_POINTERS["old_Rp"],
            },
            "current_Rp": {
                "node": CURRENT_RP,
                "positive": current_log_row.get("source_positive_certificate")
                == "actual native Rp>0",
                "source_pointer": EXPECTED_POSITIVE_RADIUS_POINTERS["current_Rp"],
            },
        }

        correction = {
            "current_physical_m15": _pointer_value(
                physical, EXPECTED_SOURCE_POINTERS["physical"]["pressure_m15"]
            ),
            "frame_only_m15": terminal_native_ports["old_frame_only"],
            "negative_complete_m15_target": _pointer_value(
                physical,
                EXPECTED_SOURCE_POINTERS["physical"]["negative_complete_m15_target"],
            ),
            "geometry_correction_m15": _pointer_value(
                physical,
                EXPECTED_SOURCE_POINTERS["physical"][
                    "pulse_geometry_correction_m15"
                ],
            ),
            "source_identity": _pointer_value(
                physical, EXPECTED_CORRECTION_SOURCE_POINTERS["source_identity"]["pointer"]
            ),
            "source_pointers": EXPECTED_CORRECTION_SOURCE_POINTERS,
        }

        source_records = {
            "actual_reports_loaded_once": True,
            "loader_attestation": _LOADER_ATTESTATION,
            "raw_hashes": raw_hashes,
            "decoded_hashes": decoded_hashes,
            "graph_digests": {
                "base": physical_definition_graph["base_nodes_digest"],
                "physical": physical_definition_graph["full_nodes_digest"],
            },
            "old_rows": rows["old_rows"],
            "current_rows": rows["current_rows"],
            "physical_rows": rows["physical_rows"],
            "geometry_substitution": _pointer_value(
                physical,
                EXPECTED_SOURCE_POINTERS["physical"]["geometry_substitution"],
            ),
            "source_pointers": EXPECTED_SOURCE_POINTERS,
            "chart_metadata": chart_metadata,
            "positive_radius": positive_radius,
            "global_node_alias": False,
            "namespace_receipt": {
                "ports_are_digest_and_handle_namespaced": terminal_interface[
                    "ports_are_digest_and_handle_namespaced"
                ],
                "source_values_or_midpoints_used": terminal_interface[
                    "source_values_or_midpoints_used"
                ],
            },
            "master_frame_only_m15": master_pressure,
            "terminal_current_physical_m15": terminal_physical_ports["pressure_m15"],
            "scale_physical_pulse_ports": scale_pulse_ports,
            "report_schemas": {
                "physical": physical_component["schema"],
                "scale": scale_component["schema"],
                "terminal": terminal["schema"],
            },
            "N_source_receipt": {
                "source_handle": FINITE_N_SOURCE_HANDLE,
                "master": {
                    "path": MASTER_RAW,
                    "pointer": EXPECTED_SOURCE_POINTERS["master"]["N_row"],
                    "definition": master_N_definition,
                },
                "own": {
                    "path": OWN_RAW,
                    "pointer": EXPECTED_SOURCE_POINTERS["own"]["N_row"],
                    "definition": own_N_definition,
                },
            },
            "physical_coordinates": physical_coordinates,
            "source_interface_pointer": EXPECTED_SOURCE_POINTERS["terminal"][
                "source_interface"
            ],
        }

        result = build_geometry_bridge(
            source_records, radius_stage, source_identity, correction
        )
        result["loader_receipt"] = {
            "called_once": True,
            "workspace_root": os.fspath(workspace_root),
            "raw_hashes": raw_hashes,
            "decoded_hashes": decoded_hashes,
            "source_pointers": EXPECTED_SOURCE_POINTERS,
            "correction_source_pointers": EXPECTED_CORRECTION_SOURCE_POINTERS,
            "positive_radius_source_pointers": EXPECTED_POSITIVE_RADIUS_POINTERS,
        }
        return result
    except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError) as error:
        return _loader_blocked(
            "authenticated source extraction failed: " + str(error),
            {"called_once": True},
        )


def build_geometry_bridge(source_records, radius_stage, source_identity, correction):
    """Return a fail-closed geometry-only bridge certificate.

    ``source_records`` must have this shape, with no implicit file loading::

        {
            "raw_hashes": {path: sha256, ...},
            "decoded_hashes": {path: sha256, ...},
            "graph_digests": {"base": BASE_GRAPH, "physical": PHYSICAL_GRAPH},
            "old_rows": {integer_handle: literal_row, ...},
            "current_rows": {integer_handle: literal_row, ...},
            "physical_rows": {integer_handle: literal_row, ...},
            "geometry_substitution": {"31161": 280182, ...},
            "audit_flags": {"numeric_values_computed": False, ...},
        }

    The function does not read or decode any record.  It only compares the
    caller's literal rows with the frozen rows above.  Every missing or
    mismatched premise returns ``status = 'blocked'`` and leaves the original
    correction ports in the returned record.
    """
    preserved = _preserve_correction(correction)
    checks = {}
    failures = []
    records = source_records if isinstance(source_records, dict) else {}
    old_rows = records.get("old_rows")
    current_rows = records.get("current_rows")
    physical_rows = records.get("physical_rows")
    radius = radius_stage if isinstance(radius_stage, dict) else {}

    _require(
        checks,
        failures,
        "actual_frozen_reports_loaded_once",
        records.get("actual_reports_loaded_once") is True,
        "caller-supplied dictionaries alone cannot authenticate this bridge",
    )
    _require(
        checks,
        failures,
        "one_shot_loader_attestation",
        records.get("loader_attestation") is _LOADER_ATTESTATION,
        "only the one-shot loader may provide the authenticated record bundle",
    )
    _require(
        checks,
        failures,
        "every_row_and_port_has_pinned_json_pointer",
        records.get("source_pointers") == EXPECTED_SOURCE_POINTERS,
        EXPECTED_SOURCE_POINTERS,
    )

    _require(
        checks,
        failures,
        "raw_hashes_pin_authenticated_inputs",
        records.get("raw_hashes") == EXPECTED_RAW_HASHES,
        "hashes pin provenance only",
    )
    _require(
        checks,
        failures,
        "decoded_hashes_pin_only_authorized_gzip_records",
        records.get("decoded_hashes") == EXPECTED_DECODED_HASHES,
        "raw SHA must be authenticated before decoding",
    )
    digests = records.get("graph_digests")
    _require(
        checks,
        failures,
        "authenticated_base_graph_digest",
        isinstance(digests, dict) and digests.get("base") == BASE_GRAPH,
        BASE_GRAPH,
    )
    _require(
        checks,
        failures,
        "authenticated_physical_graph_digest",
        isinstance(digests, dict) and digests.get("physical") == PHYSICAL_GRAPH,
        PHYSICAL_GRAPH,
    )

    for handle, expected in EXPECTED_OLD_ROWS.items():
        _require(
            checks,
            failures,
            "old_literal_row_" + str(handle),
            _row(old_rows, handle) == expected,
            expected,
        )
    for handle, expected in EXPECTED_CURRENT_ROWS.items():
        _require(
            checks,
            failures,
            "current_literal_row_" + str(handle),
            _row(current_rows, handle) == expected,
            expected,
        )
    for handle, expected in EXPECTED_PHYSICAL_ROWS.items():
        _require(
            checks,
            failures,
            "physical_literal_row_" + str(handle),
            _row(physical_rows, handle) == expected,
            expected,
        )

    _require(
        checks,
        failures,
        "geometry_substitution_is_exact",
        records.get("geometry_substitution") == GEOMETRY_SUBSTITUTION,
        GEOMETRY_SUBSTITUTION,
    )

    radius_expected = _expected_radius_stage()
    _require(
        checks,
        failures,
        "radius_stage_literal_source_evidence",
        _match_required(radius_stage, radius_expected),
        radius_expected,
    )
    _require(
        checks,
        failures,
        "radius_stage_scoped_alias_only",
        radius.get("identity_scope") == "absolute-radius log recipe only"
        and radius.get("adapter_alias") is False,
        "no global pulse alias",
    )

    identity = source_identity if isinstance(source_identity, dict) else {}
    _require(
        checks,
        failures,
        "common_Pstar_source_port",
        identity.get("current_Pstar") == PSTAR
        and identity.get("old_Pstar") == PSTAR
        and identity.get("underlying_formula_node") == PSTAR
        and identity.get("provider_node") == 25203
        and identity.get("raw_quantity_path") == "outer.owner.owner(ends).Pstar"
        and identity.get("materialized_value_alias") is False
        and identity.get("extra_Pstar_factor_applied") is False,
        PSTAR,
    )
    _require(
        checks,
        failures,
        "radius_stage_acceptance_is_explicit",
        identity.get("radius_stage_accepted") is True,
        "accepted source identity is a premise, not inferred from hashes",
    )
    _require(
        checks,
        failures,
        "finite_N_and_ordinary_Z_scope",
        identity.get("finite_N") is True
        and identity.get("finite_N_source_handle") == FINITE_N_SOURCE_HANDLE
        and identity.get("finite_N_definition") == EXPECTED_CURRENT_ROWS[FINITE_N_SOURCE_HANDLE]
        and isinstance(records.get("N_source_receipt"), dict)
        and records["N_source_receipt"].get("source_handle")
        == FINITE_N_SOURCE_HANDLE
        and records["N_source_receipt"].get("master", {}).get("definition")
        == EXPECTED_CURRENT_ROWS[FINITE_N_SOURCE_HANDLE]
        and records["N_source_receipt"].get("own", {}).get("definition")
        == EXPECTED_CURRENT_ROWS[FINITE_N_SOURCE_HANDLE]
        and identity.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS,
        {
            "finite_N": True,
            "source_handle": FINITE_N_SOURCE_HANDLE,
            "definition": EXPECTED_CURRENT_ROWS[FINITE_N_SOURCE_HANDLE],
            "ordinary_Z_orders": ORDINARY_Z_ORDERS,
        },
    )
    _require(
        checks,
        failures,
        "source_identity_premises_are_nonnumeric",
        identity.get("source_only") is True
        and identity.get("geometry_substitution_Z_independent") is True
        and identity.get("source_values_or_bounds_evaluated") is False
        and identity.get("numeric_values_computed") is False
        and identity.get("old_math_consumers_called") == 0,
        "source-only, Z-independent, no numeric consumer",
    )
    _require(
        checks,
        failures,
        "pulse_identity_not_consumed",
        identity.get("pulse_identity_consumed") is False,
        "no pulse delta identity is a premise of this bridge",
    )
    _require(
        checks,
        failures,
        "global_node_alias_is_false",
        records.get("global_node_alias") is False
        and identity.get("global_node_alias") is False,
        "only the scoped absolute-radius recipe alias is admitted",
    )

    _require(
        checks,
        failures,
        "Rp_alias_31162_to_263349_is_scoped_source_identity",
        radius.get("old_Rp") == OLD_RP
        and radius.get("native_Rp") == CURRENT_RP
        and _row(old_rows, OLD_RP) == EXPECTED_OLD_ROWS[OLD_RP]
        and _row(current_rows, CURRENT_RP) == EXPECTED_CURRENT_ROWS[CURRENT_RP],
        {"old_Rp": OLD_RP, "current_Rp": CURRENT_RP},
    )

    namespace = records.get("namespace_receipt")
    _require(
        checks,
        failures,
        "chart_ports_are_digest_and_handle_namespaced",
        isinstance(namespace, dict)
        and namespace.get("ports_are_digest_and_handle_namespaced") is True
        and namespace.get("source_values_or_midpoints_used") is False,
        "same namespace and handle ports; no source values or midpoint substitution",
    )

    positive = records.get("positive_radius")
    _require(
        checks,
        failures,
        "old_and_current_Rp_have_actual_positive_source_receipts",
        isinstance(positive, dict)
        and isinstance(positive.get("old_Rp"), dict)
        and isinstance(positive.get("current_Rp"), dict)
        and positive["old_Rp"].get("node") == OLD_RP
        and positive["old_Rp"].get("positive") is True
        and positive["old_Rp"].get("source_values_or_bounds_evaluated") is False
        and positive["old_Rp"].get("source_pointer")
        == EXPECTED_POSITIVE_RADIUS_POINTERS["old_Rp"]
        and positive["current_Rp"].get("node") == CURRENT_RP
        and positive["current_Rp"].get("positive") is True
        and positive["current_Rp"].get("source_pointer")
        == EXPECTED_POSITIVE_RADIUS_POINTERS["current_Rp"],
        EXPECTED_POSITIVE_RADIUS_POINTERS,
    )

    chart_metadata = records.get("chart_metadata")
    for chart, spec in CHARTS.items():
        physical_chart = (
            chart_metadata.get("physical", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None
        )
        master_chart = (
            chart_metadata.get("master", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None
        )
        own_chart = (
            chart_metadata.get("own", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None
        )
        metadata_ok = (
            isinstance(physical_chart, dict)
            and isinstance(master_chart, dict)
            and isinstance(own_chart, dict)
            and physical_chart.get("query") == spec["query"]
            and master_chart.get("query") == spec["query"]
            and own_chart.get("query") == spec["query"]
            and physical_chart.get("y") == spec["y"]
            and master_chart.get("y") == spec["y"]
            and own_chart.get("y") == spec["y"]
            and physical_chart.get("coordinate") == spec["coordinate"]
            and master_chart.get("coordinate") == spec["coordinate"]
            and own_chart.get("coordinate") == spec["coordinate"]
            and physical_chart.get("radius") == spec["current_radius"]
            and master_chart.get("radius") == spec["current_radius"]
            and own_chart.get("radius") == spec["old_radius"]
            and own_chart.get("old_log_radius") == spec["old_log_radius"]
            and physical_chart.get("source_pointer")
            == EXPECTED_SOURCE_POINTERS["physical"]["coordinates"]
            + "/"
            + chart
            and master_chart.get("source_pointer")
            == EXPECTED_SOURCE_POINTERS["master"]["coordinates"]
            + "/"
            + chart
            and own_chart.get("source_pointer")
            == EXPECTED_SOURCE_POINTERS["own"]["coordinates"]
            + "/"
            + chart
            and physical_chart.get("original_mu_and_actual_current_Rp_used") is True
            and master_chart.get("original_mu_and_actual_current_Rp_used") is True
        )
        _require(
            checks,
            failures,
            "chart_" + chart + "_actual_same_namespace_query_y_ports",
            metadata_ok,
            "physical/master/own receipts agree on query, y, coordinate and port",
        )

    chart_bridges = {}
    for chart in CHARTS:
        spec = CHARTS[chart]
        old_log = _row(old_rows, spec["old_log_radius"])
        old_radius = _row(old_rows, spec["old_radius"])
        current_exp = _row(current_rows, spec["current_exp_y"])
        current_radius = _row(current_rows, spec["current_radius"])
        current_log = _row(physical_rows, CURRENT_LOG_RP)
        row_ok = (
            old_log == {"operation": "sum", "arguments": [OLD_LOG_RP, spec["y"]]}
            and old_radius == {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": spec["old_log_radius"],
            }
            and current_exp == {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": spec["y"],
            }
            and current_radius == {
                "operation": "product",
                "arguments": [CURRENT_RP, spec["current_exp_y"]],
            }
            and current_log == EXPECTED_PHYSICAL_ROWS[CURRENT_LOG_RP]
            and radius.get("exact_positive_radius_log_laws_only") is True
            and spec["old_log_radius"] != spec["current_radius"]
        )
        _require(
            checks,
            failures,
            "chart_" + chart + "_same_y_positive_exp_log_law",
            row_ok,
            "same y port plus positive-radius source law",
        )
        chart_bridges[chart] = _chart_bridge(
            chart, old_rows, current_rows, physical_rows)
        chart_bridges[chart]["actual_source_receipts"] = {
            "physical": chart_metadata.get("physical", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None,
            "master": chart_metadata.get("master", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None,
            "own": chart_metadata.get("own", {}).get(chart)
            if isinstance(chart_metadata, dict)
            else None,
        }

    correction_ok = (
        isinstance(correction, dict)
        and correction.get("current_physical_m15") == SAVED_CURRENT_M15
        and correction.get("frame_only_m15") == SAVED_FRAME_ONLY_M15
        and correction.get("negative_complete_m15_target")
        == [290912, 290913, 290914, 290915]
        and correction.get("geometry_correction_m15") == SAVED_GEOMETRY_CORRECTION_M15
        and correction.get("source_identity")
        == "Delta_physical=Delta_preRp+(pulse_physical-pulse_frame_only)"
        and correction.get("source_pointers") == EXPECTED_CORRECTION_SOURCE_POINTERS
        and records.get("master_frame_only_m15") == SAVED_FRAME_ONLY_M15
        and records.get("terminal_current_physical_m15") == SAVED_CURRENT_M15
        and isinstance(records.get("scale_physical_pulse_ports"), dict)
        and records["scale_physical_pulse_ports"].get("current_pressure_m15")
        == SAVED_CURRENT_M15
        and records["scale_physical_pulse_ports"].get("old_frame_only")
        == SAVED_FRAME_ONLY_M15
        and records["scale_physical_pulse_ports"].get("negative_complete_m15_target")
        == [290912, 290913, 290914, 290915]
        and records["scale_physical_pulse_ports"].get("pulse_geometry_correction_m15")
        == SAVED_GEOMETRY_CORRECTION_M15
    )
    _require(
        checks,
        failures,
        "saved_pressure_correction_retained_without_cancellation",
        correction_ok,
        SAVED_GEOMETRY_CORRECTION_M15,
    )

    if failures:
        return {
            "schema": SCHEMA,
            "status": "blocked",
            "passed": False,
            "checks": checks,
            "failures": failures,
            "geometry_bridge": None,
            "preserved_correction": preserved,
            "scope_limits": {
                "complete_pressure_m15": False,
                "chart_measures": False,
                "chart_endpoints": False,
                "ordinary_Z_0_to_Z_3_cancellation": False,
            },
        }

    return {
        "schema": SCHEMA,
        "status": "geometry-only",
        "passed": True,
        "checks": checks,
        "failures": [],
        "geometry_bridge": {
            "Rp_alias": {
                "old_node": OLD_RP,
                "current_node": CURRENT_RP,
                "proved": True,
                "scope": "absolute-radius log recipe only",
                "source_level_only": True,
                "global_node_alias": False,
                "source_pointers": {
                    "radius_stage": {
                        "path": TERMINAL_REPORT,
                        "pointer": EXPECTED_SOURCE_POINTERS["terminal"]["radius_stage"],
                    },
                    "alias": {
                        "path": SCALE_REPORT,
                        "pointer": EXPECTED_SOURCE_POINTERS["scale"]["radius_alias"],
                    },
                    "old_row": EXPECTED_SOURCE_POINTERS["rows"]["old_rows"][
                        str(OLD_RP)
                    ],
                    "current_row": EXPECTED_SOURCE_POINTERS["rows"]["current_rows"][
                        str(CURRENT_RP)
                    ],
                },
            },
            "Rp2_multiplicative_congruence": {
                "old_node": OLD_RP2,
                "old_row": EXPECTED_OLD_ROWS[OLD_RP2],
                "current_node": CURRENT_RP2,
                "current_row": EXPECTED_PHYSICAL_ROWS[CURRENT_RP2],
                "rewrite": "product([31162,31162]) -> product([263349,263349])",
                "proved_under_scoped_Rp_alias": True,
                "numeric_evaluation": False,
                "source_pointers": {
                    "old_row": EXPECTED_SOURCE_POINTERS["rows"]["old_rows"][
                        str(OLD_RP2)
                    ],
                    "current_row": EXPECTED_SOURCE_POINTERS["rows"]["physical_rows"][
                        str(CURRENT_RP2)
                    ],
                },
            },
            "chart_radius_exp_log_laws": chart_bridges,
            "Pstar": {
                "underlying_formula_node": PSTAR,
                "current_node": PSTAR,
                "old_node": PSTAR,
                "same_source_port": True,
                "source_pointers": {
                    "physical_geometry": {
                        "path": PHYSICAL_REPORT,
                        "pointer": EXPECTED_SOURCE_POINTERS["physical"]["geometry"],
                    },
                    "physical_scale_receipt": {
                        "path": PHYSICAL_REPORT,
                        "pointer": EXPECTED_SOURCE_POINTERS["physical"][
                            "physical_scale_receipt"
                        ],
                    },
                    "alias": {
                        "path": SCALE_REPORT,
                        "pointer": EXPECTED_SOURCE_POINTERS["scale"]["pstar_alias"],
                    },
                },
            },
            "finite_N": {
                "source_handle": FINITE_N_SOURCE_HANDLE,
                "definition": EXPECTED_CURRENT_ROWS[FINITE_N_SOURCE_HANDLE],
                "finite": True,
                "not_materialized": True,
                "source_pointers": {
                    "master": {
                        "path": MASTER_RAW,
                        "pointer": EXPECTED_SOURCE_POINTERS["master"]["N_row"],
                    },
                    "own": {
                        "path": OWN_RAW,
                        "pointer": EXPECTED_SOURCE_POINTERS["own"]["N_row"],
                    },
                },
            },
            "ordinary_Z": {
                "orders": list(ORDINARY_Z_ORDERS),
                "geometry_substitution_Z_independent": True,
            },
            "global_node_alias": False,
        },
        "preserved_correction": preserved,
        "scope_limits": {
            "complete_pressure_m15": False,
            "chart_measures": False,
            "chart_endpoints": False,
            "ordinary_Z_0_to_Z_3_cancellation": False,
            "source_identity_pulseC3congruence_consumed": False,
        },
    }

