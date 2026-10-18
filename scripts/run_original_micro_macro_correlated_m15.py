"""Save or inspect the merged micro/macro source and numeric packets.

The correlated graph packet consumes accepted source receipts and appends a
relocated R_Q_y seam plus a frozen-macro six-I1 source; that graph packet stays
source-only.  The companion packet computes signed finite-N micro six-I1 C1
enclosures and the reduced micro-volume D C1 enclosure in normalized W units.
Macro signed I1 values, full pressure closure, and global closure remain open.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments/root_st073"
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction import original_micro_macro_correlated_m15 as current
from openai_ns_reconstruction import original_micro_background_I1_numeric as numeric_current


build = current.build
numeric_build = numeric_current.build
PARENT = current.PARENT
PARENT_NAME = Path(PARENT).name
PARENT_SHA256 = current.PARENT_SHA256
CORRELATED = current.CORRELATED
CORRELATED_NAME = Path(CORRELATED).name
CORRELATED_SHA256 = current.CORRELATED_SHA256
Q4_ARCHIVE = current.Q4_ARCHIVE
Q4_SOURCE = current.Q4_SOURCE
Q4_CHECK = current.Q4_CHECK
GRAPH_HELPER = current.GRAPH_HELPER
JOIN_CHECK = current.JOIN_CHECK
JOIN_SOURCE = current.JOIN_SOURCE
JOIN_ARCHIVE = current.JOIN_ARCHIVE
Q4_SHA256 = current.Q4_SHA256
Q4_SOURCE_SHA256 = current.Q4_SOURCE_SHA256
Q4_CHECK_SHA256 = current.Q4_CHECK_SHA256
GRAPH_HELPER_SHA256 = current.GRAPH_HELPER_SHA256
JOIN_PINS = dict(current.JOIN_PINS)

REVIEWED_MODULE_SHA256 = "04c671bc6899dbd2578f4519d843828b2f9ab3a70571fbcc2ba477f1ca05752e"
NUMERIC_REVIEWED_MODULE_SHA256 = "6c25842ef2e6c6f65e07f5180dcefeec226e1e9d1994c5a376b80ef1b6f1af0e"

NUMERIC = "src/openai_ns_reconstruction/original_micro_background_I1_numeric.py"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_macro_correlated_m15.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_macro_correlated_m15.py",
    NUMERIC,
    "scripts/run_original_micro_macro_correlated_m15.py",
)

BASE_NODE_COUNT = current.BASE_NODE_COUNT
CORRELATED_SUFFIX_COUNT = 624
RQY_SUFFIX_COUNT = 102
MACRO_COORDINATE = current.MACRO_COORDINATE
MACRO_RADIUS = current.MACRO_RADIUS
MACRO_LENGTH = current.MACRO_LENGTH
MACRO_U = list(current.MACRO_U)
MACRO_Q = list(current.MACRO_Q)
MACRO_Q_Y = list(current.MACRO_Q_Y)
MACRO_LOWER_RQY = list(current.MACRO_LOWER_RQY)
MACRO_UPPER_RQY = list(current.MACRO_UPPER_RQY)
NAMES = list(current.NAMES)

_REF_FIELDS = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "variable", "lower", "upper"),
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path):
    return json.loads(Path(path).read_bytes())


def _digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def _seq(value, label):
    require(isinstance(value, (list, tuple)), label + " must be a sequence")
    return list(value)


def _pair(value, label):
    value = _seq(value, label)
    require(len(value) == 2 and all(type(item) is int for item in value),
            label + " must be an integer pair")
    return value


def _passed(value, label):
    require(isinstance(value, dict) and value.get("passed") is True, label)


def _mapping(value, label):
    require(isinstance(value, dict), label + " must be an object")
    return value


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256,
            "Accepted R_Q_y seam parent report changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True,
            "Checked R_Q_y seam parent required")
    return report


def pins(report):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, expected in report[key].items():
            require(sha(base / name) == expected, "Changed inherited source pin: " + name)
            count += 1
    return count


def _refs(row):
    operation = row.get("operation")
    result = []
    for field in _REF_FIELDS.get(operation, ()):
        if field not in row:
            raise AssertionError(operation + " row is missing " + field)
        value = row[field]
        if operation == "definite_integral" and field == "variable" and isinstance(value, str):
            continue
        values = value if field == "arguments" else [value]
        require(isinstance(values, (list, tuple)), operation + "." + field + " must be a sequence")
        for offset, child in enumerate(values):
            require(type(child) is int and child >= 0,
                    operation + "." + field + " reference is not a handle")
            result.append((field + ("[" + str(offset) + "]" if field == "arguments" else ""), child))
    return result


def _validate_suffix_rows(rows, base, label):
    rows = _seq(rows, label + " rows")
    for offset, row in enumerate(rows):
        require(isinstance(row, dict), label + " row is malformed")
        operation = row.get("operation")
        require(isinstance(operation, str), label + " row operation is missing")
        for _, child in _refs(row):
            require(child < base + offset,
                    label + " row references a future or foreign suffix node")
    return rows


def _validate_parent_seam(accepted):
    seam = accepted.get("actual_micro_macro_RQy_seam")
    _passed(seam, "Accepted R_Q_y seam packet missing")
    require(seam["source_family"] == accepted["source_family"],
            "R_Q_y seam source family changed")
    patch = seam["graph_patch"]
    require(patch["base_node_count"] == BASE_NODE_COUNT
            and patch["appended_node_count"] == RQY_SUFFIX_COUNT
            and patch["base_prefix_unchanged"] is True,
            "R_Q_y seam graph prefix changed")
    contract = seam["contract"]
    require(contract["cancellation_only_in_adjacent_joined_expression"] is True
            and contract["no_outer_cap_deletion"] is True
            and contract["whole_first_to_macro_pair_not_claimed"] is True,
            "R_Q_y seam cancellation scope changed")
    guards = seam["source_guards"]
    require(guards["passed"] is True
            and guards["accepted_join_common_physical_source"] is True
            and guards["accepted_join_full_Q_y_rows"]
            and guards["accepted_join_physical_R_and_native_J_distinct"] is True,
            "R_Q_y seam source guards changed")
    return seam


def _validate_source_pins(packet, accepted):
    expected = {
        PARENT: PARENT_SHA256,
        CORRELATED: CORRELATED_SHA256,
        Q4_ARCHIVE: Q4_SHA256,
        Q4_SOURCE: Q4_SOURCE_SHA256,
        Q4_CHECK: Q4_CHECK_SHA256,
        GRAPH_HELPER: GRAPH_HELPER_SHA256,
        JOIN_CHECK: JOIN_PINS[Path(JOIN_CHECK).name],
        JOIN_SOURCE: JOIN_PINS[Path(JOIN_SOURCE).name],
        JOIN_ARCHIVE: JOIN_PINS[Path(JOIN_ARCHIVE).name],
    }
    source_pins = packet["source_pins"]
    require(set(source_pins) == set(expected), "Correlated source pin set changed")
    for path, digest in expected.items():
        item = source_pins[path]
        require(set(item) == {"path", "expected_sha256", "passed"}
                and item["path"] == path
                and item["expected_sha256"] == digest
                and item["passed"] is True,
                "Correlated source pin receipt changed: " + path)
        require(sha(ROOT / path) == digest, "Correlated source changed: " + path)
    input_pins = packet["input_pins"]
    require(set(input_pins) == {PARENT, CORRELATED},
            "Correlated input pin set changed")
    for path, digest in ((PARENT, PARENT_SHA256), (CORRELATED, CORRELATED_SHA256)):
        item = input_pins[path]
        require(item["path"] == path and item["expected_sha256"] == digest
                and item["passed"] is True,
                "Correlated input pin changed: " + path)
    # These are inherited source pins, so the new packet must not silently
    # replace a parent input or implementation digest.
    for path, digest in expected.items():
        if path == PARENT:
            continue
        key = Path(path).name if path.startswith("experiments/root_st073/") else path
        target = accepted["input_hashes"] if path.startswith("experiments/root_st073/") else accepted["implementation_hashes"]
        if key in target:
            require(target[key] == digest, "Conflicting inherited pin: " + key)


def _validate_source_guards(packet, correlated, seam):
    guards = packet["source_guards"]
    require(set(guards) == {
        "accepted_correlated_source", "accepted_rqy_join",
        "accepted_rqy_proof_rows", "join_physical_source_binding",
        "operation_specific_relocation", "macro_radius_and_length",
        "second_micro_endpoint_wrapper_identity", "rqy_boundary_contraction",
        "reduced_B_seam",
    }, "Correlated source guard fields changed")

    corr = guards["accepted_correlated_source"]
    require(set(corr) == {"passed", "standalone_equals_latest_parent_packet",
                          "packet_digest", "exact_algebra_receipt", "physical_measure",
                          "accepted_producer_executed"}
            and corr["passed"] is True
            and corr["standalone_equals_latest_parent_packet"] is True
            and corr["packet_digest"] == _digest(correlated)
            and corr["exact_algebra_receipt"] == correlated["exact_algebra"]
            and corr["physical_measure"] == "R^2*hb*ds exactly once"
            and corr["accepted_producer_executed"] is False,
            "Correlated source guard changed")

    require(guards["accepted_rqy_join"] == seam["source_guards"],
            "Accepted R_Q_y seam guard was replaced or truncated")

    proof_rows = _seq(guards["accepted_rqy_proof_rows"], "accepted R_Q_y proof rows")
    require(len(proof_rows) == 4, "Accepted R_Q_y proof row count changed")
    for order, row in enumerate(proof_rows):
        require(isinstance(row, dict)
                and row["ordinary_Z_order"] == order
                and row["conditional_boundary_identity"]["passed"] is True
                and row["conditional_boundary_identity"]["ordinary_Z_order"] == order
                and row["conditional_boundary_identity"]["common_Q_y_source_rows"] is True
                and row["conditional_boundary_identity"]["common_physical_radius_cutpoint"] is True
                and row["conditional_boundary_identity"]["one_common_N_and_global_phase"] is True
                and row["conditional_boundary_identity"]["mean_ODE_and_complete_correction_premises"] is True,
                "Accepted R_Q_y conditional proof row changed")
    binding = guards["join_physical_source_binding"]
    require(isinstance(binding, dict)
            and binding.get("passed") is True
            and binding.get("archive_prefix_equal_to_q4") is True
            and binding.get("common_cutpoints", {}).get("common_N") == 18
            and binding.get("common_cutpoints", {}).get("micro_coordinate") == 66988
            and binding.get("common_cutpoints", {}).get("macro_coordinate") == MACRO_COORDINATE,
            "Joined physical source binding changed")

    wrapper = guards["second_micro_endpoint_wrapper_identity"]
    require(set(wrapper) == {"passed", "rows",
                             "same_expression_variable_value_and_Z_independence"}
            and wrapper["passed"] is True
            and wrapper["same_expression_variable_value_and_Z_independence"] is True,
            "Second-micro endpoint wrapper guard changed")
    wrapper_rows = _seq(wrapper["rows"], "endpoint wrapper rows")
    require(len(wrapper_rows) == 2, "Endpoint wrapper row count changed")
    correlated_upper = _seq(
        correlated["graph_patch"]["endpoint_fields"]["second_micro"]["upper"]["R_Q_y"],
        "correlated second-micro upper R_Q_y",
    )[:2]
    for order, row in enumerate(wrapper_rows):
        require(set(row) == {"ordinary_Z_order", "correlated_handle", "seam_handle",
                              "row", "same_expression_variable_value_and_Z_independence"}
                and row["ordinary_Z_order"] == order
                and row["correlated_handle"] == correlated_upper[order]
                and row["seam_handle"] == [133745, 133746][order]
                and isinstance(row["row"], dict)
                and row["same_expression_variable_value_and_Z_independence"] is True,
                "Second-micro endpoint wrapper row changed")

    contraction = guards["rqy_boundary_contraction"]
    require(set(contraction) == {"passed", "rows", "opposite_FTC_signs"}
            and contraction["passed"] is True
            and contraction["opposite_FTC_signs"] is True,
            "RQy boundary contraction guard changed")
    contraction_rows = _seq(contraction["rows"], "RQy contraction rows")
    require(len(contraction_rows) == 4, "RQy contraction row count changed")
    for order, row in enumerate(contraction_rows):
        require(set(row) == {"ordinary_Z_order", "handle",
                              "negative_second_micro_upper", "frozen_macro_lower", "row"}
                and row["ordinary_Z_order"] == order
                and row["handle"] == [184365, 184383, 184401, 184419][order]
                and row["negative_second_micro_upper"] == [133745, 133746, 133747, 133748][order]
                and row["frozen_macro_lower"] == [184328, 184329, 184330, 184331][order]
                and row["row"] == {
                    "operation": "sum",
                    "arguments": [184324 + 40 + 18 * order, 184324 + 4 + order],
                }, "RQy contraction row changed")

    reduced = guards["reduced_B_seam"]
    require(set(reduced) == {"passed", "micro_B_plus_D_rows", "macro_B_plus_D_rows",
                             "shared_physical_U_Q_source", "internal_endpoint_identity",
                             "outer_B_boundary_fields", "outer_endpoint_wrappers",
                             "ordinary_Z_orders"}
            and reduced["passed"] is True
            and len(_seq(reduced["micro_B_plus_D_rows"], "micro B+D rows")) == 2
            and len(_seq(reduced["macro_B_plus_D_rows"], "macro B+D rows")) == 2
            and reduced["shared_physical_U_Q_source"] is True
            and reduced["internal_endpoint_identity"]["passed"] is True
            and reduced["internal_endpoint_identity"]["row_count"] == 6
            and set(reduced["outer_B_boundary_fields"]) == {
                "R_squared_q", "R_squared_q_squared", "R_squared_U_squared",
            }
            and reduced["ordinary_Z_orders"] == [0, 1],
            "Reduced B seam guard changed")

    relocation = guards["operation_specific_relocation"]
    require(set(relocation) == {
        "passed", "reference_fields", "exact_rational_numerator_denominator_are_literals",
        "definite_integral_variable_strings_preserved", "graph_metadata_integers_not_relocated",
        "source_graph_prefix_unchanged",
    } and relocation["passed"] is True
            and relocation["exact_rational_numerator_denominator_are_literals"] is True
            and relocation["definite_integral_variable_strings_preserved"] is True
            and relocation["graph_metadata_integers_not_relocated"] is True
            and relocation["source_graph_prefix_unchanged"] is True,
            "Graph relocation guard changed")

    macro = guards["macro_radius_and_length"]
    require(set(macro) == {
        "passed", "radius_Z_independent", "length_Z_independent",
        "physical_R_squared_length_applied_once",
    } and macro["passed"] is True
            and macro["radius_Z_independent"] is True
            and macro["length_Z_independent"] is True
            and macro["physical_R_squared_length_applied_once"] is True,
            "Macro radius/length guard changed")


def _validate_graph_namespace(packet):
    namespace = packet["graph_namespace"]
    require(set(namespace) == {
        "base_archive", "base_node_count", "accepted_correlated_suffix_count",
        "relocated_rqy_suffix_count", "new_macro_suffix_count",
        "joined_expression_suffix_count", "merged_node_count", "old_graph_unchanged",
        "namespace_digest", "serialized_suffix_count", "segment_starts",
        "serialized_suffix_covers_every_new_node",
    }, "Merged graph namespace fields changed")
    require(namespace["base_archive"] == Q4_ARCHIVE
            and namespace["base_node_count"] == BASE_NODE_COUNT
            and namespace["accepted_correlated_suffix_count"] == CORRELATED_SUFFIX_COUNT
            and namespace["relocated_rqy_suffix_count"] == RQY_SUFFIX_COUNT
            and type(namespace["new_macro_suffix_count"]) is int
            and namespace["new_macro_suffix_count"] > 0
            and type(namespace["joined_expression_suffix_count"]) is int
            and namespace["joined_expression_suffix_count"] > 0
            and namespace["merged_node_count"] == (
                BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT
                + namespace["new_macro_suffix_count"]
                + namespace["joined_expression_suffix_count"])
            and namespace["old_graph_unchanged"] is True
            and isinstance(namespace["namespace_digest"], str)
            and len(namespace["namespace_digest"]) == 64,
            "Merged graph namespace receipt changed")
    expected_starts = {
        "accepted_correlated": BASE_NODE_COUNT,
        "relocated_rqy": BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT,
        "new_macro": BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT,
        "joined_expression": (
            BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT
            + namespace["new_macro_suffix_count"]
        ),
    }
    require(namespace["serialized_suffix_count"] == (
        CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT
        + namespace["new_macro_suffix_count"]
        + namespace["joined_expression_suffix_count"]
    ) and namespace["segment_starts"] == expected_starts
            and namespace["serialized_suffix_covers_every_new_node"] is True,
            "Merged graph segment receipt changed")


def _validate_safe_graphs(packet, ports):
    graphs = packet["safe_graphs"]
    require(set(graphs) == {"base_q4", "accepted_correlated", "relocated_rqy", "new_macro",
                             "joined_expression"},
            "Safe graph namespace fields changed")

    base = graphs["base_q4"]
    require(base == {"archive": Q4_ARCHIVE, "sha256": Q4_SHA256,
                     "node_count": BASE_NODE_COUNT}, "Base Q4 graph receipt changed")

    corr = graphs["accepted_correlated"]
    require(set(corr) == {"base_node_count", "appended_node_count", "appended_nodes", "source_handles"}
            and corr["base_node_count"] == BASE_NODE_COUNT
            and corr["appended_node_count"] == CORRELATED_SUFFIX_COUNT
            and len(_seq(corr["appended_nodes"], "accepted correlated suffix")) == CORRELATED_SUFFIX_COUNT,
            "Accepted correlated graph suffix changed")
    _validate_suffix_rows(corr["appended_nodes"], BASE_NODE_COUNT,
                          "accepted correlated")
    pair = corr["source_handles"]["actual_micro_pair_increment"]
    _pair(pair, "accepted micro pair increment")

    relocated = graphs["relocated_rqy"]
    require(set(relocated) == {
        "old_base_node_count", "new_base_node_count", "appended_node_count",
        "appended_nodes", "relocation_receipt", "target_ports",
    }, "Relocated R_Q_y graph fields changed")
    relocated_base = BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT
    require(relocated["old_base_node_count"] == BASE_NODE_COUNT
            and relocated["new_base_node_count"] == relocated_base
            and relocated["appended_node_count"] == RQY_SUFFIX_COUNT,
            "Relocated R_Q_y graph bases changed")
    rows = _validate_suffix_rows(relocated["appended_nodes"], relocated_base,
                                 "relocated R_Q_y")
    require(len(rows) == RQY_SUFFIX_COUNT
            and len(_seq(relocated["relocation_receipt"], "R_Q_y relocation receipt")) == RQY_SUFFIX_COUNT,
            "Relocated R_Q_y suffix count changed")
    expected_lower = [relocated_base + i for i in range(4, 8)]
    expected_cancel = [relocated_base + i for i in (41, 59, 77, 95)]
    target = relocated["target_ports"]
    require(set(target) == {"second_micro_upper", "frozen_macro_lower", "signed_cancellation"}
            and target["second_micro_upper"] == [133745, 133746, 133747, 133748]
            and target["frozen_macro_lower"] == expected_lower
            and target["signed_cancellation"] == expected_cancel,
            "Relocated R_Q_y target ports changed")
    require(ports == target, "Remapped R_Q_y ports disagree with safe graph")

    receipt = relocated["relocation_receipt"]
    for offset, item in enumerate(receipt):
        require(set(item) == {
            "old_handle", "new_handle", "operation", "reference_fields",
            "references_before", "references_after", "literal_metadata_preserved",
        }, "R_Q_y relocation receipt schema changed")
        require(item["old_handle"] == BASE_NODE_COUNT + offset
                and item["new_handle"] == relocated_base + offset,
                "R_Q_y relocation handle sequence changed")
        require(isinstance(item["operation"], str)
                and isinstance(item["reference_fields"], list),
                "R_Q_y relocation operation receipt changed")
        require(all(value is True for value in item["literal_metadata_preserved"].values()),
                "R_Q_y relocation changed literal metadata")
        before = item["references_before"]
        after = item["references_after"]
        require(set(before) == set(after), "R_Q_y relocation reference fields changed")
        for key, value in before.items():
            if value < BASE_NODE_COUNT:
                expected = value
            else:
                require(value < BASE_NODE_COUNT + RQY_SUFFIX_COUNT,
                        "R_Q_y relocation reference escapes source suffix")
                expected = relocated_base + value - BASE_NODE_COUNT
            require(after[key] == expected, "R_Q_y reference was not relocated exactly")

    macro = graphs["new_macro"]
    require(set(macro) == {
        "base_node_count", "appended_node_count", "appended_nodes", "source_handles",
    }, "New macro graph fields changed")
    macro_base = relocated_base + RQY_SUFFIX_COUNT
    require(macro["base_node_count"] == macro_base
            and type(macro["appended_node_count"]) is int
            and macro["appended_node_count"] > 0,
            "New macro graph base changed")
    macro_rows = _validate_suffix_rows(macro["appended_nodes"], macro_base, "new macro")
    require(len(macro_rows) == macro["appended_node_count"],
            "New macro suffix count changed")
    handles = macro["source_handles"]
    require(set(handles) == {"six_I1", "increment"}
            and set(handles["six_I1"]) == set(NAMES),
            "New macro source handle set changed")
    for name in NAMES:
        _pair(handles["six_I1"][name], "new macro six-I1 " + name)
    _pair(handles["increment"], "new macro increment")

    joined = graphs["joined_expression"]
    require(set(joined) == {
        "base_node_count", "appended_node_count", "appended_nodes", "source_handles",
    }, "Joined-expression graph fields changed")
    joined_base = macro_base + macro["appended_node_count"]
    require(joined["base_node_count"] == joined_base
            and type(joined["appended_node_count"]) is int
            and joined["appended_node_count"] > 0,
            "Joined-expression graph base changed")
    joined_rows = _validate_suffix_rows(joined["appended_nodes"], joined_base,
                                         "joined expression")
    require(len(joined_rows) == joined["appended_node_count"],
            "Joined-expression suffix count changed")
    joined_handles = joined["source_handles"]
    require(set(joined_handles) == {"uncontracted", "no_seam", "shifted_tail", "pressure_increment"},
            "Joined-expression source handles changed")
    _pair(joined_handles["uncontracted"], "uncontracted joined expression")
    _pair(joined_handles["no_seam"], "no-seam joined expression")
    _pair(joined_handles["shifted_tail"], "shifted-tail joined expression")
    _pair(joined_handles["pressure_increment"], "pressure increment joined expression")


def _validate_macro(packet, ports):
    macro = packet["macro_six_I1_source"]
    require(packet["actual_macro_six_I1_source"] == macro,
            "Macro six-I1 source was duplicated inconsistently")
    require(set(macro) == {
        "passed", "ordinary_Z_orders", "six_I1_names", "physical_fields",
        "physical_measure", "coordinate", "domain", "source", "integrals",
        "endpoint_fields", "formula_context", "reduced_boundary_B", "reduced_volume_D", "reduced_S1",
        "macro_increment", "outer_RQy", "graph_suffix_start", "graph_suffix_append_index",
        "formula", "exact_pressure_source_receipt", "appended_nodes",
    }, "Macro six-I1 source fields changed")
    require(macro["passed"] is True
            and macro["ordinary_Z_orders"] == [0, 1]
            and macro["six_I1_names"] == NAMES
            and macro["coordinate"] == MACRO_COORDINATE
            and macro["domain"] == [0, 1],
            "Macro six-I1 source domain changed")

    fields = macro["physical_fields"]
    require(fields == {
        "U": MACRO_U[:3], "Q": MACRO_Q[:3],
        "U_Z": MACRO_U[1:3], "Q_Z": MACRO_Q[1:3],
    }, "Macro physical U/Q source rows changed")
    measure = macro["physical_measure"]
    require(set(measure) == {
        "radius", "radius_power", "length", "weight",
        "R_squared_length_applied_once", "Z_independent",
    } and measure["radius"] == MACRO_RADIUS
            and measure["radius_power"] == 2
            and measure["length"] == MACRO_LENGTH
            and type(measure["weight"]) is int
            and measure["R_squared_length_applied_once"] is True
            and measure["Z_independent"] is True,
            "Macro physical R-squared-length measure changed")

    require(macro["formula_context"] == {
        "z": [20, 1], "delta": [25220, 0], "L": [184439, 184438],
        "L2": [184440, 184443], "two_id": 48,
    }, "Macro formula context changed")

    require(set(macro["source"]) == set(NAMES)
            and set(macro["integrals"]) == set(NAMES),
            "Macro six-I1 source/integral set changed")
    for name in NAMES:
        _pair(macro["source"][name], "macro source " + name)
        _pair(macro["integrals"][name], "macro integral " + name)

    endpoint_fields = macro["endpoint_fields"]
    require(set(endpoint_fields) == {"lower", "upper"},
            "Macro endpoint field sides changed")
    for side in ("lower", "upper"):
        row = endpoint_fields[side]
        require(set(row) == {"R_squared_q", "R_squared_q_squared",
                              "R_squared_U_squared", "R_Q_y"},
                "Macro endpoint fields changed")
        for key in ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared"):
            _pair(row[key], "macro endpoint " + side + " " + key)
        _pair(row["R_Q_y"], "macro endpoint " + side + " R_Q_y")
    require(_seq(endpoint_fields["lower"]["R_Q_y"], "macro lower R_Q_y")
            == ports["frozen_macro_lower"][:2]
            and _seq(endpoint_fields["upper"]["R_Q_y"], "macro upper R_Q_y")
            == MACRO_UPPER_RQY[:2],
            "Macro outer R_Q_y endpoints changed")
    for key in ("reduced_boundary_B", "reduced_volume_D", "reduced_S1", "macro_increment"):
        _pair(macro[key], "macro " + key)

    outer = macro["outer_RQy"]
    require(set(outer) == {"lower", "upper"}
            and outer["lower"] == ports["frozen_macro_lower"]
            and outer["upper"] == MACRO_UPPER_RQY,
            "Macro outer R_Q_y retention changed")

    macro_base = BASE_NODE_COUNT + CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT
    require(macro["graph_suffix_start"] == macro_base
            and macro["graph_suffix_append_index"] == CORRELATED_SUFFIX_COUNT + RQY_SUFFIX_COUNT,
            "Macro suffix placement changed")
    _validate_suffix_rows(macro["appended_nodes"], macro_base, "macro source")

    require(macro["formula"] == {
        "B": "R^2 Q/L + R^2 Q^2/2 - 2 Z^2 R^2 U^2/L^2",
        "D": "-q/L + h Z q_Z/(2L) - q^2/2 + d U q_Z/L + 2 h Z^2 U^2/L^2 + 2 Z d U U_Z/L^2",
        "increment": "(B+D)/2 - R_Q_y(upper) + R_Q_y(lower)",
        "L": "1-delta*Z^2",
        "coefficients": {
            "U_squared": "2*h*Z^2/L^2",
            "U_U_Z": "2*Z*d/L^2",
            "U_q_Z": "d/L",
        },
    }, "Macro exact source formula changed")
    require(macro["exact_pressure_source_receipt"] == {
        "source_formula": "physical reduced_pressure_S with R^2*length measure",
        "no_extra_Z_on_U_q_Z": True,
        "Pstar_not_added_to_I1_weight": True,
    }, "Macro exact pressure receipt changed")


def _validate_merged_namespace(packet):
    merged = packet["merged_namespace"]
    require(set(merged) == {
        "base_archive", "base_node_count", "appended_nodes", "segment_starts",
        "merged_node_count", "roundtrip_verified",
    }, "Merged namespace fields changed")
    namespace = packet["graph_namespace"]
    graphs = packet["safe_graphs"]
    expected_suffix = (
        _seq(graphs["accepted_correlated"]["appended_nodes"], "accepted suffix")
        + _seq(graphs["relocated_rqy"]["appended_nodes"], "relocated suffix")
        + _seq(graphs["new_macro"]["appended_nodes"], "macro suffix")
        + _seq(graphs["joined_expression"]["appended_nodes"], "joined suffix")
    )
    require(merged["base_archive"] == Q4_ARCHIVE
            and merged["base_node_count"] == BASE_NODE_COUNT
            and merged["appended_nodes"] == expected_suffix
            and merged["segment_starts"] == namespace["segment_starts"]
            and merged["merged_node_count"] == namespace["merged_node_count"]
            and merged["roundtrip_verified"] is True,
            "Merged namespace roundtrip changed")


def _validate_join(packet, macro, ports, correlated):
    require(packet["remapped_accepted_seam_target_ports"] == ports,
            "Remapped accepted seam ports changed")
    joined = packet["actual_joined_pressure_function"]
    require(joined == packet["joined_pressure_function"]
            and set(joined) == {
                "ordinary_Z_orders", "micro_whole_pair", "macro_increment",
                "uncontracted_joined", "joined", "no_seam_joined", "formula",
                "contracted_formula", "reduced_boundary_B", "reduced_volume_D",
                "shared_U_Q_B_seam", "uncontracted_equals_contracted",
                "rqy_boundary_contraction",
                "rqy_boundary_contraction_rows", "rqy_signed_terms",
                "contraction_only_in_joined_expression", "latest_second_micro_increment_not_added",
                "first_micro_lower_endpoint_retained", "macro_outer_RQy_retained",
                "joined_shifted_tail", "shifted_tail_formula", "pressure_increment",
                "absolute_M15_formula", "pressure_increment_scope",
                "pressure_increment_source_defined", "pressure_increment_numerical_enclosure",
                "original_first_lower_M15_retained", "shifted_H_anchor_retained",
                "shifted_H_anchor_source", "shifted_H_anchor_digest", "endpoint_wrapper_identity",
            }, "Joined pressure function fields changed")
    require(joined["ordinary_Z_orders"] == [0, 1]
            and len(_pair(joined["micro_whole_pair"], "micro whole pair")) == 2
            and _seq(joined["macro_increment"], "joined macro increment")
            == _seq(macro["macro_increment"], "macro increment")
            and len(_pair(joined["uncontracted_joined"], "uncontracted joined")) == 2
            and len(_pair(joined["joined"], "joined pressure pair")) == 2
            and joined["no_seam_joined"] == joined["joined"]
            and joined["formula"] == "(B_outer+D_micro+D_macro)/2 - R_Q_y(m=1) + R_Q_y(first_micro lower)"
            and joined["contracted_formula"] == "(outer_B + D_micro + D_macro)/2 - R_Q_y(m=1) + R_Q_y(first_micro lower)"
            and joined["rqy_boundary_contraction"] == ports["signed_cancellation"]
            and joined["rqy_signed_terms"] == "-R_Q_y(second_micro,s=2)+R_Q_y(frozen_macro,m=0)"
            and joined["contraction_only_in_joined_expression"] is True
            and joined["latest_second_micro_increment_not_added"] is True
            and joined["first_micro_lower_endpoint_retained"] is True
            and joined["macro_outer_RQy_retained"] is True
            and len(_pair(joined["joined_shifted_tail"], "joined shifted tail")) == 2
            and joined["shifted_tail_formula"] == ".5*(B_outer+D_micro+D_macro)-R_Q_y(m=1)"
            and joined["pressure_increment"] == joined["joined"]
            and joined["absolute_M15_formula"] == "M15(m=1)=accepted_H+joined_shifted_tail"
            and joined["pressure_increment_scope"] == "first_micro_sc/2_to_macro_m=1"
            and joined["pressure_increment_source_defined"] is True
            and joined["pressure_increment_numerical_enclosure"] is False
            and joined["original_first_lower_M15_retained"] is True
            and joined["shifted_H_anchor_retained"] is True
            and joined["shifted_H_anchor_source"] == "accepted_correlated.shifted_inlet_source_digest"
            and joined["shifted_H_anchor_digest"] == correlated["shifted_inlet_source_digest"],
            "Joined pressure expression changed")

    decomposition = joined["uncontracted_equals_contracted"]
    require(set(decomposition) == {
        "passed", "uncontracted", "contracted", "accepted_micro_pair_decomposition",
        "new_macro_increment_decomposition", "RQy_identity", "B_internal_identity", "algebra",
    } and decomposition["passed"] is True
            and decomposition["uncontracted"] == joined["uncontracted_joined"]
            and decomposition["contracted"] == joined["joined"]
            and decomposition["accepted_micro_pair_decomposition"]["passed"] is True
            and decomposition["new_macro_increment_decomposition"]["passed"] is True
            and decomposition["RQy_identity"]["passed"] is True
            and decomposition["algebra"] ==
            "micro_pair + macro_increment = half(B_outer+D_micro+D_macro)-R_Q_y(m=1)+R_Q_y(first_micro lower)",
            "Joined uncontracted/contracted proof changed")

    for key in ("micro", "macro", "outer"):
        _pair(joined["reduced_boundary_B"][key], "joined B " + key)
        _pair(joined["reduced_volume_D"][key], "joined D " + key)
    seam = joined["shared_U_Q_B_seam"]
    require(set(seam) == {"micro_B_plus_D", "macro_B_plus_D", "outer_B_plus_D",
                          "same_source_physical_U_Q", "internal_endpoint_identity",
                          "ordinary_Z_orders"}
            and all(len(_pair(seam[key], "joined seam " + key)) == 2
                    for key in ("micro_B_plus_D", "macro_B_plus_D", "outer_B_plus_D"))
            and seam["same_source_physical_U_Q"] is True
            and seam["internal_endpoint_identity"]["passed"] is True
            and seam["ordinary_Z_orders"] == [0, 1],
            "Shared U/Q/B seam changed")
    reduced_guard = packet["source_guards"]["reduced_B_seam"]
    require(reduced_guard["micro_B_plus_D_rows"] == seam["micro_B_plus_D"]
            and reduced_guard["macro_B_plus_D_rows"] == seam["macro_B_plus_D"],
            "Reduced B seam guard does not bind joined source rows")
    contraction_rows = _seq(joined["rqy_boundary_contraction_rows"], "joined R_Q_y rows")
    require(len(contraction_rows) == 4
            and joined["endpoint_wrapper_identity"] == packet["source_guards"]["second_micro_endpoint_wrapper_identity"]["rows"],
            "Joined endpoint receipts changed")
    for order, row in enumerate(contraction_rows):
        require(set(row) == {"ordinary_Z_order", "handle", "negative_second_micro_upper",
                              "frozen_macro_lower", "row"}
                and row["ordinary_Z_order"] == order
                and row["handle"] == ports["signed_cancellation"][order]
                and row["negative_second_micro_upper"] == ports["second_micro_upper"][order]
                and row["frozen_macro_lower"] == ports["frozen_macro_lower"][order]
                and row["row"] == {
                    "operation": "sum",
                    "arguments": [ports["signed_cancellation"][order] - 1,
                                   ports["frozen_macro_lower"][order]],
                }, "Joined R_Q_y contraction row changed")

    identity = packet["accepted_rqy_seam_identity_consumed"]
    require(set(identity) == {
        "signed_cancellation_rows", "signed_cancellation_row_checks",
        "accepted_full_Q_y_rows", "parent_proof_rows", "second_micro_wrapper_identity",
        "common_physical_source", "R_Z_independent", "join_physical_source_binding",
        "identity_used_only_for_joined_expression",
    } and identity["signed_cancellation_rows"] == ports["signed_cancellation"]
            and len(_seq(identity["signed_cancellation_row_checks"], "identity contraction checks")) == 4
            and identity["second_micro_wrapper_identity"] == packet["source_guards"]["second_micro_endpoint_wrapper_identity"]["rows"]
            and identity["common_physical_source"] is True
            and identity["R_Z_independent"] is True
            and identity["join_physical_source_binding"] == packet["source_guards"]["join_physical_source_binding"]
            and identity["identity_used_only_for_joined_expression"] is True,
            "Accepted R_Q_y seam identity scope changed")
    require(identity["accepted_full_Q_y_rows"]
            == packet["source_guards"]["accepted_rqy_join"]["accepted_join_full_Q_y_rows"]
            and identity["accepted_full_Q_y_rows"],
            "Accepted full Q_y rows were dropped")
    proof_rows = _seq(identity["parent_proof_rows"], "joined parent proof rows")
    require(proof_rows == packet["source_guards"]["accepted_rqy_proof_rows"]
            and len(proof_rows) == 4,
            "Joined parent R_Q_y proof rows changed")


def _numeric_record(value, label):
    """Validate and decode one saved directed MPF interval record.

    The conversion is only for comparing serialized intervals; it does not
    evaluate a producer or reconstruct any source quantity.
    """
    require(isinstance(value, dict), label + " must be an interval record")
    require(set(value) == {
        "lower_exact_mpf_tuple", "upper_exact_mpf_tuple",
        "lower_display", "upper_display", "directed",
    }, label + " interval record fields changed")
    require(value["directed"] is True
            and isinstance(value["lower_display"], str)
            and isinstance(value["upper_display"], str),
            label + " interval direction/display changed")
    lower = value["lower_exact_mpf_tuple"]
    upper = value["upper_exact_mpf_tuple"]
    require(isinstance(lower, list) and isinstance(upper, list)
            and len(lower) == 4 and len(upper) == 4
            and all(type(item) is int for item in lower + upper),
            label + " exact MPF tuple changed")
    require(lower[0] in (0, 1) and upper[0] in (0, 1)
            and lower[1] >= 0 and upper[1] >= 0,
            label + " exact MPF tuple sign/mantissa changed")
    from mpmath import make_mpf
    lo = make_mpf(tuple(lower))
    hi = make_mpf(tuple(upper))
    require(lo <= hi, label + " lower endpoint exceeds upper endpoint")
    return lo, hi


def _numeric_rows(rows, label, count=2):
    rows = _seq(rows, label + " rows")
    require(len(rows) == count, label + " row count changed")
    for order, row in enumerate(rows):
        require(isinstance(row, dict)
                and set(row) == {"ordinary_Z_order", "value"}
                and row["ordinary_Z_order"] == order,
                label + " row schema changed")
        _numeric_record(row["value"], label + " Z" + str(order))
    return rows


def _validate_numeric_source_pins(packet, accepted):
    paths = (
        numeric_current.PARENT, numeric_current.B2_REPORT,
        numeric_current.I1_REPORT, numeric_current.UNIFORM,
        numeric_current.UNIFORM_SOURCE, numeric_current.LOWER_GRAPH,
        numeric_current.LOWER_SOURCE, numeric_current.LOOP_GRAPH,
        numeric_current.LOOP_SOURCE, numeric_current.BRIDGE_SOURCE,
        numeric_current.Q4_SOURCE, numeric_current.CORRELATED_SOURCE,
        numeric_current.JET_SOURCE, numeric_current.EA_SOURCE,
        numeric_current.CORE_SOURCE, numeric_current.CORE_CHECK,
        numeric_current.INLET_SOURCE, numeric_current.BRIDGE_RATES_SOURCE,
        numeric_current.BRIDGE_RATES_CHECK, numeric_current.MEAN_SOURCE,
    )
    expected = {}
    for path in paths:
        if path == numeric_current.PARENT:
            expected[path] = numeric_current.PARENT_SHA256
            continue
        group = "input_hashes" if path.startswith("experiments/root_st073/") else "implementation_hashes"
        key = Path(path).name if group == "input_hashes" else path
        digest = accepted[group].get(key)
        require(digest is not None, "Latest report has no numeric source hash for " + path)
        expected[path] = digest
    source_pins = packet["source_pins"]
    require(set(source_pins) == set(expected), "Numeric source pin set changed")
    for path, digest in expected.items():
        item = source_pins[path]
        require(set(item) == {"path", "expected_sha256"}
                and item["path"] == path
                and item["expected_sha256"] == digest,
                "Numeric source pin receipt changed: " + path)
        require(sha(ROOT / path) == digest, "Numeric source changed: " + path)
        if path != numeric_current.PARENT:
            key = Path(path).name if path.startswith("experiments/root_st073/") else path
            target = accepted["input_hashes"] if path.startswith("experiments/root_st073/") else accepted["implementation_hashes"]
            if key in target:
                require(target[key] == digest, "Conflicting inherited numeric pin: " + key)


def _validate_numeric_source_ast(source_ast):
    expected = {
        "axis": (numeric_current.EA_SOURCE, "source_jets", [
            '4*z[0]+j+sym(cap("Uz_minus_U0_absolute_bound",0,0))',
            '4+sym(cap("Uz_minus_U0_absolute_bound",0,1))',
        ]),
        "bridge_densities": (numeric_current.BRIDGE_SOURCE, "densities", [
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
        ]),
        "bridge_histories": (numeric_current.BRIDGE_SOURCE, "histories", [
            "for name,rate in RATES.items()",
            "original_incoming_memory_retained=True", "source_measure_applied_once=True",
        ]),
        "loop_transport": (numeric_current.LOOP_SOURCE, "transport", [
            "for name,rate in RATES.items():",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "predecessor_memory_never_reset=True",
        ]),
        "mean_phase": (numeric_current.MEAN_SOURCE, "_phase_rows", [
            '"phase_primitive_is_direct_ordinary_Z_row": True',
            '"source_operation": "phase_integral_of_slow_B_y" if source_name == "Cy" else "phase_integral_of_B"',
        ]),
        "mean_moving": (numeric_current.MEAN_SOURCE, "_moving_rows", [
            "exponent = graph.neg(graph.mul(hb, distance))",
            "integrand = graph.mul(kernel, primitive_at_t, hb)",
            '"physical_hb_weight_applied_once": True',
        ]),
        "mean_joined": (numeric_current.MEAN_SOURCE, "build_micro_mean", [
            'formula="exp(-hb*(s-1))"',
            "_scale_rows(g, first_JC_at_one, second_memory)",
            "_scale_rows(g, first_JCy_at_one, second_memory)",
            "second_complete = _add_rows(g, second_known, second_unresolved)",
        ]),
        "six_source": (numeric_current.CORRELATED_SOURCE, "build", [
            "source=dict(q=Q,q_Z=Q_Z,q_squared=mul(Q,Q),U_q_Z=mul(U,Q_Z),",
            'weight=fixed(g.mul(radius,radius,chart["physical_measure"]))',
        ]),
        "reduced_D": (numeric_current.CORRELATED_SOURCE, "build", [
            'D=add(neg(divide(integrals["q"])),half(divide(mul(h,z,integrals["q_Z"])))',
            'g.pscale(divide(mul(z,d,integrals["U_U_Z"]),L2),two))',
        ]),
        "physical_q": (numeric_current.Q4_SOURCE, "physical_pressure_fields", [
            "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))",
        ]),
        "core_inlet_average": (numeric_current.INLET_SOURCE, "__init__", [
            "sample=g.mul(rho,t)",
            "A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi)",
            "self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,",
        ]),
        "uniform_transfer": (numeric_current.UNIFORM_SOURCE, "run", [
            "derivative_rows=[]",
            "Uz_minus_U0_absolute_bound=eps*psi_norm*factor",
            "correction_C2=sum((eps*psi_norm*embedding(c,h,radius,0,k)",
        ]),
        "current_core_source": (numeric_current.CORE_SOURCE, "__init__", [
            "sigma=g.quotient(self.parameters['j'],c(500),'positive fixed scalar 500')",
            "logLambda=g.add(g.mul(c(4),logP),c(1000))",
            "Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))",
            "U=a.add(U0,a.scale(Psi,epsilon));mean=a.add(U0,a.scale(M,epsilon)",
        ]),
    }
    require(set(source_ast) == set(expected) | {"bridge_rates"},
            "Numeric source AST guard set changed")
    for key, (path, function, snippets) in expected.items():
        receipt = source_ast[key]
        require(isinstance(receipt, dict)
                and receipt.get("passed") is True
                and receipt.get("path") == path
                and receipt.get("function") == function
                and receipt.get("source_executed") is False
                and receipt.get("required_text") == snippets,
                "Numeric source AST guard changed: " + key)
    rates = source_ast.get("bridge_rates")
    require(isinstance(rates, dict)
            and rates.get("passed") is True
            and rates.get("path") == numeric_current.BRIDGE_RATES_SOURCE
            and rates.get("rates") == {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
            and rates.get("source_executed") is False,
            "Numeric bridge RATES source receipt changed")


def _validate_numeric_reduced_D(value):
    require(isinstance(value, dict)
            and set(value) == {
                "passed", "ordinary_Z_orders", "normalized_values", "half_normalized_values",
                "coefficientsC1rows", "same_W", "formula", "coefficient_definitions",
                "source_integrals", "finite_N_error_consumed_once",
                "background_and_finite_rows_are_signed_enclosures", "no_outer_boundary_claim",
            }
            and value["passed"] is True
            and value["ordinary_Z_orders"] == [0, 1]
            and value["same_W"] is True
            and value["formula"] == (
                "D=-I[q]/L+hZ*I[q_Z]/(2L)-I[q_squared]/2+"
                "d*I[U_q_Z]/L+2*hZ^2*I[U_squared]/L^2+"
                "2*Z*d*I[U_U_Z]/L^2"
            )
            and value["coefficient_definitions"] == {
                "h": "1-delta", "d": "1-Z^2", "L": "1-delta*Z^2",
                "hZ": "(1-delta)*Z", "L2": "L^2",
            }
            and value["source_integrals"] == "joint_I1.normalized_full_finite_N_integrals"
            and value["finite_N_error_consumed_once"] is True
            and value["background_and_finite_rows_are_signed_enclosures"] is True
            and value["no_outer_boundary_claim"] is True,
            "Numeric reduced-volume D receipt changed")
    _numeric_rows(value["normalized_values"], "Numeric reduced D normalized", 2)
    _numeric_rows(value["half_normalized_values"], "Numeric reduced D half normalized", 2)
    coefficients = value["coefficientsC1rows"]
    require(set(coefficients) == {
        "minus_inv_L", "hZ_over_2L", "minus_half", "d_over_L",
        "two_hZ2_over_L2", "two_Zd_over_L2",
    }, "Numeric reduced D coefficient set changed")
    for name in coefficients:
        _numeric_rows(coefficients[name], "Numeric reduced D coefficient " + name, 2)


def _validate_numeric_packet(packet, accepted):
    """Validate the signed six-I1 background packet and its joint radius."""
    _passed(packet, "Numeric signed background packet missing")
    require(set(packet) == {
        "passed", "source_family", "source_pins", "source_contract", "source_AST",
        "symbolic_controls", "domain", "source_axis", "signed_source_rows",
        "physical_scale", "geometry", "charts", "joint_I1", "finite_N_error",
        "reduced_micro_volume_D",
        "actual_signed_background_I1", "actual_signed_finite_N_I1_enclosures",
        "normalized_by_W", "W_definition", "finite_N_error_added_once_as_normalized_radius",
        "accepted_ancestors_replayed", "owner_builder_fixture_recurrence_norm_quadrature_replayed",
    }, "Numeric packet fields changed")
    require(packet["source_family"] == accepted["source_family"],
            "Numeric source family changed")
    _validate_numeric_source_pins(packet, accepted)

    contract = packet["source_contract"]
    mean_packet = accepted["actual_micro_mean_Z3_bounds"]
    history = mean_packet.get("new_history_controls")
    require(isinstance(history, dict)
            and history.get("passed") is True
            and history.get("physical_hb_applied_once") is True
            and history.get("second_memory_times_first_kernel_equals_joined_kernel") is True
            and history.get("ordinary_Z_commutes_with_phase_and_moving_integrals") is True,
            "Accepted explicit nonreset mean-history controls are required")
    require(set(contract) == {
        "passed", "actual_native_U_rows_are_signed_enclosures", "source_axis",
        "mean_equation", "exact_Q", "same_physical_measure", "sc_source_tuple",
        "sc_Z_independent", "finite_N_defect_packet_same_W",
        "actual_micro_I1_integrals_key", "I1_weight_provenance",
        "native_background_unit_receipt", "Q_increment_radius_receipt",
        "phase_mean_receipt", "saved_hb_scale_receipt", "bridge_rates_receipt",
        "joint_background_source_receipt", "current_core_direction_receipt",
        "current_core_fixed_point_receipt", "source_executed",
    } and contract["passed"] is True
            and contract["actual_native_U_rows_are_signed_enclosures"] is True
            and contract["source_axis"] == "U0=4*Z+j; physical U=U0+epsilon*Psi; j retained in the baseline"
            and contract["mean_equation"] == "M_y=U-M with positive exp(-h_b*(s-v)) Duhamel kernel"
            and contract["exact_Q"] == "Q=(2*Z*U-(1-delta)*Z*M-(1-Z^2)*M_Z)/(1-delta*Z^2)"
            and contract["same_physical_measure"] ==
            "R^2*h_b*ds exactly once (source receipt: R^2*hb*ds exactly once)"
            and contract["sc_Z_independent"] is True
            and contract["finite_N_defect_packet_same_W"] is True
            and contract["actual_micro_I1_integrals_key"] == "actual_micro_I1_defect_bounds"
            and contract["current_core_direction_receipt"] is True
            and contract["current_core_fixed_point_receipt"] is True
            and contract["source_executed"] is False,
            "Numeric source contract changed")
    weight = contract["I1_weight_provenance"]
    require(set(weight) == {
        "passed", "physical_measure_handle", "physical_measure_is_original_hb_not_integer_eight",
        "loghb_exact_tuples_equal_accepted_source_basis", "same_actual_radius_handle",
        "native_chart_widths_at_most_one", "actual_second_upper_number",
        "source_first_lower_is_sc_half", "positive_sc_upper_receipt_consumed",
        "Z_independent_measure_and_native_endpoints", "no_geometry_owner_or_hb_producer_executed",
        "source_ledger_sha256",
    } and weight["passed"] is True
            and weight["physical_measure_handle"] == 8
            and weight["physical_measure_is_original_hb_not_integer_eight"] is True
            and weight["loghb_exact_tuples_equal_accepted_source_basis"] is True
            and weight["same_actual_radius_handle"] == 68644
            and weight["native_chart_widths_at_most_one"] is True
            and weight["actual_second_upper_number"] == 2
            and weight["source_first_lower_is_sc_half"] is True
            and weight["positive_sc_upper_receipt_consumed"] is True
            and weight["Z_independent_measure_and_native_endpoints"] is True
            and weight["no_geometry_owner_or_hb_producer_executed"] is True
            and isinstance(weight["source_ledger_sha256"], str)
            and len(weight["source_ledger_sha256"]) == 64,
            "Numeric I1 weight provenance changed")
    require(contract["native_background_unit_receipt"] == {
        "saved_U_M_already_physical": True,
        "no_extra_Pstar_on_saved_U_M": True,
        "native_mean_is_M_not_R_times_M": True,
    }, "Numeric native background unit receipt changed")
    require(contract["Q_increment_radius_receipt"] == {
        "actual_micro_radius_handle": 68644,
        "exact_logLambda": "4*logPstar+1000",
        "physical_R_applied_once": True,
        "physical_R_Z_independent": True,
    }, "Numeric Q-increment radius receipt changed")
    require(contract["phase_mean_receipt"] == {
        "second_memory_formula": "exp(-h_b*(s-1))",
        "second_memory_formula_source": "exp(-hb*(s-1))",
        "second_history_not_reset": True,
        "new_history_controls_passed": True,
        "physical_hb_weight_applied_once": True,
    }, "Numeric phase-mean receipt changed")
    hb_receipt = contract["saved_hb_scale_receipt"]
    require(set(hb_receipt) == {
        "source_key", "saved_tuple_is_loghb", "physical_hb", "log_W",
    } and hb_receipt["source_key"] ==
            "actual_micro_mean_Z3_bounds.physical_scale_logs.hb"
            and hb_receipt["saved_tuple_is_loghb"] is True
            and hb_receipt["physical_hb"] in {"exp(loghb)", "h_b=exp(loghb)"}
            and hb_receipt["log_W"] == "2*log_Ra+loghb",
            "Numeric hb scale receipt changed")
    require(contract["bridge_rates_receipt"] == packet["source_AST"]["bridge_rates"],
            "Numeric bridge rates receipt changed")
    sc = contract["sc_source_tuple"]
    require(isinstance(sc, list) and len(sc) == 4 and all(type(item) is int for item in sc),
            "Numeric sc source tuple changed")
    receipt = contract["joint_background_source_receipt"]
    require(set(receipt) == {
        "same_native_background_on_both_micro_charts", "source_graph_functions_not_cap_values",
        "six_signed_defect_handles", "joint_endpoint_weights",
    } and receipt["same_native_background_on_both_micro_charts"] is True
            and receipt["source_graph_functions_not_cap_values"] is True
            and receipt["joint_endpoint_weights"] == {"lower": 1, "upper": -1},
            "Numeric joint background source receipt changed")
    handles = receipt["six_signed_defect_handles"]
    require(set(handles) == {"first_micro", "second_micro"},
            "Numeric six-source chart handles changed")
    for chart in handles.values():
        require(set(chart) == set(numeric_current.NAMES),
                "Numeric six-source names changed")
        for name in numeric_current.NAMES:
            _pair(chart[name], "Numeric source handle " + name)

    _validate_numeric_source_ast(packet["source_AST"])
    symbolic = packet["symbolic_controls"]
    require(set(symbolic) == {
        "passed", "Q_source_identity", "mean_source_identity",
        "positive_Duhamel_mass", "saved_hb_field", "ordinary_factorial_normalization",
    } and symbolic["passed"] is True
            and symbolic["Q_source_identity"] is True
            and symbolic["mean_source_identity"] == "M_y=U-M"
            and symbolic["positive_Duhamel_mass"] == "1-exp(-h_b*s), h_b=exp(loghb)"
            and symbolic["saved_hb_field"] == "loghb"
            and symbolic["ordinary_factorial_normalization"] is False,
            "Numeric symbolic controls changed")

    domain = packet["domain"]
    require(domain == {
        "s": ["sc/2", "1", "2"],
        "charts": ["first_micro", "second_micro"],
        "Z": ["37/100", "93/250"],
        "ordinary_Z_orders": [0, 1],
        "physical_y": "log(R/Ra)=h_b*s; h_b=exp(loghb)",
    }, "Numeric domain changed")
    axis = packet["source_axis"]
    require(set(axis) == {
        "U0", "physical_U", "baseline", "j", "ordinary_Z_orders",
        "source_rows_are_signed_enclosures", "axis_reference_not_selected_as_actual_function",
    } and axis["U0"] == "4*Z+j"
            and axis["physical_U"] == "U0+epsilon*Psi"
            and axis["ordinary_Z_orders"] == [0, 1, 2, 3]
            and axis["source_rows_are_signed_enclosures"] is True
            and axis["axis_reference_not_selected_as_actual_function"] is True,
            "Numeric source axis changed")
    _numeric_rows(axis["baseline"], "Numeric baseline", 4)
    _numeric_record(axis["j"], "Numeric j")

    signed = packet["signed_source_rows"]
    require(set(signed) == {
        "U", "M", "Q", "integrands", "core_correction", "native_y_error",
        "M_error_radius", "M3_norm_radius", "M3_center_is_zero",
        "M0_M1_M2_are_axis_plus_symmetric_error",
        "signed_function_values_are_enclosures_not_cap_centers",
    } and signed["M3_center_is_zero"] is True
            and signed["M0_M1_M2_are_axis_plus_symmetric_error"] is True
            and signed["signed_function_values_are_enclosures_not_cap_centers"] is True,
            "Numeric signed source rows changed")
    _numeric_rows(signed["U"], "Numeric U", 4)
    _numeric_rows(signed["M"], "Numeric M", 4)
    _numeric_rows(signed["Q"], "Numeric Q", 3)
    _numeric_rows(signed["core_correction"], "Numeric core correction", 4)
    _numeric_rows(signed["native_y_error"], "Numeric native y error", 4)
    m_errors = _seq(signed["M_error_radius"], "Numeric M error radius")
    require(len(m_errors) == 4, "Numeric M error radius row count changed")
    for order, value in enumerate(m_errors):
        _numeric_record(value, "Numeric M error radius Z" + str(order))
    _numeric_record(signed["M3_norm_radius"], "Numeric M3 norm radius")
    integrands = signed["integrands"]
    require(set(integrands) == set(numeric_current.NAMES), "Numeric integrand names changed")
    for name in numeric_current.NAMES:
        _numeric_rows(integrands[name], "Numeric integrand " + name, 2)

    scale = packet["physical_scale"]
    require(set(scale) == {
        "logPstar", "loghb", "log_Ra", "log_W", "W_definition", "radius_definition",
        "epsilon_definition", "saved_hb_is_loghb", "physical_hb_definition",
        "log_W_definition", "no_astronomical_exp_materialized",
    } and scale["W_definition"] == "Ra^2*h_b"
            and scale["radius_definition"] == "R=Ra*exp(h_b*s)"
            and scale["epsilon_definition"] == "exp(-(4*logPstar+1000))"
            and scale["saved_hb_is_loghb"] is True
            and scale["physical_hb_definition"] in {"exp(loghb)", "h_b=exp(loghb)"}
            and scale["log_W_definition"] == "2*log_Ra+loghb"
            and scale["no_astronomical_exp_materialized"] is True,
            "Numeric physical scale changed")
    for key in ("logPstar", "loghb", "log_Ra", "log_W"):
        _numeric_record(scale[key], "Numeric physical scale " + key)

    geometry = packet["geometry"]
    require(set(geometry) == {
        "passed", "loghb", "saved_hb_field", "physical_hb", "geometry_factor",
        "first_length", "second_length", "first_lower", "seam", "upper",
        "geometry_remainder_proof", "physical_measure_reduced_once",
    } and geometry["passed"] is True
            and geometry["first_length"] == "1-sc/2"
            and geometry["second_length"] == "1"
            and geometry["first_lower"] == "sc/2"
            and geometry["seam"] == "s=1"
            and geometry["upper"] == "2"
            and geometry["geometry_remainder_proof"] ==
            "loghb<=-1024*log(2) => h_b<=2^-1024; exp(2*h_b*s)<=1+2^-1021 for 0<=s<=2"
            and geometry["saved_hb_field"] == "loghb"
            and geometry["physical_hb"] == "h_b=exp(loghb)"
            and geometry["physical_measure_reduced_once"] is True,
            "Numeric geometry receipt changed")
    _numeric_record(geometry["loghb"], "Numeric geometry hb log")
    _numeric_record(geometry["geometry_factor"], "Numeric geometry factor")

    charts = packet["charts"]
    require(set(charts) == {"first_micro", "second_micro"}, "Numeric chart set changed")
    for chart_name in ("first_micro", "second_micro"):
        chart = charts[chart_name]
        require(set(chart) == {
            "normalized_integrals", "chart_length", "geometry_factor",
            "normalized_measure", "background_is_signed_enclosure",
            "no_midpoint_or_cap_selected_as_function_value",
        } and chart["normalized_measure"] == "integral exp(2*h_b*s) ds; h_b=exp(loghb)"
                and chart["background_is_signed_enclosure"] is True
                and chart["no_midpoint_or_cap_selected_as_function_value"] is True,
                "Numeric chart receipt changed: " + chart_name)
        integrals = chart["normalized_integrals"]
        require(set(integrals) == set(numeric_current.NAMES),
                "Numeric chart integrand set changed: " + chart_name)
        for name in numeric_current.NAMES:
            _numeric_rows(integrals[name], chart_name + " " + name, 2)
        _numeric_record(chart["chart_length"], chart_name + " length")
        _numeric_record(chart["geometry_factor"], chart_name + " geometry factor")

    joint = packet["joint_I1"]
    require(set(joint) == {
        "normalized_background_integrals", "normalized_full_finite_N_integrals",
        "background_is_signed_enclosure", "finite_N_error_radius_is_symmetric",
        "finite_N_error_added_once", "same_physical_measure",
    } and joint["background_is_signed_enclosure"] is True
            and joint["finite_N_error_radius_is_symmetric"] is True
            and joint["finite_N_error_added_once"] is True
            and joint["same_physical_measure"] ==
            "R^2*h_b*ds exactly once (source receipt: R^2*hb*ds exactly once)",
            "Numeric joint I1 receipt changed")
    background = joint["normalized_background_integrals"]
    full = joint["normalized_full_finite_N_integrals"]
    require(set(background) == set(numeric_current.NAMES)
            and set(full) == set(numeric_current.NAMES),
            "Numeric joint six-source set changed")
    for name in numeric_current.NAMES:
        bg_rows = _numeric_rows(background[name], "Numeric joint background " + name, 2)
        full_rows = _numeric_rows(full[name], "Numeric joint finite N " + name, 2)
        for order, (bg_row, full_row) in enumerate(zip(bg_rows, full_rows)):
            bg_lo, bg_hi = _numeric_record(bg_row["value"], "Numeric joint background " + name + " Z" + str(order))
            full_lo, full_hi = _numeric_record(full_row["value"], "Numeric joint finite N " + name + " Z" + str(order))
            require(full_lo <= bg_lo and bg_hi <= full_hi,
                    "Numeric finite N interval does not contain background: " + name)

    _validate_numeric_reduced_D(packet["reduced_micro_volume_D"])

    finite = packet["finite_N_error"]
    require(set(finite) == {
        "passed", "normalized_error_radius", "normalized_error_radius_definition",
        "source_log_cap_minus_log_W", "all_log_caps_below_threshold", "threshold",
        "same_source_function_error_not_a_selected_value", "physical_error_packet_not_recomputed",
    } and finite["passed"] is True
            and finite["normalized_error_radius_definition"] ==
            "2^-1024 in W=Ra^2*h_b units; saved hb field is loghb"
            and finite["all_log_caps_below_threshold"] is True
            and finite["same_source_function_error_not_a_selected_value"] is True
            and finite["physical_error_packet_not_recomputed"] is True,
            "Numeric finite N error receipt changed")
    radius = finite["normalized_error_radius"]
    _numeric_record(radius, "Numeric normalized error radius")
    require(radius["lower_exact_mpf_tuple"] == [0, 1, -1024, 1]
            and radius["upper_exact_mpf_tuple"] == [0, 1, -1024, 1],
            "Numeric normalized error radius is not exactly 2^-1024")
    _numeric_record(finite["threshold"], "Numeric finite N threshold")
    ratios = finite["source_log_cap_minus_log_W"]
    require(set(ratios) == set(numeric_current.NAMES), "Numeric finite N cap names changed")
    for name in numeric_current.NAMES:
        rows = _seq(ratios[name], "Numeric log cap minus log W " + name)
        require(len(rows) == 2, "Numeric log cap minus log W row count changed")
        for order, value in enumerate(rows):
            _numeric_record(value, "Numeric log cap minus log W " + name + " Z" + str(order))

    require(packet["actual_signed_background_I1"] is True
            and packet["actual_signed_finite_N_I1_enclosures"] is True
            and packet["normalized_by_W"] is True
            and packet["W_definition"] ==
            "W=Ra^2*h_b with h_b=exp(loghb); physical I1=W*normalized_integral"
            and packet["finite_N_error_added_once_as_normalized_radius"] is True
            and packet["accepted_ancestors_replayed"] is False
            and packet["owner_builder_fixture_recurrence_norm_quadrature_replayed"] is False,
            "Numeric signed/finite N scope changed")
    return dict(
        passed=True,
        actual_signed_background_six_I1_C1=True,
        actual_signed_finite_N_joint_radius=True,
        normalized_by_W=True,
        finite_N_radius_added_once=True,
        inherited_parent_math_payload_preserved=True,
        full_macro_pressure_global_closure_open=True,
        accepted_ancestors_replayed=False,
    )


def validate_packet(packet, accepted):
    _passed(packet, "Correlated merged packet missing")
    require(set(packet) == {
        "passed", "parent_report", "parent_report_sha256", "source_pins", "input_pins",
        "source_guards", "graph_namespace", "safe_graphs", "remapped_accepted_seam_target_ports",
        "merged_namespace",
        "macro_six_I1_source", "actual_macro_six_I1_source", "remaining_outer_macro_RQy",
        "actual_joined_pressure_function", "joined_pressure_function",
        "accepted_rqy_seam_identity_consumed", "old_212_terms_and_errors_unchanged",
        "old_six_numeric_terms_and_all_error_caps_unchanged", "signed_function_values_evaluated",
        "signed_six_integrals_evaluated", "global_goal_complete", "accepted_ancestors_replayed",
        "numerical_phase_or_quadrature_evaluated",
    }, "Correlated packet fields changed")
    require(packet["parent_report"] == PARENT
            and packet["parent_report_sha256"] == PARENT_SHA256,
            "Correlated packet parent binding changed")
    seam = _validate_parent_seam(accepted)
    correlated = _load(ROOT / CORRELATED)
    correlated_packet = correlated["actual_micro_correlated_m15"]
    _validate_source_pins(packet, accepted)
    _validate_source_guards(packet, correlated_packet, seam)
    _validate_graph_namespace(packet)
    ports = packet["remapped_accepted_seam_target_ports"]
    require(set(ports) == {"second_micro_upper", "frozen_macro_lower", "signed_cancellation"},
            "Correlated seam target ports changed")
    for key in ports:
        require(isinstance(ports[key], list) and all(type(v) is int for v in ports[key]),
                "Correlated seam target port malformed: " + key)
    _validate_safe_graphs(packet, ports)
    _validate_merged_namespace(packet)
    macro = packet["macro_six_I1_source"]
    _validate_macro(packet, ports)
    _validate_join(packet, macro, ports, correlated_packet)
    require(packet["remaining_outer_macro_RQy"] == macro["outer_RQy"]
            and packet["old_212_terms_and_errors_unchanged"] is True
            and packet["old_six_numeric_terms_and_all_error_caps_unchanged"] is True
            and packet["signed_function_values_evaluated"] is False
            and packet["signed_six_integrals_evaluated"] is False
            and packet["global_goal_complete"] is False
            and packet["accepted_ancestors_replayed"] is False
            and packet["numerical_phase_or_quadrature_evaluated"] is False,
            "Correlated source-only scope changed")
    return dict(
        passed=True,
        merged_graph_namespace=True,
        relocated_rqy_suffix=True,
        exact_macro_six_I1_source=True,
        whole_micro_pair_consumed_once=True,
        macro_outer_RQy_retained=True,
        inherited_parent_math_payload_preserved=True,
        quantitative_signed_values_and_global_closure_open=True,
        accepted_ancestors_replayed=False,
    )


def _packet_fields(report, accepted):
    require(report.get("parent_report") == PARENT_NAME
            and report.get("parent_report_sha256") == PARENT_SHA256,
            "Same R_Q_y seam predecessor required")
    mutable = {
        "stage", "parent_report", "parent_report_sha256", "input_hashes",
        "implementation_hashes", "contract", "checks", "execution_seconds",
    }
    for key, value in accepted.items():
        if key not in mutable:
            require(key in report and report[key] == value,
                    "Inherited parent packet or error changed: " + key)
    return dict(
        correlated=validate_packet(report.get("actual_micro_macro_correlated_m15"), accepted),
        numeric=_validate_numeric_packet(report.get("actual_micro_background_I1_numeric"), accepted),
    )


def _merge_pins(accepted, *packets):
    inputs = dict(accepted["input_hashes"])
    implementation = dict(accepted["implementation_hashes"])
    for packet in packets:
        for path, item in packet["source_pins"].items():
            require(item["path"] == path
                    and ("passed" not in item or item["passed"] is True),
                    "Malformed new source pin: " + path)
            digest = item["expected_sha256"]
            if path.startswith("experiments/root_st073/"):
                key, target = Path(path).name, inputs
            else:
                key, target = path, implementation
            require(target.get(key, digest) == digest,
                    "Conflicting inherited source pin: " + key)
            target[key] = digest
    implementation.update({path: sha(ROOT / path) for path in CODE})
    return inputs, implementation


def _execution_audit():
    """Record the one excluded review-time legacy symbolic check explicitly."""
    return dict(
        new_production_accepted_ancestors_replayed=False,
        new_production_old_module_executed=False,
        review_accepted_symbolic_check_replayed_once=True,
        review_accepted_symbolic_check_function=(
            "src/openai_ns_reconstruction/original_n1_micro_correlated_m15.py:exact_algebra"
        ),
        review_accepted_symbolic_check_result="6 true",
        review_accepted_symbolic_check_result_disregarded=True,
        old_data_modified=False,
        accepted_producers_norms_recurrences_quadrature_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING"
            and NUMERIC_REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed correlated module may produce this report")
    require(sha(ROOT / NUMERIC) == NUMERIC_REVIEWED_MODULE_SHA256,
            "Only the independently reviewed numeric module may produce this report")
    accepted = parent()
    pins(accepted)
    packet = build(ROOT, accepted)
    require(accepted == parent(), "Correlated construction mutated its accepted parent")
    validate_packet(packet, accepted)
    numeric_packet = numeric_build(ROOT, accepted)
    require(accepted == parent(), "Numeric construction mutated its accepted parent")
    _validate_numeric_packet(numeric_packet, accepted)
    inputs, implementation = _merge_pins(accepted, packet, numeric_packet)
    report = copy.deepcopy(accepted)
    report.update(
        stage=(
            "true_micro_macro_correlated_six_I1_source_plus_signed_micro_"
            "background_finite_N_and_reduced_D_numeric"),
        parent_report=PARENT_NAME,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_macro_correlated_m15=packet,
        actual_micro_background_I1_numeric=numeric_packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
        execution_audit=_execution_audit(),
    )
    report["contract"] = dict(
        true_merged_first_second_micro_and_macro_source=True,
        relocated_rqy_suffix_reference_safe=True,
        exact_macro_six_I1_R_squared_length_once=True,
        whole_micro_pair_consumed_once=True,
        first_micro_lower_and_macro_outer_RQy_retained=True,
        old_212_terms_and_all_error_caps_retained=True,
        actual_signed_micro_background_six_I1_C1=True,
        actual_signed_micro_finite_N_joint_radius=True,
        actual_reduced_micro_volume_D_C1=True,
        full_macro_pressure_global_closure=False,
        graph_packet_numerical=False,
        macro_signed_I1_and_full_pressure_global_closure_open=True,
        new_production_accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Merged correlated first/second-micro source, relocated R_Q_y seam, "
            "macro six-I1 R-squared-length source, and retained outer R_Q_y; "
            "plus actual signed six-I1 background and joint finite-N radius in W units; "
            "plus reduced micro-volume D C1 in the same normalized units; "
            "full macro/pressure/global closure and graph-packet numerics remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_numeric_module_sha256=sha(ROOT / NUMERIC),
        independent_runtime_worker_model_metadata_available=False,
        execution_audit=_execution_audit(),
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report = _load(path)
    accepted = parent()
    checks = report["checks"]
    require(checks["passed"] is True
            and checks.get("independent_new_math_review_completed") is True
            and checks.get("independent_review_completed") is True
            and checks.get("independently_reviewed_module_sha256")
            == report["implementation_hashes"][CODE[0]]
            == REVIEWED_MODULE_SHA256
            and checks.get("independently_reviewed_numeric_module_sha256")
            == report["implementation_hashes"][NUMERIC]
            == NUMERIC_REVIEWED_MODULE_SHA256
            and report["implementation_hashes"][CODE[2]] == sha(ROOT / CODE[2]),
            "Same independently reviewed correlated and numeric implementations required")
    require(report.get("execution_audit") == _execution_audit(),
            "Execution audit is missing or changed")
    contract = report.get("contract", {})
    require(contract.get("new_production_accepted_ancestors_replayed") is False
            and contract.get("full_macro_pressure_global_closure") is False
            and contract.get("graph_packet_numerical") is False
            and contract.get("actual_signed_micro_background_six_I1_C1") is True
            and contract.get("actual_signed_micro_finite_N_joint_radius") is True
            and contract.get("actual_reduced_micro_volume_D_C1") is True
            and contract.get("macro_signed_I1_and_full_pressure_global_closure_open") is True,
            "Top-level contract scope changed")
    return dict(
        passed=True,
        pins=pins(report),
        packet_invariants=_packet_fields(report, accepted),
        execution_audit=report["execution_audit"],
        serialized_only=True,
        new_production_accepted_owner_producer_builder_norm_recurrence_or_quadrature_replayed=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        report = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        temporary.replace(path)
        print(json.dumps(dict(
            passed=True,
            execution_seconds=report["execution_seconds"],
            bytes=path.stat().st_size,
            sha256=sha(path),
        )))
    else:
        print(json.dumps(inspect(path), indent=2))
