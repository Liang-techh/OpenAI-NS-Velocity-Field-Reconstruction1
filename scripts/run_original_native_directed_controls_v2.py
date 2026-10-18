"""Publish and inspect the standalone corrected native-control packet.

Inspection intentionally uses only JSON, hashes, AST-free structural checks,
and exact integer comparisons of serialized MPF tuples.  It does not import a
mathematical provider, execute a builder, run a fixture, or replay an older
stage.  Production is refused until the reviewed source hashes are filled in
and the coordinator's fail-closed packet checks can run.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
import sys
import time

from mpmath.libmp import mpf_cmp


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PARENT_NAME = (
    "lei_ren_part1_paper_compliant_current_original_native_control_directed_repair.json"
)
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_SHA256 = "2ba16bd7c1bfd2fae90a5f13c69736d2764cccc049f8737e653fd10f9f8dbf3c"
PARENT_BYTES = 39425302
NAME = "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
REPORT_PATH = BASE / NAME
RUNNER = "scripts/run_original_native_directed_controls_v2.py"
MODULE = "src/openai_ns_reconstruction/original_native_directed_controls_v2.py"
COMMON = "src/openai_ns_reconstruction/original_native_directed_common.py"
LONG = "src/openai_ns_reconstruction/original_long_reference_directed_controls.py"
RESTORATION = "src/openai_ns_reconstruction/original_restoration_directed_controls.py"
KEY = "actual_native_directed_controls_v2"
SCHEMA = "original_native_directed_controls_v2_v1"
ORDERS = tuple(range(5))
PARENT_INPUT_ALIASES = {
    "../../src/openai_ns_reconstruction/schedule_pressure.py",
}

# Deliberately fail closed while the source review is pending.  The root agent
# may replace these values with the reviewed SHA-256 values before the one
# authorized production run.
REVIEWED_SOURCE_HASHES = {
    COMMON: "5aa14af159faba82da6e9a73c21c49f83ec9003bf56cd40474301908c4882cb6",
    LONG: "29b5d50c4ea558a08974eb42e3a5f2746b121f9d2565685f80b5d726d7861116",
    RESTORATION: "6255fd90689e2557f5c1493087cdcc643d2ce5e1fec4ad61c128f989d94e4d2e",
    MODULE: "92dbef6efc083f1642278bf1d6881696fac0021b272775f1b0ddfa78d7a83ae2",
}


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def _safe_relative(path):
    candidate = Path(path)
    require(not candidate.is_absolute(), "Unsafe source pin path")
    if ".." in candidate.parts:
        require(candidate.as_posix() in PARENT_INPUT_ALIASES,
                "Unregistered parent traversal source pin")
    return candidate.as_posix()


def _input_path(name):
    candidate = Path(name)
    if candidate.parts and candidate.parts[0] == "experiments":
        return candidate.as_posix()
    if ".." in candidate.parts:
        return candidate.as_posix()
    return (Path("experiments/root_st073") / candidate.name).as_posix()


def _check_parent(root, parent):
    path = root / PARENT_PATH
    require(path.is_file(), "Directed-controls parent report is missing")
    raw = path.read_bytes()
    require(len(raw) == PARENT_BYTES, "Directed-controls parent byte count changed")
    require(hashlib.sha256(raw).hexdigest() == PARENT_SHA256,
            "Directed-controls parent bytes changed")
    require(json.loads(raw) == parent, "Parent mapping is not canonical")
    count = 0
    for name, expected in parent.get("input_hashes", {}).items():
        relative = _input_path(name)
        source = root / Path("experiments/root_st073") / name \
            if ".." in Path(name).parts else root / relative
        require(source.is_file(), "Pinned input is missing: " + relative)
        require(sha(source) == expected, "Pinned input changed: " + relative)
        count += 1
    for relative, expected in parent.get("implementation_hashes", {}).items():
        relative = _safe_relative(relative)
        source = root / relative
        require(source.is_file(), "Pinned implementation is missing: " + relative)
        require(sha(source) == expected, "Pinned implementation changed: " + relative)
        count += 1
    require(count == len(parent.get("input_hashes", {}))
            + len(parent.get("implementation_hashes", {})),
            "Not all parent pins were checked")
    integrity = parent.get("actual_abs_log_integrity", {})
    gate = integrity.get("self_hash_fixture_receipt", {}).get("gate_fixture", {})
    require(integrity.get("inherited_log_cap_certification_verified") is False
            and integrity.get("revalidation_required") is True
            and gate.get("no_revalidation_receipt_admitted_yet") is True,
            "Historical abs-log integrity guard was weakened")
    return count


def _canonical_tuple(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4,
            "Malformed MPF endpoint tuple")
    require(all(type(part) is int for part in value),
            "MPF endpoint components must be integers")
    sign, mantissa, exponent, bit_count = value
    require(sign in (0, 1) and mantissa >= 0,
            "Malformed MPF sign or mantissa")
    require(bit_count == mantissa.bit_length(), "Malformed MPF bit count")
    if mantissa == 0:
        require(tuple(value) == (0, 0, 0, 0), "Noncanonical MPF zero tuple")
    else:
        require(mantissa & 1, "Noncanonical MPF mantissa")
    return tuple(value)


def _mpf_compare(left, right):
    """Compare finite MPF tuples exactly, with no decimal conversion."""
    left = _canonical_tuple(left)
    right = _canonical_tuple(right)
    # Compare top binary exponents before any bounded mantissa alignment.
    # Mantissa precision does not bound the distance between the exponents.
    return mpf_cmp(left, right)


def _record(value, label="interval"):
    require(isinstance(value, Mapping) and value.get("directed") is True,
            "Directed interval record required: " + label)
    lower = _canonical_tuple(value.get("lower_exact_mpf_tuple"))
    upper = _canonical_tuple(value.get("upper_exact_mpf_tuple"))
    require(_mpf_compare(lower, upper) <= 0,
            "Reversed directed interval: " + label)
    return 1


def _record_count(value, label="packet"):
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            return _record(value, label)
        return sum(_record_count(item, label + "." + str(key))
                   for key, item in value.items())
    if isinstance(value, list):
        return sum(_record_count(item, label + "[" + str(index) + "]")
                   for index, item in enumerate(value))
    return 0


def _validate_ordered(value, label="packet"):
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = label + "." + key_text
            if key_text in {"ordinary_Z_orders", "ordinary_Z_order_list"}:
                require(isinstance(item, list)
                        and item == list(range(len(item))),
                        "Unordered ordinary-Z list at " + child)
            if key_text == "p2_ordinary_Z_orders":
                require(item == [0, 1, 2, 3],
                        "p2 ordinary-Z list must have four rows at " + child)
            if isinstance(item, list) and item and all(
                    isinstance(row, Mapping) and "ordinary_Z_order" in row
                    for row in item
            ):
                orders = [row["ordinary_Z_order"] for row in item]
                require(orders == list(range(len(item))),
                        "Unordered ordinary-Z records at " + child)
            _validate_ordered(item, child)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_ordered(item, label + "[" + str(index) + "]")


def _validate_execution_flags(value, label="packet"):
    if isinstance(value, Mapping):
        for key, item in value.items():
            name = str(key).lower()
            if name in {
                "source_executed", "accepted_ancestors_replayed",
                "accepted_builders_called", "accepted_builders_or_owners_replayed",
                "accepted_checkers_called", "accepted_checkers_replayed",
                "accepted_self_checks_called", "accepted_generic_self_checks_called",
                "accepted_recurrences_replayed", "accepted_quadratures_replayed",
                "old_stage_calls", "old_stage_called", "old_builders_called",
            }:
                require(item is False,
                        "Forbidden execution flag at " + label + "." + str(key))
            _validate_execution_flags(item, label + "." + str(key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_execution_flags(item, label + "[" + str(index) + "]")


def _validate_source_pin(root, relative, pin):
    require(isinstance(pin, Mapping), "Malformed source pin: " + str(relative))
    require(set(pin) == {
        "path", "expected_sha256", "actual_sha256", "passed", "source_executed"
    }, "Noncanonical source pin fields: " + str(relative))
    path = _safe_relative(str(relative))
    require(pin["path"] == path and pin["passed"] is True
            and pin["source_executed"] is False,
            "Noncanonical source pin: " + path)
    parent_alias = ".." in Path(path).parts
    actual_path = root / Path("experiments/root_st073") / path if parent_alias else root / path
    require(not parent_alias or path in PARENT_INPUT_ALIASES,
            "Unregistered parent traversal source pin: " + path)
    actual = sha(actual_path)
    require(actual == pin["expected_sha256"] == pin["actual_sha256"],
            "Source pin changed: " + path)
    return path, actual


def _merge_pins(root, parent, packet):
    inputs = dict(parent["input_hashes"])
    implementations = dict(parent["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    source_pins = packet.get("source_pins")
    require(isinstance(source_pins, Mapping), "Packet source pins are missing")
    for relative, pin in source_pins.items():
        path, actual = _validate_source_pin(root, relative, pin)
        candidate = Path(path)
        if path in parent["input_hashes"] and ".." in candidate.parts:
            target = inputs
            name = path
        elif candidate.parts and candidate.parts[0] == "experiments":
            target = inputs
            name = candidate.name
        else:
            target = implementations
            name = path
        if name in target:
            require(target[name] == actual, "Historical source pin conflicts: " + path)
        else:
            target[name] = actual
    # The packet must pin every accepted parent source exactly.  This catches
    # accidental omission of a safe helper or old source proof.
    for name, expected in parent["input_hashes"].items():
        path = _input_path(name)
        pin = source_pins.get(path)
        require(pin is not None and pin["actual_sha256"] == expected,
                "Accepted input pin was not frozen: " + path)
    for path, expected in parent["implementation_hashes"].items():
        pin = source_pins.get(path)
        require(pin is not None and pin["actual_sha256"] == expected,
                "Accepted implementation pin was not frozen: " + path)
    return inputs, implementations


def _strict_negative(value, label):
    if not isinstance(value, Mapping) or "upper_exact_mpf_tuple" not in value:
        return False
    _record(value, label)
    return _mpf_compare(value["upper_exact_mpf_tuple"], (0, 0, 0, 0)) < 0


def _strict_controls(packet):
    receipt = packet.get("strict_control_receipt")
    require(isinstance(receipt, Mapping) and receipt.get("passed") is True,
            "Strict-control receipt is missing")
    require(receipt.get("exact_tuple_comparisons") is True,
            "Strict-control receipt did not use exact tuple comparisons")
    slopes = receipt.get("strict_slope_records")
    margins = receipt.get("negative_b_margin_records")
    require(isinstance(slopes, list) and slopes,
            "Strict slope records are missing")
    require(isinstance(margins, list) and margins,
            "Negative b-margin records are missing")
    for path in slopes:
        require(str(path).startswith("raw_interval_tree."),
                "Strict slope path is outside the raw tree")
        value = packet
        for part in str(path).replace("[", ".").replace("]", "").split("."):
            if part:
                value = value[int(part)] if isinstance(value, list) else value[part]
        require(_strict_negative(value, path), "Strict slope is not negative: " + path)
    for path in margins:
        require(str(path).startswith("raw_interval_tree."),
                "b-margin path is outside the raw tree")
        value = packet
        for part in str(path).replace("[", ".").replace("]", "").split("."):
            if part:
                value = value[int(part)] if isinstance(value, list) else value[part]
        require(_strict_negative(value, path), "b margin is not negative: " + path)


def _serialized_native_cap_rows(rows, label):
    require(isinstance(rows, list) and len(rows) == 5,
            "Serialized cap is not five rows: " + label)
    raw = all(row is None or (
        isinstance(row, Mapping)
        and row.get("directed") is True
        and "lower_exact_mpf_tuple" in row
        and "upper_exact_mpf_tuple" in row
    ) for row in rows)
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    packed = all(isinstance(row, Mapping) and set(row) == required for row in rows)
    require(raw or packed, "Serialized cap is neither raw nor explicitly packed: " + label)
    if raw:
        for index, row in enumerate(rows):
            if row is not None:
                _record(row, label + "[" + str(index) + "]")
        return "raw"
    for index, row in enumerate(rows):
        require(row["ordinary_Z_order"] == index
                and row["exact_zero"] in (True, False),
                "Serialized packed cap order changed: " + label)
        bound = row["log_absolute_upper_bound"]
        if row["exact_zero"]:
            require(bound is None,
                    "Serialized exact-zero cap has a bound: " + label)
        else:
            _record(bound, label + "[" + str(index) + "]")
    return "packed"


def _serialized_packed_rows(rows, count, label):
    require(isinstance(rows, list) and len(rows) == count,
            "Serialized public rows changed: " + label)
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping) and set(row) == required
                and row["ordinary_Z_order"] == index
                and row["exact_zero"] in (True, False),
                "Serialized public rows are not canonical: " + label)
        bound = row["log_absolute_upper_bound"]
        if row["exact_zero"]:
            require(bound is None, "Serialized exact-zero row has a bound: " + label)
        else:
            _record(bound, label + "[" + str(index) + "]")


def _positive_record(value, label):
    _record(value, label)
    require(_mpf_compare(value["lower_exact_mpf_tuple"], (0, 0, 0, 0)) > 0,
            "Positive lower endpoint is not strict: " + label)


def _negative_packed_row(row, label):
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    require(isinstance(row, Mapping) and set(row) == required
            and row.get("exact_zero") is False,
            "Nonzero packed row required: " + label)
    bound = row["log_absolute_upper_bound"]
    _record(bound, label)
    require(_mpf_compare(bound["upper_exact_mpf_tuple"], (0, 0, 0, 0)) < 0,
            "Packed log cap is not negative: " + label)


def _validate_native_shapes(raw):
    sections = (
        ("long_reference", ("long_reshape", "reference")),
        ("restoration", ("restoration", "post_restore")),
    )
    for section_name, chart_names in sections:
        section = raw.get(section_name)
        require(isinstance(section, Mapping),
                "Serialized section is missing: " + section_name)
        charts = section.get("charts")
        require(isinstance(charts, Mapping),
                "Serialized charts are missing: " + section_name)
        for chart_name in chart_names:
            chart = charts.get(chart_name)
            require(isinstance(chart, Mapping),
                    "Serialized chart is missing: " + section_name + "." + chart_name)
            require(chart.get("ordinary_Z_orders") == [0, 1, 2, 3, 4],
                    "Serialized chart full rows changed: " + section_name + "." + chart_name)
            native = chart.get("native_caps")
            require(isinstance(native, Mapping)
                    and set(native) == {"phi", "U", "M", "A", "B", "C"},
                    "Serialized native cap names changed: " + section_name + "." + chart_name)
            for field, rows in native.items():
                _serialized_native_cap_rows(
                    rows, section_name + "." + chart_name + ".native_caps." + field
                )
            for field in ("phi_relative", "E_relative", "inverse_Utheta"):
                _serialized_packed_rows(
                    chart.get(field), 5,
                    section_name + "." + chart_name + "." + field
                )
            p2 = chart.get("p2")
            _serialized_packed_rows(
                p2, 4, section_name + "." + chart_name + ".p2"
            )
            root_caps = chart.get("root_caps")
            root_names = {"q_squared", "q", "u", "h", "r", "alpha", "K"}
            require(isinstance(root_caps, Mapping) and set(root_caps) >= root_names,
                    "Serialized root caps are incomplete: " + section_name + "." + chart_name)
            require("kappa" in root_caps,
                    "Serialized kappa cap is missing: " + section_name + "." + chart_name)
            for root_name, rows in root_caps.items():
                _serialized_packed_rows(
                    rows, 4,
                    section_name + "." + chart_name + ".root_caps." + root_name
                )
            _negative_packed_row(
                root_caps["kappa"][0],
                section_name + "." + chart_name + ".root_caps.kappa[0]",
            )
            conditioning = chart.get("conditioning")
            required_conditioning = {
                "log_D_lower", "log_lambda_lower", "positive_a_lower",
                "positive_gamma_lower", "positive_q_squared_lower_log",
                "Delta_negative",
            }
            require(isinstance(conditioning, Mapping)
                    and set(conditioning) >= required_conditioning
                    and conditioning.get("Delta_negative") is True,
                    "Serialized conditioning is incomplete: "
                    + section_name + "." + chart_name)
            for field in ("log_D_lower", "log_lambda_lower", "positive_a_lower",
                          "positive_gamma_lower", "positive_q_squared_lower_log"):
                _record(conditioning[field],
                        section_name + "." + chart_name + ".conditioning." + field)
            for field in ("positive_a_lower", "positive_gamma_lower"):
                _positive_record(
                    conditioning[field],
                    section_name + "." + chart_name + ".conditioning." + field,
                )
            if chart_name == "post_restore":
                for field in ("b", "t0"):
                    rows = root_caps[field]
                    require(all(
                        row.get("exact_zero") is True
                        and row.get("log_absolute_upper_bound") is None
                        for row in rows
                    ), "Serialized post-restore " + field + " rows are not exact zero")
                detail = chart.get("source_detail")
                require(isinstance(detail, Mapping)
                        and detail.get("physical_U_y") is True
                        and detail.get("b_formula") == "b=2*U_y/Utheta"
                        and detail.get("b_logs") == [None] * 5,
                        "Serialized post-restore U_y=0 receipt is missing")


def _validate_memory_rows(raw):
    long_reference = raw["long_reference"]["charts"]["reference"]
    restoration = raw["restoration"]["charts"]["restoration"]
    post_restore = raw["restoration"]["charts"]["post_restore"]
    expected = long_reference.get("native_caps")
    actual = restoration.get("incoming_memory")
    post_expected = restoration.get("native_caps")
    post_actual = post_restore.get("incoming_memory")
    require(isinstance(expected, Mapping) and isinstance(actual, Mapping)
            and isinstance(post_expected, Mapping) and isinstance(post_actual, Mapping),
            "Serialized memory rows are missing")
    require(all(expected.get(key) == actual.get(key) for key in ("A", "B", "C")),
            "Serialized restoration memory rows changed")
    require(all(post_expected.get(key) == post_actual.get(key)
                for key in ("A", "B", "C")),
            "Serialized post-restore memory rows changed")


def _lookup(value, path, label="value"):
    parts = path.split(".") if isinstance(path, str) else list(path)
    current = value
    for part in parts:
        require(isinstance(current, Mapping) and part in current,
                "Missing " + label + " path: " + ".".join(parts))
        current = current[part]
    return current


def _same_serialized_record(left, right, label):
    _record(left, label + ".actual")
    _record(right, label + ".saved")
    return (tuple(left["lower_exact_mpf_tuple"])
            == tuple(right["lower_exact_mpf_tuple"])
            and tuple(left["upper_exact_mpf_tuple"])
            == tuple(right["upper_exact_mpf_tuple"]))


def _parent_guard(parent, section, source_name):
    group = parent.get(section, {}).get("source_guards", {})
    if isinstance(group, Mapping):
        return group.get(source_name)
    if isinstance(group, list):
        matches = [row for row in group
                   if isinstance(row, Mapping) and row.get("function") == source_name]
        return matches[0] if len(matches) == 1 else None
    return None


def _validate_physical_geometry(parent, packet, guards):
    receipt = packet.get("physical_geometry")
    require(isinstance(receipt, Mapping)
            and receipt.get("passed") is True
            and receipt.get("unaffected_physical_geometry_reused") is True
            and receipt.get("saved_radii_exactly_equal") is True
            and receipt.get("unsafe_abs_native_cap_reused") is False
            and receipt.get("source_guard") == "physical_geometry"
            and receipt.get("source_executed") is False,
            "Physical geometry receipt is missing or changed")
    require(receipt.get("formula_source") == {
        "long_reshape": "logR=ln(110)+T",
        "reference": "logR=ln(110)+logref-8",
    }, "Physical geometry formulas changed")
    expected = (
        ("long_reshape", "long_reference.charts.long_reshape.scale.logR",
         "actual_long_reference_background_I1_numeric.charts.long_reshape.log_R_upper"),
        ("reference", "long_reference.charts.reference.scale.logR",
         "actual_long_reference_background_I1_numeric.charts.reference.log_R_upper"),
        ("restoration", "restoration.charts.restoration.scales.log_R_upper",
         "actual_restoration_native_controls.charts.restoration.log_R_upper"),
        ("post_restore", "restoration.charts.post_restore.scales.log_R_upper",
         "actual_restoration_native_controls.charts.post_restore.log_R_upper"),
    )
    rows = receipt.get("comparisons")
    require(isinstance(rows, list) and len(rows) == len(expected),
            "Physical geometry comparison receipt is incomplete")
    raw = packet.get("raw_interval_tree")
    for row, (name, provider_path, saved_path) in zip(rows, expected):
        require(isinstance(row, Mapping)
                and row.get("chart") == name
                and row.get("provider_path") == provider_path
                and row.get("saved_path") == saved_path
                and row.get("exact_endpoint_tuples_equal") is True,
                "Physical geometry comparison receipt changed: " + name)
        actual = _lookup(raw, provider_path, "physical geometry")
        saved = _lookup(parent, saved_path, "saved physical geometry")
        require(_same_serialized_record(actual, saved, "physical_geometry." + name),
                "Saved physical geometry changed: " + name)
    guard = guards.get("physical_geometry")
    require(isinstance(guard, Mapping) and guard.get("passed") is True
            and guard.get("source_executed") is False
            and guard.get("parent_guard_reused_as_data") is True,
            "Physical geometry source guard is missing")


def _validate_cut_numeric(packet):
    receipt = packet.get("cut_numeric_receipt")
    require(isinstance(receipt, Mapping)
            and receipt.get("passed") is True
            and receipt.get("formula") == "cut=1-Delta/eta"
            and receipt.get("numeric_lower_bound_above_one") is True
            and receipt.get("source_executed") is False,
            "Numeric cut receipt is missing or changed")
    rows = receipt.get("rows")
    require(isinstance(rows, list) and len(rows) == 4,
            "Numeric cut receipt does not cover four charts")
    expected = {
        ("long_reference", "long_reshape"),
        ("long_reference", "reference"),
        ("restoration", "restoration"),
        ("restoration", "post_restore"),
    }
    found = set()
    for row in rows:
        require(isinstance(row, Mapping)
                and row.get("formula") == "cut=1-Delta/eta"
                and row.get("derived_from") ==
                "kappa_upper<1, eta_upper<=1, Delta=-(2-kappa)"
                and row.get("Delta_negative") is True
                and row.get("source_executed") is False,
                "Numeric cut row changed")
        record = row.get("numeric_lower_bound")
        _record(record, "cut_numeric_receipt.numeric_lower_bound")
        require(_mpf_compare(record["lower_exact_mpf_tuple"], (0, 0, 0, 0)) > 0,
                "Numeric cut lower bound is not positive")
        require(_mpf_compare(record["lower_exact_mpf_tuple"], (0, 1, 0, 1)) > 0,
                "Numeric cut lower bound is not above one")
        eta_log, kappa_log = row.get("eta_source_log"), row.get("kappa_log_cap")
        _record(eta_log, "cut eta source log")
        _record(kappa_log, "cut kappa log cap")
        raw = packet["raw_interval_tree"]
        source_eta = raw["shared"]["logeta"]
        source_kappa = raw[row["section"]]["charts"][row["chart"]]["root_caps"]["kappa"][0]["log_absolute_upper_bound"]
        require(_same_serialized_record(eta_log, source_eta, "cut eta")
                and _same_serialized_record(kappa_log, source_kappa, "cut kappa")
                and _mpf_compare(eta_log["upper_exact_mpf_tuple"], (0, 0, 0, 0)) <= 0
                and _mpf_compare(kappa_log["upper_exact_mpf_tuple"], (0, 0, 0, 0)) < 0
                and row.get("stored_logarithms_exponentiated") is False
                and tuple(record["lower_exact_mpf_tuple"]) == (0, 1, 1, 1)
                and tuple(record["upper_exact_mpf_tuple"]) == (0, 1, 1, 1),
                "Algebraic cut lower bound or its logarithmic premises changed")
        found.add((row.get("section"), row.get("chart")))
    require(found == expected, "Numeric cut chart coverage changed")


def _validate_source_guards(root, parent, guards):
    specs = (
        ("full_p2", "actual_long_reference_native_controls", "full_p2"),
        ("functions_source", "actual_long_reference_native_controls", "actual_constructor"),
        ("bridge_fields", "actual_long_reference_native_controls", "generic_fields"),
        ("roots_source", "actual_long_reference_native_controls", "roots"),
        ("loop_source", "actual_long_reference_native_controls", "loop"),
        ("anchored_gradient", "actual_long_reference_native_controls", "anchored_gradient"),
        ("history_densities", "actual_long_reference_native_controls", "history_densities"),
        ("complete_histories", "actual_long_reference_native_controls", "complete_histories"),
        ("sigma_formula", "actual_long_reference_native_controls", "sigma_formula"),
        ("physical_geometry", "actual_long_reference_background_I1_numeric", "__init__"),
        ("restoration_functions", "actual_restoration_native_controls", "c6_chartmaps"),
        ("restoration_p2", "actual_restoration_native_controls", "q4_full_formula"),
    )
    for output, section, source_name in specs:
        guard = guards.get(output)
        old = _parent_guard(parent, section, source_name)
        require(isinstance(guard, Mapping) and isinstance(old, Mapping),
                "Source guard evidence is missing: " + output)
        for field in ("path", "function", "AST_sha256"):
            require(guard.get(field) == old.get(field),
                    "Source guard evidence changed: " + output + "." + field)
        require(guard.get("required_text") == list(old.get("required_text", [])),
                "Source guard snippets changed: " + output)
        require(guard.get("passed") is True
                and guard.get("source_executed") is False
                and guard.get("parent_guard_reused_as_data") is True,
                "Source guard execution receipt changed: " + output)
        source = root / _safe_relative(guard["path"])
        require(guard.get("source_sha256") == sha(source),
                "Source guard bytes changed: " + output)
    proof = guards.get("source_proof")
    require(isinstance(proof, Mapping)
            and proof.get("guarded_functions") == [output for output, _, _ in specs]
            and proof.get("passed") is True
            and proof.get("parent_guards_reused_as_data") is True
            and proof.get("source_executed") is False
            and proof.get("accepted_source_builders_called") is False,
            "Source proof receipt changed")


def _validate_packet(root, parent, packet):
    require(isinstance(packet, Mapping), "Directed-controls packet is missing")
    require(packet.get("schema") == SCHEMA and packet.get("passed") is True,
            "Directed-controls packet schema or pass flag changed")
    scope = packet.get("scope")
    require(isinstance(scope, Mapping), "Directed-controls scope is missing")
    required_true = (
        "full_corrected_native_controls_standalone",
        "corrected_native_controls_verified",
        "historical_abs_log_integrity_guard_preserved",
    )
    for key in required_true:
        require(scope.get(key) is True, "Required scope flag is missing: " + key)
    required_false = (
        "inherited_loop_revalidated", "current_N_revalidated",
        "full_pressure_revalidated", "global_five_moment_closure",
        "full_goal_complete", "inherited_log_caps_revalidated",
    )
    for key in required_false:
        require(scope.get(key) is False, "Scope was promoted: " + key)
    parent_scope = packet.get("parent_scope")
    require(isinstance(parent_scope, Mapping)
            and parent_scope.get("report") == PARENT_PATH
            and parent_scope.get("report_sha256") == PARENT_SHA256
            and parent_scope.get("report_bytes") == PARENT_BYTES
            and parent_scope.get("exact_parent_mapping_checked") is True
            and parent_scope.get("all_parent_pins_checked") is True
            and parent_scope.get("parent_pin_count") == (
                len(parent.get("input_hashes", {}))
                + len(parent.get("implementation_hashes", {}))
            ),
            "Parent scope receipt changed")
    require(packet.get("source_family") == parent.get("source_family"),
            "Source family changed")
    guards = packet.get("source_guards")
    require(isinstance(guards, Mapping), "Source guards are missing")
    for key in (
        "full_p2", "functions_source", "bridge_fields", "roots_source",
        "loop_source", "anchored_gradient", "history_densities",
        "complete_histories", "sigma_formula", "physical_geometry",
        "restoration_functions",
        "restoration_p2",
    ):
        guard = guards.get(key)
        require(isinstance(guard, Mapping) and guard.get("passed") is True
                and isinstance(guard.get("required_text"), list)
                and guard.get("source_executed") is False
                and guard.get("parent_guard_reused_as_data") is True,
                "Source AST guard is missing: " + key)
    proof = guards.get("source_proof")
    require(isinstance(proof, Mapping) and proof.get("passed") is True
            and proof.get("parent_guards_reused_as_data") is True
            and proof.get("source_executed") is False,
            "Source proof guard is missing")
    _validate_source_guards(root, parent, guards)
    premises = packet.get("source_premises")
    require(isinstance(premises, Mapping)
            and premises.get("passed") is True
            and premises.get("source_executed") is False,
            "Source premise receipt is missing")
    gradient_premise = premises.get("anchored_gradient", {})
    require(gradient_premise.get("passed") is True
            and gradient_premise.get("ordinary_not_taylor") is True
            and gradient_premise.get("row_count") == 4
            and gradient_premise.get("source_executed") is False,
            "Anchored-gradient source premise changed")
    sigma_premise = premises.get("sigma_derivative_bound", {})
    require(sigma_premise.get("passed") is True
            and sigma_premise.get("value_range") == [0, 1]
            and sigma_premise.get("derivative_range") == [0, 8]
            and sigma_premise.get("argument") ==
                "t=(1-2q)^2, x=8*sqrt(t)/(1-t)^2, sigma'=8*(1+3t)/(1-t)^3/cosh(x)^2"
            and sigma_premise.get("inequality") ==
                "cosh(x)^2>=1+x^2; 64-(6-3t+t^2)*(1-t)=58+9t-4t^2+t^3>0 for 0<=t<1"
            and sigma_premise.get("positive_polynomial") == "t**3 - 4*t**2 + 9*t + 58"
            and sigma_premise.get("flat_ends_by_definition") is True
            and sigma_premise.get("source_executed") is False,
            "Sigma derivative source premise changed")
    history_premise = premises.get("history_rate_source", {})
    expected_rates = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
    require(history_premise.get("passed") is True
            and history_premise.get("source") ==
                "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
            and history_premise.get("bridge_alias") == "RATES=bridge.RATES"
            and history_premise.get("rates") == expected_rates
            and history_premise.get("active_native_rates") == {"A": 1, "B": 2, "C": 1}
            and history_premise.get("source_executed") is False,
            "History-rate source premise changed")
    components = packet.get("components")
    require(isinstance(components, Mapping), "Component receipts are missing")
    common = components.get("common_inputs", {})
    long_component = components.get("long_reference", {})
    restoration_component = components.get("restoration", {})
    require(common.get("passed") is True and common.get("build_calls") == 1
            and common.get("fixture_calls") == 1
            and common.get("source_executed") is False,
            "Common component receipt is not single-call")
    require(long_component.get("passed") is True
            and long_component.get("compute_calls") == 1
            and long_component.get("source_executed") is False,
            "Long component receipt is not single-call")
    require(restoration_component.get("passed") is True
            and restoration_component.get("compute_calls") == 1
            and restoration_component.get("source_executed") is False,
            "Restoration component receipt is not single-call")
    fixture = packet.get("fixture")
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True
            and fixture.get("closed_log_fixtures") == 2
            and fixture.get("independent_pressure_sectors") == 6
            and fixture.get("nonzero_pressure_rows") == 24
            and fixture.get("historical_inputs_used") is False
            and fixture.get("accepted_self_checks_called") is False,
            "Common fixture receipt is missing or changed")
    packet_checks = packet.get("checks")
    require(isinstance(packet_checks, Mapping)
            and packet_checks.get("passed") is True
            and packet_checks.get("parent_bytes_checked") is True
            and packet_checks.get("parent_mapping_checked") is True
            and packet_checks.get("all_parent_pins_checked") is True
            and packet_checks.get("source_pins_checked") is True
            and packet_checks.get("source_guards_checked") is True
            and packet_checks.get("source_premises_checked") is True
            and packet_checks.get("context_precision_bits") == 1600
            and packet_checks.get("common_built_once") is True
            and packet_checks.get("long_computed_after_common") is True
            and packet_checks.get("restoration_computed_after_long") is True
            and packet_checks.get("raw_interval_tree_serialized") is True
            and packet_checks.get("full_five_row_shapes") is True
            and packet_checks.get("four_row_p2_shapes") is True
            and packet_checks.get("strict_slope_exact_tuple_check") is True
            and packet_checks.get("negative_b_margin_exact_tuple_check") is True
            and packet_checks.get("source_memory_input_checked") is True
            and packet_checks.get("physical_geometry_checked") is True
            and packet_checks.get("numeric_cut_checked") is True
            and packet_checks.get("old_errors_and_212_retained") is True
            and packet_checks.get("historical_directed_records") == 256
            and packet_checks.get("accepted_ancestors_replayed") is False
            and packet_checks.get("source_executed") is False,
            "Coordinator check receipt is missing or changed")
    _validate_execution_flags(packet)
    _validate_ordered(packet, "packet")
    raw = packet.get("raw_interval_tree")
    require(isinstance(raw, Mapping), "Serialized raw interval tree is missing")
    count = _record_count(raw, "raw_interval_tree")
    require(count > 0, "Serialized packet has no interval records")
    json.dumps(raw, allow_nan=False)
    _validate_native_shapes(raw)
    _validate_memory_rows(raw)
    _validate_physical_geometry(parent, packet, guards)
    _validate_cut_numeric(packet)
    shape = packet.get("shape_receipt")
    require(isinstance(shape, Mapping) and shape.get("passed") is True
            and shape.get("chart_full_rows") is True
            and shape.get("chart_p2_four_rows") is True
            and shape.get("native_cap_rows_raw_or_packed") is True
            and shape.get("root_names_complete") is True
            and shape.get("conditioning_checked") is True,
            "Full-five/four-row shape receipt is missing")
    memory = packet.get("memory_receipt")
    require(isinstance(memory, Mapping) and memory.get("passed") is True
            and memory.get("source_executed") is False
            and memory.get("numeric_actual_rows_equal") is True
            and memory.get("post_restore_rows_equal") is True
            and memory.get("fields") == ["A", "B", "C"],
            "Source memory input receipt is missing")
    preservation = packet.get("preservation")
    require(isinstance(preservation, Mapping) and preservation.get("passed") is True
            and preservation.get("parent_math_data_retained") is True,
            "Historical preservation receipt is missing")
    flags = preservation.get("flags", {})
    for key in ("old_212_terms_and_errors_unchanged",
                "six_signed_integrals_retained", "incoming_histories_not_reset"):
        require(flags.get(key) is True, "Historical preservation flag is missing: " + key)
    _strict_controls(packet)
    inputs, implementations = _merge_pins(root, parent, packet)
    return count, inputs, implementations


def _reviewed_check(root):
    for relative, expected in REVIEWED_SOURCE_HASHES.items():
        require(expected != "REVIEW_PENDING",
                "Independent source review is pending: " + relative)
        require(sha(root / relative) == expected,
                "Reviewed source bytes changed: " + relative)


def _historical_equal(parent, report):
    mutable = {
        "stage", "parent_report", "parent_report_sha256", "input_hashes",
        "implementation_hashes", "contract", "checks", "execution_seconds",
        "execution_audit", KEY,
    }
    require(set(report) == set(parent) | {KEY},
            "Report has unexpected or missing top-level fields")
    for key, value in parent.items():
        if key not in mutable:
            require(report.get(key) == value,
                    "Historical mathematical payload changed: " + key)


def _validate_report(root, parent, report):
    _historical_equal(parent, report)
    require(report.get("parent_report") == PARENT_NAME
            and report.get("parent_report_sha256") == PARENT_SHA256,
            "Report predecessor changed")
    packet = report.get(KEY)
    count, inputs, implementations = _validate_packet(root, parent, packet)
    require(report.get("input_hashes") == inputs,
            "Report input pin map is not the exact parent-plus-new map")
    require(report.get("implementation_hashes") == implementations,
            "Report implementation pin map is not the exact parent-plus-new map")
    checks = report.get("checks")
    require(isinstance(checks, Mapping) and checks.get("passed") is True
            and checks.get("scope") == "full_corrected_native_controls_standalone"
            and checks.get("inherited_caps_revalidated") is False
            and checks.get("directed_records") == 256
            and checks.get("inherited_loop_revalidated") is False
            and checks.get("current_N_revalidated") is False
            and checks.get("full_pressure_revalidated") is False
            and checks.get("global_five_moment_closure") is False
            and checks.get("full_goal_complete") is False,
            "Report checks promoted beyond corrected native controls")
    contract = report.get("contract")
    require(isinstance(contract, Mapping)
            and contract.get("corrected_native_controls_verified") is True
            and contract.get("inherited_log_cap_certification_verified") is False
            and contract.get("inherited_log_caps_revalidated") is False
            and contract.get("global_five_moment_closure") is False
            and contract.get("full_goal_complete") is False,
            "Report contract changed historical certification scope")
    return count, len(inputs) + len(implementations)


def produce(parent):
    _reviewed_check(ROOT)
    before = _check_parent(ROOT, parent)
    started = time.monotonic()
    sys.path.insert(0, str(ROOT / "src"))
    # Import the new coordinator only in production.  Inspection above and
    # below remains serialized-only and never imports this module.
    from openai_ns_reconstruction import original_native_directed_controls_v2 as coordinator

    packet = coordinator.build(ROOT, parent)
    report = copy.deepcopy(parent)
    report[KEY] = packet
    report["stage"] = "current_original_native_directed_controls_v2"
    report["parent_report"] = PARENT_NAME
    report["parent_report_sha256"] = PARENT_SHA256
    report["input_hashes"], report["implementation_hashes"] = _merge_pins(
        ROOT, parent, packet
    )
    contract = dict(report["contract"])
    contract.update(
        corrected_native_controls_verified=True,
        inherited_log_cap_certification_verified=False,
        inherited_log_caps_revalidated=False,
        full_pressure_revalidated=False,
        global_five_moment_closure=False,
        full_goal_complete=False,
    )
    report["contract"] = contract
    report["checks"] = {
        "passed": True,
        "scope": "full_corrected_native_controls_standalone",
        "corrected_native_controls_verified": True,
        "inherited_caps_revalidated": False,
        "directed_records": parent["checks"]["directed_records"],
        "inherited_loop_revalidated": False,
        "current_N_revalidated": False,
        "full_pressure_revalidated": False,
        "global_five_moment_closure": False,
        "full_goal_complete": False,
        "previous_packets_unaltered": True,
        "canonical_records": packet["checks"]["canonical_records"],
        "accepted_ancestors_replayed": False,
    }
    report["execution_audit"] = {
        "new_builder": "original_native_directed_controls_v2.build",
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_replayed": False,
        "accepted_self_checks_replayed": False,
        "old_stage_calls": False,
        "source_executed": False,
        "parent_payload_deepcopied": True,
        "parent_pins_checked": before,
    }
    report["execution_seconds"] = time.monotonic() - started
    _validate_report(ROOT, parent, report)
    return report


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-report", action="store_true")
    mode.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    parent = load(BASE / PARENT_NAME)
    _check_parent(ROOT, parent)
    if args.write_report:
        require(not REPORT_PATH.exists(),
                "Directed-controls report already exists; production replay refused")
        report = produce(parent)
        REPORT_PATH.write_text(
            json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    else:
        report = load(REPORT_PATH)
        records, pin_count = _validate_report(ROOT, parent, report)
        print(json.dumps({
            "passed": True,
            "report": str(REPORT_PATH),
            "bytes": REPORT_PATH.stat().st_size,
            "sha256": sha(REPORT_PATH),
            "canonical_records": records,
            "pins": pin_count,
            "scope": "full_corrected_native_controls_standalone",
            "inherited_loop_revalidated": False,
            "current_N_revalidated": False,
            "full_pressure_revalidated": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        }, ensure_ascii=False))


if __name__ == "__main__":
    main()
