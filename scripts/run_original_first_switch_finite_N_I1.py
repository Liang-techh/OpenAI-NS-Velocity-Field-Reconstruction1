"""Produce or inspect the first-switch finite-N I1 source packet.

The runner keeps the accepted macro finite-N report byte-for-byte as the
predecessor payload and appends the five reviewed first-switch packets.  The
new stages are imported lazily so that a static fixture can run while the
independent review pins are still pending.  Production is deliberately
review-gated and, once released, executes only the five new builders in their
declared order.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"

PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
PARENT_SHA256 = "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c"
NAME = "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"

BACKGROUND_KEY = "actual_first_switch_background_I1_numeric"
INPUT_KEY = "actual_first_switch_loop_Z3_inputs"
LOOP_KEY = "actual_first_switch_loop_Z3_bounds"
NUMERIC_KEY = "actual_first_switch_finite_N_I1_numeric"
PRESSURE_KEY = "actual_macro_first_switch_pressure_source"

BACKGROUND_MODULE = "src/openai_ns_reconstruction/original_first_switch_background_I1_numeric.py"
INPUT_MODULE = "src/openai_ns_reconstruction/original_first_switch_loop_Z3_inputs.py"
LOOP_MODULE = "src/openai_ns_reconstruction/original_first_switch_loop_Z3_bounds.py"
NUMERIC_MODULE = "src/openai_ns_reconstruction/original_first_switch_finite_N_I1_numeric.py"
# The source module carries the repository's established ``original_macro``
# prefix; the report key above names the actual first-switch stage.
PRESSURE_MODULE = "src/openai_ns_reconstruction/original_macro_first_switch_pressure_source.py"
RUNNER_MODULE = "scripts/run_original_first_switch_finite_N_I1.py"

CODE = (
    BACKGROUND_MODULE,
    INPUT_MODULE,
    LOOP_MODULE,
    NUMERIC_MODULE,
    PRESSURE_MODULE,
    RUNNER_MODULE,
)

# These are intentionally review gates rather than guessed hashes.  The root
# reviewer replaces the five values after independent mathematical review.
REVIEWED = {
    BACKGROUND_MODULE: "0c60e4bfdbd458c857363b119dcb61714a73d640230a3095c4c5f8583e74127d",
    INPUT_MODULE: "d372ccf26e226dac56aae49c6c3b2f4e383551c60c6927e68b1638f552f95f6d",
    LOOP_MODULE: "76d82c7661f52f5a9315d3aa4932747599298df8b537fa9d1107ee935b21ac88",
    NUMERIC_MODULE: "094dc8d966b77159953aab6552684c4c1ef498bdd319fd917c71df3a18b27baf",
    PRESSURE_MODULE: "24f5cbae0098ff9f0d7b4e146d0242f797f45f86343744e235e4f02d68115d0f",
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
    raw = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def parent():
    path = BASE / PARENT_NAME
    require(path.is_file(), "Accepted macro finite-N parent is missing")
    require(sha(path) == PARENT_SHA256, "Accepted macro finite-N parent bytes changed")
    result = load(path)
    require(result.get("checks", {}).get("passed") is True,
            "Accepted macro finite-N parent did not pass")
    return result


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]),
                          (ROOT, report["implementation_hashes"])):
        require(isinstance(group, dict), "Report source pin table is missing")
        for name, expected in group.items():
            require(isinstance(expected, str), "Source pin is not a SHA256: " + name)
            require(sha(folder / name) == expected,
                    "Source bytes changed: " + name)
            count += 1
    return count


def merge_pins(accepted, *packets):
    saved = dict(accepted["input_hashes"])
    implementations = dict(accepted["implementation_hashes"])
    saved[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        pins = packet.get("source_pins")
        require(isinstance(pins, dict) and pins,
                "Every first-switch packet must carry source_pins")
        for key, item in pins.items():
            path = Path(item.get("path", key) if isinstance(item, dict) else key)
            expected = item.get("expected_sha256") if isinstance(item, dict) else item
            require(isinstance(expected, str), "Missing source pin SHA256: " + str(path))
            if len(path.parts) == 1:
                actual = BASE / path
                name, target = path.name, saved
            elif path.as_posix().startswith("experiments/root_st073/"):
                actual = ROOT / path
                name, target = path.name, saved
            else:
                actual = ROOT / path
                name, target = path.as_posix(), implementations
            require(actual.is_file(), "Pinned source is missing: " + path.as_posix())
            require(sha(actual) == expected,
                    "New packet source bytes changed: " + path.as_posix())
            require(target.get(name, expected) == expected,
                    "Conflicting source pin: " + name)
            target[name] = expected
    for path in CODE:
        require((ROOT / path).is_file(), "New implementation is missing: " + path)
        implementations[path] = sha(ROOT / path)
    return saved, implementations


def check_records(value, path=""):
    """Check directed binary ordering without expanding large powers."""
    from mpmath.libmp import mpf_cmp

    count = 0
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lo, hi = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            require(
                all(
                    isinstance(item, list)
                    and len(item) == 4
                    and all(type(number) is int for number in item)
                    and item[0] in (0, 1)
                    and item[1] >= 0
                    and item[3] >= 0
                    for item in (lo, hi)
                ),
                "Malformed or nonfinite directed interval",
            )
            directed = value.get("directed") is True
            legacy_exact_source = (
                path.endswith(".R100_geometry_certificate.actual_T_singleton_tuple")
                and "directed" not in value and lo == hi
            )
            require(mpf_cmp(tuple(lo), tuple(hi)) <= 0 and (directed or legacy_exact_source),
                    "Reversed or nondirected interval: " + path)
            count += int(directed)
        count += sum(check_records(item, path + "." + str(key))
                     for key, item in value.items())
    elif isinstance(value, (list, tuple)):
        count += sum(check_records(item, path + "." + str(key))
                     for key, item in enumerate(value))
    return count


def _orders(rows, label, expected):
    require(isinstance(rows, list), label + " must be a list")
    require([row.get("ordinary_Z_order") for row in rows
             if isinstance(row, dict)] == list(expected),
            "Missing or unordered ordinary Z rows: " + label)
    return rows


def _record(row, field="value"):
    require(isinstance(row, dict), "Directed row is not a mapping")
    value = row.get(field)
    require(isinstance(value, dict), "Missing directed row field: " + field)
    require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
            "Incomplete directed row field: " + field)
    return value


def _compare_mpf(left, right):
    from mpmath.libmp import mpf_cmp
    return mpf_cmp(tuple(left), tuple(right))


def _require_parent_packet(packet, label):
    require(isinstance(packet, dict) and packet.get("passed") is True,
            label + " packet did not pass")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            label + " packet has the wrong parent SHA")
    require(isinstance(packet.get("source_pins"), dict)
            and packet["source_pins"],
            label + " source_pins are required")
    require(packet.get("physical_y_derivatives_claimed",
                       packet.get("limits", {}).get("physical_y_derivatives_claimed")) is False,
            label + " made a physical-y derivative claim")


def _validate_background(packet, accepted):
    _require_parent_packet(packet, "First-switch background")
    require(packet.get("ordinary_Z_orders") == [0, 1, 2, 3, 4],
            "First-switch background must retain ordinary Z rows 0..4")
    require(packet.get("actual_phase_Z_independent") is True
            and packet.get("finite_N_added") is False
            and packet.get("global_five_moment_closure") is False,
            "First-switch background scope changed")
    fourth = accepted["actual_macro_loop_Z3_inputs"]["native_Z4_extension"]["comparison_fifth_radius_proof"]
    frozen = _record(packet["comparison"]["frozen_phi"][0])
    source = fourth["comparison_phi"][0]
    floor = fourth["comparison_phi_floor"]
    require(frozen["lower_exact_mpf_tuple"] == source["lower_exact_mpf_tuple"]
            and _compare_mpf(floor["lower_exact_mpf_tuple"], [0, 0, 0, 0]) > 0
            and _compare_mpf(floor["lower_exact_mpf_tuple"],
                             source["lower_exact_mpf_tuple"]) <= 0,
            "Comparison reciprocal must use its own directed positive floor")
    geometry = packet.get("geometry_measure")
    require(isinstance(geometry, dict), "First-switch geometry receipt is missing")
    require(geometry.get("physical_measure") == "R^2*hbS*ds"
            and geometry.get("Wswitch_definition") == "10000*hbS"
            and packet["source_contract"].get("normalized_by_Wswitch") is True,
            "First-switch Wswitch source measure changed")
    for name in ("normalized_K", "loghbS", "logWswitch", "R_upper"):
        require(isinstance(geometry.get(name), dict),
                "Missing directed first-switch geometry record: " + name)
    proof = geometry.get("normalized_K_proof")
    require(isinstance(proof, dict) and all(proof.get(name) is True for name in (
        "log_width_plus_log4_lt_log_2minus1024",
        "two_hbS_le_half_2minus1024",
        "exp_half_width_lt_one_plus_2minus1024",
        "K_is_enclosure_not_selected_value",
        "hbS_not_materialized",
    )), "First-switch normalized measure proof is incomplete")
    fields = packet.get("native_fields")
    require(isinstance(fields, dict), "First-switch native fields are missing")
    for name in ("phi", "U", "M", "A", "B", "C"):
        _orders(fields.get(name), "background native " + name, range(5))
    i1 = packet.get("background_I1", packet.get("I1"))
    require(isinstance(i1, dict) and i1.get("ordinary_Z_orders") == [0, 1],
            "First-switch background ordinary C1 I1 packet is missing")
    _orders(i1.get("D_rows"), "background D_rows", range(2))
    contract = packet.get("source_contract", {})
    require(contract.get("first_switch_measure") == "R^2*hbS*ds"
            and contract.get("normalized_by_Wswitch") is True
            and contract.get("incoming_memory_retained") is True,
            "First-switch background source contract changed")
    return dict(passed=True, ordinary_C1=True, native_rows_0_to_4=True,
                own_Wswitch_enclosure=True, physical_y_claimed=False)


def _validate_inputs(packet):
    _require_parent_packet(packet, "First-switch loop input")
    require(packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
            "First-switch loop input must have four ordinary Z rows")
    proof = packet.get("proofs", {})
    require(isinstance(proof, dict) and all(proof.get(name) is True for name in (
        "kappa_less_than_one_cut_one",
        "current_first_switch_kappa_less_than_one_cut_one",
        "eta_positive_le_one",
        "eta_positive",
        "eta_at_most_one",
        "a_positive",
        "actual_N18",
        "actual_N18_is_graph_handle",
        "actual_N_not_materialized",
        "actual_phase_Z_independent",
        "real_first_switch_source_parameters",
        "first_switch_background_packet_passed",
        "F2_fifth_cap_proved",
    )) and proof.get("physical_y_derivatives_claimed") is False,
        "First-switch input premise certificate changed")
    kappa = packet.get("kappa_log_upper")
    require(isinstance(kappa, dict), "First-switch kappa log bound is missing")
    require(_compare_mpf(kappa["upper_exact_mpf_tuple"], [0, 0, 0, 0]) < 0,
            "First-switch kappa<1 certificate failed")
    conditioning = packet.get("conditioning")
    require(isinstance(conditioning, dict)
            and isinstance(conditioning.get("log_D_lower"), dict)
            and isinstance(conditioning.get("log_lambda_lower"), dict),
            "First-switch conditioning floors are missing")
    caps = packet.get("log_caps")
    require(isinstance(caps, dict), "First-switch logarithmic input caps are missing")
    for name in ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K"):
        _orders(caps.get(name), "first-switch log cap " + name, range(4))
    return dict(passed=True, kappa_cut_one=True, eta_in_unit_interval=True,
                actual_phase_Z_independent=True, four_log_cap_rows=True)


def _validate_loop(packet):
    _require_parent_packet(packet, "First-switch loop")
    require(packet.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and packet.get("universal_loop_controls_consumed") is False
            and packet.get("universal_B_over_E_value_cap_tightened") is False
            and packet.get("actual_phase_Z_independent") is True
            and packet.get("finite_N_phase_not_selected") is True,
            "First-switch loop scope changed")
    rows = packet.get("B_over_Pstar_over_frozen_E0_log_caps")
    _orders(rows, "B over Pstar over frozen E0", range(4))
    for row in rows:
        if row.get("exact_zero"):
            require(row.get("log_absolute_upper_bound") is None,
                    "Exact-zero B row carries a logarithmic record")
        else:
            _record(row, "log_absolute_upper_bound")
    proof = packet.get("inverse_proof", {})
    require(isinstance(proof, dict) and proof.get("passed") is True
            and proof.get("Z_only") is True
            and proof.get("physical_y_derivatives") is False
            and all(proof.get(name) is True for name in (
                "ordinary_Z_derivatives",
                "phase_terms_divided_by_Z_and_angle_factorials",
                "returned_Z_factorials_restored",
                "delta_zero_is_exact",
            )), "First-switch implicit inverse certificate changed")
    source = packet.get("source_formulas", {})
    require(source.get("B") == "E*M/2"
            and "compose_inverse" in source.get("T1hat", "")
            and "Phi[k,j>=1]" in source.get("Phi", ""),
            "First-switch phase source formulas are incomplete")
    return dict(passed=True, ordinary_Z_orders=[0, 1, 2, 3],
                inverse_Z_only=True, generic_C0_tightening_used=False)


def _validate_numeric(packet):
    _require_parent_packet(packet, "First-switch finite-N numeric")
    require(packet.get("accepted_ancestors_replayed") is False
            and packet.get("full_goal_complete") is False
            and packet.get("all_normalized_I1_errors_below_2_power_minus1024") is True,
            "First-switch finite-N numeric scope changed")
    contract = packet.get("source_contract", {})
    require(contract.get("passed") is True
            and contract.get("chart") == "first_switch"
            and contract.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and contract.get("N_handle") == 18
            and contract.get("physical_J_handle") == 9
            and contract.get("phase_Z_independent") is True
            and contract.get("physical_y_derivatives_claimed") is False
            and contract.get("nonreset_mean_positive_kernel") is True
            and contract.get("finite_N_radius_added_once_to_background") is True,
            "First-switch finite-N source contract changed")
    fields = packet.get("physical_error_log_bounds", {})
    for name, length in (("delta_U", 4), ("delta_M", 4), ("delta_Q", 3)):
        _orders(fields.get(name), "first-switch physical error " + name, range(length))
        for row in fields[name]:
            if not row.get("exact_zero"):
                _record(row, "log_absolute_upper_bound")
    i1 = packet.get("I1_full", {})
    d_full = packet.get("D_full", {})
    require(i1.get("passed") is True and d_full.get("passed") is True
            and i1.get("signed_background_unchanged") is True
            and i1.get("error_radius_is_not_a_selected_value") is True
            and d_full.get("computed_from_full_I1_once") is True
            and d_full.get("separate_D_error_not_added") is True,
            "First-switch finite-N I1/D enclosure changed")
    for rows in list(i1.get("six_integrals", {}).values()) + [d_full.get("D_rows")]:
        _orders(rows, "first-switch signed C1 rows", range(2))
    radius = packet.get("signed_I1_radius")
    require(isinstance(radius, dict)
            and radius.get("lower_exact_mpf_tuple") == [0, 1, -1024, 1]
            and radius.get("upper_exact_mpf_tuple") == [0, 1, -1024, 1],
            "Wrong first-switch finite-N I1 error radius")
    errors = packet.get("six_I1_difference_error_log_bounds", {})
    require(set(errors) == {"q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z"},
            "Incomplete first-switch six-I1 error certificate")
    threshold = packet["I1_error_threshold"]["lower_exact_mpf_tuple"]
    for name, rows in errors.items():
        _orders(rows, "first-switch error " + name, range(2))
        for row in rows:
            if row.get("exact_zero"):
                continue
            cap = _record(row, "log_absolute_upper_bound")
            require(_compare_mpf(cap["upper_exact_mpf_tuple"], threshold) < 0,
                    "First-switch finite-N error exceeds 2^-1024: " + name)
    limits = packet.get("limits", {})
    require(limits.get("first_switch_full_finite_N_six_I1_C1_enclosed") is True
            and limits.get("outer_RQy_retained") is True
            and limits.get("outer_RQy_evaluated") is False
            and limits.get("complete_pressure_global_closure") is False
            and limits.get("global_five_moment_closure") is False
            and limits.get("physical_y_derivatives_claimed") is False,
            "First-switch finite-N limits changed")
    joined = packet.get("joined_macro_first_switch_volume", {})
    require(joined.get("passed") is True
            and joined.get("complete_pressure_enclosed") is False
            and joined.get("outer_RQy_retained") is True,
            "Joined Wmacro first-switch enclosure was promoted")
    return dict(passed=True, ordinary_C1=True, finite_N_error_below_2minus1024=True,
                own_Wswitch_and_joined_Wmacro=True, pressure_outer_endpoint_retained=True,
                closure=False)


def _validate_pressure(packet):
    _require_parent_packet(packet, "First-switch pressure source")
    require(packet.get("accepted_ancestors_replayed") is False
            and packet.get("signed_relation_not_global_closure") is True,
            "First-switch pressure source scope changed")
    contract = packet.get("source_contract", {})
    require(contract.get("passed") is True
            and contract.get("chart") == "first_switch"
            and contract.get("ordinary_Z_orders") == [0, 1]
            and contract.get("phase_Z_independent") is True
            and contract.get("physical_y_derivatives_claimed") is False
            and contract.get("physical_Pstar_applied_once") is True
            and contract.get("local_S1_added_once") is True,
            "First-switch pressure source contract changed")
    limits = packet.get("limits", {})
    require(limits.get("outer_RQy_retained") is True
            and limits.get("outer_RQy_evaluated") is False
            and limits.get("complete_pressure_global_closure") is False
            and limits.get("global_closure") is False,
            "First-switch pressure endpoint or closure scope changed")
    guards = packet.get("source_guards")
    require(isinstance(guards, dict) and guards.get("passed") is True,
            "First-switch pressure source guards are missing")
    graph = packet.get("graph_patch")
    require(isinstance(graph, dict)
            and graph.get("base_prefix_unchanged") is True
            and graph.get("base_node_count") == 184703
            and graph.get("new_node_start") == 184703
            and graph.get("appended_node_count") == len(graph.get("appended_nodes", []))
            and graph.get("full_node_count") == 184703 + graph["appended_node_count"]
            and graph.get("base_prefix_canonical_nodes_sha256") ==
            "84f458ff131da40c3cdd6ed6dd0014fc661dad6352e2e24b01c9b46039293ac9",
            "First-switch pressure graph receipt is incomplete")
    joined = packet["joined_pressure_source"]
    require(joined["parent_joined"] == [184685, 184686]
            and joined["first_switch_local_S1"] == [148458, 148459]
            and joined["first_switch_upper_RQy"] == [148514, 148515]
            and joined["first_switch_lower_RQy"] == [184687, 184688]
            and joined["conditional_identity_only_in_joined_expression"] is True
            and contract["six_I1_D_unit"] == "raw physical graph; numeric packet Wswitch"
            and contract["numeric_D_is_not_full_S1"] is True,
            "Source join duplicates or mis-scales a pressure segment")
    return dict(passed=True, pressure_source_graph=True,
                outer_endpoint_retained=True, global_closure=False,
                source_graph_receipt="full prefix, local S1 once, retained external endpoints")


def validate_report(report, accepted):
    mutable = {
        "stage", "parent_report", "parent_report_sha256", "input_hashes",
        "implementation_hashes", "contract", "checks", "execution_seconds",
        "execution_audit",
    }
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Accepted payload changed: " + key)
    require(report.get("parent_report") == PARENT_NAME
            and report.get("parent_report_sha256") == PARENT_SHA256,
            "Wrong accepted predecessor")
    audit = report.get("execution_audit", {})
    require(audit.get("current_turn_accepted_ancestors_replayed") is False
            and audit.get("current_turn_accepted_producer_owner_norm_recurrence_quadrature_replayed") is False
            and audit.get("inherited_parent_execution_audit") == accepted["execution_audit"],
            "Current or inherited execution audit changed")
    packets = {
        BACKGROUND_KEY: report.get(BACKGROUND_KEY),
        INPUT_KEY: report.get(INPUT_KEY),
        LOOP_KEY: report.get(LOOP_KEY),
        NUMERIC_KEY: report.get(NUMERIC_KEY),
        PRESSURE_KEY: report.get(PRESSURE_KEY),
    }
    counts = {}
    for key, packet in packets.items():
        _require_parent_packet(packet, key)
        counts[key] = check_records(packet)
        require(counts[key] > 0, "Missing directed numeric records: " + key)
    return dict(
        background=_validate_background(packets[BACKGROUND_KEY], accepted),
        inputs=_validate_inputs(packets[INPUT_KEY]),
        loop=_validate_loop(packets[LOOP_KEY]),
        numeric=_validate_numeric(packets[NUMERIC_KEY]),
        pressure=_validate_pressure(packets[PRESSURE_KEY]),
        directed_records=counts,
        accepted_math_payload_preserved=True,
    )


def _new_builders():
    """Import only the five new first-switch builders, after review gating."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from openai_ns_reconstruction import (
        original_first_switch_background_I1_numeric as background,
        original_first_switch_finite_N_I1_numeric as numeric,
        original_first_switch_loop_Z3_bounds as loop,
        original_first_switch_loop_Z3_inputs as inputs,
        original_macro_first_switch_pressure_source as pressure,
    )
    return background, inputs, loop, numeric, pressure


def run():
    started = time.perf_counter()
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING",
                "Independent final mathematical review is pending: " + path)
        require(sha(ROOT / path) == expected,
                "Reviewed source bytes changed: " + path)
    accepted = parent()
    check_pins(accepted)
    background, inputs, loop, numeric, pressure = _new_builders()

    # This is the required production order.  Each stage consumes only the
    # packet produced immediately before it, while the accepted parent remains
    # a read-only predecessor.
    background_packet = background.build(ROOT, accepted)
    input_packet = inputs.build(ROOT, accepted, background_packet)
    loop_packet = loop.build(ROOT, accepted, input_packet)
    numeric_packet = numeric.build(ROOT, accepted, background_packet, loop_packet)
    pressure_packet = pressure.build(ROOT, accepted, numeric_packet)
    require(accepted == parent(), "New first-switch builders mutated accepted data")

    report = copy.deepcopy(accepted)
    saved, implementations = merge_pins(
        accepted, background_packet, input_packet, loop_packet,
        numeric_packet, pressure_packet,
    )
    report.update(
        stage="first_switch_finite_N_I1_C1_with_joined_Wmacro_enclosure",
        parent_report=PARENT_NAME,
        parent_report_sha256=PARENT_SHA256,
        input_hashes=saved,
        implementation_hashes=implementations,
    )
    report[BACKGROUND_KEY] = background_packet
    report[INPUT_KEY] = input_packet
    report[LOOP_KEY] = loop_packet
    report[NUMERIC_KEY] = numeric_packet
    report[PRESSURE_KEY] = pressure_packet
    report["execution_audit"] = dict(
        current_turn_accepted_ancestors_replayed=False,
        current_turn_accepted_producer_owner_norm_recurrence_quadrature_replayed=False,
        current_turn_new_source_builders_only=True,
        pre_acceptance_repairs=[
            "New D coefficient uses 2*h*Z^2, not 2*h^2*Z^2",
            "Width, smoothing, inlet and physical-unit source receipts checked",
            "First unaccepted integration exited 1 on a renamed coefficient key; new/parent spellings now supported",
            "Second unaccepted integration exited 1 on a legacy exact T singleton without a directed flag; validated separately, not counted as a directed enclosure",
        ],
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]),
    )
    report["contract"] = dict(
        first_switch_background_I1_C1_enclosed=True,
        first_switch_loop_Z3_enclosed=True,
        first_switch_finite_N_I1_C1_enclosed=True,
        first_switch_pressure_source_appended=True,
        own_Wswitch_enclosure=True,
        joined_Wmacro_enclosure=True,
        pressure_outer_endpoint_retained=True,
        complete_pressure_global_closure=False,
        terminal_five_moment_closure=False,
        physical_y_derivatives_claimed=False,
        full_goal_complete=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True,
        independent_runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=REVIEWED,
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    json.dumps(report, allow_nan=False)
    return report


def inspect(path):
    report, accepted = load(path), parent()
    require(report.get("checks", {}).get("passed") is True
            and report["checks"].get("independent_new_math_review_completed") is True
            and report["checks"].get("reviewed_implementation_sha256") == REVIEWED,
            "Saved first-switch report review binding changed")
    for path_name, expected in REVIEWED.items():
        require(report["implementation_hashes"].get(path_name) == expected
                == sha(ROOT / path_name),
                "Saved reviewed source changed: " + path_name)
    require(report["implementation_hashes"].get(RUNNER_MODULE) == sha(ROOT / RUNNER_MODULE),
            "Runner source pin changed")
    return dict(
        passed=True,
        pins=check_pins(report),
        packet_invariants=validate_report(report, accepted),
        serialized_only=True,
        accepted_builders_producers_or_checkers_replayed=False,
    )


def fixture_bound_check():
    """Pure static preflight; it never imports or executes a production stage."""
    require(CODE[-1] == RUNNER_MODULE, "Runner must be the final source pin")
    require(len(REVIEWED) == 5 and all(value == "REVIEW_PENDING" or
            (len(value) == 64 and all(ch in "0123456789abcdef" for ch in value))
            for value in REVIEWED.values()), "Five well-formed first-switch review gates required")
    require(PARENT_SHA256 ==
            "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c",
            "Parent source pin changed")
    return dict(
        passed=True,
        production_blocked_until_review=any(v == "REVIEW_PENDING" for v in REVIEWED.values()),
        production_order=[BACKGROUND_KEY, INPUT_KEY, LOOP_KEY, NUMERIC_KEY, PRESSURE_KEY],
        inherited_parent_payload=True,
        current_accepted_replay=False,
        own_Wswitch_and_joined_Wmacro_enclosure=True,
        finite_N_error_threshold="2^-1024",
        pressure_outer_endpoint_retained=True,
        complete_pressure_global_closure=False,
        source_graph_receipt="full prefix and local-source join checked during production",
    )


def self_check():
    return fixture_bound_check()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        result = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        temporary.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=result["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
