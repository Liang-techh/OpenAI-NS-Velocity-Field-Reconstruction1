"""Save or inspect the exact second-micro/frozen-macro Q4 R_Q_y seam.

The packet appends source graph rows for the adjacent endpoint identity.  It
preserves the accepted I1 report, incoming and cumulative S1 histories, and
all old signed/error packets.  It does not evaluate phase values, quadrature,
finite-N norms, or global closure.
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

from openai_ns_reconstruction import original_micro_macro_RQy_seam as current


build = current.build
PARENT_NAME = current.PARENT_NAME
PARENT = f"experiments/root_st073/{PARENT_NAME}"
PARENT_SHA256 = current.PARENT_SHA256
PINNED_HASHES = dict(current.PINNED_HASHES)
JOIN_CHECK_NAME = current.JOIN_CHECK_NAME
JOIN_SOURCE_NAME = current.JOIN_SOURCE_NAME
JOIN_ARCHIVE_NAME = current.JOIN_ARCHIVE_NAME
Q4_SOURCE_NAME = current.Q4_SOURCE_NAME
Q4_ARCHIVE_NAME = current.Q4_ARCHIVE_NAME
Q4_CHECK_NAME = current.Q4_CHECK_NAME

REVIEWED_MODULE_SHA256 = "06bd81597832f7308d2c9c5209b4d7dbb02ed9605c6bc12effc49ba7ccafd80f"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_macro_RQy_seam.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_macro_RQy_seam.py",
    "scripts/run_original_micro_macro_RQy_seam.py",
)

Q4_BASE_NODE_COUNT = 183700
FORBIDDEN_MACRO_ENDPOINTS = [141114, 141115, 141116, 141117]
MICRO_ENDPOINT = [133745, 133746, 133747, 133748]
MACRO_FINITE = [141027, 141028, 141029, 141030]
MICRO_QY = [133501, 133505, 133511, 133517]
MACRO_QY = [140868, 140872, 140878, 140884]
MICRO_QY_FULL = [133501, 133505, 133511, 133517, 133523]
MACRO_QY_FULL = [140868, 140872, 140878, 140884, 140890]
MICRO_INCOMING_S1 = [129638, 129639, 129640, 129641]
MICRO_CUMULATIVE_S1 = [133693, 133694, 133695, 133696]
MICRO_CUMULATIVE_S1_ENDPOINT = [133729, 133730, 133731, 133732]
FIRST_MICRO_ENDPOINT_S1 = [129638, 129639, 129640, 129641]
MACRO_INCOMING_S1 = [133729, 133730, 133731, 133732]
MACRO_CUMULATIVE_S1 = [141062, 141063, 141064, 141065]
MACRO_CUMULATIVE_S1_ENDPOINT = [141098, 141099, 141100, 141101]
GLOBAL_INITIAL_RQY = [129567, 129568, 129569, 129570]
MICRO_GLOBAL_M15_FINITE = [133721, 133722, 133723, 133724]
MICRO_GLOBAL_M15_ENDPOINT = [133737, 133738, 133739, 133740]
MACRO_GLOBAL_M15_FINITE = [141090, 141091, 141092, 141093]
MACRO_GLOBAL_M15_ENDPOINT = [141106, 141107, 141108, 141109]
Q4_APPENDED_NODE_COUNT = 102


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


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256,
            "Accepted I1 parent report changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True,
            "Checked I1 parent required")
    return report


def pins(report):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, expected in report[key].items():
            require(sha(base / name) == expected, "Changed source pin: " + name)
            count += 1
    return count


def _passed(value, label):
    require(isinstance(value, dict) and value.get("passed") is True, label)


def _int_handle(value, label, minimum=0):
    require(isinstance(value, int) and not isinstance(value, bool)
            and value >= minimum, label + " must be a graph handle")


def _handle_rows(value, label, count=4, minimum=Q4_BASE_NODE_COUNT):
    require(isinstance(value, list) and len(value) == count,
            label + " row count changed")
    for item in value:
        _int_handle(item, label, minimum)


def _exact_handles(value, expected, label):
    require(value == expected, label + " handles changed")


def _source_ast(receipt, label, path, functions):
    _passed(receipt, label + " source AST receipt missing")
    require(receipt.get("path") == Path(path).name,
            label + " source AST path changed")
    actual = receipt.get("functions")
    require(isinstance(actual, dict) and set(actual) == set(functions),
            label + " source AST function set changed")
    for name, snippets in functions.items():
        row = actual[name]
        require(isinstance(row, dict)
                and row.get("present") is True
                and row.get("matched") is True
                and row.get("tokens") == list(snippets),
                label + " source AST guard changed: " + name)


def _source_guards(packet):
    guards = packet.get("source_guards")
    require(isinstance(guards, dict), "Q4 seam source guards missing")
    require(set(guards) == {
        "passed", "accepted_join_checks", "join_source_ast", "q4_source_ast",
        "q4_rows", "join_archive_exact_Q4_prefix",
        "join_archive_full_Q_y_bindings",
        "forbidden_macro_endpoint_handles_consumed",
        "accepted_join_common_physical_source", "accepted_join_common_Q_y_rows",
        "accepted_join_full_Q_y_rows",
        "full_Q4_physical_Q_y_rows_bound", "physical_Q_y_Z0_Z3_rows",
        "accepted_join_common_U_mean_recovery",
        "accepted_join_physical_R_and_native_J_distinct",
        "accepted_smooth_source_premises", "accepted_Q4_Z_independent_binders",
        "radius_Z_independent_transitive",
    }, "Q4 seam source guard fields changed")
    require(guards["passed"] is True
            and guards["join_archive_exact_Q4_prefix"] is True
            and guards["join_archive_full_Q_y_bindings"] is True
            and guards["forbidden_macro_endpoint_handles_consumed"] is False
            and guards["accepted_join_common_physical_source"] is True
            and guards["accepted_join_full_Q_y_rows"] == {
                "second_micro": MICRO_QY_FULL,
                "frozen_macro": MACRO_QY_FULL,
            }
            and guards["full_Q4_physical_Q_y_rows_bound"] is True
            and guards["accepted_join_common_U_mean_recovery"] is True
            and guards["accepted_join_physical_R_and_native_J_distinct"] is True
            and guards["accepted_Q4_Z_independent_binders"] is True
            and guards["radius_Z_independent_transitive"] is True,
            "Q4 seam source guard controls changed")

    checks = guards["accepted_join_checks"]
    require(isinstance(checks, dict) and set(checks) == {
        "source_contract", "actual_background_recipes",
        "independent_physical_transport_seam",
        "actual_leading_recovery_and_common_N", "actual_macro_pressure_and_FTC",
        "actual_own_mixed_seam",
    }, "Accepted join check set changed")
    for key, value in checks.items():
        _passed(value, "Accepted join check missing: " + key)

    _source_ast(guards["join_source_ast"], "Join", JOIN_SOURCE_NAME, {
        "__init__": (
            "self.background=self.bridge['macro']",
            "self.radius=self.ref(self.background['radius'])",
            "self.phase=self.ref(self.loop['frozen_macro']['actual_phase'])",
            "self.N=self.base.N",
            "micro_hb_and_macro_length_native_Jacobians_remain_distinct=True",
        ),
        "_second_micro_packet": (
            "coordinate=self.ref(bg['phase'])",
            "radius=self.ref(bg['radius'])",
            "phase=self.ref(self.loop['second_micro']['actual_phase'])",
            "physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])]",
        ),
    })
    _source_ast(guards["q4_source_ast"], "Q4", Q4_SOURCE_NAME, {
        "physical_pressure_fields": (
            "U=a.scale(fields['V'],Pstar)",
            "Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar)",
            "my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        ),
        "__init__": (
            "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",
            "physical_measures='S0:R*Jacobian*dcoordinate; S1:R^2*Jacobian*dcoordinate'",
            "m15_minus_endpoint_sign_applied_once=True",
        ),
    })

    rows = guards["q4_rows"]
    require(rows == {
        "micro_endpoint_rows": MICRO_ENDPOINT,
        "macro_lower_source_rows": MACRO_FINITE,
        "cumulative_S1_rows_verified": True,
        "physical_Q_y_rows_verified": {
            "second_micro": MICRO_QY_FULL,
            "frozen_macro": MACRO_QY_FULL,
        },
        "radius_Z_independent_transitive": True,
        "forbidden_macro_endpoint_rows": FORBIDDEN_MACRO_ENDPOINTS,
    }, "Accepted Q4 source rows changed")
    require(guards["accepted_join_common_Q_y_rows"] == {
        "second_micro": MICRO_QY,
        "frozen_macro": MACRO_QY,
    }, "Accepted common physical Q_y rows changed")
    require(guards["physical_Q_y_Z0_Z3_rows"] == {
        "second_micro": MICRO_QY,
        "frozen_macro": MACRO_QY,
    }, "Full Q4 physical Q_y rows changed")
    smooth = guards["accepted_smooth_source_premises"]
    require(isinstance(smooth, dict) and set(smooth) == {
        "common_global_log_radius_phase_at_micro2_macro0_and_one_N",
        "original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows",
        "actual_micro_complete_correction_is_macro_incoming_so_mean_values_match",
        "original_mean_ODE_and_one_minus_delta_recovery_give_common_y4_Z4_Q",
    } and all(value is True for value in smooth.values()),
            "Accepted smooth source seam premises changed")


def _source_contract(packet, accepted):
    contract = packet.get("contract")
    require(isinstance(contract, dict), "Q4 seam contract missing")
    require(set(contract) == {
        "ordinary_Z_orders", "second_micro_upper_handles",
        "scope", "second_micro_lower_coordinate",
        "frozen_macro_lower_is_new_finite_boundary_substitution",
        "macro_upper_accepted_endpoint_handles_not_used",
        "signed_boundary_terms", "joined_pressure_expression",
        "global_initial_and_incoming_histories_retained",
        "first_micro_history_retained_as_incoming_source",
        "whole_first_to_macro_pair_not_claimed",
        "old_six_signed_terms_and_caps_retained",
        "cancellation_only_in_adjacent_joined_expression",
        "no_numeric_phase_or_quadrature", "no_outer_cap_deletion",
        "quantitative_mixed_norms_remain_open",
        "residual_constructed_without_cancellation_pair",
    }, "Q4 seam contract fields changed")
    require(
        contract["ordinary_Z_orders"] == [0, 1, 2, 3]
        and contract["second_micro_upper_handles"] == MICRO_ENDPOINT
        and contract["scope"] == "adjacent second_micro s=1..2 followed by frozen_macro m=0..1"
        and contract["second_micro_lower_coordinate"] == "s=1"
        and contract["frozen_macro_lower_is_new_finite_boundary_substitution"] is True
        and contract["macro_upper_accepted_endpoint_handles_not_used"] is True
        and contract["signed_boundary_terms"] == "-R_Q_y(second_micro,s=2)+R_Q_y(frozen_macro,m=0)"
        and contract["joined_pressure_expression"] == (
            "1/2*(S1_micro_at_s2-S1_micro_incoming)-R_Q_y(s=2)+R_Q_y(s=1) + "
            "1/2*(S1_macro_at_m1-S1_macro_incoming)-R_Q_y(m=1)+R_Q_y(m=0)"
        )
        and contract["global_initial_and_incoming_histories_retained"] is True
        and contract["first_micro_history_retained_as_incoming_source"] is True
        and contract["whole_first_to_macro_pair_not_claimed"] is True
        and contract["old_six_signed_terms_and_caps_retained"] is True
        and contract["cancellation_only_in_adjacent_joined_expression"] is True
        and contract["no_numeric_phase_or_quadrature"] is True
        and contract["no_outer_cap_deletion"] is True
        and contract["quantitative_mixed_norms_remain_open"] is True
        and contract["residual_constructed_without_cancellation_pair"] is True,
        "Q4 seam contract changed",
    )
    require(packet.get("source_family") == accepted["source_family"],
            "Q4 seam source family changed")


def _graph_patch(packet):
    patch = packet.get("graph_patch")
    require(isinstance(patch, dict) and set(patch) == {
        "base_archive", "base_node_count", "base_prefix_unchanged",
        "appended_nodes", "appended_node_count", "endpoint_handles",
    }, "Q4 seam graph patch fields changed")
    require(patch["base_archive"] == Q4_ARCHIVE_NAME
            and patch["base_node_count"] == Q4_BASE_NODE_COUNT
            and patch["base_prefix_unchanged"] is True,
            "Q4 graph prefix provenance changed")
    appended = patch["appended_nodes"]
    require(isinstance(appended, list)
            and isinstance(patch["appended_node_count"], int)
            and patch["appended_node_count"] == len(appended)
            and len(appended) == Q4_APPENDED_NODE_COUNT,
            "Q4 seam appended-node count changed")
    endpoints = patch["endpoint_handles"]
    require(isinstance(endpoints, dict) and set(endpoints) == {
        "second_micro_upper", "frozen_macro_lower", "signed_cancellation",
    }, "Q4 seam endpoint handle fields changed")
    _exact_handles(endpoints["second_micro_upper"], MICRO_ENDPOINT,
                   "Second-micro upper")
    for key in ("frozen_macro_lower", "signed_cancellation"):
        _handle_rows(endpoints[key], "Graph patch " + key)
    for row in appended:
        require(isinstance(row, dict) and isinstance(row.get("operation"), str),
                "Q4 appended graph row malformed")
        refs = []
        for key in ("arguments", "argument", "numerator", "denominator",
                    "expression", "variable", "value", "integrand",
                    "lower", "upper"):
            value = row.get(key)
            if isinstance(value, list):
                refs.extend(value)
            elif isinstance(value, int):
                refs.append(value)
        require(not any(value in FORBIDDEN_MACRO_ENDPOINTS for value in refs),
                "Forbidden accepted macro endpoint consumed")


def _proof_rows(packet):
    rows = packet.get("proof_rows")
    require(isinstance(rows, list) and len(rows) == 4,
            "Q4 seam proof rows must cover ordinary Z0..Z3")
    patch = packet["graph_patch"]
    base = patch["base_node_count"]
    appended = patch["appended_nodes"]

    def node(handle):
        _int_handle(handle, "Proof graph reference", base)
        offset = handle - base
        require(offset < len(appended), "Proof graph reference outside appended suffix")
        return appended[offset]

    def refs(row):
        result = []
        for key in ("arguments", "argument", "numerator", "denominator",
                    "expression", "variable", "value", "integrand",
                    "lower", "upper"):
            value = row.get(key)
            if isinstance(value, list):
                result.extend(value)
            elif isinstance(value, int):
                result.append(value)
        return result

    for order, row in enumerate(rows):
        require(isinstance(row, dict) and row.get("ordinary_Z_order") == order,
                "Q4 seam ordinary-Z proof order changed")
        require(set(row) == {
            "ordinary_Z_order", "second_micro", "frozen_macro",
            "signed_boundary_cancellation", "joined_adjacent_increment",
            "residual_after_boundary_cancellation",
            "reconstructed_joined_increment", "factorized_micro_endpoint",
            "factorized_macro_lower", "factorized_micro_endpoint_identity",
            "factorized_macro_lower_identity", "conditional_boundary_identity",
        }, "Q4 seam proof row fields changed")
        second = row["second_micro"]
        require(isinstance(second, dict) and set(second) == {
            "incoming_S1", "first_micro_history_endpoint",
            "cumulative_S1_finite_source",
            "cumulative_S1_at_s_2", "local_S1_difference",
            "R_Q_y_lower_at_s_1", "R_Q_y_upper_at_s_2",
            "global_initial_R_Q_y", "global_m15_finite_source",
            "global_m15_endpoint", "adjacent_increment",
        }, "Second-micro proof fields changed")
        require(second["incoming_S1"] == MICRO_INCOMING_S1[order],
                       "Second-micro incoming")
        require(second["first_micro_history_endpoint"]
                == FIRST_MICRO_ENDPOINT_S1[order],
                "First-micro incoming endpoint changed")
        require(second["cumulative_S1_finite_source"] == MICRO_CUMULATIVE_S1[order],
                       "Second-micro cumulative")
        require(second["cumulative_S1_at_s_2"] == MICRO_CUMULATIVE_S1_ENDPOINT[order],
                "Second-micro upper cumulative endpoint changed")
        _int_handle(second["local_S1_difference"], "Second-micro local S1", Q4_BASE_NODE_COUNT)
        _int_handle(second["R_Q_y_lower_at_s_1"], "Second-micro lower", Q4_BASE_NODE_COUNT)
        require(second["R_Q_y_upper_at_s_2"] == MICRO_ENDPOINT[order],
                "Second-micro upper endpoint changed")
        require(second["global_initial_R_Q_y"] == GLOBAL_INITIAL_RQY[order]
                and second["global_m15_finite_source"] == MICRO_GLOBAL_M15_FINITE[order]
                and second["global_m15_endpoint"] == MICRO_GLOBAL_M15_ENDPOINT[order],
                "Second-micro retained global source changed")
        _int_handle(second["adjacent_increment"], "Second-micro increment", Q4_BASE_NODE_COUNT)
        local_node = node(second["local_S1_difference"])
        require(local_node.get("operation") == "sum"
                and local_node.get("arguments", [None, None])[0]
                == second["cumulative_S1_at_s_2"],
                "Second-micro local S1 does not use the s=2 endpoint")
        local_args = local_node["arguments"]
        require(len(local_args) == 2
                and node(local_args[1]).get("operation") == "negative"
                and node(local_args[1]).get("argument") == second["incoming_S1"],
                "Second-micro incoming S1 subtraction changed")

        macro = row["frozen_macro"]
        require(isinstance(macro, dict) and set(macro) == {
            "incoming_S1", "cumulative_S1_finite_source",
            "cumulative_S1_at_m_1", "local_S1_difference",
            "R_Q_y_lower_at_m_0", "R_Q_y_upper_derived_at_m_1",
            "global_m15_finite_source", "global_m15_endpoint",
            "adjacent_increment",
        }, "Frozen-macro proof fields changed")
        require(macro["incoming_S1"] == MACRO_INCOMING_S1[order]
                and macro["cumulative_S1_finite_source"] == MACRO_CUMULATIVE_S1[order]
                and macro["cumulative_S1_at_m_1"] == MACRO_CUMULATIVE_S1_ENDPOINT[order]
                and macro["global_m15_finite_source"] == MACRO_GLOBAL_M15_FINITE[order]
                and macro["global_m15_endpoint"] == MACRO_GLOBAL_M15_ENDPOINT[order],
                "Frozen-macro retained source changed")
        for key in ("local_S1_difference", "R_Q_y_lower_at_m_0",
                    "R_Q_y_upper_derived_at_m_1", "adjacent_increment"):
            _int_handle(macro[key], "Frozen-macro " + key, Q4_BASE_NODE_COUNT)
        macro_local_node = node(macro["local_S1_difference"])
        require(macro_local_node.get("operation") == "sum"
                and macro_local_node.get("arguments", [None, None])[0]
                == macro["cumulative_S1_at_m_1"],
                "Frozen-macro local S1 does not use the m=1 endpoint")
        macro_local_args = macro_local_node["arguments"]
        require(len(macro_local_args) == 2
                and node(macro_local_args[1]).get("operation") == "negative"
                and node(macro_local_args[1]).get("argument") == macro["incoming_S1"],
                "Frozen-macro incoming S1 subtraction changed")
        for key in ("signed_boundary_cancellation", "joined_adjacent_increment",
                    "residual_after_boundary_cancellation",
                    "reconstructed_joined_increment", "factorized_micro_endpoint",
                    "factorized_macro_lower"):
            _int_handle(row[key], "Joined proof " + key, Q4_BASE_NODE_COUNT)
        residual_node = node(row["residual_after_boundary_cancellation"])
        require(row["signed_boundary_cancellation"] not in refs(residual_node),
                "Residual graph references the cancellation pair")
        for key, finite, wrapper, qy, radius, variable, value in (
            ("factorized_micro_endpoint_identity", [133660, 133661, 133662, 133663][order],
             MICRO_ENDPOINT[order], MICRO_QY[order], 68644, 66988, 48),
            ("factorized_macro_lower_identity", MACRO_FINITE[order],
             row["frozen_macro"]["R_Q_y_lower_at_m_0"], MACRO_QY[order],
             71448, 70265, 0),
        ):
            identity = row[key]
            require(isinstance(identity, dict) and set(identity) == {
                "passed", "finite_product", "actual_wrapper",
                "factorized_wrapper_product", "substituted_Q_y",
                "substituted_radius", "variable", "value",
                "Z_independent_substitution", "identity",
            } and identity["passed"] is True
                    and identity["finite_product"] == finite
                    and identity["actual_wrapper"] == wrapper
                    and identity["factorized_wrapper_product"] >= Q4_BASE_NODE_COUNT
                    and identity["substituted_Q_y"] >= Q4_BASE_NODE_COUNT
                    and identity["substituted_radius"] >= Q4_BASE_NODE_COUNT
                    and identity["variable"] == variable
                    and identity["value"] == value
                    and identity["Z_independent_substitution"] is True
                    and identity["identity"] ==
                    "sub(product(Q_y,R),x,a)=product(sub(Q_y,x,a),sub(R,x,a))",
                    "Factorized source-wrapper identity changed")
        boundary = row["conditional_boundary_identity"]
        require(isinstance(boundary, dict) and set(boundary) == {
            "passed", "ordinary_Z_order", "common_Q_y_source_rows",
            "common_physical_radius_cutpoint", "one_common_N_and_global_phase",
            "mean_ODE_and_complete_correction_premises", "identity",
        } and boundary["passed"] is True
                and boundary["ordinary_Z_order"] == order
                and boundary["common_Q_y_source_rows"] is True
                and boundary["common_physical_radius_cutpoint"] is True
                and boundary["one_common_N_and_global_phase"] is True
                and boundary["mean_ODE_and_complete_correction_premises"] is True
                and boundary["identity"] ==
                "R_micro(2)*Q_y_micro(2)=R_macro(0)*Q_y_macro(0)",
                "Conditional boundary identity changed")

    endpoints = patch["endpoint_handles"]
    require(
        endpoints["frozen_macro_lower"]
        == [row["frozen_macro"]["R_Q_y_lower_at_m_0"] for row in rows]
        and endpoints["signed_cancellation"]
        == [row["signed_boundary_cancellation"] for row in rows],
        "Q4 seam endpoint handle rows do not match proof rows",
    )


def _symbolic_controls(packet):
    controls = packet.get("symbolic_controls")
    _passed(controls, "Q4 seam symbolic controls missing")
    require(set(controls) == {
        "passed", "radius_cutpoint", "radius_identity", "pressure_increment",
        "signed_FTC_identity", "conditional_boundary_pair_zero",
        "Z_independent_substitution_product_identity", "Q_y_operator",
        "ordinary_Z_orders", "Q_y_common_source_identity",
        "source_bound_full_Q_y_arrays", "source_bound_smooth_seam_premises",
    }, "Q4 seam symbolic control fields changed")
    require(
        controls["radius_identity"] is True
        and controls["signed_FTC_identity"] is True
        and controls["conditional_boundary_pair_zero"] is True
        and controls["Z_independent_substitution_product_identity"] is True
        and controls["ordinary_Z_orders"] == [0, 1, 2, 3]
        and controls["Q_y_common_source_identity"] is True
        and controls["source_bound_full_Q_y_arrays"] is True
        and controls["source_bound_smooth_seam_premises"] is True
        and controls["radius_cutpoint"] == (
            "R_micro(s)=R_a*exp(hb*s), R_macro(m)=R_a*exp(2*hb+length*m); s=2,m=0"
        )
        and controls["pressure_increment"] == (
            "Delta m15 = (S1_cumulative-S1_incoming)/2 - R_Q_y(upper) + R_Q_y(lower)"
        )
        and controls["Q_y_operator"] == (
            "[2 Z U_y - (1-delta) Z m_y - (1-Z^2) (m_y)_Z]/(1-delta Z^2)"
        ),
        "Q4 seam symbolic controls changed",
    )


def _packet_fields(report, accepted):
    require(report.get("parent_report") == Path(PARENT).name
            and report.get("parent_report_sha256") == PARENT_SHA256,
            "Same I1 predecessor required")
    mutable = {
        "stage", "parent_report", "parent_report_sha256", "input_hashes",
        "implementation_hashes", "contract", "checks", "execution_seconds",
    }
    for key, value in accepted.items():
        if key not in mutable:
            require(key in report and report[key] == value,
                    "Inherited I1 packet or error changed: " + key)
    packet = report.get("actual_micro_macro_RQy_seam")
    _passed(packet, "Q4 seam packet missing")
    require(set(packet) == {
        "passed", "source_family", "source_pins", "input_pins", "source_guards",
        "proof_rows", "graph_patch", "symbolic_controls", "contract",
    }, "Q4 seam packet fields changed")
    _source_contract(packet, accepted)
    _source_guards(packet)
    _graph_patch(packet)
    _proof_rows(packet)
    _symbolic_controls(packet)

    expected_pins = {PARENT_NAME: PARENT_SHA256, **PINNED_HASHES}
    source_pins = packet.get("source_pins")
    input_pins = packet.get("input_pins")
    require(source_pins == expected_pins and input_pins == expected_pins,
            "Q4 seam source/input pin set changed")
    for name, expected in expected_pins.items():
        require(sha(DIRECTORY / name) == expected,
                "Changed Q4 seam input: " + name)
        require(name == PARENT_NAME
                or accepted["input_hashes"].get(name) == expected,
                "Q4 seam input is not inherited from the accepted parent: " + name)
    return dict(
        passed=True,
        exact_second_micro_upper_to_frozen_macro_lower_RQy_identity=True,
        ordinary_Z0_Z3_endpoint_and_history_rows=True,
        accepted_Q4_physical_Q_y_rows_bound=True,
        distinct_native_phase_history_and_physical_y_sources=True,
        signed_boundary_cancellation_only_local_join=True,
        first_micro_history_retained_as_incoming_source=True,
        whole_first_to_macro_pair_not_claimed=True,
        old_212_six_signed_and_error_packets_retained=True,
        quantitative_norms_and_global_closure_open=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed seam module may produce this report")
    accepted = parent()
    pins(accepted)
    packet = build(ROOT, accepted)
    require(accepted == parent(), "New seam mutated its accepted parent")
    _packet_fields({
        **accepted,
        "parent_report": Path(PARENT).name,
        "parent_report_sha256": PARENT_SHA256,
        "actual_micro_macro_RQy_seam": packet,
    }, accepted)

    inputs = dict(accepted["input_hashes"])
    implementation = dict(accepted["implementation_hashes"])
    for name, expected in packet["input_pins"].items():
        require(sha(DIRECTORY / name) == expected,
                "Changed new input: " + name)
        if name != PARENT_NAME:
            require(inputs.get(name) == expected,
                    "Conflicting inherited source pin: " + name)
        inputs[name] = expected
    implementation.update({path: sha(ROOT / path) for path in CODE})
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_micro_macro_RQy_seam",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_macro_RQy_seam=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        exact_second_micro_upper_to_frozen_macro_lower_RQy_identity=True,
        ordinary_Z0_to_Z3_endpoint_and_history_rows=True,
        accepted_Q4_physical_Q_y_rows_bound=True,
        distinct_native_phase_history_and_physical_y_sources=True,
        signed_boundary_cancellation_only_local_join=True,
        first_micro_history_retained_as_incoming_source=True,
        whole_first_to_macro_pair_not_claimed=True,
        old_212_six_signed_and_error_packets_retained=True,
        quantitative_norms_and_global_closure_open=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Exact source-linked Q4 second-micro upper/frozen-macro lower R_Q_y "
            "endpoint identity, ordinary Z0..Z3 graph rows, retained incoming and "
            "cumulative S1 histories, and local signed boundary cancellation; "
            "quantitative finite-N norms, outer caps, and global closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independent_runtime_worker_model_metadata_available=False,
        accepted_ancestors_replayed=False,
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report = _load(path)
    accepted = parent()
    checks = report["checks"]
    require(
        checks["passed"]
        and checks.get("independent_new_math_review_completed") is True
        and checks.get("independent_review_completed") is True
        and checks["independently_reviewed_module_sha256"]
        == report["implementation_hashes"][CODE[0]]
        == REVIEWED_MODULE_SHA256
        and report["implementation_hashes"][CODE[1]] == sha(ROOT / CODE[1]),
        "Same independently reviewed seam implementation required",
    )
    return dict(
        passed=True,
        pins=pins(report),
        packet_invariants=_packet_fields(report, accepted),
        serialized_only=True,
        accepted_owner_norm_recurrence_or_quadrature_replayed=False,
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
