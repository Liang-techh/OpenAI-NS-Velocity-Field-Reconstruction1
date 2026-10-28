"""UNEXECUTED M-channel formula-level source-linkage consumer draft.

This is a separate loader for the remaining M port.  It reads only saved JSON
publications and the saved result JSON embedded in the accepted D report.  It
does not import a reconstruction module, call a source expression, decode
gzip, compare source graphs, evaluate a numerical value, or enter the failed
generic M/Z0 atom family.

The loader uses the source-owned symbolic definition

    M_C = T_C_m / A_C
    M_N = T_N_m / A_N
    delta_m = T_N_m - T_C_m
    a = A_N - A_C
    M_gap = (A_C * delta_m - a * T_C_m) / (A_C * A_N)

It binds the four saved D ``m`` rows and validates their source-role
namespace, field, order, and source-node handles against the joint role
vectors.  The saved base source recipe connects the C3/C5 M numerator,
amplitude, and typed endpoint only at ordinary Z0.  Z1--Z3 endpoint positions
are listed as reads, not claimed as bindings; their true derivative relations
remain open.  The Z0 result is only a saved positive base-recipe handle
relation, not a cross-archive source-cut identity.
It deliberately does not access the joint report's literal M numerator,
denominator, or M atom-term subtree.  The scoped C3/native amplitude certificate
is loaded as a premise, but its exact binding to both joint C3/C5 amplitude
vectors is not established here; therefore ``a`` is retained symbolically and
is never replaced by zero.

Endpoint equality to ``T_m / A`` remains a conditional open obligation, as do
common physical points, units, finite-N error transport, the M-to-m15 map, and
the independent reference offset.  There is intentionally no ``__main__``
entry point.  A future root-owned review may invoke
``load_M_formula_source_linkage``; this draft has not been invoked.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path


SCHEMA = "compact_M_formula_level_source_linkage_loader_draft_v4"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
FORMULA = "M_gap = (A_C * delta_m - a * T_C_m) / (A_C * A_N)"

D_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_compact_D_km_source_binding.json"
)
# This is the raw SHA of the publication report itself.
D_REPORT_SHA256 = (
    "9af596c029d7eadd78b25a7895ffc8f4f7f24aefcea0a3b6f4f07f6c28669587"
)
# This is the SHA of the embedded D result.json bytes, not the report SHA.
D_EMBEDDED_RESULT_SHA256 = (
    "0ebea65173f200afcee2a917804d3674861c66b24bd6120c17612c993a845e24"
)

JOINT_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "joint_target_numerator_transport.json"
)
JOINT_REPORT_SHA256 = (
    "cff727e06c1af18eebeba0bbf027a488128c2021f9feb839bafdfc368b18f133"
)

AMPLITUDE_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "C3_outer_amplitude_identity.json"
)
AMPLITUDE_REPORT_SHA256 = (
    "987e1aba8d33a22b325e1f0386a00e522778d296b65f98ac7d826da228292975"
)

RAW_ARCHIVE_SHA256 = {
    "C3": "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "C5": "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
}

PRIOR_GENERIC_M_Z0_FAILURE = {
    "report_path": (
        "docs/CURRENT_COMPACT_JOINT_HISTORY_GAP_BINDING_FAILURE_2026_10_28.md"
    ),
    "result_sha256": "50257d1c62728ba554228f69e68d50d4c4be673ff0b0bf13e79b3c6b75b3be3c",
    "exception": "KeyError: 0",
    "reentered": False,
    "proof_body_reused": False,
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha256_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_authenticated_json(root, relative, expected_sha256):
    """Read saved JSON after a byte hash; never imports or evaluates source."""
    path = Path(root) / relative
    _require(path.is_file(), f"missing saved JSON: {relative}")
    _require(
        _sha256_file(path) == expected_sha256,
        f"saved JSON SHA changed: {relative}",
    )
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"saved JSON is not an object: {relative}")
    return value


def _read_embedded_D_result(report):
    """Decode only the saved result JSON, not source or gzip archive bytes."""
    checkpoint = report["raw_checkpoints"]["result.json"]
    raw = base64.b64decode(checkpoint["base64"], validate=True)
    _require(_sha256_bytes(raw) == D_EMBEDDED_RESULT_SHA256,
             "embedded D result SHA changed")
    _require(checkpoint["sha256"] == D_EMBEDDED_RESULT_SHA256,
             "D result checkpoint digest disagrees")
    value = json.loads(raw.decode("utf-8"))
    _require(isinstance(value, dict), "embedded D result is not an object")
    _require(value.get("result") is not None, "embedded D result has no result")
    return value["result"]


def _assert_publication_state(D_report, D_result, joint_report, amplitude_report):
    """Check saved publication flags without accessing the M atom subtree."""
    _require(
        D_report["schema"] == "compact_D_km_source_binding_publication_v1",
        "unexpected D publication schema",
    )
    _require(D_report["source_only"] is True, "D publication is not source-only")
    _require(D_report["numeric_values_computed"] is False,
             "D publication contains numeric values")
    _require(D_report["saved_generic_failure_reentered"] is False,
             "D publication reentered the failed generic family")
    _require(D_result["schema"] == "compact_D_km_reference_offset_source_binding_draft_v1",
             "unexpected embedded D result schema")
    _require(D_result["source_only"] is True, "embedded D result is not source-only")
    _require(D_result["signed_D_history_row_count"] == 8,
             "embedded D result does not retain eight signed rows")

    _require(joint_report["schema"] ==
             "original_compact_joint_target_numerator_transport_publication_v1",
             "unexpected joint publication schema")
    _require(joint_report["source_only"] is True,
             "joint publication is not source-only")
    _require(joint_report["source_packet_definition_accepted"] is True,
             "joint source-packet definition is not accepted")
    _require(joint_report["source_packet_defined"] is True,
             "joint source packet is not defined")
    _require(joint_report["joint_target_numerator_source_defined"] is True,
             "joint target numerator source is not defined")
    _require(joint_report["five_joint_numerator_identities_proved"] is True,
             "five joint numerator identities are not accepted")
    _require(joint_report["ordinary_Z_quotient_transport_proved"] is True,
             "joint quotient transport is not accepted")
    _require(joint_report["signed_target_gaps_evaluated"] is False,
             "joint publication evaluated signed gaps")
    _require(joint_report["numeric_values_computed"] is False,
             "joint publication contains numeric values")

    _require(
        amplitude_report["schema"] ==
        "original_compact_C3_outer_amplitude_identity_publication_v1",
        "unexpected amplitude publication schema",
    )
    _require(amplitude_report["source_only"] is True,
             "amplitude publication is not source-only")
    _require(amplitude_report["original_source_evaluations"] == 0,
             "amplitude publication evaluated an original source")
    _require(amplitude_report["numeric_values_computed"] is False,
             "amplitude publication contains numeric values")
    _require(amplitude_report["C3_outer_amplitude_source_semantics_bound"] is True,
             "C3 amplitude source semantics are not bound")
    _require(amplitude_report["actual_native_C3_amplitude_identity_proved"] is True,
             "scoped native/C3 amplitude identity is not present")
    _require(amplitude_report["native_outer_to_C5_amplitude_source_alias_proved"] is True,
             "native-to-C5 amplitude alias premise is not present")


def _handle_tuple(value, label):
    """Parse saved role handles only; no source expression is interpreted."""
    if isinstance(value, str):
        tokens = value.split()
    else:
        tokens = list(value)
    try:
        handles = tuple(int(token) for token in tokens)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not a saved handle vector") from exc
    _require(len(handles) == len(ORDINARY_Z_ORDERS),
             f"{label} does not have four ordinary-Z handles")
    return handles


def _read_typed_roles(joint_report):
    """Read typed endpoint/amplitude/history roles, never literal M terms."""
    roles = joint_report["source_recipe_contract"]["ordinary_Z_roles"]
    _require(tuple(roles["orders"]) == ORDINARY_Z_ORDERS,
             "ordinary-Z role order changed")
    _require(roles["orders_are_source_roles"] is True,
             "ordinary-Z order vector is not marked as source roles")
    _require(roles["range_values_are_not_function_values"] is True,
             "ordinary-Z role ranges are marked as function values")
    C3 = roles["C3"]
    C5 = roles["C5"]
    _require(C3["namespace"] == "repair_C3_band",
             "C3 source-role namespace changed")
    _require(C5["namespace"] == "true_compact_C5",
             "C5 source-role namespace changed")
    return {
        "C3": {
            "manifest_namespace": "C3",
            "role_namespace": C3["namespace"],
            "M_endpoint_handles": _handle_tuple(
                C3["target_ports_Z0_to_Z3"]["M"], "C3 M typed endpoints"
            ),
            "M_endpoint_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C3."
                "target_ports_Z0_to_Z3.M"
            ),
            "m_history_handles": _handle_tuple(
                C3["terminal_history_ports_Z0_to_Z3"]["m"],
                "C3 m history roles",
            ),
            "m_history_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C3."
                "terminal_history_ports_Z0_to_Z3.m"
            ),
            "m_history_source_role": (
                "source_recipe_contract.ordinary_Z_roles.C3."
                "terminal_history_ports_Z0_to_Z3.m"
            ),
            "amplitude_handles": _handle_tuple(
                C3["amplitude_ports_Z0_to_Z3"], "C3 amplitude roles"
            ),
            "amplitude_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C3."
                "amplitude_ports_Z0_to_Z3"
            ),
        },
        "C5": {
            "manifest_namespace": "C5",
            "role_namespace": C5["namespace"],
            "M_endpoint_handles": _handle_tuple(
                C5["target_ports_Z0_to_Z3"]["M"], "C5 M typed endpoints"
            ),
            "M_endpoint_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C5."
                "target_ports_Z0_to_Z3.M"
            ),
            "m_history_handles": _handle_tuple(
                C5["terminal_scaled_history_ports_Z0_to_Z3"]["m"],
                "C5 scaled m history roles",
            ),
            "m_history_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C5."
                "terminal_scaled_history_ports_Z0_to_Z3.m"
            ),
            "m_history_source_role": (
                "source_recipe_contract.ordinary_Z_roles.C5."
                "terminal_scaled_history_ports_Z0_to_Z3.m"
            ),
            "raw_correction_handles": _handle_tuple(
                C5["raw_correction_ports_Z0_to_Z3"]["m"],
                "C5 raw correction roles",
            ),
            "amplitude_handles": _handle_tuple(
                C5["amplitude_ports_Z0_to_Z3"], "C5 amplitude roles"
            ),
            "amplitude_pointer": (
                "$.source_recipe_contract.ordinary_Z_roles.C5."
                "amplitude_ports_Z0_to_Z3"
            ),
        },
        "orders_pointer": "$.source_recipe_contract.ordinary_Z_roles.orders",
    }


def _read_D_m_rows(D_result, roles):
    """Read D ``m`` bindings and cross-check saved source-role fields."""
    rows = D_result["signed_D_history_rows"]["m"]
    _require(set(rows) == {"0", "1", "2", "3"},
             "D saved m rows are not exactly ordinary Z0 to Z3")
    selected = {}
    for order in ORDINARY_Z_ORDERS:
        row = rows[str(order)]
        for side in ("C3", "C5"):
            source = row[f"{side}_source"]
            role = roles[side]
            _require(source["namespace"] == role["manifest_namespace"],
                     f"D {side} saved row {order} namespace changed")
            _require(source["field"] == "m",
                     f"D {side} saved row {order} is not m")
            _require(source["ordinary_Z_order"] == order,
                     f"D {side} saved row {order} order changed")
            _require(source["source_node"] == role["m_history_handles"][order],
                     f"D {side} saved row {order} source node disagrees with role")
            _require(source["source_role"] == role["m_history_source_role"],
                     f"D {side} saved row {order} source role changed")
            _require(source["archive_raw_sha256"] == RAW_ARCHIVE_SHA256[side],
                     f"D {side} saved row {order} archive digest changed")
        _require(row["N_factor"] == 18,
                 f"D saved m row {order} changed its N factor")
        _require(row["C5_role_field"] == "native_N_scaled_terminal_history",
                 f"D C5 saved row {order} lost its N-scaled role")
        _require(row["C5_raw_correction_node"] ==
                 roles["C5"]["raw_correction_handles"][order],
                 f"D C5 saved row {order} raw correction node changed")
        _require(row["gap_value_selected"] is False,
                 f"D saved m row {order} selected a gap value")
        _require(row["common_physical_domain_claimed"] is False,
                 f"D saved m row {order} claimed a common domain")
        _require(row["physical_Pstar_lift_applied"] is False,
                 f"D saved m row {order} applied a physical lift")
        selected[order] = row
    return selected


def _read_M_base_source_recipes(joint_report, D_rows, roles):
    """Bind only the saved positive base-recipe handles at ordinary Z0.

    This accesses source-recipe records, never the joint report's M atom
    numerator/denominator or literal atom-term subtree.  It establishes a
    saved handle relation, not a cross-archive source-cut identity, function
    value, or physical function equality.
    """
    contract = joint_report["source_recipe_contract"]
    _require(contract["passed"] is True,
             "joint source-recipe contract did not pass")
    _require(contract["source_positive_certificates_retained"] is True,
             "joint positive-source premise is not retained")
    _require(contract["source_cuts_not_evaluated_or_identified_across_archives"]
             is True,
             "joint source-cut contract unexpectedly evaluated a source cut")
    recipes = contract["base_source_recipes"]
    _require(contract["base_source_recipe_count"] == 10,
             "joint base source-recipe count is not ten")
    _require(len(recipes) == contract["base_source_recipe_count"],
             "joint base source-recipe list length disagrees")
    selected = [recipe for recipe in recipes
                if recipe.get("target") == "M"]
    _require(len(selected) == 2,
             "saved M source recipes are not exactly C3 and C5")
    by_namespace = {}
    for recipe in selected:
        namespace = recipe.get("namespace")
        _require(namespace in ("C3", "C5"),
                 "saved M source recipe has an unexpected namespace")
        _require(namespace not in by_namespace,
                 f"duplicate saved M source recipe for {namespace}")
        _require(recipe["base_recipe_checked"] is True,
                 f"saved {namespace} M source recipe was not checked")
        by_namespace[namespace] = recipe
    _require(set(by_namespace) == {"C3", "C5"},
             "saved M source recipes do not cover C3 and C5")

    result = {
        "scope": "ordinary_Z0_only",
        "recipe_pointer":
            "$.source_recipe_contract.base_source_recipes[*]"
            " (target=M, namespace=C3|C5)",
        "base_recipe_checked": True,
        "positive_source_premise_retained": True,
        "cross_archive_source_cut_identity_proved": False,
        "saved_positive_base_source_recipe_handle_relation": True,
        "function_value_evaluated": False,
        "physical_function_equality_proved": False,
        "M_literal_atom_subtree_accessed": False,
    }
    for side in ("C3", "C5"):
        recipe = by_namespace[side]
        role = roles[side]
        D_m_source = D_rows[0][f"{side}_source"]
        _require(recipe["numerator"] == D_m_source["source_node"],
                 f"{side} Z0 recipe numerator disagrees with D m source")
        _require(recipe["numerator"] == role["m_history_handles"][0],
                 f"{side} Z0 recipe numerator disagrees with role")
        _require(recipe["denominator"] == role["amplitude_handles"][0],
                 f"{side} Z0 recipe denominator disagrees with amplitude role")
        _require(recipe["source_node"] == role["M_endpoint_handles"][0],
                 f"{side} Z0 recipe endpoint disagrees with M role")
        result[side] = {
            "namespace": side,
            "numerator_m_source_node": recipe["numerator"],
            "denominator_amplitude_source_node": recipe["denominator"],
            "M_endpoint_source_node": recipe["source_node"],
            "saved_positive_base_source_recipe_handle_relation": True,
        }
    return result


def _read_M_quotient_certificate(joint_report):
    """Return four quotient summaries without accessing the M atom subtree."""
    certificate = joint_report["joint_target_numerator_certificate"]
    _require(certificate["passed"] is True,
             "joint numerator certificate did not pass")
    _require(certificate["conditional"] is True,
             "joint quotient certificate is not marked conditional")
    denominator_contract = certificate["conditional_denominator_contract"]
    for key in ("A", "C", "N", "mu"):
        _require(isinstance(denominator_contract[key], str) and
                 denominator_contract[key],
                 f"joint quotient conditional prerequisite {key} is absent")
    _require(denominator_contract["numeric_intervals_used"] is False,
             "joint quotient conditional contract used numeric intervals")
    _require(certificate["quotient_product_certificate_passed"] is True,
             "joint quotient product certificate did not pass")
    jet_rows = certificate["quotient_jet_relations"]["M"]
    product_rows = certificate["quotient_product_certificate"]["M"]
    _require(len(jet_rows) == 4 and len(product_rows) == 4,
             "saved M quotient certificate is not four ordinary-Z rows")
    selected = {}
    for order in ORDINARY_Z_ORDERS:
        jet = jet_rows[order]
        product = product_rows[order]
        _require(jet["channel"] == "M" and
                 jet["ordinary_Z_order"] == order,
                 f"saved M quotient jet row {order} changed")
        _require(jet["triangular"] is True,
                 f"saved M quotient jet row {order} is not triangular")
        _require(jet["conditional_on_source_positive_denominator"] is True,
                 f"saved M quotient jet row {order} lacks its conditional premise")
        _require(product["channel"] == "M" and
                 product["ordinary_Z_order"] == order,
                 f"saved M product row {order} changed")
        _require(product["conditional_on_source_positive_denominator"] is True,
                 f"saved M product row {order} lacks its conditional premise")
        _require(product["passed"] is True,
                 f"saved M product certificate row {order} did not pass")
        _require(product["division_performed"] is False,
                 f"saved M product row {order} performed division")
        _require(product["quotient_jet_values_evaluated"] is False,
                 f"saved M product row {order} evaluated a quotient value")
        selected[order] = {
            "ordinary_Z_order": order,
            "quotient_jet_pointer": (
                "$.joint_target_numerator_certificate."
                f"quotient_jet_relations.M[{order}]"
            ),
            "quotient_product_pointer": (
                "$.joint_target_numerator_certificate."
                f"quotient_product_certificate.M[{order}]"
            ),
            "conditional_on_source_positive_denominator": True,
            "conditional_prerequisites_present": True,
            "positive_denominator_reproved_here": False,
            "triangular": True,
            "product_passed": True,
            "division_performed": False,
            "quotient_jet_values_evaluated": False,
            "M_atom_subtree_accessed": False,
            "source_recipe_scope": (
                "ordinary_Z0_positive_base_recipe"
                if order == 0 else "read_only_role_metadata"
            ),
            "endpoint_derivative_relation_proved": False,
        }
    return selected


def _read_scoped_amplitude_certificate(amplitude_report):
    """Read the scoped amplitude premise without turning it into a C5 claim."""
    consequence = amplitude_report["C3_Rc_amplitude_consequence"]
    source_binding = amplitude_report["source_binding"]
    certificate = amplitude_report["C3_outer_amplitude_identity_certificate"]
    _require(certificate["passed"] is True,
             "scoped C3 amplitude certificate did not pass")
    _require(consequence["ordinary_Z_orders"] == list(ORDINARY_Z_ORDERS),
             "amplitude certificate does not cover ordinary Z0 to Z3")
    _require(consequence["terminal_zero_not_claimed"] is True,
             "amplitude certificate unexpectedly claims a terminal zero")
    _require(source_binding["physical_Pstar_lift_claimed"] is False,
             "amplitude certificate unexpectedly claims a physical lift")
    return {
        "report_scope": consequence["scope"],
        "function_identity_literal": consequence["function_identity"],
        "ordinary_Z_orders": tuple(consequence["ordinary_Z_orders"]),
        "certificate_scope": certificate["scope"],
        "certificate_passed": True,
        "C3_rebuilt_source_leaf_ports": tuple(
            source_binding["rebuilt_source_leaf_ports"]
        ),
        "C3_primary_source_ports": tuple(source_binding["primary_source_ports"]),
        "C3_A3_rebuilt_source_leaf": consequence["rebuilt_A3_source_leaf_port"],
        "C3_A3_primary_source": consequence["primary_A3_source_port"],
        "rebuilt_expression_24177_used": (
            source_binding["rebuilt_expression_A3_24177_used"]
        ),
        "native_outer_to_C5_alias_flag": True,
        "exact_joint_C3_C5_vector_binding_proved": False,
        "a_zero_substitution_used": False,
        "a_retained_symbolically": True,
        "global_identity_admitted": False,
    }


def _build_order_bindings(
    D_rows, roles, quotient_rows, amplitude_scope, z0_source_recipe
):
    """Join D/source-role handles without binding typed M endpoint values."""
    bindings = []
    for order in ORDINARY_Z_ORDERS:
        row = D_rows[order]
        C3 = roles["C3"]
        C5 = roles["C5"]
        bindings.append({
            "ordinary_Z_order": order,
            "formula": FORMULA,
            "T_C_m_source": {
                "D_signed_row_pointer": row["signed_row_pointer"],
                "source_node_handle": row["C3_source"]["source_node"],
                "manifest_pointer": row["C3_source"]["manifest_pointer"],
                "source_role": row["C3_source"]["source_role"],
                "normalization_pointer": row["C3_source"]["normalization_pointer"],
            },
            "T_N_m_source": {
                "D_signed_row_pointer": row["signed_row_pointer"],
                "source_node_handle": row["C5_source"]["source_node"],
                "manifest_pointer": row["C5_source"]["manifest_pointer"],
                "source_role": row["C5_source"]["source_role"],
                "raw_correction_node_handle": row["C5_raw_correction_node"],
                "normalization_pointer": row["C5_source"]["normalization_pointer"],
                "N_factor_node": row["N_factor"],
            },
            "delta_m_source": {
                "formal_definition": "delta_m = T_N_m - T_C_m",
                "C3_m_source_row": row["C3_source"]["source_node"],
                "C5_m_source_row": row["C5_source"]["source_node"],
                "D_k_m_gap_role_field_present": False,
                "D_gap_role_verified": False,
                "source_binding_verified": False,
                "selected_value": False,
            },
            "A_C_source": {
                "source_node_handle": C3["amplitude_handles"][order],
                "source_role_pointer": C3["amplitude_pointer"],
            },
            "A_N_source": {
                "source_node_handle": C5["amplitude_handles"][order],
                "source_role_pointer": C5["amplitude_pointer"],
            },
            "C3_typed_M_endpoint_read_position": {
                "manifest_namespace": C3["manifest_namespace"],
                "role_namespace": C3["role_namespace"],
                "ordinary_Z_order": order,
                "source_node": C3["M_endpoint_handles"][order],
                "source_role_pointer": C3["M_endpoint_pointer"],
                "accessed_as": "read_position_only",
                "endpoint_source_bound": False,
            },
            "C5_typed_M_endpoint_read_position": {
                "manifest_namespace": C5["manifest_namespace"],
                "role_namespace": C5["role_namespace"],
                "ordinary_Z_order": order,
                "source_node": C5["M_endpoint_handles"][order],
                "source_role_pointer": C5["M_endpoint_pointer"],
                "accessed_as": "read_position_only",
                "endpoint_source_bound": False,
            },
            "saved_positive_base_source_recipe": (
                z0_source_recipe if order == 0 else None
            ),
            "saved_positive_base_source_recipe_handle_relation": order == 0,
            "quotient_certificate": quotient_rows[order],
            "scoped_amplitude_identity": {
                "certificate_report_scope": amplitude_scope["report_scope"],
                "a_zero_substitution_used": False,
                "a_retained_symbolically": True,
            },
            "endpoint_equals_T_over_A_proved": False,
            "endpoint_derivative_relation_proved": False,
            "common_point_proved": False,
            "same_units_proved": False,
        })
    return bindings


def load_M_formula_source_linkage(repo_root):
    """Load authenticated saved ports and return a non-numerical M packet.

    This function is intentionally not called in this draft.  If reviewed and
    invoked later, it reads saved JSON and embedded JSON only.  It never opens
    either gzip archive and never imports or evaluates a source expression.
    """
    root = Path(repo_root)
    D_report = _read_authenticated_json(root, D_REPORT, D_REPORT_SHA256)
    joint_report = _read_authenticated_json(
        root, JOINT_REPORT, JOINT_REPORT_SHA256
    )
    amplitude_report = _read_authenticated_json(
        root, AMPLITUDE_REPORT, AMPLITUDE_REPORT_SHA256
    )
    D_result = _read_embedded_D_result(D_report)
    _assert_publication_state(D_report, D_result, joint_report, amplitude_report)
    roles = _read_typed_roles(joint_report)
    D_rows = _read_D_m_rows(D_result, roles)
    z0_source_recipe = _read_M_base_source_recipes(
        joint_report, D_rows, roles
    )
    quotient_rows = _read_M_quotient_certificate(joint_report)
    amplitude_scope = _read_scoped_amplitude_certificate(amplitude_report)
    bindings = _build_order_bindings(
        D_rows, roles, quotient_rows, amplitude_scope, z0_source_recipe
    )

    return {
        "schema": SCHEMA,
        "source_only": True,
        "binding_kind": (
            "saved_D_m_plus_Z0_source_recipe_plus_quotient_certificate_metadata"
        ),
        "formula": FORMULA,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "inputs": {
            "D_publication": {
                "path": str(D_REPORT).replace("\\", "/"),
                "report_sha256": D_REPORT_SHA256,
                "embedded_result_sha256": D_EMBEDDED_RESULT_SHA256,
                "m_rows_read": 4,
                "D_atom_subtree_accessed": False,
            },
            "joint_publication": {
                "path": str(JOINT_REPORT).replace("\\", "/"),
                "report_sha256": JOINT_REPORT_SHA256,
                "M_quotient_rows_read": 4,
                "M_Z0_base_recipe_accessed": True,
                "M_endpoint_role_positions_listed": 8,
                "M_literal_atom_subtree_accessed": False,
            },
            "scoped_amplitude_publication": {
                "path": str(AMPLITUDE_REPORT).replace("\\", "/"),
                "report_sha256": AMPLITUDE_REPORT_SHA256,
                "scope": amplitude_scope["report_scope"],
                "a_zero_substitution_used": False,
                "joint_C3_C5_vector_binding_proved": False,
            },
            "raw_archive_sha256": RAW_ARCHIVE_SHA256,
            "gzip_decoded": False,
        },
        "source_port_binding": {
            "four_D_m_rows_bound": True,
            "D_m_namespace_field_order_source_node_checked": True,
            "D_N_factor_18_checked": True,
            "saved_positive_base_source_recipe_handle_relation": True,
            "Z0_saved_base_recipe_links_D_m_amplitude_and_M_endpoint_handles": True,
            "cross_archive_source_cut_identity_proved": False,
            "Z0_source_recipe_is_not_a_function_value": True,
            "Z0_physical_function_equality_proved": False,
            "four_typed_C3_C5_M_endpoint_rows_bound": False,
            "eight_typed_M_endpoint_positions_listed": True,
            "four_non_atom_quotient_certificate_rows_bound": False,
            "four_non_atom_quotient_certificate_rows_read_as_metadata": True,
            "Z1_Z3_endpoint_derivative_relations_proved": False,
            "scoped_amplitude_certificate_loaded": True,
            "formula_source_identity_proved": False,
            "endpoint_equals_T_over_A_proved": False,
        },
        "zero_order_saved_positive_base_source_recipe_handle_relation": (
            z0_source_recipe
        ),
        "ordinary_Z_bindings": bindings,
        "N_policy": {
            "N_node": 18,
            "invN_node": 19,
            "N_materialized": False,
            "extra_N_divisions": 0,
            "D_N_factor_all_orders": 18,
            "D_N_factor_consistency_checked": True,
            "finite_N_M_error_proved": False,
        },
        "normalization_and_units": {
            "literal": D_rows[0]["normalization_literal"],
            "common_physical_domain_proved": False,
            "same_units_proved": False,
            "physical_Pstar_lift_applied": False,
        },
        "M_vs_m15": {
            "M_equals_m15_proved": False,
            "physical_m15_map_entered": False,
            "distinct_symbols_retained": True,
        },
        "reference_offset": {
            "status": "inherited_unresolved_symbolic_marker_only",
            "independent": True,
            "numeric_offset_selected": False,
            "cancelled_in_M_formula": False,
        },
        "prior_generic_failure": PRIOR_GENERIC_M_Z0_FAILURE,
        "retained_legacy_scope": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "J_terms": 16,
            "K_terms": 17,
            "outer_B_terms": 2,
            "all_old_error_data_retained": True,
            "cap_transfer": False,
        },
        "nonclaims": {
            "signed_gap_value": False,
            "M_gap_zero": False,
            "M_Z0_function_value": False,
            "M_Z0_physical_function_equality": False,
            "M_Z1_Z3_endpoint_derivative_relations": False,
            "common_physical_point": False,
            "common_physical_domain": False,
            "same_units": False,
            "finite_N_error_bound": False,
            "physical_Pstar_lift": False,
            "M_equals_m15": False,
            "reference_offset_cancellation": False,
            "cap_transfer": False,
            "a_zero": False,
        },
        "audit": {
            "generation_scope":
                "call_time_output_of_load_M_formula_source_linkage",
            "authoring_draft_invoked": False,
            "source_binding_is_not_an_authoring_time_observation": True,
            "source_modules_imported": 0,
            "source_functions_executed": 0,
            "source_graph_compared": 0,
            "saved_role_namespace_field_order_node_consistency_checked": True,
            "numeric_values_evaluated": 0,
            "formula_identity_replayed": 0,
            "M_Z0_literal_atom_catalog_entered": False,
            "M_Z0_literal_numerator_reentered": False,
            "M_Z0_literal_denominator_reentered": False,
            "generic_M_Z0_family_reentered": False,
            "prior_failure_reentered": False,
            "gzip_archives_decoded": 0,
        },
        "remaining_scope": [
            "Treat the Z0 saved positive base-recipe handle relation only; it is not a cross-archive source-cut identity, function value, or physical function equality.",
            "Keep Z1, Z2 and Z3 typed M endpoint positions as read-only metadata and prove their true derivative relations separately.",
            "Keep endpoint equality to T_m/A conditional at every order until a source function certificate is available.",
            "Bind the scoped amplitude identity to the joint C3/C5 vectors before any a=0 use.",
            "Prove common physical point, coordinate map, and unit/scale convention.",
            "Derive a separate finite-N M error bound; do not reuse a cap as a value.",
            "Keep M distinct from m15 and retain the independent reference offset.",
        ],
    }


# No __main__ entry point: intentionally unexecuted draft.
