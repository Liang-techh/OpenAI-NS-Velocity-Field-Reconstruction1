"""UNEXECUTED I/S/Cp-only source-port draft.

This is a new, deliberately narrow source consumer design for the remaining
``h/e/p`` signed-history rows and the saved ordinary-Z ``I/S/Cp`` numerator /
denominator atom occurrences.  It never imports a reconstruction module,
executes a source function, evaluates a gap, reconstructs a quotient, or
replays a formula proof.  It does not iterate the ``M`` or ``D`` channels.

The prior D-only publication is inherited as immutable metadata only.  Its
eight accepted ``k/m`` rows are not read, recomputed, or rebound here.  The
independent reference offset, common physical chart, units, Pstar lift,
physical ratio, and cap transfer remain explicit nonclaims.

There is intentionally no ``__main__`` entry point.  This file is a draft and
must remain unexecuted until the root agent separately reviews and freezes it.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


SCHEMA = "compact_ISCp_gap_atoms_new_draft_v1"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
ISCP_FIELDS = ("h", "e", "p")
ISCP_CHANNELS = ("I", "S", "Cp")
CHANNEL_FOR_FIELD = {"h": "I", "e": "S", "p": "Cp"}

SIGNED_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "signed_history_defect_transport.json"
)
JOINT_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "joint_target_numerator_transport.json"
)
RAW_ARCHIVES = {
    "C3": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz"
    ),
    "C5": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
    ),
}

SIGNED_REPORT_SHA256 = (
    "dc367908bd44991ef6f4a973d688b4bd27bb0cab43464e0e396b3abb1201c861"
)
JOINT_REPORT_SHA256 = (
    "cff727e06c1af18eebeba0bbf027a488128c2021f9feb839bafdfc368b18f133"
)
RAW_ARCHIVE_SHA256 = {
    "C3": "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "C5": "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
}

# These are source-node identifiers, not source-function values.  The row
# layout is m, h, k, e, p in the saved twenty-row packet; only h/e/p offsets
# are entered by this draft.
HISTORY_ROW_START = {"h": 4, "e": 12, "p": 16}
EXPECTED_C3_HISTORY = {
    "h": (4928, 4929, 12043, 24178),
    "e": (4932, 4933, 12045, 24180),
    "p": (4888, 4889, 12041, 24176),
}
EXPECTED_C5_HISTORY = {
    "h": (258889, 258890, 258891, 258892),
    "e": (258901, 258902, 258903, 258904),
    "p": (258907, 258908, 258909, 258910),
}
EXPECTED_C5_RAW_CORRECTION = {
    "h": (258652, 258653, 258654, 258655),
    "e": (258664, 258665, 258666, 258667),
    "p": (258670, 258671, 258672, 258673),
}
EXPECTED_C3_AMPLITUDE = (4934, 4935, 12046, 24181)
EXPECTED_C5_AMPLITUDE = (251009, 251010, 251011, 251012)

SIGNED_NORMALIZATION = (
    "generic common E=Utheta/Pstar, V=Uz/Pstar, "
    "m/k scaled1/Pstar, h/e/p retained per original moment definitions"
)

EXPECTED_RETAINED_CONSTRUCTION_DATA = {
    "all_old_error_data_retained": True,
    "cumulative_terms": 212,
    "inherited_errors_deleted": False,
    "raw_signed_terms": 213,
    "reduced_terms": 183,
    "reserved_contribution_registered_once": True,
    "six_signed_channels_retained": True,
}
EXPECTED_RETAINED_INTEGRAL_INVENTORY = {
    "J_terms": 16,
    "J_terms_retained": True,
    "K_terms": 17,
    "K_terms_retained": True,
    "outer_B_terms": 2,
    "outer_B_terms_retained": True,
    "reserved_contribution_registered_once": True,
}
OLD_ERROR_FLAGS = (
    "old212_adapter_identity",
    "old212_quantclosure",
    "old_errors_removed",
    "old_error_terms_removed",
)
REFERENCE_OFFSET_MARKER = "independent old reference offset"

# Inherited report only.  No helper, source file, or D atom is imported or
# called; this is a provenance pointer to the already accepted D-only packet.
INHERITED_D_ONLY_REPORT = {
    "path": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_compact_D_km_source_binding.json"
    ),
    "sha256": "9af596c029d7eadd78b25a7895ffc8f4f7f24aefcea0a3b6f4f07f6c28669587",
    "saved_result_path": "tmp/compact_D_km_source_binding_v1/result.json",
    "saved_result_sha256": "0ebea65173f200afcee2a917804d3674861c66b24bd6120c17612c993a845e24",
    "status": "accepted_limited_source_port_binding",
    "accepted_history_rows": 8,
    "accepted_channel": "D",
    "metadata_only_inheritance": True,
    "reentered": False,
}

# Static failure evidence only.  The failed generic v1 source consumer is not
# loaded, copied, or entered by this file.
PRIOR_GENERIC_V1_FAILURE = {
    "result_path": "tmp/compact_joint_history_gap_binding_v1/result.json",
    "result_sha256": "50257d1c62728ba554228f69e68d50d4c4be673ff0b0bf13e79b3c6b75b3be3c",
    "status": "failed",
    "exception": "KeyError: 0",
    "reentered": False,
    "proof_body_reused": False,
    "reason_not_reentered": (
        "the saved generic family indexed field-keyed history nodes before "
        "selecting the field; this draft has a separate h/e/p-only entry path"
    ),
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha256(path):
    """Hash bytes only; this helper does not decode a gzip archive."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _authenticate_inherited_D_bytes(root):
    """Pin saved D report/result bytes only; never decode or replay D."""
    report = Path(root) / INHERITED_D_ONLY_REPORT["path"]
    result = Path(root) / INHERITED_D_ONLY_REPORT["saved_result_path"]
    _require(_sha256(report) == INHERITED_D_ONLY_REPORT["sha256"], "inherited D publication bytes changed")
    _require(_sha256(result) == INHERITED_D_ONLY_REPORT["saved_result_sha256"], "inherited D result bytes changed")
    return {**INHERITED_D_ONLY_REPORT, "report_and_result_raw_sha256_authenticated": True,
            "D_report_or_result_decoded": False}


def _read_authenticated_json(root, relative, expected_sha256):
    path = Path(root) / relative
    _require(path.is_file(), f"missing authenticated JSON: {relative}")
    _require(_sha256(path) == expected_sha256,
             f"authenticated JSON SHA changed: {relative}")
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"saved JSON is not an object: {relative}")
    return value


def _authenticate_raw_archives(root):
    """Check authoritative raw SHA values before any possible gzip decode."""
    checked = {}
    for namespace in ("C3", "C5"):
        relative = RAW_ARCHIVES[namespace]
        path = Path(root) / relative
        _require(path.is_file(), f"missing raw {namespace} archive: {relative}")
        digest = _sha256(path)
        _require(digest == RAW_ARCHIVE_SHA256[namespace],
                 f"raw {namespace} archive SHA changed before decode")
        checked[namespace] = {
            "relative_path": str(relative).replace("\\", "/"),
            "raw_sha256": digest,
            "raw_sha256_checked_before_gzip_decode": True,
            "gzip_decoded": False,
        }
    return checked


def _audit_saved_flags(signed, joint):
    """Retain source-only guards without evaluating any saved source node."""
    for label, report in (
        ("signed.top", signed),
        ("signed.full", signed["full_result"]),
        ("signed.component", signed["component_result"]),
    ):
        _require(report["source_only"] is True, f"{label} is not source-only")
        _require(report["original_source_evaluations"] == 0,
                 f"{label} evaluated an original source")
        _require(report["signed_history_source_gap_functions_defined"] is True,
                 f"{label} lacks saved signed-history rows")
        _require(report["native_N_scaling_applied_once_bound"] is True,
                 f"{label} lost one-time N scaling")

    for label, report in (
        ("joint.top", joint),
        ("joint.full", joint["full_result"]),
        ("joint.component", joint["component_result"]),
    ):
        _require(report["source_only"] is True, f"{label} is not source-only")
        _require(report["signed_target_gaps_evaluated"] is False,
                 f"{label} contains evaluated signed gaps")
        _require(report["numeric_values_computed"] is False,
                 f"{label} contains numeric values")
        _require(report["joint_target_numerator_source_defined"] is True,
                 f"{label} lacks saved joint source packet")
        _require(report["actual_C3_C5_incoming_source_identity_proved"] is False,
                 f"{label} upgrades incoming source identity")
        binding = report["source_binding"]
        _require(binding["actual_target_gap_zero_proved"] is False,
                 f"{label} claims a zero target gap")
        _require(binding["physical_Pstar_lift_applied"] is False,
                 f"{label} claims a physical Pstar lift")
        _require(report["cap_to_source_identity_proved"] is False,
                 f"{label} claims cap transfer")


def _retained_legacy_snapshot(signed, joint):
    """Return the inherited inventory and unresolved offset marker verbatim."""
    reports = (
        ("signed.top", signed),
        ("signed.full", signed["full_result"]),
        ("signed.component", signed["component_result"]),
        ("joint.top", joint),
        ("joint.full", joint["full_result"]),
        ("joint.component", joint["component_result"]),
    )
    for label, report in reports:
        _require(report["retained_construction_data"] ==
                 EXPECTED_RETAINED_CONSTRUCTION_DATA,
                 f"{label} changed retained 212/213/183 inventory")
        _require(report["retained_integral_term_inventory"] ==
                 EXPECTED_RETAINED_INTEGRAL_INVENTORY,
                 f"{label} changed retained integral inventory")
        for flag in OLD_ERROR_FLAGS:
            _require(report[flag] is False, f"{label} changed old error {flag}")
        caps = report["retained_cap_provenance"]
        _require(caps["all_original_caps_retained"] is True and
                 caps["cap_to_source_bridge_proved"] is False and
                 caps["numeric_centers_used_as_source_values"] is False,
                 f"{label} changed retained cap provenance")

    markers = []
    for label, report in (
        ("joint.top", joint),
        ("joint.full", joint["full_result"]),
        ("joint.component", joint["component_result"]),
    ):
        scope = report["remaining_scope"]
        _require(isinstance(scope, list), f"{label} remaining scope is not a list")
        matches = [
            item for item in scope
            if isinstance(item, str) and REFERENCE_OFFSET_MARKER in item
        ]
        _require(matches, f"{label} lost independent reference-offset marker")
        markers.append({
            "report": label,
            "json_pointer": "$.remaining_scope[*]",
            "marker": matches[0],
        })
    return {
        "retained_construction_data": EXPECTED_RETAINED_CONSTRUCTION_DATA,
        "retained_integral_term_inventory": EXPECTED_RETAINED_INTEGRAL_INVENTORY,
        "old_error_flags": {name: False for name in OLD_ERROR_FLAGS},
        "reference_offset": {
            "status": "inherited_unresolved_symbolic_marker_only",
            "numeric_offset_selected": False,
            "markers": markers,
        },
        "cap_transfer_claimed": False,
    }


def _read_N_handles(signed):
    """Preserve N node handles; never materialize N or evaluate 1/N."""
    full = signed["full_result"]
    terminal = full["terminal_history_common_N_binding"]
    literal = full["source_binding"]["common_N_literal_binding"]
    _require(terminal["N_node"] == 18 and terminal["invN_node"] == 19,
             "saved N handles changed")
    _require(terminal["J"] == 1358356628656378313 and
             terminal["N_integer_materialized"] is False,
             "saved unmaterialized N descriptor changed")
    _require(terminal["C3_archive_sha256"] == RAW_ARCHIVE_SHA256["C3"] and
             terminal["native_archive_sha256"] == RAW_ARCHIVE_SHA256["C5"],
             "saved N archive binding changed")
    _require(literal["C3_node_18"] == 18 and literal["native_node_18"] == 18 and
             literal["C3_invN_node_19"] == 19 and
             literal["native_invN_node_19"] == 19 and
             literal["N_integer_materialized"] is False,
             "saved common-N literal handles changed")
    return {
        "N_node": 18,
        "invN_node": 19,
        "J": 1358356628656378313,
        "definition": "N=2^(2^J)",
        "type": "DoubleDyadicInteger",
        "N_Z_independent": True,
        "N_integer_materialized": False,
        "invN_operation": "positive_quotient(numerator=1, denominator_node=18)",
        "source_pointers": {
            "terminal": "$.full_result.terminal_history_common_N_binding",
            "literal": "$.full_result.source_binding.common_N_literal_binding",
        },
    }


def _roles_and_manifest(joint):
    contract = joint["source_recipe_contract"]
    roles = contract["ordinary_Z_roles"]
    _require(roles["orders"] == list(ORDINARY_Z_ORDERS),
             "ordinary-Z role orders changed")
    _require(roles["orders_are_source_roles"] is True and
             roles["range_values_are_not_function_values"] is True,
             "ordinary-Z role semantics changed")
    _require(roles["C3"]["namespace"] == "repair_C3_band" and
             roles["C5"]["namespace"] == "true_compact_C5",
             "C3/C5 source namespaces changed")

    manifest = joint["joint_target_numerator_certificate"]["source_gap_manifest"]
    _require(isinstance(manifest, list) and len(manifest) == 2,
             "joint source-gap manifest shape changed")
    by_namespace = {}
    for index, entry in enumerate(manifest):
        namespace = entry["namespace"]
        _require(namespace in ("C3", "C5") and namespace not in by_namespace,
                 "joint manifest namespace changed")
        _require(entry["archive_raw_sha256"] == RAW_ARCHIVE_SHA256[namespace],
                 f"{namespace} manifest archive changed")
        role = roles[namespace]
        history_key = (
            "terminal_history_ports_Z0_to_Z3"
            if namespace == "C3"
            else "terminal_scaled_history_ports_Z0_to_Z3"
        )
        _require(entry["amplitude_source_nodes_Z0_to_Z3"] ==
                 role["amplitude_ports_Z0_to_Z3"],
                 f"{namespace} amplitude role changed")
        _require({key: tuple(entry["history_source_nodes_Z0_to_Z3"][key]) for key in ISCP_FIELDS}
                 == {key: tuple(role[history_key][key]) for key in ISCP_FIELDS},
                 f"{namespace} history role changed")
        by_namespace[namespace] = {"index": index, "entry": entry}
    _require(set(by_namespace) == {"C3", "C5"},
             "joint source-gap manifest incomplete")
    return roles, by_namespace


def _port_witness(roles, manifest, namespace, kind, field, order, source_node):
    """Attach a manifest witness to one literal source atom."""
    _require(namespace in ("C3", "C5"), "invalid source namespace")
    _require(order in ORDINARY_Z_ORDERS, "invalid ordinary-Z order")
    role = roles[namespace]
    if kind == "amplitude":
        _require(field == "A", "amplitude field must be A")
        declared = role["amplitude_ports_Z0_to_Z3"]
        pointer_tail = f"amplitude_source_nodes_Z0_to_Z3[{order}]"
    elif kind == "history":
        _require(field in ISCP_FIELDS, "history field outside I/S/Cp scope")
        history_key = (
            "terminal_history_ports_Z0_to_Z3"
            if namespace == "C3"
            else "terminal_scaled_history_ports_Z0_to_Z3"
        )
        declared = role[history_key][field]
        pointer_tail = f"history_source_nodes_Z0_to_Z3.{field}[{order}]"
    else:
        raise ValueError("unsupported source witness kind")
    _require(declared[order] == source_node,
             "saved atom node does not match ordinary-Z source role")
    manifest_index = manifest[namespace]["index"]
    return {
        "namespace": namespace,
        "contract_namespace": role["namespace"],
        "archive_raw_sha256": RAW_ARCHIVE_SHA256[namespace],
        "kind": kind,
        "field": field,
        "ordinary_Z_order": order,
        "source_node": source_node,
        "manifest_pointer": (
            "$.joint_target_numerator_certificate.source_gap_manifest["
            f"{manifest_index}].{pointer_tail}"
        ),
    }


def _make_iscp_atom_index(roles, manifest, history):
    """Index only A and h/e/p tokens used by selected I/S/Cp rows."""
    index = {}
    for namespace in ("C3", "C5"):
        role = roles[namespace]
        for order, source_node in enumerate(role["amplitude_ports_Z0_to_Z3"]):
            token = f"tn_{namespace}_A_node_{source_node}"
            index[token] = _port_witness(
                roles, manifest, namespace, "amplitude", "A", order, source_node
            )
        history_key = (
            "terminal_history_ports_Z0_to_Z3"
            if namespace == "C3"
            else "terminal_scaled_history_ports_Z0_to_Z3"
        )
        for field in ISCP_FIELDS:
            for order, source_node in enumerate(role[history_key][field]):
                token = f"tn_{namespace}_{field}_node_{source_node}"
                witness = _port_witness(
                    roles, manifest, namespace, "history", field, order,
                    source_node,
                )
                row = next(item for item in history if item["field"] == field and item["ordinary_Z_order"] == order)
                source_key = "C3_source" if namespace == "C3" else "C5_source"
                _require(row[source_key]["source_node"] == source_node,
                         "signed source node does not match joint history atom")
                witness["signed_row_pointer"] = row["signed_row_pointer"]
                witness["normalization_pointer"] = row["normalization_pointer"]
                index[token] = witness
    # mu is a saved symbolic positive parameter, not a selected value.
    index["tn_mu"] = {
        "namespace": "joint-common",
        "kind": "positive_common_parameter",
        "field": "mu",
        "ordinary_Z_order": None,
        "source_node": None,
        "source_function_value_selected": False,
        "json_pointer": (
            "$.joint_target_numerator_certificate."
            "conditional_denominator_contract.mu"
        ),
    }
    return index


def _bind_iscp_history_rows(signed, roles, manifest):
    """Bind only the twelve saved h/e/p rows, preserving signed correction data."""
    rows = signed["full_result"]["signed_history_gap_function_nodes"]
    _require(isinstance(rows, list) and len(rows) == 20,
             "saved signed-history packet is not twenty rows")
    bound = []
    for field in ISCP_FIELDS:
        start = HISTORY_ROW_START[field]
        for order in ORDINARY_Z_ORDERS:
            row_index = start + order
            row = rows[row_index]
            c3_node = EXPECTED_C3_HISTORY[field][order]
            c5_node = EXPECTED_C5_HISTORY[field][order]
            raw_node = EXPECTED_C5_RAW_CORRECTION[field][order]
            _require(row["gap_index"] == row_index and
                     row["field"] == field and
                     row["ordinary_Z_order"] == order and
                     row["operation"] == "signed_difference",
                     "selected h/e/p signed-history row key changed")
            left = row["left"]
            right = row["right"]
            _require(left["namespace"] == "true_compact_C5" and
                     left["archive_sha256"] == RAW_ARCHIVE_SHA256["C5"] and
                     left["role"] == "native_N_scaled_terminal_history" and
                     left["source_node"] == c5_node and
                     left["raw_correction_node"] == raw_node and
                     left["raw_correction_role"] == "native_function_substitution_at_coordinate_2" and
                     left["N_factor"] == 18 and
                     left["N_scaling_applied_once"] is True,
                     "selected C5 signed-history source row changed")
            _require(right["namespace"] == "repair_C3_band" and
                     right["archive_sha256"] == RAW_ARCHIVE_SHA256["C3"] and
                     right["role"] == "C3_rebuilt_terminal_N_scaled_history" and
                     right["source_node"] == c3_node and
                     right["already_N_scaled"] is True and
                     right["source_role"] ==
                     "source_recipe_contract.ordinary_Z_roles.C3."
                     f"terminal_history_ports_Z0_to_Z3.{field}",
                     "selected C3 signed-history source row changed")
            _require(roles["C3"]["terminal_history_ports_Z0_to_Z3"][field][order] == right["source_node"] and
                     roles["C5"]["terminal_scaled_history_ports_Z0_to_Z3"][field][order] == left["source_node"] and
                     roles["C5"]["raw_correction_ports_Z0_to_Z3"][field][order] == left["raw_correction_node"] and
                     manifest["C3"]["entry"]["history_source_nodes_Z0_to_Z3"][field][order] == right["source_node"] and
                     manifest["C5"]["entry"]["history_source_nodes_Z0_to_Z3"][field][order] == left["source_node"],
                     "signed history source does not match joint manifest")
            _require(row["normalization"] == SIGNED_NORMALIZATION and
                     row["values_uncomputed"] is True and
                     row["gap_zero_claimed"] is False and
                     row["same_units_claimed"] is False and
                     row["conditional_common_point_pullback_unbound"] is True,
                     "selected h/e/p signed-history guard changed")
            n_binding = row["common_N_literal_binding"]
            _require(n_binding == {"C3_invN_node_19": 19, "C3_node_18": 18,
                                   "J": 1358356628656378313, "N_integer_materialized": False,
                                   "native_invN_node_19": 19, "native_node_18": 18},
                     "selected signed row has a different common N source")
            bound.append({
                "field": field,
                "channel": CHANNEL_FOR_FIELD[field],
                "ordinary_Z_order": order,
                "signed_row_pointer": (
                    "$.full_result.signed_history_gap_function_nodes["
                    f"{row_index}]"
                ),
                "signed_gap_operation": "saved signed difference only",
                "C5_source": {
                    "namespace": left["namespace"],
                    "archive_raw_sha256": left["archive_sha256"],
                    "source_node": left["source_node"],
                    "raw_correction_node": left["raw_correction_node"],
                    "raw_correction_role": left["raw_correction_role"],
                    "role": left["role"],
                    "source_role_pointer": "$.source_recipe_contract.ordinary_Z_roles.C5.terminal_scaled_history_ports_Z0_to_Z3." + field,
                    "N_factor": left["N_factor"],
                    "N_scaling_applied_once": left["N_scaling_applied_once"],
                },
                "C3_source": {
                    "namespace": right["namespace"],
                    "archive_raw_sha256": right["archive_sha256"],
                    "source_node": right["source_node"],
                    "already_N_scaled": right["already_N_scaled"],
                    "role": right["role"],
                    "source_role": right["source_role"],
                },
                "normalization_literal": row["normalization"],
                "normalization_pointer": (
                    "$.full_result.signed_history_gap_function_nodes["
                    f"{row_index}].normalization"
                ),
                "physical_Pstar_lift_applied": False,
                "extra_physical_ratio_claimed": False,
                "common_physical_domain_claimed": False,
                "same_units_claimed": False,
                "gap_zero_claimed": False,
            })
    return bound


def _atom_bindings_for_row(row, channel, order, atom_index, side):
    """Bind actual saved monomial tokens without rewriting their formula."""
    terms = row["terms"]
    _require(isinstance(terms, list),
             f"saved {channel} {side} terms are not a list")
    bound_terms = []
    for term_index, term in enumerate(terms):
        monomial = term["monomial"]
        _require(isinstance(monomial, list),
                 f"saved {channel} {side} monomial is not a list")
        atoms = []
        for atom_index_in_term, token in enumerate(monomial):
            _require(token in atom_index,
                     f"unscoped saved atom in selected {channel} row: {token}")
            atoms.append({
                "atom_index": atom_index_in_term,
                "symbolic_atom": token,
                "occurrence_pointer": (
                    "$.joint_target_numerator_certificate."
                    "numerators_and_denominators."
                    f"{channel}.{side}[{order}].terms[{term_index}]"
                    f".monomial[{atom_index_in_term}]"
                ),
                "source_witness": atom_index[token],
            })
        bound_terms.append({
            "term_index": term_index,
            "coefficient_literal": term["coefficient"],
            "monomial_literal": monomial,
            "atom_bindings": atoms,
        })
    return bound_terms


def _bind_actual_iscp_atoms(joint, atom_index):
    """Read only I/S/Cp ordinary-Z rows from the saved joint packet."""
    formulas = joint["joint_target_numerator_certificate"]
    rows_by_channel = formulas["numerators_and_denominators"]
    bindings = {}
    for channel in ISCP_CHANNELS:
        channel_rows = rows_by_channel[channel]
        field = {"I": "h", "S": "e", "Cp": "p"}[channel]
        _require(channel_rows["raw_gap_role"] == "delta_" + field,
                 f"{channel} raw gap role does not match signed field")
        channel_binding = {}
        for order in ORDINARY_Z_ORDERS:
            numerator = channel_rows["numerator"][order]
            denominator = channel_rows["denominator"][order]
            _require(numerator["ordinary_Z_order"] == order and
                     denominator["ordinary_Z_order"] == order,
                     f"saved {channel} ordinary-Z row changed")
            numerator_tokens = {token for term in numerator["terms"] for token in term["monomial"]}
            denominator_tokens = {token for term in denominator["terms"] for token in term["monomial"]}
            _require({f"tn_C3_{field}_node_{EXPECTED_C3_HISTORY[field][order]}",
                      f"tn_C5_{field}_node_{EXPECTED_C5_HISTORY[field][order]}"}.issubset(numerator_tokens),
                     f"{channel} numerator lacks current-order C3/C5 history atoms")
            _require({f"tn_C3_A_node_{EXPECTED_C3_AMPLITUDE[order]}",
                      f"tn_C5_A_node_{EXPECTED_C5_AMPLITUDE[order]}"}.issubset(denominator_tokens),
                     f"{channel} denominator lacks current-order C3/C5 amplitude atoms")
            channel_binding[str(order)] = {
                "channel": channel,
                "ordinary_Z_order": order,
                "raw_gap_role": channel_rows["raw_gap_role"],
                "raw_gap_role_pointer": "$.joint_target_numerator_certificate.numerators_and_denominators." + channel + ".raw_gap_role",
                "numerator_pointer": (
                    "$.joint_target_numerator_certificate."
                    "numerators_and_denominators."
                    f"{channel}.numerator[{order}]"
                ),
                "denominator_pointer": (
                    "$.joint_target_numerator_certificate."
                    "numerators_and_denominators."
                    f"{channel}.denominator[{order}]"
                ),
                "numerator_terms": _atom_bindings_for_row(
                    numerator, channel, order, atom_index, "numerator"
                ),
                "denominator_terms": _atom_bindings_for_row(
                    denominator, channel, order, atom_index, "denominator"
                ),
                "source_values_selected": False,
                "numeric_gap_evaluated": False,
                "formula_identity_replayed": False,
                "physical_Pstar_lift_applied": False,
                "extra_physical_ratio_claimed": False,
                "cap_transfer_claimed": False,
            }
        bindings[channel] = channel_binding
    return bindings


def bind_iscp_saved_ports(repo_root):
    """Construct the unexecuted design's h/e/p and I/S/Cp literal packet."""
    root = Path(repo_root)
    inherited_D = _authenticate_inherited_D_bytes(root)
    raw_archives = _authenticate_raw_archives(root)
    signed = _read_authenticated_json(
        root, SIGNED_REPORT, SIGNED_REPORT_SHA256
    )
    joint = _read_authenticated_json(
        root, JOINT_REPORT, JOINT_REPORT_SHA256
    )
    _audit_saved_flags(signed, joint)
    retained = _retained_legacy_snapshot(signed, joint)
    n_handles = _read_N_handles(signed)
    roles, manifest = _roles_and_manifest(joint)
    history = _bind_iscp_history_rows(signed, roles, manifest)
    atom_index = _make_iscp_atom_index(roles, manifest, history)
    atoms = _bind_actual_iscp_atoms(joint, atom_index)
    return {
        "schema": SCHEMA,
        "source_only": True,
        "binding_kind": "iscp_h_e_p_direct_saved_json_ports_only",
        "scope": {
            "history_fields": list(ISCP_FIELDS),
            "channels": list(ISCP_CHANNELS),
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "signed_history_rows_bound": len(history),
            "numerator_rows_bound": len(ISCP_CHANNELS) * len(ORDINARY_Z_ORDERS),
            "denominator_rows_bound": len(ISCP_CHANNELS) * len(ORDINARY_Z_ORDERS),
            "M_or_D_traversed": False,
            "not_entered_channels": ["M", "D"],
        },
        "inputs": {
            "signed_history_report": {
                "path": str(SIGNED_REPORT).replace("\\", "/"),
                "sha256": SIGNED_REPORT_SHA256,
            },
            "joint_numerator_report": {
                "path": str(JOINT_REPORT).replace("\\", "/"),
                "sha256": JOINT_REPORT_SHA256,
            },
            "raw_archives": raw_archives,
            "inherited_D_only_report": inherited_D,
        },
        "prior_generic_v1_failure": PRIOR_GENERIC_V1_FAILURE,
        "N_handles": n_handles,
        "h_e_p_history_bindings": history,
        "actual_I_S_Cp_numerator_denominator_atom_bindings": atoms,
        "retained_legacy": retained,
        "nonclaims": {
            "signed_target_gaps_evaluated": False,
            "actual_gap_zero": False,
            "common_physical_domain": False,
            "same_units": False,
            "physical_Pstar_lift": False,
            "extra_physical_ratio": False,
            "cap_transfer": False,
            "reference_offset_selected": False,
            "old_errors_removed": False,
        },
        "scope_audit": {
            "source_modules_imported": 0,
            "source_functions_executed": 0,
            "numeric_values_evaluated": 0,
            "formula_identity_replayed": 0,
            "prior_generic_v1_reentered": False,
            "D_only_leaf_reentered": False,
            "gzip_archives_decoded": 0,
        },
        "json_pointers": {
            "selected_signed_rows": [row["signed_row_pointer"] for row in history],
            "joint_selected_atoms": [atoms[channel][str(order)][side + "_pointer"]
                                     for channel in ISCP_CHANNELS
                                     for order in ORDINARY_Z_ORDERS
                                     for side in ("numerator", "denominator")],
            "N_handles": "$.full_result.terminal_history_common_N_binding",
            "reference_offset": "$.remaining_scope[3]",
        },
        "remaining_scope": [
            "Bind no common physical point or unit convention in this packet.",
            "Do not select source-function values, gap values, or extra physical ratios.",
            "Keep the independent reference offset as an unresolved symbolic marker.",
            "Retain all old 212/213/183 data, errors, and signed corrections.",
        ],
    }


# No __main__ entry point: this draft is intentionally unexecuted.
