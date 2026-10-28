"""UNEXECUTED D-channel-only source-bound successor draft.

This is intentionally narrower and independent of the rejected history/atom
families.  It binds only the saved ``k`` and ``m`` signed-history rows at
ordinary ``Z0..Z3`` to the already saved D-channel numerator/denominator atom
strings, and records the independent reference-offset marker as an unresolved
source pointer.  It does not inspect M/I/S/Cp rows, does not build a shared atom
catalog, does not replay a numerator identity, and does not evaluate a source.

Only standard-library JSON and raw-byte hashing are present.  There is no
``__main__`` entry point; this file must remain unexecuted until root review.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


SCHEMA = "compact_D_km_reference_offset_source_binding_draft_v1"
ORDERS = (0, 1, 2, 3)
D_FIELDS = ("k", "m")

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
ARCHIVE_PATHS = {
    "C3": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz"
    ),
    "C5": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
    ),
}
REPORT_SHA256 = {
    "signed": "dc367908bd44991ef6f4a973d688b4bd27bb0cab43464e0e396b3abb1201c861",
    "joint": "cff727e06c1af18eebeba0bbf027a488128c2021f9feb839bafdfc368b18f133",
}
ARCHIVE_SHA256 = {
    "C3": "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "C5": "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
}

EXPECTED_C3 = {
    "amplitude": (4934, 4935, 12046, 24181),
    "k": (4930, 4931, 12044, 24179),
    "m": (4920, 4921, 12042, 24177),
}
EXPECTED_C5 = {
    "amplitude": (251009, 251010, 251011, 251012),
    "k": (258895, 258896, 258897, 258898),
    "m": (258883, 258884, 258885, 258886),
}
EXPECTED_C5_RAW = {
    "k": (258658, 258659, 258660, 258661),
    "m": (258646, 258647, 258648, 258649),
}
NORMALIZATION_LITERAL = (
    "generic common E=Utheta/Pstar, V=Uz/Pstar, "
    "m/k scaled1/Pstar, h/e/p retained per original moment definitions"
)
REFERENCE_OFFSET_MARKER = "independent old reference offset"
OLD_ERROR_FLAGS = (
    "old212_adapter_identity",
    "old212_quantclosure",
    "old_errors_removed",
    "old_error_terms_removed",
)
EXPECTED_OLD_DATA = {
    "all_old_error_data_retained": True,
    "cumulative_terms": 212,
    "inherited_errors_deleted": False,
    "raw_signed_terms": 213,
    "reduced_terms": 183,
    "reserved_contribution_registered_once": True,
    "six_signed_channels_retained": True,
}
EXPECTED_INTEGRAL_INVENTORY = {
    "J_terms": 16,
    "J_terms_retained": True,
    "K_terms": 17,
    "K_terms_retained": True,
    "outer_B_terms": 2,
    "outer_B_terms_retained": True,
    "reserved_contribution_registered_once": True,
}


def _source_role_pointer(namespace, field):
    """Return the saved source-role path for one C3/C5 history port."""
    _require(namespace in ("C3", "C5") and field in D_FIELDS,
             "invalid source-role pointer request")
    if namespace == "C3":
        return (
            "source_recipe_contract.ordinary_Z_roles.C3."
            f"terminal_history_ports_Z0_to_Z3.{field}"
        )
    return (
        "source_recipe_contract.ordinary_Z_roles.C5."
        f"terminal_scaled_history_ports_Z0_to_Z3.{field}"
    )

PRIOR_FAILED_FAMILY = {
    "result_path": "tmp/compact_joint_history_gap_binding_v1/result.json",
    "result_sha256": "50257d1c62728ba554228f69e68d50d4c4be673ff0b0bf13e79b3c6b75b3be3c",
    "status": "failed",
    "exception": "KeyError: 0",
    "reentered": False,
    "proof_body_reused": False,
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _raw_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _raw_archive_receipts(root):
    """Authenticate gzip bytes before any possible decode; no decode is done."""
    receipts = {}
    for namespace in ("C3", "C5"):
        relative = ARCHIVE_PATHS[namespace]
        path = Path(root) / relative
        _require(path.is_file(), f"missing {namespace} raw archive")
        digest = _raw_sha256(path)
        _require(digest == ARCHIVE_SHA256[namespace],
                 f"{namespace} raw archive SHA changed")
        receipts[namespace] = {
            "relative_path": str(relative).replace("\\", "/"),
            "raw_sha256": digest,
            "raw_checked_before_gzip_decode": True,
            "gzip_decoded": False,
        }
    return receipts


def _load_report(root, relative, expected_sha256):
    path = Path(root) / relative
    _require(path.is_file(), f"missing saved report: {relative}")
    _require(_raw_sha256(path) == expected_sha256,
             f"saved report SHA changed: {relative}")
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"saved report is not an object: {relative}")
    return value


def _audit_non_numeric_saved_state(signed, joint):
    for label, report in (
        ("signed", signed),
        ("signed.full", signed["full_result"]),
        ("signed.component", signed["component_result"]),
    ):
        _require(report["source_only"] is True, f"{label} is not source-only")
        _require(report["original_source_evaluations"] == 0,
                 f"{label} evaluated an original source")
        _require(report["native_N_scaling_applied_once_bound"] is True,
                 f"{label} lost one-time N scaling")
    for label, report in (
        ("joint", joint),
        ("joint.full", joint["full_result"]),
        ("joint.component", joint["component_result"]),
    ):
        _require(report["source_only"] is True, f"{label} is not source-only")
        _require(report["signed_target_gaps_evaluated"] is False and
                 report["numeric_values_computed"] is False,
                 f"{label} contains evaluated gaps or values")
        _require(report["actual_C3_C5_incoming_source_identity_proved"] is False,
                 f"{label} upgrades incoming source equality")
        binding = report["source_binding"]
        _require(binding["actual_target_gap_zero_proved"] is False and
                 binding["physical_Pstar_lift_applied"] is False,
                 f"{label} upgrades a nonclaim")
        _require(report["cap_to_source_identity_proved"] is False,
                 f"{label} upgrades cap transfer")


def _audit_old_data_and_offset(signed, joint):
    reports = (
        signed,
        signed["full_result"],
        signed["component_result"],
        joint,
        joint["full_result"],
        joint["component_result"],
    )
    for report in reports:
        _require(report["retained_construction_data"] == EXPECTED_OLD_DATA,
                 "old 212/213/183 inventory changed")
        _require(report["retained_integral_term_inventory"] == EXPECTED_INTEGRAL_INVENTORY,
                 "retained integral inventory changed")
        for flag in OLD_ERROR_FLAGS:
            _require(report[flag] is False, f"old error flag changed: {flag}")
        provenance = report["retained_cap_provenance"]
        _require(provenance["all_original_caps_retained"] is True and
                 provenance["cap_to_source_bridge_proved"] is False and
                 provenance["numeric_centers_used_as_source_values"] is False,
                 "retained cap provenance changed")

    markers = []
    for index, report in enumerate((joint, joint["full_result"], joint["component_result"])):
        scope = report["remaining_scope"]
        _require(isinstance(scope, list) and len(scope) > 3,
                 "joint remaining scope is not indexable at [3]")
        marker = scope[3]
        _require(isinstance(marker, str) and REFERENCE_OFFSET_MARKER in marker,
                 "independent reference-offset marker is not at $.remaining_scope[3]")
        markers.append({
            "report_index": index,
            "pointer": "$.remaining_scope[3]",
            "marker": marker,
            "inherited_symbolic": True,
        })
    return {
        "old_data_retained": EXPECTED_OLD_DATA,
        "integral_inventory_retained": EXPECTED_INTEGRAL_INVENTORY,
        "old_error_flags": {flag: False for flag in OLD_ERROR_FLAGS},
        "reference_offset": {
            "status": "inherited_symbolic_marker_only",
            "inherited_symbolic": True,
            "marker_pointer": "$.remaining_scope[3]",
            "numeric_offset_selected": False,
            "markers": markers,
        },
    }


def _read_common_N_descriptor(signed):
    terminal = signed["full_result"]["terminal_history_common_N_binding"]
    literal = signed["full_result"]["source_binding"]["common_N_literal_binding"]
    _require(terminal["N_node"] == 18 and terminal["invN_node"] == 19 and
             terminal["J"] == 1358356628656378313 and
             terminal["N_integer_materialized"] is False,
             "N handle or materialization state changed")
    _require(terminal["C3_archive_sha256"] == ARCHIVE_SHA256["C3"] and
             terminal["native_archive_sha256"] == ARCHIVE_SHA256["C5"],
             "N archive namespaces changed")
    _require(literal["C3_node_18"] == 18 and literal["native_node_18"] == 18 and
             literal["C3_invN_node_19"] == 19 and literal["native_invN_node_19"] == 19 and
             literal["N_integer_materialized"] is False,
             "common N literal handles changed")
    return {
        "N_node": 18,
        "invN_node": 19,
        "J": 1358356628656378313,
        "definition": "N=2^(2^J)",
        "type": "DoubleDyadicInteger",
        "N_Z_independent": True,
        "N_integer_materialized": False,
        "source_pointers": {
            "terminal": "$.full_result.terminal_history_common_N_binding",
            "literal": "$.full_result.source_binding.common_N_literal_binding",
        },
    }


def _read_D_roles(joint):
    roles = joint["source_recipe_contract"]["ordinary_Z_roles"]
    _require(roles["orders"] == [0, 1, 2, 3] and
             roles["orders_are_source_roles"] is True and
             roles["range_values_are_not_function_values"] is True,
             "ordinary-Z role contract changed")
    c3 = roles["C3"]
    c5 = roles["C5"]
    _require(c3["namespace"] == "repair_C3_band" and
             c5["namespace"] == "true_compact_C5",
             "D source namespaces changed")
    _require(tuple(c3["amplitude_ports_Z0_to_Z3"]) == EXPECTED_C3["amplitude"] and
             tuple(c5["amplitude_ports_Z0_to_Z3"]) == EXPECTED_C5["amplitude"] and
             tuple(c3["terminal_history_ports_Z0_to_Z3"]["k"]) == EXPECTED_C3["k"] and
             tuple(c3["terminal_history_ports_Z0_to_Z3"]["m"]) == EXPECTED_C3["m"] and
             tuple(c5["terminal_scaled_history_ports_Z0_to_Z3"]["k"]) == EXPECTED_C5["k"] and
             tuple(c5["terminal_scaled_history_ports_Z0_to_Z3"]["m"]) == EXPECTED_C5["m"],
             "D source ports changed")

    manifest = joint["joint_target_numerator_certificate"]["source_gap_manifest"]
    _require(isinstance(manifest, list) and len(manifest) == 2,
             "joint source manifest cardinality changed")
    manifest_by_namespace = {}
    for index, entry in enumerate(manifest):
        namespace = entry["namespace"]
        _require(namespace in ("C3", "C5") and namespace not in manifest_by_namespace,
                 "joint source manifest namespace changed")
        _require(entry["archive_raw_sha256"] == ARCHIVE_SHA256[namespace],
                 "joint source manifest archive changed")
        history_key = ("history_source_nodes_Z0_to_Z3")
        expected_history = (c3["terminal_history_ports_Z0_to_Z3"] if namespace == "C3"
                            else c5["terminal_scaled_history_ports_Z0_to_Z3"])
        _require(entry["amplitude_source_nodes_Z0_to_Z3"] ==
                 (c3["amplitude_ports_Z0_to_Z3"] if namespace == "C3"
                  else c5["amplitude_ports_Z0_to_Z3"]) and
                 tuple(entry[history_key]["k"]) == tuple(expected_history["k"]) and
                 tuple(entry[history_key]["m"]) == tuple(expected_history["m"]),
                 "D manifest ports do not match D source roles")
        manifest_by_namespace[namespace] = {"index": index, "entry": entry}
    _require(set(manifest_by_namespace) == {"C3", "C5"},
             "D manifest lacks C3 or C5")
    return roles, manifest_by_namespace


def _history_source_witness(
    namespace, field, order, source_node, manifest, source_role,
    normalization_pointer,
):
    _require(namespace in ("C3", "C5") and field in D_FIELDS and order in ORDERS,
             "invalid D history witness")
    manifest_index = manifest[namespace]["index"]
    return {
        "namespace": namespace,
        "archive_raw_sha256": ARCHIVE_SHA256[namespace],
        "field": field,
        "ordinary_Z_order": order,
        "source_node": source_node,
        "source_role": source_role,
        "normalization_pointer": normalization_pointer,
        "manifest_pointer": (
            "$.joint_target_numerator_certificate.source_gap_manifest["
            f"{manifest_index}].history_source_nodes_Z0_to_Z3.{field}[{order}]"
        ),
    }


def _bind_signed_D_history(signed, roles, manifest):
    rows = signed["full_result"]["signed_history_gap_function_nodes"]
    _require(isinstance(rows, list) and len(rows) == 20,
             "saved signed history does not retain 20 rows")
    output = {}
    row_offsets = {"m": 0, "k": 8}
    for field in D_FIELDS:
        output[field] = {}
        for order in ORDERS:
            row_index = row_offsets[field] + order
            row = rows[row_index]
            left = row["left"]
            right = row["right"]
            _require(row["gap_index"] == row_index and row["field"] == field and
                     row["ordinary_Z_order"] == order and
                     row["operation"] == "signed_difference",
                     "D signed-history row order changed")
            _require(left["namespace"] == "true_compact_C5" and
                     left["archive_sha256"] == ARCHIVE_SHA256["C5"] and
                     left["source_node"] == EXPECTED_C5[field][order] and
                     left["raw_correction_node"] == EXPECTED_C5_RAW[field][order] and
                     left["raw_correction_role"] == "native_function_substitution_at_coordinate_2" and
                     left["N_factor"] == 18 and
                     left["N_scaling_applied_once"] is True and
                     left["role"] == "native_N_scaled_terminal_history",
                     "C5 D signed-history source changed")
            _require(right["namespace"] == "repair_C3_band" and
                     right["archive_sha256"] == ARCHIVE_SHA256["C3"] and
                     right["source_node"] == EXPECTED_C3[field][order] and
                     right["already_N_scaled"] is True and
                     right["role"] == "C3_rebuilt_terminal_N_scaled_history" and
                     right.get("source_role") == _source_role_pointer("C3", field),
                     "C3 D signed-history source changed")
            _require(
                roles["C5"]["terminal_scaled_history_ports_Z0_to_Z3"][field][order] ==
                left["source_node"] and
                roles["C5"]["raw_correction_ports_Z0_to_Z3"][field][order] ==
                left["raw_correction_node"] and
                roles["C3"]["terminal_history_ports_Z0_to_Z3"][field][order] ==
                right["source_node"] and
                manifest["C5"]["entry"]["namespace"] == "C5" and
                manifest["C3"]["entry"]["namespace"] == "C3" and
                manifest["C5"]["entry"]["archive_raw_sha256"] == left["archive_sha256"] and
                manifest["C3"]["entry"]["archive_raw_sha256"] == right["archive_sha256"] and
                manifest["C5"]["entry"]["history_source_nodes_Z0_to_Z3"][field][order] ==
                left["source_node"] and
                manifest["C3"]["entry"]["history_source_nodes_Z0_to_Z3"][field][order] ==
                right["source_node"],
                "signed D source node disagrees with joint source manifest",
            )
            _require(row["normalization"] == NORMALIZATION_LITERAL and
                     row["values_uncomputed"] is True and
                     row["gap_zero_claimed"] is False and
                     row["same_units_claimed"] is False and
                     row["conditional_common_point_pullback_unbound"] is True,
                     "D signed-history nonclaims changed")
            normalization_pointer = (
                "$.full_result.signed_history_gap_function_nodes["
                f"{row_index}].normalization"
            )
            c5_source_role = _source_role_pointer("C5", field)
            output[field][str(order)] = {
                "signed_row_pointer": (
                    "$.full_result.signed_history_gap_function_nodes["
                    f"{row_index}]"
                ),
                "C5_source": _history_source_witness(
                    "C5", field, order, left["source_node"], manifest,
                    c5_source_role, normalization_pointer,
                ),
                "C3_source": _history_source_witness(
                    "C3", field, order, right["source_node"], manifest,
                    right["source_role"], normalization_pointer,
                ),
                "C5_role_field": left["role"],
                "C5_raw_correction_node": left["raw_correction_node"],
                "N_factor": left["N_factor"],
                "normalization_literal": row["normalization"],
                "normalization_pointer": normalization_pointer,
                "physical_Pstar_lift_applied": False,
                "common_physical_domain_claimed": False,
                "gap_value_selected": False,
            }
    return output


def _normalization_context_pointers(signed_D, order):
    _require(order in ORDERS, "invalid normalization context order")
    return {
        "k": signed_D["k"][str(order)]["normalization_pointer"],
        "m": signed_D["m"][str(order)]["normalization_pointer"],
    }


def _d_atom_witnesses(roles, manifest, signed_D, order):
    _require(order in ORDERS, "invalid D ordinary-Z order")
    c3 = roles["C3"]
    c5 = roles["C5"]
    direct = {}
    for atom_order in ORDERS:
        c3_a = c3["amplitude_ports_Z0_to_Z3"][atom_order]
        c5_a = c5["amplitude_ports_Z0_to_Z3"][atom_order]
        direct[f"tn_C3_A_node_{c3_a}"] = {
            "namespace": "repair_C3_band",
            "archive_raw_sha256": ARCHIVE_SHA256["C3"],
            "kind": "amplitude",
            "ordinary_Z_order": atom_order,
            "source_node": c3_a,
            "normalization_context_pointers": _normalization_context_pointers(
                signed_D, atom_order
            ),
            "normalization_applies_directly": False,
            "manifest_pointer": (
                "$.joint_target_numerator_certificate.source_gap_manifest["
                f"{manifest['C3']['index']}].amplitude_source_nodes_Z0_to_Z3[{atom_order}]"
            ),
        }
        direct[f"tn_C5_A_node_{c5_a}"] = {
            "namespace": "true_compact_C5",
            "archive_raw_sha256": ARCHIVE_SHA256["C5"],
            "kind": "amplitude",
            "ordinary_Z_order": atom_order,
            "source_node": c5_a,
            "normalization_context_pointers": _normalization_context_pointers(
                signed_D, atom_order
            ),
            "normalization_applies_directly": False,
            "manifest_pointer": (
                "$.joint_target_numerator_certificate.source_gap_manifest["
                f"{manifest['C5']['index']}].amplitude_source_nodes_Z0_to_Z3[{atom_order}]"
            ),
        }
        for field in ("k", "m"):
            c3_node = c3["terminal_history_ports_Z0_to_Z3"][field][atom_order]
            c5_node = c5["terminal_scaled_history_ports_Z0_to_Z3"][field][atom_order]
            direct[f"tn_C3_{field}_node_{c3_node}"] = {
                "namespace": "repair_C3_band",
                "archive_raw_sha256": ARCHIVE_SHA256["C3"],
                "kind": "history",
                "field": field,
                "ordinary_Z_order": atom_order,
                "source_node": c3_node,
                "signed_row_pointer": signed_D[field][str(atom_order)]["signed_row_pointer"],
                "normalization_pointer": signed_D[field][str(atom_order)]["normalization_pointer"],
                "normalization_context_pointers": _normalization_context_pointers(
                    signed_D, atom_order
                ),
                "normalization_applies_directly": True,
                "manifest_pointer": (
                    "$.joint_target_numerator_certificate.source_gap_manifest["
                    f"{manifest['C3']['index']}].history_source_nodes_Z0_to_Z3.{field}[{atom_order}]"
                ),
            }
            direct[f"tn_C5_{field}_node_{c5_node}"] = {
                "namespace": "true_compact_C5",
                "archive_raw_sha256": ARCHIVE_SHA256["C5"],
                "kind": "history",
                "field": field,
                "ordinary_Z_order": atom_order,
                "source_node": c5_node,
                "signed_row_pointer": signed_D[field][str(atom_order)]["signed_row_pointer"],
                "normalization_pointer": signed_D[field][str(atom_order)]["normalization_pointer"],
                "normalization_context_pointers": _normalization_context_pointers(
                    signed_D, atom_order
                ),
                "normalization_applies_directly": True,
                "manifest_pointer": (
                    "$.joint_target_numerator_certificate.source_gap_manifest["
                    f"{manifest['C5']['index']}].history_source_nodes_Z0_to_Z3.{field}[{atom_order}]"
                ),
            }
    direct["tn_mu"] = {
        "namespace": "joint-common",
        "kind": "positive_common_parameter",
        "source_value_selected": False,
        "normalization_context_pointers": _normalization_context_pointers(
            signed_D, order
        ),
        "normalization_applies_directly": False,
        "json_pointer": (
            "$.joint_target_numerator_certificate."
            "conditional_denominator_contract.mu"
        ),
    }
    return direct


def _literal_row_tokens(row):
    _require(isinstance(row, dict) and isinstance(row.get("terms"), list),
             "D formula row terms are not a list")
    tokens = set()
    for term in row["terms"]:
        _require(isinstance(term, dict), "D formula term is not an object")
        monomial = term.get("monomial")
        _require(isinstance(monomial, list),
                 "D formula term monomial is not a list")
        tokens.update(monomial)
    return tokens


def _require_D_row_atom_roles(order, numerator, denominator):
    """Require only literal atom presence, not a polynomial identity proof."""
    numerator_tokens = _literal_row_tokens(numerator)
    denominator_tokens = _literal_row_tokens(denominator)

    required_history = {
        f"tn_C3_{field}_node_{EXPECTED_C3[field][order]}"
        for field in D_FIELDS
    }
    required_history.update(
        f"tn_C5_{field}_node_{EXPECTED_C5[field][order]}"
        for field in D_FIELDS
    )
    _require(required_history.issubset(numerator_tokens),
             "D numerator lacks a current-order C3/C5 k or m atom")

    required_amplitude = {
        f"tn_C3_A_node_{EXPECTED_C3['amplitude'][order]}",
        f"tn_C5_A_node_{EXPECTED_C5['amplitude'][order]}",
        "tn_mu",
    }
    _require(required_amplitude.issubset(denominator_tokens),
             "D denominator lacks current-order amplitude or tn_mu")

    if order > 0:
        # The saved derivative rows may also contain lower-order atoms.  They
        # must remain inside the saved source-role ports; only current-order
        # k/m and amplitude ports are mandatory here.
        allowed_history = {
            f"tn_C3_{field}_node_{node}"
            for field in D_FIELDS
            for node in EXPECTED_C3[field]
        }
        allowed_history.update(
            f"tn_C5_{field}_node_{node}"
            for field in D_FIELDS
            for node in EXPECTED_C5[field]
        )
        for token in numerator_tokens:
            if token.startswith(("tn_C3_k_node_", "tn_C3_m_node_",
                                 "tn_C5_k_node_", "tn_C5_m_node_")):
                _require(token in allowed_history,
                         "D derivative numerator has an unbound lower-order k/m atom")
    return {
        "required_numerator_current_order_k_m_atoms": sorted(required_history),
        "required_denominator_current_order_amplitude_atoms": sorted(
            required_amplitude - {"tn_mu"}
        ),
        "required_denominator_tn_mu": True,
        "order_gt_zero_lower_order_k_m_atoms_allowed": order > 0,
        "formula_identity_reproved": False,
    }


def _bind_saved_D_rows(joint, roles, manifest, signed_D):
    certificate = joint["joint_target_numerator_certificate"]
    formula = certificate["numerators_and_denominators"]["D"]
    _require(formula["raw_gap_role"] == "joint_delta_k_and_delta_m" and
             formula["joint_numerator"] is True,
             "saved D channel metadata changed")
    bound = []
    for order in ORDERS:
        numerator = formula["numerator"][order]
        denominator = formula["denominator"][order]
        _require(numerator["ordinary_Z_order"] == order and
                 denominator["ordinary_Z_order"] == order,
                 "saved D ordinary-Z row changed")
        atom_presence_contract = _require_D_row_atom_roles(
            order, numerator, denominator
        )
        direct = _d_atom_witnesses(roles, manifest, signed_D, order)
        side_rows = []
        for side, row in (("numerator", numerator), ("denominator", denominator)):
            _require(isinstance(row, dict) and isinstance(row.get("terms"), list),
                     f"saved D {side} terms are not a list")
            terms = []
            for term_index, term in enumerate(row["terms"]):
                _require(isinstance(term, dict),
                         f"saved D {side} term is not an object")
                literal_atoms = term["monomial"]
                _require(isinstance(literal_atoms, list),
                         f"saved D {side} term monomial is not a list")
                atom_bindings = []
                for atom_index, token in enumerate(literal_atoms):
                    _require(token in direct,
                             f"D {side} atom has no direct source witness: {token}")
                    atom_bindings.append({
                        "atom_index": atom_index,
                        "symbolic_atom": token,
                        "occurrence_pointer": (
                            "$.joint_target_numerator_certificate."
                            "numerators_and_denominators.D."
                            f"{side}[{order}].terms[{term_index}].monomial[{atom_index}]"
                        ),
                        "source_witness": direct[token],
                    })
                terms.append({
                    "term_index": term_index,
                    "literal_term": term,
                    "atom_bindings": atom_bindings,
                })
            side_rows.append({
                "side": side,
                "row_pointer": (
                    "$.joint_target_numerator_certificate."
                    "numerators_and_denominators.D."
                    f"{side}[{order}]"
                ),
                "terms": terms,
            })
        bound.append({
            "ordinary_Z_order": order,
            "signed_k_row_pointer": signed_D["k"][str(order)]["signed_row_pointer"],
            "signed_m_row_pointer": signed_D["m"][str(order)]["signed_row_pointer"],
            "D_source_rows": side_rows,
            "literal_atom_presence_contract": atom_presence_contract,
            "D_k_m_gap_role": "joint_delta_k_and_delta_m",
            "extra_N_divisions": 0,
            "gap_value_selected": False,
            "formula_identity_reproved": False,
        })
    return bound


def bind_D_km_and_reference_offset(repo_root):
    """Build only D/k/m and reference-offset source bindings; never called here."""
    root = Path(repo_root)
    raw_archives = _raw_archive_receipts(root)
    signed = _load_report(root, SIGNED_REPORT, REPORT_SHA256["signed"])
    joint = _load_report(root, JOINT_REPORT, REPORT_SHA256["joint"])
    _audit_non_numeric_saved_state(signed, joint)
    retained = _audit_old_data_and_offset(signed, joint)
    N = _read_common_N_descriptor(signed)
    roles, manifest = _read_D_roles(joint)
    signed_D = _bind_signed_D_history(signed, roles, manifest)
    d_rows = _bind_saved_D_rows(joint, roles, manifest, signed_D)
    return {
        "schema": SCHEMA,
        "source_only": True,
        "scope": "D-channel signed k/m source binding plus reference-offset marker only",
        "prior_failed_family": PRIOR_FAILED_FAMILY,
        "report_sha256": REPORT_SHA256,
        "raw_archive_receipts": raw_archives,
        "N_descriptor": N,
        "signed_D_history_rows": signed_D,
        "signed_D_history_row_count": 8,
        "actual_D_source_rows": d_rows,
        "reference_offset_source_binding": retained["reference_offset"],
        "retained_legacy_data": retained,
        "explicitly_unbound_history_fields": ["h", "e", "p"],
        "explicitly_unentered_channels": ["M", "I", "S", "Cp"],
        "nonclaims": {
            "signed_gap_value": False,
            "D_gap_zero": False,
            "common_physical_domain": False,
            "same_units": False,
            "physical_Pstar_lift": False,
            "cap_transfer": False,
            "old_error_reduction": False,
        },
        "audit": {
            "source_modules_imported": 0,
            "source_functions_executed": 0,
            "numeric_values_evaluated": 0,
            "formula_proof_replayed": 0,
            "M_Z0_term_family_reentered": False,
            "prior_v1_helper_called": False,
            "gzip_archives_decoded": 0,
        },
        "json_pointers": {
            "signed_k_rows": [
                "$.full_result.signed_history_gap_function_nodes[8]",
                "$.full_result.signed_history_gap_function_nodes[9]",
                "$.full_result.signed_history_gap_function_nodes[10]",
                "$.full_result.signed_history_gap_function_nodes[11]",
            ],
            "signed_m_rows": [
                "$.full_result.signed_history_gap_function_nodes[0]",
                "$.full_result.signed_history_gap_function_nodes[1]",
                "$.full_result.signed_history_gap_function_nodes[2]",
                "$.full_result.signed_history_gap_function_nodes[3]",
            ],
            "D_numerator": [
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.numerator[0]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.numerator[1]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.numerator[2]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.numerator[3]",
            ],
            "D_denominator": [
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.denominator[0]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.denominator[1]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.denominator[2]",
                "$.joint_target_numerator_certificate.numerators_and_denominators.D.denominator[3]",
            ],
            "reference_offset": "$.remaining_scope[3]",
        },
        "remaining_scope": [
            "The h/e/p signed-history rows remain unbound by this D-only packet.",
            "M/I/S/Cp numerator atoms remain intentionally unentered.",
            "The reference offset has only a saved unresolved marker; no numeric offset is selected.",
            "No common physical chart, Pstar lift, gap equality, or cap transfer is admitted.",
        ],
    }


# No __main__ entry point: this draft is intentionally unexecuted.
