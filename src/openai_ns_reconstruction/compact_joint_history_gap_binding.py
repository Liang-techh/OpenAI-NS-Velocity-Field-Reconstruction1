"""UNEXECUTED draft: bind saved signed-history ports to joint symbolic gaps.

This is a source-only literal/namespace consumer design.  It deliberately does
not import reconstruction modules, evaluate source functions, build a jet,
perform polynomial arithmetic, or read the compact volume/error packet.

If later executed by the root agent, archive authentication is raw-byte SHA-256
only; this draft never decodes a gzip archive.  The returned packet remains a
literal source-port binding and carries no common-domain, value-equality, or cap
transfer claim.

The intended first execution must be performed by the root agent after this
draft is reviewed and frozen.  Until then this file is only a scoped design
artifact under W/.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


SCHEMA = "compact_joint_history_gap_binding_source_consumer_draft_v1"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
HISTORY_FIELDS = ("m", "h", "k", "e", "p")
CHANNEL_FOR_FIELD = {
    "m": "M",
    "h": "I",
    "e": "S",
    "p": "Cp",
    "k": "D",
}
RAW_GAP_ROLE = {
    "M": "delta_m",
    "I": "delta_h",
    "S": "delta_e",
    "Cp": "delta_p",
    "D": "joint_delta_k_and_delta_m",
}
JOINT_HISTORY_FIELDS = {
    "M": ("m",),
    "I": ("h",),
    "S": ("e",),
    "Cp": ("p",),
    "D": ("k", "m"),
}
JOINT_NUMERATOR_FORMULAS = {
    "M": "C*delta_m-a*S",
    "I": "C*delta_h-a*H",
    "S": "C^2*delta_e-a*(2*C+a)*E",
    "Cp": "C^2*delta_p-a*(2*C+a)*P",
    "D": "C^2*delta_k-C^2*A*delta_m-a*(2*C+a)*K+a*C*A*S",
}

COMMON_N_LITERAL_DESCRIPTOR = {
    "N": {
        "node_id": 18,
        "operation": "exact_positive_integer_expression",
        "definition": "N=2^(2^J)",
        "type": "DoubleDyadicInteger",
        "J": 1358356628656378313,
        "Z_independent": True,
        "outer_integer_and_inner_power_not_materialized": True,
    },
    "invN": {
        "node_id": 19,
        "operation": "positive_quotient",
        "numerator": 1,
        "denominator_node_id": 18,
        "Z_independent": True,
    },
    "endpoint": {
        "node_id": 48,
        "operation": "exact_rational",
        "numerator": 2,
        "denominator": 1,
    },
    "source_guard": (
        "src/openai_ns_reconstruction/"
        "original_compact_signed_history_defect_transport.py:shd_c3_binding"
    ),
}

# These are the authenticated raw archive digests used by both saved packets.
C3_ARCHIVE_RAW_SHA256 = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
C5_ARCHIVE_RAW_SHA256 = (
    "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
)
SIGNED_REPORT_SHA256 = (
    "dc367908bd44991ef6f4a973d688b4bd27bb0cab43464e0e396b3abb1201c861"
)
JOINT_REPORT_SHA256 = (
    "cff727e06c1af18eebeba0bbf027a488128c2021f9feb839bafdfc368b18f133"
)

SIGNED_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "signed_history_defect_transport.json"
)
JOINT_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "joint_target_numerator_transport.json"
)
ARCHIVE_RELATIVE = {
    "C3": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz"
    ),
    "C5": Path(
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
    ),
}

SIGNED_NORMALIZATION_LITERAL = (
    "generic common E=Utheta/Pstar, V=Uz/Pstar, "
    "m/k scaled1/Pstar, h/e/p retained per original moment definitions"
)
REFERENCE_OFFSET_MARKER = "independent old reference offset"
OLD_ERROR_FLAG_NAMES = (
    "old212_adapter_identity",
    "old212_quantclosure",
    "old_errors_removed",
    "old_error_terms_removed",
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
EXPECTED_RETAINED_INTEGRAL_TERM_INVENTORY = {
    "J_terms": 16,
    "J_terms_retained": True,
    "K_terms": 17,
    "K_terms_retained": True,
    "outer_B_terms": 2,
    "outer_B_terms_retained": True,
    "reserved_contribution_registered_once": True,
}


EXPECTED_C5_SCALED_PORTS = {
    "m": (258883, 258884, 258885, 258886),
    "h": (258889, 258890, 258891, 258892),
    "k": (258895, 258896, 258897, 258898),
    "e": (258901, 258902, 258903, 258904),
    "p": (258907, 258908, 258909, 258910),
}
EXPECTED_C5_RAW_PORTS = {
    "m": (258646, 258647, 258648, 258649),
    "h": (258652, 258653, 258654, 258655),
    "k": (258658, 258659, 258660, 258661),
    "e": (258664, 258665, 258666, 258667),
    "p": (258670, 258671, 258672, 258673),
}
EXPECTED_C3_TERMINAL_PORTS = {
    "m": (4920, 4921, 12042, 24177),
    "h": (4928, 4929, 12043, 24178),
    "k": (4930, 4931, 12044, 24179),
    "e": (4932, 4933, 12045, 24180),
    "p": (4888, 4889, 12041, 24176),
}
EXPECTED_C5_AMPLITUDE_PORTS = (251009, 251010, 251011, 251012)
EXPECTED_C3_AMPLITUDE_PORTS = (4934, 4935, 12046, 24181)

MANIFEST_ALIAS = {
    "C3": {
        "manifest_namespace": "C3",
        "contract_namespace": "repair_C3_band",
        "archive_sha256": C3_ARCHIVE_RAW_SHA256,
    },
    "C5": {
        "manifest_namespace": "C5",
        "contract_namespace": "true_compact_C5",
        "archive_sha256": C5_ARCHIVE_RAW_SHA256,
    },
}

REPORT_POINTERS = {
    "signed_rows_full": "$.full_result.signed_history_gap_function_nodes",
    "signed_rows_component": "$.component_result.signed_history_gap_function_nodes",
    "signed_source_binding_top": "$.source_binding",
    "signed_source_binding_full": "$.full_result.source_binding",
    "signed_source_binding_component": "$.component_result.source_binding",
    "signed_common_N_full": "$.full_result.terminal_history_common_N_binding",
    "signed_common_N_component": "$.component_result.terminal_history_common_N_binding",
    "joint_contract_top": "$.source_recipe_contract",
    "joint_contract_full": "$.full_result.source_recipe_contract",
    "joint_contract_component": "$.component_result.source_recipe_contract",
    "joint_certificate_top": "$.joint_target_numerator_certificate",
    "joint_certificate_full": "$.full_result.joint_target_numerator_certificate",
    "joint_certificate_component": "$.component_result.joint_target_numerator_certificate",
    "joint_manifest": "$.joint_target_numerator_certificate.source_gap_manifest",
    "joint_formulas": "$.joint_target_numerator_certificate.numerators_and_denominators",
    "archive_C3_raw": str(ARCHIVE_RELATIVE["C3"]).replace("\\", "/"),
    "archive_C5_raw": str(ARCHIVE_RELATIVE["C5"]).replace("\\", "/"),
    "signed_retained_construction_data": "$.retained_construction_data",
    "joint_retained_construction_data": "$.retained_construction_data",
    "joint_retained_cap_provenance": "$.retained_cap_provenance",
    "joint_reference_offset_scope": "$.remaining_scope[*]",
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_authenticated_json(path, expected_sha256):
    _require(_sha256(path) == expected_sha256, "saved report SHA changed")
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def _check_authenticated_archives(repo_root):
    """Authenticate gzip files by raw bytes without decoding either archive."""
    checked = {}
    expected = {
        "C3": C3_ARCHIVE_RAW_SHA256,
        "C5": C5_ARCHIVE_RAW_SHA256,
    }
    for namespace, relative in ARCHIVE_RELATIVE.items():
        path = repo_root / relative
        _require(path.is_file(), f"missing {namespace} raw archive")
        raw_sha256 = _sha256(path)
        _require(raw_sha256 == expected[namespace],
                 f"{namespace} raw archive SHA changed before gzip decode")
        checked[namespace] = {
            "relative_path": str(relative).replace("\\", "/"),
            "raw_sha256": raw_sha256,
            "raw_sha256_checked_before_gzip_decode": True,
            "gzip_decoded": False,
        }
    return checked


def _same_json(left, right):
    return json.dumps(left, sort_keys=True, separators=(",", ":")) == json.dumps(
        right, sort_keys=True, separators=(",", ":")
    )


def _signed_slices(signed_report):
    full = signed_report["full_result"]
    component = signed_report["component_result"]
    _require(
        _same_json(full["signed_history_gap_function_nodes"],
                   component["signed_history_gap_function_nodes"]),
        "full/component signed-history row pointers diverged",
    )
    _require(
        _same_json(signed_report["source_binding"], full["source_binding"]) and
        _same_json(full["source_binding"], component["source_binding"]),
        "signed source-binding pointers diverged",
    )
    _require(
        _same_json(full["terminal_history_common_N_binding"],
                   component["terminal_history_common_N_binding"]),
        "full/component common-N pointers diverged",
    )
    return full, component


def _joint_slices(joint_report):
    full = joint_report["full_result"]
    component = joint_report["component_result"]
    _require(
        _same_json(joint_report["source_recipe_contract"],
                   full["source_recipe_contract"]) and
        _same_json(full["source_recipe_contract"],
                   component["source_recipe_contract"]),
        "joint source-recipe pointers diverged",
    )
    _require(
        _same_json(joint_report["joint_target_numerator_certificate"],
                   full["joint_target_numerator_certificate"]) and
        _same_json(full["joint_target_numerator_certificate"],
                   component["joint_target_numerator_certificate"]),
        "joint certificate pointers diverged",
    )
    return full, component


def _ordinary_roles(joint_report):
    _joint_slices(joint_report)
    contract = joint_report["source_recipe_contract"]
    roles = contract["ordinary_Z_roles"]
    _require(roles["orders"] == list(ORDINARY_Z_ORDERS), "ordinary Z roles changed")
    _require(roles["orders_are_source_roles"] is True, "ordinary Z role flag changed")
    _require(roles["range_values_are_not_function_values"] is True,
             "range values were admitted as functions")
    return roles["C3"], roles["C5"]


def _check_saved_contract_flags(signed_report, joint_report):
    signed_full, signed_component = _signed_slices(signed_report)
    joint_full, joint_component = _joint_slices(joint_report)
    for report in (signed_report, signed_full, signed_component):
        _require(report["source_only"] is True, "signed report is not source-only")
        _require(report["original_source_evaluations"] == 0,
                 "signed report evaluated a source")
        _require(report["signed_history_source_gap_functions_defined"] is True,
                 "signed history packet is not defined")
        _require(report["native_N_scaling_applied_once_bound"] is True,
                 "native N scaling was not authenticated")
    for report in (joint_report, joint_full, joint_component):
        _require(report["source_only"] is True, "joint report is not source-only")
        _require(report["signed_target_gaps_evaluated"] is False,
                 "joint signed gaps were evaluated")
        _require(report["numeric_values_computed"] is False,
                 "joint numeric values were computed")
        _require(report["joint_target_numerator_source_defined"] is True,
                 "joint numerator source packet is not defined")
        _require(report["five_joint_numerator_identities_proved"] is True,
                 "joint numerator identities are unavailable")
        _require(report["actual_C3_C5_incoming_source_identity_proved"] is False,
                 "would silently upgrade an existing false incoming identity")
        _require(report["source_binding"]["actual_target_gap_zero_proved"] is False,
                 "would silently claim a zero target gap")
        _require(report["source_binding"]["physical_Pstar_lift_applied"] is False,
                 "would silently claim a physical Pstar lift")
        _require(report["cap_to_source_identity_proved"] is False,
                 "would silently claim cap transfer")


def _check_common_N_descriptor(signed_report):
    signed_full, signed_component = _signed_slices(signed_report)
    expected_row = {
        "C3_invN_node_19": 19,
        "C3_node_18": 18,
        "J": COMMON_N_LITERAL_DESCRIPTOR["N"]["J"],
        "N_integer_materialized": False,
        "native_invN_node_19": 19,
        "native_node_18": 18,
    }
    expected_terminal = {
        "C3_archive_sha256": C3_ARCHIVE_RAW_SHA256,
        "J": COMMON_N_LITERAL_DESCRIPTOR["N"]["J"],
        "N_integer_materialized": False,
        "N_node": 18,
        "endpoint_literal_rows_bound": {
            "C3": {
                "node_id": 48,
                "row": {
                    "denominator": 1,
                    "numerator": 2,
                    "operation": "exact_rational",
                },
            },
            "native": {
                "node_id": 48,
                "row": {
                    "denominator": 1,
                    "numerator": 2,
                    "operation": "exact_rational",
                },
            },
        },
        "inlet_cap_equality_claimed": False,
        "invN_node": 19,
        "native_archive_sha256": C5_ARCHIVE_RAW_SHA256,
        "same_units_claimed": False,
    }
    for report in (signed_report, signed_full, signed_component):
        _require(report["source_binding"]["common_N_literal_binding"] == expected_row,
                 "common-N source-binding row changed")
    for report in (signed_full, signed_component):
        _require(report["terminal_history_common_N_binding"] == expected_terminal,
                 "common-N terminal descriptor changed")
        for row in report["signed_history_gap_function_nodes"]:
            _require(row["common_N_literal_binding"] == expected_row,
                     "gap row common-N descriptor changed")
    return {
        "literal_descriptor": COMMON_N_LITERAL_DESCRIPTOR,
        "report_row_descriptor": expected_row,
        "terminal_descriptor": expected_terminal,
        "report_pointers": {
            "full": REPORT_POINTERS["signed_common_N_full"],
            "component": REPORT_POINTERS["signed_common_N_component"],
        },
    }


def _signed_rows_by_key(signed_report):
    signed_full, signed_component = _signed_slices(signed_report)
    rows = signed_full["signed_history_gap_function_nodes"]
    _require(len(rows) == 20, "signed-history row count changed")
    by_key = {}
    for expected_index, row in enumerate(rows):
        _require(row["gap_index"] == expected_index,
                 "signed-history gap_index is not 0..19 in source order")
        key = (row["field"], row["ordinary_Z_order"])
        _require(key not in by_key, "duplicate signed-history field/order")
        _require(row["ordinary_Z_order"] in ORDINARY_Z_ORDERS,
                 "unexpected ordinary Z order")
        by_key[key] = row
    _require(set(by_key) == {
        (field, order)
        for field in HISTORY_FIELDS
        for order in ORDINARY_Z_ORDERS
    }, "signed-history field/order grid changed")
    return by_key


def _check_joint_gap_definitions(joint_report):
    expected = {
        "delta_m": "T_N_m-T_C_m",
        "delta_h": "T_N_h-T_C_h",
        "delta_k": "T_N_k-T_C_k",
        "delta_e": "T_N_e-T_C_e",
        "delta_p": "T_N_p-T_C_p",
        "a": "A_N-A_C",
    }
    _joint_slices(joint_report)
    top = joint_report["joint_target_numerator_certificate"]["signed_gap_definitions"]
    top_cert = joint_report["joint_target_numerator_certificate"]
    full_cert = joint_report["full_result"]["joint_target_numerator_certificate"]
    component_cert = joint_report["component_result"]["joint_target_numerator_certificate"]
    full = full_cert["signed_gap_definitions"]
    component = component_cert["signed_gap_definitions"]
    _require(top == expected and full == expected and component == expected,
             "joint symbolic gap definitions changed")
    for certificate in (top_cert, full_cert, component_cert):
        _require(certificate["source_equality_admitted"] is False and
                 certificate["gap_zero_assumed"] is False and
                 certificate["cap_admitted"] is False,
                 "joint certificate nonclaims changed")


def _check_joint_manifest(joint_report, c3_roles, c5_roles):
    _joint_slices(joint_report)
    manifest = joint_report["joint_target_numerator_certificate"]["source_gap_manifest"]
    _require(len(manifest) == 2, "joint source-gap manifest cardinality changed")
    _require(
        [(entry["namespace"], entry["archive_raw_sha256"]) for entry in manifest]
        == [("C3", C3_ARCHIVE_RAW_SHA256), ("C5", C5_ARCHIVE_RAW_SHA256)],
        "joint C3/C5 manifest namespace or archive mapping changed",
    )
    expected = {
        "C3": (c3_roles["amplitude_ports_Z0_to_Z3"],
                c3_roles["terminal_history_ports_Z0_to_Z3"]),
        "C5": (c5_roles["amplitude_ports_Z0_to_Z3"],
                c5_roles["terminal_scaled_history_ports_Z0_to_Z3"]),
    }
    by_namespace = {}
    for index, entry in enumerate(manifest):
        namespace = entry["namespace"]
        amplitude, histories = expected[namespace]
        _require(tuple(entry["amplitude_source_nodes_Z0_to_Z3"]) == tuple(amplitude),
                 "joint manifest amplitude ports changed")
        _require(
            {key: tuple(value) for key, value in entry["history_source_nodes_Z0_to_Z3"].items()}
            == {key: tuple(value) for key, value in histories.items()},
            "joint manifest history ports changed",
        )
        alias = MANIFEST_ALIAS[namespace]
        _require(alias["archive_sha256"] == entry["archive_raw_sha256"],
                 "joint manifest archive alias changed")
        by_namespace[namespace] = {
            "manifest_index": index,
            "manifest_namespace": alias["manifest_namespace"],
            "contract_namespace": alias["contract_namespace"],
            "archive_raw_sha256": entry["archive_raw_sha256"],
            "entry": entry,
        }
    _require(c3_roles["namespace"] == MANIFEST_ALIAS["C3"]["contract_namespace"],
             "C3 contract namespace alias changed")
    _require(c5_roles["namespace"] == MANIFEST_ALIAS["C5"]["contract_namespace"],
             "C5 contract namespace alias changed")
    return by_namespace


def _check_retained_legacy_provenance(signed_report, joint_report):
    """Keep the old 212/error inventory and the unresolved reference marker."""
    signed_full, signed_component = _signed_slices(signed_report)
    joint_full, joint_component = _joint_slices(joint_report)
    signed_reports = (signed_report, signed_full, signed_component)
    joint_reports = (joint_report, joint_full, joint_component)
    for label, reports in (("signed", signed_reports), ("joint", joint_reports)):
        for index, report in enumerate(reports):
            _require(report["retained_construction_data"] ==
                     EXPECTED_RETAINED_CONSTRUCTION_DATA,
                     f"{label}[{index}] old 212 construction inventory changed")
            _require(report["retained_integral_term_inventory"] ==
                     EXPECTED_RETAINED_INTEGRAL_TERM_INVENTORY,
                     f"{label}[{index}] retained integral inventory changed")
            for flag in OLD_ERROR_FLAG_NAMES:
                _require(report[flag] is False,
                         f"{label}[{index}] old-error flag changed: {flag}")
            cap = report["retained_cap_provenance"]
            _require(cap["all_original_caps_retained"] is True and
                     cap["cap_to_source_bridge_proved"] is False and
                     cap["numeric_centers_used_as_source_values"] is False and
                     cap["same_family_is_not_a_function_identity"] is True,
                     f"{label}[{index}] retained cap provenance changed")

    marker_reports = []
    for index, report in enumerate(joint_reports):
        scope = report["remaining_scope"]
        _require(isinstance(scope, list),
                 f"joint[{index}] remaining scope is not a list")
        markers = [item for item in scope
                   if isinstance(item, str) and REFERENCE_OFFSET_MARKER in item]
        _require(markers,
                 f"joint[{index}] independent reference-offset marker missing")
        marker_reports.append({
            "report_index": index,
            "pointer": "$.remaining_scope[*]",
            "marker": markers[0],
        })
    return {
        "old_error_flags": {flag: False for flag in OLD_ERROR_FLAG_NAMES},
        "retained_construction_data": EXPECTED_RETAINED_CONSTRUCTION_DATA,
        "retained_integral_term_inventory": EXPECTED_RETAINED_INTEGRAL_TERM_INVENTORY,
        "retained_cap_provenance": {
            "all_original_caps_retained": True,
            "cap_to_source_bridge_proved": False,
            "numeric_centers_used_as_source_values": False,
            "same_family_is_not_a_function_identity": True,
        },
        "reference_offset": {
            "status": "unresolved_scope_marker_retained",
            "numeric_offset_selected": False,
            "common_physical_domain_claimed": False,
            "markers": marker_reports,
        },
    }


def _manifest_witness(manifest_by_namespace, namespace, kind, field, order, node):
    alias = manifest_by_namespace[namespace]
    manifest_index = alias["manifest_index"]
    if kind == "amplitude":
        pointer = (
            f"$.joint_target_numerator_certificate.source_gap_manifest[{manifest_index}]"
            f".amplitude_source_nodes_Z0_to_Z3[{order}]"
        )
    else:
        pointer = (
            f"$.joint_target_numerator_certificate.source_gap_manifest[{manifest_index}]"
            f".history_source_nodes_Z0_to_Z3.{field}[{order}]"
        )
    return {
        "manifest_namespace": alias["manifest_namespace"],
        "contract_namespace": alias["contract_namespace"],
        "archive_raw_sha256": alias["archive_raw_sha256"],
        "kind": kind,
        "field": field,
        "ordinary_Z_order": order,
        "source_node": node,
        "manifest_pointer": pointer,
    }


def _joint_atom_catalog(manifest_by_namespace, c3_roles, c5_roles):
    catalog = {}
    for namespace, roles, token_namespace in (
        ("C3", c3_roles, "C3"),
        ("C5", c5_roles, "C5"),
    ):
        amplitude = roles["amplitude_ports_Z0_to_Z3"]
        for order, node in enumerate(amplitude):
            token = f"tn_{token_namespace}_A_node_{node}"
            catalog[token] = _manifest_witness(
                manifest_by_namespace, namespace, "amplitude", "A", order, node
            )
        histories = (
            roles["terminal_history_ports_Z0_to_Z3"] if namespace == "C3"
            else roles["terminal_scaled_history_ports_Z0_to_Z3"]
        )
        for field, nodes in histories.items():
            for order, node in enumerate(nodes):
                token = f"tn_{token_namespace}_{field}_node_{node}"
                catalog[token] = _manifest_witness(
                    manifest_by_namespace, namespace, "history", field, order, node
                )
    catalog["tn_mu"] = {
        "kind": "common_positive_parameter",
        "literal_pointer": (
            "$.joint_target_numerator_certificate."
            "conditional_denominator_contract.mu"
        ),
        "source_function_value_selected": False,
    }
    return catalog


def _formula_atom_rows(formula, channel, order, catalog, c3_roles, c5_roles):
    _require(formula["raw_gap_role"] == RAW_GAP_ROLE[channel],
             f"{channel} raw gap role changed")
    _require(formula["joint_numerator"] is (channel == "D"),
             f"{channel} joint-numerator flag changed")
    numerator = formula["numerator"][order]
    denominator = formula["denominator"][order]
    _require(numerator["ordinary_Z_order"] == order and
             denominator["ordinary_Z_order"] == order,
             f"{channel} ordinary-Z formula row changed")
    atom_witnesses = []
    all_terms = list(numerator["terms"]) + list(denominator["terms"])
    seen = set()
    for term in all_terms:
        _require(isinstance(term["monomial"], list),
                 f"{channel} monomial is not a literal atom list")
        for token in term["monomial"]:
            _require(token in catalog,
                     f"{channel} contains an unbound symbolic atom: {token}")
            if token not in seen:
                atom_witnesses.append({
                    "symbolic_atom_literal": token,
                    "source_witness": catalog[token],
                })
                seen.add(token)
    fields = JOINT_HISTORY_FIELDS[channel]
    for field in fields:
        for token_namespace, roles in (("C3", c3_roles), ("C5", c5_roles)):
            nodes = (
                roles["terminal_history_ports_Z0_to_Z3"]
                if token_namespace == "C3"
                else roles["terminal_scaled_history_ports_Z0_to_Z3"]
            )
            token = f"tn_{token_namespace}_{field}_node_{nodes[order]}"
            _require(token in seen, f"{channel} omitted current {field} atom")
    for token_namespace, roles in (("C3", c3_roles), ("C5", c5_roles)):
        token = f"tn_{token_namespace}_A_node_{roles['amplitude_ports_Z0_to_Z3'][order]}"
        _require(token in seen, f"{channel} omitted current amplitude atom")
    if channel == "D":
        _require("tn_mu" in seen, "D denominator omitted tn_mu")
    return {
        "numerator_literal": numerator,
        "denominator_literal": denominator,
        "atom_witnesses": atom_witnesses,
        "json_pointers": {
            "numerator": (
                f"$.joint_target_numerator_certificate.numerators_and_denominators."
                f"{channel}.numerator[{order}]"
            ),
            "denominator": (
                f"$.joint_target_numerator_certificate.numerators_and_denominators."
                f"{channel}.denominator[{order}]"
            ),
        },
    }


def _bind_joint_formula_graph(joint_report, manifest_by_namespace, c3_roles, c5_roles):
    _joint_slices(joint_report)
    cert = joint_report["joint_target_numerator_certificate"]
    formulas = cert["numerators_and_denominators"]
    catalog = _joint_atom_catalog(manifest_by_namespace, c3_roles, c5_roles)
    graph = {}
    for channel in ("M", "I", "S", "Cp", "D"):
        graph[channel] = {
            str(order): _formula_atom_rows(
                formulas[channel], channel, order, catalog, c3_roles, c5_roles
            )
            for order in ORDINARY_Z_ORDERS
        }
    return {
        "channels": graph,
        "atom_catalog": catalog,
        "D_actual_graph_pointer": (
            "$.joint_target_numerator_certificate.numerators_and_denominators.D"
        ),
    }


def _check_gap_row(row, field, order, c5_node, c3_node):
    _require(row["field"] == field and row["ordinary_Z_order"] == order,
             "signed row key changed")
    _require(row["gap_index"] == HISTORY_FIELDS.index(field) * 4 + order,
             "signed row gap_index changed")
    _require(row["operation"] == "signed_difference", "signed operation changed")
    left, right = row["left"], row["right"]
    _require(left["namespace"] == "true_compact_C5" and
             left["archive_sha256"] == C5_ARCHIVE_RAW_SHA256 and
             left["role"] == "native_N_scaled_terminal_history" and
             left["source_node"] == c5_node and
             left["raw_correction_node"] == EXPECTED_C5_RAW_PORTS[field][order] and
             left["raw_correction_role"] ==
             "native_function_substitution_at_coordinate_2" and
             left["N_factor"] == 18 and
             left["N_scaling_applied_once"] is True,
             "C5 signed-history source identity changed")
    _require(right["namespace"] == "repair_C3_band" and
             right["archive_sha256"] == C3_ARCHIVE_RAW_SHA256 and
             right["role"] == "C3_rebuilt_terminal_N_scaled_history" and
             right["source_node"] == c3_node and
             right["source_role"] ==
             f"source_recipe_contract.ordinary_Z_roles.C3."
             f"terminal_history_ports_Z0_to_Z3.{field}" and
             right["already_N_scaled"] is True,
             "C3 signed-history source identity changed")
    _require(row["common_N_literal_binding"] == {
        "C3_invN_node_19": 19,
        "C3_node_18": 18,
        "J": COMMON_N_LITERAL_DESCRIPTOR["N"]["J"],
        "N_integer_materialized": False,
        "native_invN_node_19": 19,
        "native_node_18": 18,
    }, "gap row common-N literal changed")
    _require(row["values_uncomputed"] is True and
             row["gap_zero_claimed"] is False and
             row["same_units_claimed"] is False,
             "gap scope or unit claim changed")
    _require(row["normalization"] == SIGNED_NORMALIZATION_LITERAL,
             "signed normalization/Pstar literal changed")
    _require(row["conditional_common_point_pullback_unbound"] is True,
             "common physical source-domain status changed")


def _bind_history_rows(signed_report, joint_report, manifest_by_namespace, joint_graph):
    c3_roles, c5_roles = _ordinary_roles(joint_report)
    _require(c3_roles["namespace"] == "repair_C3_band", "C3 namespace changed")
    _require(c5_roles["namespace"] == "true_compact_C5", "C5 namespace changed")
    c3_ports = c3_roles["terminal_history_ports_Z0_to_Z3"]
    c5_ports = c5_roles["terminal_scaled_history_ports_Z0_to_Z3"]
    _require({key: tuple(value) for key, value in c3_ports.items()} ==
             EXPECTED_C3_TERMINAL_PORTS, "joint C3 terminal ports changed")
    _require({key: tuple(value) for key, value in c5_ports.items()} ==
             EXPECTED_C5_SCALED_PORTS, "joint C5 scaled ports changed")

    signed = _signed_rows_by_key(signed_report)
    bound_rows = []
    for field in HISTORY_FIELDS:
        channel = CHANNEL_FOR_FIELD[field]
        for order in ORDINARY_Z_ORDERS:
            c5_node = c5_ports[field][order]
            c3_node = c3_ports[field][order]
            row = signed[(field, order)]
            _check_gap_row(row, field, order, c5_node, c3_node)
            formula_binding = joint_graph["channels"][channel][str(order)]
            bound_rows.append({
                "field": field,
                "ordinary_Z_order": order,
                "channel": channel,
                "signed_history_row": row["gap_index"],
                "signed_history_row_pointer": (
                    f"$.full_result.signed_history_gap_function_nodes[{row['gap_index']}]"
                ),
                "signed_gap": {
                    "left_namespace": "true_compact_C5",
                    "left_source_node": c5_node,
                    "left_raw_correction_node": row["left"]["raw_correction_node"],
                    "left_raw_correction_role": row["left"]["raw_correction_role"],
                    "right_namespace": "repair_C3_band",
                    "right_source_node": c3_node,
                    "operation": "signed_difference",
                },
                "source_manifest_witnesses": {
                    "C5": _manifest_witness(
                        manifest_by_namespace, "C5", "history", field, order, c5_node
                    ),
                    "C3": _manifest_witness(
                        manifest_by_namespace, "C3", "history", field, order, c3_node
                    ),
                },
                "joint_formula_binding": {
                    "raw_gap_role": RAW_GAP_ROLE[channel],
                    "formula_literal": JOINT_NUMERATOR_FORMULAS[channel],
                    "actual_numerator_literal": formula_binding["numerator_literal"],
                    "actual_denominator_literal": formula_binding["denominator_literal"],
                    "actual_atom_witnesses": formula_binding["atom_witnesses"],
                    "json_pointers": formula_binding["json_pointers"],
                },
                "source_unit_contract": {
                    "normalization_literal": row["normalization"],
                    "normalization_pointer": (
                        f"$.full_result.signed_history_gap_function_nodes[{row['gap_index']}]"
                        ".normalization"
                    ),
                    "C5_N_scaling": "raw correction times N exactly once",
                    "C3_N_scaling": "terminal history already N-scaled",
                    "physical_Pstar_lift_applied": False,
                    "Pstar_status": "present only in generic normalized source literal; no physical lift",
                    "common_physical_domain_claimed": False,
                    "same_units_claimed": False,
                    "gap_zero_claimed": False,
                },
            })
    return bound_rows


def _bind_amplitude(joint_report, manifest_by_namespace):
    c3_roles, c5_roles = _ordinary_roles(joint_report)
    c3 = tuple(c3_roles["amplitude_ports_Z0_to_Z3"])
    c5 = tuple(c5_roles["amplitude_ports_Z0_to_Z3"])
    _require(c3 == EXPECTED_C3_AMPLITUDE_PORTS, "C3 amplitude ports changed")
    _require(c5 == EXPECTED_C5_AMPLITUDE_PORTS, "C5 amplitude ports changed")
    return [
        {
            "ordinary_Z_order": order,
            "C5_source_manifest_witness": _manifest_witness(
                manifest_by_namespace, "C5", "amplitude", "A", order, c5[order]
            ),
            "C3_source_manifest_witness": _manifest_witness(
                manifest_by_namespace, "C3", "amplitude", "A", order, c3[order]
            ),
            "gap_role": "a = A_N-A_C",
            "source_units": "amplitude source units only; no physical lift",
            "gap_zero_claimed": False,
            "identity_scope": "canonical leading-background terminal E source-query only",
        }
        for order in ORDINARY_Z_ORDERS
    ]


def consume_saved_reports(repo_root):
    """Build a literal binding packet; intentionally never called in this draft.

    The returned packet would be a source identity between two already saved
    source packets, not a numerical target-gap evaluation.  In particular, it
    does not claim C3/C5 equality, physical unit equality, volume binding, or
    cap transfer.
    """
    repo_root = Path(repo_root)
    archive_raw = _check_authenticated_archives(repo_root)
    signed_path = repo_root / SIGNED_REPORT_RELATIVE
    joint_path = repo_root / JOINT_REPORT_RELATIVE
    signed_report = _read_authenticated_json(signed_path, SIGNED_REPORT_SHA256)
    joint_report = _read_authenticated_json(joint_path, JOINT_REPORT_SHA256)
    _check_saved_contract_flags(signed_report, joint_report)
    _check_joint_gap_definitions(joint_report)
    common_N = _check_common_N_descriptor(signed_report)
    c3_roles, c5_roles = _ordinary_roles(joint_report)
    manifest_by_namespace = _check_joint_manifest(joint_report, c3_roles, c5_roles)
    retained_legacy = _check_retained_legacy_provenance(signed_report, joint_report)
    joint_graph = _bind_joint_formula_graph(
        joint_report, manifest_by_namespace, c3_roles, c5_roles
    )
    history = _bind_history_rows(
        signed_report, joint_report, manifest_by_namespace, joint_graph
    )
    amplitude = _bind_amplitude(joint_report, manifest_by_namespace)
    history_by_key = {(row["field"], row["ordinary_Z_order"]): row for row in history}
    d_joint = []
    for order in ORDINARY_Z_ORDERS:
        k_row = history_by_key[("k", order)]
        m_row = history_by_key[("m", order)]
        d_joint.append({
            "ordinary_Z_order": order,
            "history_row_refs": {
                "k": k_row["signed_history_row"],
                "m": m_row["signed_history_row"],
            },
            "symbolic_gap_refs": ("delta_k", "delta_m"),
            "joint_numerator_formula": JOINT_NUMERATOR_FORMULAS["D"],
            "k_formula_binding": k_row["joint_formula_binding"],
            "m_formula_binding": m_row["joint_formula_binding"],
            "actual_D_graph": joint_graph["channels"]["D"][str(order)],
            "extra_N_divisions": 0,
        })
    return {
        "schema": SCHEMA,
        "source_only": True,
        "source_binding_defined": True,
        "authenticated_report_sha256": {
            "signed_history": SIGNED_REPORT_SHA256,
            "joint_numerator": JOINT_REPORT_SHA256,
        },
        "authenticated_archive_raw_sha256": {
            "C3": C3_ARCHIVE_RAW_SHA256,
            "C5": C5_ARCHIVE_RAW_SHA256,
        },
        "archive_raw_authentication": archive_raw,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "history_row_count": len(history),
        "history_bindings": history,
        "amplitude_bindings": amplitude,
        "common_N_descriptor": common_N,
        "manifest_namespace_aliases": manifest_by_namespace,
        "joint_formula_graph": joint_graph,
        "channel_field_contract": {
            channel: list(fields) for channel, fields in JOINT_HISTORY_FIELDS.items()
        },
        "joint_numerator_formula_contract": JOINT_NUMERATOR_FORMULAS,
        "D_joint_k_m_bindings": d_joint,
        "joint_numerator_retained": True,
        "retained_legacy_provenance": retained_legacy,
        "source_equality_admitted": False,
        "signed_target_gaps_evaluated": False,
        "physical_Pstar_lift_claimed": False,
        "compact_volume_source_binding_claimed": False,
        "error_cap_transfer_claimed": False,
        "source_inspection_audit": {
            "source_modules_imported": 0,
            "source_functions_executed": 0,
            "numeric_values_evaluated": 0,
            "symbolic_math_operations": 0,
            "range_or_midpoint_values_used": False,
        },
        "report_pointers": REPORT_POINTERS,
        "remaining_scope": [
            "Bind the resulting gap nodes to a common physical source domain.",
            "Resolve physical scale conventions without adding Pstar by name.",
            "Keep all retained volume/error rows and the independent reference offset unchanged.",
            "The saved reports expose the reference offset only as an unresolved scope marker; no numeric offset is selected.",
            "No common physical domain, gap equality, or cap transfer is asserted by this packet.",
        ],
    }


# No __main__ entry point is provided: this draft must not be executed here.
