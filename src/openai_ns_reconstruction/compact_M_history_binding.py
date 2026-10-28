"""Fresh, source-only M signed-history port binding.

This leaf reads two already-saved JSON publications and selects only their
explicit ordinary-Z role vectors.  It does not import a reconstruction
module, evaluate a source, decode an archive, walk a graph, or access the
failed joint atom family.  Field selection is deliberately completed before
an ordinary-Z order is used; this avoids the old ``nodes[order]`` mistake.

The returned packet keeps the quotient formula symbolic.  It records the
source-owned C3/C5 M endpoint, incoming ``m`` history, and amplitude ports for
Z0--Z3, while retaining ``a = A_N-A_C``, the P1/m15 distinction, and the old
212/213/183 error ledger.  No function value, seam, common physical point,
unit identification, or error reduction is claimed.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


SCHEMA = "compact_M_history_binding_v1"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
FORMULA = "M_gap = (A_C * delta_m - a * T_C_m) / (A_C * A_N)"

JOINT_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_"
    "joint_target_numerator_transport.json"
)
JOINT_REPORT_SHA256 = (
    "cff727e06c1af18eebeba0bbf027a488128c2021f9feb839bafdfc368b18f133"
)

FAILURE_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_compact_"
    "joint_history_gap_binding_failure.json"
)
FAILURE_REPORT_SHA256 = (
    "5f2a57011c3a07f3ac4ef6a21d794678aa79225af98d2be815864a099b498468"
)

RAW_ARCHIVE_SHA256 = {
    "C3": "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "C5": "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_saved_json(root: Path, relative: Path, expected_sha256: str) -> dict:
    """Authenticate and read one saved publication, without source evaluation."""
    path = root / relative
    _require(path.is_file(), f"missing saved report: {relative}")
    _require(
        _sha256_file(path) == expected_sha256,
        f"saved report SHA changed: {relative}",
    )
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"saved report is not an object: {relative}")
    return value


def _ordered_field_ports(container: dict, field: str, label: str) -> tuple[int, ...]:
    """Select a named field first, then read its four ordinary-Z positions."""
    _require(isinstance(container, dict), f"{label} is not a field map")
    selected = container.get(field)
    _require(
        isinstance(selected, (list, tuple)) and len(selected) == 4,
        f"{label}.{field} is not a four-position source vector",
    )
    values = []
    for order in ORDINARY_Z_ORDERS:
        value = selected[order]
        _require(isinstance(value, int) and not isinstance(value, bool),
                 f"{label}.{field}[{order}] is not an integer source handle")
        values.append(value)
    return tuple(values)


def _read_M_roles(joint_report: dict) -> dict:
    """Read only the explicit M, incoming-m, and amplitude role vectors."""
    contract = joint_report.get("source_recipe_contract")
    _require(isinstance(contract, dict), "joint source-role contract is absent")
    roles = contract.get("ordinary_Z_roles")
    _require(isinstance(roles, dict), "ordinary-Z source roles are absent")
    _require(tuple(roles.get("orders", ())) == ORDINARY_Z_ORDERS,
             "ordinary-Z role order changed")
    _require(roles.get("orders_are_source_roles") is True,
             "ordinary-Z orders are not marked as source roles")
    _require(roles.get("range_values_are_not_function_values") is True,
             "ordinary-Z role ranges were marked as function values")

    result = {}
    for namespace, history_key, history_role in (
        ("C3", "terminal_history_ports_Z0_to_Z3", "terminal_history"),
        ("C5", "terminal_scaled_history_ports_Z0_to_Z3", "terminal_scaled_history"),
    ):
        side = roles.get(namespace)
        _require(isinstance(side, dict), f"{namespace} source roles are absent")
        expected_role_namespace = {
            "C3": "repair_C3_band",
            "C5": "true_compact_C5",
        }[namespace]
        _require(side.get("namespace") == expected_role_namespace,
                 f"{namespace} source namespace changed")
        history = side.get(history_key)
        result[namespace] = {
            "role_namespace": side["namespace"],
            "M_endpoint": _ordered_field_ports(
                side.get("target_ports_Z0_to_Z3"), "M",
                f"{namespace} target_ports_Z0_to_Z3",
            ),
            "m_history": _ordered_field_ports(
                history, "m", f"{namespace} {history_key}"
            ),
            "amplitude": tuple(
                _ordered_field_ports(
                    {"A": side.get("amplitude_ports_Z0_to_Z3")},
                    "A", f"{namespace} amplitude_ports_Z0_to_Z3",
                )
            ),
            "history_key": history_key,
            "history_role": history_role,
        }
    return result


def _read_prior_failure(failure_report: dict) -> dict:
    """Retain the failed-family boundary without touching its embedded data."""
    _require(
        failure_report.get("schema") ==
        "compact_joint_history_gap_binding_failure_publication_v1",
        "unexpected prior failure publication schema",
    )
    _require(failure_report.get("binding_accepted") is False,
             "prior failed joint family was marked accepted")
    _require(failure_report.get("exception") == "KeyError: 0",
             "prior failure exception changed")
    _require(failure_report.get("old_errors_removed") is False,
             "prior failure changed old-error ledger")
    _require(failure_report.get("same_failed_family_must_not_be_reentered") is True,
             "prior failure did not retain the no-reentry guard")
    return {
        "path": str(FAILURE_REPORT).replace("\\", "/"),
        "report_sha256": FAILURE_REPORT_SHA256,
        "status": failure_report.get("status"),
        "exception": failure_report.get("exception"),
        "failure_stage": failure_report.get("failure_stage"),
        "binding_accepted": False,
        "reentered": False,
        "old_errors_removed": False,
        "same_failed_family_must_not_be_reentered": True,
    }


def bind_M_history_binding(repo_root: str | Path) -> dict:
    """Return four M source-port rows using saved roles only.

    This is intentionally a source-reference packet.  The only order-indexed
    reads occur after selecting the named ``M``, ``m``, or ``A`` field.
    """
    root = Path(repo_root)
    joint_report = _read_saved_json(root, JOINT_REPORT, JOINT_REPORT_SHA256)
    failure_report = _read_saved_json(
        root, FAILURE_REPORT, FAILURE_REPORT_SHA256
    )
    roles = _read_M_roles(joint_report)
    prior_failure = _read_prior_failure(failure_report)

    rows = []
    for order in ORDINARY_Z_ORDERS:
        c3 = roles["C3"]
        c5 = roles["C5"]
        rows.append({
            "channel": "M",
            "ordinary_Z_order": order,
            "C3_source": {
                "namespace": "C3",
                "role_namespace": c3["role_namespace"],
                "history_field": "m",
                "history_source_node": c3["m_history"][order],
                "history_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C3."
                    "terminal_history_ports_Z0_to_Z3.m"
                ),
                "amplitude_source_node": c3["amplitude"][order],
                "amplitude_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C3."
                    "amplitude_ports_Z0_to_Z3"
                ),
                "M_endpoint_source_node": c3["M_endpoint"][order],
                "M_endpoint_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C3."
                    "target_ports_Z0_to_Z3.M"
                ),
                "archive_raw_sha256": RAW_ARCHIVE_SHA256["C3"],
            },
            "C5_source": {
                "namespace": "C5",
                "role_namespace": c5["role_namespace"],
                "history_field": "m",
                "history_source_node": c5["m_history"][order],
                "history_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C5."
                    "terminal_scaled_history_ports_Z0_to_Z3.m"
                ),
                "amplitude_source_node": c5["amplitude"][order],
                "amplitude_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C5."
                    "amplitude_ports_Z0_to_Z3"
                ),
                "M_endpoint_source_node": c5["M_endpoint"][order],
                "M_endpoint_source_role": (
                    "source_recipe_contract.ordinary_Z_roles.C5."
                    "target_ports_Z0_to_Z3.M"
                ),
                "archive_raw_sha256": RAW_ARCHIVE_SHA256["C5"],
            },
            "terminal_histories": {
                "C3": "T_C_m from the C3 terminal m history",
                "C5": "T_N_m from the C5 N-scaled terminal m history",
                "values_selected": False,
            },
            "saved_role_references_selected": True,
            "same_source_semantics_proved": False,
            "cross_archive_source_identity_proved": False,
            "formula": {
                "conditional_on_same_normalization_and_source_interpretation": True,
                "M_C": "T_C_m / A_C",
                "M_N": "T_N_m / A_N",
                "delta_m": "T_N_m-T_C_m",
                "a": "A_N-A_C",
                "M_gap": FORMULA,
                "a_zero_substitution_used": False,
            },
        })

    return {
        "schema": SCHEMA,
        "source_only": True,
        "binding_kind": "M_only_signed_history_source_port_binding",
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "inputs": {
            "joint_role_publication": {
                "path": str(JOINT_REPORT).replace("\\", "/"),
                "report_sha256": JOINT_REPORT_SHA256,
                "role_pointer": (
                    "$.source_recipe_contract.ordinary_Z_roles"
                ),
                "M_endpoint_rows_read": 4,
                "terminal_m_history_rows_read": 4,
                "amplitude_rows_read": 4,
                "M_atom_subtree_accessed": False,
                "graph_traversed": False,
            },
            "prior_failure_publication": prior_failure,
            "gzip_decoded": False,
        },
        "source_port_binding": {
            "M_C3_C5_endpoint_history_amplitude_rows_selected": 4,
            "actual_cross_packet_binding_proved": False,
            "M_field_selected_before_order_index": True,
            "incoming_histories_assumed_zero": False,
            "source_values_computed": False,
            "signed_gap_evaluated": False,
            "endpoint_equals_T_over_A_proved": False,
            "M_Z0_Z3_function_equality_proved": False,
        },
        "ordinary_Z_bindings": rows,
        "N_policy": {
            "N_node_inherited_unverified": 18,
            "invN_node_inherited_unverified": 19,
            "N_materialized": False,
            "extra_N_divisions": 0,
        },
        "symbol_separation": {
            "a_definition": "a=A_N-A_C",
            "P1_is_not_m15": True,
            "M_distinct_from_m15": True,
            "M_equals_m15_proved": False,
        },
        "retained_legacy_scope": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "J_terms": 16,
            "K_terms": 17,
            "outer_B_terms": 2,
            "reserved_contribution_registered_once": True,
            "six_signed_channels_retained": True,
            "all_old_error_data_retained": True,
            "inherited_errors_deleted": False,
            "cap_transfer": False,
        },
        "nonclaims": {
            "signed_gap_value": False,
            "M_gap_zero": False,
            "M_Z0_Z3_function_value": False,
            "common_physical_point": False,
            "common_physical_domain": False,
            "same_units": False,
            "physical_Pstar_lift": False,
            "finite_N_error_bound": False,
            "M_equals_m15": False,
            "reference_offset_cancellation": False,
            "old_error_reduction": False,
            "cap_transfer": False,
        },
        "runtime_model_metadata": None,
    }


def load_M_history_binding(repo_root: str | Path) -> dict:
    """Compatibility spelling for the one permitted source consumer call."""
    return bind_M_history_binding(repo_root)
