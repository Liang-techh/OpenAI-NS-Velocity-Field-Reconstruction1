"""Source-only conditional P-gap minus current-radius G-gap expression."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "pulse_gap_physical_current_difference_v1"
ORDERS = (0, 1)
ALL_ORDERS = (0, 1, 2, 3)
PHYSICAL_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json"
)
G_RESULT = Path("tmp/pulse_gap_current_radius_source_v1/result.json")
EXPECTED_PHYSICAL_REPORT_SHA256 = "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b"
EXPECTED_G_RESULT_SHA256 = "c6d170f6935e0dad3937e155fbbb7e67c6590ee19369291594d10479883fef05"
P_OWNER, P_NAMESPACE = "P", "physical_current_native_report"
G_OWNER, G_NAMESPACE = "G", "current_radius_saved_result"
P_BINDING = "$.component_result.physical_pulse_bindings"
P_GAP_M15 = f"{P_BINDING}.pieces.gap.m15"
P_GAP_ENDPOINTS = f"{P_BINDING}.endpoints.gap"
G_SOURCE = "$.result.G_current_radius_source"
G_GAP_ENDPOINTS = f"{G_SOURCE}.endpoints.gap"
P_PRESSURE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py"
    ":pulse_pressure_functions"
)
P_PRESSURE_MAP = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions.py"
    ":pressure_functions"
)
P_ENDPOINT_FIELDS = (
    "R_q", "R_squared_q", "R_q_squared", "R_squared_q_squared",
    "R_U_squared", "R_squared_U_squared", "q_plus_q_y", "R_q_y",
)


class SourceExpressionInputError(ValueError):
    """An authenticated saved source boundary is missing or changed."""


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise SourceExpressionInputError(message)


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _rel(path: Path) -> str:
    return path.as_posix()


def _read_json(root: Path, relative: Path, expected: str, label: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = root / relative
    _require(path.is_file(), f"missing {label}: {_rel(relative)}")
    raw = path.read_bytes()
    actual = _sha256(raw)
    _require(actual == expected, f"{label} raw SHA-256 changed: {actual} != {expected}")
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceExpressionInputError(f"invalid JSON in {label}") from exc
    _require(isinstance(value, Mapping), f"{label} top level is not an object")
    return value, {
        "path": _rel(relative),
        "expected_raw_sha256": expected,
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": False,
    }


def _map(value: Any, pointer: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"expected object at {pointer}")
    return value


def _rows(
    value: Any,
    pointer: str,
    owner: str,
    namespace: str,
    *,
    saved_ref_dicts: bool = False,
    saved_source_pointer: str | None = None,
) -> dict[str, Any]:
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes)), f"expected rows at {pointer}")
    _require(len(value) == 4, f"expected four ordinary-Z rows at {pointer}")
    result = []
    for order, row in enumerate(value):
        row_pointer = f"{pointer}[{order}]"
        if isinstance(row, int) and not isinstance(row, bool):
            _require(not saved_ref_dicts, f"G row is not a saved direct-ref dictionary at {row_pointer}")
            source_node_id = row
        else:
            ref = _map(row, row_pointer)
            _require(saved_ref_dicts, f"P row unexpectedly wrapped at {row_pointer}")
            _require(ref.get("status") == "direct_source_ref", f"row status changed at {row_pointer}")
            expected_saved_role = (
                f"{saved_source_pointer}[{order}]"
                if saved_source_pointer is not None
                else row_pointer
            )
            _require(ref.get("source_role") == expected_saved_role, f"row source role changed at {row_pointer}")
            _require(ref.get("node_id_is_not_a_numeric_value") is True, f"opaque marker missing at {row_pointer}")
            _require(ref.get("source_value_selected") is False, f"source value selected at {row_pointer}")
            source_node_id = ref.get("source_node_id")
            _require(isinstance(source_node_id, int) and not isinstance(source_node_id, bool), f"node id missing at {row_pointer}")
        result.append({
            "status": "direct_source_ref",
            "source_role": (
                f"{saved_source_pointer}[{order}]"
                if saved_source_pointer is not None
                else row_pointer
            ),
            "source_node_id": source_node_id,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
            "owner": owner,
            "namespace": namespace,
            "ordinary_Z_order": order,
        })
    return {
        "status": "direct_source_rows",
        "owner": owner,
        "namespace": namespace,
        "source_role": pointer,
        "ordinary_Z_orders": list(ALL_ORDERS),
        "ordinary_Z0_to_Z1": result[:2],
        "all_saved_rows": result,
        "source_values_selected": False,
    }


def _endpoint_roles(binding: Mapping[str, Any]) -> dict[str, Any]:
    endpoints = _map(binding.get("endpoints"), f"{P_BINDING}.endpoints")
    raw = endpoints.get("gap")
    _require(isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)) and len(raw) == 2, f"expected two P gap endpoints at {P_GAP_ENDPOINTS}")
    packet: dict[str, Any] = {"status": "direct_physical_finite_endpoint_source_roles", "owner": P_OWNER, "namespace": P_NAMESPACE, "chart": "gap", "orientation": "lower_to_upper", "source_values_selected": False}
    for index, label in enumerate(("lower", "upper")):
        endpoint = _map(raw[index], f"{P_GAP_ENDPOINTS}[{index}]")
        fields = {
            field: _rows(endpoint.get(field), f"{P_GAP_ENDPOINTS}[{index}].{field}", P_OWNER, P_NAMESPACE)
            for field in P_ENDPOINT_FIELDS
        }
        packet[label] = {"source_role": f"{P_GAP_ENDPOINTS}[{index}]", "fields": fields, "source_values_selected": False}
    return packet


def _p_rows(physical: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    _require(physical.get("schema") == "original_native_physical_pulse_source_stage_v1", "physical top-level schema changed")
    component = _map(physical.get("component_result"), "$.component_result")
    _require(component.get("schema") == "original_native_physical_pulse_source_v1", "physical component schema changed")
    _require(component.get("passed") is True, "physical component_result.passed is not true")
    _require(component.get("source_only") is True, "physical component is not source-only")
    binding = _map(component.get("physical_pulse_bindings"), P_BINDING)
    pieces = _map(binding.get("pieces"), f"{P_BINDING}.pieces")
    gap = _map(pieces.get("gap"), f"{P_BINDING}.pieces.gap")
    p_rows = _rows(gap.get("m15"), P_GAP_M15, P_OWNER, P_NAMESPACE)
    endpoint_packet = _endpoint_roles(binding)
    corrections = _rows(binding.get("pulse_geometry_correction_m15"), f"{P_BINDING}.pulse_geometry_correction_m15", P_OWNER, P_NAMESPACE)
    return p_rows, endpoint_packet, corrections


def _g_rows(document: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    _require(document.get("component") == "pulse_gap_current_radius_source_v1", "G component schema changed")
    _require(document.get("status") == "returned", "G result status changed")
    result = _map(document.get("result"), "$.result")
    _require(result.get("schema") == "pulse_gap_current_radius_source_v1", "G source schema changed")
    _require(result.get("source_bound_G") is True, "G source_bound_G is not true")
    _require(result.get("source_only") is True, "G source is not source-only")
    _require(result.get("admissible_for_numeric_use") is False, "G unexpectedly admits numeric use")
    source = _map(result.get("G_current_radius_source"), G_SOURCE)
    _require(source.get("owner") == "G", "G source owner changed")
    endpoints = _map(source.get("endpoints"), f"{G_SOURCE}.endpoints")
    gap = _map(endpoints.get("gap"), G_GAP_ENDPOINTS)
    packets: dict[str, Any] = {}
    for side in ("left", "right"):
        side_pointer = f"{G_GAP_ENDPOINTS}.{side}"
        side_data = _map(gap.get(side), side_pointer)
        increment = _map(side_data.get("signed_pressure_increment"), f"{side_pointer}.signed_pressure_increment")
        pointer = f"{side_pointer}.signed_pressure_increment.m15"
        packets[side] = _rows(
            increment.get("m15"),
            pointer,
            G_OWNER,
            G_NAMESPACE,
            saved_ref_dicts=True,
            saved_source_pointer=(
                "$.pulse_radial_pressure_C3_bindings.endpoints.gap."
                f"{side}.signed_pressure_increment.m15"
            ),
        )
    return packets["left"], packets["right"]


def _retained_scope(physical: Mapping[str, Any]) -> dict[str, Any]:
    retained = _map(physical.get("retained_construction_data"), "$.retained_construction_data")
    expected = {"cumulative_terms": 212, "raw_signed_terms": 213, "reduced_terms": 183}
    counts: dict[str, Any] = {}
    for key, wanted in expected.items():
        value = retained.get(key)
        _require(isinstance(value, int) and not isinstance(value, bool) and value == wanted, f"retained count changed at {key}")
        counts[key] = {"value": value, "kind": "saved_metadata_count_not_source_node_id", "source_role": f"$.retained_construction_data.{key}"}
    _require(retained.get("six_signed_channels_retained") is True, "six signed channels not retained")
    _require(retained.get("reserved_contribution_registered_once") is True, "reserved contribution retention changed")
    _require(retained.get("all_old_error_data_retained") is True, "old error data not retained")
    _require(retained.get("inherited_errors_deleted") is False, "inherited errors marked deleted")
    return {"status": "retained_not_added_to_expression", "range_label": "212 cumulative / 213 raw / 183 reduced", "counts": counts, "six_signed_channels_retained": True, "reserved_contribution_registered_once": True, "all_old_error_data_retained": True, "inherited_errors_deleted": False, "source_values_selected": False}


def _expr(operator: str, *args: Any) -> dict[str, Any]:
    return {"operator": operator, "arguments": list(args), "source_values_selected": False}


def _failed(provenance: Mapping[str, Any], reason: str) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "failed_closed_source_boundary",
        "source_only": True,
        "formal_source_expression_defined": False,
        "admissible_for_numeric_use": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "source_values_selected": False,
        "diagnostic": reason,
        "unresolved_gates": ["saved_input_evidence_complete"],
        "provenance": dict(provenance),
    }


def build_gap_physical_current_difference(repo_root: str | Path) -> dict[str, Any]:
    """Build P-gap minus the saved G right-minus-left increment, Z0/Z1 only."""

    root = Path(repo_root)
    provenance: dict[str, Any] = {
        "physical_report": {"path": _rel(PHYSICAL_REPORT), "expected_raw_sha256": EXPECTED_PHYSICAL_REPORT_SHA256},
        "G_result": {"path": _rel(G_RESULT), "expected_raw_sha256": EXPECTED_G_RESULT_SHA256},
    }
    try:
        physical, p_meta = _read_json(root, PHYSICAL_REPORT, EXPECTED_PHYSICAL_REPORT_SHA256, "physical report")
        g_document, g_meta = _read_json(root, G_RESULT, EXPECTED_G_RESULT_SHA256, "G saved result")
        provenance["physical_report"] = p_meta
        provenance["G_result"] = g_meta
        p_rows, p_endpoints, corrections = _p_rows(physical)
        g_left, g_right = _g_rows(g_document)
        retained = _retained_scope(physical)
    except (SourceExpressionInputError, OSError) as exc:
        return _failed(provenance, str(exc))

    expressions = []
    for order in ORDERS:
        p_ref = p_rows["all_saved_rows"][order]
        right_ref = g_right["all_saved_rows"][order]
        left_ref = g_left["all_saved_rows"][order]
        g_gap = _expr("subtract", right_ref, left_ref)
        expressions.append({
            "ordinary_Z_order": order,
            "operation_order": ["G_right", "G_left", "G_right_minus_G_left", "P_gap", "P_gap_minus_G_gap"],
            "P_gap": p_ref,
            "G_right_cumulative_endpoint_increment": right_ref,
            "G_left_cumulative_endpoint_increment": left_ref,
            "G_gap": g_gap,
            "P_gap_minus_G_gap": _expr("subtract", p_ref, g_gap),
            "source_values_selected": False,
            "numeric_difference_evaluated": False,
        })

    return {
        "schema": SCHEMA,
        "status": "conditional_source_expression_not_admissible_for_numeric_use",
        "source_only": True,
        "formal_source_expression_defined": True,
        "admissible_for_numeric_use": False,
        "source_identity_proved": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "source_values_selected": False,
        "inputs": {
            "physical_report": provenance["physical_report"],
            "G_result": provenance["G_result"],
            "physical_component_result_passed": True,
            "physical_component_schema": "original_native_physical_pulse_source_v1",
            "source_bound_G": True,
        },
        "provenance": {
            "P": {"owner": P_OWNER, "namespace": P_NAMESPACE, "binding": P_BINDING, "gap_piece": P_GAP_M15, "finite_endpoints": P_GAP_ENDPOINTS},
            "G": {"owner": G_OWNER, "namespace": G_NAMESPACE, "saved_source": G_SOURCE, "finite_endpoints": G_GAP_ENDPOINTS, "left_m15": f"{G_GAP_ENDPOINTS}.left.signed_pressure_increment.m15", "right_m15": f"{G_GAP_ENDPOINTS}.right.signed_pressure_increment.m15", "source_bound_G": True},
            "inputs": provenance,
        },
        "P_gap_local_piece": {
            "rows": p_rows,
            "finite_endpoint_source_roles": p_endpoints,
            "source_definition": {
                "status": "source_pointer_only_not_authenticated_by_this_leaf",
                "source_function": P_PRESSURE_SOURCE,
                "delegated_signed_pressure_map": P_PRESSURE_MAP,
                "chart": "gap",
                "local_piece": "pieces.gap.m15",
                "physical_measure_weights": {"I0": "(R/mu) dxi", "I1": "(R^2/mu) dxi"},
                "finite_endpoint_bracket": "upper.R_q_y - lower.R_q_y",
                "signed_m15_definition": "m15 = S1/2 - (upper.R_q_y - lower.R_q_y)",
                "original_m15_negative_weighted_integral_sign_applied_once": True,
                "source_authentication_performed": False,
                "source_values_selected": False,
            },
        },
        "G_gap_cumulative_endpoint_increments": {"left": g_left, "right": g_right, "orientation": "right_minus_left", "signed_pressure_increment": True, "source_values_selected": False},
        "ordinary_Z_orders": list(ORDERS),
        "formal_difference": {"formula": "P_gap - (G_right - G_left)", "expression_kind": "source_function_definition", "not_an_upper_bound": True, "status": "conditional_source_expression_not_evaluated", "expressions": expressions, "source_values_selected": False, "numeric_difference_evaluated": False},
        "conditional_proof_gates": {
            "status": "unresolved_conditional_gate_set",
            "all_proved": False,
            "requirements": [
                "common_physical_domain",
                "common_mu",
                "common_gap_coordinate",
                "common_Pstar",
                "physical_units_and_measure_normalization",
                "ordinary_Z_derivative_and_Z1_transport",
                "signed_pressure_map_equivalence",
                "finite_endpoint_function_equivalence",
                "right_minus_left_orientation_matches_local_gap_piece",
                "common_gap_radius_and_L_nonzero_domain",
            ],
            "source_identity_proved": False,
        },
        "historical_G_minus_F": {"status": "retained_unresolved_source_difference", "source_role": "current-radius G source minus historical frame-only F source", "preserved": True, "evaluated": False, "zero_proved": False, "source_identity_proved": False},
        "P_minus_F_decomposition": {"exact_additive_form": "(P-G)+(G-F)", "G_minus_F_retained": True, "either_term_zero_proved": False, "numeric_values_computed": False},
        "retained_prior_physical_pulse_corrections": {"status": "direct_source_roles_retained_not_subtracted", "physical_current_pulse_geometry_correction_m15": corrections, "correction_count": 4, "correction_removed": False},
        "retained_212_scope": retained,
        "independent_reference_offset": {"status": "retained_unresolved_symbolic_marker", "source_role": "independent_reference_offset", "selected": False, "numeric_value_computed": False, "included_in_P_minus_G_expression": False},
        "unresolved_gates": [
            "common_physical_domain",
            "common_mu",
            "common_gap_coordinate",
            "common_Pstar",
            "physical_units_and_measure_normalization",
            "ordinary_Z_derivative_and_Z1_transport",
            "signed_pressure_map_equivalence",
            "finite_endpoint_function_equivalence",
            "right_minus_left_orientation_matches_local_gap_piece",
            "common_gap_radius_and_L_nonzero_domain",
        ],
        "scope_audit": {"raw_hash_authentication_before_parse": True, "gzip_archives_decoded": 0, "source_graph_traversed": False, "graph_id_arithmetic": False, "source_constructors_called": 0, "source_checkers_called": 0, "old_family_replayed": False, "incident_pair_DAG_reentered": False, "root_mapping_performed": False, "quadratures": 0, "numeric_values_computed": False, "numeric_difference_evaluated": False, "source_values_selected": False, "old_errors_removed": False, "extra_terms_added": False},
    }


__all__ = ["build_gap_physical_current_difference"]
