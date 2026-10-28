"""Source-only current physical/frame-only gap expression consumer.

Only two frozen JSON artifacts are consumed.  This leaf does not import a
source constructor or checker, inspect an ``exact_graph_nodes`` array, follow
graph edges, evaluate node ids, integrate a field, or compare the two source
owners.  The gzip archive is authenticated on raw bytes before decoding.  The
returned packet contains opaque source references and conditional symbolic
expressions; it contains no source values or numerical differences.

The accepted OriginalNativePulseM15 gap helper is deliberately not called.
This module supplies the separate current-physical/frame-only source boundary
needed by a later evaluator.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "pulse_gap_actual_source_program_v1"
ORDERS = (0, 1, 2)

PHYSICAL_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json"
)
FRAME_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
RAW_SHA256 = {
    "physical": "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    "frame": "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
}


class SourceProgramInputError(ValueError):
    """A frozen input boundary is absent or has changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceProgramInputError(message)


def _raw_sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _relative(path: Path) -> str:
    return path.as_posix()


def _read_json_after_raw_auth(
    root: Path, relative: Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = root / relative
    _require(path.is_file(), f"missing saved report: {_relative(relative)}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(actual == expected_sha256, f"raw SHA changed: {_relative(relative)}")
    value = json.loads(raw.decode("utf-8"))
    _require(isinstance(value, Mapping), f"saved report is not an object: {_relative(relative)}")
    return value, {
        "path": _relative(relative),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": False,
    }


def _read_gzip_json_after_raw_auth(
    root: Path, relative: Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    """Authenticate raw gzip bytes before the first decompression call."""

    path = root / relative
    _require(path.is_file(), f"missing raw gzip archive: {_relative(relative)}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(actual == expected_sha256, f"raw gzip SHA changed: {_relative(relative)}")
    decoded = gzip.decompress(raw)
    value = json.loads(decoded.decode("utf-8"))
    _require(isinstance(value, Mapping), f"decoded archive is not an object: {_relative(relative)}")
    return value, {
        "path": _relative(relative),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": True,
    }


def _direct(value: Any, path: Sequence[str], label: str) -> Any:
    """Read a literal JSON path; this never interprets a graph node."""

    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            raise SourceProgramInputError(f"missing direct source field {label}")
        current = current[key]
    return current


def _optional(value: Any, path: Sequence[str]) -> Any | None:
    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            return None
        current = current[key]
    return current


def _node_ref(value: Any, pointer: str, label: str, gates: list[str]) -> dict[str, Any]:
    """Keep a saved node id opaque; never treat it as a numeric value."""

    if not isinstance(value, int) or isinstance(value, bool):
        gates.append(f"{label}:direct_source_node_id_missing_or_noninteger")
        return {
            "status": "unresolved",
            "source_role": pointer,
            "node_id_is_not_a_numeric_value": True,
        }
    return {
        "status": "direct_source_ref",
        "source_role": pointer,
        "source_node_id": value,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _row_refs(
    value: Any, pointer: str, label: str, gates: list[str]
) -> dict[str, Any]:
    """Expose saved ordinary-Z rows without differentiating or traversing them."""

    if not isinstance(value, list) or len(value) < len(ORDERS):
        gates.append(f"{label}:fewer_than_three_saved_ordinary_Z_rows")
        return {
            "status": "unresolved",
            "source_role": pointer,
            "ordinary_Z_orders": list(ORDERS),
            "source_value_selected": False,
        }
    all_rows = [
        _node_ref(item, f"{pointer}[{index}]", f"{label}[{index}]", gates)
        for index, item in enumerate(value)
    ]
    return {
        "status": "direct_source_rows",
        "source_role": pointer,
        "ordinary_Z_orders": list(ORDERS),
        "ordinary_Z0_to_Z2": all_rows[: len(ORDERS)],
        "all_saved_rows": all_rows,
        "source_values_selected": False,
    }


def _first_row(rows: Mapping[str, Any], order: int) -> dict[str, Any]:
    if rows.get("status") != "direct_source_rows":
        return {
            "status": "unresolved",
            "ordinary_Z_order": order,
            "source_value_selected": False,
        }
    return dict(rows["ordinary_Z0_to_Z2"][order], ordinary_Z_order=order)


def _at_xi_11(function_ref: Mapping[str, Any], coordinate_ref: Mapping[str, Any]) -> dict[str, Any]:
    """Represent evaluation at xi=11 without creating a graph node."""

    if function_ref.get("status") != "direct_source_ref" or coordinate_ref.get("status") != "direct_source_ref":
        return {
            "status": "unresolved",
            "operator": "source_function_at",
            "coordinate_literal": "xi=11",
            "source_value_selected": False,
        }
    return {
        "status": "conditional_source_expression",
        "operator": "source_function_at",
        "function_ref": dict(function_ref),
        "coordinate_owner": dict(coordinate_ref),
        "coordinate_literal": "xi=11",
        "source_value_selected": False,
    }


def _formal_symbol(name: str) -> dict[str, Any]:
    return {"kind": "formal_symbol", "name": name}


def _literal(value: str) -> dict[str, Any]:
    return {"kind": "exact_scalar_literal", "value": value}


def _expr(operator: str, *arguments: Any) -> dict[str, Any]:
    return {"operator": operator, "arguments": list(arguments)}


def _side_snapshot(
    side: str, container: Mapping[str, Any], pointer: str
) -> tuple[dict[str, Any], list[str]]:
    gates: list[str] = []

    raw_geometry = _optional(container, ("geometry",))
    geometry: dict[str, Any] = {}
    if not isinstance(raw_geometry, Mapping):
        gates.append(f"{side}:geometry_owner_missing")
        raw_geometry = {}
    for name in ("mu", "delta", "Pstar", "N", "Rp"):
        if name in raw_geometry:
            geometry[name] = _node_ref(
                raw_geometry[name],
                f"{pointer}.geometry.{name}",
                f"{side}.geometry.{name}",
                gates,
            )
        else:
            gates.append(f"{side}:geometry.{name}:direct_source_owner_missing")
            geometry[name] = {
                "status": "unresolved",
                "source_role": f"{pointer}.geometry.{name}",
                "node_id_is_not_a_numeric_value": True,
            }

    coordinate_key = "coordinates" if side == "P" else "pulse_coordinates"
    coordinates = _optional(container, (coordinate_key, "gap"))
    coordinate_pointer = f"{pointer}.{coordinate_key}.gap"
    if not isinstance(coordinates, Mapping):
        gates.append(f"{side}:gap_coordinate_owner_missing")
        coordinates = {}
    coordinate_owner: dict[str, Any] = {}
    for name in ("query", "y", "radius"):
        if side == "F" and name == "radius":
            metadata = _node_ref(
                coordinates.get(name), f"{coordinate_pointer}.{name}",
                "F.gap.current_radius_metadata", gates,
            )
            coordinate_owner[name] = {
                "status": "unresolved",
                "reason": "frame pulse_coordinates radius is current-geometry metadata; its identity with the radius in the frame-only field graph is not established",
                "current_radius_metadata": metadata,
                "source_value_selected": False,
            }
            gates.append("F:frame_only_field_radius_owner_unresolved")
            continue
        if name in coordinates:
            coordinate_owner[name] = _node_ref(
                coordinates[name],
                f"{coordinate_pointer}.{name}",
                f"{side}.gap.{name}",
                gates,
            )
        else:
            gates.append(f"{side}:gap.{name}:direct_source_owner_missing")
            coordinate_owner[name] = {
                "status": "unresolved",
                "source_role": f"{coordinate_pointer}.{name}",
                "node_id_is_not_a_numeric_value": True,
            }
    coordinate_owner["coordinate_literal"] = coordinates.get("coordinate")
    if coordinate_owner["coordinate_literal"] != "xi=mu*log(R/Rp)":
        gates.append(f"{side}:gap_coordinate_literal_missing_or_changed")
    coordinate_owner["source_value_selected"] = False

    selected = _optional(container, ("selected_pulse",))
    selected_refs: dict[str, Any] = {}
    if not isinstance(selected, Mapping):
        gates.append(f"{side}:selected_pulse_owner_missing")
        selected = {}
    for name in ("C1", "C2"):
        selected_refs[name] = _row_refs(
            selected.get(name),
            f"{pointer}.selected_pulse.{name}",
            f"{side}.selected_pulse.{name}",
            gates,
        )

    leading = _optional(container, ("leading_fields", "gap"))
    if not isinstance(leading, Mapping):
        gates.append(f"{side}:leading_fields.gap_owner_missing")
        leading = {}
    leading_refs: dict[str, Any] = {}
    for name in ("U", "U_y", "m", "m_y", "Q", "Q_y"):
        leading_refs[name] = _row_refs(
            leading.get(name),
            f"{pointer}.leading_fields.gap.{name}",
            f"{side}.leading_fields.gap.{name}",
            gates,
        )

    radius_at_11 = _at_xi_11(coordinate_owner["radius"], coordinate_owner["query"])
    q_rows = leading_refs["Q"]
    qy_rows = leading_refs["Q_y"]
    q_at_11 = _at_xi_11(_first_row(q_rows, 0), coordinate_owner["query"])
    qz_at_11 = _at_xi_11(_first_row(q_rows, 1), coordinate_owner["query"])
    qzz_at_11 = _at_xi_11(_first_row(q_rows, 2), coordinate_owner["query"])
    qy_at_11 = _at_xi_11(_first_row(qy_rows, 0), coordinate_owner["query"])
    qyz_at_11 = _at_xi_11(_first_row(qy_rows, 1), coordinate_owner["query"])

    endpoints = _optional(container, ("endpoints", "gap"))
    if isinstance(endpoints, list) and len(endpoints) >= 2:
        endpoint_fields = ("R_q", "R_squared_q", "R_q_squared", "R_squared_q_squared", "R_U_squared", "R_squared_U_squared", "q_plus_q_y", "R_q_y")
        endpoint_roles: dict[str, Any] = {
            "status": "direct_endpoint_function_roles",
            "lower": {
                "source_role": f"{pointer}.endpoints.gap[0]",
                "fields": {name: _row_refs(endpoints[0].get(name), f"{pointer}.endpoints.gap[0].{name}", f"{side}.gap.lower.{name}", gates) for name in endpoint_fields} if isinstance(endpoints[0], Mapping) else {},
                "source_value_selected": False,
            },
            "upper": {
                "source_role": f"{pointer}.endpoints.gap[1]",
                "fields": {name: _row_refs(endpoints[1].get(name), f"{pointer}.endpoints.gap[1].{name}", f"{side}.gap.upper.{name}", gates) for name in endpoint_fields} if isinstance(endpoints[1], Mapping) else {},
                "source_value_selected": False,
            },
        }
    else:
        endpoint_roles = {
            "status": "unresolved_missing_direct_endpoint_roles",
            "source_value_selected": False,
        }
        gates.append(f"{side}:gap_finite_endpoint_radius_owner_missing")

    return {
        "side": side,
        "source_container_role": pointer,
        "geometry": geometry,
        "coordinate_owner": coordinate_owner,
        "selected_pulse": selected_refs,
        "leading_fields_gap": leading_refs,
        "finite_endpoints": endpoint_roles,
        "H_source_expressions_at_xi_11": {
            "H": _expr("multiply", radius_at_11, q_at_11),
            "H_Z": _expr("multiply", radius_at_11, qz_at_11),
            "H_ZZ": _expr("multiply", radius_at_11, qzz_at_11),
        },
        "Q_y_source_expression_at_xi_11": qy_at_11,
        "Q_yZ_source_expression_at_xi_11": qyz_at_11,
        "gap_prerequisites": {
            "U_equals_zero": {
                "status": "caller_condition_not_machine_proved",
                "source_rows": leading_refs["U"],
            },
            "U_y_equals_zero": {
                "status": "caller_condition_not_machine_proved",
                "source_rows": leading_refs["U_y"],
            },
            "m_y_equals_negative_m": {
                "status": "caller_condition_not_machine_proved",
                "source_rows": {
                    "m": leading_refs["m"],
                    "m_y": leading_refs["m_y"],
                },
            },
            "Q_y_equals_negative_Q": {
                "status": "caller_condition_not_machine_proved",
                "source_rows": {"Q": q_rows, "Q_y": qy_rows},
            },
            "R_Q_equals_H": {
                "status": "conditional_on_Q_y_equals_negative_Q",
                "source_rows": {"Q": q_rows, "Q_y": qy_rows},
            },
            "radius_Z_independent": {
                "status": "caller_condition_not_machine_proved",
                "source_owner": coordinate_owner["radius"],
            },
        },
        "source_graph_traversed": False,
        "source_values_selected": False,
    }, gates


def _unresolved_input(side: str, quantity: str, endpoint_roles: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "status": "caller_supplied_source_difference_unresolved",
        "kind": "source_endpoint_difference",
        "side": side,
        "quantity": quantity,
        "endpoint_owner": dict(endpoint_roles),
        "source_value_selected": False,
        "numeric_difference_evaluated": False,
    }


def _formal_difference_formulas(
    physical: Mapping[str, Any], frame: Mapping[str, Any]
) -> dict[str, Any]:
    """Build formal J/J_Z trees while refusing to merge side metadata."""

    z = _formal_symbol("Z")
    delta = {
        "kind": "conditional_common_source_handle",
        "name": "delta_common",
        "source_value_selected": False,
    }
    one, two, four = (_literal(value) for value in ("1", "2", "4"))
    L = _expr("subtract", one, _expr("multiply", delta, _expr("power", z, two)))
    h = _expr("subtract", one, delta)
    c = _expr("divide", _expr("multiply", h, z), _expr("multiply", four, L))
    c_z = _expr(
        "divide",
        _expr("multiply", h, _expr("add", one, _expr("multiply", delta, _expr("power", z, two)))),
        _expr("multiply", four, _expr("power", L, two)),
    )

    hp = physical["H_source_expressions_at_xi_11"]
    hf = frame["H_source_expressions_at_xi_11"]
    dr_p = _unresolved_input("P", "Delta_R", physical["finite_endpoints"])
    dr_f = _unresolved_input("F", "Delta_R", frame["finite_endpoints"])
    drz_p = _unresolved_input("P", "Delta_R_Z", physical["finite_endpoints"])
    drz_f = _unresolved_input("F", "Delta_R_Z", frame["finite_endpoints"])
    dy_p = _unresolved_input("P", "Delta_y", physical["finite_endpoints"])
    dy_f = _unresolved_input("F", "Delta_y", frame["finite_endpoints"])
    dyz_p = _unresolved_input("P", "Delta_y_Z", physical["finite_endpoints"])
    dyz_f = _unresolved_input("F", "Delta_y_Z", frame["finite_endpoints"])
    dy = {
        "kind": "conditional_common_source_difference",
        "name": "Delta_y_common",
        "side_inputs": {"P": dy_p, "F": dy_f},
        "same_coordinate_and_mu_proved": False,
        "source_value_selected": False,
    }
    dy_z = {
        "kind": "conditional_common_source_difference",
        "name": "Delta_y_Z_common",
        "side_inputs": {"P": dyz_p, "F": dyz_f},
        "same_coordinate_and_mu_proved": False,
        "source_value_selected": False,
    }
    linear = _expr(
        "subtract",
        _expr("multiply", hp["H_Z"], dr_p),
        _expr("multiply", hf["H_Z"], dr_f),
    )
    linear_z = _expr(
        "subtract",
        _expr(
            "add",
            _expr("multiply", hp["H_ZZ"], dr_p),
            _expr("multiply", hp["H_Z"], drz_p),
        ),
        _expr(
            "add",
            _expr("multiply", hf["H_ZZ"], dr_f),
            _expr("multiply", hf["H_Z"], drz_f),
        ),
    )
    quadratic = _expr(
        "subtract", _expr("power", hp["H"], two), _expr("power", hf["H"], two)
    )
    quadratic_z = _expr(
        "subtract",
        _expr("multiply", two, _expr("multiply", hp["H"], hp["H_Z"])),
        _expr("multiply", two, _expr("multiply", hf["H"], hf["H_Z"])),
    )
    J = _expr(
        "add",
        _expr("multiply", c, linear),
        _expr("negate", _expr("multiply", _expr("divide", dy, four), quadratic)),
    )
    J_Z = _expr(
        "add",
        _expr("multiply", c_z, linear),
        _expr("multiply", c, linear_z),
        _expr("negate", _expr("multiply", _expr("divide", dy_z, four), quadratic)),
        _expr("negate", _expr("multiply", _expr("divide", dy, four), quadratic_z)),
    )
    def separate_side(side, snapshot, H, dR, dR_Z, dY, dY_Z):
        # These parameters are side-specific until a common source is proved.
        side_delta = {
            "kind": "source_parameter_reference",
            "side": side,
            "owner": snapshot["geometry"]["delta"],
            "source_value_selected": False,
        }
        side_L = _expr("subtract", one, _expr("multiply", side_delta, _expr("power", z, two)))
        side_h = _expr("subtract", one, side_delta)
        side_c = _expr("divide", _expr("multiply", side_h, z), _expr("multiply", four, side_L))
        side_c_Z = _expr("divide", _expr("multiply", side_h, _expr("add", one, _expr("multiply", side_delta, _expr("power", z, two)))), _expr("multiply", four, _expr("power", side_L, two)))
        linear_side = _expr("multiply", H["H_Z"], dR)
        linear_side_Z = _expr("add", _expr("multiply", H["H_ZZ"], dR), _expr("multiply", H["H_Z"], dR_Z))
        square = _expr("power", H["H"], two)
        moment = _expr("subtract", _expr("multiply", side_c, linear_side), _expr("multiply", _expr("divide", dY, four), square))
        moment_Z = _expr("add", _expr("multiply", side_c_Z, linear_side), _expr("multiply", side_c, linear_side_Z), _expr("negate", _expr("multiply", _expr("divide", dY_Z, four), square)), _expr("negate", _expr("multiply", _expr("divide", dY, two), _expr("multiply", H["H"], H["H_Z"]))))
        return {"delta": side_delta, "L": side_L, "c": side_c, "moment": moment, "ordinary_Z1": moment_Z, "delta_Z_zero_required": True, "source_value_selected": False}
    side_p = separate_side("P", physical, hp, dr_p, drz_p, dy_p, dyz_p)
    side_f = separate_side("F", frame, hf, dr_f, drz_f, dy_f, dyz_f)
    return {
        "status": "conditional_formal_expression_not_evaluated",
        "common_geometry": {
            "Z": z,
            "delta": delta,
            "L": L,
            "c": c,
            "c_Z": c_z,
            "delta_Z_zero_required": True,
            "same_mu_delta_coordinate_proved": False,
        },
        "side_delta_R": {"P": dr_p, "F": dr_f},
        "side_delta_R_Z": {"P": drz_p, "F": drz_f},
        "side_delta_y": {"P": dy_p, "F": dy_f},
        "side_delta_y_Z": {"P": dyz_p, "F": dyz_f},
        "side_terms_kept_separate": {"P": side_p, "F": side_f},
        "unmerged_gap_difference": _expr("subtract", side_p["moment"], side_f["moment"]),
        "unmerged_gap_difference_Z1": _expr("subtract", side_p["ordinary_Z1"], side_f["ordinary_Z1"]),
        "conditional_common_J": J,
        "conditional_common_J_Z": J_Z,
        "merge_gate": {
            "status": "unresolved",
            "requires": [
                "physical_and_frame_share_the_same_mu_handle",
                "physical_and_frame_share_the_same_delta_handle",
                "L_is_nonzero_on_the_common_domain",
                "both_sides_share_the_same_gap_coordinate_and_two_finite_endpoints",
                "both_side_radius_functions_are_Z_independent_for_H_Z_and_H_ZZ",
                "caller_supplies_source_owned_Delta_R_and_Delta_y_jets",
                "frame_only_field_radius_owner_is_source_bound",
                "gap_U_and_U_y_are_zero_and_Q_y_equals_negative_Q_on_each_side",
                "ordinary_Z1_Q_yZ_and_source_field_derivatives_are_valid",
                "physical_measure_and_finite_endpoint_ownership_are_bound_on_each_side",
                "physical_Pstar_and_units_are_matched_once_on_each_side",
            ],
            "source_owner_equality_proved": False,
        },
    }


def build_actual_gap_source_program(repo_root: str | Path) -> dict[str, Any]:
    """Return a JSON-serializable current physical/frame-only source packet."""

    root = Path(repo_root)
    physical_raw, physical_receipt = _read_json_after_raw_auth(
        root, PHYSICAL_REPORT, RAW_SHA256["physical"]
    )
    frame_raw, frame_receipt = _read_gzip_json_after_raw_auth(
        root, FRAME_ARCHIVE, RAW_SHA256["frame"]
    )
    physical = _direct(
        physical_raw,
        ("component_result", "physical_pulse_bindings"),
        "physical.component_result.physical_pulse_bindings",
    )
    frame = _direct(
        frame_raw,
        ("Rp_native_C5_pulse_source_bindings",),
        "frame.Rp_native_C5_pulse_source_bindings",
    )
    _require(isinstance(physical, Mapping), "physical bindings are not an object")
    _require(isinstance(frame, Mapping), "frame bindings are not an object")

    p_snapshot, p_gates = _side_snapshot(
        "P", physical, "$.component_result.physical_pulse_bindings"
    )
    f_snapshot, f_gates = _side_snapshot(
        "F", frame, "$.Rp_native_C5_pulse_source_bindings"
    )
    formulas = _formal_difference_formulas(p_snapshot, f_snapshot)

    physical_correction = _optional(physical, ("pulse_geometry_correction_m15",))
    if isinstance(physical_correction, list):
        correction_roles: dict[str, Any] = {
            "status": "direct_physical_source_correction_roles_retained",
            "source_role": "$.component_result.physical_pulse_bindings.pulse_geometry_correction_m15",
            "row_count": len(physical_correction),
            "source_node_ids": list(physical_correction),
            "source_node_ids_are_not_values": True,
            "source_values_selected": False,
            "correction_removed": False,
        }
    else:
        correction_roles = {
            "status": "unresolved_missing_physical_correction_roles",
            "source_values_selected": False,
            "correction_removed": False,
        }
        p_gates.append("P:physical_correction_roles_missing")

    return {
        "schema": SCHEMA,
        "status": "incomplete_source_boundary_not_admissible_for_numeric_use",
        "admissible_for_numeric_use": False,
        "source_gap_identity_proved": False,
        "source_only": True,
        "inputs": {
            "physical_report": physical_receipt,
            "frame_archive": frame_receipt,
            "frame_archive_decoded_only_after_raw_hash": True,
        },
        "sides": {"P": p_snapshot, "F": f_snapshot},
        "gap_coordinate_contract": {
            "status": "proposed_literals_not_a_source_identity",
            "P_coordinate_owner": p_snapshot["coordinate_owner"],
            "F_coordinate_owner": f_snapshot["coordinate_owner"],
            "coordinate_literal": "xi=mu*log(R/Rp)",
            "finite_gap_endpoint_literals": {
                "lower": "xi=11",
                "upper": "xi=13-4*mu",
                "end_anchor": "s=-4",
            },
            "literal_values_not_evaluated": True,
            "frame_radius_metadata_is_not_field_radius_proof": True,
        },
        "formal_difference": formulas,
        "independent_reference_offset": {
            "status": "unresolved_symbolic_marker",
            "selected": False,
            "source_role": "independent_reference_offset",
        },
        "retained_212_scope": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "six_signed_channels_retained": True,
            "reserved_contribution_registered_once": True,
            "old_errors_removed": False,
            "source_values_selected": False,
        },
        "physical_pulse_correction": correction_roles,
        "unresolved_gates": sorted(set(p_gates + f_gates)),
        "flags": {
            "numeric_values_computed": False,
            "numeric_difference_evaluated": False,
            "source_function_values_selected": False,
            "source_graph_traversed": False,
            "old_source_constructor_imported": False,
            "old_source_checker_imported": False,
            "incident_pair_DAG_reentered": False,
            "root_mapping_performed": False,
            "rewritten_node_mapping_consumed": False,
            "physical_frame_source_equality_proved": False,
            "common_mu_delta_L_proved": False,
            "common_coordinate_proved": False,
            "ordinary_Z0_to_Z2_source_refs_preserved": all(
                snapshot["leading_fields_gap"][field].get("status") == "direct_source_rows"
                for snapshot in (p_snapshot, f_snapshot) for field in ("Q", "Q_y")
            ),
            "Q_y_equals_negative_Q_machine_proof": False,
            "H_equals_RQ_machine_proof": False,
            "physical_pulse_correction_removed": False,
            "reference_offset_selected": False,
            "error_terms_removed": False,
        },
        "scope_audit": {
            "direct_literal_field_retrieval_only": True,
            "exact_graph_nodes_subtree_accessed": False,
            "complete_frame_archive_json_decoded": True,
            "rewritten_node_mapping_read": False,
            "source_constructors_called": 0,
            "source_checkers_called": 0,
            "graph_comparisons": 0,
            "quadratures": 0,
            "mathematical_evaluations": 0,
            "gzip_archives_decoded": 1,
            "gzip_raw_hashes_authenticated_before_decode": True,
        },
    }


# No ``__main__`` entry point: root freezes and invokes this leaf exactly once.
