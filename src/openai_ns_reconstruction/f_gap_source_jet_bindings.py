"""Bound opaque F-gap source handles for the seven-jet adapter.

This leaf consumes only the already-published
``pulse_gap_actual_source_program_v1/result.json`` packet.  The packet is
authenticated by raw SHA-256 before its JSON is parsed.  The selected rows are
source handles, not values: a saved integer node id is retained as an opaque
identifier and is never used as arithmetic input, graph index, or derivative
value.

The requested adapter input order is the twelve field handles
``U[0:3], m[0:4], U_y[0:2], m_y[0:3]`` followed by
``delta0, delta1, delta2, Z``.  Only ``delta0`` can use the producer's
opaque ``geometry.delta`` candidate; the other three adapter operands remain
unresolved because this packet does not authenticate their owners.  The
additional geometry candidates ``mu, Pstar, Rp`` remain in the geometry
section and are not adapter operands.  The producer packet labels only ordinary Z0--Z2;
``m[3]`` is admitted only as the fourth item of ``all_saved_rows`` and keeps
``ordinary_Z3_semantics_proved=False``.  Coordinate and radius metadata are
retained separately.  In particular, the radius metadata is never promoted to
the F field-radius owner.

No source constructor, source checker, graph edge, archive, numerical value,
or old leaf is imported or called here.  There is intentionally no executable
entry point; a coordinating runner may call the two public functions.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA = "f_gap_source_jet_bindings_v1"
PRIOR_COMPONENT = "pulse_gap_actual_source_program_v1"
RESULT_PATH = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
RESULT_SHA256 = "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"

F_SOURCE_NAMESPACE = "frame_Rp_native_C5_archive"
F_OWNER = "F"
F_CONTAINER_POINTER = "$.Rp_native_C5_pulse_source_bindings"
F_LEADING_POINTER = F_CONTAINER_POINTER + ".leading_fields.gap"
F_GEOMETRY_POINTER = F_CONTAINER_POINTER + ".geometry"
F_COORDINATE_POINTER = F_CONTAINER_POINTER + ".pulse_coordinates.gap"

FIELD_ORDERS: dict[str, tuple[int, ...]] = {
    "U": (0, 1, 2),
    "m": (0, 1, 2, 3),
    "U_y": (0, 1),
    "m_y": (0, 1, 2),
}
GEOMETRY_NAMES = ("mu", "delta", "Pstar", "Rp")
COORDINATE_NAMES = ("query", "y", "radius")
INPUT_ORDER = (
    "U0",
    "U1",
    "U2",
    "m0",
    "m1",
    "m2",
    "m3",
    "Uy0",
    "Uy1",
    "my0",
    "my1",
    "my2",
    "delta0",
    "delta1",
    "delta2",
    "Z",
)


class FGapSourceJetBindingInputError(ValueError):
    """Raised when the frozen raw-input authentication contract is invalid."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapSourceJetBindingInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    result = value.lower()
    _require(all(char in "0123456789abcdef" for char in result), label + " SHA malformed")
    return result


def _normalise_key(raw_path: Any, root: Path, label: str) -> str:
    """Resolve a manifest path without allowing a workspace escape."""

    _require(isinstance(raw_path, str) and raw_path, label + " path missing")
    normalised = raw_path.replace("\\", "/")
    candidate = Path(normalised)
    root_resolved = root.resolve()
    if candidate.is_absolute():
        absolute = candidate.resolve()
        try:
            relative = absolute.relative_to(root_resolved)
        except ValueError:
            suffix = "/" + RESULT_PATH.as_posix()
            _require(
                normalised.endswith(suffix) or normalised == RESULT_PATH.as_posix(),
                label + " absolute path is outside the workspace",
            )
            relative = RESULT_PATH
    else:
        relative = candidate
    _require(".." not in relative.parts, label + " path escapes the workspace")
    key = relative.as_posix()
    _require(key, label + " path is empty")
    resolved = (root_resolved / relative).resolve()
    _require(resolved.is_relative_to(root_resolved), label + " path escapes the workspace")
    return key


def _manifest_entries(
    source: Mapping[str, Any],
    hashes: Mapping[str, Any] | None,
    root: Path,
) -> dict[str, str]:
    """Collect all raw input hashes without parsing any input bytes."""

    entries: dict[str, str] = {}

    def add(raw_path: Any, expected: Any, label: str) -> None:
        key = _normalise_key(raw_path, root, label)
        value = _digest(expected, label)
        prior = entries.get(key)
        _require(prior is None or prior == value, "manifest hash conflict: " + key)
        entries[key] = value

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            add(item.get("path"), item.get("sha256"), f"source.files[{index}]")
    for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for raw_path, expected in table.items():
                add(raw_path, expected, f"source.{field}")
    if isinstance(hashes, Mapping):
        for raw_path, expected in hashes.items():
            add(raw_path, expected, f"caller hashes.{raw_path}")

    _require(entries, "caller frozen manifest contains no raw SHA entries")
    result_key = RESULT_PATH.as_posix()
    _require(result_key in entries, "published result SHA is not pinned")
    _require(entries[result_key] == RESULT_SHA256, "published result SHA pin changed")
    return entries


def _authenticate_raw_inputs(
    root: Path,
    source: Mapping[str, Any],
    hashes: Mapping[str, Any] | None,
) -> tuple[bytes, dict[str, Any]]:
    """Hash every supplied raw input before decoding the published result."""

    entries = _manifest_entries(source, hashes, root)
    actual: dict[str, str] = {}
    result_raw: bytes | None = None
    for key, expected in entries.items():
        path = (root / key).resolve()
        _require(path.is_relative_to(root), "authenticated path escapes the workspace: " + key)
        _require(path.is_file(), "authenticated input is missing: " + key)
        raw = path.read_bytes()
        observed = _sha(raw)
        _require(observed == expected, "raw SHA changed: " + key)
        actual[key] = observed
        if key == RESULT_PATH.as_posix():
            result_raw = raw

    _require(result_raw is not None, "published result bytes were not authenticated")
    return result_raw, {
        "all_manifest_raw_hashes_checked_before_parse": True,
        "manifest_entry_count": len(entries),
        "required_result": {
            "path": RESULT_PATH.as_posix(),
            "expected_sha256": RESULT_SHA256,
            "actual_sha256": actual[RESULT_PATH.as_posix()],
        },
        "authenticated_inputs": {
            key: {"expected_sha256": entries[key], "actual_sha256": actual[key]}
            for key in sorted(entries)
        },
        "parsed_paths": [RESULT_PATH.as_posix()],
        "other_inputs_parsed": False,
    }


def _json_object(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapSourceJetBindingInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _validated_prior_inner(wrapper: Mapping[str, Any]) -> tuple[Mapping[str, Any] | None, list[str]]:
    """Validate the frozen outer/inner wrapper contract before field access."""

    errors: list[str] = []
    if wrapper.get("component") != PRIOR_COMPONENT:
        errors.append("outer component is not pulse_gap_actual_source_program_v1")
    if wrapper.get("status") != "returned":
        errors.append("outer status is not returned")
    if "schema" in wrapper:
        errors.append("outer wrapper unexpectedly contains schema")
    if "error" in wrapper:
        errors.append("outer wrapper unexpectedly contains error")
    inner = wrapper.get("result")
    if not isinstance(inner, Mapping):
        errors.append("outer result is missing or not an object")
        return None, errors
    if inner.get("schema") != PRIOR_COMPONENT:
        errors.append("inner result schema does not match the producer component")
    if inner.get("status") != "incomplete_source_boundary_not_admissible_for_numeric_use":
        errors.append("inner result status is not the frozen incomplete-source status")
    if inner.get("source_only") is not True:
        errors.append("inner result source_only is not true")
    if inner.get("admissible_for_numeric_use") is not False:
        errors.append("inner result admissible_for_numeric_use is not false")
    return (inner if not errors else None), errors


def _unresolved(
    expected_role: str,
    reason: str,
    *,
    kind: str,
    field: str | None = None,
    order: int | None = None,
    original: Any = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "unresolved",
        "kind": kind,
        "reason": reason,
        "expected_source_role": expected_role,
        "source_value_selected": False,
        "node_id_is_not_a_numeric_value": True,
        "ordinary_Z_semantics_proved": False,
        "ordinary_Z3_semantics_proved": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
    }
    if field is not None:
        result["field"] = field
    if order is not None:
        result["requested_ordinary_Z_order"] = order
    if original is not None:
        result["original_source_ref"] = copy.deepcopy(original)
    return result


def _opaque_handle(
    raw: Any,
    expected_role: str,
    *,
    kind: str,
    field: str | None = None,
    order: int | None = None,
) -> dict[str, Any]:
    """Bind one direct ref while keeping its node id opaque."""

    if not isinstance(raw, Mapping):
        return _unresolved(
            expected_role,
            "direct source ref is missing or not an object",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    if raw.get("status") != "direct_source_ref":
        return _unresolved(
            expected_role,
            "direct source ref status is not direct_source_ref",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    node = raw.get("source_node_id")
    role = raw.get("source_role")
    namespace = raw.get("namespace")
    owner = raw.get("owner")
    if type(node) is not int:
        return _unresolved(
            expected_role,
            "source_node_id is absent or not an opaque integer handle",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    if role != expected_role:
        return _unresolved(
            expected_role,
            "source role does not match the selected F path",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    if namespace is not None and namespace != F_SOURCE_NAMESPACE:
        result = _unresolved(
            expected_role,
            "source namespace differs from the F archive candidate namespace",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
        result["observed_namespace"] = namespace
        result["expected_namespace"] = F_SOURCE_NAMESPACE
        return result
    if owner is not None and owner != F_OWNER:
        result = _unresolved(
            expected_role,
            "source owner differs from F",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
        result["observed_owner"] = owner
        return result
    if raw.get("node_id_is_not_a_numeric_value") is not True:
        return _unresolved(
            expected_role,
            "node id lacks the required opaque-identity guard",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    if raw.get("source_value_selected") is not False:
        return _unresolved(
            expected_role,
            "source value was selected by the input packet",
            kind=kind,
            field=field,
            order=order,
            original=raw,
        )
    return {
        "status": "bound_opaque_source_handle",
        "kind": kind,
        "field": field,
        "requested_ordinary_Z_order": order,
        "source_role": role,
        "source_node_id": node,
        "namespace": F_SOURCE_NAMESPACE,
        "namespace_observed": namespace,
        "namespace_inferred": namespace is None,
        "namespace_verified": False,
        "owner": F_OWNER,
        "owner_inferred": owner is None,
        "owner_verified": False,
        "source_ref": copy.deepcopy(raw),
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
        "ordinary_Z_semantics_proved": False,
        "ordinary_Z3_semantics_proved": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
    }


def _ref_identity(value: Any) -> tuple[Any, ...] | None:
    if not isinstance(value, Mapping):
        return None
    node = value.get("source_node_id")
    role = value.get("source_role")
    if type(node) is not int or not isinstance(role, str):
        return None
    return (node, role, value.get("namespace"), value.get("owner"))


def _bind_field_rows(
    leading: Mapping[str, Any] | None,
    field: str,
) -> dict[str, Any]:
    """Select requested rows, with m3 restricted to all_saved_rows[3]."""

    expected_record_role = F_LEADING_POINTER + "." + field
    orders = FIELD_ORDERS[field]
    record = leading.get(field) if isinstance(leading, Mapping) else None
    if not isinstance(record, Mapping):
        return {
            "status": "unresolved_field_record",
            "field": field,
            "source_role": expected_record_role,
            "ordinary_Z_orders_requested": list(orders),
            "rows": {
                str(order): _unresolved(
                    f"{expected_record_role}[{order}]",
                    "leading field record is missing",
                    kind="field_jet_handle",
                    field=field,
                    order=order,
                )
                for order in orders
            },
            "source_values_selected": False,
        }

    record_contract_errors: list[str] = []
    if record.get("status") != "direct_source_rows":
        record_contract_errors.append("record status is not direct_source_rows")
    if record.get("source_role") != expected_record_role:
        record_contract_errors.append("record source_role does not match the expected F pointer")
    if record.get("ordinary_Z_orders") != [0, 1, 2]:
        record_contract_errors.append("record ordinary_Z_orders is not exactly [0, 1, 2]")
    if record.get("source_values_selected") is not False:
        record_contract_errors.append("record source_values_selected is not false")
    if record_contract_errors:
        return {
            "status": "unresolved_field_record_contract",
            "field": field,
            "source_role": expected_record_role,
            "ordinary_Z_orders_requested": list(orders),
            "rows": {
                str(order): _unresolved(
                    f"{expected_record_role}[{order}]",
                    "; ".join(record_contract_errors),
                    kind="field_jet_handle",
                    field=field,
                    order=order,
                    original=record,
                )
                for order in orders
            },
            "record_contract_errors": record_contract_errors,
            "source_values_selected": False,
            "ordinary_Z_semantics_proved": False,
            "ordinary_Z3_semantics_proved": False,
        }

    z02 = record.get("ordinary_Z0_to_Z2")
    all_rows = record.get("all_saved_rows")
    row_map: dict[str, Any] = {}
    for order in orders:
        expected_role = f"{expected_record_role}[{order}]"
        selected: Any = None
        selected_from = ""
        if order <= 2 and isinstance(z02, list) and len(z02) > order:
            selected = z02[order]
            selected_from = f"ordinary_Z0_to_Z2[{order}]"
            if isinstance(all_rows, list) and len(all_rows) > order:
                ordinary_identity = _ref_identity(selected)
                all_identity = _ref_identity(all_rows[order])
                if ordinary_identity is not None and all_identity is not None and ordinary_identity != all_identity:
                    row_map[str(order)] = _unresolved(
                        expected_role,
                        "ordinary_Z0_to_Z2 and all_saved_rows disagree for the same order",
                        kind="field_jet_handle",
                        field=field,
                        order=order,
                        original={
                            "ordinary_Z0_to_Z2": selected,
                            "all_saved_rows": all_rows[order],
                        },
                    )
                    continue
        elif isinstance(all_rows, list) and len(all_rows) > order:
            selected = all_rows[order]
            selected_from = f"all_saved_rows[{order}]"
        else:
            row_map[str(order)] = _unresolved(
                expected_role,
                "requested source row is absent",
                kind="field_jet_handle",
                field=field,
                order=order,
            )
            continue

        binding = _opaque_handle(
            selected,
            expected_role,
            kind="field_jet_handle",
            field=field,
            order=order,
        )
        binding["record_source_role"] = record.get("source_role", expected_record_role)
        binding["selected_from"] = selected_from
        binding["ordinary_Z3_semantics_proved"] = False
        if field == "m" and order == 3:
            binding["extra_row_semantics_promoted"] = False
            binding["ordinary_Z3_semantics_proved"] = False
            binding["source_row_is_only_fourth_all_saved_row"] = True
        else:
            binding["source_row_is_only_fourth_all_saved_row"] = False
        row_map[str(order)] = binding

    return {
        "status": "direct_source_rows_selected" if all(
            item.get("status") == "bound_opaque_source_handle" for item in row_map.values()
        ) else "partially_unresolved_source_rows",
        "field": field,
        "source_role": record.get("source_role", expected_record_role),
        "ordinary_Z_orders_requested": list(orders),
        "producer_record_status": record.get("status"),
        "rows": row_map,
        "source_values_selected": False,
        "ordinary_Z_semantics_proved": False,
        "ordinary_Z3_semantics_proved": False,
    }


def _bind_geometry(f_side: Mapping[str, Any] | None) -> dict[str, Any]:
    geometry = f_side.get("geometry") if isinstance(f_side, Mapping) else None
    result: dict[str, Any] = {}
    for name in GEOMETRY_NAMES:
        role = f"{F_GEOMETRY_POINTER}.{name}"
        raw = geometry.get(name) if isinstance(geometry, Mapping) else None
        result[name] = _opaque_handle(raw, role, kind="geometry_handle")
        result[name]["geometry_name"] = name
    return result


def _bind_coordinate_metadata(f_side: Mapping[str, Any] | None) -> dict[str, Any]:
    owner = f_side.get("coordinate_owner") if isinstance(f_side, Mapping) else None
    result: dict[str, Any] = {}
    for name in COORDINATE_NAMES:
        role = f"{F_COORDINATE_POINTER}.{name}"
        raw = owner.get(name) if isinstance(owner, Mapping) else None
        if name != "radius":
            binding = _opaque_handle(raw, role, kind="gap_coordinate_metadata")
            binding["metadata_only"] = True
            result[name] = binding
            continue

        metadata = raw.get("current_radius_metadata") if isinstance(raw, Mapping) else raw
        binding = _opaque_handle(metadata, role, kind="gap_radius_metadata")
        result[name] = {
            "status": "metadata_only_unresolved_field_owner"
            if binding.get("status") == "bound_opaque_source_handle"
            else "unresolved",
            "kind": "gap_radius_metadata",
            "metadata_ref": binding,
            "metadata_radius_not_field_owner": True,
            "field_radius_owner_proved": False,
            "source_value_selected": False,
            "source_capture_proved": False,
            "physical_scale_proved": False,
            "reason": (
                "pulse_coordinates.gap.radius is retained as current-geometry "
                "metadata; it is not promoted to the frame-only field-radius owner"
            ),
        }
    literal = owner.get("coordinate_literal") if isinstance(owner, Mapping) else None
    result["coordinate_literal"] = literal if isinstance(literal, str) else None
    result["coordinate_literal_is_value"] = False
    return result


def _adapter_inputs(
    field_bindings: Mapping[str, Any], geometry: Mapping[str, Any]
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    field_names = {
        "U0": ("U", 0),
        "U1": ("U", 1),
        "U2": ("U", 2),
        "m0": ("m", 0),
        "m1": ("m", 1),
        "m2": ("m", 2),
        "m3": ("m", 3),
        "Uy0": ("U_y", 0),
        "Uy1": ("U_y", 1),
        "my0": ("m_y", 0),
        "my1": ("m_y", 1),
        "my2": ("m_y", 2),
    }
    for name in INPUT_ORDER:
        if name in field_names:
            field, order = field_names[name]
            record = field_bindings.get(field, {})
            binding = record.get("rows", {}).get(str(order))
            if not isinstance(binding, Mapping):
                binding = _unresolved(
                    f"{F_LEADING_POINTER}.{field}[{order}]",
                    "selected field binding is missing",
                    kind="field_jet_handle",
                    field=field,
                    order=order,
                )
            result.append(
                {
                    "input_name": name,
                    "kind": "field_jet_handle",
                    "field": field,
                    "requested_ordinary_Z_order": order,
                    "binding": copy.deepcopy(binding),
                }
            )
        elif name == "delta0":
            binding = geometry.get("delta")
            if not isinstance(binding, Mapping):
                binding = _unresolved(
                    f"{F_GEOMETRY_POINTER}.delta",
                    "delta0 has no geometry.delta candidate",
                    kind="adapter_delta0_handle",
                )
            else:
                binding = copy.deepcopy(binding)
                binding["kind"] = "adapter_delta0_handle"
                binding["adapter_input_name"] = "delta0"
                binding["source_geometry_role"] = f"{F_GEOMETRY_POINTER}.delta"
                binding["ordinary_Z_derivative_semantics_proved"] = False
                binding["source_value_selected"] = False
            result.append(
                {
                    "input_name": name,
                    "kind": "adapter_delta0_handle",
                    "source_geometry_name": "delta",
                    "binding": binding,
                }
            )
        else:
            unresolved_roles = {
                "delta1": f"{F_GEOMETRY_POINTER}.delta_Z",
                "delta2": f"{F_GEOMETRY_POINTER}.delta_ZZ",
                "Z": f"{F_CONTAINER_POINTER}.fixed_xi.Z",
            }
            binding = _unresolved(
                unresolved_roles[name],
                "the authenticated producer packet has no corresponding direct source handle",
                kind="adapter_input_handle",
            )
            binding["adapter_input_name"] = name
            result.append(
                {
                    "input_name": name,
                    "kind": "adapter_input_handle",
                    "source_owner_proved": False,
                    "binding": copy.deepcopy(binding),
                }
            )
    return result


def _empty_side() -> tuple[
    dict[str, Any], dict[str, Any], dict[str, Any], list[str]
]:
    fields = {
        field: _bind_field_rows(None, field)
        for field in FIELD_ORDERS
    }
    geometry = _bind_geometry(None)
    coordinates = _bind_coordinate_metadata(None)
    unresolved = [
        f"F.leading_fields_gap.{field}.Z{order}"
        for field, orders in FIELD_ORDERS.items()
        for order in orders
    ]
    unresolved.extend(f"F.geometry.{name}" for name in GEOMETRY_NAMES)
    unresolved.extend(
        f"F.coordinate_owner.{name}" for name in COORDINATE_NAMES
    )
    return fields, geometry, coordinates, unresolved


def _build_capsule_from_inner(
    inner: Mapping[str, Any] | None,
    receipt: Mapping[str, Any],
    *,
    wrapper_status: Any = None,
    parse_error: str | None = None,
) -> dict[str, Any]:
    f_side = None
    if isinstance(inner, Mapping):
        sides = inner.get("sides")
        if isinstance(sides, Mapping) and isinstance(sides.get("F"), Mapping):
            f_side = sides["F"]

    if f_side is None:
        fields, geometry, coordinates, unresolved = _empty_side()
    else:
        leading = f_side.get("leading_fields_gap")
        fields = {field: _bind_field_rows(leading, field) for field in FIELD_ORDERS}
        geometry = _bind_geometry(f_side)
        coordinates = _bind_coordinate_metadata(f_side)
        unresolved = []
        for field, orders in FIELD_ORDERS.items():
            rows = fields[field].get("rows", {})
            for order in orders:
                if rows.get(str(order), {}).get("status") != "bound_opaque_source_handle":
                    unresolved.append(f"F.leading_fields_gap.{field}.Z{order}")
        for name, binding in geometry.items():
            if binding.get("status") != "bound_opaque_source_handle":
                unresolved.append(f"F.geometry.{name}")
        for name in COORDINATE_NAMES:
            if name == "radius":
                unresolved.append("F.field_radius_owner")
            elif coordinates[name].get("status") != "bound_opaque_source_handle":
                unresolved.append(f"F.coordinate_owner.{name}")

    adapter_inputs = _adapter_inputs(fields, geometry)
    unresolved.extend(
        [
            "F.field_radius_owner",
            "F.m.Z3_semantics",
            "F.delta1",
            "F.delta2",
            "F.Z",
            "F.source_capture",
            "F.physical_scale_and_units",
            "F.ordinary_Z_derivative_semantics",
            "F.namespace_identity",
        ]
    )
    unresolved = sorted(set(unresolved))
    all_bound = not any(
        item.get("binding", {}).get("status") != "bound_opaque_source_handle"
        for item in adapter_inputs
    )
    status = "source_handles_bound_with_unresolved_gates" if all_bound else "partially_unresolved_source_handles"
    if parse_error is not None:
        status = "unresolved_result_packet"

    source_wrapper: dict[str, Any] = {
        "status": wrapper_status,
        "inner_schema": inner.get("schema") if isinstance(inner, Mapping) else None,
        "inner_status": inner.get("status") if isinstance(inner, Mapping) else None,
        "result_path": RESULT_PATH.as_posix(),
        "source_values_selected": False,
        "contract": {
            "outer_component_expected": PRIOR_COMPONENT,
            "outer_status_expected": "returned",
            "outer_schema_required_absent": True,
            "outer_error_required_absent": True,
            "inner_schema_expected": PRIOR_COMPONENT,
            "inner_status_expected": "incomplete_source_boundary_not_admissible_for_numeric_use",
            "inner_source_only_required_true": True,
            "inner_admissible_for_numeric_use_required_false": True,
            "validated": parse_error is None,
        },
    }
    if parse_error is not None:
        source_wrapper["parse_error"] = parse_error

    return {
        "schema": SCHEMA,
        "status": status,
        "source_only": True,
        "admissible_for_numeric_use": False,
        "source_identity_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_graph_traversed": False,
        "graph_edges_followed": False,
        "node_ids_are_values": False,
        "node_id_is_not_a_numeric_value": True,
        "ordinary_Z_semantics_proved": False,
        "ordinary_Z3_semantics_proved": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
        "field_radius_owner_proved": False,
        "namespace_identity_proved": False,
        "namespace_candidates_verified": False,
        "owner_candidates_verified": False,
        "raw_hash_authentication": copy.deepcopy(dict(receipt)),
        "source_wrapper": source_wrapper,
        "F": {
            "owner": F_OWNER,
            "namespace": F_SOURCE_NAMESPACE,
            "namespace_is_candidate_only": True,
            "source_container_role": F_CONTAINER_POINTER,
            "leading_fields_gap": fields,
            "geometry": geometry,
            "gap_coordinate_metadata": coordinates,
            "field_radius_owner": {
                "status": "unresolved",
                "reason": "radius metadata is not a field-radius owner",
                "source_value_selected": False,
                "field_radius_owner_proved": False,
            },
        },
        "adapter_input_order": list(INPUT_ORDER),
        "adapter_inputs": adapter_inputs,
        "adapter_input_count": len(adapter_inputs),
        "unresolved_gates": unresolved,
        "flags": {
            "source_value_selected": False,
            "numeric_values_computed": False,
            "source_graph_traversed": False,
            "source_capture_proved": False,
            "physical_scale_proved": False,
            "ordinary_Z_derivative_semantics_proved": False,
            "ordinary_Z3_semantics_proved": False,
            "m3_extra_row_promoted_to_ordinary_Z3": False,
            "m3_ordinary_Z3_semantics_proved": False,
            "delta1_binding_proved": False,
            "delta2_binding_proved": False,
            "Z_binding_proved": False,
            "radius_metadata_promoted_to_field_owner": False,
            "namespace_identity_proved": False,
            "namespace_candidates_verified": False,
            "owner_candidates_verified": False,
            "numeric_use_admissible": False,
        },
        "scope": {
            "consumed_result_path_only": True,
            "archive_read": False,
            "graph_rows_indexed": False,
            "graph_edges_followed": False,
            "old_leaf_imported": False,
            "old_leaf_called": False,
            "source_constructor_called": False,
            "source_checker_called": False,
        },
    }


def build_f_gap_source_jet_bindings(
    root: str | Path,
    source: Mapping[str, Any],
    hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return an opaque F seven-jet source-handle capsule.

    ``source`` and ``hashes`` are caller-owned raw-input manifests.  All paths
    named by either manifest are hashed before the one required result packet is
    parsed.  Missing fields in an authenticated packet are represented as
    ``unresolved`` handles; no numeric zero is substituted.
    """

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    root_path = Path(root).resolve()
    result_raw, receipt = _authenticate_raw_inputs(root_path, source, hashes)
    try:
        wrapper = _json_object(result_raw, RESULT_PATH.as_posix())
    except FGapSourceJetBindingInputError as exc:
        return _build_capsule_from_inner(
            None,
            receipt,
            wrapper_status="unresolved",
            parse_error=str(exc),
        )

    inner, wrapper_errors = _validated_prior_inner(wrapper)
    if wrapper_errors:
        return _build_capsule_from_inner(
            None,
            receipt,
            wrapper_status=wrapper.get("status"),
            parse_error="; ".join(wrapper_errors),
        )
    return _build_capsule_from_inner(
        inner,
        receipt,
        wrapper_status=wrapper.get("status"),
    )


def bind_f_gap_source_jets(
    root: str | Path,
    prior_result_path: str | Path,
    prior_result_sha256: str,
    source_manifest: Mapping[str, Any],
    input_hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Runner-compatible alias with an explicit prior-result pin.

    The coordinating runner supplies the result path and SHA separately.  They
    must agree with this leaf's fixed publication pin before the common builder
    is entered; the result is still authenticated together with every manifest
    input before JSON parsing.
    """

    supplied_path = str(prior_result_path).replace("\\", "/")
    if supplied_path.endswith("/" + RESULT_PATH.as_posix()):
        supplied_path = RESULT_PATH.as_posix()
    _require(
        supplied_path == RESULT_PATH.as_posix(),
        "runner prior-result path differs from the fixed publication",
    )
    _require(
        isinstance(prior_result_sha256, str)
        and prior_result_sha256.lower() == RESULT_SHA256,
        "runner prior-result SHA differs from the fixed publication",
    )
    return build_f_gap_source_jet_bindings(root, source_manifest, input_hashes)


def _synthetic_ref(role: str, node: int, *, namespace: str | None = None) -> dict[str, Any]:
    value: dict[str, Any] = {
        "status": "direct_source_ref",
        "source_role": role,
        "source_node_id": node,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }
    if namespace is not None:
        value["namespace"] = namespace
    return value


def _synthetic_rows(field: str, count: int, *, namespace: str | None = None) -> dict[str, Any]:
    role = f"{F_LEADING_POINTER}.{field}"
    rows = [_synthetic_ref(f"{role}[{index}]", 1000 + count * 10 + index, namespace=namespace) for index in range(count)]
    return {
        "status": "direct_source_rows",
        "source_role": role,
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z0_to_Z2": rows[:3],
        "all_saved_rows": rows,
        "source_values_selected": False,
    }


def _synthetic_side(*, missing_m3: bool = False, wrong_namespace: bool = False) -> dict[str, Any]:
    leading = {
        "U": _synthetic_rows("U", 3),
        "m": _synthetic_rows("m", 3 if missing_m3 else 4),
        "U_y": _synthetic_rows("U_y", 2),
        "m_y": _synthetic_rows("m_y", 3),
    }
    if wrong_namespace:
        leading["U"]["ordinary_Z0_to_Z2"][1]["namespace"] = "wrong_namespace"
        leading["U"]["all_saved_rows"][1]["namespace"] = "wrong_namespace"
    geometry = {
        name: _synthetic_ref(f"{F_GEOMETRY_POINTER}.{name}", 2000 + index)
        for index, name in enumerate(GEOMETRY_NAMES)
    }
    coordinates = {
        "query": _synthetic_ref(f"{F_COORDINATE_POINTER}.query", 3000),
        "y": _synthetic_ref(f"{F_COORDINATE_POINTER}.y", 3001),
        "radius": {
            "status": "unresolved",
            "reason": "current geometry metadata",
            "current_radius_metadata": _synthetic_ref(f"{F_COORDINATE_POINTER}.radius", 3002),
            "source_value_selected": False,
        },
        "coordinate_literal": "xi=mu*log(R/Rp)",
    }
    return {
        "leading_fields_gap": leading,
        "geometry": geometry,
        "coordinate_owner": coordinates,
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Exercise wrapper, field-record, m3, namespace, and radius guards."""

    complete = _build_capsule_from_inner(
        {"sides": {"F": _synthetic_side()}},
        {"all_manifest_raw_hashes_checked_before_parse": False},
        wrapper_status="synthetic",
    )
    missing_m3 = _build_capsule_from_inner(
        {"sides": {"F": _synthetic_side(missing_m3=True)}},
        {"all_manifest_raw_hashes_checked_before_parse": False},
        wrapper_status="synthetic",
    )
    wrong_namespace = _build_capsule_from_inner(
        {"sides": {"F": _synthetic_side(wrong_namespace=True)}},
        {"all_manifest_raw_hashes_checked_before_parse": False},
        wrapper_status="synthetic",
    )
    bad_wrapper = {
        "schema": "unexpected_outer_schema",
        "component": PRIOR_COMPONENT,
        "status": "returned",
        "result": {
            "schema": PRIOR_COMPONENT,
            "status": "incomplete_source_boundary_not_admissible_for_numeric_use",
            "source_only": True,
            "admissible_for_numeric_use": False,
            "sides": {"F": _synthetic_side()},
        },
    }
    bad_wrapper_inner, bad_wrapper_errors = _validated_prior_inner(bad_wrapper)
    bad_wrapper_capsule = _build_capsule_from_inner(
        bad_wrapper_inner,
        {"all_manifest_raw_hashes_checked_before_parse": False},
        wrapper_status=bad_wrapper.get("status"),
        parse_error="; ".join(bad_wrapper_errors),
    )
    bad_field_side = _synthetic_side()
    bad_field_side["leading_fields_gap"]["m"]["status"] = "unresolved"
    bad_field_inner = {"sides": {"F": bad_field_side}}
    bad_field_record = _build_capsule_from_inner(
        bad_field_inner,
        {"all_manifest_raw_hashes_checked_before_parse": False},
        wrapper_status="synthetic",
    )

    complete_inputs = {item["input_name"]: item for item in complete["adapter_inputs"]}
    missing_inputs = {item["input_name"]: item for item in missing_m3["adapter_inputs"]}
    wrong_inputs = {item["input_name"]: item for item in wrong_namespace["adapter_inputs"]}
    bad_wrapper_inputs = {item["input_name"]: item for item in bad_wrapper_capsule["adapter_inputs"]}
    bad_field_inputs = {item["input_name"]: item for item in bad_field_record["adapter_inputs"]}
    checks = {
        "sixteen_adapter_inputs": complete["adapter_input_count"] == 16
        and complete["adapter_input_order"] == list(INPUT_ORDER),
        "complete_m3_uses_fourth_all_saved_row": (
            complete_inputs["m3"]["binding"].get("selected_from") == "all_saved_rows[3]"
            and complete_inputs["m3"]["binding"].get("extra_row_semantics_promoted") is False
            and complete_inputs["m3"]["binding"].get("ordinary_Z3_semantics_proved") is False
        ),
        "delta0_is_only_geometry_delta_candidate": (
            complete_inputs["delta0"]["binding"].get("status") == "bound_opaque_source_handle"
            and complete_inputs["delta0"]["binding"].get("source_geometry_role")
            == f"{F_GEOMETRY_POINTER}.delta"
        ),
        "delta1_delta2_and_Z_stay_unresolved": all(
            complete_inputs[name]["binding"].get("status") == "unresolved"
            for name in ("delta1", "delta2", "Z")
        ),
        "bad_wrapper_is_unresolved": (
            bool(bad_wrapper_errors)
            and bad_wrapper_capsule["status"] == "unresolved_result_packet"
            and all(item["binding"].get("status") == "unresolved" for item in bad_wrapper_inputs.values())
        ),
        "bad_field_record_is_unresolved": (
            bad_field_record["F"]["leading_fields_gap"]["m"].get("status")
            == "unresolved_field_record_contract"
            and all(
                bad_field_inputs[name]["binding"].get("status") == "unresolved"
                for name in ("m0", "m1", "m2", "m3")
            )
        ),
        "missing_m3_is_unresolved_not_zero": (
            missing_inputs["m3"]["binding"].get("status") == "unresolved"
            and "source_node_id" not in missing_inputs["m3"]["binding"]
        ),
        "namespace_difference_is_rejected": (
            wrong_inputs["U1"]["binding"].get("status") == "unresolved"
            and wrong_inputs["U1"]["binding"].get("observed_namespace") == "wrong_namespace"
        ),
        "radius_metadata_is_not_field_owner": (
            complete["F"]["gap_coordinate_metadata"]["radius"].get("metadata_radius_not_field_owner") is True
            and complete["F"]["gap_coordinate_metadata"]["radius"].get("field_radius_owner_proved") is False
            and complete["F"]["field_radius_owner"].get("field_radius_owner_proved") is False
        ),
        "all_handles_remain_opaque": all(
            item["binding"].get("node_id_is_not_a_numeric_value") is True
            for item in complete["adapter_inputs"]
            if item["binding"].get("status") == "bound_opaque_source_handle"
        ),
        "semantic_and_capture_flags_stay_false": all(
            complete["flags"].get(name) is False
            for name in (
                "source_capture_proved",
                "physical_scale_proved",
                "ordinary_Z_derivative_semantics_proved",
                "ordinary_Z3_semantics_proved",
                "m3_ordinary_Z3_semantics_proved",
                "delta1_binding_proved",
                "delta2_binding_proved",
                "Z_binding_proved",
                "namespace_identity_proved",
                "namespace_candidates_verified",
                "owner_candidates_verified",
            )
        ),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": SCHEMA,
        "status": "synthetic_checks_passed" if all(checks.values()) else "synthetic_checks_failed",
        "passed": all(checks.values()),
        "checks_passed": all(checks.values()),
        "synthetic_only": True,
        "checks": checks,
        "check_count": len(checks),
        "source_values_selected": False,
        "source_values_computed": False,
        "source_graph_traversed": False,
        "source_identity_proved": False,
        "archive_read": False,
        "real_inputs_read": False,
    }


__all__ = [
    "SCHEMA",
    "RESULT_PATH",
    "RESULT_SHA256",
    "FGapSourceJetBindingInputError",
    "build_f_gap_source_jet_bindings",
    "bind_f_gap_source_jets",
    "verify_synthetic_fixture",
]
