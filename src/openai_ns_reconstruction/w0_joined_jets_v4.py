"""Standalone source-bound join of the W0 ordinary ``Z0/Z1/Z2`` sectors.

The caller supplies two already authenticated packets: the old W0 packet
owns ordinary ``Z0`` and ``Z1`` and the new W0 packet owns ordinary ``Z2``.
Only ``F1``, ``Uz1`` and ``K1`` are admitted from the new packet.  The old
packet's ``d_x_F1`` port is retained at orders zero and one for the product
consumer; its Z2 port, ``P1``, ``d_x_Uz1`` and all feedback remain unknown.

This is a pure record join.  It does not import a producer, call an ancestor,
evaluate a source, add intervals, materialise a formal atom, add ``g5`` or
append an m15 term.  ``state`` contains the canonical term lists separately
under ``core`` and ``local`` so a downstream product may consume them without
this leaf silently summing the two sectors.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from fractions import Fraction
from typing import Any


SCHEMA = "w0_joined_jets_v4"
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
CELL_COUNT = 4
ORDERS = (0, 1, 2)
OLD_ORDERS = (0, 1)
Z2_ORDER = 2
FIELDS = ("F1", "Uz1", "K1")
DERIVATIVE_PORT = "d_x_F1"
OLD_FIELDS = FIELDS + (DERIVATIVE_PORT,)
UNKNOWN_FIELDS = ("P1", "d_x_Uz1")
ALL_COMPONENTS = OLD_FIELDS + UNKNOWN_FIELDS
ROWS = {"F1": "5", "Uz1": "6", "K1": "6", "d_x_F1": "5"}
COMPONENT_IDS = {
    "F1": 1,
    "Uz1": 2,
    "K1": 3,
    "P1": 4,
    "d_x_F1": 5,
    "d_x_Uz1": 6,
}


def _fail(label: str, message: str) -> None:
    raise ValueError(f"{label}: {message}")


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(label, "mapping required")
    return value


def _lookup(value: Any, key: Any, label: str) -> Any:
    if isinstance(value, Mapping):
        if key in value:
            return value[key]
        text = str(key)
        if text in value:
            return value[text]
    _fail(label, f"missing key {key!r}")


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _first(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _raw_interval(value: Any, label: str) -> tuple[Any, Any]:
    raw = getattr(value, "_mpi_", None)
    if raw is None and isinstance(value, Mapping):
        lower = _first(value, "lower_exact_mpf_tuple", "lower_exact_tuple")
        upper = _first(value, "upper_exact_mpf_tuple", "upper_exact_tuple")
        if lower is not None and upper is not None:
            raw = (tuple(lower), tuple(upper))
    if not isinstance(raw, (tuple, list)) or len(raw) != 2:
        _fail(label, "directed interval required")
    return raw[0], raw[1]


def _endpoint_fraction(endpoint: Any, label: str) -> Fraction:
    if not isinstance(endpoint, (tuple, list)) or len(endpoint) != 4:
        _fail(label, "exact MPF endpoint tuple required")
    sign, mantissa, exponent, bit_count = endpoint
    if any(type(item) is not int for item in (sign, mantissa, exponent, bit_count)):
        _fail(label, "MPF endpoint tuple must contain integers")
    if sign not in (0, 1) or mantissa < 0:
        _fail(label, "invalid MPF endpoint tuple")
    if mantissa == 0:
        return Fraction(0)
    numerator = -mantissa if sign else mantissa
    return Fraction(numerator * (2**exponent), 1) if exponent >= 0 else Fraction(numerator, 2 ** (-exponent))


def _point(value: Any, label: str) -> Fraction:
    lower, upper = _raw_interval(value, label)
    if lower != upper:
        _fail(label, "q endpoint must be a singleton")
    return _endpoint_fraction(lower, label)


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        try:
            return _raw_interval(left, "left") == _raw_interval(right, "right")
        except ValueError:
            return left == right


def _q_cells(cells: Any, label: str) -> list[dict[str, Any]]:
    if not isinstance(cells, (list, tuple)) or len(cells) != CELL_COUNT:
        _fail(label, "four q cells required")
    result: list[dict[str, Any]] = []
    for index, raw in enumerate(cells):
        cell = _mapping(raw, f"{label}[{index}]")
        if cell.get("index", index) != index:
            _fail(f"{label}[{index}]", "cell index mismatch")
        for key in ("q_lower", "q_upper", "s_lower", "s_upper"):
            if key not in cell:
                _fail(f"{label}[{index}]", f"missing {key}")
        lower = _point(cell["q_lower"], f"{label}[{index}].q_lower")
        upper = _point(cell["q_upper"], f"{label}[{index}].q_upper")
        if lower != Fraction(index, 4) or upper != Fraction(index + 1, 4):
            _fail(f"{label}[{index}]", "q endpoints do not identify the declared cell")
        if lower >= upper:
            _fail(f"{label}[{index}]", "q interval is not ordered")
        _raw_interval(cell["s_lower"], f"{label}[{index}].s_lower")
        _raw_interval(cell["s_upper"], f"{label}[{index}].s_upper")
        result.append(deepcopy(dict(cell)))
    return result


def _term_ids(term: Mapping[str, Any], label: str) -> list[Any]:
    origin = term.get("origin")
    if not isinstance(origin, Mapping):
        _fail(label, "origin metadata required")
    ids: list[Any] = []
    for value in (term.get("source_ids"), origin.get("source_ids"), term.get("input_id"), origin.get("input_id")):
        for item in _as_list(value):
            if item is not None and item not in ids:
                ids.append(deepcopy(item))
    if not ids:
        _fail(label, "source ID required")
    return ids


def _validate_term(
    term: Any,
    label: str,
    *,
    role: str,
    order: int,
    field: str,
    cell: int | None,
) -> dict[str, Any]:
    if not isinstance(term, Mapping):
        _fail(label, "term mapping required")
    powers = term.get("powers")
    if not isinstance(powers, (list, tuple)) or len(powers) != len(ATOMS):
        _fail(label, "four formal powers required")
    if any(type(power) is not int for power in powers):
        _fail(label, "formal powers must be integers")
    _raw_interval(term.get("coefficient"), label + ".coefficient")
    origin = term.get("origin")
    if not isinstance(origin, Mapping):
        _fail(label, "origin metadata required")
    _term_ids(term, label)
    if "error_classes" not in term and "error_classes" not in origin:
        _fail(label, "error metadata required")

    direct = origin.get("operation") in {
        "local_green_cell",
        "w0_z2_four_cell_v1_local",
        "w0_z2_four_cell_v2_local",
        "w0_z2_four_cell_v2_join_local",
    }
    for key in ("packet_role", "source_packet_role"):
        marker = term.get(key)
        if marker is not None and marker != role:
            _fail(label, f"{key} claims {marker!r}, expected {role!r}")
    marker_order = term.get("ordinary_Z_order")
    if marker_order is not None and marker_order != order:
        _fail(label, "ordinary Z order metadata disagrees with packet order")
    marker_field = term.get("field")
    if marker_field is not None and marker_field != field:
        _fail(label, "field metadata disagrees with packet field")
    marker_row = term.get("source_row")
    if marker_row is None and direct:
        marker_row = origin.get("source_row")
    if marker_row is not None and str(marker_row) != ROWS[field]:
        _fail(label, "source row metadata disagrees with packet field")

    ordinary_derivative = term.get("source_z_derivative_is_ordinary")
    factorial_divided = term.get("ordinary_Z_factorial_divided")
    if ordinary_derivative is not None:
        if order != Z2_ORDER or ordinary_derivative is not True:
            _fail(label, "Z2 source must be marked as an ordinary derivative")
    if factorial_divided is not None:
        if order != Z2_ORDER or factorial_divided is not False:
            _fail(label, "Z2 source must not be factorial-divided")
    if cell is not None:
        marker_cell = term.get("phase_cell")
        if marker_cell is None and direct:
            marker_cell = origin.get("phase_cell")
        if type(marker_cell) is not int or marker_cell != cell:
            _fail(label, "phase cell metadata disagrees with local container")
        refs = term.get("source_node_refs", origin.get("source_node_refs"))
        for ref in _as_list(refs):
            if not isinstance(ref, Mapping):
                continue
            ref_cell = _first(ref, "phase_cell", "cell")
            if ref_cell is not None and ref_cell != cell:
                _fail(label, "source node has the wrong q-cell")
            ref_role = _first(ref, "packet_role", "source_packet_role")
            if ref_role is not None and ref_role != role:
                _fail(label, "source node has the wrong packet role")
    return deepcopy(dict(term))


def _terms(
    value: Any,
    label: str,
    *,
    role: str,
    order: int,
    field: str,
    cell: int | None,
) -> list[dict[str, Any]]:
    if not isinstance(value, (list, tuple)) or not value:
        _fail(label, "non-empty source term list required")
    return [
        _validate_term(term, f"{label}[{index}]", role=role, order=order, field=field, cell=cell)
        for index, term in enumerate(value)
    ]


def _role_ok(value: Any, role: str) -> bool:
    aliases = {
        "z0_z1": {"z0_z1", "z01", "old", "old_w0", "old_W0_Z0_Z1"},
        "z2": {"z2", "new", "new_w0", "w0_z2"},
    }
    return value in aliases[role]


def _unwrap(packet: Any, label: str, role: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    wrapper = _mapping(packet, label)
    verified = _first(
        wrapper,
        "source_verified",
        "verified_source_package",
        "authenticated_source_package",
        "source_package_verified",
    )
    payload = _first(wrapper, "payload", "packet", "source_payload", "result")
    if payload is None:
        payload = wrapper
    payload = _mapping(payload, label + ".payload")
    payload_verified = _first(
        payload,
        "source_verified",
        "verified_source_package",
        "authenticated_source_package",
        "source_package_verified",
    )
    if verified is None:
        verified = payload_verified
    if verified is not True:
        _fail(label, "caller-authenticated source package required")
    if payload_verified is not None and payload_verified is not True:
        _fail(label, "payload authentication marker is false")
    packet_role = _first(wrapper, "packet_role", "jet_role", "source_packet_role")
    payload_role = _first(payload, "packet_role", "jet_role", "source_packet_role")
    if packet_role is None:
        packet_role = payload_role
    if not _role_ok(packet_role, role):
        _fail(label, f"packet role {packet_role!r} is not {role!r}")
    if payload_role is not None and not _role_ok(payload_role, role):
        _fail(label, f"payload packet role {payload_role!r} is not {role!r}")
    family = _first(wrapper, "source_family", "source_family_id")
    payload_family = _first(payload, "source_family", "source_family_id")
    if family is None:
        family = payload_family
    if family is not None and payload_family is not None and family != payload_family:
        _fail(label, "wrapper and payload source families differ")
    atoms = _first(wrapper, "source_atoms", "formal_atoms")
    payload_atoms = _first(payload, "source_atoms", "formal_atoms")
    if atoms is None:
        atoms = payload_atoms
    if atoms is not None and payload_atoms is not None and atoms != payload_atoms:
        _fail(label, "wrapper and payload source atoms differ")
    basis = _first(wrapper, "source_power_basis", "formal_power_basis")
    payload_basis = _first(payload, "source_power_basis", "formal_power_basis")
    if basis is None:
        basis = payload_basis
    if basis is not None and payload_basis is not None and tuple(basis) != tuple(payload_basis):
        _fail(label, "wrapper and payload source power bases differ")
    return payload, {
        "source_family": deepcopy(family),
        "source_atoms": deepcopy(atoms),
        "source_power_basis": deepcopy(basis),
    }


def _packet(packet: Any, label: str, role: str, orders: tuple[int, ...]) -> dict[str, Any]:
    payload, metadata = _unwrap(packet, label, role)
    if tuple(metadata.get("source_power_basis") or ()) != ATOMS:
        _fail(label, "source power basis differs from ATOMS")
    if payload.get("source_power_basis") is not None and tuple(payload["source_power_basis"]) != ATOMS:
        _fail(label, "payload source power basis differs from ATOMS")
    if role == "z2" and (
        payload.get("source_z2_is_ordinary_derivative") is not True
        or payload.get("ordinary_Z_factorial_divided") is not False
    ):
        _fail(label, "Z2 must be an ordinary derivative without factorial division")
    if tuple(payload.get("ordinary_Z_orders") or ()) != orders:
        _fail(label, f"ordinary orders must be {orders}")
    atoms = metadata.get("source_atoms")
    if atoms is not None and (not isinstance(atoms, Mapping) or set(atoms) != set(ATOMS)):
        _fail(label, "source atom names are incomplete")
    q_cells = _q_cells(payload.get("q_cells"), label + ".q_cells")

    core_raw = payload.get("core_transport")
    if not isinstance(core_raw, Mapping):
        _fail(label, "core_transport missing")
    required = OLD_FIELDS if role == "z0_z1" else FIELDS
    if role == "z2" and all(field in core_raw for field in FIELDS):
        by_order: Mapping[Any, Any] = {Z2_ORDER: core_raw}
    else:
        by_order = core_raw
    core: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in orders:
        row = _mapping(_lookup(by_order, order, label + ".core_transport"), f"{label}.core.Z{order}")
        core[order] = {
            field: _terms(
                _lookup(row, field, f"{label}.core.Z{order}"),
                f"{label}.core.Z{order}.{field}",
                role=role,
                order=order,
                field=field,
                cell=None,
            )
            for field in required
        }

    local_raw = payload.get("local")
    if not isinstance(local_raw, Mapping):
        _fail(label, "local rows missing")
    local: dict[int, dict[int, dict[str, list[dict[str, Any]]]]] = {order: {} for order in orders}

    def local_cell(raw_cell: Any, order: int, cell: int) -> dict[str, list[dict[str, Any]]]:
        row = _mapping(raw_cell, f"{label}.local.Z{order}.cell{cell}")
        return {
            field: _terms(
                _lookup(row, field, f"{label}.local.Z{order}.cell{cell}"),
                f"{label}.local.Z{order}.cell{cell}.{field}",
                role=role,
                order=order,
                field=field,
                cell=cell,
            )
            for field in required
        }

    if role == "z0_z1":
        for order in orders:
            order_raw = _mapping(_lookup(local_raw, order, label + ".local"), f"{label}.local.Z{order}")
            if set(str(key) for key in order_raw) != {"0", "1", "2", "3"}:
                _fail(label + ".local", "four cells required")
            for cell in range(CELL_COUNT):
                local[order][cell] = local_cell(_lookup(order_raw, cell, label + ".local"), order, cell)
    else:
        first = _mapping(_lookup(local_raw, 0, label + ".local"), label + ".local[0]")
        if any(field in first for field in required):
            if set(str(key) for key in local_raw) != {"0", "1", "2", "3"}:
                _fail(label + ".local", "four cells required")
            for cell in range(CELL_COUNT):
                local[Z2_ORDER][cell] = local_cell(_lookup(local_raw, cell, label + ".local"), Z2_ORDER, cell)
        else:
            order_raw = _mapping(_lookup(local_raw, Z2_ORDER, label + ".local"), label + ".local.Z2")
            if set(str(key) for key in order_raw) != {"0", "1", "2", "3"}:
                _fail(label + ".local.Z2", "four cells required")
            for cell in range(CELL_COUNT):
                local[Z2_ORDER][cell] = local_cell(_lookup(order_raw, cell, label + ".local.Z2"), Z2_ORDER, cell)
    return {
        "payload": payload,
        "metadata": metadata,
        "source_family": metadata.get("source_family"),
        "source_atoms": atoms,
        "q_cells": q_cells,
        "core": core,
        "local": local,
    }


def _same_q_cells(left: list[dict[str, Any]], right: list[dict[str, Any]]) -> None:
    for index, (old_cell, new_cell) in enumerate(zip(left, right)):
        for key in ("q_lower", "q_upper", "s_lower", "s_upper"):
            if not _same_interval(old_cell[key], new_cell[key]):
                _fail(f"q_cells[{index}]", f"old/new {key} differ")


def _ref(container: str, order: int, field: str, index: int, cell: int | None = None) -> dict[str, Any]:
    result = {
        "kind": "source_term_ref",
        "container": container,
        "ordinary_Z_order": order,
        "field": field,
        "term_index": index,
    }
    if cell is not None:
        result["phase_cell"] = cell
    return result


def _unknown(component: str, order: int, reason: str) -> dict[str, Any]:
    return {
        "status": "unknown",
        "component": component,
        "ordinary_Z_order": order,
        "core": [],
        "local": [],
        "reason": reason,
        "zero_padded": False,
    }


def build(
    old_packet: Any = None,
    z2_packet: Any = None,
    *,
    z01_packet: Any = None,
    new_packet: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
) -> dict[str, Any]:
    """Return a source-only, product-ready W0 three-order jet join."""

    if old_packet is None:
        old_packet = z01_packet
    elif z01_packet is not None and z01_packet is not old_packet:
        _fail("old_packet", "conflicting z01_packet alias")
    if z2_packet is None:
        z2_packet = new_packet
    elif new_packet is not None and new_packet is not z2_packet:
        _fail("z2_packet", "conflicting new_packet alias")
    if old_packet is None or z2_packet is None:
        _fail("build", "old Z0/Z1 and new Z2 packets are both required")

    old = _packet(old_packet, "old_packet", "z0_z1", OLD_ORDERS)
    new = _packet(z2_packet, "z2_packet", "z2", (Z2_ORDER,))
    family = old["source_family"]
    if family is None or new["source_family"] is None:
        _fail("source_family", "both packages must provide source_family")
    if family != new["source_family"]:
        _fail("source_family", "old and Z2 packages differ")
    if source_family is not None and source_family != family:
        _fail("source_family", "caller family differs from packages")
    old_atoms = old.get("source_atoms")
    new_atoms = new.get("source_atoms")
    if old_atoms is not None and new_atoms is not None and old_atoms != new_atoms:
        _fail("source_atoms", "old and Z2 atom maps differ")
    if source_atoms is not None:
        if not isinstance(source_atoms, Mapping) or set(source_atoms) != set(ATOMS):
            _fail("source_atoms", "caller source atoms must contain ATOMS")
        if old_atoms is not None and old_atoms != source_atoms:
            _fail("source_atoms", "caller atoms differ from old package")
        if new_atoms is not None and new_atoms != source_atoms:
            _fail("source_atoms", "caller atoms differ from Z2 package")
    _same_q_cells(old["q_cells"], new["q_cells"])

    core_transport = {
        0: {field: deepcopy(old["core"][0][field]) for field in OLD_FIELDS},
        1: {field: deepcopy(old["core"][1][field]) for field in OLD_FIELDS},
        2: {field: deepcopy(new["core"][2][field]) for field in FIELDS},
    }
    local = {
        order: {
            cell: {
                field: deepcopy((old if order < 2 else new)["local"][order][cell][field])
                for field in (OLD_FIELDS if order < 2 else FIELDS)
            }
            for cell in range(CELL_COUNT)
        }
        for order in ORDERS
    }

    state: dict[int, dict[tuple[int, int], dict[str, Any]]] = {}
    unknown_state: dict[int, dict[tuple[int, int], dict[str, Any]]] = {}
    for cell in range(CELL_COUNT):
        state[cell] = {}
        unknown_state[cell] = {}
        for order in ORDERS:
            fields = OLD_FIELDS if order < 2 else FIELDS
            for field in fields:
                component = COMPONENT_IDS[field]
                core_terms = core_transport[order][field]
                local_terms = local[order][cell][field]
                state[cell][(component, order)] = {
                    "status": "known",
                    "component": field,
                    "component_id": component,
                    "ordinary_Z_order": order,
                    "core": core_terms,
                    "local": local_terms,
                    "core_refs": [_ref("core_transport", order, field, i) for i in range(len(core_terms))],
                    "local_refs": [_ref("local", order, field, i, cell) for i in range(len(local_terms))],
                    "implicit_addition": False,
                    "core_history_retained_once": True,
                }
            for field in UNKNOWN_FIELDS:
                unknown_state[cell][(COMPONENT_IDS[field], order)] = _unknown(
                    field, order, "component is outside the bounded product state; no zero padding"
                )
        unknown_state[cell][(COMPONENT_IDS[DERIVATIVE_PORT], Z2_ORDER)] = _unknown(
            DERIVATIVE_PORT, Z2_ORDER, "new Z2 packet has no d_x_F1 row"
        )

    feedback = {
        "status": "unknown",
        "included": False,
        "eliminated": False,
        "cancelled": False,
        "zero_padded": False,
        "source_ids": [],
        "reason": "higher feedback remainder was not supplied to this join",
    }
    return {
        "schema": SCHEMA,
        "source_bound": True,
        "source_only": True,
        "source_family": deepcopy(family),
        "source_atoms": deepcopy(source_atoms if source_atoms is not None else (old_atoms if old_atoms is not None else new_atoms)),
        "source_power_basis": list(ATOMS),
        "source_packages": {
            "z0_z1": {
                "verified": True,
                "packet_role": "z0_z1",
                "source_family": deepcopy(old["source_family"]),
                "ordinary_Z_orders": list(OLD_ORDERS),
                "q_cells": {"container": "q_cells", "packet": "old"},
            },
            "z2": {
                "verified": True,
                "packet_role": "z2",
                "source_family": deepcopy(new["source_family"]),
                "ordinary_Z_orders": [Z2_ORDER],
                "q_cells": {"container": "q_cells", "packet": "new"},
            },
        },
        "ordinary_Z_orders": list(ORDERS),
        "ordinary_Z_order": Z2_ORDER,
        "fields": list(FIELDS),
        "component_ids": dict(COMPONENT_IDS),
        "rows": dict(ROWS),
        "q_cells": deepcopy(old["q_cells"]),
        "core_transport": core_transport,
        "core_history": {"container": "core_transport", "retained_once": True, "per_cell_repetition": False},
        "local": local,
        "state": state,
        "product_state": state,
        "unknown_state": unknown_state,
        "joined_jet": state,
        "state_shape": "state[cell][(component, ordinary_Z_order)] = {'core': [terms], 'local': [terms]}",
        "state_semantics": "state[cell][(component_id, ordinary_Z_order)] has separate core/local source-term records; no implicit sum",
        "known_components_by_order": {
            "0": [COMPONENT_IDS[field] for field in OLD_FIELDS],
            "1": [COMPONENT_IDS[field] for field in OLD_FIELDS],
            "2": [COMPONENT_IDS[field] for field in FIELDS],
        },
        "unknown_components": {
            field: {"component_id": COMPONENT_IDS[field], "orders": list(ORDERS), "status": "unknown", "zero_padded": False}
            for field in UNKNOWN_FIELDS
        },
        "unknown_Z2_d_x_F1": {"component_id": COMPONENT_IDS[DERIVATIVE_PORT], "ordinary_Z_order": 2, "status": "unknown", "zero_padded": False},
        "feedback_remainder": feedback,
        "core_history_reset": False,
        "core_retained_once": True,
        "local_added_once_in_definition": False,
        "local_terms_are_separate": True,
        "source_ids_preserved": True,
        "intervals_preserved": True,
        "formal_powers_preserved": True,
        "error_classes_preserved": True,
        "q_cell_ownership_checked": True,
        "source_package_ownership_checked": True,
        "source_terms_are_values": False,
        "source_values_computed": False,
        "actual_forced_evaluation": False,
        "formal_powers_materialized": False,
        "source_z2_is_ordinary_derivative": True,
        "ordinary_Z_factorial_divided": False,
        "geometry_Z_independent": True,
        "hb_assumption": "0<hb<1/2",
        "hb_materialized": False,
        "formal_hb2_materialized": False,
        "source_z2_whole_collar_enclosure": True,
        "source_endpoint_not_selected": True,
        "kernel_bound_includes_formal_hb2": False,
        "physical_or_cartesian_measure_included": False,
        "no_physical_cartesian_measure": True,
        "g_uniform_on_entire_first_collar_required": True,
        "old_W1_consumed": False,
        "uncertainty_preserved": True,
        "live_owner_claimed": False,
        "physical_cartesian_frame_claimed": False,
        "uncertainty_correlation_proved": False,
        "external_W0_Z0_Z1_recomputed": False,
        "g5_added": False,
        "g5_recomputed": False,
        "g5_reference_only": True,
        "m15_computed": False,
        "m15_ledger_changed": False,
        "m15_terms_added": 0,
        "m15_ledger_size_before": 212,
        "m15_ledger_size_after": 212,
        "ledger_appended": False,
        "allword_tail_added": False,
        "feedback_remainder_eliminated": False,
        "ancestor_imports": False,
        "ancestor_called": False,
        "receipt_reads": False,
        "quadrature_used": False,
        "midpoint_values_used": False,
        "physical_cartesian_owner_proved": False,
        "same_source_cancellation_claimed": False,
    }


__all__ = [
    "ALL_COMPONENTS",
    "ATOMS",
    "CELL_COUNT",
    "COMPONENT_IDS",
    "DERIVATIVE_PORT",
    "FIELDS",
    "OLD_FIELDS",
    "OLD_ORDERS",
    "ORDERS",
    "ROWS",
    "SCHEMA",
    "build",
]
