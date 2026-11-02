"""Conditional target-tail error transport to ordinary m15 swirl."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_kchi_p_tail_m15_v1"
CELL_COUNT = 8
ORDERS = ("Z0", "Z1")
FIELD = "F1"
BRANCHES = ("G0_P_Rg", "D_P_Rg", "R_P_G0g", "R_P_Dg", "R_P_Rg")
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_POWER_BASIS = (
    "hb", "F0_anchor", "Pstar_squared", "Lambda", "F0_local", "L_local"
)
Z1_ROUTES = ("Phi1_E0", "Phi0_E1", "Phi0_ell_E0")


class PTailM15InputError(ValueError):
    """Raised when a saved packet cannot be bound safely."""


def _copy(value: Any) -> Any:
    return deepcopy(value)


def _mapping(value: Any, label: str) -> Mapping[Any, Any]:
    if not isinstance(value, Mapping):
        raise PTailM15InputError(label + " must be a mapping")
    return value


def _sequence(value: Any, label: str, *, allow_empty: bool = True) -> list[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PTailM15InputError(label + " must be a sequence")
    result = list(value)
    if not allow_empty and not result:
        raise PTailM15InputError(label + " must not be empty")
    return result


def _lookup(mapping: Mapping[Any, Any], key: Any, label: str) -> Any:
    candidates: list[Any] = [key]
    if isinstance(key, str):
        if key.startswith("Z") and key[1:].isdigit():
            number = int(key[1:])
            candidates.extend((number, str(number)))
        elif key.isdigit():
            number = int(key)
            candidates.extend((number, "Z" + key))
    elif isinstance(key, int):
        candidates.extend((str(key), "Z" + str(key)))
    for candidate in candidates:
        if candidate in mapping:
            return mapping[candidate]
    raise PTailM15InputError(label + " is missing " + str(key))


def _payload(value: Any, label: str, *names: str) -> Mapping[Any, Any]:
    """Accept a direct value or the small receipt/result wrapper."""

    root = _mapping(value, label)
    for name in names:
        candidate = root.get(name)
        if isinstance(candidate, Mapping):
            return candidate
    result = root.get("result")
    if isinstance(result, Mapping):
        for name in names:
            candidate = result.get(name)
            if isinstance(candidate, Mapping):
                return candidate
        output = result.get("output")
        if isinstance(output, Mapping):
            for name in names:
                candidate = output.get(name)
                if isinstance(candidate, Mapping):
                    return candidate
            return output
        return result
    output = root.get("output")
    if isinstance(output, Mapping):
        for name in names:
            candidate = output.get(name)
            if isinstance(candidate, Mapping):
                return candidate
        return output
    return root


def _source_family_logs(packet: Mapping[Any, Any], label: str) -> tuple[Any, Any]:
    family = packet.get("source_family")
    logs = packet.get("source_logs")
    if family is None:
        raise PTailM15InputError(label + ".source_family is required")
    if logs is None:
        raise PTailM15InputError(label + ".source_logs is required")
    return _copy(family), _copy(logs)


def _bind_source_labels(
    tail: Mapping[Any, Any],
    geometry: Mapping[Any, Any],
    symbolic: Mapping[Any, Any],
    connection: Mapping[Any, Any],
) -> tuple[Any, Any]:
    packets = ((tail, "tail"), (geometry, "geometry"),
               (symbolic, "symbolic"), (connection, "connection"))
    family, logs = _source_family_logs(*packets[0])
    for packet, label in packets[1:]:
        other_family, other_logs = _source_family_logs(packet, label)
        if other_family != family:
            raise PTailM15InputError("source_family differs between saved packets")
        if other_logs != logs:
            raise PTailM15InputError("source_logs differs between saved packets")
    return family, logs


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _endpoint_record(value: Any, label: str) -> tuple[Any, str]:
    if _is_interval(value):
        return value, "direct_interval"
    record = _mapping(value, label)
    if "lower_exact_mpf_tuple" in record or "upper_exact_mpf_tuple" in record:
        return record, "direct"
    for key in ("interval", "bounds", "range", "value", "coefficient"):
        nested = record.get(key)
        if isinstance(nested, Mapping) and (
            "lower_exact_mpf_tuple" in nested
            or "upper_exact_mpf_tuple" in nested
        ):
            return nested, key
    raise PTailM15InputError(label + " must retain a binary interval record")


def _check_endpoint_shape(record: Any, label: str) -> None:
    if _is_interval(record):
        return
    endpoint = _mapping(record, label)
    for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        raw = endpoint.get(key)
        if (not isinstance(raw, (list, tuple)) or len(raw) != 4
                or any(type(part) is not int for part in raw)):
            raise PTailM15InputError(label + " must retain a four-integer endpoint")
    # Old saved logs omit this marker.  An explicit false marker is unsafe.
    if "directed" in endpoint and endpoint.get("directed") is not True:
        raise PTailM15InputError(label + ".directed must remain true when present")


def _interval(c: Any, value: Any, label: str) -> Any:
    if _is_interval(value):
        return value
    record, _path = _endpoint_record(value, label)
    if _is_interval(record):
        return record
    _check_endpoint_shape(record, label)
    lower = tuple(record["lower_exact_mpf_tuple"])
    upper = tuple(record["upper_exact_mpf_tuple"])
    interval = c.make_mpf((lower, upper))
    if not interval.a <= interval.b:
        raise PTailM15InputError(label + " has reversed endpoints")
    return interval


def _interval_provenance(value: Any, label: str) -> dict[str, Any]:
    record, path = _endpoint_record(value, label)
    if not _is_interval(record):
        _check_endpoint_shape(record, label)
    return {
        "wrapper": _copy(value),
        "endpoint_record": _copy(record),
        "endpoint_path": path,
        "directed_marker_present": bool(
            isinstance(record, Mapping) and "directed" in record
        ),
    }


def _upper_abs_point(c: Any, value: Any, label: str) -> Any:
    magnitude = abs(_interval(c, value, label))
    raw = magnitude._mpi_[1]
    return c.make_mpf((raw, raw))


def _symmetric(c: Any, value: Any) -> Any:
    radius = _upper_abs_point(c, value, "radius")
    negative = -radius
    return c.make_mpf((negative._mpi_[0], radius._mpi_[1]))


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if value is None or value == "":
            continue
        if not any(value == old for old in result):
            result.append(_copy(value))
    return result


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        return [value]
    return list(value)


def _ids(value: Any) -> list[Any]:
    if not isinstance(value, Mapping):
        return []
    result: list[Any] = []
    for key in ("source_ids", "source_id", "source_names", "source_name",
                "source_paths", "source_path", "term_ids", "term_id", "id"):
        result.extend(_items(value.get(key)))
    origin = value.get("origin")
    if isinstance(origin, Mapping):
        result.extend(_ids(origin))
    return _unique(result)


def _powers(term: Mapping[Any, Any], label: str) -> list[int]:
    raw = term.get("powers")
    if raw is None:
        raw = term.get("formal_source_powers")
    if raw is None:
        raw = term.get("source_powers")
    if (isinstance(raw, (str, bytes, bytearray)) or not isinstance(raw, Sequence)
            or len(raw) != 4 or any(type(power) is not int for power in raw)):
        raise PTailM15InputError(label + " must contain four integer powers")
    return [int(power) for power in raw]


def _output_powers(source: Sequence[int], *, ell: bool) -> list[int]:
    powers = [int(item) for item in source] + [0, 0]
    powers[3] -= 2
    powers[4] = 1
    if ell:
        powers[3] += 1
    powers[5] = 0
    return powers


def _term_interval_value(term: Mapping[Any, Any], c: Any, label: str) -> Any:
    for key in ("coefficient", "interval", "error", "bound"):
        if key in term:
            return _interval(c, term[key], label + "." + key)
    raise PTailM15InputError(label + " is missing its error interval")


def _term_records(value: Any, label: str) -> list[Mapping[Any, Any]]:
    if isinstance(value, Mapping):
        return [value]
    rows = _sequence(value, label)
    return [_mapping(row, label + "[" + str(index) + "]")
            for index, row in enumerate(rows)]


def _cell_map(value: Any, label: str) -> dict[int, Mapping[Any, Any]]:
    if isinstance(value, Mapping):
        rows = [_lookup(value, index, label) for index in range(CELL_COUNT)]
    else:
        rows = _sequence(value, label, allow_empty=False)
    if len(rows) != CELL_COUNT:
        raise PTailM15InputError(label + " must contain eight cells")
    result: dict[int, Mapping[Any, Any]] = {}
    for position, raw in enumerate(rows):
        cell = _mapping(raw, label + "[" + str(position) + "]")
        index = cell.get("index", cell.get("target_index",
                                             cell.get("target_cell", position)))
        if type(index) is not int or index != position:
            raise PTailM15InputError(label + " has noncanonical cell order")
        result[position] = cell
    return result


def _tail_cells(tail: Mapping[Any, Any], branch: str) -> list[Mapping[Any, Any]]:
    branches = _mapping(tail.get("branches"), "tail.branches")
    packet = _mapping(_lookup(branches, branch, "tail branch"), "tail branch packet")
    return list(_cell_map(packet.get("target_cells"),
                          "tail." + branch + ".target_cells").values())


def _tail_terms(tail_cell: Mapping[Any, Any], target: int, source: int,
                order: str, branch: str) -> list[Mapping[Any, Any]]:
    by_source = _mapping(tail_cell.get("by_source_cell"),
                         "tail." + branch + ".by_source_cell")
    source_packet = _mapping(_lookup(by_source, source, "tail source cell"),
                             "tail source packet")
    fields = _mapping(source_packet.get("fields"), "tail source packet.fields")
    order_packet = _mapping(_lookup(fields, order, "tail fields order"),
                            "tail " + order + " fields")
    if FIELD not in order_packet:
        raise PTailM15InputError("tail source packet is missing F1")
    return _term_records(order_packet[FIELD], "tail F1 terms")


def _validate_error_term(term: Mapping[Any, Any], *, target: int, source: int,
                         order: str, branch: str, position: int,
                         c: Any) -> tuple[Any, list[int], list[Any]]:
    label = ("tail." + branch + ".target_cells[" + str(target)
             + "].by_source_cell[" + str(source) + "].fields." + order
             + ".F1[" + str(position) + "]")
    if term.get("error_only") is not True:
        raise PTailM15InputError(label + " must have error_only=True")
    if term.get("symmetric_error", term.get("symmetric")) is not True:
        raise PTailM15InputError(label + " must have symmetric_error=True")
    if term.get("signed_field") is not False:
        raise PTailM15InputError(label + " must have signed_field=False")
    powers = _powers(term, label)
    coefficient = _term_interval_value(term, c, label)
    lower, upper = coefficient._mpi_
    if not (lower[1] == 0 and upper[1] == 0):
        if lower[0] != 1 or upper[0] != 0 or lower[1:] != upper[1:]:
            raise PTailM15InputError(label + " must retain symmetric binary endpoints")
    if "target_cell" in term and term.get("target_cell") != target:
        raise PTailM15InputError(label + " target_cell differs from its container")
    if "source_cell" in term and term.get("source_cell") != source:
        raise PTailM15InputError(label + " source_cell differs from its container")
    return coefficient, powers, _ids(term)


def _background_terms(cell: Mapping[Any, Any], name: str,
                      index: int) -> list[Mapping[Any, Any]]:
    terms_root = _mapping(cell.get("terms"),
                          "symbolic.cells[" + str(index) + "].terms")
    swirl = _mapping(terms_root.get("swirl"),
                     "symbolic.cells[" + str(index) + "].terms.swirl")
    found: list[Mapping[Any, Any]] = []
    for order in ORDERS:
        rows = _term_records(_lookup(swirl, order, "symbolic swirl " + order),
                             "symbolic swirl " + order)
        for row in rows:
            if row.get("background_factor_name") == name:
                background = row.get("background_factor")
                _endpoint_record(background, "symbolic " + name + " background")
                found.append(_mapping(background,
                                      "symbolic " + name + " background"))
    if not found:
        raise PTailM15InputError("symbolic cell " + str(index)
                                + " has no " + name + " background factor")
    return found


def _background_source_marker(value: Mapping[Any, Any], fallback: str) -> list[Any]:
    marker = _ids(value)
    return marker or [fallback]


def _collect_backgrounds(symbolic: Mapping[Any, Any]) -> dict[str, dict[str, Any]]:
    cells = _cell_map(symbolic.get("cells"), "symbolic.cells")
    output: dict[str, dict[str, Any]] = {}
    for name in ("Phi_Z0", "Phi_Z1"):
        base_wrapper: Mapping[Any, Any] | None = None
        base_record: Any = None
        base_marker: list[Any] | None = None
        records_by_cell: dict[str, Any] = {}
        for index in range(CELL_COUNT):
            candidates = _background_terms(cells[index], name, index)
            for candidate in candidates:
                marker = _background_source_marker(candidate, name)
                endpoint, _path = _endpoint_record(candidate,
                                                    "symbolic " + name)
                if base_wrapper is None:
                    base_wrapper, base_record, base_marker = candidate, endpoint, marker
                else:
                    if marker != base_marker:
                        raise PTailM15InputError(name + " source name differs")
                    if endpoint != base_record:
                        raise PTailM15InputError(name + " interval record differs")
            # The producer repeats one background record on every symbolic
            # term.  Keep a count and one representative per cell instead
            # of multiplying that redundant data by every new error term.
            records_by_cell[str(index)] = {
                "candidate_count": len(candidates), "record": _copy(candidates[0])
            }
        if base_wrapper is None or base_record is None or base_marker is None:
            raise PTailM15InputError("missing symbolic " + name + " background")
        output[name] = {
            "name": name,
            "wrapper": _copy(base_wrapper),
            "endpoint_record": _copy(base_record),
            "source_ids": _copy(base_marker),
            "records_by_cell": records_by_cell,
            "source_identity_proved": False,
            "source_identity_status": "same source name/record only; owner unverified",
        }
    return output


def _geometry_cells(geometry: Mapping[Any, Any]) -> dict[int, Mapping[Any, Any]]:
    verified = (geometry.get("coordinate_verified") is True
                or geometry.get("coordinate_correspondence_verified") is True)
    if not verified:
        raise PTailM15InputError("geometry coordinate verification is required")
    for key in ("cutoff_owner_verified", "same_live_cutoff_owner_verified",
                "cutoff_source_identity_proved"):
        if geometry.get(key) is True:
            raise PTailM15InputError("geometry cutoff same-source ownership is not admitted")
    return _cell_map(geometry.get("cells"), "geometry.cells")


def _swirl_weight(c: Any, geometry_cell: Mapping[Any, Any],
                  index: int) -> tuple[Any, dict[str, Any]]:
    weights = geometry_cell.get("weights")
    if not isinstance(weights, Mapping):
        raise PTailM15InputError("geometry cell " + str(index) + " is missing weights")
    record = weights.get("swirl", weights.get("swirl_weight"))
    if not isinstance(record, Mapping):
        raise PTailM15InputError("geometry cell " + str(index) + " is missing swirl weight")
    numerator, denominator = record.get("numerator"), record.get("denominator")
    if type(numerator) is not int or type(denominator) is not int or denominator == 0:
        raise PTailM15InputError("geometry cell " + str(index) + " has malformed swirl weight")
    if numerator == 0:
        positive_numerator, positive_denominator = 0, 1
    elif denominator < 0:
        positive_numerator, positive_denominator = abs(numerator), abs(denominator)
    else:
        positive_numerator, positive_denominator = abs(numerator), denominator
    rational = {
        "numerator": positive_numerator,
        "denominator": positive_denominator,
        "signed_input": _copy(record),
        "absolute": True,
        "source": "geometry.weights.swirl",
    }
    point = c.mpf(positive_numerator) / positive_denominator
    raw = point._mpi_
    return c.make_mpf((raw[0], raw[1])), rational


def _connection_payload(connection: Mapping[Any, Any]) -> tuple[Any, Any, dict[str, Any]]:
    f0 = connection.get("F0_connection")
    if not isinstance(f0, Mapping):
        f0 = connection
    ratio = f0.get("ell_over_lambda")
    if ratio is None and f0 is not connection:
        ratio = connection.get("ell_over_lambda")
    if ratio is None:
        raise PTailM15InputError("connection.F0_connection.ell_over_lambda is required")
    _endpoint_record(ratio, "connection.F0_connection.ell_over_lambda")
    provenance = {
        "path": "F0_connection.ell_over_lambda",
        "record": _copy(ratio),
        "F0_connection": {"ell_over_lambda": _copy(ratio)},
        "full_connection_record_location": "authenticated first-word receipt result.output.F0_connection",
        "connection_source": {
            "source_family": _copy(connection.get("source_family")),
            "source_logs": _copy(connection.get("source_logs")),
        },
        "source_ids": _unique(_ids(f0) + _ids(ratio) + ["saved_F0_connection:ell_over_lambda"]),
        "added_source_id_is_a_receipt_path_label_not_owner_proof": True,
        "source_identity_proved": False,
        "ratio_is_enclosure_not_value": True,
        "connection_recomputed": False,
    }
    return f0, ratio, provenance


def _radius(c: Any, factors: Sequence[Any]) -> Any:
    value = c.mpf(1)
    for factor in factors:
        value *= _upper_abs_point(c, factor, "product factor")
    return _symmetric(c, value)


def _make_output_term(
    c: Any, *, error: Mapping[Any, Any], error_value: Any,
    source_powers: Sequence[int], branch: str, target: int, source: int,
    order: str, route: str, position: int, phi_name: str,
    phi_meta: Mapping[str, Any], weight_value: Any,
    weight_meta: Mapping[str, Any], ell_value: Any | None,
    ell_meta: Mapping[str, Any] | None, coordinates: Any,
) -> dict[str, Any]:
    uses_ell = route == "Phi0_ell_E0"
    output_powers = _output_powers(source_powers, ell=uses_ell)
    input_ids = _ids(error)
    phi_ids = list(phi_meta.get("source_ids", []))
    ell_ids = list(ell_meta.get("source_ids", [])) if ell_meta is not None else []
    source_ids = _unique(input_ids + phi_ids + ell_ids)
    factors = ["abs(swirl_weight)", "abs(" + phi_name + ")"]
    if uses_ell:
        factors.append("abs(ell_over_lambda)")
    product_factors: list[Any] = [weight_value, phi_meta["endpoint_record"]]
    if uses_ell and ell_value is not None:
        product_factors.append(ell_value)
    product_factors.append(error_value)
    radius = _radius(c, product_factors)
    ordinary_order = 0 if order == "Z0" else 1
    input_order = "Z1" if route == "Phi0_E1" else "Z0"
    path = ("branches." + branch + ".target_cells[" + str(target)
            + "].by_source_cell[" + str(source) + "].fields." + input_order
            + ".F1[" + str(position) + "]")
    return {
        "coefficient": radius,
        "powers": _copy(output_powers),
        "formal_source_powers": _copy(output_powers),
        "formal_source_basis": list(OUTPUT_POWER_BASIS),
        "saved_formal_source_powers": _copy(list(source_powers)),
        "source_formal_source_powers": _copy(list(source_powers)),
        "original_error": _copy(error),
        "original_error_record": _copy(error),
        "input_error": _copy(error),
        "input_error_record": _copy(error),
        "input_term_position": position,
        "input_term_index": position,
        "input_order": input_order,
        "input_field": FIELD,
        "input_path": path,
        "input_source_ids": _copy(input_ids),
        "parent_source_ids": _copy(input_ids),
        "source_ids": source_ids,
        "source_ids_by_factor": {
            "tail_error": _copy(input_ids),
            "phi": _copy(phi_ids),
            "ell_over_lambda": _copy(ell_ids),
        },
        "error_classes": _copy(error.get("error_classes", [])),
        "input_error_classes": _copy(error.get("error_classes", [])),
        "target_cell": target,
        "source_cell": source,
        "target_index": target,
        "source_index": source,
        "target_source_ids": {
            "target_cell": target,
            "source_cell": source,
            "input_source_ids": _copy(input_ids),
        },
        "coordinates": _copy(coordinates),
        "ordinary_Z": order,
        "ordinary_Z_order": ordinary_order,
        "branch": branch,
        "route": route,
        "swirl_error_route": route,
        "factor_formula": "abs(swirl_weight)*abs(" + phi_name + ")*E" + ("1" if route == "Phi0_E1" else "0") + ("*abs(ell_over_lambda)" if uses_ell else ""),
        "factors": factors,
        "error_only": True,
        "symmetric_error": True,
        "signed_field": False,
        "uncertainty": True,
        "coefficient_is_error_enclosure": True,
        "bound_is_not_a_field_value": True,
        "no_summation": True,
        "no_cancellation": True,
        "no_known_signed_terms": True,
        "weight_applied_once": True,
        "no_extra_minus_two": True,
        "swirl_weight_includes_minus_16_delta_t_four": True,
        "geometry_weight": _copy(weight_meta),
        "positive_swirl_weight_factor": weight_value,
        "phi_provenance": _copy(phi_meta),
        "ell_provenance": _copy(ell_meta) if ell_meta is not None else None,
        "conditional_L_local_identity": uses_ell,
        "conditional_identity": "L_local = Lambda * ell_over_lambda" if uses_ell else None,
        "F0_anchor_local_independent": True,
        "source_identity_status": "conditional prior; source owner unverified",
        "source_identity_proved": False,
        "input_term_preserved": True,
    }


def build(c: Any, tail: Any, geometry: Any, symbolic: Any, connection: Any) -> dict[str, Any]:
    """Return separate conditional Z0/Z1 swirl error records for five routes."""

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide mpf and make_mpf")
    tail_packet = _payload(tail, "tail", "target_tail")
    geometry_packet = _payload(geometry, "geometry", "geometry")
    symbolic_packet = _payload(symbolic, "symbolic", "real_source_expressions")
    connection_root = _mapping(connection, "connection")
    # The historical output contains F0_connection and result as siblings;
    # selecting the latter would silently discard the saved connection.
    connection_packet = (connection_root if "F0_connection" in connection_root
                         else _payload(connection_root, "connection", "output"))
    family, logs = _bind_source_labels(tail_packet, geometry_packet,
                                       symbolic_packet, connection_packet)
    raw_branches = _mapping(tail_packet.get("branches"), "tail.branches")
    if set(raw_branches) != set(BRANCHES) or len(raw_branches) != len(BRANCHES):
        raise PTailM15InputError("tail must retain exactly the five canonical branches")
    geometry_cells = _geometry_cells(geometry_packet)
    backgrounds = _collect_backgrounds(symbolic_packet)
    _connection_record, ell_record, ell_provenance = _connection_payload(connection_packet)
    ell_value = _interval(c, ell_record, "connection.F0_connection.ell_over_lambda")
    output_branches: dict[str, Any] = {}
    total_counts = {order: 0 for order in ORDERS}
    for branch in BRANCHES:
        tail_cells = _tail_cells(tail_packet, branch)
        branch_cells: list[dict[str, Any]] = []
        branch_counts = {order: 0 for order in ORDERS}
        for target in range(CELL_COUNT):
            geometry_cell = geometry_cells[target]
            weight_value, weight_meta = _swirl_weight(c, geometry_cell, target)
            coordinate = geometry_cell.get("coordinates", geometry_cell.get("coordinate"))
            if coordinate is None:
                raise PTailM15InputError("geometry cell " + str(target) + " is missing coordinates")
            target_row: dict[str, Any] = {
                "target_index": target, "target_cell": target,
                "coordinates": _copy(coordinate), "geometry_cell_index": target,
                "by_source_cell": {},
            }
            for source in range(target + 1):
                source_row: dict[str, Any] = {
                    "target_index": target, "source_index": source,
                    "target_cell": target, "source_cell": source,
                    "coordinates": _copy(coordinate),
                    "fields": {order: {FIELD: []} for order in ORDERS},
                    "no_aggregation": True, "no_cancellation": True,
                    "known_signed_terms_separate": True,
                }
                z0_terms = _tail_terms(tail_cells[target], target, source, "Z0", branch)
                z1_terms = _tail_terms(tail_cells[target], target, source, "Z1", branch)
                for position, term in enumerate(z0_terms):
                    error_value, powers, _ids_from_error = _validate_error_term(
                        term, target=target, source=source, order="Z0",
                        branch=branch, position=position, c=c)
                    source_row["fields"]["Z0"][FIELD].append(_make_output_term(
                        c, error=term, error_value=error_value,
                        source_powers=powers, branch=branch, target=target,
                        source=source, order="Z0", route="Phi0_E0",
                        position=position, phi_name="Phi0",
                        phi_meta=backgrounds["Phi_Z0"], weight_value=weight_value,
                        weight_meta=weight_meta, ell_value=None, ell_meta=None,
                        coordinates=coordinate))
                for position, term in enumerate(z1_terms):
                    error_value, powers, _ids_from_error = _validate_error_term(
                        term, target=target, source=source, order="Z1",
                        branch=branch, position=position, c=c)
                    source_row["fields"]["Z1"][FIELD].append(_make_output_term(
                        c, error=term, error_value=error_value,
                        source_powers=powers, branch=branch, target=target,
                        source=source, order="Z1", route="Phi0_E1",
                        position=position, phi_name="Phi0",
                        phi_meta=backgrounds["Phi_Z0"], weight_value=weight_value,
                        weight_meta=weight_meta, ell_value=None, ell_meta=None,
                        coordinates=coordinate))
                # Each Z1 route is separate.  The E0 stream is reused only as
                # an input record; no route is added to another route.
                for position, term in enumerate(z0_terms):
                    error_value, powers, _ids_from_error = _validate_error_term(
                        term, target=target, source=source, order="Z0",
                        branch=branch, position=position, c=c)
                    source_row["fields"]["Z1"][FIELD].append(_make_output_term(
                        c, error=term, error_value=error_value,
                        source_powers=powers, branch=branch, target=target,
                        source=source, order="Z1", route="Phi1_E0",
                        position=position, phi_name="Phi1",
                        phi_meta=backgrounds["Phi_Z1"], weight_value=weight_value,
                        weight_meta=weight_meta, ell_value=None, ell_meta=None,
                        coordinates=coordinate))
                    source_row["fields"]["Z1"][FIELD].append(_make_output_term(
                        c, error=term, error_value=error_value,
                        source_powers=powers, branch=branch, target=target,
                        source=source, order="Z1", route="Phi0_ell_E0",
                        position=position, phi_name="Phi0",
                        phi_meta=backgrounds["Phi_Z0"], weight_value=weight_value,
                        weight_meta=weight_meta, ell_value=ell_value,
                        ell_meta=ell_provenance, coordinates=coordinate))
                branch_counts["Z0"] += len(source_row["fields"]["Z0"][FIELD])
                branch_counts["Z1"] += len(source_row["fields"]["Z1"][FIELD])
                target_row["by_source_cell"][str(source)] = source_row
            branch_cells.append(target_row)
        total_counts["Z0"] += branch_counts["Z0"]
        total_counts["Z1"] += branch_counts["Z1"]
        output_branches[branch] = {
            "branch": branch, "target_cells": branch_cells,
            "term_counts": branch_counts, "error_only": True,
            "symmetric_error": True, "signed_field": False,
            "no_aggregation": True, "no_cancellation": True,
            "known_signed_terms_separate": True, "d_x_F1_not_consumed": True,
        }
    return {
        "schema": SCHEMA, "status": "conditional_m15_swirl_tail_error",
        "source_family": _copy(family), "source_logs": _copy(logs),
        "source_identity_status": "source_family/source_logs match only; owner unverified",
        "source_identity_verified": False, "source_owner_identity_proved": False,
        "owner_verified": False,
        "geometry": {
            "coordinate_verified": True,
            "coordinate_correspondence_verified": geometry_packet.get("coordinate_correspondence_verified") is True,
            "cutoff_same_source_owner_proved": False,
            "cutoff_owner_verified": False, "source_identity_proved": False,
            "source_family_logs_match_only": True,
        },
        "background_provenance": {
            "Phi_Z0": _copy(backgrounds["Phi_Z0"]),
            "Phi_Z1": _copy(backgrounds["Phi_Z1"]),
            "same_source_name_records_checked": True,
            "source_identity_proved": False,
        },
        "connection_provenance": _copy(ell_provenance),
        "ell_over_lambda": _copy(ell_record),
        "source_power_basis": list(OUTPUT_POWER_BASIS),
        "formal_source_basis": list(OUTPUT_POWER_BASIS),
        "ordinary_Z_orders": list(ORDERS), "branches": output_branches,
        "branch_names": list(BRANCHES), "term_counts": total_counts,
        "weight_formula": "abs(geometry.weights.swirl.numerator/denominator)",
        "swirl_weight_includes_minus_16_delta_t_four": True,
        "no_extra_minus_two": True, "weights_applied_once": True,
        "Z0_formula": "abs(swirl_weight)*abs(Phi0)*E0",
        "Z1_formula": "abs(swirl_weight)*abs(Phi1)*E0; abs(swirl_weight)*abs(Phi0)*E1; abs(swirl_weight)*abs(Phi0)*abs(ell_over_lambda)*E0",
        "Z1_routes_separate": list(Z1_ROUTES),
        "coefficients_are_symmetric_error_intervals": True,
        "all_coefficients_symmetric_error_only": True,
        "coefficients_are_not_field_values": True, "original_errors_retained": True,
        "known_signed_terms_separate": True, "ledger_separate": True,
        "appended_to_ledger": False, "ledger_changed": False,
        "pressure_included": False, "collar_included": False,
        "full_word_included": False, "allword_tail_included": False,
        "complete_W1": False, "complete_m15": False,
        "fixed_x2_exit_computed": False, "terminal_closure_proved": False,
        "source_values_computed": False, "point_value_selected": False,
        "integration_performed": True,
        "integration_scope": "target-cell upper error enclosure using saved exact measure weight",
        "signed_field_integrated": False, "old_kernel_helper_called": False,
        "connection_recomputed": False, "conditional_prior_no_new_identity": True,
        "F0_connection_reused_without_recomputation": True,
        "F0_anchor_local_independent": True, "formal_powers_not_expanded": True,
        "directed_endpoint_policy": "saved endpoint tuples use c.make_mpf; symmetric radii use upper abs",
    }


propagate = build


__all__ = ["SCHEMA", "CELL_COUNT", "ORDERS", "FIELD", "BRANCHES",
           "SOURCE_POWER_BASIS", "OUTPUT_POWER_BASIS", "Z1_ROUTES",
           "PTailM15InputError", "build", "propagate"]
