"""Strict coordinate bridge for the saved Kchi and K0 cell packets.

This leaf is intentionally narrower than either producer.  It accepts the
already-bound ``original_kchi_first_saved_cells_v2`` packet and the
``original_k0_symbolic_m15_cells_v1`` packet, and checks only their coordinate
metadata.  It does not inspect coefficient records, import an old helper,
perform interval arithmetic, integrate a field, or make a live-owner claim.

The checked canonical map is

``Kchi x_i = [i/4, (i+1)/4]``;
``K0 t_i = [i/8, (i+1)/8]``;
``K0 x_i = 2*t_i``;
``K0 rho_i = x_i**2 = 4*t_i**2``.

The endpoint decoder is deliberately bounded.  It converts only small binary
coordinate endpoints to exact :class:`fractions.Fraction` values; an endpoint
with a large exponent or mantissa is rejected before any large shift is made.
The resulting weights are exact rational metadata for a later integration
leaf, not integrals or field values.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from fractions import Fraction
from typing import Any


SCHEMA = "original_kchi_k0_geometry_bridge_v1"
KCHI_SCHEMA = "original_kchi_first_saved_cells_v2"
K0_SCHEMA = "original_k0_symbolic_m15_cells_v1"
CELL_COUNT = 8

# Coordinate endpoints are small dyadics.  These caps make malformed or
# coefficient-sized tuples fail closed without allocating an exponent-sized
# integer.  No coefficient endpoint is ever passed to this decoder.
MAX_ABS_EXPONENT = 256
MAX_MANTISSA_BITS = 256
MAX_COUNT = 256

AXIAL_WEIGHT_FORMULA = "8*(t_hi^2-t_lo^2)"
SWIRL_WEIGHT_FORMULA = "-16*(t_hi^4-t_lo^4)"
RADIAL_MEASURE_FORMULA = "d_rho=8*t*dt=2*x*dx"
RHO_FORMULA = "rho=x^2=4*t^2"


class GeometryBridgeInputError(ValueError):
    """Raised when a saved coordinate packet is not safe to bridge."""


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise GeometryBridgeInputError(label + " must be a mapping")
    return value


def _sequence(value: Any, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise GeometryBridgeInputError(label + " must be a sequence")
    return value


def _binary_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    """Validate one bounded exact-mpf tuple without expanding its exponent."""

    if (
        not isinstance(value, list)
        or len(value) != 4
        or any(type(item) is not int for item in value)
    ):
        raise GeometryBridgeInputError(label + " must be a four-integer binary tuple")
    sign, mantissa, exponent, count = value
    if sign not in (0, 1) or mantissa < 0 or count < 0:
        raise GeometryBridgeInputError(label + " has invalid sign, mantissa, or count")
    if abs(exponent) > MAX_ABS_EXPONENT:
        raise GeometryBridgeInputError(label + " exponent exceeds the coordinate bound")
    if mantissa == 0:
        if sign != 0 or exponent != 0 or count != 0:
            raise GeometryBridgeInputError(label + " has an invalid zero tuple")
        return 0, 0, 0, 0
    if mantissa.bit_length() > MAX_MANTISSA_BITS:
        raise GeometryBridgeInputError(label + " mantissa exceeds the coordinate bound")
    if count != mantissa.bit_length() or count > MAX_COUNT:
        raise GeometryBridgeInputError(label + " has an invalid mantissa count")
    return sign, mantissa, exponent, count


def _tuple_fraction(value: tuple[int, int, int, int]) -> Fraction:
    sign, mantissa, exponent, _count = value
    if mantissa == 0:
        return Fraction(0)
    magnitude = Fraction(mantissa, 1)
    if exponent >= 0:
        magnitude *= 1 << exponent
    else:
        magnitude /= 1 << (-exponent)
    return -magnitude if sign else magnitude


def _endpoint_record(
    value: Any,
    label: str,
    *,
    require_directed: bool = True,
) -> tuple[dict[str, Any], Fraction, Fraction]:
    """Return a copied endpoint record and its bounded exact values."""

    record = _mapping(value, label)
    if require_directed and record.get("directed") is not True:
        raise GeometryBridgeInputError(label + " must retain directed=True")
    lower_tuple = _binary_tuple(record.get("lower_exact_mpf_tuple"), label + ".lower")
    upper_tuple = _binary_tuple(record.get("upper_exact_mpf_tuple"), label + ".upper")
    lower = _tuple_fraction(lower_tuple)
    upper = _tuple_fraction(upper_tuple)
    if lower > upper:
        raise GeometryBridgeInputError(label + " has reversed endpoints")
    return deepcopy(dict(record)), lower, upper


def _require_interval(
    value: Any,
    expected_lower: Fraction,
    expected_upper: Fraction,
    label: str,
    *,
    require_directed: bool = True,
) -> dict[str, Any]:
    record, lower, upper = _endpoint_record(
        value, label, require_directed=require_directed
    )
    if lower != expected_lower or upper != expected_upper:
        raise GeometryBridgeInputError(
            label + " is not the expected canonical interval"
        )
    return record


def _fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": int(value.numerator), "denominator": int(value.denominator)}


def _expected_t(index: int) -> tuple[Fraction, Fraction]:
    return Fraction(index, CELL_COUNT), Fraction(index + 1, CELL_COUNT)


def _expected_x(index: int) -> tuple[Fraction, Fraction]:
    lower, upper = _expected_t(index)
    return 2 * lower, 2 * upper


def _expected_rho(index: int) -> tuple[Fraction, Fraction]:
    lower, upper = _expected_t(index)
    return 4 * lower * lower, 4 * upper * upper


def _expected_source_indices(index: int) -> list[int]:
    return list(range(index + 1))


def _packet_header(
    packet: Any,
    *,
    schema: str,
    label: str,
) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Sequence[Any]]:
    root = _mapping(packet, label)
    if root.get("schema") != schema:
        raise GeometryBridgeInputError(label + " has the wrong schema")
    family = _mapping(root.get("source_family"), label + ".source_family")
    logs = _mapping(root.get("source_logs"), label + ".source_logs")
    cells = _sequence(root.get("cells"), label + ".cells")
    if len(cells) != CELL_COUNT:
        raise GeometryBridgeInputError(
            label + ".cells must contain exactly eight cells"
        )
    return root, family, logs, cells


def _check_kchi_cell(cell_value: Any, index: int) -> dict[str, Any]:
    label = "kchi_saved_cells.cells[{}]".format(index)
    cell = _mapping(cell_value, label)
    if (
        cell.get("index") != index
        or cell.get("target_cell") != index
        or cell.get("target_index") != index
    ):
        raise GeometryBridgeInputError(label + " has a noncanonical index")
    expected_sources = _expected_source_indices(index)
    if cell.get("source_indices") != expected_sources:
        raise GeometryBridgeInputError(label + " has the wrong source_indices")
    if cell.get("triangular_diagonal_source") != index:
        raise GeometryBridgeInputError(label + " has the wrong triangular diagonal")
    if (
        "source_j_le_target_i" in cell
        and cell.get("source_j_le_target_i") is not True
    ):
        raise GeometryBridgeInputError(label + " must retain source_j_le_target_i=True")

    x_lower, x_upper = _expected_x(index)
    x_record = _require_interval(cell.get("x"), x_lower, x_upper, label + ".x")
    target_x_record = _require_interval(
        cell.get("target_x"), x_lower, x_upper, label + ".target_x"
    )

    source_values = _sequence(cell.get("source_cells"), label + ".source_cells")
    if len(source_values) != index + 1:
        raise GeometryBridgeInputError(label + " has the wrong number of source cells")
    source_records: list[dict[str, Any]] = []
    for source_index, source_value in enumerate(source_values):
        source_label = label + ".source_cells[{}]".format(source_index)
        source = _mapping(source_value, source_label)
        if source.get("source_index") != source_index:
            raise GeometryBridgeInputError(source_label + " has a noncanonical source index")
        if source.get("target_index") != index:
            raise GeometryBridgeInputError(source_label + " has the wrong target index")
        if source.get("triangular") is not (source_index == index):
            raise GeometryBridgeInputError(source_label + " has the wrong triangular flag")
        source_lower, source_upper = _expected_x(source_index)
        source_x = _require_interval(
            source.get("source_x"),
            source_lower,
            source_upper,
            source_label + ".source_x",
        )
        source_target_x = _require_interval(
            source.get("target_x"),
            x_lower,
            x_upper,
            source_label + ".target_x",
        )
        source_records.append(
            {
                "source_index": source_index,
                "target_index": index,
                "triangular": source_index == index,
                "source_x": source_x,
                "target_x": source_target_x,
            }
        )

    return {
        "index": index,
        "x": x_record,
        "target_x": target_x_record,
        "source_indices": deepcopy(expected_sources),
        "triangular_diagonal_source": index,
        "source_j_le_target_i": True,
        "source_cells": source_records,
    }


def _check_k0_cell(cell_value: Any, index: int) -> dict[str, Any]:
    label = "k0_symbolic_output.cells[{}]".format(index)
    cell = _mapping(cell_value, label)
    if cell.get("index") != index or cell.get("target_cell") != index:
        raise GeometryBridgeInputError(label + " has a noncanonical index")
    expected_sources = _expected_source_indices(index)
    if cell.get("source_indices") != expected_sources:
        raise GeometryBridgeInputError(label + " has the wrong source_indices")
    if cell.get("triangular_diagonal_source") != index:
        raise GeometryBridgeInputError(label + " has the wrong triangular diagonal")

    t_lower, t_upper = _expected_t(index)
    x_lower, x_upper = _expected_x(index)
    rho_lower, rho_upper = _expected_rho(index)
    t_record = _require_interval(
        cell.get("t"), t_lower, t_upper, label + ".t", require_directed=False
    )
    x_record = _require_interval(
        cell.get("x"), x_lower, x_upper, label + ".x", require_directed=False
    )
    rho_record = _require_interval(
        cell.get("rho"),
        rho_lower,
        rho_upper,
        label + ".rho",
        require_directed=False,
    )

    if (x_lower, x_upper) != (2 * t_lower, 2 * t_upper):
        raise GeometryBridgeInputError(label + " violates x=2*t")
    if (rho_lower, rho_upper) != (x_lower * x_lower, x_upper * x_upper):
        raise GeometryBridgeInputError(label + " violates rho=x^2")

    return {
        "index": index,
        "t": t_record,
        "x": x_record,
        "rho": rho_record,
        "source_indices": deepcopy(expected_sources),
        "triangular_diagonal_source": index,
    }


def _weights(index: int) -> dict[str, Any]:
    t_lower, t_upper = _expected_t(index)
    rho_lower, rho_upper = _expected_rho(index)
    rho_delta = rho_upper - rho_lower
    axial = 8 * (t_upper * t_upper - t_lower * t_lower)
    swirl = -16 * (t_upper**4 - t_lower**4)
    return {
        "rho_delta": _fraction_record(rho_delta),
        "axial": _fraction_record(axial),
        "swirl": _fraction_record(swirl),
        "axial_formula": AXIAL_WEIGHT_FORMULA,
        "swirl_formula": SWIRL_WEIGHT_FORMULA,
        "radial_measure": RADIAL_MEASURE_FORMULA,
        "is_predefined_rational_weight": True,
        "is_integral": False,
    }


def build(kchi_saved_cells: Any, k0_symbolic_output: Any) -> dict[str, Any]:
    """Verify the eight-cell coordinate bridge and return safe metadata.

    The two packets are treated as already decoded caller inputs.  Only schema,
    source-family/log equality, cell indices, coordinate endpoint records and
    triangular source ownership are read.  In particular, no ``terms`` or
    coefficient field is accessed.
    """

    _kchi_root, kchi_family, kchi_logs, kchi_cells = _packet_header(
        kchi_saved_cells, schema=KCHI_SCHEMA, label="kchi_saved_cells"
    )
    _k0_root, k0_family, k0_logs, k0_cells = _packet_header(
        k0_symbolic_output, schema=K0_SCHEMA, label="k0_symbolic_output"
    )
    if kchi_family != k0_family:
        raise GeometryBridgeInputError("source_family records differ")
    if kchi_logs != k0_logs:
        raise GeometryBridgeInputError("source_logs records differ")

    cells: list[dict[str, Any]] = []
    weights_by_cell: dict[str, Any] = {}
    for index in range(CELL_COUNT):
        kchi = _check_kchi_cell(kchi_cells[index], index)
        k0 = _check_k0_cell(k0_cells[index], index)
        weight = _weights(index)
        weights_by_cell[str(index)] = deepcopy(weight)
        cells.append(
            {
                "index": index,
                "kchi": kchi,
                "k0": k0,
                "coordinates": {
                    "kchi_x": deepcopy(kchi["x"]),
                    "kchi_target_x": deepcopy(kchi["target_x"]),
                    "k0_t": deepcopy(k0["t"]),
                    "k0_x": deepcopy(k0["x"]),
                    "k0_rho": deepcopy(k0["rho"]),
                },
                "weights": weight,
                "provenance": {
                    "kchi_schema": KCHI_SCHEMA,
                    "k0_schema": K0_SCHEMA,
                    "kchi_cell_index": index,
                    "k0_cell_index": index,
                    "source_family_logs_equal": True,
                    "live_owner_verified": False,
                },
            }
        )

    return {
        "schema": SCHEMA,
        "status": "coordinate_weight_verified",
        "cell_count": CELL_COUNT,
        "cells": cells,
        "weights_by_cell": weights_by_cell,
        "weight_formulas": {
            "axial": AXIAL_WEIGHT_FORMULA,
            "swirl": SWIRL_WEIGHT_FORMULA,
        },
        "rho_formula": RHO_FORMULA,
        "radial_measure": RADIAL_MEASURE_FORMULA,
        "coordinate_correspondence_verified": True,
        "canonical_x_partition_verified": True,
        "canonical_t_partition_verified": True,
        "rho_equals_x_squared_verified": True,
        "source_indices_are_forcing_cells": True,
        "target_cells_are_rho_measure_cells": True,
        "source_family": deepcopy(dict(kchi_family)),
        "source_logs": deepcopy(dict(kchi_logs)),
        "source_family_logs_equal": True,
        "source_identity_status": "source_family/source_logs match; live owner unverified",
        "source_identity_verified": False,
        "live_owner_verified": False,
        "cutoff_owner_verified": False,
        "coefficient_records_read": False,
        "interval_arithmetic_performed": False,
        "integration_performed": False,
        "integrals_claimed": False,
        "ledger_consulted": False,
        "ledger_changed": False,
        "W1_complete": False,
        "m15_complete": False,
        "complete_W1": False,
        "complete_m15": False,
        "scope": "coordinate and predefined target-rho weight verification only",
    }


bind = build
verify = build


__all__ = [
    "SCHEMA",
    "KCHI_SCHEMA",
    "K0_SCHEMA",
    "CELL_COUNT",
    "MAX_ABS_EXPONENT",
    "MAX_MANTISSA_BITS",
    "MAX_COUNT",
    "GeometryBridgeInputError",
    "build",
    "bind",
    "verify",
]
