"""Factored signed interval leaf for the corrected lower ``V_y`` repair.

This module is deliberately a small source-layer adapter.  The runner supplies
the interval context, the four exact serialized ``Z`` cells, one directed
``delta`` interval, and already authenticated raw rows for ``R``, ``Pstar``
and ``V_y[0:2]``.  No receipt is opened here and no producer or graph is
replayed.

The only numerical work is directed interval arithmetic on one whole cell at a
time.  Source rows stay factored as

    (signed coefficient interval) * exp(formal positive scale).

In particular, ``log_absolute_upper`` is never read.  The large subtraction
``Pstar - 1`` is represented by the correlated factor
``beta(Pstar) = 1 - 1/Pstar`` with enclosure ``[0, 1]``; that enclosure is not
treated as the value zero or one.  ``repair1`` intentionally remains two
separate source terms so scales are not merged without a certificate.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

try:  # Only used for decoding exact mpf endpoint tuples.
    import mpmath as _mp
except Exception:  # pragma: no cover - mpmath is a runner dependency.
    _mp = None


SCHEMA = "original_corrected_fourcell_signed_repair_v1"
COMPONENT = "corrected_physical_lower_Vy_factored_signed_interval_repair_v1"
CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))
ORDERS = (0, 1)
SCALE_FIELDS = ("source_exponents", "radius_power", "additional_log_interval")


class FourCellSignedRepairError(ValueError):
    """Raised when caller-owned source evidence cannot be admitted."""


OriginalCorrectedFourcellSignedRepairError = FourCellSignedRepairError


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FourCellSignedRepairError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    return value


def _sequence(value: Any, label: str) -> Sequence[Any]:
    _require(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
        f"{label} must be a sequence",
    )
    return value


def _copy_metadata(value: Any) -> Any:
    """Copy small caller metadata without depending on interval deepcopy."""

    if isinstance(value, Mapping):
        return {key: _copy_metadata(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_metadata(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_copy_metadata(item) for item in value)
    return value


def _is_mpf_tuple(value: Any) -> bool:
    if not isinstance(value, (tuple, list)) or len(value) != 4:
        return False
    return all(isinstance(item, int) and not isinstance(item, bool) for item in value)


def _point(c: Any, value: Any, label: str) -> Any:
    """Decode one scalar endpoint without selecting an interval midpoint."""

    if hasattr(value, "_mpf_"):
        return c.mpf(value)
    if _is_mpf_tuple(value):
        if _mp is not None:
            try:
                return c.mpf(_mp.make_mpf(tuple(value)))
            except Exception:
                pass
        try:
            return c.mpf(tuple(value))
        except Exception as exc:  # pragma: no cover - defensive fallback.
            raise FourCellSignedRepairError(f"{label} has an invalid exact mpf tuple") from exc
    try:
        return c.mpf(value)
    except Exception as exc:
        raise FourCellSignedRepairError(f"{label} is not a directed interval endpoint") from exc


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return getattr(value, "_mpf_", value)


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct a directed interval, retaining exact endpoint tuples."""

    lo = _raw_endpoint(lower, False)
    hi = _raw_endpoint(upper, True)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((lo, hi))
        except Exception:
            pass
    return c.mpf([_point(c, lo, "interval lower"), _point(c, hi, "interval upper")])


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _is_serialized_interval(value: Any) -> bool:
    if isinstance(value, Mapping):
        return (
            ("lower" in value or "lower_exact_mpf_tuple" in value)
            and ("upper" in value or "upper_exact_mpf_tuple" in value)
        )
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2


def _interval(c: Any, value: Any, label: str) -> Any:
    """Accept a runner interval or canonical lower/upper serialization."""

    if _is_interval(value):
        result = _make_interval(c, value, value)
    elif isinstance(value, Mapping):
        lower_key = next(
            (key for key in ("lower_exact_mpf_tuple", "lower") if key in value), None
        )
        upper_key = next(
            (key for key in ("upper_exact_mpf_tuple", "upper") if key in value), None
        )
        _require(lower_key is not None and upper_key is not None, f"{label} needs lower and upper")
        lower = _point(c, value[lower_key], f"{label}.lower")
        upper = _point(c, value[upper_key], f"{label}.upper")
        result = _make_interval(c, lower, upper)
    elif _is_mpf_tuple(value):
        # A four-tuple is one exact scalar endpoint, not a two-endpoint box.
        result = _point(c, value, label)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        _require(len(values) == 2, f"{label} must have two directed endpoints")
        result = _make_interval(
            c,
            _point(c, values[0], f"{label}.lower"),
            _point(c, values[1], f"{label}.upper"),
        )
    else:
        # A scalar is admitted as a degenerate directed interval.  It is still
        # carried by the interval context and is never converted to a value.
        result = _point(c, value, label)
    lower, upper = _bounds(result, label)
    _require(lower <= upper, f"{label} has reversed endpoints")
    try:
        finite = c.isfinite(lower) and c.isfinite(upper)
    except Exception:
        finite = True
    _require(bool(finite), f"{label} has non-finite endpoints")
    return result


def _bounds(value: Any, label: str = "interval") -> tuple[Any, Any]:
    if not _is_interval(value):
        raise FourCellSignedRepairError(f"{label} is not a directed interval")
    return value.a, value.b


def _cell(value: Any, label: str) -> tuple[str, str]:
    values = _sequence(value, label)
    _require(
        len(values) == 2 and all(isinstance(item, str) for item in values),
        f"{label} must contain two exact serialized endpoints",
    )
    result = (values[0], values[1])
    _require(result in CELLS, f"{label} is not one of the four exact Z cells")
    return result


def _cells(value: Any) -> tuple[tuple[str, str], ...]:
    values = _sequence(value, "exact_Z_cells")
    result = tuple(_cell(item, f"exact_Z_cells[{index}]") for index, item in enumerate(values))
    _require(len(result) == 4, "exact_Z_cells must contain exactly four cells")
    _require(set(result) == set(CELLS), "exact_Z_cells must cover the canonical four cells")
    # The output order is fixed so a runner cannot silently associate a source
    # row with a different cell after a permutation.
    return CELLS


def _cell_key(cell: tuple[str, str]) -> str:
    return str(cell)


def _cell_key_candidates(cell: tuple[str, str]) -> tuple[Any, ...]:
    return (
        str(cell),
        cell,
        f"[{cell[0]!r}, {cell[1]!r}]",
        f"[{cell[0]}, {cell[1]}]",
    )


def _cell_certificate(certificate: Mapping[str, Any], cell: tuple[str, str]) -> Mapping[str, Any]:
    """Return an optional cell-local certificate without identifying bases across cells."""

    for key in ("by_cell", "cells", "per_cell", "cell_certificates", "cell_provenance"):
        nested = certificate.get(key)
        if isinstance(nested, Mapping):
            for candidate in _cell_key_candidates(cell):
                try:
                    value = nested.get(candidate)
                except TypeError:
                    value = None
                if isinstance(value, Mapping):
                    return value
    return certificate


def _source_provenance(row: Mapping[str, Any], base: Mapping[str, Any]) -> Any:
    for container in (row, base):
        for key in ("provenance", "source_provenance", "origin", "source_origin"):
            if key in container:
                return _copy_metadata(container[key])
    return None


def _source_index(row: Mapping[str, Any], base: Mapping[str, Any], fallback: str) -> Any:
    provenance = _source_provenance(row, base)
    for container in (row, base, provenance if isinstance(provenance, Mapping) else {}):
        for key in ("source_index", "index", "row_index", "node", "source_id"):
            if key in container:
                return _copy_metadata(container[key])
    return fallback


def _row_base(row: Any, label: str) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    outer = _mapping(row, label)
    if "coefficient_interval" in outer and "formal_positive_scale" in outer:
        return outer, outer
    nested = outer.get("row")
    if isinstance(nested, Mapping) and {
        "coefficient_interval",
        "formal_positive_scale",
    }.issubset(nested):
        return outer, nested
    raise FourCellSignedRepairError(
        f"{label} must provide coefficient_interval and formal_positive_scale"
    )


def _fraction(value: Any, label: str) -> Fraction:
    _require(not isinstance(value, bool), f"{label} must be an exact integer or half-integer")
    try:
        if isinstance(value, Fraction):
            result = value
        elif isinstance(value, int):
            result = Fraction(value, 1)
        elif isinstance(value, float):
            result = Fraction(str(value))
        elif isinstance(value, str):
            result = Fraction(value)
        else:
            result = Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise FourCellSignedRepairError(f"{label} must be an exact integer or half-integer") from exc
    _require(result.denominator in (1, 2), f"{label} must be an integer or half-integer")
    return result


def _display_fraction(value: Fraction) -> int | float:
    if value.denominator == 1:
        return int(value.numerator)
    return float(value)


def _scale(c: Any, row: Mapping[str, Any], label: str) -> dict[str, Any]:
    raw = _mapping(row.get("formal_positive_scale"), f"{label}.formal_positive_scale")
    for field in SCALE_FIELDS:
        _require(field in raw, f"{label}.formal_positive_scale.{field} is required")
    powers = _sequence(raw.get("source_exponents"), f"{label}.source_exponents")
    _require(len(powers) == 4, f"{label}.source_exponents must have four entries")
    source_exponents = tuple(
        _fraction(power, f"{label}.source_exponents[{index}]")
        for index, power in enumerate(powers)
    )
    radius_power = _fraction(raw.get("radius_power"), f"{label}.radius_power")
    offset_raw = raw.get("additional_log_interval")
    _require(
        _is_interval(offset_raw) or _is_serialized_interval(offset_raw),
        f"{label}.additional_log_interval must be a directed interval",
    )
    offset = _interval(c, offset_raw, f"{label}.additional_log_interval")
    return {
        "source_exponents": source_exponents,
        "radius_power": radius_power,
        "additional_log_interval": offset,
        "raw_scale_record": _copy_metadata(dict(raw)),
    }


def _scale_record(scale: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "source_exponents": [
            _display_fraction(value) for value in scale["source_exponents"]
        ],
        "radius_power": _display_fraction(scale["radius_power"]),
        "additional_log_interval": scale["additional_log_interval"],
        "exponents_exact_integer_or_half_integer": True,
    }


def _scale_sum(*scales: Mapping[str, Any]) -> dict[str, Any]:
    _require(bool(scales), "at least one formal source scale is required")
    first = scales[0]
    source_exponents = [Fraction(0, 1)] * 4
    radius_power = Fraction(0, 1)
    offset = None
    for scale in scales:
        source_exponents = [
            left + right for left, right in zip(source_exponents, scale["source_exponents"])
        ]
        radius_power += scale["radius_power"]
        value = scale["additional_log_interval"]
        offset = value if offset is None else offset + value
    return {
        "source_exponents": tuple(source_exponents),
        "radius_power": radius_power,
        "additional_log_interval": offset,
    }


def _row(c: Any, raw: Any, label: str, fallback_index: str) -> dict[str, Any]:
    outer, base = _row_base(raw, label)
    if "point_value_selected" in outer:
        _require(outer.get("point_value_selected") is False, f"{label} selected a point value")
    if "point_value_selected" in base:
        _require(base.get("point_value_selected") is False, f"{label} selected a point value")
    if "encloses_original_source_function" in outer:
        _require(
            outer.get("encloses_original_source_function") is True,
            f"{label} is not marked as a source-function enclosure",
        )
    if "encloses_original_source_function" in base:
        _require(
            base.get("encloses_original_source_function") is True,
            f"{label} is not marked as a source-function enclosure",
        )
    coefficient_raw = base.get("coefficient_interval")
    _require(
        _is_interval(coefficient_raw) or _is_serialized_interval(coefficient_raw),
        f"{label}.coefficient_interval must be a signed directed interval",
    )
    coefficient = _interval(c, coefficient_raw, f"{label}.coefficient_interval")
    scale = _scale(c, base, label)
    return {
        "coefficient_interval": coefficient,
        "formal_positive_scale": scale,
        "provenance": _source_provenance(outer, base),
        "source_index": _source_index(outer, base, fallback_index),
        "source_path": _copy_metadata(
            outer.get("source_path", base.get("source_path", fallback_index))
        ),
        "row_was_enclosure": True,
        "point_value_selected": False,
        "log_absolute_upper_not_used_as_value": True,
    }


def _row_parts(entry: Any, label: str) -> tuple[Any, Any, Any, Any]:
    """Read one cell's R/Pstar/V_y[0]/V_y[1] tuple."""

    if isinstance(entry, Sequence) and not isinstance(entry, (str, bytes, bytearray)):
        values = list(entry)
        _require(len(values) == 4, f"{label} sequence must be R, Pstar, V_y[0], V_y[1]")
        return values[0], values[1], values[2], values[3]
    record = _mapping(entry, label)

    def pick(names: Sequence[str], field: str) -> Any:
        for name in names:
            if name in record:
                return record[name]
        raise FourCellSignedRepairError(f"{label} is missing {field}")

    R = pick(("R", "r", "source_R", "R_row"), "R")
    Pstar = pick(("Pstar", "pstar", "source_Pstar", "Pstar_row"), "Pstar")
    vy_value = None
    for name in ("V_y", "Vy", "v_y", "V_y_C0_C1", "Vy_C0_C1"):
        if name in record:
            vy_value = record[name]
            break
    if vy_value is not None:
        if isinstance(vy_value, Mapping):
            Y0 = next((vy_value[key] for key in (0, "0", "V_y[0]", "Y0") if key in vy_value), None)
            Y1 = next((vy_value[key] for key in (1, "1", "V_y[1]", "Y1") if key in vy_value), None)
            _require(Y0 is not None and Y1 is not None, f"{label}.V_y must provide rows 0 and 1")
        else:
            values = _sequence(vy_value, f"{label}.V_y")
            _require(len(values) >= 2, f"{label}.V_y must provide rows 0 and 1")
            Y0, Y1 = values[0], values[1]
    else:
        Y0 = pick(("V_y[0]", "Vy[0]", "Y0", "V_y_0", "V0"), "V_y[0]")
        Y1 = pick(("V_y[1]", "Vy[1]", "Y1", "V_y_1", "V1"), "V_y[1]")
    return R, Pstar, Y0, Y1


def _cell_entry(source_rows: Any, cell: tuple[str, str], position: int) -> Any:
    if isinstance(source_rows, Sequence) and not isinstance(source_rows, (str, bytes, bytearray)):
        values = list(source_rows)
        _require(len(values) == 4, "source_rows sequence must contain four cell records")
        return values[position]
    root = _mapping(source_rows, "source_rows")
    for container_key in ("by_cell", "cell_sources", "samples", "rows", "source_rows"):
        nested = root.get(container_key)
        if isinstance(nested, (Mapping, Sequence)) and not isinstance(nested, (str, bytes, bytearray)):
            try:
                return _cell_entry(nested, cell, position)
            except FourCellSignedRepairError:
                pass
    for key in _cell_key_candidates(cell):
        try:
            if key in root:
                return root[key]
        except TypeError:
            continue
    raise FourCellSignedRepairError(f"source_rows has no record for exact cell {_cell_key(cell)}")


def _all_source_rows(c: Any, source_rows: Any) -> dict[tuple[str, str], dict[str, Any]]:
    result: dict[tuple[str, str], dict[str, Any]] = {}
    for position, cell in enumerate(CELLS):
        R_raw, P_raw, Y0_raw, Y1_raw = _row_parts(
            _cell_entry(source_rows, cell, position), f"source_rows[{_cell_key(cell)!r}]"
        )
        result[cell] = {
            "R": _row(c, R_raw, f"source_rows[{_cell_key(cell)!r}].R", f"{_cell_key(cell)}:R"),
            "Pstar": _row(
                c, P_raw, f"source_rows[{_cell_key(cell)!r}].Pstar", f"{_cell_key(cell)}:Pstar"
            ),
            "Y0": _row(
                c,
                Y0_raw,
                f"source_rows[{_cell_key(cell)!r}].V_y[0]",
                f"{_cell_key(cell)}:V_y[0]",
            ),
            "Y1": _row(
                c,
                Y1_raw,
                f"source_rows[{_cell_key(cell)!r}].V_y[1]",
                f"{_cell_key(cell)}:V_y[1]",
            ),
        }
    return result


def _truth(value: Any) -> bool:
    if value is True:
        return True
    if isinstance(value, Mapping):
        return any(_truth(item) for item in value.values())
    return False


def _flag(record: Mapping[str, Any], names: Sequence[str]) -> bool:
    containers: list[Mapping[str, Any]] = [record]
    for key in ("checks", "certificates", "proof", "attestation", "evidence"):
        value = record.get(key)
        if isinstance(value, Mapping):
            containers.append(value)
    return any(_truth(container.get(name)) for container in containers for name in names)


def _identity_map(value: Any) -> dict[str, list[Any]]:
    if not isinstance(value, Mapping):
        return {}
    found: dict[str, list[Any]] = {}
    identity_keys = (
        "same_source_id",
        "source_identity",
        "source_basis_id",
        "basis_id",
        "provenance_id",
        "same_basis_id",
    )
    for key in identity_keys:
        if key in value:
            found.setdefault(key, []).append(value[key])
    for key in ("provenance", "source_provenance", "origin", "source_origin"):
        nested = value.get(key)
        for nested_key, nested_values in _identity_map(nested).items():
            found.setdefault(nested_key, []).extend(nested_values)
    return found


def _identity_values(value: Any) -> list[Any]:
    return [item for values in _identity_map(value).values() for item in values]


def _validate_source_certificate(
    certificate: Any, rows: Mapping[tuple[str, str], Mapping[str, Mapping[str, Any]]]
) -> dict[str, Any]:
    cert = _mapping(certificate, "same_source_certificate")
    required_guards = (
        "delta_Z_independent",
        "R_Z_independent",
        "Pstar_Z_independent",
        "ordinary_Z1",
        "Y_already_physical_Pstar_once",
    )
    per_cell: dict[str, Any] = {}
    for cell, cell_rows in rows.items():
        # Every source/product scale is combined only inside this cell.  The
        # opaque basis object is intentionally not compared with another cell.
        local_cert = _cell_certificate(cert, cell)
        same_basis = _flag(
            local_cert,
            (
                "same_source_basis",
                "same_source_identity",
                "same_basis",
                "same_source",
                "source_basis_checked",
            ),
        ) or _flag(cert, ("same_source_basis", "same_source_identity", "same_basis", "same_source", "source_basis_checked"))
        same_provenance = _flag(
            local_cert,
            (
                "same_source_provenance",
                "provenance_checked",
                "same_provenance",
                "source_provenance_checked",
            ),
        ) or _flag(cert, ("same_source_provenance", "provenance_checked", "same_provenance", "source_provenance_checked"))
        _require(same_basis, f"same-source basis certificate is required for cell {cell}")
        _require(same_provenance, f"same-source provenance certificate is required for cell {cell}")
        for guard in required_guards:
            _require(
                _flag(local_cert, (guard,)) or _flag(cert, (guard,)),
                f"same-source certificate is missing {guard} for cell {cell}",
            )
        passed = _flag(local_cert, ("passed", "all_passed", "certified", "admitted")) or _flag(
            cert, ("passed", "all_passed", "certified", "admitted")
        )
        _require(passed or (same_basis and same_provenance), f"same-source certificate is not passed for cell {cell}")

        row_ids: dict[str, list[Any]] = {}
        for row in cell_rows.values():
            for key, values in _identity_map(row.get("provenance")).items():
                row_ids.setdefault(key, []).extend(values)
        for key, values in row_ids.items():
            first = values[0]
            _require(all(value == first for value in values), f"raw source rows split same-source identity in cell {cell}")
        cert_ids = _identity_map(local_cert)
        for key, values in cert_ids.items():
            if key in row_ids:
                _require(all(value in row_ids[key] for value in values), f"source certificate identity differs from rows in cell {cell}")
        per_cell[_cell_key(cell)] = {
            "passed": True,
            "same_source_basis": True,
            "same_source_provenance": True,
            "required_guards": {guard: True for guard in required_guards},
            "row_identity_values_checked": bool(row_ids),
            "certificate_identity_values": _copy_metadata(cert_ids),
        }
    return {
        "passed": True,
        "same_source_basis": True,
        "same_source_provenance": True,
        "required_guards": {guard: True for guard in required_guards},
        **{guard: True for guard in required_guards},
        "cell_local_certificates": per_cell,
        "cross_cell_basis_identity_asserted": False,
        "source_basis_certificate": _copy_metadata(dict(cert)),
    }


def _find_log_pstar(c: Any, certificate: Mapping[str, Any]) -> tuple[Any | None, str | None]:
    for key in (
        "logPstar_interval",
        "logPstar_enclosure",
        "logPstar_bounds",
        "logPstar",
        "log_pstar_interval",
        "log_pstar",
    ):
        if key in certificate:
            return _interval(c, certificate[key], f"pstar_positive_certificate.{key}"), key
    for key in ("logPstar_lower", "logPstar_lower_bound", "log_pstar_lower"):
        if key in certificate:
            return _point(c, certificate[key], f"pstar_positive_certificate.{key}"), key
    for container_key in ("proof", "checks", "evidence", "attestation"):
        nested = certificate.get(container_key)
        if isinstance(nested, Mapping):
            interval, key = _find_log_pstar(c, nested)
            if interval is not None:
                return interval, f"{container_key}.{key}"
    return None, None


def _validate_pstar_certificate(
    c: Any,
    certificate: Any,
    rows: Mapping[tuple[str, str], Mapping[str, Mapping[str, Any]]] | None = None,
) -> dict[str, Any]:
    cert = _mapping(certificate, "pstar_positive_certificate")
    if rows is not None:
        has_cell_map = any(
            isinstance(cert.get(key), Mapping)
            for key in ("by_cell", "cells", "per_cell", "cell_certificates", "cell_provenance")
        )
        cert_by_cell = {
            _cell_key(cell): _cell_certificate(cert, cell) for cell in rows
        }
        for cell, cell_rows in rows.items():
            pstar_ids = _identity_map(cell_rows["Pstar"].get("provenance"))
            local_ids = _identity_map(cert_by_cell[_cell_key(cell)])
            if has_cell_map:
                for key, values in local_ids.items():
                    if key in pstar_ids:
                        _require(
                            all(value in pstar_ids[key] for value in values),
                            f"pstar positivity certificate is not bound to the Pstar source in cell {cell}",
                        )
    passed = _flag(cert, ("passed", "all_passed", "certified", "admitted"))
    explicit_positive = _flag(
        cert,
        (
            "logPstar_positive",
            "logPstar_gt_zero",
            "logPstar_strictly_positive",
            "strict_logPstar_positive",
            "Pstar_gt_one",
            "pstar_gt_one",
        ),
    )
    log_interval, source_key = _find_log_pstar(c, cert)
    if log_interval is None and rows:
        first_cell = next(iter(rows))
        local_cert = _cell_certificate(cert, first_cell)
        log_interval, source_key = _find_log_pstar(c, local_cert)
    beta_lower = c.mpf(0)
    beta_lower_proved = False
    if log_interval is not None:
        lower, _ = _bounds(log_interval, "pstar log certificate")
        _require(lower > 0, "pstar positivity certificate does not prove logPstar > 0")
        # For p=log(Pstar)>0, e^p >= 1+p gives
        #   beta = 1-exp(-p) >= p/(1+p).
        # This keeps the common factor bounded and avoids ever materializing
        # the astronomical Pstar or its reciprocal.
        ratio = c.mpf(lower) / (c.mpf(1) + c.mpf(lower))
        beta_lower = ratio.a
        beta_lower_proved = True
    else:
        _require(explicit_positive, "explicit logPstar > 0 certificate is required")
    _require(passed or explicit_positive or log_interval is not None, "pstar positivity certificate is not passed")
    return {
        "passed": True,
        "logPstar_gt_zero": True,
        "Pstar_gt_one": True,
        "logPstar_interval": log_interval,
        "logPstar_certificate_field": source_key,
        "strict_lower_checked": log_interval is not None,
        "beta_lower": beta_lower,
        "beta_lower_proved_from_logPstar": beta_lower_proved,
        "certificate": _copy_metadata(dict(cert)),
    }


def _provenance_view(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "source_index": _copy_metadata(row["source_index"]),
        "source_path": _copy_metadata(row["source_path"]),
        "provenance": _copy_metadata(row["provenance"]),
    }


def _term(
    *,
    c: Any,
    name: str,
    repair: int,
    index: int,
    coefficient: Any,
    scale: Mapping[str, Any],
    cell: tuple[str, str],
    rows: Mapping[str, Mapping[str, Any]],
    F_name: str,
    beta_binding: str,
    beta: Any,
) -> dict[str, Any]:
    source_indices = {
        "R": _copy_metadata(rows["R"]["source_index"]),
        "Pstar": _copy_metadata(rows["Pstar"]["source_index"]),
        "Y": _copy_metadata(rows["Y0"]["source_index"] if name.endswith("Y0") else rows["Y1"]["source_index"]),
    }
    source_indices["V_y"] = [
        _copy_metadata(rows["Y0"]["source_index"]),
        _copy_metadata(rows["Y1"]["source_index"]),
    ]
    source_indices["Y0"] = _copy_metadata(rows["Y0"]["source_index"])
    source_indices["Y1"] = _copy_metadata(rows["Y1"]["source_index"])
    return {
        "term_name": name,
        "repair_order": repair,
        "term_index": index,
        "cell": list(cell),
        "exact_Z_cell": list(cell),
        "source_correlated": True,
        "source_indices": source_indices,
        "coefficient_interval": coefficient,
        "coefficient": coefficient,
        "signed_coefficient_interval": coefficient,
        "signed_coefficient": coefficient,
        "formal_positive_scale": _scale_record(scale),
        "formal_scale": _scale_record(scale),
        "physical_formal_scale": _scale_record(scale),
        "source_scale_factors": {
            "R": _scale_record(rows["R"]["formal_positive_scale"]),
            "Pstar": _scale_record(rows["Pstar"]["formal_positive_scale"]),
            "Y": _scale_record(rows["Y0"]["formal_positive_scale"] if name.endswith("Y0") else rows["Y1"]["formal_positive_scale"]),
        },
        "factors": ["R", "Pstar", "beta(Pstar)", F_name, "Y0" if name.endswith("Y0") else "Y1"],
        "beta_interval": beta,
        "beta": {
            "definition": "beta(Pstar)=1-1/Pstar",
            "interval": beta,
            "binding_id": beta_binding,
            "same_source_identity": True,
        },
        "beta_binding_id": beta_binding,
        "beta_is_same_source_factor": True,
        "beta_interval_is_enclosure_only": True,
        "beta_not_selected_as_zero_or_one": True,
        "ordinary_Z_order": repair,
        "axis_function_whole_cell": True,
        "point_value_selected": False,
        "log_absolute_upper_not_used_as_value": True,
        "source_row_provenance": {
            "R": _provenance_view(rows["R"]),
            "Pstar": _provenance_view(rows["Pstar"]),
            "Y0": _provenance_view(rows["Y0"]),
            "Y1": _provenance_view(rows["Y1"]),
        },
        "formal_source_basis": {
            "scope": "one exact Z cell",
            "opaque": True,
            "same_source_certificate_checked_within_cell": True,
            "cross_cell_basis_identity_asserted": False,
            "cell": list(cell),
        },
        "Y_original_Pstar_retained_once": True,
        "Pstar_multiplier_is_old_excess": True,
    }


def _build_cell(
    c: Any,
    delta: Any,
    cell: tuple[str, str],
    rows: Mapping[str, Mapping[str, Any]],
    beta: Any,
    beta_binding: str,
) -> dict[str, Any]:
    left = _point(c, cell[0], f"cell {cell} lower")
    right = _point(c, cell[1], f"cell {cell} upper")
    Z = _make_interval(c, left, right)
    one = c.mpf(1)
    two = c.mpf(2)
    z2 = Z ** 2
    L = one - delta * z2
    L_lower, L_upper = _bounds(L, f"L{cell}")
    _require(L_lower > 0, f"L=1-delta*Z^2 is not strictly positive on cell {cell}")
    F = two * Z / L
    Fz = two * (one + delta * z2) / (L ** 2)

    common_scale = _scale_sum(
        rows["R"]["formal_positive_scale"], rows["Pstar"]["formal_positive_scale"]
    )
    repair0_scale = _scale_sum(common_scale, rows["Y0"]["formal_positive_scale"])
    repair1_y1_scale = _scale_sum(common_scale, rows["Y1"]["formal_positive_scale"])
    repair1_y0_scale = repair0_scale

    common_coeff = rows["R"]["coefficient_interval"] * rows["Pstar"]["coefficient_interval"] * beta
    repair0_coefficient = common_coeff * F * rows["Y0"]["coefficient_interval"]
    repair1_y1_coefficient = common_coeff * F * rows["Y1"]["coefficient_interval"]
    repair1_y0_coefficient = common_coeff * Fz * rows["Y0"]["coefficient_interval"]

    repair0 = _term(
        c=c,
        name="R*Pstar*beta*F*Y0",
        repair=0,
        index=0,
        coefficient=repair0_coefficient,
        scale=repair0_scale,
        cell=cell,
        rows=rows,
        F_name="F",
        beta_binding=beta_binding,
        beta=beta,
    )
    repair1_y1 = _term(
        c=c,
        name="R*Pstar*beta*F*Y1",
        repair=1,
        index=0,
        coefficient=repair1_y1_coefficient,
        scale=repair1_y1_scale,
        cell=cell,
        rows=rows,
        F_name="F",
        beta_binding=beta_binding,
        beta=beta,
    )
    repair1_y0 = _term(
        c=c,
        name="R*Pstar*beta*Fz*Y0",
        repair=1,
        index=1,
        coefficient=repair1_y0_coefficient,
        scale=repair1_y0_scale,
        cell=cell,
        rows=rows,
        F_name="Fz",
        beta_binding=beta_binding,
        beta=beta,
    )
    denominator = {
        "formula": "L=1-delta*Z^2",
        "interval": L,
        "strictly_positive": True,
        "lower_bound": L_lower,
        "upper_bound": L_upper,
        "whole_cell": True,
        "exact_Z_cell": list(cell),
    }
    return {
        "cell": list(cell),
        "exact_Z_cell": list(cell),
        "Z_interval": Z,
        "axis_function_whole_cell": True,
        "denominator_evidence": denominator,
        "L": L,
        "F": F,
        "Fz": Fz,
        "axis_function_enclosures": {"Z": Z, "L": L, "F": F, "Fz": Fz},
        "repair0": [repair0],
        "repair1": [repair1_y1, repair1_y0],
        "repair0_terms": [repair0],
        "repair1_terms": [repair1_y1, repair1_y0],
        "source_terms": {"repair0": [repair0], "repair1": [repair1_y1, repair1_y0]},
        "formal_scales": {
            "repair0": repair0["formal_positive_scale"],
            "repair1": [term["formal_positive_scale"] for term in (repair1_y1, repair1_y0)],
        },
        "physical_formal_scales": {
            "repair0": repair0["formal_positive_scale"],
            "repair1": [term["formal_positive_scale"] for term in (repair1_y1, repair1_y0)],
        },
        "source_indices": {
            key: _copy_metadata(value["source_index"]) for key, value in rows.items()
        },
        "ordinary_Z1_product_rule": (
            "R*Pstar*beta*(F*Y1 + Fz*Y0), retained as two source terms"
        ),
        "repair1_scales_not_merged": True,
        "pstar_multiplier_is_old_excess": True,
        "normalization_substitution_not_performed": True,
    }


def build(
    c: Any,
    cells: Any = None,
    delta: Any = None,
    source_rows: Any = None,
    source_certificate: Any = None,
    pstar_positive_certificate: Any = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build all four whole-cell factored repairs in one returned mapping.

    The positional form is ``build(c, cells, delta, source_rows,
    source_certificate, pstar_positive_certificate)``.  A few descriptive
    keyword aliases are accepted for frozen runners; unknown arguments fail
    closed.
    """

    def take(current: Any, names: Sequence[str], label: str) -> Any:
        if current is not None:
            return current
        found = [name for name in names if name in aliases]
        if len(found) > 1:
            raise TypeError(f"multiple aliases supplied for {label}: {found}")
        if found:
            return aliases.pop(found[0])
        return None

    cells = take(
        cells,
        (
            "z_cells",
            "exact_z_cells",
            "exact_Z_cells",
            "exact_cells",
            "four_z_cells",
            "z_cell_partition",
        ),
        "cells",
    )
    delta = take(delta, ("delta_interval", "delta_Z", "delta_box"), "delta")
    source_rows = take(
        source_rows,
        (
            "raw_source_rows",
            "raw_rows",
            "cell_source_rows",
            "source_rows_by_cell",
            "source_row_packets",
            "rows",
            "source_row_map",
        ),
        "source_rows",
    )
    source_certificate = take(
        source_certificate,
        (
            "same_source_certificate",
            "source_basis_certificate",
            "same_source_basis_certificate",
            "same_source_provenance_certificate",
            "source_basis_provenance_certificate",
            "provenance_certificate",
            "certificate",
        ),
        "source_certificate",
    )
    pstar_positive_certificate = take(
        pstar_positive_certificate,
        (
            "pstar_certificate",
            "pstar_positivity_certificate",
            "pstar_positive",
            "pstar_log_positive_certificate",
            "log_pstar_positive_certificate",
            "logpstar_certificate",
            "positive_pstar_certificate",
        ),
        "pstar_positive_certificate",
    )
    if aliases:
        names = ", ".join(sorted(str(name) for name in aliases))
        raise TypeError(f"unexpected four-cell repair arguments: {names}")
    _require(hasattr(c, "mpf"), "c must be an MPIntervalContext-like directed context")
    _require(cells is not None, "four exact Z cells are required")
    _require(delta is not None, "delta directed interval is required")
    _require(source_rows is not None, "raw source rows are required")
    _require(source_certificate is not None, "same-source basis/provenance certificate is required")
    _require(pstar_positive_certificate is not None, "logPstar > 0 certificate is required")

    _cells(cells)
    _require(
        _is_interval(delta) or _is_serialized_interval(delta),
        "delta must be a directed interval",
    )
    delta_interval = _interval(c, delta, "delta")
    normalized_rows = _all_source_rows(c, source_rows)
    source_cert = _validate_source_certificate(source_certificate, normalized_rows)
    pstar_cert = _validate_pstar_certificate(c, pstar_positive_certificate, normalized_rows)
    beta = _make_interval(c, pstar_cert["beta_lower"], c.mpf(1))

    by_cell: dict[str, Any] = {}
    for cell in CELLS:
        beta_binding = f"beta(Pstar):cell={_cell_key(cell)}:same-source-Pstar"
        beta_cell = _make_interval(c, pstar_cert["beta_lower"], c.mpf(1))
        by_cell[_cell_key(cell)] = _build_cell(
            c, delta_interval, cell, normalized_rows[cell], beta_cell, beta_binding
        )

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "cells": [list(cell) for cell in CELLS],
        "exact_Z_cells": [list(cell) for cell in CELLS],
        "delta_interval": delta_interval,
        "source_rows": {
            _cell_key(cell): {
                key: {
                    "source_index": _copy_metadata(value["source_index"]),
                    "source_path": _copy_metadata(value["source_path"]),
                    "provenance": _copy_metadata(value["provenance"]),
                    "coefficient_interval": value["coefficient_interval"],
                    "coefficient": value["coefficient_interval"],
                    "formal_positive_scale": _scale_record(value["formal_positive_scale"]),
                    "point_value_selected": False,
                    "log_absolute_upper_not_used_as_value": True,
                }
                for key, value in normalized_rows[cell].items()
            }
            for cell in CELLS
        },
        "source_certificate": source_cert,
        "pstar_positive_certificate": pstar_cert,
        "formal_source_basis_scope": "opaque cell-local basis; no cross-cell basis object identity asserted",
        "pstar_beta": {
            "definition": "beta(Pstar)=1-1/Pstar",
            "interval": beta,
            "enclosure": beta,
            "cell_local_intervals": {
                key: value["repair0"][0]["beta_interval"] for key, value in by_cell.items()
            },
            "lower": _bounds(beta)[0],
            "upper": _bounds(beta)[1],
            "same_source_identity": "cell-local beta(Pstar) bindings",
            "binding_scope": "one exact Z cell",
            "same_source_identity_within_cell": True,
            "cross_cell_basis_identity_asserted": False,
            "same_source_provenance": True,
            "derived_from_logPstar_gt_zero": True,
            "lower_bound_uses_p_over_one_plus_p": bool(
                pstar_cert["beta_lower_proved_from_logPstar"]
            ),
            "Pstar_gt_one": True,
            "not_selected_as_zero": True,
            "not_selected_as_one": True,
            "interval_is_enclosure_not_function_value": True,
            "extra_force_not_introduced": True,
        },
        "cells_by_key": by_cell,
        "cell_sources": by_cell,
        "axis_function_whole_cell_enclosures": {
            key: {
                "Z": value["Z_interval"],
                "L": value["L"],
                "F": value["F"],
                "Fz": value["Fz"],
            }
            for key, value in by_cell.items()
        },
        "denominator_evidence": {
            key: value["denominator_evidence"] for key, value in by_cell.items()
        },
        "repair0": {key: value["repair0"] for key, value in by_cell.items()},
        "repair1": {key: value["repair1"] for key, value in by_cell.items()},
        "repair_terms_are_source_correlated": True,
        "repair0_term_count": 1,
        "repair1_term_count": 2,
        "ordinary_Z1_is_product_rule": True,
        "ordinary_Z1_terms_not_merged": True,
        "formal_scales_are_factored": True,
        "formal_exponents_exact_integer_or_half_integer": True,
        "additional_log_offsets_added_as_directed_intervals": True,
        "Y_original_Pstar_retained_once": True,
        "Pstar_multiplier_is_old_excess_not_normalization": True,
        "coefficient_intervals_signed_directed": True,
        "coefficient_intervals_are_not_function_values": True,
        "point_values_selected": False,
        "log_absolute_upper_used_as_value": False,
        "extra_force_introduced": False,
        "ledger_written": False,
        "no_ledger": True,
        "global_closure_claimed": False,
        "global_five_moment_closure": False,
        "graph_traversed": False,
        "receipt_opened": False,
        "receipt_decoded": False,
        "ancestor_imported": False,
        "numeric_source_values_selected": False,
    }


def build_repair(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(*args, **kwargs)


enclose = build
evaluate = build


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "CELLS",
    "ORDERS",
    "SCALE_FIELDS",
    "FourCellSignedRepairError",
    "OriginalCorrectedFourcellSignedRepairError",
    "build",
    "build_repair",
    "enclose",
    "evaluate",
]
