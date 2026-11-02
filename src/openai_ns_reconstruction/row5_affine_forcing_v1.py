"""Source-bound affine row-5 forcing and its fixed-:math:`x` ``Z1`` jet.

The raw row-5 relation used here is

``W1_x = W5`` and
``L0 W1 = g5 + C0_51 W1 + C0_52 W2 + C0_53 W3 + C0_55 W5``
``                 + B1_51 W1_Z + B1_52 W2_Z + B1_53 W3_Z``.

The ``Kchi`` right hand side is that row plus ``2*chi*W1``.  This module is
only a source-bound expression builder.  It does not solve for ``W``, assign
missing state entries the value zero, evaluate a physical ``h5`` field, or
perform a Navier--Stokes validation.  In particular, all source sectors stay
as separate additive records.  The four formal powers
``(hb, F0_anchor, Pstar_squared, Lambda)`` are copied from the source term;
they are never evaluated while building the expression.

The ordinary derivative is at fixed ``x``.  Therefore the output keeps the
four product-rule pieces for every matrix entry and the two pieces from
``2*chi*W1``::

    C0_Z1*W + C0_Z0*W_Z
    B1_Z1*W_Z + B1_Z0*W_ZZ
    2*chi_Z*W1 + 2*chi*W1_Z

``evaluate`` is deliberately a separate caller-owned operation.  It accepts
complete state jets and (when needed by nonzero formal powers) the four
formal atom intervals, multiplies interval factors with directed arithmetic,
and only then forms a total enclosure.  A missing state or formal atom raises
instead of becoming a silent zero.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "row5_affine_forcing_v1"
COMPONENT = "source_bound_row5_affine_forcing"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
STATE_COMPONENTS = (1, 2, 3, 5)
STATE_Z_ORDERS = (0, 1, 2)
C0_KEYS = ("5,1", "5,2", "5,3", "5,5")
B1_KEYS = ("5,1", "5,2", "5,3")
MATRIX_KEYS = {"C0": C0_KEYS, "B1": B1_KEYS}

_MISSING = object()
_COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "raw_interval",
    "directed_interval",
    "interval",
    "enclosure",
    "box",
    "value",
)
_POWER_KEYS = ("powers", "formal_powers", "source_powers")


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS",
    "STATE_COMPONENTS",
    "STATE_Z_ORDERS",
    "C0_KEYS",
    "B1_KEYS",
    "build",
    "build_row5_affine_forcing",
    "row5_affine_forcing",
    "enclose",
    "evaluate",
    "evaluate_expression",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None or value is _MISSING:
        return []
    return list(value) if _is_sequence(value) else [value]


def _clone(value: Any) -> Any:
    """Clone ordinary receipt-shaped metadata without importing a producer."""

    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None:
            return nested[1 if upper else 0]
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return nested
        return endpoint
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _validate_context(c: Any) -> None:
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide mpf and make_mpf directed interval operations")


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label} must be a four-entry MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} must contain integers only")
    sign, mantissa, exponent, bitcount = (int(item) for item in value)
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(f"{label} has invalid MPF tuple fields")
    if mantissa == 0 and bitcount != 0:
        raise ValueError(f"{label} is a non-finite MPF tuple")
    return sign, mantissa, exponent, bitcount


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode a directed interval without midpoint or float conversion."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        result = _make_interval(c, value, value)
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_tuple(value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_tuple(value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            result = c.make_mpf((lower, upper))
        else:
            nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
            if nested is not _MISSING and nested is not value:
                result = _interval(c, nested, label + ".coefficient")
            else:
                lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
                upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
                if lower is _MISSING or upper is _MISSING:
                    raise TypeError(f"{label} is not interval-compatible")
                result = _make_interval(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make_interval(c, value[0], value[1])
    else:
        try:
            point = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
        result = point if getattr(point, "_mpi_", None) is not None else _make_interval(c, point, point)

    raw = getattr(result, "_mpi_", None)
    if raw is None:
        raise TypeError(f"{label} did not produce a directed interval")
    lower = c.make_mpf((raw[0], raw[0]))
    upper = c.make_mpf((raw[1], raw[1]))
    if upper < lower:
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _power_values(raw: Mapping[Any, Any], label: str) -> list[int]:
    powers: Any = _MISSING
    for key in _POWER_KEYS:
        if key in raw:
            powers = raw[key]
            break
    if powers is _MISSING:
        nested = raw.get("origin")
        if isinstance(nested, Mapping):
            for key in _POWER_KEYS:
                if key in nested:
                    powers = nested[key]
                    break
    if powers is _MISSING:
        powers = (
            raw.get("hb_power", _MISSING),
            raw.get("F0_power", raw.get("F0_anchor_power", _MISSING)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", _MISSING)),
            raw.get("Lambda_power", raw.get("lambda_power", _MISSING)),
        )
        if any(item is _MISSING for item in powers):
            raise KeyError(f"{label}.powers is required; formal powers may not be inferred")
    if isinstance(powers, Mapping):
        powers = [powers.get(atom, _MISSING) for atom in SOURCE_POWER_BASIS]
    if not _is_sequence(powers) or len(powers) != 4:
        raise TypeError(f"{label}.powers must contain four integer powers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise TypeError(f"{label}.powers must contain four integer powers")
    return [int(item) for item in powers]


def _metadata_values(raw: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(_items(raw[name]))
    for container_name in ("origin", "metadata", "original_term"):
        container = raw.get(container_name)
        if isinstance(container, Mapping):
            for name in names:
                if name in container:
                    values.extend(_items(container[name]))
    return values


def _paths(raw: Mapping[Any, Any]) -> list[Any]:
    paths: list[Any] = []
    # A singular source_path may itself be a list of path segments.  Retain
    # that list as one opaque path instead of flattening it.
    if "source_path" in raw:
        paths.append(_clone(raw["source_path"]))
    elif "path" in raw:
        paths.append(_clone(raw["path"]))
    for key in ("source_paths", "paths"):
        if key in raw:
            paths.extend(_clone(item) for item in _items(raw[key]))
    for container_name in ("origin", "metadata", "original_term"):
        container = raw.get(container_name)
        if not isinstance(container, Mapping):
            continue
        if "source_path" in container:
            paths.append(_clone(container["source_path"]))
        elif "path" in container:
            paths.append(_clone(container["path"]))
        for key in ("source_paths", "paths"):
            if key in container:
                paths.extend(_clone(item) for item in _items(container[key]))
    return paths


def _bool_metadata(raw: Mapping[Any, Any], name: str, default: bool = False) -> bool:
    if name in raw:
        value = raw[name]
    elif name == "uncertainty" and "is_uncertainty" in raw:
        value = raw["is_uncertainty"]
    else:
        value = _MISSING
    if value is _MISSING:
        for container_name in ("origin", "metadata", "original_term"):
            container = raw.get(container_name)
            if isinstance(container, Mapping) and name in container:
                value = container[name]
                break
    if value is _MISSING:
        return bool(default)
    if not isinstance(value, bool):
        raise TypeError(f"term.{name} must be boolean")
    return value


def _normalise_term(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    coefficient = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), _MISSING)
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    coefficient_interval = _interval(c, coefficient, label + ".coefficient")
    powers = _power_values(raw, label)
    source_ids = _metadata_values(raw, ("source_ids", "source_id", "sources", "source", "input_id"))
    if not source_ids:
        raise ValueError(f"{label} must retain at least one source_id")
    source_paths = _paths(raw)
    error_classes = _metadata_values(raw, ("error_classes", "error_class", "uncertainty_class"))
    source_error_classes = _metadata_values(raw, ("source_error_classes", "source_error_class"))
    if not source_error_classes:
        source_error_classes = list(error_classes)
    uncertainty = _bool_metadata(raw, "uncertainty", False)
    error_only = _bool_metadata(raw, "error_only", False)
    symmetric_error = _bool_metadata(raw, "symmetric_error", error_only)
    signed_field = _bool_metadata(raw, "signed_field", False)
    if error_only and not symmetric_error:
        raise ValueError(f"{label}.error_only requires symmetric_error")
    if error_only and signed_field:
        raise ValueError(f"{label}.error_only may not be signed_field")
    if error_only:
        lower = _interval(c, coefficient_interval, label + ".coefficient")._mpi_[0]
        upper = _interval(c, coefficient_interval, label + ".coefficient")._mpi_[1]
        lower_value = c.make_mpf((lower, lower))
        upper_value = c.make_mpf((upper, upper))
        if lower_value != -upper_value:
            raise ValueError(f"{label}.error_only coefficient must be symmetric")
    original = _clone(dict(raw))
    source_path = _clone(raw.get("source_path", raw.get("path")))
    result = _clone(dict(raw))
    result.update(
        {
            "coefficient": coefficient_interval,
            "coefficient_interval": coefficient_interval,
            "powers": list(powers),
            "formal_powers": list(powers),
            "power_basis": list(SOURCE_POWER_BASIS),
            "source_ids": list(source_ids),
            "source_paths": source_paths,
            "source_path": source_path,
            "error_classes": list(error_classes),
            "source_error_classes": list(source_error_classes),
            "uncertainty": uncertainty,
            "error_only": error_only,
            "symmetric_error": symmetric_error,
            "signed_field": signed_field,
            "original_term": original,
            "original_source_term": _clone(original),
            "midpoint_values_used": False,
        }
    )
    return result


def _order_key(key: Any) -> int | None:
    if key == 0 or (isinstance(key, str) and key.lower() in {"0", "z0", "ordinary_z0", "value"}):
        return 0
    if key == 1 or (isinstance(key, str) and key.lower() in {"1", "z1", "ordinary_z1", "derivative"}):
        return 1
    return None


def _term_list(c: Any, raw: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and any(key in raw for key in _COEFFICIENT_KEYS + _POWER_KEYS):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a term sequence")
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_normalise_term(c, item, f"{label}[{index}]") for index, item in enumerate(values)]


def _unwrap_rows(raw: Any, wrappers: Sequence[str], label: str) -> Mapping[Any, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a mapping")
    current = raw
    for wrapper in wrappers:
        nested = current.get(wrapper, _MISSING)
        if isinstance(nested, Mapping) and nested is not current:
            current = nested
            break
    return current


def _normalise_g5(c: Any, raw: Any) -> dict[str, list[dict[str, Any]]]:
    current = _unwrap_rows(raw, ("g5", "g", "ordinary_Z", "source_Z", "rows", "output"), "source_packet.g5")
    # Actual saved output may be g[order]["5"].
    if set(current) and all(_order_key(key) is not None for key in current):
        found = {}
        for key, value in current.items():
            order = _order_key(key)
            if isinstance(value, Mapping) and "5" in value:
                value = value["5"]
            elif isinstance(value, Mapping) and 5 in value:
                value = value[5]
            found[f"Z{order}"] = value
    else:
        found = {}
        for key, value in current.items():
            order = _order_key(key)
            if order is not None:
                found[f"Z{order}"] = value
            elif str(key) in {"0", "1"} and isinstance(value, Mapping) and "5" in value:
                found[f"Z{key}"] = value["5"]
    # When the outer map is g["0"]["5"], the first branch above does not
    # see Z rows because the values are row dictionaries; recover them here.
    if set(found) != set(ORDINARY_Z_LABELS):
        found = {}
        for key, value in current.items():
            order = _order_key(key)
            if order is not None:
                if isinstance(value, Mapping) and "5" in value:
                    value = value["5"]
                found[f"Z{order}"] = value
    if set(found) != set(ORDINARY_Z_LABELS):
        raise ValueError("source_packet.g5 requires explicit Z0 and Z1 rows")
    return {label: _term_list(c, found[label], "source_packet.g5." + label) for label in ORDINARY_Z_LABELS}


def _matrix_section(raw: Any, name: str) -> Mapping[Any, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"source_packet.{name} must be a mapping")
    current = raw
    for wrapper in (name, "ordinary_Z", "source_Z", "rows", "jets", "terms"):
        nested = current.get(wrapper, _MISSING)
        if isinstance(nested, Mapping) and nested is not current:
            current = nested
            break
    return current


def _normalise_matrix(c: Any, raw: Any, name: str) -> dict[str, dict[str, list[dict[str, Any]]]]:
    current = _matrix_section(raw, name)
    rows: dict[str, Any] = {}
    for key, value in current.items():
        order = _order_key(key)
        if order is not None:
            if f"Z{order}" in rows:
                raise ValueError(f"source_packet.{name} supplies Z{order} twice")
            rows[f"Z{order}"] = value
        elif isinstance(key, str) and key.lower() in {f"{name.lower()}_z0", f"{name.lower()}_z1"}:
            rows[f"Z{key[-1]}"] = value
    if set(rows) != set(ORDINARY_Z_LABELS):
        raise ValueError(f"source_packet.{name} requires explicit Z0 and Z1 rows")
    result: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for label in ORDINARY_Z_LABELS:
        row = rows[label]
        if not isinstance(row, Mapping):
            raise TypeError(f"source_packet.{name}.{label} must be a matrix row mapping")
        by_key: dict[str, Any] = {}
        for key, value in row.items():
            canonical = f"{int(key[0])},{int(key[1])}" if isinstance(key, (tuple, list)) and len(key) == 2 else str(key)
            by_key[canonical] = value
        required = set(MATRIX_KEYS[name])
        if not required.issubset(by_key):
            missing = sorted(required - set(by_key))
            raise ValueError(f"source_packet.{name}.{label} missing matrix entries {missing}")
        result[label] = {
            key: _term_list(c, by_key[key], f"source_packet.{name}.{label}.{key}")
            for key in MATRIX_KEYS[name]
        }
    return result


def _normalise_axes(c: Any, raw: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping):
        current = raw
        # A direct term wins over wrappers named ``chi`` or ``axis``.
        if not any(key in current for key in _COEFFICIENT_KEYS + _POWER_KEYS):
            for wrapper in (label, "axis", "terms", "rows", "value"):
                nested = current.get(wrapper, _MISSING)
                if nested is not _MISSING and nested is not current:
                    current = nested
                    break
        values = _term_list(c, current, label)
    elif _is_sequence(raw):
        values = _term_list(c, raw, label)
    else:
        raise TypeError(f"{label} must be an axis term or term sequence")
    result: list[dict[str, Any]] = []
    for index, term in enumerate(values):
        if term["powers"] != [0, 0, 0, 0]:
            raise ValueError(f"{label}[{index}] formal powers must be zero")
        if term["error_only"]:
            raise ValueError(f"{label}[{index}] may not be error-only")
        # An explicit false certificate must fail closed; absence is a
        # caller condition and is not silently promoted to proof.
        raw_term = term["original_term"]
        if isinstance(raw_term, Mapping):
            for key in ("whole_Z", "whole_z", "certified_whole_Z", "whole_Z_certified", "Z_independent"):
                if raw_term.get(key) is False:
                    raise ValueError(f"{label}[{index}] does not certify whole-Z scope")
        term["axis_whole_Z_scope"] = True
        term["axis_name"] = label
        result.append(term)
    return result


def _packet_sections(source_packet: Any) -> tuple[Any, Any, Any, Any]:
    if not isinstance(source_packet, Mapping):
        raise TypeError("source_packet must be a mapping")
    current = source_packet
    for wrapper in ("source_packet", "output", "result", "forcing"):
        nested = current.get(wrapper, _MISSING)
        if isinstance(nested, Mapping) and nested is not current:
            current = nested
            break
    g5 = current.get("g5", current.get("g", _MISSING))
    matrix = current.get("matrix", current)
    if g5 is _MISSING:
        raise KeyError("source_packet.g5 is required")
    if not isinstance(matrix, Mapping):
        raise TypeError("source_packet.matrix must be a mapping")
    c0 = matrix.get("C0", _MISSING)
    b1 = matrix.get("B1", _MISSING)
    if c0 is _MISSING or b1 is _MISSING:
        raise KeyError("source_packet.C0 and source_packet.B1 are required")
    return g5, c0, b1, current


def _state_reference(component: int, order: int) -> dict[str, int]:
    return {"component": int(component), "ordinary_Z_order": int(order)}


def _emit_term(
    source: Mapping[str, Any],
    *,
    output_label: str,
    source_family: str,
    source_key: str,
    source_order: int,
    scalar: int,
    state_reference: Mapping[str, int] | None,
    role: str,
    matrix_name: str | None = None,
) -> dict[str, Any]:
    original = _clone(source["original_term"])
    result = {
        "coefficient": source["coefficient"],
        "coefficient_interval": source["coefficient_interval"],
        "powers": list(source["powers"]),
        "formal_powers": list(source["formal_powers"]),
        "power_basis": list(SOURCE_POWER_BASIS),
        "scalar": int(scalar),
        "state_reference": _clone(state_reference),
        "state_ref": _clone(state_reference),
        "state": _clone(state_reference),
        "reference": _clone(state_reference),
        "source_only": state_reference is None,
        "source_only_g5": state_reference is None and source_family == "g5",
        "output": output_label,
        "ordinary_Z_order": 0 if output_label == "h5" else 1,
        "source_ordinary_Z_order": int(source_order),
        "source_family": source_family,
        "source_key": source_key,
        "matrix_name": matrix_name,
        "matrix_key": source_key if matrix_name is not None else None,
        "role": role,
        "term_role": role,
        "source_ids": list(source["source_ids"]),
        "source_paths": _clone(source["source_paths"]),
        "source_path": _clone(source.get("source_path")),
        "uncertainty": bool(source["uncertainty"]),
        "error_only": bool(source["error_only"]),
        "symmetric_error": bool(source["symmetric_error"]),
        "signed_field": bool(source["signed_field"]),
        "error_classes": list(source["error_classes"]),
        "source_error_classes": list(source["source_error_classes"]),
        "original_term": original,
        "original_source_term": _clone(original),
        "source_term": _clone(original),
        "midpoint_values_used": False,
        "formal_powers_preserved": True,
        "source_provenance_preserved": True,
    }
    if state_reference is not None:
        result["state_component"] = int(state_reference["component"])
        result["state_ordinary_Z_order"] = int(state_reference["ordinary_Z_order"])
    else:
        result["state_component"] = None
        result["state_ordinary_Z_order"] = None
    return result


def _build_rows(
    g5: Mapping[str, Sequence[Mapping[str, Any]]],
    c0: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    b1: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    chi: Sequence[Mapping[str, Any]],
    chi_z: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    h5: list[dict[str, Any]] = []
    h5_z1: list[dict[str, Any]] = []

    for source in g5["Z0"]:
        h5.append(_emit_term(source, output_label="h5", source_family="g5", source_key="g5", source_order=0,
                              scalar=1, state_reference=None, role="source_only_g5"))
    for matrix_name, entries, state_order in (("C0", c0, 0), ("B1", b1, 1)):
        for key in MATRIX_KEYS[matrix_name]:
            component = int(key.split(",")[1])
            for source in entries["Z0"][key]:
                h5.append(_emit_term(source, output_label="h5", source_family=matrix_name, source_key=key,
                                      source_order=0, scalar=1,
                                      state_reference=_state_reference(component, state_order),
                                      role=f"{matrix_name}_Z0_times_state" , matrix_name=matrix_name))
    for source in chi:
        h5.append(_emit_term(source, output_label="h5", source_family="axis", source_key="chi", source_order=0,
                             scalar=2, state_reference=_state_reference(1, 0), role="two_chi_times_W1"))

    for source in g5["Z1"]:
        h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="g5", source_key="g5", source_order=1,
                                 scalar=1, state_reference=None, role="source_only_g5_Z1"))

    # C0_Z1*W_j + C0_Z0*W_j_Z.
    for key in C0_KEYS:
        component = int(key.split(",")[1])
        for source in c0["Z1"][key]:
            h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="C0", source_key=key,
                                     source_order=1, scalar=1, state_reference=_state_reference(component, 0),
                                     role="C0_Z1_times_state", matrix_name="C0"))
        for source in c0["Z0"][key]:
            h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="C0", source_key=key,
                                     source_order=0, scalar=1, state_reference=_state_reference(component, 1),
                                     role="C0_Z0_times_state_Z1", matrix_name="C0"))

    # B1_Z1*W_j_Z + B1_Z0*W_j_ZZ.
    for key in B1_KEYS:
        component = int(key.split(",")[1])
        for source in b1["Z1"][key]:
            h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="B1", source_key=key,
                                     source_order=1, scalar=1, state_reference=_state_reference(component, 1),
                                     role="B1_Z1_times_state_Z1", matrix_name="B1"))
        for source in b1["Z0"][key]:
            h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="B1", source_key=key,
                                     source_order=0, scalar=1, state_reference=_state_reference(component, 2),
                                     role="B1_Z0_times_state_Z2", matrix_name="B1"))

    for source in chi_z:
        h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="axis", source_key="chi_Z", source_order=1,
                                scalar=2, state_reference=_state_reference(1, 0), role="two_chi_Z_times_W1"))
    for source in chi:
        h5_z1.append(_emit_term(source, output_label="h5_Z1", source_family="axis", source_key="chi", source_order=0,
                                scalar=2, state_reference=_state_reference(1, 1), role="two_chi_times_W1_Z"))
    return h5, h5_z1


def build(c: Any, source_packet: Any, chi: Any, chi_z: Any) -> dict[str, Any]:
    """Build the source-bound row-5 ``h5`` and fixed-``x`` ``h5_Z1`` rows.

    ``source_packet`` may be the explicit shape described in the module
    docstring or a saved-output-shaped mapping where ``g[Z]["5"]`` and
    ``matrix[C0|B1][Z]["5,j"]`` are nested under ``output``/``result``.
    Every required row and matrix entry is required.  No interval arithmetic
    is performed between source terms here.
    """

    _validate_context(c)
    g5_raw, c0_raw, b1_raw, packet_current = _packet_sections(source_packet)
    g5 = _normalise_g5(c, g5_raw)
    c0 = _normalise_matrix(c, c0_raw, "C0")
    b1 = _normalise_matrix(c, b1_raw, "B1")
    chi_rows = _normalise_axes(c, chi, "chi")
    chi_z_rows = _normalise_axes(c, chi_z, "chi_z")
    h5, h5_z1 = _build_rows(g5, c0, b1, chi_rows, chi_z_rows)
    normal_packet = {"g5": g5, "C0": c0, "B1": b1}
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_source_bound_row5_affine_forcing",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "state_components": list(STATE_COMPONENTS),
        "state_Z_orders": list(STATE_Z_ORDERS),
        "source_packet": normal_packet,
        "g5": g5,
        "matrix": {"C0": c0, "B1": b1},
        "C0": c0,
        "B1": b1,
        "chi": chi_rows,
        "chi_z": chi_z_rows,
        "chi_Z": chi_z_rows,
        "axes": {"chi": chi_rows, "chi_Z": chi_z_rows},
        "h5": h5,
        "h5_Z1": h5_z1,
        "h5_terms": h5,
        "h5_Z1_terms": h5_z1,
        "rows": {"h5": h5, "h5_Z1": h5_z1},
        "terms": {"h5": h5, "h5_Z1": h5_z1},
        "additive_terms": {"h5": h5, "h5_Z1": h5_z1},
        "formula": (
            "h5 = g5 + C0_51*W1 + C0_52*W2 + C0_53*W3 + C0_55*W5 "
            "+ B1_51*W1_Z + B1_52*W2_Z + B1_53*W3_Z + 2*chi*W1"
        ),
        "h5_formula": (
            "g5 + sum_j(C0_5j*Wj) + sum_j(B1_5j*Wj_Z) + 2*chi*W1, "
            "j in {1,2,3,5} for C0 and j in {1,2,3} for B1"
        ),
        "h5_Z1_formula": (
            "g5_Z1 + sum_j(C0_5j_Z1*Wj + C0_5j_Z0*Wj_Z) "
            "+ sum_j(B1_5j_Z1*Wj_Z + B1_5j_Z0*Wj_ZZ) "
            "+ 2*chi_Z*W1 + 2*chi*W1_Z"
        ),
        "row5_identity": "W1_x = W5",
        "ordinary_Z_semantics": "fixed x; ordinary Z1 product rule",
        "fixed_x_ordinary_Z": True,
        "fixed_x_is_caller_semantics": True,
        "rho_equals_x_squared": True,
        "rho_equals_x_squared_is_premise": True,
        "rho_equals_x_squared_proved": False,
        "geometry_Z_independent": True,
        "geometry_Z_independent_is_premise": True,
        "geometry_Z_independent_proved": False,
        "physical_dx_applied": False,
        "rho_x2_chain_rule_applied": False,
        "source_terms_individual": True,
        "interval_sum_in_build": False,
        "source_duplicates_preserved": True,
        "source_ids_preserved": True,
        "source_paths_preserved": True,
        "formal_powers_preserved": True,
        "formal_powers_materialized": False,
        "uncertainty_preserved": True,
        "error_classes_preserved_without_mixing": True,
        "state_references_explicit": True,
        "missing_state_is_not_zero": True,
        "axis_same_owner_required": True,
        "axis_same_owner_proved": False,
        "axis_same_owner": False,
        "g5_matrix_same_receipt_path_assumed": True,
        "g5_matrix_same_input_jets_assumed": True,
        "g5_matrix_same_receipt_and_input_jets_is_premise": True,
        "g5_matrix_same_receipt_proved": False,
        # A runner may attach and authenticate the historical source receipt,
        # but this direct adapter receives only a packet and does not verify
        # that external evidence.  Keep these proof flags false.
        "same_receipt_matrix_and_g": False,
        "same_production_path": False,
        "same_receipt_matrix_and_g_condition": "caller/runner receipt evidence required; not verified by this leaf",
        "same_production_path_condition": "caller/runner receipt evidence required; not verified by this leaf",
        "runner_receipt_verification_required": True,
        "same_W_owner": False,
        "c0_row5_raw_B0_interpretation_is_premise": True,
        "c0_row5_raw_B0_owner_proved": False,
        "c0_kinematic_subtraction_rows": ["(1,5)", "(2,6)", "(3,6)"],
        "same_state_owner_proved": False,
        "source_owner_proved": False,
        "same_owner_proved": False,
        "same_source_identity_proved": False,
        "source_bound_only": True,
        "physical_h5_materialized": False,
        "physical_field_materialized": False,
        "h5_is_physical_function": False,
        "complete_W1": False,
        "state_solved": False,
        "closed_state": False,
        "full_state": False,
        "global_NS": False,
        "m15": False,
        "ns_validation": False,
        "scope": {
            "source_bound_affine_expression": True,
            "raw_g5_consumed": True,
            "raw_C0_B1_consumed": True,
            "axis_chi_and_chi_Z_consumed": True,
            "axis_scope_same_owner_is_condition": True,
            "fixed_x": True,
            "rho_equals_x_squared_is_premise": True,
            "geometry_Z_independent_is_premise": True,
            "C0_row5_raw_B0_is_source_premise": True,
            "kinematic_subtraction_rows": ["(1,5)", "(2,6)", "(3,6)"],
            "state_W_consumed": False,
            "state_W_owner_proved": False,
            "state_materialized": False,
            "physical_h5": False,
            "physical_dx": False,
            "complete_collar": False,
            "complete_six_field": False,
            "m15": False,
            "ns_validation": False,
            "receipt_reads": False,
            "ancestor_imports": False,
            "ledger_written": False,
            "midpoint_values_used": False,
        },
        "conditions": {
            "source_packet_has_explicit_g5_Z0_Z1": True,
            "source_packet_has_explicit_C0_B1_row5_entries": True,
            "formal_power_basis": list(SOURCE_POWER_BASIS),
            "axis_chi_whole_Z_scope_required": True,
            "axis_chi_Z_whole_Z_scope_required": True,
            "axis_same_owner": "caller condition; not proved by this leaf",
            "W1_x_equals_W5": "raw-system identity supplied by caller/source",
            "fixed_x_ordinary_Z": True,
            "rho_equals_x_squared": "caller premise; no coordinate conversion performed here",
            "geometry_Z_independent": "caller premise; not proved here",
            "g5_matrix_same_receipt_and_input_jets": "caller/runner premise; not verified by this direct leaf",
            "C0_row5_equals_raw_B0": "source-code premise; only raw B0 kinematic subtractions at (1,5),(2,6),(3,6)",
            "missing_W_or_W_Z_or_W_ZZ": "reject; never infer zero",
        },
        "resource_limits": {
            "input_term_count": sum(
                len(g5[label]) for label in ORDINARY_Z_LABELS
            ) + sum(
                len(section[label][key])
                for section in (c0, b1)
                for label in ORDINARY_Z_LABELS
                for key in MATRIX_KEYS["C0" if section is c0 else "B1"]
            ) + len(chi_rows) + len(chi_z_rows),
            "h5_term_count": len(h5),
            "h5_Z1_term_count": len(h5_z1),
        },
        "caller_packet_wrapper_seen": packet_current is not source_packet,
    }


def _formal_atoms(state: Mapping[str, Any]) -> Mapping[Any, Any]:
    for key in ("formal_atoms", "atom_values", "formal_atom_values", "atoms"):
        value = state.get(key, _MISSING)
        if isinstance(value, Mapping):
            return value
    return state


def _atom_value(c: Any, atoms: Mapping[Any, Any], name: str) -> Any:
    aliases = {
        "hb": ("hb", "H_b", "h_b"),
        "F0_anchor": ("F0_anchor", "F0", "F0Anchor", "F0_anchor_value"),
        "Pstar_squared": ("Pstar_squared", "Pstar^2", "Pstar2", "Pstar"),
        "Lambda": ("Lambda", "lambda", "Λ", "Lambda_value"),
    }
    for key in aliases[name]:
        if key in atoms:
            return _interval(c, atoms[key], "formal atom " + name)
    raise KeyError(f"caller state is missing formal atom {name!r}")


def _formal_factor(c: Any, atoms: Mapping[Any, Any], powers: Sequence[int]) -> Any:
    value = c.mpf(1)
    for name, power in zip(SOURCE_POWER_BASIS, powers):
        power = int(power)
        if power == 0:
            continue
        atom = _atom_value(c, atoms, name)
        if power < 0:
            lower = c.make_mpf((atom._mpi_[0], atom._mpi_[0]))
            upper = c.make_mpf((atom._mpi_[1], atom._mpi_[1]))
            if lower <= 0 <= upper:
                raise ZeroDivisionError(f"formal atom {name!r} contains zero for negative power")
            for _ in range(-power):
                value = value / atom
        else:
            for _ in range(power):
                value = value * atom
    return value


def _state_container(state: Mapping[str, Any]) -> Mapping[Any, Any]:
    for key in ("W", "state", "components", "jets", "state_jets"):
        value = state.get(key, _MISSING)
        if isinstance(value, Mapping):
            return value
    return state


def _component_entry(container: Mapping[Any, Any], component: int) -> Any:
    for key in (component, str(component), f"W{component}", f"W_{component}"):
        if key in container:
            return container[key]
    raise KeyError(f"caller state is missing W{component}")


def _state_order_entry(entry: Any, order: int) -> Any:
    if isinstance(entry, Mapping):
        for wrapper in ("ordinary_Z", "jets", "values", "state"):
            nested = entry.get(wrapper, _MISSING)
            if isinstance(nested, Mapping) and nested is not entry:
                try:
                    return _state_order_entry(nested, order)
                except KeyError:
                    pass
        for key in (order, str(order), f"Z{order}", f"ordinary_Z{order}", f"W_Z{order}"):
            if key in entry:
                return entry[key]
        if order == 0:
            for key in ("value", "coefficient", "coefficient_interval", "interval", "enclosure"):
                if key in entry:
                    return entry
        raise KeyError(f"caller state is missing ordinary Z{order}")
    if order == 0:
        return entry
    raise KeyError(f"caller state is missing ordinary Z{order}")


def _state_value(c: Any, value: Any, atoms: Mapping[Any, Any], label: str) -> Any:
    """Decode one caller-owned state jet, retaining any formal powers."""

    if isinstance(value, Mapping):
        if "terms" in value and _is_sequence(value["terms"]):
            total = c.mpf(0)
            for index, item in enumerate(value["terms"]):
                total = total + _state_value(c, item, atoms, f"{label}.terms[{index}]")
            return total
        coefficient = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
        if coefficient is not _MISSING:
            powers = value.get("powers", value.get("formal_powers", _MISSING))
            if powers is _MISSING:
                raise KeyError(f"{label}.powers is required for a formal state mapping")
            if not _is_sequence(powers) or len(powers) != 4:
                raise TypeError(f"{label}.powers must contain four integer powers")
            if any(isinstance(item, bool) or not isinstance(item, int) for item in powers):
                raise TypeError(f"{label}.powers must contain four integer powers")
            return _interval(c, coefficient, label + ".coefficient") * _formal_factor(c, atoms, powers)
    if _is_sequence(value) and not (len(value) == 2 and not any(isinstance(item, Mapping) for item in value)):
        total = c.mpf(0)
        for index, item in enumerate(value):
            total = total + _state_value(c, item, atoms, f"{label}[{index}]")
        return total
    return _interval(c, value, label)


def _state_factor(c: Any, state: Mapping[str, Any], reference: Mapping[str, Any]) -> Any:
    if not isinstance(state, Mapping):
        raise TypeError("state must be a mapping of caller-owned W jets")
    component = reference.get("component")
    order = reference.get("ordinary_Z_order")
    if isinstance(component, bool) or not isinstance(component, int) or component not in STATE_COMPONENTS:
        raise ValueError("invalid state component reference")
    if isinstance(order, bool) or not isinstance(order, int) or order not in STATE_Z_ORDERS:
        raise ValueError("invalid state ordinary Z order reference")
    # Flat caller maps are useful when a source receipt already names jets as
    # W1_Z0, W2_Z1, ... .  They remain caller-owned and are never synthesized.
    for key in (
        f"W{component}_Z{order}",
        f"W_{component}_Z{order}",
        f"W{component}Z{order}",
        (component, order),
    ):
        if key in state:
            return _state_value(c, state[key], _formal_atoms(state), f"W{component}_Z{order}")
    container = _state_container(state)
    entry = _component_entry(container, component)
    selected = _state_order_entry(entry, order)
    return _state_value(c, selected, _formal_atoms(state), f"W{component}_Z{order}")


def _term_value(c: Any, term: Mapping[str, Any], state: Mapping[str, Any]) -> Any:
    value = _interval(c, term["coefficient"], "expression coefficient")
    scalar = term.get("scalar", 1)
    if isinstance(scalar, bool) or not isinstance(scalar, int):
        raise TypeError("expression scalar must be an integer")
    value = value * scalar
    value = value * _formal_factor(c, _formal_atoms(state), term.get("powers", [0, 0, 0, 0]))
    reference = term.get("state_reference", term.get("state_ref"))
    if reference is not None:
        value = value * _state_factor(c, state, reference)
    return value


def _row_terms(expression: Mapping[str, Any], label: str) -> Sequence[Mapping[str, Any]]:
    rows = expression.get("rows")
    if isinstance(rows, Mapping) and label in rows:
        values = rows[label]
    else:
        values = expression.get(label, _MISSING)
    if values is _MISSING or not _is_sequence(values):
        raise KeyError(f"expression is missing additive row {label}")
    return values


def evaluate(c: Any, expression: Mapping[str, Any], state: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate a built expression against complete caller-owned intervals.

    The returned ``h5`` and ``h5_Z1`` are directed enclosures.  The
    ``term_values`` rows retain one interval per source-linked additive term,
    so the final interval sum does not erase provenance.  ``state`` must
    provide W components 1, 2, 3 and 5 wherever referenced, through ordinary
    Z order 2, plus formal atom intervals for every nonzero term power.
    """

    _validate_context(c)
    if not isinstance(expression, Mapping):
        raise TypeError("expression must be the mapping returned by build")
    if not isinstance(state, Mapping):
        raise TypeError("state must be a mapping of caller-owned intervals")
    outputs: dict[str, Any] = {}
    term_values: dict[str, list[Any]] = {}
    term_enclosures: dict[str, list[dict[str, Any]]] = {}
    for label in ("h5", "h5_Z1"):
        values: list[Any] = []
        records: list[dict[str, Any]] = []
        for index, term in enumerate(_row_terms(expression, label)):
            if not isinstance(term, Mapping):
                raise TypeError(f"{label}[{index}] must be a term mapping")
            value = _term_value(c, term, state)
            values.append(value)
            records.append({"index": index, "value": value, "term": term})
        total = c.mpf(0)
        for value in values:
            total = total + value
        outputs[label] = total
        term_values[label] = values
        term_enclosures[label] = records
    return {
        "schema": SCHEMA,
        "status": "caller_evaluated_directed_row5_affine_forcing",
        "h5": outputs["h5"],
        "h5_Z1": outputs["h5_Z1"],
        "rows": outputs,
        "term_values": term_values,
        "term_enclosures": term_enclosures,
        "directed": True,
        "interval_sum_performed": True,
        "source_terms_individual": True,
        "formal_atoms_caller_owned": True,
        "state_caller_owned": True,
        "missing_state_is_not_zero": True,
        "physical_h5_materialized": False,
        "complete_W1": False,
        "ns_validation": False,
    }


build_row5_affine_forcing = build
row5_affine_forcing = build
enclose = build
evaluate_expression = evaluate

