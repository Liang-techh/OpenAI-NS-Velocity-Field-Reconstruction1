"""Independent v2 finite physical ``m15`` sectors for the first ``Kchi`` feedback words.

The input to this leaf is deliberately small and source-bound: eight target
cells ``t in [i/8, (i+1)/8]`` with ``rho = 4*t**2``, already physical leading
``U``/``F`` rows, and four finite sector branches.  The branches are the
signed first ``W0`` chi word and the three signed ``W1`` routes
(``E01``, ``E10``, and ``E11``).  No producer, archive,
tail estimate, or pressure/g4 term is consulted here.

For a target cell the exact directed radial weights are

``2*epsilon*integral(Ulead*Uz1*drho)``

and

``-2*epsilon**2*integral(rho*Flead*F1*drho)``.

Thus the cell factors are ``8*(t_hi**2-t_lo**2)`` and
``-16*(t_hi**4-t_lo**4)``; the external formal powers are respectively
``Lambda**-1`` and ``Lambda**-2``.  Ordinary ``Z1`` is the two-term product
rule.  Input fields may be absent, in which case that field is the natural
zero.  A symmetric error or a tail bound is retained in a separate error
stream and is never inserted into the signed field stream.

This is a finite sector calculation only.  It does not reconstruct the
already physical leading fields or their ``F0`` connection, does not add
``g4``/pressure, and does not claim a complete ``m15`` response.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


COMPONENT = "n1_kchi_feedback_m15"
SCHEMA = "n1_kchi_feedback_m15_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
OUTER_CELL_COUNT = 8
FIELDS = ("F1", "Uz1")
# Canonical finite graph: the W0 chi-first word plus three independently
# source-bound W1 products.  Historical chi/B0/B1 spellings remain aliases;
# E01/E10/E11 are never merged.
BRANCHES = ("W0_chi", "E01", "E10", "E11")
REQUIRED_CERTIFICATES = ()

_FIELD_ALIASES: dict[str, tuple[Any, ...]] = {
    "F1": ("F1", "f1", 1, "1"),
    "Uz1": ("Uz1", "UZ1", "U_z1", "uz1", "Uz", 2, "2"),
}
_ORDER_ALIASES: dict[int, tuple[Any, ...]] = {
    0: (0, "0", "Z0", "z0", "ordinary_Z0"),
    1: (1, "1", "Z1", "z1", "ordinary_Z1"),
}
_BRANCH_ALIASES: dict[str, tuple[Any, ...]] = {
    "W0_chi": (
        "W0_chi",
        "w0_chi",
        "chi",
        "χ",
        "chi_first",
        "W0_chi_first",
        "w0_chi_first",
        "initial_function",
        "computed_prefix",
    ),
    # Historical B0/B1 spellings are aliases for E01/E10 only.  E11 is
    # always a distinct route when supplied.
    "E01": (
        "E01",
        "e01",
        "feedback_E01",
        "w1_E01",
        "W1_E01",
        "B0_feedback",
        "b0_feedback",
        "feedback_B0",
        "feedback0",
        "B0",
        "w1_B0",
        "W1_B0",
    ),
    "E10": (
        "E10",
        "e10",
        "feedback_E10",
        "w1_E10",
        "W1_E10",
        "B1_feedback",
        "b1_feedback",
        "feedback_B1",
        "feedback1",
        "B1",
        "w1_B1",
        "W1_B1",
        "B1_derivative_feedback",
        "b1_derivative_feedback",
    ),
    "E11": (
        "E11",
        "e11",
        "feedback_E11",
        "w1_E11",
        "W1_E11",
        "B1_derivative_delta_feedback",
    ),
}
_ERROR_KEYS = (
    "error",
    "errors",
    "error_sectors",
    "symmetric_error",
    "symmetric_errors",
    "tail",
    "tail_bound",
    "remainder",
    "uncertainty_sectors",
)


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "OUTER_CELL_COUNT",
    "FIELDS",
    "BRANCHES",
    "REQUIRED_CERTIFICATES",
    "build",
    "build_kchi_feedback_m15",
    "adapt",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    """Return ordinary scalar bounds for validation only."""

    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    from mpmath import mp

    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _freeze(value: Any) -> Any:
    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Iterable[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Collect source IDs, classes, and uncertainty without midpointing."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {"source_ids": [], "error_classes": [], "uncertainty": False}

    if isinstance(value, Mapping):
        _seen.add(marker)
        ids = list(_as_items(value.get("source_ids", value.get("sources"))))
        if not ids:
            ids = list(_as_items(value.get("source_id", value.get("source"))))
        classes = list(
            _as_items(
                value.get(
                    "error_classes",
                    value.get("uncertainty_class", value.get("error_class")),
                )
            )
        )
        if value.get("error_class") is not None:
            classes.extend(_as_items(value.get("error_class")))
        operation = str(value.get("operation", "")).lower()
        kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        bound_type = str(value.get("bound_type", "")).lower()
        uncertain = bool(
            value.get(
                "uncertainty",
                value.get("is_uncertainty", value.get("symmetric_error", False)),
            )
        )
        uncertain = uncertain or bool(value.get("error_only", False))
        uncertain = uncertain or any(
            token in operation or token in kind or token in bound_type
            for token in ("cauchy", "symmetric", "uncertainty", "tail", "remainder")
        )
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            children = (
                child
                if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray))
                else (child,)
            )
            for item in children:
                info = _origin_info(item, _seen=_seen)
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
        }

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = False
        for item in value:
            info = _origin_info(item, _seen=_seen)
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": uncertain,
        }

    return {"source_ids": [], "error_classes": [], "uncertainty": False}


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
    }


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    for item in inputs:
        info = _origin_info(item)
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or bool(info["uncertainty"])
    ids.extend(_as_items(metadata.get("source_ids")))
    classes.extend(_as_items(metadata.get("error_classes")))
    origin: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "symbolic_cancellation_claimed": False,
    }
    if inputs:
        origin["parents"] = tuple(_origin_ref(item) for item in inputs)
    for key, value in metadata.items():
        if key not in {"source_ids", "error_classes", "uncertainty"}:
            origin[key] = value
    return origin


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not isinstance(power, int) or isinstance(power, bool) for power in value)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(power) for power in value)  # type: ignore[return-value]


def _term_powers(value: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    if "powers" in value:
        return _check_powers(value["powers"], f"{label}.powers")
    if any(
        key in value
        for key in ("hb_power", "F0_power", "Pstar_squared_power", "Pstar_power", "Lambda_power")
    ):
        return _check_powers(
            (
                value.get("hb_power", 0),
                value.get("F0_power", 0),
                value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
                value.get("Lambda_power", 0),
            ),
            f"{label}.formal powers",
        )
    # Zero is the neutral formal sector and is useful for synthetic callers.
    return (0, 0, 0, 0)


def _term_marked_error(value: Mapping[str, Any]) -> bool:
    """Recognize explicit bound markers only.

    kind=source_interval_uncertainty and similar provenance labels describe an
    uncertain signed source interval; they are not kernel-tail records.
    """
    if any(
        bool(value.get(key, False))
        for key in ("tail_bound", "remainder", "error_only", "symmetric_error", "is_tail_bound")
    ):
        return True
    bound_type = str(value.get("bound_type", "")).lower()
    kind = str(value.get("kind", "")).lower()
    explicit = {
        "tail",
        "tail_bound",
        "kernel_tail",
        "kchi_tail",
        "remainder",
        "symmetric",
        "symmetric_error",
        "error",
    }
    return bound_type in explicit or kind in explicit


def _check_term(
    c: Any,
    value: Any,
    label: str,
    *,
    signed_field: bool,
    error_kind: str | None = None,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    coefficient = value.get("coefficient", value.get("value", value.get("interval")))
    if coefficient is None:
        raise ValueError(f"{label} requires coefficient")
    if not _is_interval(coefficient):
        coefficient = c.mpf(coefficient)
    powers = _term_powers(value, label)
    raw_origin = value.get("origin", value.get("metadata", {}))
    info = _origin_info(raw_origin)
    ids = list(info["source_ids"])
    for key in ("source_ids", "source_id", "input_id", "source"):
        if value.get(key) is not None:
            ids.extend(_as_items(value[key]))
    ids = _unique(ids)
    if not ids:
        raise ValueError(f"{label}.origin must retain source_ids")
    classes = list(info["error_classes"])
    classes.extend(_as_items(value.get("error_classes")))
    classes.extend(_as_items(value.get("error_class")))
    classes = _unique(classes)
    uncertain = bool(info["uncertainty"])
    uncertain = uncertain or bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    if not signed_field:
        uncertain = True
    origin = dict(raw_origin) if isinstance(raw_origin, Mapping) else {}
    for marker in ("tail_bound", "remainder", "error_only", "symmetric_error", "bound_type", "kind"):
        if marker in value:
            origin[marker] = value[marker]
    origin.update(
        {
            "source_ids": list(ids),
            "error_classes": list(classes),
            "uncertainty": bool(uncertain),
            "signed_field": bool(signed_field),
            "symbolic_cancellation_claimed": False,
        }
    )
    if error_kind is not None:
        origin["bound_type"] = error_kind
        origin["error_only"] = True
    return {
        "powers": powers,
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": ids,
        "error_classes": classes,
        "uncertainty": bool(uncertain),
        "signed_field": bool(signed_field),
    }


def _is_term_mapping(value: Any) -> bool:
    return isinstance(value, Mapping) and any(
        key in value
        for key in ("coefficient", "value", "interval", "powers", "F0_power", "Lambda_power")
    )


def _term_list(
    c: Any,
    value: Any,
    label: str,
    *,
    signed_field: bool,
    error_kind: str | None = None,
) -> list[dict[str, Any]]:
    if value is None:
        return []
    if _is_term_mapping(value):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    result: list[dict[str, Any]] = []
    for index, item in enumerate(values):
        if not isinstance(item, Mapping):
            raise TypeError(f"{label}[{index}] must be a sector mapping")
        result.append(
            _check_term(
                c,
                item,
                f"{label}[{index}]",
                signed_field=signed_field,
                error_kind=error_kind,
            )
        )
    return result


def _order_key(value: Any, label: str) -> int | None:
    if isinstance(value, bool):
        return None
    for order, aliases in _ORDER_ALIASES.items():
        if value in aliases:
            return order
    if isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError:
            return None
        return parsed if parsed in ORDINARY_Z_ORDERS else None
    return value if isinstance(value, int) and value in ORDINARY_Z_ORDERS else None


def _field_key(value: Any) -> str | None:
    for field, aliases in _FIELD_ALIASES.items():
        if value in aliases:
            return field
    return None


def _empty_fields() -> dict[int, dict[str, list[dict[str, Any]]]]:
    return {
        order: {field: [] for field in FIELDS}
        for order in ORDINARY_Z_ORDERS
    }


def _merge_field_rows(
    c: Any,
    target: dict[int, dict[str, list[dict[str, Any]]]],
    value: Any,
    label: str,
    *,
    signed_field: bool,
    error_kind: str | None,
) -> bool:
    """Parse one branch payload into order/field rows.

    The parser accepts the two common forms used by saved sectors:
    ``{Z0: {F1: ..., Uz1: ...}, Z1: ...}`` and
    ``{F1: {Z0: ..., Z1: ...}, Uz1: ...}``.  Numeric columns ``1`` and ``2``
    are accepted as the legacy F1/Uz1 names.
    """

    if value is None:
        return False
    if isinstance(value, Mapping) and _is_term_mapping(value):
        # A direct term cannot identify a field.  The caller only invokes this
        # path for a field-specific payload; keep it as Z0 in that case.
        return False
    if isinstance(value, Mapping):
        for wrapper in ("ordinary_Z", "ordinary", "rows", "fields", "sector_fields"):
            if wrapper in value and isinstance(value[wrapper], (Mapping, Sequence)):
                return _merge_field_rows(
                    c,
                    target,
                    value[wrapper],
                    f"{label}.{wrapper}",
                    signed_field=signed_field,
                    error_kind=error_kind,
                )

        seen = False
        # First recognize rows keyed by ordinary Z order.
        order_items: list[tuple[int, Any]] = []
        for raw_key, row in value.items():
            order = _order_key(raw_key, label)
            if order is not None:
                order_items.append((order, row))
        if order_items:
            for order, row in order_items:
                if isinstance(row, Mapping):
                    for raw_field, sectors in row.items():
                        field = _field_key(raw_field)
                        if field is None:
                            continue
                        field_signed, field_errors = (
                            _split_error_wrappers(sectors)
                            if isinstance(sectors, Mapping)
                            else (sectors, [])
                        )
                        target[order][field].extend(
                            _term_list(
                                c,
                                field_signed,
                                f"{label}[{order}].{field}",
                                signed_field=signed_field,
                                error_kind=error_kind,
                            )
                        )
                        for error_key, error_payload in field_errors:
                            target[order][field].extend(
                                _term_list(
                                    c,
                                    error_payload,
                                    f"{label}[{order}].{field}.{error_key}",
                                    signed_field=False,
                                    error_kind=(
                                        "symmetric_error"
                                        if "symmetric" in error_key
                                        else error_key
                                    ),
                                )
                            )
                        seen = True
                elif isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                    # A row with no field label is not assigned: silently
                    # assigning it to both fields would duplicate a source.
                    continue
            if seen:
                return True

        # Then recognize fields keyed at the outer level.
        for raw_field, payload in value.items():
            field = _field_key(raw_field)
            if field is None:
                continue
            seen = True
            if isinstance(payload, Mapping) and not _is_term_mapping(payload):
                nested_orders: list[tuple[int, Any]] = []
                for raw_order, sectors in payload.items():
                    order = _order_key(raw_order, label)
                    if order is not None:
                        nested_orders.append((order, sectors))
                if nested_orders:
                    for order, sectors in nested_orders:
                        target[order][field].extend(
                            _term_list(
                                c,
                                sectors,
                                f"{label}.{field}[{order}]",
                                signed_field=signed_field,
                                error_kind=error_kind,
                            )
                        )
                else:
                    target[0][field].extend(
                        _term_list(
                            c,
                            payload,
                            f"{label}.{field}",
                            signed_field=signed_field,
                            error_kind=error_kind,
                        )
                    )
            else:
                target[0][field].extend(
                    _term_list(
                        c,
                        payload,
                        f"{label}.{field}",
                        signed_field=signed_field,
                        error_kind=error_kind,
                    )
                )
        if seen:
            return True
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        # A two-row sequence is an ordinary-Z sequence if each row is a field
        # mapping.  Otherwise it is a Z0 sector list for the field-specific
        # caller, handled by _merge_field_payload below.
        if len(value) == 2 and all(isinstance(row, Mapping) for row in value):
            for order, row in enumerate(value):
                _merge_field_rows(
                    c,
                    target,
                    row,
                    f"{label}[{order}]",
                    signed_field=signed_field,
                    error_kind=error_kind,
                )
            return True
    return False


def _merge_field_payload(
    c: Any,
    target: dict[int, dict[str, list[dict[str, Any]]]],
    field: str,
    payload: Any,
    label: str,
    *,
    signed_field: bool,
    error_kind: str | None,
) -> None:
    """Parse a field-specific payload, including a direct Z0 term list."""

    if payload is None:
        return
    if isinstance(payload, Mapping) and _is_term_mapping(payload):
        target[0][field].extend(
            _term_list(c, payload, label, signed_field=signed_field, error_kind=error_kind)
        )
        return
    if isinstance(payload, Mapping):
        nested_orders: list[tuple[int, Any]] = []
        for raw_order, sectors in payload.items():
            order = _order_key(raw_order, label)
            if order is not None:
                nested_orders.append((order, sectors))
        if nested_orders:
            for order, sectors in nested_orders:
                target[order][field].extend(
                    _term_list(
                        c,
                        sectors,
                        f"{label}[{order}]",
                        signed_field=signed_field,
                        error_kind=error_kind,
                    )
                )
            return
    if isinstance(payload, Sequence) and not isinstance(payload, (str, bytes, bytearray)):
        if len(payload) == 2 and all(
            isinstance(row, Mapping) and not _is_term_mapping(row) for row in payload
        ):
            # Explicit two-row ordinary-Z field payload.
            for order, sectors in enumerate(payload):
                target[order][field].extend(
                    _term_list(
                        c,
                        sectors,
                        f"{label}[{order}]",
                        signed_field=signed_field,
                        error_kind=error_kind,
                    )
                )
            return
        target[0][field].extend(
            _term_list(c, payload, label, signed_field=signed_field, error_kind=error_kind)
        )
        return
    raise TypeError(f"{label} is not a field sector payload")


def _split_error_wrappers(value: Mapping[str, Any]) -> tuple[Any, list[tuple[str, Any]]]:
    signed: Any = value
    errors: list[tuple[str, Any]] = []
    for key in _ERROR_KEYS:
        if key in value:
            errors.append((key, value[key]))
    if any(key in value for key in ("field", "signed", "signed_field", "sectors", "terms", "value_sectors")):
        for key in ("field", "signed", "signed_field", "sectors", "terms", "value_sectors"):
            if key in value:
                signed = value[key]
                break
    elif errors and not any(_field_key(key) is not None for key in value):
        signed = None
    return signed, errors


def _normalize_branch(c: Any, value: Any, label: str) -> dict[str, Any]:
    fields = _empty_fields()
    errors = _empty_fields()
    if value is None:
        return {"fields": fields, "errors": errors, "raw": value}
    signed_value = value
    error_values: list[tuple[str, Any]] = []
    if isinstance(value, Mapping):
        signed_value, error_values = _split_error_wrappers(value)

        # A direct generic W0/W1 ordinary-Z wrapper is the signed payload.
        for wrapper in (
            "W0_ordinary_Z",
            "W1_ordinary_Z",
            "ordinary_Z",
            "ordinary",
            "sectors",
            "fields",
        ):
            if wrapper in value and signed_value is value:
                signed_value = value[wrapper]
                break

        # Parse direct field names, order rows, or legacy numeric columns.
        parsed = False
        if isinstance(signed_value, Mapping):
            for raw_field, payload in signed_value.items():
                field = _field_key(raw_field)
                if field is not None:
                    field_signed, field_errors = (
                        _split_error_wrappers(payload)
                        if isinstance(payload, Mapping)
                        else (payload, [])
                    )
                    _merge_field_payload(
                        c,
                        fields,
                        field,
                        field_signed,
                        f"{label}.{field}",
                        signed_field=True,
                        error_kind=None,
                    )
                    for error_key, error_payload in field_errors:
                        _merge_field_payload(
                            c,
                            errors,
                            field,
                            error_payload,
                            f"{label}.{field}.{error_key}",
                            signed_field=False,
                            error_kind=(
                                "symmetric_error" if "symmetric" in error_key else error_key
                            ),
                        )
                    parsed = True
        if not parsed:
            _merge_field_rows(
                c,
                fields,
                signed_value,
                label,
                signed_field=True,
                error_kind=None,
            )

        for key, payload in error_values:
            error_kind = "symmetric_error" if "symmetric" in key else key
            # Field-specific error wrappers are preferred.  A direct error
            # sector has no safe field inference and is therefore kept in the
            # explicit untyped error row for auditability.
            if isinstance(payload, Mapping):
                parsed_error = False
                for raw_field, field_payload in payload.items():
                    field = _field_key(raw_field)
                    if field is None:
                        continue
                    _merge_field_payload(
                        c,
                        errors,
                        field,
                        field_payload,
                        f"{label}.{key}.{field}",
                        signed_field=False,
                        error_kind=error_kind,
                    )
                    parsed_error = True
                if not parsed_error:
                    _merge_field_rows(
                        c,
                        errors,
                        payload,
                        f"{label}.{key}",
                        signed_field=False,
                        error_kind=error_kind,
                    )
            else:
                _merge_field_rows(
                    c,
                    errors,
                    payload,
                    f"{label}.{key}",
                    signed_field=False,
                    error_kind=error_kind,
                )

        # Terms marked as tails/errors inside a signed field are moved out of
        # the signed stream.  This is intentionally conservative: no signed
        # field is synthesized from a bound-only term.
        for order in ORDINARY_Z_ORDERS:
            for field in FIELDS:
                kept: list[dict[str, Any]] = []
                for term in fields[order][field]:
                    if _term_marked_error(term["origin"]):
                        term["origin"]["signed_field"] = False
                        term["origin"]["error_only"] = True
                        term["origin"].setdefault("bound_type", "tail_bound")
                        term["uncertainty"] = True
                        errors[order][field].append(term)
                    else:
                        kept.append(term)
                fields[order][field] = kept
    else:
        _merge_field_rows(c, fields, signed_value, label, signed_field=True, error_kind=None)
    return {"fields": fields, "errors": errors, "raw": value}


def _lookup_branch(mapping: Mapping[Any, Any], branch: str) -> Any:
    for key in _BRANCH_ALIASES.get(branch, (branch,)):
        if key in mapping:
            return mapping[key]
    return None


def _cell_branch_value(cell: Mapping[Any, Any], branch: str) -> Any:
    containers: list[Any] = []
    for key in ("branches", "feedback", "W1_feedback", "w1_feedback", "finite_sectors"):
        item = cell.get(key)
        if isinstance(item, Mapping):
            containers.append(item)
    for container in containers:
        found = _lookup_branch(container, branch)
        if found is not None:
            return found
    found = _lookup_branch(cell, branch)
    if found is not None:
        return found
    if branch in {"W0_chi", "chi", "chi_first"}:
        for key in ("W0_ordinary_Z", "W0_chi_ordinary_Z", "chi_ordinary_Z", "W0_chi"):
            if key in cell:
                return cell[key]
    return None


def _cell_input_value(value: Any, branch: str, index: int) -> Any:
    """Resolve an optional separate eight-cell branch input."""

    if value is None:
        return None
    if isinstance(value, Mapping):
        direct = _lookup_branch(value, branch)
        if direct is not None and direct is not value:
            return _cell_input_value(direct, branch, index)
        for key in ("cells", "outer_cells", "profile", "target_cells"):
            if key in value:
                return _cell_input_value(value[key], branch, index)
        for key in (index, str(index)):
            if key in value:
                return value[key]
        return value
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == OUTER_CELL_COUNT:
            return value[index]
    return value


def _canonical_branch_name(value: Any) -> str:
    for branch, aliases in _BRANCH_ALIASES.items():
        if value in aliases:
            return branch
    return str(value)


def _discover_branches(profile: Any, separate_branches: Mapping[str, Any]) -> tuple[str, ...]:
    """Discover explicit route names while retaining a stable canonical order."""

    if isinstance(profile, Mapping):
        for wrapper in ("outer_cells", "profile", "target_cells", "cells"):
            if wrapper in profile:
                profile = profile[wrapper]
                break
    names: list[str] = []

    def add(raw: Any) -> None:
        name = _canonical_branch_name(raw)
        if name not in names:
            names.append(name)

    def scan_mapping(mapping: Mapping[Any, Any]) -> None:
        for container_key in ("branches", "feedback", "W1_feedback", "w1_feedback", "finite_sectors"):
            container = mapping.get(container_key)
            if isinstance(container, Mapping):
                for raw in container:
                    add(raw)
        for raw in mapping:
            if _canonical_branch_name(raw) in BRANCHES or raw in {
                "W0_ordinary_Z",
                "W0_chi_ordinary_Z",
                "chi_ordinary_Z",
            }:
                if raw not in {"W0_ordinary_Z", "W0_chi_ordinary_Z", "chi_ordinary_Z"}:
                    add(raw)
            elif raw not in {"t", "rho", "branches", "feedback", "W1_feedback", "w1_feedback", "finite_sectors"}:
                # A direct unknown route is accepted if its value resembles a
                # sector packet; geometry/metadata keys are ignored.
                if isinstance(mapping[raw], (Mapping, Sequence)) and not isinstance(mapping[raw], (str, bytes, bytearray)):
                    add(raw)

    if isinstance(profile, Sequence) and not isinstance(profile, (str, bytes, bytearray)):
        for cell in profile:
            if isinstance(cell, Mapping):
                scan_mapping(cell)
    for raw, value in separate_branches.items():
        if value is not None:
            add(raw)
            if raw in {"w1_feedback", "feedback_branches"} and isinstance(value, Mapping):
                for nested in value:
                    add(nested)
    # A generic W0 ordinary-Z cell is the chi-first route even when no explicit
    # branch wrapper was supplied.
    if isinstance(profile, Sequence) and any(
        isinstance(cell, Mapping)
        and any(key in cell for key in ("W0_ordinary_Z", "W0_chi_ordinary_Z", "chi_ordinary_Z"))
        for cell in profile
    ):
        add("W0_chi")
    if not names:
        return BRANCHES
    # Keep the W0 chi route first, followed by the conventional E routes, then
    # any caller-defined names in their source order.
    preferred = [name for name in BRANCHES if name in names]
    return tuple(preferred + [name for name in names if name not in preferred])


def _normalize_profile(
    c: Any,
    profile: Any,
    separate_branches: Mapping[str, Any],
    branch_names: Sequence[str],
) -> list[dict[str, Any]]:
    if isinstance(profile, Mapping):
        for wrapper in ("outer_cells", "profile", "target_cells", "cells"):
            if wrapper in profile:
                profile = profile[wrapper]
                break
    if not isinstance(profile, Sequence) or isinstance(profile, (str, bytes, bytearray)):
        raise TypeError("profile must be a sequence of eight target cells")
    if len(profile) != OUTER_CELL_COUNT:
        raise ValueError("profile must contain exactly eight target cells")

    normalized: list[dict[str, Any]] = []
    for index, cell in enumerate(profile):
        label = f"profile[{index}]"
        if not isinstance(cell, Mapping):
            raise TypeError(f"{label} must be a cell mapping")
        if "t" not in cell or "rho" not in cell:
            raise ValueError(f"{label} requires t and rho")
        t, rho = cell["t"], cell["rho"]
        if not _is_interval(t) or not _is_interval(rho):
            raise TypeError(f"{label}.t and .rho must be directed intervals")
        left = c.mpf(index) / c.mpf(OUTER_CELL_COUNT)
        right = c.mpf(index + 1) / c.mpf(OUTER_CELL_COUNT)
        expected_t = _make_interval(c, left, right)
        expected_rho = c.mpf(4) * expected_t * expected_t
        if not _same_interval(t, expected_t) or not _same_interval(rho, expected_rho):
            raise ValueError(f"{label} is not the directed rho=4*t^2 target cell")

        packets: dict[str, dict[str, Any]] = {}
        for branch in branch_names:
            value = _cell_branch_value(cell, branch)
            if value is None:
                value = _cell_input_value(separate_branches.get(branch), branch, index)
            packets[branch] = _normalize_branch(c, value, f"{label}.{branch}")
        normalized.append(
            {
                "index": index,
                "t": t,
                "rho": rho,
                "branches": packets,
                "raw": cell,
            }
        )
    return normalized


def _normalize_leading(c: Any, value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    rows = {0: [], 1: []}
    if value is None:
        return rows
    if isinstance(value, Mapping):
        for wrapper in ("ordinary_Z", "ordinary", "leading", "rows"):
            if wrapper in value and isinstance(value[wrapper], (Mapping, Sequence)):
                return _normalize_leading(c, value[wrapper], f"{label}.{wrapper}")
        if _is_term_mapping(value):
            rows[0] = _term_list(c, value, f"{label}[0]", signed_field=True)
            return rows
        order_items = []
        for raw_order, sectors in value.items():
            order = _order_key(raw_order, label)
            if order is not None:
                order_items.append((order, sectors))
        if order_items:
            for order, sectors in order_items:
                rows[order] = _term_list(c, sectors, f"{label}[{order}]", signed_field=True)
            return rows
        # Accept a full leading U/F row when callers pass it directly.
        for key in ("U", "Uz", "Uz0", "velocity", "field", "sectors", "terms"):
            if key in value:
                return _normalize_leading(c, value[key], f"{label}.{key}")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2 and all(isinstance(item, (Mapping, Sequence)) for item in value):
            # A two-row ordinary-Z list is distinguishable from a two-term list
            # because each row is not a sector mapping with a coefficient.
            if not all(_is_term_mapping(item) for item in value):
                rows[0] = _term_list(c, value[0], f"{label}[0]", signed_field=True)
                rows[1] = _term_list(c, value[1], f"{label}[1]", signed_field=True)
                return rows
        rows[0] = _term_list(c, value, f"{label}[0]", signed_field=True)
        return rows
    raise TypeError(f"{label} must provide ordinary-Z sectors")


def _validate_source_atoms(source_atoms: Any) -> None:
    names: set[Any] = set()
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        for item in source_atoms:
            if isinstance(item, Mapping):
                name = item.get("name", item.get("atom", item.get("power")))
                if name is not None:
                    names.add(name)
            else:
                names.add(item)
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError(
            "source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda"
        )


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, tuple[Any, ...], bool, Any, Any]:
    info = _origin_info(term["origin"])
    classes = tuple(sorted((_freeze(item) for item in info["error_classes"]), key=repr))
    origin = term["origin"]
    return (
        tuple(term["powers"]),
        bool(info["uncertainty"]),
        classes,
        bool(origin.get("signed_field", term.get("signed_field", True))),
        origin.get("bound_type"),
        origin.get("derivative_piece"),
    )


def _aggregate(terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Outward-add only compatible sectors; never claim cancellation."""

    groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for term in terms:
        info = _origin_info(term["origin"])
        key = _group_key(term)
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(term["powers"]),
                "coefficient": term["coefficient"],
                "ids": {},
                "classes": {},
                "uncertainty": False,
                "signed_field": bool(term.get("signed_field", True)),
                "bound_type": term["origin"].get("bound_type"),
                "derivative_piece": term["origin"].get("derivative_piece"),
                "merged_term_count": 0,
            }
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for item in info["source_ids"]:
            group["ids"][_freeze(item)] = item
        for item in info["error_classes"]:
            group["classes"][_freeze(item)] = item
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["merged_term_count"] += 1

    result: list[dict[str, Any]] = []
    for group in groups.values():
        origin: dict[str, Any] = {
            "operation": "finite_m15_aggregate",
            "source_ids": list(group["ids"].values()),
            "error_classes": list(group["classes"].values()),
            "uncertainty": bool(group["uncertainty"]),
            "signed_field": bool(group["signed_field"]),
            "symbolic_cancellation_claimed": False,
            "merged_term_count": group["merged_term_count"],
        }
        if group["derivative_piece"] is not None:
            origin["derivative_piece"] = group["derivative_piece"]
        if group["bound_type"] is not None:
            origin["bound_type"] = group["bound_type"]
            origin["error_only"] = True
        term = {
            "powers": tuple(group["powers"]),
            "coefficient": group["coefficient"],
            "origin": origin,
            "signed_field": bool(group["signed_field"]),
        }
        # Legacy consumers inspect these two projections.
        term["F0_power"] = term["powers"][1]
        term["Lambda_power"] = term["powers"][3]
        result.append(term)
    return result


def _cell_weight(c: Any, index: int, channel: str) -> Any:
    left = c.mpf(index) / c.mpf(OUTER_CELL_COUNT)
    right = c.mpf(index + 1) / c.mpf(OUTER_CELL_COUNT)
    if channel == "axial":
        return c.mpf(8) * (right**2 - left**2)
    if channel == "swirl":
        return -c.mpf(16) * (right**4 - left**4)
    raise ValueError(f"unknown m15 channel {channel!r}")


def _product(
    leading: Mapping[str, Any],
    word: Mapping[str, Any],
    factor: Any,
    *,
    branch: str,
    channel: str,
    field: str,
    order: int,
    cell: int,
    derivative_piece: str,
    external_powers: tuple[int, int, int, int],
    signed_field: bool,
    error_kind: str | None = None,
) -> dict[str, Any]:
    powers = tuple(
        a + b + offset
        for a, b, offset in zip(leading["powers"], word["powers"], external_powers)
    )
    origin = _new_origin(
        "kchi_feedback_m15_product",
        leading["origin"],
        word["origin"],
        branch=branch,
        channel=channel,
        field=field,
        ordinary_Z_order=order,
        target_cell=cell,
        derivative_piece=derivative_piece,
        external_formal_powers=external_powers,
        signed_field=bool(signed_field),
        symbolic_cancellation_claimed=False,
    )
    if error_kind is not None:
        origin["bound_type"] = error_kind
        origin["error_only"] = True
        origin["uncertainty"] = True
    return {
        "powers": powers,
        "coefficient": factor * leading["coefficient"] * word["coefficient"],
        "origin": origin,
        "signed_field": bool(signed_field),
    }


def _assemble_branch(
    c: Any,
    cells: Sequence[Mapping[str, Any]],
    leading: Mapping[int, list[dict[str, Any]]],
    *,
    branch: str,
    field: str,
    channel: str,
    external_powers: tuple[int, int, int, int],
) -> tuple[dict[int, list[dict[str, Any]]], dict[int, list[dict[str, Any]]], dict[int, dict[int, list[dict[str, Any]]]], dict[int, dict[int, list[dict[str, Any]]]]]:
    signed_total: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    error_total: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    signed_cells: dict[int, dict[int, list[dict[str, Any]]]] = {
        order: {index: [] for index in range(OUTER_CELL_COUNT)}
        for order in ORDINARY_Z_ORDERS
    }
    error_cells: dict[int, dict[int, list[dict[str, Any]]]] = {
        order: {index: [] for index in range(OUTER_CELL_COUNT)}
        for order in ORDINARY_Z_ORDERS
    }
    for cell in cells:
        index = int(cell["index"])
        factor = _cell_weight(c, index, channel)
        packet = cell["branches"][branch]
        words = packet["fields"]
        error_words = packet["errors"]
        raw: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
        raw_errors: dict[int, list[dict[str, Any]]] = {0: [], 1: []}

        for left in leading[0]:
            for right in words[0][field]:
                raw[0].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=0,
                        cell=index,
                        derivative_piece="leading0*word0",
                        external_powers=external_powers,
                        signed_field=True,
                    )
                )
            for right in error_words[0][field]:
                raw_errors[0].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=0,
                        cell=index,
                        derivative_piece="leading0*error0",
                        external_powers=external_powers,
                        signed_field=False,
                        error_kind=right["origin"].get("bound_type", "symmetric_error"),
                    )
                )

        for left in leading[1]:
            for right in words[0][field]:
                raw[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=1,
                        cell=index,
                        derivative_piece="leading1*word0",
                        external_powers=external_powers,
                        signed_field=True,
                    )
                )
            for right in error_words[0][field]:
                raw_errors[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=1,
                        cell=index,
                        derivative_piece="leading1*error0",
                        external_powers=external_powers,
                        signed_field=False,
                        error_kind=right["origin"].get("bound_type", "symmetric_error"),
                    )
                )
        for left in leading[0]:
            for right in words[1][field]:
                raw[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=1,
                        cell=index,
                        derivative_piece="leading0*word1",
                        external_powers=external_powers,
                        signed_field=True,
                    )
                )
            for right in error_words[1][field]:
                raw_errors[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        branch=branch,
                        channel=channel,
                        field=field,
                        order=1,
                        cell=index,
                        derivative_piece="leading0*error1",
                        external_powers=external_powers,
                        signed_field=False,
                        error_kind=right["origin"].get("bound_type", "symmetric_error"),
                    )
                )

        signed_cells[0][index] = _aggregate(raw[0])
        signed_cells[1][index] = _aggregate(raw[1])
        error_cells[0][index] = _aggregate(raw_errors[0])
        error_cells[1][index] = _aggregate(raw_errors[1])
        signed_total[0].extend(raw[0])
        signed_total[1].extend(raw[1])
        error_total[0].extend(raw_errors[0])
        error_total[1].extend(raw_errors[1])

    return (
        {order: _aggregate(signed_total[order]) for order in ORDINARY_Z_ORDERS},
        {order: _aggregate(error_total[order]) for order in ORDINARY_Z_ORDERS},
        signed_cells,
        error_cells,
    )


def _coalesce(value: Any, aliases: Mapping[str, Any], names: Sequence[str], label: str) -> Any:
    candidates = [] if value is None else [(label, value)]
    candidates.extend((name, aliases[name]) for name in names if aliases.get(name) is not None)
    if not candidates:
        return None
    first_name, first = candidates[0]
    for name, candidate in candidates[1:]:
        if candidate is not first:
            raise TypeError(f"{label} was supplied more than once ({first_name}, {name})")
    return first


def build(
    c: Any,
    *positional: Any,
    profile: Any | None = None,
    leading_U: Any | None = None,
    leading_F: Any | None = None,
    w0_chi: Any | None = None,
    B0_feedback: Any | None = None,
    B1_derivative_feedback: Any | None = None,
    E01: Any | None = None,
    E10: Any | None = None,
    E11: Any | None = None,
    w1_feedback: Any | None = None,
    feedback_branches: Any | None = None,
    source_family: Any | None = None,
    source_atoms: Any | None = None,
    certificate: Any | None = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build finite physical axial/swirl sectors through ordinary ``Z1``.

    ``profile`` may carry each branch in a cell, or the optional
    ``w0_chi``, ``E01``, ``E10`` and ``E11`` (or a ``w1_feedback`` mapping)
    may provide separate eight-cell arrays.  Missing branch/order/field rows
    are zeros.  Explicit route names in a feedback mapping are retained as
    independent output branches.  The supplied ``leading_U`` and
    ``leading_F`` rows are used as already physical fields; their ``F0``
    connection is not reconstructed.
    """

    if positional:
        if len(positional) > 3:
            raise TypeError("build accepts at most profile, leading_U, leading_F positionally")
        positional_values = list(positional) + [None] * (3 - len(positional))
        profile = _coalesce(profile, {"positional_profile": positional_values[0]}, ("positional_profile",), "profile")
        leading_U = _coalesce(leading_U, {"positional_leading_U": positional_values[1]}, ("positional_leading_U",), "leading_U")
        leading_F = _coalesce(leading_F, {"positional_leading_F": positional_values[2]}, ("positional_leading_F",), "leading_F")

    if "g4" in aliases or "pressure" in aliases:
        raise ValueError("g4/pressure is not part of the finite kchi feedback m15 adapter")
    profile = _coalesce(
        profile,
        aliases,
        ("outer_cells", "cells", "target_cells", "w1_profile"),
        "profile",
    )
    leading_U = _coalesce(
        leading_U,
        aliases,
        ("leading_u", "U_lead", "Ulead", "leading_velocity"),
        "leading_U",
    )
    leading_F = _coalesce(
        leading_F,
        aliases,
        ("leading_f", "F_lead", "Flead", "leading_swirl"),
        "leading_F",
    )
    w0_chi = _coalesce(
        w0_chi,
        aliases,
        ("chi", "chi_first", "W0_chi", "w0_chi_first", "initial_function"),
        "w0_chi",
    )
    B0_feedback = _coalesce(
        B0_feedback,
        aliases,
        ("feedback_B0", "feedback0", "B0", "w1_B0", "W1_B0"),
        "B0_feedback",
    )
    B1_derivative_feedback = _coalesce(
        B1_derivative_feedback,
        aliases,
        ("feedback_B1", "feedback1", "B1", "w1_B1", "W1_B1"),
        "B1_derivative_feedback",
    )
    E01 = _coalesce(
        E01,
        aliases,
        ("e01", "feedback_E01", "w1_E01", "W1_E01"),
        "E01",
    )
    E10 = _coalesce(
        E10,
        aliases,
        ("e10", "feedback_E10", "w1_E10", "W1_E10"),
        "E10",
    )
    E11 = _coalesce(
        E11,
        aliases,
        ("e11", "feedback_E11", "w1_E11", "W1_E11"),
        "E11",
    )
    w1_feedback = _coalesce(
        w1_feedback,
        aliases,
        ("W1_feedback", "feedback", "w1_feedback_branches"),
        "w1_feedback",
    )
    feedback_branches = _coalesce(
        feedback_branches,
        aliases,
        ("feedback_routes",),
        "feedback_branches",
    )
    unexpected = set(aliases) - {
        "outer_cells",
        "cells",
        "target_cells",
        "w1_profile",
        "leading_u",
        "U_lead",
        "Ulead",
        "leading_velocity",
        "leading_f",
        "F_lead",
        "Flead",
        "leading_swirl",
        "chi",
        "chi_first",
        "W0_chi",
        "w0_chi_first",
        "initial_function",
        "feedback_B0",
        "feedback0",
        "B0",
        "w1_B0",
        "W1_B0",
        "feedback_B1",
        "feedback1",
        "B1",
        "w1_B1",
        "W1_B1",
        "e01",
        "feedback_E01",
        "w1_E01",
        "W1_E01",
        "e10",
        "feedback_E10",
        "w1_E10",
        "W1_E10",
        "e11",
        "feedback_E11",
        "w1_E11",
        "W1_E11",
        "W1_feedback",
        "feedback",
        "w1_feedback_branches",
        "feedback_routes",
    }
    if unexpected:
        names = ", ".join(sorted(str(name) for name in unexpected))
        raise TypeError(f"unexpected kchi feedback arguments: {names}")
    if profile is None:
        raise TypeError("profile is required")
    if leading_U is None or leading_F is None:
        raise TypeError("leading_U and leading_F are required")
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is not None and isinstance(source_family, str) and not source_family.strip():
        raise ValueError("source_family must not be empty")
    if source_atoms is not None:
        _validate_source_atoms(source_atoms)
    if certificate is not None and not isinstance(certificate, Mapping):
        raise TypeError("certificate must be a mapping when supplied")

    leading_u = _normalize_leading(c, leading_U, "leading_U")
    leading_f = _normalize_leading(c, leading_F, "leading_F")
    separate_branches = {
        "W0_chi": w0_chi,
        "E01": E01 if E01 is not None else B0_feedback,
        "E10": E10 if E10 is not None else B1_derivative_feedback,
        "E11": E11,
    }
    if isinstance(w1_feedback, Mapping):
        separate_branches.update(w1_feedback)
    if isinstance(feedback_branches, Mapping):
        separate_branches.update(feedback_branches)
    branch_names = _discover_branches(profile, separate_branches)
    cells = _normalize_profile(c, profile, separate_branches, branch_names)

    branch_results: dict[str, dict[int, dict[str, Any]]] = {}
    axial_by_order: dict[int, dict[str, list[dict[str, Any]]]] = {
        order: {} for order in ORDINARY_Z_ORDERS
    }
    swirl_by_order: dict[int, dict[str, list[dict[str, Any]]]] = {
        order: {} for order in ORDINARY_Z_ORDERS
    }
    axial_errors_by_order: dict[int, dict[str, list[dict[str, Any]]]] = {
        order: {} for order in ORDINARY_Z_ORDERS
    }
    swirl_errors_by_order: dict[int, dict[str, list[dict[str, Any]]]] = {
        order: {} for order in ORDINARY_Z_ORDERS
    }
    for branch in branch_names:
        axial, axial_errors, axial_cells, axial_error_cells = _assemble_branch(
            c,
            cells,
            leading_u,
            branch=branch,
            field="Uz1",
            channel="axial",
            external_powers=(0, 0, 0, -1),
        )
        swirl, swirl_errors, swirl_cells, swirl_error_cells = _assemble_branch(
            c,
            cells,
            leading_f,
            branch=branch,
            field="F1",
            channel="swirl",
            external_powers=(0, 0, 0, -2),
        )
        branch_results[branch] = {
            order: {
                "axial": axial[order],
                "swirl": swirl[order],
                "axial_error": axial_errors[order],
                "swirl_error": swirl_errors[order],
                "symmetric_error_axial": axial_errors[order],
                "symmetric_error_swirl": swirl_errors[order],
                "axial_cells": axial_cells[order],
                "swirl_cells": swirl_cells[order],
                "axial_error_cells": axial_error_cells[order],
                "swirl_error_cells": swirl_error_cells[order],
            }
            for order in ORDINARY_Z_ORDERS
        }
        for order in ORDINARY_Z_ORDERS:
            axial_by_order[order][branch] = axial[order]
            swirl_by_order[order][branch] = swirl[order]
            axial_errors_by_order[order][branch] = axial_errors[order]
            swirl_errors_by_order[order][branch] = swirl_errors[order]

    target_cells: list[dict[str, Any]] = []
    for cell in cells:
        index = int(cell["index"])
        target_cells.append(
            {
                "index": index,
                "t": cell["t"],
                "rho": cell["rho"],
                "branches": {
                    branch: {
                        order: {
                            "axial": branch_results[branch][order]["axial_cells"][index],
                            "swirl": branch_results[branch][order]["swirl_cells"][index],
                            "axial_error": branch_results[branch][order]["axial_error_cells"][index],
                            "swirl_error": branch_results[branch][order]["swirl_error_cells"][index],
                        }
                        for order in ORDINARY_Z_ORDERS
                    }
                    for branch in branch_names
                },
            }
        )

    # A branch-oriented view is useful to callers assembling one source route
    # at a time.  The order-oriented views above preserve the existing adapter
    # convention and make ordinary Z0/Z1 explicit.
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": ["Z0", "Z1"],
        "outer_cell_count": OUTER_CELL_COUNT,
        "branches": branch_results,
        "branch_order": list(branch_names),
        "axial": axial_by_order,
        "swirl": swirl_by_order,
        "axial_by_branch": {
            branch: {order: axial_by_order[order][branch] for order in ORDINARY_Z_ORDERS}
            for branch in branch_names
        },
        "swirl_by_branch": {
            branch: {order: swirl_by_order[order][branch] for order in ORDINARY_Z_ORDERS}
            for branch in branch_names
        },
        "symmetric_error_axial": axial_errors_by_order,
        "symmetric_error_swirl": swirl_errors_by_order,
        "error_sectors": {
            branch: {
                order: {
                    "axial": axial_errors_by_order[order][branch],
                    "swirl": swirl_errors_by_order[order][branch],
                }
                for order in ORDINARY_Z_ORDERS
            }
            for branch in branch_names
        },
        "target_cells": target_cells,
        "profile": profile,
        "axial_weight_formula": "8*(t_hi^2-t_lo^2)",
        "swirl_weight_formula": "-16*(t_hi^4-t_lo^4)",
        "axial_formula": "2*epsilon*integral(Ulead*Uz1*drho)",
        "swirl_formula": "-2*epsilon^2*integral(rho*Flead*F1*drho)",
        "rho_formula": "4*t^2",
        "formal_axial_external_powers": (0, 0, 0, -1),
        "formal_swirl_external_powers": (0, 0, 0, -2),
        "ordinary_Z1_product_rule": True,
        "leading_U_is_physical": True,
        "leading_F_is_physical": True,
        "leading_F_includes_F0_connection": True,
        "F0_connection_added_here": False,
        "F0_recomputed": False,
        "g4_used": False,
        "g4_pressure_added": False,
        "pressure_included": False,
        "pressure_added": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "symmetric_errors_separate": True,
        "tail_bounds_as_signed_fields": False,
        "aggregation_key": "(powers, uncertainty, error_classes, signed_field, bound_type, derivative_piece)",
        "aggregation_is_outward_interval_addition_only": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "tail_added": False,
        "tail_included": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "m15_computed": False,
        "total_m15_computed": False,
        "finite_sector_only": True,
        "finite_kchi_feedback_m15_computed": True,
        "geometry_assumptions_caller_certified": certificate is not None,
        "geometry_certificate_independently_verified": False,
        "scope": "finite W0 chi-first plus independent E01/E10/E11 W1 feedback physical axial/swirl m15 sectors only; no g4/pressure or tail",
    }


build_kchi_feedback_m15 = build
adapt = build
