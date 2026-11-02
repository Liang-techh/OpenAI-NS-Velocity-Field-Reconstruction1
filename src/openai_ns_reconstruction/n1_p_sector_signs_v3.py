"""Exact sign-of-enclosure accounting for the saved finite ``P`` sectors.

The caller supplies the already authenticated finite-P packet.  Its shape is
``branch -> ordinary-Z order -> target cell -> raw term records`` with four
branches, two ordinary-Z rows, and eight cells.  This leaf deliberately does
not evaluate a coefficient, re-run a cell, add records, or choose a midpoint.
It reads only the two exact binary endpoint tuples of each record and emits
integer counts for the four possible closed-interval sign classes.

An endpoint is accepted either as a finite mpmath ``mpf``/interval object
(``_mpi_``/``_mpf_``) or as a saved exact four-integer mpf tuple.  The endpoint
tuples are retained in every output record.  Comparisons use integer binary
ordering only; no mpmath arithmetic is performed.  Non-finite endpoints and
reversed intervals are rejected.

``error_only``/``symmetric_error`` records and records carrying
``signed_field=False`` are never counted as signed records.  They remain in
separate streams, with their original metadata and endpoint classification.
The result therefore states only a conditional sign of the supplied
enclosure under the caller's upstream assumptions.  It proves no field sign,
source owner, gauge identity, terminal closure, or ledger statement.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_p_sector_signs_v3"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = (
    "G0_P_W0",
    "G0_P_DeltaW0",
    "D_P_W0",
    "D_P_DeltaW0",
)
ORDINARY_Z_ORDERS = (0, 1)
TARGET_CELLS = tuple(range(8))
SIGN_CATEGORIES = (
    "strictly_positive",
    "strictly_negative",
    "exact_zero",
    "contains_zero",
)
COUNT_KEYS = SIGN_CATEGORIES

_MISSING = object()
_ZERO = (0, 0, 0, 0)

__all__ = [
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "BRANCHES",
    "ORDINARY_Z_ORDERS",
    "TARGET_CELLS",
    "SIGN_CATEGORIES",
    "COUNT_KEYS",
    "build",
    "classify_interval",
    "classify_sign",
    "classify",
    "signs",
    "count_signs",
    "build_sign_report",
    "aggregate",
    "combine",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _is_context(value: Any) -> bool:
    """Recognize a caller interval context without invoking it."""

    return value is not None and hasattr(value, "make_mpf") and hasattr(value, "mpf")


def _lookup(mapping: Mapping[Any, Any], *names: Any) -> Any:
    """Return an exact key match, then a case-insensitive spelling match."""

    for name in names:
        try:
            if name in mapping:
                return mapping[name]
        except (TypeError, ValueError):
            continue
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        result = lowered.get(str(name).lower(), _MISSING)
        if result is not _MISSING:
            return result
    return _MISSING


def _map_key(mapping: Mapping[Any, Any], raw: Any, label: str) -> Any:
    """Find a branch/order/cell key without coercing arbitrary values."""

    try:
        if raw in mapping:
            return mapping[raw]
    except (TypeError, ValueError):
        pass
    text = str(raw)
    for key, value in mapping.items():
        if str(key) == text:
            return value
    raise KeyError(f"{label} is missing {raw!r}")


def _int_order(raw: Any, label: str) -> int:
    if type(raw) is int:
        value = raw
    elif isinstance(raw, str) and raw.strip().upper() in {"0", "Z0", "1", "Z1"}:
        text = raw.strip().upper()
        value = int(text.removeprefix("Z"))
    else:
        raise TypeError(f"{label} must be ordinary-Z order 0 or 1")
    if value not in ORDINARY_Z_ORDERS:
        raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    return value


def _int_cell(raw: Any, label: str) -> int:
    if type(raw) is int:
        value = raw
    elif isinstance(raw, str) and raw.strip().isdigit():
        value = int(raw.strip())
    else:
        raise TypeError(f"{label} must be target cell 0 through 7")
    if value not in TARGET_CELLS:
        raise ValueError(f"{label} must be target cell 0 through 7")
    return value


def _typed_int(raw: Any, label: str) -> int:
    if type(raw) is not int:
        raise TypeError(f"{label} must be an integer")
    return raw


def _endpoint_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    """Extract and validate one finite exact mpmath mpf tuple.

    mpmath stores a finite endpoint as ``(sign, mantissa, exponent, bc)``.
    A non-integer exponent/bit-count is how the special inf/NaN tuples are
    represented; rejecting those fields also rejects every non-finite value.
    """

    raw = getattr(value, "_mpf_", _MISSING)
    if raw is not _MISSING:
        value = raw
    if isinstance(value, Mapping):
        nested = _lookup(value, "exact_mpf_tuple", "mpf_tuple", "__exact_mpf_tuple__")
        if nested is not _MISSING and nested is not value:
            value = nested
        elif all(name in value for name in ("sign", "man", "exp", "bc")):
            value = (value["sign"], value["man"], value["exp"], value["bc"])
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label} must be a four-integer exact mpf tuple")
    if any(type(item) is not int for item in value):
        raise ValueError(f"{label} must be finite; special mpf endpoints are rejected")
    sign, mantissa, exponent, bit_count = tuple(value)
    if sign not in (0, 1) or mantissa < 0 or bit_count < 0:
        raise ValueError(f"{label} is not a valid finite mpf tuple")
    return sign, mantissa, exponent, bit_count


def _pair_from_mapping(value: Mapping[Any, Any], label: str) -> tuple[Any, Any] | None:
    lower = _lookup(
        value,
        "lower_exact_mpf_tuple",
        "lower_mpf_tuple",
        "lower_exact",
        "lower",
        "lo",
        "a",
    )
    upper = _lookup(
        value,
        "upper_exact_mpf_tuple",
        "upper_mpf_tuple",
        "upper_exact",
        "upper",
        "hi",
        "b",
    )
    if lower is not _MISSING or upper is not _MISSING:
        if lower is _MISSING or upper is _MISSING:
            raise ValueError(f"{label} has only one interval endpoint")
        return lower, upper
    nested = _lookup(
        value,
        "__exact_interval__",
        "exact_mpf_endpoints",
        "exact_endpoints",
        "endpoints",
        "interval_endpoints",
    )
    if nested is not _MISSING:
        if not _is_sequence(nested) or len(nested) != 2:
            raise TypeError(f"{label} exact endpoints must be a pair")
        return nested[0], nested[1]
    return None


def _interval_endpoints(value: Any, label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    """Read directed endpoints without hydrating or evaluating the interval."""

    raw_interval = getattr(value, "_mpi_", _MISSING)
    if raw_interval is not _MISSING:
        if not _is_sequence(raw_interval) or len(raw_interval) != 2:
            raise TypeError(f"{label} does not expose a directed mpf interval")
        lower_raw, upper_raw = raw_interval[0], raw_interval[1]
    elif isinstance(value, Mapping):
        pair = _pair_from_mapping(value, label)
        if pair is None:
            raise TypeError(f"{label} must contain exact binary interval endpoints")
        lower_raw, upper_raw = pair
    elif _is_sequence(value) and len(value) == 2:
        lower_raw, upper_raw = value[0], value[1]
    else:
        raise TypeError(f"{label} must be a directed interval or exact endpoint pair")

    lower = _endpoint_tuple(lower_raw, f"{label}.lower")
    upper = _endpoint_tuple(upper_raw, f"{label}.upper")
    if _compare_mpf(lower, upper) > 0:
        raise ValueError(f"{label} has reversed endpoints")
    return lower, upper


def _compare_magnitude(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    """Compare positive finite mpf magnitudes using integer binary order."""

    # The top occupied binary position decides first.  At equal top position
    # the exponent difference is bounded by the two mantissa bit lengths, so
    # these shifts never depend on an unbounded exponent-sized allocation.
    left_top = left[2] + left[1].bit_length()
    right_top = right[2] + right[1].bit_length()
    if left_top < right_top:
        return -1
    if left_top > right_top:
        return 1
    left_exp, right_exp = left[2], right[2]
    if left_exp >= right_exp:
        left_value = left[1] << (left_exp - right_exp)
        right_value = right[1]
    else:
        left_value = left[1]
        right_value = right[1] << (right_exp - left_exp)
    return (left_value > right_value) - (left_value < right_value)


def _compare_mpf(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    """Full exact ordering for finite mpf tuples, with no mpmath arithmetic."""

    left_zero = left[1] == 0
    right_zero = right[1] == 0
    if left_zero and right_zero:
        return 0
    if left_zero:
        return -1 if right[0] == 0 else 1
    if right_zero:
        return 1 if left[0] == 0 else -1
    if left[0] != right[0]:
        return -1 if left[0] == 1 else 1
    magnitude = _compare_magnitude(left, right)
    return -magnitude if left[0] == 1 else magnitude


def classify_interval(value: Any, label: str = "interval") -> dict[str, Any]:
    """Classify one closed directed interval by exact endpoint signs.

    The returned endpoint tuples are the original binary endpoints in a
    canonical tuple form.  No point value is selected.
    """

    lower, upper = _interval_endpoints(value, label)
    lower_sign = _compare_mpf(lower, _ZERO)
    upper_sign = _compare_mpf(upper, _ZERO)
    if lower_sign > 0:
        category = "strictly_positive"
    elif upper_sign < 0:
        category = "strictly_negative"
    elif lower_sign == 0 and upper_sign == 0:
        category = "exact_zero"
    else:
        category = "contains_zero"
    return {
        "classification": category,
        "sign_class": category,
        "lower_exact_mpf_tuple": lower,
        "upper_exact_mpf_tuple": upper,
        "exact_endpoint_tuple": (lower, upper),
        "midpoint_used": False,
    }


classify_sign = classify_interval
classify = classify_interval


def _serialise(value: Any) -> Any:
    """Make metadata JSON-safe while retaining exact endpoint tuples."""

    if isinstance(value, Mapping):
        result: dict[Any, Any] = {}
        for key, item in value.items():
            safe_key = key if isinstance(key, (str, int, bool)) or key is None else repr(key)
            result[safe_key] = _serialise(item)
        return result
    raw_interval = getattr(value, "_mpi_", _MISSING)
    if raw_interval is not _MISSING:
        lower, upper = _interval_endpoints(value, "metadata interval")
        return {
            "lower_exact_mpf_tuple": lower,
            "upper_exact_mpf_tuple": upper,
        }
    raw_mpf = getattr(value, "_mpf_", _MISSING)
    if raw_mpf is not _MISSING:
        return _endpoint_tuple(value, "metadata mpf")
    if isinstance(value, (list, tuple)):
        return [_serialise(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [_serialise(item) for item in sorted(value, key=repr)]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    # Provenance is opaque to this leaf.  A repr is preferable to returning a
    # non-serializable producer object and does not participate in sign logic.
    return repr(value)


def _sequence_metadata(value: Any, label: str) -> list[Any]:
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a sequence")
    return [_serialise(item) for item in value]


def _field(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any:
    return _lookup(mapping, *names)


def _flag(term: Mapping[Any, Any], origin: Mapping[Any, Any], names: Sequence[str], label: str) -> bool | None:
    values: list[bool] = []
    for mapping in (term, origin):
        raw = _field(mapping, names)
        if raw is _MISSING:
            continue
        if not isinstance(raw, bool):
            raise TypeError(f"{label} must be boolean")
        values.append(raw)
    if values and any(value != values[0] for value in values[1:]):
        raise ValueError(f"{label} has conflicting term/origin values")
    return values[0] if values else None


def _metadata(term: Mapping[Any, Any], label: str, branch: str, order: int, cell: int) -> dict[str, Any]:
    origin_raw = term.get("origin", _MISSING)
    if origin_raw is _MISSING:
        origin: Mapping[Any, Any] = {}
    elif isinstance(origin_raw, Mapping):
        origin = origin_raw
    else:
        raise TypeError(f"{label}.origin must be a provenance mapping")

    powers_raw = _field(term, ("powers", "formal_powers", "power_tuple"))
    if powers_raw is _MISSING:
        powers_raw = _field(origin, ("powers", "formal_powers", "power_tuple"))
    if powers_raw is _MISSING or not _is_sequence(powers_raw) or len(powers_raw) != 4:
        raise ValueError(f"{label}.powers must contain four integers")
    powers = tuple(_typed_int(value, f"{label}.powers[{index}]") for index, value in enumerate(powers_raw))

    target_raw = _field(term, ("target_cell", "cell"))
    if target_raw is _MISSING:
        target_raw = _field(origin, ("target_cell", "cell"))
    target = cell if target_raw is _MISSING else _int_cell(target_raw, f"{label}.target_cell")
    if target != cell:
        raise ValueError(f"{label}.target_cell does not match cell {cell}")

    lead_raw = _field(term, ("lead_order", "lead"))
    if lead_raw is _MISSING:
        lead_raw = _field(origin, ("lead_order", "lead"))
    word_raw = _field(term, ("word_order", "word"))
    if word_raw is _MISSING:
        word_raw = _field(origin, ("word_order", "word"))
    if lead_raw is _MISSING or word_raw is _MISSING:
        raise ValueError(f"{label} must retain lead_order and word_order")
    lead = _typed_int(lead_raw, f"{label}.lead_order")
    word = _typed_int(word_raw, f"{label}.word_order")
    if lead not in (0, 1) or word not in (0, 1) or lead + word != order:
        raise ValueError(f"{label} has incompatible lead_order/word_order for Z{order}")

    z_raw = _field(term, ("ordinary_Z_order", "z_order", "ordinary_z_order"))
    if z_raw is _MISSING:
        z_raw = _field(origin, ("ordinary_Z_order", "z_order", "ordinary_z_order"))
    if z_raw is not _MISSING and _int_order(z_raw, f"{label}.ordinary_Z_order") != order:
        raise ValueError(f"{label}.ordinary_Z_order does not match row")

    classes_raw = _field(term, ("error_classes", "error_class", "uncertainty_class"))
    if classes_raw is _MISSING:
        classes_raw = _field(origin, ("error_classes", "error_class", "uncertainty_class"))
    if classes_raw is _MISSING:
        classes: list[Any] = []
    elif isinstance(classes_raw, (str, bytes, bytearray)):
        classes = [_serialise(classes_raw)]
    else:
        classes = _sequence_metadata(classes_raw, f"{label}.error_classes")

    uncertainty_raw = _field(term, ("uncertainty", "is_uncertainty"))
    if uncertainty_raw is _MISSING:
        uncertainty_raw = _field(origin, ("uncertainty", "is_uncertainty"))
    uncertainty = False if uncertainty_raw is _MISSING else uncertainty_raw
    if not isinstance(uncertainty, bool):
        raise TypeError(f"{label}.uncertainty must be boolean")

    error_only = _flag(term, origin, ("error_only",), f"{label}.error_only")
    symmetric = _flag(term, origin, ("symmetric_error", "symmetric_tail_error"), f"{label}.symmetric_error")
    signed_field = _flag(term, origin, ("signed_field", "signed"), f"{label}.signed_field")

    path = _field(term, ("path", "source_path", "source_paths", "input_path", "input_paths", "receipt_path"))
    if path is _MISSING:
        path = _field(origin, ("path", "source_path", "source_paths", "input_path", "input_paths", "receipt_path"))
    path_value = None if path is _MISSING else _serialise(path)

    coefficient = _field(
        term,
        (
            "coefficient",
            "coefficient_interval",
            "interval",
            "enclosure",
            "box",
            "value",
            "exact_interval",
            "exact_endpoints",
            "coefficient_endpoints",
        ),
    )
    if coefficient is _MISSING:
        # A saved receipt may put the exact endpoint pair directly on the
        # term instead of wrapping it under ``coefficient``.
        direct_pair = _pair_from_mapping(term, f"{label}.coefficient")
        if direct_pair is not None:
            coefficient = direct_pair
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient interval")
    endpoint_info = classify_interval(coefficient, f"{label}.coefficient")

    excluded: list[str] = []
    if error_only is True:
        excluded.append("error_only")
    if symmetric is True:
        excluded.append("symmetric_error")
    if signed_field is False:
        excluded.append("signed_field_false")
    is_signed = not excluded

    normalized: dict[Any, Any] = {}
    for key, value in term.items():
        if str(key).lower() in {
            "coefficient",
            "coefficient_interval",
            "interval",
            "enclosure",
            "box",
            "value",
        }:
            normalized[key] = {
                "lower_exact_mpf_tuple": endpoint_info["lower_exact_mpf_tuple"],
                "upper_exact_mpf_tuple": endpoint_info["upper_exact_mpf_tuple"],
            }
        else:
            normalized[key] = _serialise(value)
    normalized.update(
        {
            "branch": branch,
            "ordinary_Z_order": order,
            "target_cell": cell,
            "lead_order": lead,
            "word_order": word,
            "powers": powers,
            "origin": _serialise(origin),
            "path": path_value,
            "source_path": path_value,
            "error_classes": classes,
            "uncertainty": uncertainty,
            "error_only": bool(error_only),
            "symmetric_error": bool(symmetric),
            "signed_field": signed_field,
            "classification": endpoint_info["classification"],
            "sign_class": endpoint_info["classification"],
            "lower_exact_mpf_tuple": endpoint_info["lower_exact_mpf_tuple"],
            "upper_exact_mpf_tuple": endpoint_info["upper_exact_mpf_tuple"],
            "endpoint_tuple": endpoint_info["exact_endpoint_tuple"],
            "excluded_from_signed_counts": not is_signed,
            "excluded_reason": tuple(excluded),
        }
    )
    return {
        "branch": branch,
        "ordinary_Z_order": order,
        "target_cell": cell,
        "lead_order": lead,
        "word_order": word,
        "powers": powers,
        "origin": _serialise(origin),
        "path": path_value,
        "source_path": path_value,
        "error_classes": classes,
        "uncertainty": uncertainty,
        "error_only": bool(error_only),
        "symmetric_error": bool(symmetric),
        "signed_field": signed_field,
        "is_signed": is_signed,
        "excluded_reason": tuple(excluded),
        "classification": endpoint_info["classification"],
        "lower_exact_mpf_tuple": endpoint_info["lower_exact_mpf_tuple"],
        "upper_exact_mpf_tuple": endpoint_info["upper_exact_mpf_tuple"],
        "record": normalized,
    }


def _mapping_orders(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        nested = _lookup(value, "ordinary_Z", "ordinary_z", "orders", "by_order")
        if nested is not _MISSING and nested is not value:
            value = nested
        rows: dict[int, Any] = {}
        for raw_order, row in value.items():
            order = _int_order(raw_order, f"{label} key")
            if order in rows:
                raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
            rows[order] = row
        if set(rows) != set(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label} must provide ordinary-Z orders 0 and 1")
        return rows
    if _is_sequence(value) and len(value) == len(ORDINARY_Z_ORDERS):
        return {order: value[order] for order in ORDINARY_Z_ORDERS}
    raise TypeError(f"{label} must provide ordinary-Z orders 0 and 1")


def _mapping_cells(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        nested = _lookup(value, "target_cells", "cells", "by_cell")
        if nested is not _MISSING and nested is not value:
            value = nested
        cells: dict[int, Any] = {}
        for raw_cell, terms in value.items():
            cell = _int_cell(raw_cell, f"{label} key")
            if cell in cells:
                raise ValueError(f"{label} contains duplicate target cell {cell}")
            cells[cell] = terms
        if set(cells) != set(TARGET_CELLS):
            raise ValueError(f"{label} must provide target cells 0 through 7")
        return cells
    if _is_sequence(value) and len(value) == len(TARGET_CELLS):
        return {cell: value[cell] for cell in TARGET_CELLS}
    raise TypeError(f"{label} must provide target cells 0 through 7")


def _branch_rows(branch: Any, label: str) -> dict[int, dict[int, Sequence[Any]]]:
    if not isinstance(branch, Mapping):
        raise TypeError(f"{label} must be an ordinary-Z/cell mapping")
    orders = _mapping_orders(branch, label)
    result: dict[int, dict[int, Sequence[Any]]] = {}
    for order in ORDINARY_Z_ORDERS:
        cells = _mapping_cells(orders[order], f"{label}.Z{order}")
        result[order] = {}
        for cell in TARGET_CELLS:
            terms = cells[cell]
            if not _is_sequence(terms):
                raise TypeError(f"{label}.Z{order}.cell{cell} must be a term list")
            result[order][cell] = terms
    return result


def _branches(packet: Any) -> tuple[Mapping[Any, Any], dict[str, dict[int, dict[int, Sequence[Any]]]]]:
    if not isinstance(packet, Mapping):
        raise TypeError("branches must be a branch mapping")
    original = packet
    nested = _lookup(packet, "branches")
    if isinstance(nested, Mapping):
        packet = nested
    selected: dict[str, Any] = {}
    for branch in BRANCHES:
        selected[branch] = _map_key(packet, branch, "branches")
    if len(packet) != len(BRANCHES) or {str(key) for key in packet} != set(BRANCHES):
        raise ValueError("branches must contain exactly the four finite P branches")
    rows = {
        branch: _branch_rows(selected[branch], f"branches[{branch!r}]")
        for branch in BRANCHES
    }
    return original, rows


def _obligation(value: Any) -> dict[str, Any]:
    if value is None:
        return {
            "provided": False,
            "authenticated_branches": (),
            "upstream_source_authenticated": None,
            "same_interval_context": None,
            "condition": (
                "sign classes are conditional on the caller's prior-P receipt "
                "authentication and interval-context assumptions"
            ),
        }
    if not isinstance(value, Mapping):
        raise TypeError("caller_obligation must be a mapping when supplied")
    raw_branches = _field(
        value,
        (
            "authenticated_branches",
            "authenticated_finite_P_branches",
            "finite_P_branches",
            "branch_roster",
        ),
    )
    if isinstance(raw_branches, Mapping):
        selected = [branch for branch in BRANCHES if raw_branches.get(branch) is True]
    elif raw_branches is None or raw_branches is _MISSING:
        selected = list(BRANCHES)
    elif _is_sequence(raw_branches):
        selected = list(raw_branches)
    else:
        raise TypeError("caller_obligation authenticated branches must be a sequence or mapping")
    if set(selected) != set(BRANCHES):
        raise ValueError("caller_obligation must authenticate all four finite P branches")
    auth = _field(value, ("upstream_source_authenticated", "source_authenticated", "authenticated"))
    if auth is not _MISSING and auth is not None and auth is not True:
        raise ValueError("caller_obligation source authentication must be true")
    same = _field(value, ("same_interval_context", "same_mpmath_interval_context", "interval_context_same"))
    if same is not _MISSING and same is not None and same is not True:
        raise ValueError("caller_obligation must state one common interval context")
    return {
        "provided": True,
        "authenticated_branches": tuple(selected),
        "upstream_source_authenticated": None if auth is _MISSING else auth,
        "same_interval_context": None if same is _MISSING else same,
        "condition": _serialise(value.get(
            "condition",
            "caller-authenticated prior P physical receipt and one common exact interval context",
        )),
        "record": _serialise(value),
    }


def _counts() -> dict[str, int]:
    return {name: 0 for name in SIGN_CATEGORIES}


def _group(branch: str, order: int, cell: int, lead: int, word: int) -> dict[str, Any]:
    counts = _counts()
    return {
        "branch": branch,
        "ordinary_Z_order": order,
        "target_cell": cell,
        "lead_order": lead,
        "word_order": word,
        "counts": counts,
        "category_counts": dict(counts),
        "records": [],
        "signed_records": [],
        "symmetric_tail_errors": [],
        "error_only_records": [],
        "unsigned_records": [],
        "record_count": 0,
        "signed_record_count": 0,
        "symmetric_tail_error_count": 0,
        "error_only_count": 0,
        "unsigned_record_count": 0,
        "coefficient_arithmetic_performed": False,
        "same_power_terms_merged": False,
        "interval_sum_performed": False,
    }


def _finalise_group(group: dict[str, Any]) -> dict[str, Any]:
    for category in SIGN_CATEGORIES:
        group[category] = group["counts"][category]
    return group


def _cell_result(
    branch: str,
    order: int,
    cell: int,
    terms: Sequence[Any],
) -> dict[str, Any]:
    groups: dict[tuple[int, int], dict[str, Any]] = {}
    all_records: list[dict[str, Any]] = []
    signed_records: list[dict[str, Any]] = []
    symmetric_errors: list[dict[str, Any]] = []
    unsigned_records: list[dict[str, Any]] = []

    for index, term in enumerate(terms):
        label = f"{branch}.Z{order}.cell{cell}[{index}]"
        if not isinstance(term, Mapping):
            raise TypeError(f"{label} must be a term mapping")
        info = _metadata(term, label, branch, order, cell)
        key = (info["lead_order"], info["word_order"])
        group = groups.setdefault(key, _group(branch, order, cell, *key))
        record = dict(info["record"])
        record.update(
            {
                "input_index": index,
                "index": index,
                "is_signed": info["is_signed"],
                "excluded_from_signed_counts": not info["is_signed"],
                "excluded_reason": info["excluded_reason"],
            }
        )
        all_records.append(record)
        group["records"].append(record)
        group["record_count"] += 1
        if info["is_signed"]:
            signed_records.append(record)
            group["signed_records"].append(record)
            group["signed_record_count"] += 1
            category = info["classification"]
            group["counts"][category] += 1
            group["category_counts"][category] += 1
        elif info["error_only"] or info["symmetric_error"]:
            symmetric_errors.append(record)
            group["symmetric_tail_errors"].append(record)
            group["symmetric_tail_error_count"] += 1
            if info["error_only"]:
                group["error_only_records"].append(record)
                group["error_only_count"] += 1
        else:
            unsigned_records.append(record)
            group["unsigned_records"].append(record)
            group["unsigned_record_count"] += 1

    final_groups = [_finalise_group(groups[key]) for key in sorted(groups)]
    by_lead_word: dict[int, dict[int, dict[str, Any]]] = {}
    counts_by_lead_word: dict[int, dict[int, dict[str, int]]] = {}
    for group in final_groups:
        lead, word = group["lead_order"], group["word_order"]
        by_lead_word.setdefault(lead, {})[word] = group
        counts_by_lead_word.setdefault(lead, {})[word] = dict(group["counts"])
    return {
        "branch": branch,
        "ordinary_Z_order": order,
        "target_cell": cell,
        "records": all_records,
        "raw_records": all_records,
        "input_records": all_records,
        "groups": final_groups,
        "aggregates": final_groups,
        "by_lead_word": by_lead_word,
        "counts": counts_by_lead_word,
        "category_counts": counts_by_lead_word,
        "signed_records": signed_records,
        "symmetric_tail_errors": symmetric_errors,
        "error_only_records": [record for record in symmetric_errors if record["error_only"]],
        "unsigned_records": unsigned_records,
        "input_record_count": len(all_records),
        "signed_record_count": len(signed_records),
        "symmetric_tail_error_count": len(symmetric_errors),
        "unsigned_record_count": len(unsigned_records),
        "records_are_endpoint_serialized": True,
        "coefficient_arithmetic_performed": False,
        "same_power_terms_merged": False,
        "interval_sum_performed": False,
        "midpoint_used": False,
        "cap_used": False,
    }


def build(
    context: Any = None,
    branches: Any = None,
    *,
    caller_obligation: Any = None,
    obligation: Any = None,
    caller_condition: Any = None,
    source_certificate: Any = None,
    prior_receipt: Any = None,
) -> dict[str, Any]:
    """Classify caller-supplied finite-P interval records.

    ``context`` is accepted for API compatibility with the prior cell leaf,
    but is never called.  It may be omitted when every coefficient already
    carries exact endpoint tuples.  ``branches`` must contain exactly the
    four finite-P branches, both ordinary-Z rows, and all eight cells.
    """

    if not _is_context(context) and _is_context(branches):
        context, branches = branches, context
    if branches is None and isinstance(context, Mapping):
        try:
            if all(_lookup(context, branch) is not _MISSING for branch in BRANCHES) or isinstance(_lookup(context, "branches"), Mapping):
                branches = context
                context = None
        except (TypeError, ValueError):
            pass
    if branches is None:
        raise TypeError("branches are required")
    if context is not None and not _is_context(context):
        raise TypeError("context must be the caller's interval context or None")

    for alias, value in (
        ("obligation", obligation),
        ("caller_condition", caller_condition),
        ("source_certificate", source_certificate),
        ("prior_receipt", prior_receipt),
    ):
        if value is None:
            continue
        if caller_obligation is not None and caller_obligation is not value:
            raise ValueError(f"caller_obligation and {alias} disagree")
        caller_obligation = value

    _input_packet, rows = _branches(branches)
    obligation_record = _obligation(caller_obligation)

    branch_outputs: dict[str, Any] = {}
    all_counts: dict[str, Any] = {}
    all_records: dict[str, Any] = {}
    all_symmetric: dict[str, Any] = {}
    all_error_only: dict[str, Any] = {}
    all_unsigned: dict[str, Any] = {}
    total_input = 0
    total_signed = 0
    total_symmetric = 0
    total_unsigned = 0

    for branch in BRANCHES:
        order_outputs: dict[int, Any] = {}
        branch_counts: dict[int, Any] = {}
        branch_records: dict[int, Any] = {}
        branch_symmetric: dict[int, Any] = {}
        branch_error_only: dict[int, Any] = {}
        branch_unsigned: dict[int, Any] = {}
        for order in ORDINARY_Z_ORDERS:
            cell_outputs: dict[int, Any] = {}
            order_counts: dict[int, Any] = {}
            order_records: dict[int, Any] = {}
            order_symmetric: dict[int, Any] = {}
            order_error_only: dict[int, Any] = {}
            order_unsigned: dict[int, Any] = {}
            for cell in TARGET_CELLS:
                result = _cell_result(branch, order, cell, rows[branch][order][cell])
                cell_outputs[cell] = result
                order_counts[cell] = result["counts"]
                order_records[cell] = result["records"]
                order_symmetric[cell] = result["symmetric_tail_errors"]
                order_error_only[cell] = result["error_only_records"]
                order_unsigned[cell] = result["unsigned_records"]
                total_input += result["input_record_count"]
                total_signed += result["signed_record_count"]
                total_symmetric += result["symmetric_tail_error_count"]
                total_unsigned += result["unsigned_record_count"]
            order_outputs[order] = {
                "ordinary_Z_order": order,
                "cells": cell_outputs,
                "target_cells": cell_outputs,
                "counts": order_counts,
                "category_counts": order_counts,
                "records": order_records,
                "raw_records": order_records,
                "symmetric_tail_errors": order_symmetric,
                "error_only_records": order_error_only,
                "unsigned_records": order_unsigned,
            }
            branch_counts[order] = order_counts
            branch_records[order] = order_records
            branch_symmetric[order] = order_symmetric
            branch_error_only[order] = order_error_only
            branch_unsigned[order] = order_unsigned
        branch_outputs[branch] = {
            "branch": branch,
            "ordinary_Z": order_outputs,
            "orders": order_outputs,
            "by_order": order_outputs,
            "counts": branch_counts,
            "category_counts": branch_counts,
            "records": branch_records,
            "raw_records": branch_records,
            "symmetric_tail_errors": branch_symmetric,
            "error_only_records": branch_error_only,
            "unsigned_records": branch_unsigned,
            "branch_input_shape": "ordinary_Z -> target_cell -> raw term records",
        }
        all_counts[branch] = branch_counts
        all_records[branch] = branch_records
        all_symmetric[branch] = branch_symmetric
        all_error_only[branch] = branch_error_only
        all_unsigned[branch] = branch_unsigned

    return {
        "schema": SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "branch_names": list(BRANCHES),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "target_cells": list(TARGET_CELLS),
        "sign_categories": list(SIGN_CATEGORIES),
        "branches": branch_outputs,
        "counts": all_counts,
        "category_counts": all_counts,
        "sign_counts": all_counts,
        "records": all_records,
        "raw_records": all_records,
        "symmetric_tail_errors": all_symmetric,
        "error_only_records": all_error_only,
        "unsigned_records": all_unsigned,
        "input_packet_shape": "4 branches x Z0/Z1 x 8 cells x raw term records",
        "input_packet_retained": False,
        "input_endpoint_records_retained": True,
        "input_record_count": total_input,
        "signed_record_count": total_signed,
        "symmetric_tail_error_count": total_symmetric,
        "unsigned_record_count": total_unsigned,
        "caller_obligation": obligation_record.get("record", None),
        "caller_obligation_record": obligation_record,
        "authenticated_branches_condition": obligation_record["authenticated_branches"],
        "upstream_source_authenticated": obligation_record["upstream_source_authenticated"],
        "same_interval_context": obligation_record["same_interval_context"],
        "sign_of_enclosure_only": True,
        "sign_of_enclosure_scope": (
            "conditional on the caller-authenticated prior-P physical receipt, "
            "its signed-field contract, and the supplied directed interval endpoints"
        ),
        "field_sign_claimed": False,
        "physical_field_claimed": False,
        "physical_measure_reapplied": False,
        "source_owner_proved": False,
        "same_source_claimed": False,
        "same_gauge_claimed": False,
        "terminal_closure_proved": False,
        "ledger_appended": False,
        "prior_receipt_reexecuted": False,
        "same_cell_interval_reexecuted": False,
        "coefficient_arithmetic_performed": False,
        "interval_addition_performed": False,
        "interval_subtraction_performed": False,
        "interval_multiplication_performed": False,
        "coefficient_selected": False,
        "same_power_terms_merged": False,
        "integer_count_aggregation_only": True,
        "raw_endpoint_tuples_used": True,
        "exact_binary_endpoint_comparison": True,
        "nonfinite_endpoints_rejected": True,
        "reversed_intervals_rejected": True,
        "midpoint_used": False,
        "cap_used": False,
        "error_only_as_signed": False,
        "signed_field_false_as_signed": False,
        "symmetric_tail_errors_separate": True,
        "records_endpoint_serialized": True,
        "context_supplied_but_not_called": context is not None,
        "input": {
            "branch_names": list(BRANCHES),
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "target_cells": list(TARGET_CELLS),
        },
    }


signs = build
count_signs = build
build_sign_report = build
aggregate = build
combine = build
