"""Finite first-word residual at the native ``rho=4`` inlet seam.

This is a source-only, pure-function leaf.  It does not open files, import a
saved producer, evaluate a source, or add a tail estimate.  The caller owns
the conversion of saved records to interval-valued terms and is expected to
have checked the common source family, common ``Z`` tube, anchor, and scale.

The canonical input data format is deliberately small and can be built by a
caller without importing this module's historical producers.  A sector term
is a mapping with exactly these semantic fields (additional provenance is
allowed)::

    {
        "source_exponents": (e0, e1, e2, e3),
        "combined_source_log_interval": {"lower": L0, "upper": L1},
        "coefficient_interval": {"lower": c0, "upper": c1},
    }

The four exponent entries are integers.  The two intervals may instead be a
two-item sequence or an :class:`Interval` instance.  Endpoint values are
never converted to a midpoint or to ``float``.  Signed coefficient intervals
are multiplied by the four-corner interval rule and are negated when a raw
residual contribution is formed.  Certified callers must supply exact
endpoints or directed interval endpoint objects; ordinary floating-point
scalars alone do not provide outward rounding.

Matrices use one-based field indices and ordinary-``Z`` jets::

    {
        "Z0": {(row, column): [sector_term, ...]},
        "Z1": {(row, column): [sector_term, ...]},
    }

String keys such as ``"4,1"`` are accepted for JSON-shaped in-memory
records.  A vector uses the analogous column map::

    {
        "Z0": {column: [sector_term, ...]},
        "Z1": {column: [sector_term, ...]},
    }

The primary function, :func:`compute_inlet_word_residual`, takes ``B0``,
``B1`` and the saved homogeneous ``W1``.  It removes exactly the kinematic
entries ``(1,5)``, ``(2,6)``, and ``(3,6)`` from ``B0`` to form ``C0`` and
assembles only dynamics rows ``4, 5, 6``.  For ``S1 = W0 + W1`` the finite
word identity is

``L_raw S1 - g = -(C0 W1 + B1 * d_Z W1)``.

The returned ``Z0`` sectors are fully assembled.  The returned
``Z1_known`` sectors contain only
``-(C0_Z W1 + C0 W1_Z + B1_Z W1_Z)``.  Since the input carries only ``W1_Z0``
and ``W1_Z1``, the missing ``-(B1 W1_ZZ)`` is returned separately, one
symbolic term per supplied ``B1`` sector.  Consequently no output claims a
fully evaluated first ``Z`` derivative.

Terms with equal source-exponent tuples are grouped with interval addition
and multiplication.  If their source-log intervals differ, the result keeps
each original log/coeff pair in ``source_log_terms`` and exposes a hull in
``combined_source_log_interval``.  The hull is taken from the optional
caller-supplied ``source_log_hulls`` mapping when present; otherwise it is a
conservative hull of the supplied intervals.  The grouped coefficient is a
signed interval sum, but ``collapsed_as_single_exponential`` is false for a
multi-log group, so consumers cannot mistake the diagnostic hull for a
midpoint repair or erase source-log correlation.

The omitted infinite-word tail is intentionally not part of this residual.
The result records that it remains independent and unpropagated.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias


COMPONENT = "n1_inlet_word_residual"
SCHEMA = "n1_inlet_word_residual_v1"
INPUT_SCHEMA = "n1_inlet_word_residual_input_v1"
RESULT_SCHEMA = "n1_inlet_word_residual_result_v1"

DYNAMIC_ROWS = (4, 5, 6)
KINEMATIC_ENTRIES = ((1, 5), (2, 6), (3, 6))
ORDINARY_Z_ORDERS = (0, 1)
SOURCE_DIMENSION = 4

Number: TypeAlias = Any
ExponentTuple: TypeAlias = tuple[int, int, int, int]


class ResidualInputError(ValueError):
    """Raised when a caller-supplied sector/matrix contract is malformed."""


@dataclass(frozen=True)
class Interval:
    """A closed interval retaining caller endpoint objects exactly.

    The class is intentionally agnostic about endpoint type.  ``Decimal``,
    ``Fraction``, arbitrary-precision scalars, and caller interval backends
    that expose ordered endpoints can all be retained without a conversion.
    """

    lower: Number
    upper: Number

    def as_mapping(self) -> dict[str, Number]:
        return {"lower": self.lower, "upper": self.upper}


@dataclass(frozen=True)
class SectorTerm:
    """Canonical immutable representation of one source-log sector."""

    source_exponents: ExponentTuple
    combined_source_log_interval: Interval
    coefficient_interval: Interval

    def as_mapping(self) -> dict[str, Any]:
        return {
            "source_exponents": self.source_exponents,
            "combined_source_log_interval": self.combined_source_log_interval.as_mapping(),
            "coefficient_interval": self.coefficient_interval.as_mapping(),
        }


def _fail(message: str) -> None:
    raise ResidualInputError(message)


def _as_interval(value: Any, label: str) -> Interval:
    """Normalize a two-endpoint value without changing endpoint scale."""

    if isinstance(value, Interval):
        interval = value
    elif isinstance(value, Mapping):
        lower = value.get("lower", value.get("lower_bound"))
        upper = value.get("upper", value.get("upper_bound"))
        if lower is None or upper is None:
            # Exact-endpoint records are a useful caller-facing spelling.
            lower = value.get("lower_exact")
            upper = value.get("upper_exact")
        if lower is None or upper is None:
            _fail(f"{label} must contain lower and upper endpoints")
        interval = Interval(lower, upper)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            _fail(f"{label} must be a two-endpoint interval")
        interval = Interval(value[0], value[1])
    else:
        # A scalar is a point interval.  This is useful for exact zero terms
        # and does not create a midpoint of a nontrivial caller interval.
        interval = Interval(value, value)
    try:
        reversed_interval = interval.upper < interval.lower
        if bool(reversed_interval):
            _fail(f"{label} has reversed endpoints")
    except TypeError:
        # Some arbitrary-precision interval endpoint wrappers deliberately do
        # not expose scalar ordering.  Their arithmetic backend remains the
        # caller's responsibility; retaining the endpoints is preferable to
        # coercing them through float or midpoint arithmetic.
        pass
    return interval


def interval_add(left: Any, right: Any) -> Interval:
    """Return the outward endpoint sum ``left + right``."""

    a, b = _as_interval(left, "left interval"), _as_interval(right, "right interval")
    lower, upper = a.lower + b.lower, a.upper + b.upper
    return Interval(getattr(lower, "a", lower), getattr(upper, "b", upper))


def interval_negate(value: Any) -> Interval:
    """Negate an interval while preserving signed coefficient orientation."""

    interval = _as_interval(value, "interval")
    lower, upper = -interval.upper, -interval.lower
    return Interval(getattr(lower, "a", lower), getattr(upper, "b", upper))


def interval_multiply(left: Any, right: Any) -> Interval:
    """Multiply two intervals using all four endpoint products."""

    a, b = _as_interval(left, "left interval"), _as_interval(right, "right interval")
    products = (
        a.lower * b.lower,
        a.lower * b.upper,
        a.upper * b.lower,
        a.upper * b.upper,
    )
    try:
        lower = min(getattr(value, "a", value) for value in products)
        upper = max(getattr(value, "b", value) for value in products)
    except TypeError as exc:
        raise ResidualInputError("interval endpoints must be orderable after multiplication") from exc
    return Interval(lower, upper)


def interval_hull(values: Iterable[Any]) -> Interval:
    """Return the endpoint hull; no midpoint or sign repair is performed."""

    intervals = tuple(_as_interval(value, "hull member") for value in values)
    if not intervals:
        _fail("cannot form an empty interval hull")
    try:
        lower = min(getattr(interval.lower, "a", interval.lower) for interval in intervals)
        upper = max(getattr(interval.upper, "b", interval.upper) for interval in intervals)
    except TypeError as exc:
        raise ResidualInputError("interval endpoints must be orderable for a hull") from exc
    return Interval(lower, upper)


# Short aliases are useful to interval callbacks that use ``add``/``mul``
# vocabulary.  They remain pure wrappers around the explicit operations.
interval_sum = interval_add
interval_product = interval_multiply


def _interval_mapping(value: Interval) -> dict[str, Number]:
    return value.as_mapping()


def _exponents(value: Any, label: str) -> ExponentTuple:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        _fail(f"{label} must be a four-entry exponent tuple")
    if len(value) != SOURCE_DIMENSION or any(type(item) is not int for item in value):
        _fail(f"{label} must contain exactly four integers")
    return tuple(int(item) for item in value)  # type: ignore[return-value]


def _term(value: Any, label: str) -> SectorTerm:
    if isinstance(value, SectorTerm):
        return value
    if not isinstance(value, Mapping):
        _fail(f"{label} must be a sector mapping")
    exponents = _exponents(value.get("source_exponents"), f"{label}.source_exponents")
    log_value = value.get("combined_source_log_interval", value.get("source_log_interval"))
    coeff_value = value.get("coefficient_interval", value.get("coefficient"))
    if log_value is None:
        _fail(f"{label}.combined_source_log_interval is missing")
    if coeff_value is None:
        _fail(f"{label}.coefficient_interval is missing")
    return SectorTerm(
        source_exponents=exponents,
        combined_source_log_interval=_as_interval(log_value, f"{label}.combined_source_log_interval"),
        coefficient_interval=_as_interval(coeff_value, f"{label}.coefficient_interval"),
    )


def _term_output(term: SectorTerm, *, extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    out = term.as_mapping()
    if extra:
        out.update(dict(extra))
    return out


def add_source_exponents(left: Any, right: Any) -> ExponentTuple:
    """Add two four-dimensional source exponent tuples componentwise."""

    a = _exponents(left, "left source_exponents")
    b = _exponents(right, "right source_exponents")
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def multiply_sector_terms(left: Any, right: Any) -> dict[str, Any]:
    """Multiply two source sectors, retaining their log and signed factors."""

    a, b = _term(left, "left sector"), _term(right, "right sector")
    product = SectorTerm(
        source_exponents=add_source_exponents(a.source_exponents, b.source_exponents),
        combined_source_log_interval=interval_add(
            a.combined_source_log_interval,
            b.combined_source_log_interval,
        ),
        coefficient_interval=interval_multiply(a.coefficient_interval, b.coefficient_interval),
    )
    return _term_output(
        product,
        extra={
            "factor_source_exponents": (a.source_exponents, b.source_exponents),
            "factor_source_log_intervals": (
                _interval_mapping(a.combined_source_log_interval),
                _interval_mapping(b.combined_source_log_interval),
            ),
            "factor_coefficient_intervals": (
                _interval_mapping(a.coefficient_interval),
                _interval_mapping(b.coefficient_interval),
            ),
        },
    )


def _negated_sector(value: Any) -> dict[str, Any]:
    term = _term(value, "sector")
    return _term_output(
        SectorTerm(
            source_exponents=term.source_exponents,
            combined_source_log_interval=term.combined_source_log_interval,
            coefficient_interval=interval_negate(term.coefficient_interval),
        ),
        extra={"signed_residual_factor": -1},
    )


def _same_interval(left: Interval, right: Interval) -> bool:
    return left.lower == right.lower and left.upper == right.upper


def _hull_for_exponents(
    exponent: ExponentTuple,
    observed: tuple[Interval, ...],
    source_log_hulls: Mapping[Any, Any] | None,
) -> tuple[Interval, str]:
    """Resolve a log hull and identify whether it came from the caller."""

    if source_log_hulls is not None:
        candidate: Any = None
        if exponent in source_log_hulls:
            candidate = source_log_hulls[exponent]
        else:
            text_key = ",".join(str(item) for item in exponent)
            if text_key in source_log_hulls:
                candidate = source_log_hulls[text_key]
        if candidate is not None:
            hull = _as_interval(candidate, f"source_log_hulls[{exponent!r}]")
            try:
                covered = all(hull.lower <= item.lower and item.upper <= hull.upper for item in observed)
            except TypeError as exc:
                raise ResidualInputError("source-log hull endpoints must be orderable") from exc
            if not covered:
                _fail(f"caller source-log hull does not cover exponent sector {exponent!r}")
            return hull, "caller_hull"
    return interval_hull(observed), "input_hull"


def merge_sector_terms(
    terms: Iterable[Any],
    *,
    source_log_hulls: Mapping[Any, Any] | None = None,
) -> list[dict[str, Any]]:
    """Merge sectors by exponent tuple using interval arithmetic.

    Equal exponent tuples and equal source-log intervals are summed directly.
    For one exponent tuple carrying multiple distinct source logs, all source
    log/coeff pairs are retained in ``source_log_terms`` and a caller hull (or
    input hull) is exposed for inspection.  The result explicitly refuses to
    call that multi-log group a single exponential.
    """

    canonical = tuple(_term(value, f"sector[{index}]") for index, value in enumerate(terms))
    groups: dict[ExponentTuple, list[SectorTerm]] = {}
    for item in canonical:
        groups.setdefault(item.source_exponents, []).append(item)

    merged: list[dict[str, Any]] = []
    for exponent in sorted(groups):
        members = groups[exponent]
        log_groups: list[list[SectorTerm]] = []
        for member in members:
            for group in log_groups:
                if _same_interval(member.combined_source_log_interval, group[0].combined_source_log_interval):
                    group.append(member)
                    break
            else:
                log_groups.append([member])

        if len(log_groups) == 1:
            group = log_groups[0]
            coefficient = group[0].coefficient_interval
            for member in group[1:]:
                coefficient = interval_add(coefficient, member.coefficient_interval)
            merged.append(
                _term_output(
                    SectorTerm(exponent, group[0].combined_source_log_interval, coefficient),
                    extra={
                        "merged_term_count": len(group),
                        "source_log_terms": [
                            _term_output(member) for member in group
                        ],
                        "source_log_correlation_preserved": True,
                        "collapsed_as_single_exponential": True,
                        "source_log_hull_source": "exact_equal_log",
                    },
                )
            )
            continue

        logs = tuple(member.combined_source_log_interval for member in members)
        hull, hull_source = _hull_for_exponents(exponent, logs, source_log_hulls)
        coefficient = members[0].coefficient_interval
        for member in members[1:]:
            coefficient = interval_add(coefficient, member.coefficient_interval)
        merged.append(
            _term_output(
                SectorTerm(exponent, hull, coefficient),
                extra={
                    "merged_term_count": len(members),
                    "source_log_terms": [
                        _term_output(member) for member in members
                    ],
                    "source_log_correlation_preserved": True,
                    "collapsed_as_single_exponential": False,
                    "source_log_hull_source": hull_source,
                    "source_log_mismatch": True,
                },
            )
        )
    return merged


def _pair_key(value: Any, label: str) -> tuple[int, int]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        row, column = value
    elif isinstance(value, str):
        pieces = value.replace(":", ",").replace("/", ",").split(",")
        if len(pieces) != 2:
            _fail(f"{label} must be a row,column key")
        row, column = pieces
    else:
        _fail(f"{label} must be a row,column key")
    if isinstance(row, bool) or isinstance(column, bool):
        _fail(f"{label} row and column must be integers")
    try:
        row_int, col_int = int(row), int(column)
    except (TypeError, ValueError) as exc:
        raise ResidualInputError(f"{label} row and column must be integers") from exc
    if row_int <= 0 or col_int <= 0:
        _fail(f"{label} row and column must be positive one-based indices")
    return row_int, col_int


def _term_list(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, Mapping):
        for key in ("sectors", "sector_terms", "terms"):
            if key in value:
                return _term_list(value[key], f"{label}.{key}")
        return [value]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    _fail(f"{label} must be a sector or a sequence of sectors")


def _order_key(value: Any, label: str) -> int:
    text = str(value).replace("ordinary_Z", "Z")
    if text in {"Z0", "0"}:
        return 0
    if text in {"Z1", "1"}:
        return 1
    _fail(f"{label} must be Z0 or Z1")


def _matrix_order_payload(payload: Any, order: int, label: str) -> dict[tuple[int, int], list[Any]]:
    out: dict[tuple[int, int], list[Any]] = {}
    if payload is None:
        return out
    if isinstance(payload, Mapping):
        # Canonical pair-key map.
        pair_like = True
        for raw_key, value in payload.items():
            try:
                pair = _pair_key(raw_key, f"{label} key")
            except ResidualInputError:
                pair_like = False
                break
            out[pair] = _term_list(value, f"{label}[{raw_key!r}]")
        if pair_like and (payload or not out):
            return out
        # Nested row -> column map is convenient for in-memory callers.
        out.clear()
        for raw_row, columns in payload.items():
            try:
                row = int(raw_row)
            except (TypeError, ValueError) as exc:
                raise ResidualInputError(f"{label} row key is not an integer") from exc
            if not isinstance(columns, Mapping):
                _fail(f"{label}[{raw_row!r}] must map columns to sectors")
            for raw_column, terms in columns.items():
                pair = _pair_key((row, raw_column), f"{label} entry")
                out[pair] = _term_list(terms, f"{label}[{raw_row!r}][{raw_column!r}]")
        return out
    if isinstance(payload, Sequence) and not isinstance(payload, (str, bytes, bytearray)):
        for index, entry in enumerate(payload):
            if not isinstance(entry, Mapping):
                _fail(f"{label}[{index}] must be a mapping")
            pair = _pair_key((entry.get("row"), entry.get("column", entry.get("col"))), f"{label}[{index}]")
            terms = entry.get("sectors", entry.get("sector_terms", entry.get("terms")))
            out.setdefault(pair, []).extend(_term_list(terms, f"{label}[{index}].terms"))
        return out
    _fail(f"{label} must be a matrix order mapping or sequence")


def _normalize_matrix(matrix: Any, label: str) -> dict[int, dict[tuple[int, int], list[Any]]]:
    out = {0: {}, 1: {}}
    if isinstance(matrix, Mapping):
        order_keys = []
        for raw_key, payload in matrix.items():
            try:
                order = _order_key(raw_key, label)
            except ResidualInputError:
                order_keys = []
                break
            order_keys.append(order)
            out[order] = _matrix_order_payload(payload, order, f"{label}.{raw_key}")
        if order_keys:
            return out
        # Pair-key map whose values contain Z0/Z1.
        for raw_pair, value in matrix.items():
            pair = _pair_key(raw_pair, f"{label} key")
            if isinstance(value, Mapping) and any(key in value for key in ("Z0", "Z1", "ordinary_Z0", "ordinary_Z1")):
                for raw_order, payload in value.items():
                    try:
                        order = _order_key(raw_order, label)
                    except ResidualInputError:
                        continue
                    out[order].setdefault(pair, []).extend(_term_list(payload, f"{label}[{raw_pair!r}].{raw_order}"))
                continue
            _fail(f"{label}[{raw_pair!r}] must provide Z0 and Z1 terms")
        return out
    if isinstance(matrix, Sequence) and not isinstance(matrix, (str, bytes, bytearray)):
        for index, entry in enumerate(matrix):
            if not isinstance(entry, Mapping):
                _fail(f"{label}[{index}] must be a mapping")
            pair = _pair_key((entry.get("row"), entry.get("column", entry.get("col"))), f"{label}[{index}]")
            for raw_order, payload in entry.items():
                try:
                    order = _order_key(raw_order, label)
                except ResidualInputError:
                    continue
                out[order].setdefault(pair, []).extend(_term_list(payload, f"{label}[{index}].{raw_order}"))
        return out
    _fail(f"{label} must be a Z0/Z1 matrix mapping")


def _vector_order_payload(payload: Any, order: int, label: str) -> dict[int, list[Any]]:
    out: dict[int, list[Any]] = {}
    if payload is None:
        return out
    if isinstance(payload, Mapping):
        if any(key in payload for key in ("sectors", "sector_terms", "terms")):
            column = payload.get("column", payload.get("col"))
            if column is None:
                _fail(f"{label} sector entry has no column")
            out[int(column)] = _term_list(payload, label)
            return out
        for raw_column, terms in payload.items():
            try:
                column = int(raw_column)
            except (TypeError, ValueError) as exc:
                raise ResidualInputError(f"{label} column key is not an integer") from exc
            if column <= 0:
                _fail(f"{label} column must be positive one-based")
            out[column] = _term_list(terms, f"{label}[{raw_column!r}]")
        return out
    if isinstance(payload, Sequence) and not isinstance(payload, (str, bytes, bytearray)):
        for index, entry in enumerate(payload):
            if not isinstance(entry, Mapping):
                _fail(f"{label}[{index}] must be a mapping")
            try:
                column = int(entry.get("column", entry.get("col")))
            except (TypeError, ValueError) as exc:
                raise ResidualInputError(f"{label}[{index}] column is not an integer") from exc
            out[column] = _term_list(entry.get("sectors", entry.get("sector_terms", entry.get("terms"))), f"{label}[{index}]")
        return out
    _fail(f"{label} must be a vector order mapping or sequence")


def _normalize_vector(vector: Any, label: str) -> dict[int, dict[int, list[Any]]]:
    out = {0: {}, 1: {}}
    if isinstance(vector, Mapping):
        orders_present = []
        for raw_key, payload in vector.items():
            try:
                order = _order_key(raw_key, label)
            except ResidualInputError:
                orders_present = []
                break
            orders_present.append(order)
            out[order] = _vector_order_payload(payload, order, f"{label}.{raw_key}")
        if orders_present:
            return out
        # Column-key map whose values contain Z0/Z1.
        for raw_column, value in vector.items():
            try:
                column = int(raw_column)
            except (TypeError, ValueError) as exc:
                raise ResidualInputError(f"{label} column key is not an integer") from exc
            if not isinstance(value, Mapping):
                _fail(f"{label}[{raw_column!r}] must provide Z0 and Z1 terms")
            for raw_order, payload in value.items():
                try:
                    order = _order_key(raw_order, label)
                except ResidualInputError:
                    continue
                out[order][column] = _term_list(payload, f"{label}[{raw_column!r}].{raw_order}")
        return out
    if isinstance(vector, Sequence) and not isinstance(vector, (str, bytes, bytearray)):
        for index, entry in enumerate(vector):
            if not isinstance(entry, Mapping):
                _fail(f"{label}[{index}] must be a mapping")
            try:
                column = int(entry.get("column", entry.get("col")))
            except (TypeError, ValueError) as exc:
                raise ResidualInputError(f"{label}[{index}] column is not an integer") from exc
            for raw_order, payload in entry.items():
                try:
                    order = _order_key(raw_order, label)
                except ResidualInputError:
                    continue
                out[order][column] = _term_list(payload, f"{label}[{index}].{raw_order}")
        return out
    _fail(f"{label} must be a Z0/Z1 vector mapping")


def remove_kinematic_entries(B0: Any) -> dict[str, Any]:
    """Return a copy-shaped normalized matrix with the three kinematic entries removed."""

    matrix = _normalize_matrix(B0, "B0")
    out: dict[str, Any] = {"Z0": {}, "Z1": {}}
    for order, order_name in ((0, "Z0"), (1, "Z1")):
        for pair, terms in matrix[order].items():
            if pair in KINEMATIC_ENTRIES:
                continue
            out[order_name][pair] = list(terms)
    return out


def _products(
    matrix_terms: Iterable[Any],
    vector_terms: Iterable[Any],
    *,
    source_log_hulls: Mapping[Any, Any] | None,
) -> list[dict[str, Any]]:
    products: list[dict[str, Any]] = []
    for matrix_term in matrix_terms:
        for vector_term in vector_terms:
            products.append(_negated_sector(multiply_sector_terms(matrix_term, vector_term)))
    # Group only once, after all columns and contributions have been joined.
    # Re-grouping an already merged multi-log sector would discard its
    # source_log_terms and falsely collapse distinct logarithm enclosures.
    return products


def _unknown_second_Z_terms(
    row: int,
    B1_order0: Mapping[tuple[int, int], list[Any]],
    *,
    source_log_hulls: Mapping[Any, Any] | None,
) -> list[dict[str, Any]]:
    """List ``-B1[row,col] * W1_ZZ[col]`` symbolically, term by term."""

    # The caller does not provide W1_ZZ in this leaf.  Each B1 sector is
    # retained with its source log and signed coefficient; the symbolic factor
    # below prevents it from being interpreted as a computed sector value.
    entries: list[dict[str, Any]] = []
    by_column: dict[int, list[dict[str, Any]]] = {}
    for (matrix_row, column), matrix_terms in B1_order0.items():
        if matrix_row != row:
            continue
        for index, matrix_term in enumerate(matrix_terms):
            negated = _negated_sector(matrix_term)
            by_column.setdefault(column, []).append(
                _term_output(
                    _term(negated, "negated B1 sector"),
                    extra={
                        "row": row,
                        "column": column,
                        "matrix_term_index": index,
                        "symbolic_factor": f"W1_ZZ[{column}]",
                        "symbolic_formula": f"-B1[{row},{column}] * W1_ZZ[{column}]",
                        "cannot_omit": True,
                        "fully_evaluated": False,
                        "source_log_hulls_supplied": source_log_hulls is not None,
                    },
                )
            )
    for column in sorted(by_column):
        entries.extend(by_column[column])
    return entries


def _contribution_bundle(
    row: int,
    *,
    C0: dict[int, dict[tuple[int, int], list[Any]]],
    B1: dict[int, dict[tuple[int, int], list[Any]]],
    W1: dict[int, dict[int, list[Any]]],
    source_log_hulls: Mapping[Any, Any] | None,
) -> dict[str, Any]:
    def matrix_row(matrix: dict[tuple[int, int], list[Any]], column: int) -> list[Any]:
        return matrix.get((row, column), [])

    columns = sorted({column for (matrix_row_index, column) in C0[0] if matrix_row_index == row}
                     | {column for (matrix_row_index, column) in C0[1] if matrix_row_index == row}
                     | {column for (matrix_row_index, column) in B1[0] if matrix_row_index == row}
                     | {column for (matrix_row_index, column) in B1[1] if matrix_row_index == row})
    # A missing W1 column means no saved field is available for a nonzero
    # matrix sector.  Empty lists are valid explicit zero sectors.
    used_columns = [
        column
        for column in columns
        if any(matrix_row(C0[order], column) or matrix_row(B1[order], column) for order in ORDINARY_Z_ORDERS)
    ]
    for column in used_columns:
        if column not in W1[0] or column not in W1[1]:
            _fail(f"W1 column {column} must provide both Z0 and Z1 sector lists")

    c0_w1 = []
    b1_w1z = []
    c0z_w1 = []
    c0_w1z = []
    b1z_w1z = []
    for column in used_columns:
        c0_w1.extend(_products(matrix_row(C0[0], column), W1[0].get(column, []), source_log_hulls=source_log_hulls))
        b1_w1z.extend(_products(matrix_row(B1[0], column), W1[1].get(column, []), source_log_hulls=source_log_hulls))
        c0z_w1.extend(_products(matrix_row(C0[1], column), W1[0].get(column, []), source_log_hulls=source_log_hulls))
        c0_w1z.extend(_products(matrix_row(C0[0], column), W1[1].get(column, []), source_log_hulls=source_log_hulls))
        b1z_w1z.extend(_products(matrix_row(B1[1], column), W1[1].get(column, []), source_log_hulls=source_log_hulls))

    known_z0 = merge_sector_terms(c0_w1 + b1_w1z, source_log_hulls=source_log_hulls)
    known_z1 = merge_sector_terms(c0z_w1 + c0_w1z + b1z_w1z, source_log_hulls=source_log_hulls)
    return {
        "row": row,
        "Z0": {
            "sectors": known_z0,
            "fully_evaluated": True,
            "formula": "-(C0 W1 + B1 W1_Z)",
            "contributions": {
                "-C0_W1": c0_w1,
                "-B1_W1_Z": b1_w1z,
            },
        },
        "Z1_known": {
            "sectors": known_z1,
            "fully_evaluated": False,
            "formula": "-(C0_Z W1 + C0 W1_Z + B1_Z W1_Z)",
            "contributions": {
                "-C0_Z_W1": c0z_w1,
                "-C0_W1_Z": c0_w1z,
                "-B1_Z_W1_Z": b1z_w1z,
            },
        },
        "Z1_omitted_remainder": {
            "sectors": _unknown_second_Z_terms(row, B1[0], source_log_hulls=source_log_hulls),
            "fully_evaluated": False,
            "cannot_omit": True,
            "formula": "-(B1 W1_ZZ)",
            "status": "symbolic_uncomputed_remainder",
        },
    }


def compute_inlet_word_residual(
    *,
    B0: Any,
    B1: Any,
    W1: Any,
    source_log_hulls: Mapping[Any, Any] | None = None,
    rows: Sequence[int] = DYNAMIC_ROWS,
    location: str = "rho=4 core/first seam",
) -> dict[str, Any]:
    """Assemble the finite-word raw residual for dynamics rows ``4,5,6``.

    ``B0``, ``B1`` and ``W1`` are already-normalized interval records; this
    function never calls a source callback.  ``source_log_hulls`` is optional
    because callers may have no sector collision.  When supplied, a hull is
    selected by the exact four-exponent tuple (or its comma-separated string)
    and must contain every observed log interval for that tuple.  No midpoint
    is accepted or constructed.
    """

    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes, bytearray)):
        _fail("rows must be a sequence of one-based dynamic row indices")
    selected_rows = tuple(int(row) for row in rows)
    if selected_rows != DYNAMIC_ROWS:
        _fail("the raw inlet residual is defined for dynamics rows exactly (4, 5, 6)")

    b0 = _normalize_matrix(B0, "B0")
    b1 = _normalize_matrix(B1, "B1")
    c0_public = remove_kinematic_entries(B0)
    c0 = _normalize_matrix(c0_public, "C0")
    w1 = _normalize_vector(W1, "W1")

    row_results = {
        row: _contribution_bundle(
            row,
            C0=c0,
            B1=b1,
            W1=w1,
            source_log_hulls=source_log_hulls,
        )
        for row in DYNAMIC_ROWS
    }
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "input_schema": INPUT_SCHEMA,
        "location": location,
        "finite_word": 1,
        "finite_partial_sum": "S1 = W0 + W1",
        "equation": "L_raw S1 - g = -(C0 W1 + B1 W1_Z)",
        "kinematic_entries_removed_from_B0": [list(pair) for pair in KINEMATIC_ENTRIES],
        "matrices": {
            "C0": c0_public,
            "B1": {
                "Z0": dict(b1[0]),
                "Z1": dict(b1[1]),
            },
        },
        "dynamic_rows": list(DYNAMIC_ROWS),
        "raw_signed_residual": {str(row): row_results[row]["Z0"]["sectors"] for row in DYNAMIC_ROWS},
        "ordinary_Z": {
            "Z0": {str(row): row_results[row]["Z0"] for row in DYNAMIC_ROWS},
            "Z1_known": {str(row): row_results[row]["Z1_known"] for row in DYNAMIC_ROWS},
            "Z1_omitted_remainder": {
                str(row): row_results[row]["Z1_omitted_remainder"] for row in DYNAMIC_ROWS
            },
            "Z1_fully_evaluated": False,
            "Z1_formula_with_uncomputed_term": "-(C0_Z W1 + C0 W1_Z + B1_Z W1_Z + B1 W1_ZZ)",
        },
        "rows": {str(row): row_results[row] for row in DYNAMIC_ROWS},
        "source_log_policy": {
            "source_dimension": SOURCE_DIMENSION,
            "hull_used_instead_of_midpoint": True,
            "caller_hulls_supplied": source_log_hulls is not None,
            "different_source_logs_preserved": True,
            "signed_coefficients_preserved": True,
        },
        "omitted_word_tail": {
            "included_in_residual": False,
            "status": "independent_unpropagated",
            "statement": "The all-word tail is not added to the finite S1 residual.",
        },
        "limitations": {
            "W1_ZZ": "not supplied; represented by explicit non-omittable symbolic remainder",
            "Z1_fully_evaluated": False,
            "source_identity_checked_here": False,
            "same_family_checked_by_caller": True,
            "Z_tube_checked_by_caller": True,
            "anchor_checked_by_caller": True,
            "scale_checked_by_caller": True,
            "interval_callbacks_invoked": False,
            "files_read": False,
            "old_modules_imported": False,
        },
    }


# Descriptive aliases make the leaf easy to discover without introducing a
# second implementation or changing the canonical result contract.
assemble_inlet_word_residual = compute_inlet_word_residual
build_inlet_word_residual = compute_inlet_word_residual
n1_inlet_word_residual = compute_inlet_word_residual


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "INPUT_SCHEMA",
    "RESULT_SCHEMA",
    "DYNAMIC_ROWS",
    "KINEMATIC_ENTRIES",
    "Interval",
    "SectorTerm",
    "ResidualInputError",
    "interval_add",
    "interval_negate",
    "interval_multiply",
    "interval_hull",
    "interval_sum",
    "interval_product",
    "add_source_exponents",
    "multiply_sector_terms",
    "merge_sector_terms",
    "remove_kinematic_entries",
    "compute_inlet_word_residual",
    "assemble_inlet_word_residual",
    "build_inlet_word_residual",
    "n1_inlet_word_residual",
]
