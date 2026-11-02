"""Independent finite-K0 profile for the first-word ordinary ``Z2`` packet.

The input is the already saved ``feedback_Z2`` packet from the accepted
``n1_core_word1_z2_v2`` result.  That packet has one term list for each of
the rows ``4``, ``5`` and ``6`` and already contains all six Leibniz
branches.  This leaf only performs the radial Green integration.  It does
not rebuild the source, differentiate a matrix, call an ancestor, choose a
point from an interval, or materialize ``F0``/``Lambda``.

The normalized radius is ``t=x/2`` and the eight canonical cells are
``[i/8,(i+1)/8]``.  A target cell ``i`` receives complete source cells
``j<i`` and the moving-endpoint triangular part of ``j=i``.  The first
diagonal cell uses the analytic axis continuation; no ``log(T/0)`` or
division by the axis endpoint is attempted.  Every product remains a
separate term, so this is a finite K0 condition package, not a complete
first-word field or a moment closure.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "core_firstword_z2_profile_v1"
INPUT_SCHEMA = "n1_core_word1_z2_v2"
CELL_COUNT = 8
ROWS = ("4", "5", "6")
FIELDS = (
    ("F1", "5"),
    ("Uz1", "6"),
    ("K1", "6"),
    ("P1", "4"),
    ("d_x_F1", "5"),
    ("d_x_Uz1", "6"),
)
FIELD_ROW = dict(FIELDS)
ORDERS = ("Z2",)
B1_SEAM_KEYS = ("5,1", "5,2", "5,3")
STATE_COLUMN = {"1": "F1", "2": "Uz1", "3": "K1"}
SEAM_FIELDS = ("F1", "d_x_F1")

MATERIALIZED_KEYS = {
    "F0",
    "Lambda",
    "F0_value",
    "Lambda_value",
    "formal_atoms",
    "atom_values",
}


__all__ = [
    "SCHEMA",
    "INPUT_SCHEMA",
    "CELL_COUNT",
    "ROWS",
    "FIELDS",
    "FIELD_ROW",
    "ORDERS",
    "B1_SEAM_KEYS",
    "STATE_COLUMN",
    "SEAM_FIELDS",
    "cell_kernel_masses",
    "target_radius_kernel_masses",
    "green_cell_kernel_masses",
    "kernel_masses",
    "build",
    "build_missing_seam",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _lookup(mapping: Mapping[Any, Any], key: Any, label: str = "mapping") -> Any:
    if key in mapping:
        return mapping[key]
    alternate = int(key) if isinstance(key, str) and key.isdigit() else str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(f"{label} is missing key {key!r}")


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label} must be a four-entry MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} must contain exact integers")
    sign, mantissa, exponent, bits = (int(item) for item in value)
    if sign not in (0, 1) or mantissa < 0 or bits < 0:
        raise ValueError(f"{label} is not a finite MPF tuple")
    if mantissa == 0 and bits != 0:
        raise ValueError(f"{label} has a noncanonical zero")
    return sign, mantissa, exponent, bits


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None and len(raw) == 2:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None and len(nested) == 2:
            return nested[1 if upper else 0]
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return nested
        return endpoint
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    if hasattr(c, "make_mpf"):
        return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))
    return c.mpf([_raw_endpoint(lower, False), _raw_endpoint(upper, True)])


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode a directed interval without replacing it by a representative."""

    if isinstance(value, float):
        raise TypeError(f"{label} may not contain a Python float")
    if getattr(value, "_mpi_", None) is not None:
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return _make_interval(c, value.a, value.b)
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            if not has_lower or not has_upper:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_tuple(value["lower_exact_mpf_tuple"], label + ".lower")
            upper = _exact_tuple(value["upper_exact_mpf_tuple"], label + ".upper")
            result = _make_interval(c, lower, upper)
        elif "lower" in value or "upper" in value or "lo" in value or "hi" in value:
            lower = value.get("lower", value.get("lo", value.get("left")))
            upper = value.get("upper", value.get("hi", value.get("right")))
            if lower is None or upper is None:
                raise TypeError(f"{label} requires lower and upper endpoints")
            result = _make_interval(c, lower, upper)
        else:
            nested = value.get("coefficient", value.get("directed_interval"))
            if nested is None or nested is value:
                raise TypeError(f"{label} is not interval-compatible")
            return _interval(c, nested, label + ".coefficient")
    elif _is_sequence(value) and len(value) == 2:
        result = _make_interval(c, value[0], value[1])
    else:
        if isinstance(value, bool):
            raise TypeError(f"{label} may not be boolean")
        try:
            result = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
    if result.a > result.b:
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _endpoint_bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    return _make_interval(c, value._mpi_[0], value._mpi_[0]), _make_interval(
        c, value._mpi_[1], value._mpi_[1]
    )


def _exact_cell(c: Any, index: int) -> Any:
    left = c.mpf(index) / CELL_COUNT
    right = c.mpf(index + 1) / CELL_COUNT
    return _make_interval(c, left, right)


def _degenerate_endpoint(c: Any, value: Any, side: int) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[side]
    return _make_interval(c, raw, raw)


def _positive(c: Any, value: Any, label: str) -> None:
    lo, hi = _endpoint_bounds(c, value)
    if hi <= 0:
        raise ValueError(f"{label} must be positive")
    if lo <= 0:
        raise ValueError(f"{label} interval contains the axis")


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value, label)
    lower, upper = value._mpi_
    if _make_interval(c, upper, upper) < 0:
        raise ValueError(f"{label} has a negative upper bound")
    if _make_interval(c, lower, lower) < 0:
        lower = c.mpf(0)._mpi_[0]
    return _make_interval(c, lower, upper)


def _primitive_uz(c: Any, target_t: Any, source_t: Any) -> Any:
    """Primitive of ``2*s^2*log(T/s)+s^2`` in normalized radius."""

    source_t = _interval(c, source_t)
    target_t = _interval(c, target_t)
    source_lower, _ = _endpoint_bounds(c, source_t)
    if source_lower == 0:
        return c.mpf(0)
    _positive(c, source_t, "positive source endpoint")
    _positive(c, target_t, "positive target radius")
    return 2 * source_t * source_t * c.log(target_t / source_t) + source_t * source_t


def cell_kernel_masses(
    c: Any,
    target_t: Any,
    source_left_t: Any,
    source_right_t: Any,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Return directed K0 Green masses for one normalized source cell.

    The formulas are the exact cell integrals for
    ``F1=s(1-s^2/x^2)/2``, ``Uz1=s*log(x/s)``, ``K1=-F1``, ``P1=1``,
    ``d_x_F1=s^3/x^3`` and ``d_x_Uz1=s/x`` after ``x=2*t``.  For a
    diagonal cell ``triangular=True`` replaces the source upper endpoint by
    the moving target interval.  At the axis the removable polynomial limits
    are used explicitly.
    """

    T = _interval(c, target_t, "target_t")
    a = _interval(c, source_left_t, "source_left_t")
    b = T if triangular else _interval(c, source_right_t, "source_right_t")
    Tlo, Thi = _endpoint_bounds(c, T)
    alo, ahi = _endpoint_bounds(c, a)
    blo, bhi = _endpoint_bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo:
        raise ValueError("invalid nonnegative source interval")
    if not triangular and bhi > Thi:
        raise ValueError("source cell extends beyond target radius")
    if triangular and not _same_interval(b, T):
        raise ValueError("triangular source endpoint must be the target interval")

    # This covers the complete first diagonal target cell [0,1/8].
    if triangular and alo == 0 and Tlo == 0:
        f1 = T * T / 2
        uz1 = T * T
        p1 = 2 * T
        dxf1 = T / 2
        dxuz1 = T
    else:
        _positive(c, T, "target radius")
        T2 = T * T
        a2 = a * a
        b2 = b * b
        a4 = a2 * a2
        b4 = b2 * b2
        f1 = (b2 - a2) - (b4 - a4) / (2 * T2)
        uz1 = _primitive_uz(c, T, b) - _primitive_uz(c, T, a)
        p1 = 2 * (b - a)
        dxf1 = (b4 - a4) / (2 * T * T2)
        dxuz1 = (b2 - a2) / T

    f1_nonnegative = _nonnegative(c, f1, "F1 Green mass")
    return {
        "F1": f1_nonnegative,
        "Uz1": _nonnegative(c, uz1, "Uz1 Green mass"),
        "K1": -f1_nonnegative,
        "P1": _nonnegative(c, p1, "P1 Green mass"),
        "d_x_F1": _nonnegative(c, dxf1, "d_x_F1 Green mass"),
        "d_x_Uz1": _nonnegative(c, dxuz1, "d_x_Uz1 Green mass"),
    }


def target_radius_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Physical-radius spelling of :func:`cell_kernel_masses`."""

    return cell_kernel_masses(
        c,
        _interval(c, target_x, "target_x") / 2,
        _interval(c, source_left_x, "source_left_x") / 2,
        _interval(c, source_right_x, "source_right_x") / 2,
        triangular=triangular,
    )


green_cell_kernel_masses = cell_kernel_masses
kernel_masses = cell_kernel_masses


def _formal(term: Mapping[str, Any], label: str) -> tuple[int, int]:
    if "F0_power" not in term or "Lambda_power" not in term:
        raise KeyError(f"{label} requires F0_power and Lambda_power")
    f0 = term["F0_power"]
    lam = term["Lambda_power"]
    if isinstance(f0, bool) or not isinstance(f0, int):
        raise TypeError(f"{label}.F0_power must be an integer")
    if isinstance(lam, bool) or not isinstance(lam, int):
        raise TypeError(f"{label}.Lambda_power must be an integer")
    if any(key in term for key in MATERIALIZED_KEYS):
        raise ValueError(f"{label} may not materialize F0 or Lambda")
    if not any(
        key in term
        for key in (
            "coefficient",
            "coefficient_interval",
            "coefficient_enclosure",
            "directed_interval",
            "enclosure",
            "interval",
        )
    ):
        raise KeyError(f"{label} requires a directed coefficient")
    return int(f0), int(lam)


def _coefficient(c: Any, term: Mapping[str, Any], label: str) -> tuple[Any, Any]:
    for key in (
        "coefficient",
        "coefficient_interval",
        "coefficient_enclosure",
        "directed_interval",
        "enclosure",
        "interval",
    ):
        if key in term:
            raw = term[key]
            return _interval(c, raw, label + ".coefficient"), raw
    raise KeyError(f"{label} requires a coefficient")


def _items(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        return [value]
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a term list or mapping")
    result = list(value)
    if any(not isinstance(item, Mapping) for item in result):
        raise TypeError(f"{label} contains a non-mapping term")
    return result


def _collect_error_classes(term: Mapping[str, Any]) -> list[Any]:
    result: list[Any] = []

    def collect(value: Any) -> None:
        if isinstance(value, str):
            if value not in result:
                result.append(value)
        elif _is_sequence(value):
            for nested in value:
                collect(nested)

    for key in (
        "error_classes",
        "source_error_classes",
        "state_error_classes",
        "matrix_error_classes",
        "unknown_error_classes",
    ):
        if key in term:
            collect(term[key])
    for wrapper in ("provenance", "metadata", "source", "state", "matrix", "input"):
        nested = term.get(wrapper)
        if isinstance(nested, Mapping):
            for key in (
                "error_classes",
                "source_error_classes",
                "state_error_classes",
                "matrix_error_classes",
            ):
                if key in nested:
                    collect(nested[key])
    return result


def _ids(term: Mapping[str, Any]) -> list[Any]:
    result: list[Any] = []

    def add(value: Any) -> None:
        if isinstance(value, str):
            if value not in result:
                result.append(value)
        elif _is_sequence(value):
            for nested in value:
                add(nested)

    for key in (
        "input_id",
        "source_id",
        "term_id",
        "node_id",
        "source_cell_id",
        "matrix_source_id",
        "W0_source_id",
        "input_ids",
        "source_ids",
    ):
        if key in term:
            add(term[key])
    for wrapper in ("provenance", "metadata", "source", "input"):
        nested = term.get(wrapper)
        if isinstance(nested, Mapping):
            for key in ("input_id", "source_id", "term_id", "node_id", "source_ids"):
                if key in nested:
                    add(nested[key])
    return result


def _as_id_list(value: Any, fallback: list[Any]) -> list[Any]:
    if isinstance(value, str):
        return [value]
    if _is_sequence(value):
        return list(value)
    return list(fallback)


def _cell_t(c: Any, cell: Mapping[str, Any], index: int) -> Any:
    if "t" not in cell:
        raise KeyError(f"source cell {index} requires t")
    t = _interval(c, cell["t"], f"source cell {index}.t")
    expected = _exact_cell(c, index)
    if not _same_interval(t, expected):
        raise ValueError(f"source cell {index} has noncanonical t interval")
    if "rho" not in cell:
        raise KeyError(f"source cell {index} requires rho")
    rho = _interval(c, cell["rho"], f"source cell {index}.rho")
    if not _same_interval(rho, 4 * t * t):
        raise ValueError(f"source cell {index} has inconsistent rho=4*t^2")
    return t


def _feedback_terms(
    c: Any, cell: Mapping[str, Any], index: int
) -> dict[str, list[Mapping[str, Any]]]:
    if "feedback_Z2" not in cell:
        raise KeyError(f"source cell {index} is missing feedback_Z2")
    feedback = cell["feedback_Z2"]
    if not isinstance(feedback, Mapping):
        raise TypeError(f"source cell {index}.feedback_Z2 must be a mapping")
    # A few saved wrappers put the order label around the row map.  The
    # accepted v2 receipt uses the direct row map, and both forms are read
    # without reconstructing any branch.
    if "Z2" in feedback and isinstance(feedback["Z2"], Mapping):
        feedback = feedback["Z2"]
    rows: dict[str, list[Mapping[str, Any]]] = {}
    for row in ROWS:
        values = _items(_lookup(feedback, row, f"source cell {index}.feedback_Z2"), row)
        if not values:
            raise ValueError(f"source cell {index} feedback row {row} is empty")
        for term_index, term in enumerate(values):
            _formal(term, f"source cell {index} row {row} term {term_index}")
            _coefficient(c, term, f"source cell {index} row {row} term {term_index}")
        rows[row] = values
    return rows


def _scaled_term(
    c: Any,
    term: Mapping[str, Any],
    mass: Any,
    *,
    source_index: int,
    target_index: int,
    field: str,
    row: str,
    source_term_index: int,
    source_t: Any,
    target_t: Any,
    source_left_t: Any,
    source_right_t: Any,
    triangular: bool,
) -> dict[str, Any]:
    f0, lam = _formal(term, "feedback term")
    coefficient, raw_coefficient = _coefficient(c, term, "feedback term")
    identifiers = _ids(term)
    input_id = identifiers[0] if identifiers else (
        f"feedback_Z2:source:{source_index}:row:{row}:term:{source_term_index}"
    )
    errors = _collect_error_classes(term)
    result = dict(term)
    # Keep the original directed object alongside the scaled coefficient;
    # this is useful when auditing a saved receipt without another full term
    # copy.
    result["input_coefficient"] = raw_coefficient
    result["coefficient"] = coefficient * mass
    result["coefficient_is_directed"] = True
    result["F0_power"] = f0
    result["Lambda_power"] = lam
    result.setdefault("error_classes", list(errors))
    result.setdefault("source_error_classes", list(errors))
    result["input_id"] = term.get("input_id", input_id)
    result["input_ids"] = _as_id_list(
        term.get("input_ids"), identifiers or [input_id]
    )
    result["input_source_cell"] = term.get(
        "source_cell", term.get("source_cell_index", term.get("W0_source_cell"))
    )
    result["source_cell"] = source_index
    result["source_cell_index"] = source_index
    result["target_cell"] = target_index
    result["target_cell_index"] = target_index
    result["source_term_index"] = source_term_index
    if "term_id" in term:
        result["input_term_id"] = term["term_id"]
    result["term_id"] = (
        f"firstword_Z2:target:{target_index}:source:{source_index}:{field}:{source_term_index}"
    )
    result["source_term_id"] = input_id
    result["source_row"] = row
    result["target_field"] = field
    result["ordinary_Z_order"] = 2
    result["W1_order"] = "Z2"
    result["kernel_mass"] = mass
    result["kernel_mass_field"] = field
    result["kernel_mass_exact"] = True
    result["target_t"] = target_t
    result["source_t"] = source_t
    result["source_left_t"] = source_left_t
    result["source_right_t"] = source_right_t
    result["triangular_diagonal"] = triangular
    result["no_cancellation"] = True
    result["formal_atoms_materialized"] = False
    result["F0_materialized"] = False
    result["Lambda_materialized"] = False
    return result


def _resolve_cells(cells: Any) -> list[Mapping[str, Any]]:
    if isinstance(cells, Mapping):
        for key in ("cells", "target_cells", "outer_cells", "saved_outer_cells"):
            if key in cells:
                cells = cells[key]
                break
        else:
            raise KeyError("cells requires an eight-cell sequence")
    if not _is_sequence(cells):
        raise TypeError("cells must be an eight-cell sequence")
    result = list(cells)
    if len(result) != CELL_COUNT:
        raise ValueError("exactly eight radial source cells are required")
    if any(not isinstance(cell, Mapping) for cell in result):
        raise TypeError("radial source cells must be mappings")
    return result


def build(c: Any, cells: Any) -> dict[str, Any]:
    """Build eight target-cell ``W1_Z2`` Green condition packages.

    ``cells`` is the saved ``result.outer_cells`` sequence.  Each cell must
    carry canonical ``t``/``rho`` and ``feedback_Z2`` rows.  The result keeps
    one canonical term table per target cell under ``target_cells[*].W1_Z2``;
    envelope-level ``terms`` views are compact source-cell/term-ID references.
    """

    source_cells = _resolve_cells(cells)
    source_t = [_cell_t(c, cell, index) for index, cell in enumerate(source_cells)]
    feedback = [
        _feedback_terms(c, cell, index)
        for index, cell in enumerate(source_cells)
    ]

    target_cells: list[dict[str, Any]] = []
    # ``target_cells[*].W1_Z2`` is the one canonical full term table.  All
    # envelope-level ``terms``/by-source views below contain references only;
    # this avoids serializing the same large directed intervals repeatedly.
    target_terms: dict[str, dict[str, list[dict[str, Any]]]] = {}
    target_by_source: dict[str, dict[str, dict[str, list[dict[str, Any]]]]] = {}
    target_kernels: dict[str, dict[str, dict[str, Any]]] = {}

    for target_index in range(CELL_COUNT):
        target_t = source_t[target_index]
        aggregate = {field: [] for field, _ in FIELDS}
        aggregate_refs = {field: [] for field, _ in FIELDS}
        by_source: dict[str, dict[str, list[dict[str, Any]]]] = {}
        kernels_by_source: dict[str, dict[str, Any]] = {}
        for source_index in range(target_index + 1):
            # Strict causality: complete source cells use source_index <
            # target_index.  The only source with source_index == target_index
            # is the moving-upper-endpoint triangular diagonal.
            diagonal = source_index == target_index
            source_cell_t = source_t[source_index]
            source_left_t = _degenerate_endpoint(c, source_cell_t, 0)
            source_right_t = (
                target_t
                if diagonal
                else _degenerate_endpoint(c, source_cell_t, 1)
            )
            masses = cell_kernel_masses(
                c,
                target_t,
                source_left_t,
                source_right_t,
                triangular=diagonal,
            )
            kernels_by_source[str(source_index)] = masses
            source_refs = {field: [] for field, _ in FIELDS}
            for field, row in FIELDS:
                for term_index, term in enumerate(feedback[source_index][row]):
                    emitted = _scaled_term(
                        c,
                        term,
                        masses[field],
                        source_index=source_index,
                        target_index=target_index,
                        field=field,
                        row=row,
                        source_term_index=term_index,
                        source_t=source_cell_t,
                        target_t=target_t,
                        source_left_t=source_left_t,
                        source_right_t=source_right_t,
                        triangular=diagonal,
                    )
                    source_refs[field].append(
                        {
                            "term_id": emitted["term_id"],
                            "input_id": emitted.get("input_id"),
                            "source_cell": source_index,
                            "source_term_index": term_index,
                            "target_cell": target_index,
                            "field": field,
                        }
                    )
                    aggregate_refs[field].append(source_refs[field][-1])
                    aggregate[field].append(emitted)
            by_source[str(source_index)] = source_refs

        target_terms[str(target_index)] = aggregate_refs
        target_by_source[str(target_index)] = by_source
        target_kernels[str(target_index)] = kernels_by_source
        target_cells.append(
            {
                "index": target_index,
                "t": target_t,
                "x": 2 * target_t,
                "rho": 4 * target_t * target_t,
                "source_indices": list(range(target_index + 1)),
                "source_j_lt_target_i_for_full_cells": True,
                "source_j_le_target_i": True,
                "source_j_gt_target_i_excluded": True,
                "strict_causality": True,
                "triangular_diagonal_source": target_index,
                "axis_endpoint_continuation": target_index == 0,
                "terms": {"Z2": aggregate_refs},
                "terms_Z2": aggregate_refs,
                "W1_Z2": aggregate,
                "W1_by_order": {"Z2": {"canonical": "W1_Z2", "term_refs": aggregate_refs}},
                "W1_by_order_by_source_cell": {"Z2": by_source},
                "terms_by_source_cell": {"Z2": by_source},
                "kernel_masses": kernels_by_source,
                "kernel_mass_scope": "K0 source-cell integral; diagonal upper endpoint moves with target t",
                "input_feedback": "feedback_Z2",
            }
        )

    return {
        "schema": SCHEMA,
        "input_schema": INPUT_SCHEMA,
        "target_cells": target_cells,
        "outer_cells": [
            {"index": cell["index"], "target_cell_ref": cell["index"]}
            for cell in target_cells
        ],
        "terms": {"Z2": target_terms},
        "terms_by_order": {"Z2": target_terms},
        "terms_by_target_cell": target_terms,
        "W1_Z2_target_cells": target_terms,
        "W1_Z2_target_cells_by_source": target_by_source,
        "kernel_masses_by_target_source": target_kernels,
        "saved_source_cells": [
            {
                "index": index,
                "source_cell_ref": index,
                "input_feedback": "feedback_Z2",
            }
            for index in range(CELL_COUNT)
        ],
        "saved_outer_cells_reused": True,
        "feedback_Z2_reused": True,
        "feedback_Z2_reconstructed": False,
        "all_six_leibniz_branches_input": True,
        "all_six_leibniz_branches_reconstructed": False,
        "no_source_cells_collapsed": True,
        "source_indices_are_causal_prefixes": True,
        "strict_j_less_than_i_for_complete_cells": True,
        "diagonal_source_is_triangular": True,
        "axis_endpoint_continuation": True,
        "axis_log_division_avoided": True,
        "source_ids_retained": True,
        "provenance_retained": True,
        "error_classes_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "no_cancellation": True,
        "cancellation_performed": False,
        "formal_atoms_materialized": False,
        "F0_materialized": False,
        "Lambda_materialized": False,
        "source_exponentials_materialized": False,
        "green_kernel_mass_definition": "exact directed source-cell integral; diagonal source cell is triangular",
        "green_kernel_masses_exact": True,
        "principal_kernel": "K0 only; no Kchi or higher kernel words",
        "ordinary_Z_orders": [2],
        "input_feedback_order": 2,
        "output_Z_order": 2,
        "coordinate": "rho=4*t^2, normalized t=x/2, real Z in [-1,1]",
        "normalized_radius": "t=x/2",
        "finite_K0_condition_package": True,
        "complete_firstword_Z2": False,
        "full_firstword_claimed": False,
        "physical_field_claimed": False,
        "global_five_moment_closure": False,
        "first_collar_added": False,
        "allword_tail_added": False,
        "m15_computed": False,
    }


def _resolve_profile_cells(profile: Any) -> tuple[list[Mapping[str, Any]], Mapping[str, Any] | None]:
    envelope: Mapping[str, Any] | None = profile if isinstance(profile, Mapping) else None
    if isinstance(profile, Mapping):
        values = profile.get("target_cells", profile.get("outer_cells"))
        if values is None:
            raise KeyError("profile requires target_cells")
    elif _is_sequence(profile):
        values = profile
    else:
        raise TypeError("profile must be a mapping or target-cell sequence")
    values = list(values)
    if len(values) != CELL_COUNT or any(not isinstance(cell, Mapping) for cell in values):
        raise ValueError("exactly eight firstword target cells are required")
    return values, envelope


def _resolve_matrix_cells(
    profile_cells: list[Mapping[str, Any]],
    matrix_cells: Any,
    envelope: Mapping[str, Any] | None,
) -> list[Mapping[str, Any]]:
    if matrix_cells is None and envelope is not None:
        matrix_cells = envelope.get("matrix_cells", envelope.get("raw_matrix_cells"))
    if matrix_cells is None:
        extracted: list[Mapping[str, Any]] = []
        for cell in profile_cells:
            nested = cell.get("matrix_cell", cell.get("raw_matrix"))
            if nested is None and "B1" in cell:
                nested = cell
            if not isinstance(nested, Mapping):
                extracted = []
                break
            extracted.append(nested)
        matrix_cells = extracted or None
    if matrix_cells is None or not _is_sequence(matrix_cells):
        raise KeyError("raw B1 matrix_cells are required")
    result = list(matrix_cells)
    if len(result) != CELL_COUNT or any(not isinstance(cell, Mapping) for cell in result):
        raise ValueError("exactly eight raw matrix cells are required")
    return result


def _matrix_terms(cell: Mapping[str, Any], order: str, key: str) -> list[Mapping[str, Any]]:
    try:
        b1 = _lookup(cell, "B1", "matrix cell")
        ordered = _lookup(b1, order, "B1")
        return _items(_lookup(ordered, key, f"B1[{order}]"), f"B1[{order}][{key}]")
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell is missing B1[{order}][{key}]") from exc


def _state_terms(cell: Mapping[str, Any], field: str) -> list[Mapping[str, Any]]:
    # ``W1_Z2`` is the canonical full term table emitted by this leaf;
    # envelope aliases deliberately contain only term references.
    family: Any = cell.get("W1_Z2")
    if family is None:
        family = cell.get("W1_by_order")
    if family is not None and isinstance(family, Mapping):
        if "Z2" in family:
            family = family["Z2"]
        elif "2" in family:
            family = family["2"]
    if family is None:
        family = cell.get("terms", cell.get("terms_Z2", cell.get("W1_Z2")))
        if isinstance(family, Mapping) and "Z2" in family:
            family = family["Z2"]
    if not isinstance(family, Mapping):
        raise KeyError(f"target cell is missing W1 Z2 terms for {field}")
    return _items(_lookup(family, field, "W1 Z2"), f"W1 Z2 {field}")


def _seam_term(
    c: Any,
    matrix_term: Mapping[str, Any],
    state_term: Mapping[str, Any],
    mass: Any,
    *,
    source_index: int,
    source_t: Any,
    matrix_key: str,
    matrix_term_index: int,
    state_field: str,
    state_term_index: int,
    target_field: str,
    target_x: Any,
    target_t: Any,
) -> dict[str, Any]:
    matrix_f0, matrix_lam = _formal(matrix_term, "B1[0] matrix term")
    state_f0, state_lam = _formal(state_term, "firstword Z2 state term")
    matrix_coefficient, matrix_raw = _coefficient(c, matrix_term, "B1[0] matrix term")
    state_coefficient, state_raw = _coefficient(c, state_term, "firstword Z2 state term")
    matrix_ids = _ids(matrix_term)
    state_ids = _ids(state_term)
    matrix_id = matrix_ids[0] if matrix_ids else f"matrix:B1:cell:{source_index}:0:{matrix_key}:{matrix_term_index}"
    state_id = state_ids[0] if state_ids else f"state:W1:cell:{source_index}:Z2:{state_field}:{state_term_index}"
    matrix_ids = matrix_ids or [matrix_id]
    state_ids = state_ids or [state_id]
    classes = _collect_error_classes(matrix_term)
    for item in _collect_error_classes(state_term):
        if item not in classes:
            classes.append(item)
    result = {
        "F0_power": matrix_f0 + state_f0,
        "Lambda_power": matrix_lam + state_lam,
        "coefficient": matrix_coefficient * state_coefficient * mass,
        "coefficient_is_directed": True,
        "ordinary_Z_order": 1,
        "origin": "B1[0]*W1_Z2 (missing Z1 branch)",
        "missing_product_rule": "B1[0]*W1_Z2",
        "matrix_name": "B1",
        "matrix_order": "0",
        "matrix_key": matrix_key,
        "row": 5,
        "column": int(matrix_key.split(",")[1]),
        "state_field": state_field,
        "state_order": "2",
        "target_field": target_field,
        "matrix_term_index": matrix_term_index,
        "state_term_index": state_term_index,
        "source_cell": source_index,
        "source_cell_index": source_index,
        "source_radius_cell": source_index,
        "target_cell": "x=2 seam",
        "target_cell_index": -1,
        "target_radius_cell": -1,
        "source_t": source_t,
        "target_x": target_x,
        "target_t": target_t,
        "kernel_mass": mass,
        "kernel_mass_field": target_field,
        "kernel_mass_exact": True,
        "kernel_mass_positive": target_field != "K1",
        "matrix_input_id": matrix_id,
        "state_input_id": state_id,
        "state_profile_term_id": state_term.get("term_id"),
        "matrix_receipt_input_path": matrix_term.get("receipt_input_path"),
        "state_receipt_input_path": state_term.get("receipt_input_path"),
        "input_ids": matrix_ids + [item for item in state_ids if item not in matrix_ids],
        "source_ids": matrix_ids + [item for item in state_ids if item not in matrix_ids],
        "input_origins": {
            "matrix": matrix_term.get("origin"),
            "state": state_term.get("origin"),
        },
        "matrix_source_cell": matrix_term.get("source_cell", source_index),
        "matrix_source_cell_index": matrix_term.get("source_cell_index", source_index),
        "state_source_cell": state_term.get("source_cell", source_index),
        "state_source_cell_index": state_term.get("source_cell_index", source_index),
        "matrix_identity_opaque": True,
        "state_identity_opaque": True,
        "error_classes": classes,
        "source_error_classes": _collect_error_classes(matrix_term),
        "state_error_classes": _collect_error_classes(state_term),
        "uncertainty": bool(classes),
        "signed": True,
        "signed_field": True,
        "input_matrix_coefficient": matrix_raw,
        "input_state_coefficient": state_raw,
        "matrix_provenance": matrix_term.get("provenance", matrix_term.get("origin")),
        "state_provenance": state_term.get("provenance", state_term.get("origin")),
        "term_id": f"firstword_Z2_missing_seam:cell:{source_index}:Z1:{target_field}:{matrix_key}:{matrix_term_index}:{state_term_index}",
        "no_cancellation": True,
        "cancellation_performed": False,
        "formal_atoms_materialized": False,
        "F0_materialized": False,
        "Lambda_materialized": False,
    }
    return result


def build_missing_seam(c: Any, profile: Any, matrix_cells: Any = None) -> dict[str, Any]:
    """Emit the finite missing ``B1[0]*W1_Z2`` seam branch.

    The branch is the ordinary-Z1 product omitted by the saved Z0/Z1
    firstword profile.  Columns ``1,2,3`` of each raw ``B1[0]`` cell are
    paired with the corresponding ``F1``, ``Uz1`` and ``K1`` state terms in
    that same target/source radial cell.  Each pair is integrated separately
    to the fixed ``x=2`` seam with exact directed K0 masses (including both
    ``F1`` and ``d_x_F1`` / ``Fx`` scales).  No prior term is changed,
    cancelled, or folded into this result.
    """

    # The frozen one-entry consumer may probe an optional seam method with
    # the eight raw B1 cells alone.  Such a packet deliberately contains no
    # firstword Z2 state, so fail closed with an explicit non-materialized
    # record rather than inventing state terms.  Supplying ``profile`` and
    # ``matrix_cells`` together takes the actual product path below.
    if (
        matrix_cells is None
        and _is_sequence(profile)
        and len(profile) == CELL_COUNT
        and all(
            isinstance(cell, Mapping)
            and "B1" in cell
            and not any(key in cell for key in ("W1_Z2", "terms_Z2", "W1_by_order", "terms"))
            for cell in profile
        )
    ):
        return {
            "schema": SCHEMA + ":missing_seam_v1",
            "input_schema": SCHEMA,
            "status": "not_materialized",
            "reason": "firstword Z2 profile cells are required alongside raw B1 seam cells",
            "missing_product_rule": "B1[0]*W1_Z2",
            "missing_branch_included": False,
            "missing_branch_materialized": False,
            "finite_K0_condition_package": True,
            "complete_secondword": False,
            "physical_field_claimed": False,
            "formal_powers_retained": True,
            "signed_coefficients_retained": True,
            "provenance_retained": True,
            "error_classes_retained": True,
            "no_cancellation": True,
            "cancellation_performed": False,
            "formal_atoms_materialized": False,
            "F0_materialized": False,
            "Lambda_materialized": False,
        }

    profile_cells, envelope = _resolve_profile_cells(profile)
    matrices = _resolve_matrix_cells(profile_cells, matrix_cells, envelope)
    target_x = c.mpf(2)
    target_t = c.mpf(1)
    output = {field: [] for field in SEAM_FIELDS}
    by_source: dict[str, dict[str, list[dict[str, Any]]]] = {}
    kernels: dict[str, dict[str, Any]] = {}
    source_records: list[dict[str, Any]] = []

    for source_index, (target_cell, matrix_cell) in enumerate(zip(profile_cells, matrices)):
        source_cell_t = _cell_t(c, target_cell, source_index)
        matrix_t = _cell_t(c, matrix_cell, source_index)
        if not _same_interval(source_cell_t, matrix_t):
            raise ValueError(f"target and matrix cell {source_index} have different t intervals")
        source_left_t = _degenerate_endpoint(c, source_cell_t, 0)
        source_right_t = _degenerate_endpoint(c, source_cell_t, 1)
        masses = cell_kernel_masses(
            c,
            target_t,
            source_left_t,
            source_right_t,
            triangular=False,
        )
        kernels[str(source_index)] = {
            "source_cell": source_index,
            "source_cell_index": source_index,
            "t": source_cell_t,
            "a": source_left_t,
            "b": source_right_t,
            "target_x": target_x,
            "target_t": target_t,
            "F1": masses["F1"],
            "Fx": masses["d_x_F1"],
            "d_x_F1": masses["d_x_F1"],
            "exact": True,
            "positive": True,
            "independent_polynomial": True,
        }
        cell_refs = {field: [] for field in SEAM_FIELDS}
        term_ids: list[str] = []
        for matrix_key in B1_SEAM_KEYS:
            column = matrix_key.split(",")[1]
            state_field = STATE_COLUMN[column]
            matrix_values = _matrix_terms(matrix_cell, "0", matrix_key)
            state_values = _state_terms(target_cell, state_field)
            for matrix_term_index, matrix_term in enumerate(matrix_values):
                _formal(matrix_term, "B1[0] matrix term")
                _coefficient(c, matrix_term, "B1[0] matrix term")
                for state_term_index, state_term in enumerate(state_values):
                    _formal(state_term, "firstword Z2 state term")
                    _coefficient(c, state_term, "firstword Z2 state term")
                    for target_field in SEAM_FIELDS:
                        emitted = _seam_term(
                            c,
                            matrix_term,
                            state_term,
                            masses[target_field],
                            source_index=source_index,
                            source_t=source_cell_t,
                            matrix_key=matrix_key,
                            matrix_term_index=matrix_term_index,
                            state_field=state_field,
                            state_term_index=state_term_index,
                            target_field=target_field,
                            target_x=target_x,
                            target_t=target_t,
                        )
                        output[target_field].append(emitted)
                        cell_refs[target_field].append(
                            {
                                "term_id": emitted["term_id"],
                                "matrix_input_id": emitted.get("matrix_input_id"),
                                "state_input_id": emitted.get("state_input_id"),
                                "source_cell": source_index,
                                "target_field": target_field,
                            }
                        )
                        term_ids.append(emitted["term_id"])
        by_source[str(source_index)] = cell_refs
        source_records.append(
            {
                "index": source_index,
                "source_cell": source_index,
                "source_cell_index": source_index,
                "t": source_cell_t,
                "rho": _interval(c, target_cell["rho"], f"target cell {source_index}.rho"),
                "matrix_source_cell": source_index,
                "state_source_cell": source_index,
                "columns": [1, 2, 3],
                "term_ids": term_ids,
                "term_counts": {field: len(values) for field, values in cell_refs.items()},
            }
        )

    return {
        "schema": SCHEMA + ":missing_seam_v1",
        "input_schema": SCHEMA,
        "source_cells": source_records,
        "target_seam": {
            "x": target_x,
            "t": target_t,
            "rho": c.mpf(4),
            "fixed": True,
            "target_cell_index": -1,
            "kernel_mass_scope": "source cell [i/8,(i+1)/8] at x=2",
        },
        "kernels": kernels,
        "terms_by_order": {"Z1": output},
        "missing_Z1": {
            field: [
                {
                    "term_id": term["term_id"],
                    "source_cell": term.get("source_cell"),
                    "target_field": field,
                }
                for term in values
            ]
            for field, values in output.items()
        },
        "terms": {"Z1": {
            field: [
                {
                    "term_id": term["term_id"],
                    "source_cell": term.get("source_cell"),
                    "target_field": field,
                }
                for term in values
            ]
            for field, values in output.items()
        }},
        "by_source_cell": {"Z1": by_source},
        "missing_product_rule": "B1[0]*W1_Z2",
        "missing_branch_included": True,
        "known_Z1_records_untouched": True,
        "old_records_merged": False,
        "cancellation_performed": False,
        "no_cancellation": True,
        "B1_zero_columns": [1, 2, 3],
        "state_columns": {"1": "F1", "2": "Uz1", "3": "K1"},
        "output_fields": list(SEAM_FIELDS),
        "ordinary_Z_orders": [1],
        "state_Z_order": 2,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "provenance_retained": True,
        "error_classes_retained": True,
        "source_ids_retained": True,
        "matrix_sector_identity_proved": False,
        "state_matrix_same_product_proved": False,
        "kernel_mass_exact": True,
        "kernel_mass_positive": True,
        "F_and_Fx_exact": True,
        "finite_K0_condition_package": True,
        "principal_kernel": "K0 only; missing seam branch only",
        "formal_atoms_materialized": False,
        "F0_materialized": False,
        "Lambda_materialized": False,
        "complete_secondword": False,
        "physical_field_claimed": False,
        "m15_computed": False,
    }

