"""Finite raw-K0 second-feedback source leaf for the core row 5 sector.

The input to this module is the already saved first-word profile and the raw
``C0``/``B1`` matrix sectors saved beside its eight target cells.  It is an
input-only adapter: it never calls an ancestor, reconstructs a first-word
state, or materializes ``F0`` or ``Lambda``.

The row-5 columns ``1`` (``F1``), ``2`` (``Uz1``), ``3`` (``K1``), and ``5``
(``d_x_F1``) are admitted wherever the saved raw matrix has that sector.
Thus C0 uses columns ``1,2,3,5`` and B1 uses columns ``1,2,3``.  The known
ordinary-Z rules are

``Z0: C0[0] W[0] + B1[0] W[1]``
``Z1: C0[1] W[0] + C0[0] W[1] + B1[1] W[1]``.

The omitted ``B1[0] W[2]`` term is returned as explicit unknown records for
the available B1 columns.  No absent B1 column 5 sector is invented.  Thus
the result is a finite source subset, not a complete second word or a
physical field.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "core_secondword_partial_v1"
INPUT_SCHEMA = "n1_core_w1_profile_v2"
OUTPUT_FIELDS = ("F1", "d_x_F1")
STATE_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
STATE_COLUMN = {"1": "F1", "2": "Uz1", "3": "K1", "4": "P1", "5": "d_x_F1", "6": "d_x_Uz1"}
USED_COLUMNS = (1, 2, 3, 5)
RAW_C0_COLUMNS = (1, 2, 3, 5)
RAW_B1_COLUMNS = (1, 2, 3)
UNKNOWN_COLUMNS = ()
MATRIX_NAMES = ("C0", "B1")
ORDERS = ("0", "1")
CELL_COUNT = 8
C0_MATRIX_KEYS = ("5,1", "5,2", "5,3", "5,5")
B1_MATRIX_KEYS = ("5,1", "5,2", "5,3")

PRODUCT_RULE_Z0 = (
    ("C0", "0", "0", "C0[0]*W[0]"),
    ("B1", "0", "1", "B1[0]*W[1]"),
)
PRODUCT_RULE_Z1 = (
    ("C0", "1", "0", "C0[1]*W[0]"),
    ("C0", "0", "1", "C0[0]*W[1]"),
    ("B1", "1", "1", "B1[1]*W[1]"),
)
PRODUCT_RULES = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}

_MISSING = object()
_COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "coefficient_enclosure",
    "directed_interval",
    "enclosure",
    "interval",
)

__all__ = [
    "SCHEMA",
    "INPUT_SCHEMA",
    "OUTPUT_FIELDS",
    "STATE_FIELDS",
    "USED_COLUMNS",
    "RAW_C0_COLUMNS",
    "RAW_B1_COLUMNS",
    "C0_MATRIX_KEYS",
    "B1_MATRIX_KEYS",
    "PRODUCT_RULE_Z0",
    "PRODUCT_RULE_Z1",
    "PRODUCT_RULES",
    "target_seam_kernel_masses",
    "cell_kernel_masses",
    "build",
    "build_secondword_partial",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


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


def _raw_endpoints(value: Any) -> tuple[Any, Any]:
    return _raw_endpoint(value, False), _raw_endpoint(value, True)


def _is_interval(value: Any) -> bool:
    return getattr(value, "_mpi_", None) is not None and hasattr(value, "a") and hasattr(value, "b")


def _reject_float(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(f"{label} may not contain a Python float")
    if isinstance(value, Mapping):
        for key, nested in value.items():
            _reject_float(nested, f"{label}.{key}")
    elif _is_sequence(value):
        for index, nested in enumerate(value):
            _reject_float(nested, f"{label}[{index}]")


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4 or any(type(item) is not int for item in value):
        raise TypeError(f"{label} must be a four-entry exact MPF tuple")
    result = tuple(int(item) for item in value)
    sign, mantissa, _, bits = result
    if sign not in (0, 1) or mantissa < 0 or bits < 0 or (mantissa == 0 and bits != 0):
        raise ValueError(f"{label} is not a finite exact MPF tuple")
    return result


def _interval(c: Any, value: Any, label: str) -> Any:
    """Decode an input as a directed interval in the caller's context."""

    _reject_float(value, label)
    if _is_interval(value):
        result = value
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
                return _interval(c, nested, label + ".coefficient")
            lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
            upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
            if lower is _MISSING or upper is _MISSING:
                raise TypeError(f"{label} is not interval-compatible")
            result = c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))
    elif _is_sequence(value) and len(value) == 2:
        result = c.make_mpf((_raw_endpoint(value[0], False), _raw_endpoint(value[1], True)))
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


def _canonical_cell(c: Any, index: int) -> Any:
    left = c.mpf(index) / CELL_COUNT
    right = c.mpf(index + 1) / CELL_COUNT
    return c.make_mpf((_raw_endpoint(left, False), _raw_endpoint(right, True)))


def target_seam_kernel_masses(c: Any, index: int) -> dict[str, Any]:
    """Return exact positive K0 seam masses for ``x=2`` and cell ``index``.

    The endpoints are canonical zero-width rationals.  Consequently these
    masses are independent of any saved target-cell interval and use only the
    two polynomial primitives requested by the core seam:

    ``F = (b**2-a**2) - (b**4-a**4)/2`` and
    ``Fx = (b**4-a**4)/2``.
    """

    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < CELL_COUNT:
        raise ValueError("target seam cell index must be an integer in [0, 7]")
    a = c.mpf(index) / CELL_COUNT
    b = c.mpf(index + 1) / CELL_COUNT
    a2, b2 = a * a, b * b
    a4, b4 = a2 * a2, b2 * b2
    f_mass = (b2 - a2) - (b4 - a4) / 2
    fx_mass = (b4 - a4) / 2
    if f_mass.a < 0 or fx_mass.a < 0 or f_mass.b <= 0 or fx_mass.b <= 0:
        raise ValueError("target seam K0 masses must be positive")
    return {
        "F1": f_mass,
        "d_x_F1": fx_mass,
        "dxF1": fx_mass,
        "F": f_mass,
        "Fx": fx_mass,
        "a": a,
        "b": b,
        "formula": {
            "F1": "(b^2-a^2)-(b^4-a^4)/2",
            "d_x_F1": "(b^4-a^4)/2",
        },
        "positive": True,
        "exact": True,
        "independent_polynomial": True,
    }


cell_kernel_masses = target_seam_kernel_masses


def _lookup(mapping: Mapping[Any, Any], key: Any, label: str = "mapping") -> Any:
    if key in mapping:
        return mapping[key]
    alternate = int(key) if isinstance(key, str) and key.isdigit() else str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(f"{label} is missing key {key!r}")


def _items(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        return [value]
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a term list or mapping")
    result = list(value)
    if any(not isinstance(item, Mapping) for item in result):
        raise TypeError(f"{label} contains a non-mapping term")
    return result


def _formal(term: Mapping[str, Any], label: str) -> tuple[int, int]:
    if "F0_power" not in term or "Lambda_power" not in term:
        raise KeyError(f"{label} requires F0_power and Lambda_power")
    f0 = term["F0_power"]
    lam = term["Lambda_power"]
    if isinstance(f0, bool) or not isinstance(f0, int) or isinstance(lam, bool) or not isinstance(lam, int):
        raise TypeError(f"{label} powers must be integers")
    if not any(key in term for key in _COEFFICIENT_KEYS):
        raise KeyError(f"{label} requires a directed coefficient")
    if any(key in term for key in ("F0", "Lambda", "F0_value", "Lambda_value", "formal_atoms", "atom_values")):
        raise ValueError(f"{label} may not materialize F0 or Lambda")
    return int(f0), int(lam)


def _coefficient(c: Any, term: Mapping[str, Any], label: str) -> Any:
    raw = next((term[key] for key in _COEFFICIENT_KEYS if key in term), _MISSING)
    if raw is _MISSING:
        raise KeyError(f"{label} requires a coefficient")
    return _interval(c, raw, label + ".coefficient")


def _term_classes(term: Mapping[str, Any]) -> list[str]:
    classes: list[str] = []

    def collect(value: Any) -> None:
        if isinstance(value, str):
            if value not in classes:
                classes.append(value)
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
    for wrapper in ("origin", "metadata", "provenance", "source", "state", "matrix"):
        nested = term.get(wrapper)
        if isinstance(nested, Mapping):
            for key in ("error_classes", "source_error_classes", "state_error_classes", "matrix_error_classes"):
                if key in nested:
                    collect(nested[key])
    return classes


def _values(term: Mapping[str, Any], names: tuple[str, ...]) -> list[Any]:
    result: list[Any] = []
    for name in names:
        if name in term:
            value = term[name]
            if isinstance(value, str):
                result.append(value)
            elif _is_sequence(value):
                result.extend(_clone(item) for item in value)
            else:
                result.append(_clone(value))
    return result


def _paths(term: Mapping[str, Any]) -> list[Any]:
    result: list[Any] = []
    for name in (
        "source_path",
        "state_path",
        "matrix_path",
        "path",
        "parent_receipt_path",
        "receipt_input_path",
        "parent_input_path",
        "matrix_sector_identity",
        "W0_sector_identity",
    ):
        if name in term:
            result.append(_clone(term[name]))
    for name in ("source_paths", "state_paths", "matrix_paths", "paths"):
        if name in term:
            value = term[name]
            result.extend(_clone(item) for item in (value if _is_sequence(value) else [value]))
    return result


def _node_ref(term: Mapping[str, Any], *, role: str, cell_index: int, term_index: int, order: str, field: str | None = None) -> dict[str, Any]:
    ids = _values(
        term,
        (
            "node_id",
            "source_id",
            "state_id",
            "term_id",
            "input_id",
            "source_cell_id",
            "matrix_source_id",
            "W0_source_id",
            "source_ids",
            "state_ids",
        ),
    )
    ids = [item for item in ids if item is not None]
    if not ids:
        ids = [f"{role}:cell:{cell_index}:Z{order}:{field or 'matrix'}:{term_index}"]
    paths = _paths(term)
    return {
        "role": role,
        "node_id": next((item for item in ids if isinstance(item, str)), f"{role}:cell:{cell_index}:Z{order}:{field or 'matrix'}:{term_index}"),
        "ids": ids,
        "source_ids": ids,
        "paths": paths,
        "source_paths": paths,
        "error_classes": _term_classes(term),
        "source_error_classes": _term_classes(term),
        "source_cell": term.get("source_cell", cell_index),
        "source_cell_index": term.get("source_cell_index", cell_index),
        "source_cell_id": term.get("source_cell_id", term.get("source_cell", cell_index)),
        "matrix_source_id": term.get("matrix_source_id"),
        "W0_source_id": term.get("W0_source_id"),
        "sector_identity": term.get("matrix_sector_identity", term.get("W0_sector_identity")),
        "term_index": term_index,
        "ordinary_Z_order": int(order),
        "field": field,
    }


def _bool_value(term: Mapping[str, Any], key: str, default: bool = False) -> bool:
    value = term.get(key, default)
    if not isinstance(value, bool):
        raise TypeError(f"{key} must be boolean when supplied")
    return value


def _matrix_terms(cell: Mapping[str, Any], matrix_name: str, order: str, key: str) -> list[Mapping[str, Any]]:
    sector = _lookup(cell, matrix_name, "matrix cell")
    if not isinstance(sector, Mapping):
        raise TypeError(f"matrix cell {matrix_name} sector must be a mapping")
    ordered = _lookup(sector, order, matrix_name)
    if not isinstance(ordered, Mapping):
        raise TypeError(f"matrix cell {matrix_name} order {order} must be a mapping")
    return _items(_lookup(ordered, key, f"{matrix_name}[{order}]"), f"{matrix_name}[{order}][{key}]")


def _optional_matrix_terms(cell: Mapping[str, Any], matrix_name: str, order: str, key: str) -> list[Mapping[str, Any]]:
    try:
        return _matrix_terms(cell, matrix_name, order, key)
    except (KeyError, TypeError):
        return []


def _state_terms(cell: Mapping[str, Any], order: str, field: str) -> list[Mapping[str, Any]]:
    family = cell.get("W1_by_order")
    if family is None:
        family = cell.get("W1_by_ordinary_Z")
    if family is None:
        family = cell.get("W1_Z" + order)
    if not isinstance(family, Mapping):
        raise KeyError(f"target cell is missing W1 order {order}")
    if "Z" + order in family and isinstance(family["Z" + order], Mapping):
        ordered = family["Z" + order]
    else:
        ordered = _lookup(family, order, "W1_by_order")
    if not isinstance(ordered, Mapping):
        raise TypeError(f"W1 order {order} must be a mapping")
    return _items(_lookup(ordered, field, f"W1 Z{order}"), f"W1 Z{order} {field}")


def _cell_t(c: Any, cell: Mapping[str, Any], index: int, label: str) -> Any:
    if "t" not in cell or "rho" not in cell:
        raise KeyError(f"{label} requires t and rho")
    t = _interval(c, cell["t"], label + ".t")
    expected = _canonical_cell(c, index)
    if not _same_interval(t, expected):
        raise ValueError(f"{label} has noncanonical t cell {index}")
    rho = _interval(c, cell["rho"], label + ".rho")
    if not _same_interval(rho, 4 * t * t):
        raise ValueError(f"{label} has inconsistent rho=4*t^2")
    return t


def _validate_profile_cell(c: Any, cell: Mapping[str, Any], index: int) -> Any:
    t = _cell_t(c, cell, index, "firstword target cell")
    for order in ORDERS:
        for field in STATE_FIELDS:
            terms = _state_terms(cell, order, field)
            if not terms:
                raise ValueError(f"target cell {index} W1 Z{order} {field} has no source terms")
            for term_index, term in enumerate(terms):
                _formal(term, f"target cell {index} W1 Z{order} {field}[{term_index}]")
                _coefficient(c, term, f"target cell {index} W1 Z{order} {field}[{term_index}]")
    return t


def _validate_matrix_cell(c: Any, cell: Mapping[str, Any], index: int) -> Any:
    t = _cell_t(c, cell, index, "raw matrix cell")
    for matrix_name in MATRIX_NAMES:
        keys = C0_MATRIX_KEYS if matrix_name == "C0" else B1_MATRIX_KEYS
        for order in ORDERS:
            for key in keys:
                terms = _matrix_terms(cell, matrix_name, order, key)
                if not terms:
                    raise ValueError(f"matrix cell {index} {matrix_name}[{order}][{key}] has no source terms")
                for term_index, term in enumerate(terms):
                    _formal(term, f"matrix cell {index} {matrix_name}[{order}][{key}][{term_index}]")
                    _coefficient(c, term, f"matrix cell {index} {matrix_name}[{order}][{key}][{term_index}]")
    return t


def _resolve_inputs(profile: Any, matrix_cells: Any) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], Mapping[str, Any] | None]:
    if _is_sequence(profile):
        target_cells = list(profile)
        envelope: Mapping[str, Any] | None = None
    elif isinstance(profile, Mapping):
        envelope = profile
        target_cells = profile.get("target_cells", profile.get("outer_cells"))
        if target_cells is None and "W1_target_cells" in profile:
            target_cells = profile["W1_target_cells"]
        if target_cells is None:
            raise KeyError("firstword profile requires target_cells")
        target_cells = list(target_cells)
    else:
        raise TypeError("firstword profile must be a mapping or target-cell sequence")
    if matrix_cells is None and envelope is not None:
        matrix_cells = envelope.get("matrix_cells", envelope.get("raw_matrix_cells"))
        if matrix_cells is None:
            extracted = []
            for cell in target_cells:
                if not isinstance(cell, Mapping):
                    raise TypeError("target cells must be mappings")
                nested = cell.get("matrix_cell", cell.get("raw_matrix"))
                if nested is None and "C0" in cell and "B1" in cell:
                    nested = cell
                if nested is None:
                    extracted = []
                    break
                extracted.append(nested)
            if extracted:
                matrix_cells = extracted
    if matrix_cells is None:
        raise KeyError("raw C0/B1 matrix_cells are required")
    matrix_cells = list(matrix_cells)
    if len(target_cells) != CELL_COUNT or len(matrix_cells) != CELL_COUNT:
        raise ValueError("exactly eight firstword target cells and eight matrix cells are required")
    if any(not isinstance(cell, Mapping) for cell in target_cells + matrix_cells):
        raise TypeError("target and matrix cells must be mappings")
    return target_cells, matrix_cells, envelope


def _extract_input_schema(envelope: Mapping[str, Any] | None) -> Any:
    if envelope is None:
        return INPUT_SCHEMA
    return envelope.get("schema", envelope.get("input_schema", INPUT_SCHEMA))


def _combine(
    c: Any,
    matrix_term: Mapping[str, Any],
    state_term: Mapping[str, Any],
    mass: Any,
    *,
    source_cell: int,
    source_t: Any,
    order: str,
    matrix_name: str,
    matrix_order: str,
    state_order: str,
    matrix_key: str,
    matrix_term_index: int,
    state_field: str,
    state_term_index: int,
    origin: str,
    target_x: Any,
    target_t: Any,
) -> dict[str, Any]:
    matrix_f0, matrix_lam = _formal(matrix_term, "matrix term")
    state_f0, state_lam = _formal(state_term, "firstword state term")
    matrix_coefficient = _coefficient(c, matrix_term, "matrix term")
    state_coefficient = _coefficient(c, state_term, "firstword state term")
    coefficient = matrix_coefficient * state_coefficient * mass
    matrix_classes = _term_classes(matrix_term) or ["raw_C0_B1_coefficient_interval"]
    state_classes = _term_classes(state_term) or ["saved_firstword_profile_coefficient_interval"]
    classes = matrix_classes + [item for item in state_classes if item not in matrix_classes]
    matrix_ref = _node_ref(matrix_term, role="matrix_source_node", cell_index=source_cell, term_index=matrix_term_index, order=matrix_order)
    state_ref = _node_ref(state_term, role="firstword_state_node", cell_index=source_cell, term_index=state_term_index, order=state_order, field=state_field)
    matrix_source_cell = matrix_term.get("source_cell", source_cell)
    state_source_cell = state_term.get("source_cell", source_cell)
    result = {
        "F0_power": matrix_f0 + state_f0,
        "Lambda_power": matrix_lam + state_lam,
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "directed": True,
        "function_value": False,
        "signed": True,
        "signed_field": _bool_value(matrix_term, "signed_field", True) and _bool_value(state_term, "signed_field", True),
        "ordinary_Z_order": int(order),
        "origin": origin,
        "input_origins": {
            "matrix": _clone(matrix_term.get("origin")),
            "state": _clone(state_term.get("origin")),
        },
        "matrix_name": matrix_name,
        "matrix_order": matrix_order,
        "matrix_key": matrix_key,
        "row": 5,
        "column": int(matrix_key.split(",")[1]),
        "state_field": state_field,
        "state_order": state_order,
        "matrix_term_index": matrix_term_index,
        "state_term_index": state_term_index,
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "source_radius_cell": source_cell,
        "source_t": source_t,
        "target_cell": "x=2 seam",
        "target_cell_index": -1,
        "target_radius_cell": -1,
        "target_x": target_x,
        "target_t": target_t,
        "kernel_mass": mass,
        "kernel_mass_field": None,
        "kernel_mass_positive": True,
        "kernel_mass_exact": True,
        "matrix_source_cell": matrix_source_cell,
        "matrix_source_cell_index": matrix_term.get("source_cell_index", source_cell),
        "state_source_cell": state_source_cell,
        "state_source_cell_index": state_term.get("source_cell_index", source_cell),
        "matrix_node_ref": matrix_ref,
        "source_node_ref": matrix_ref,
        "state_node_ref": state_ref,
        "source_ids": matrix_ref["source_ids"] + [item for item in state_ref["source_ids"] if item not in matrix_ref["source_ids"]],
        "state_ids": state_ref["source_ids"],
        "source_paths": matrix_ref["source_paths"] + [item for item in state_ref["source_paths"] if item not in matrix_ref["source_paths"]],
        "state_paths": state_ref["source_paths"],
        "error_classes": classes,
        "source_error_classes": matrix_classes,
        "state_error_classes": state_classes,
        "uncertainty": _bool_value(matrix_term, "uncertainty", False) or _bool_value(state_term, "uncertainty", False),
        "error_only": _bool_value(matrix_term, "error_only", False) or _bool_value(state_term, "error_only", False),
        "symmetric_error": _bool_value(matrix_term, "symmetric_error", False) or _bool_value(state_term, "symmetric_error", False),
        "scalar_multiplier": 1,
        "scalar_applied_once": True,
        "leibniz_multiplier": 1,
        "formal_atoms_materialized": False,
        "cancellation_performed": False,
        "no_cancellation": True,
    }
    return result


def _unknown_b1_z2(
    matrix_term: Mapping[str, Any],
    *,
    source_cell: int,
    matrix_term_index: int,
    matrix_key: str,
    state_field: str,
    source_t: Any,
    target_x: Any,
    target_t: Any,
) -> dict[str, Any]:
    matrix_ref = _node_ref(matrix_term, role="matrix_source_node", cell_index=source_cell, term_index=matrix_term_index, order="0")
    classes = _term_classes(matrix_term) or ["raw_B1_coefficient_interval"]
    if "unknown:B1[0]*W[2]" not in classes:
        classes.append("unknown:B1[0]*W[2]")
    unknown_ref = {
        "role": "missing_firstword_Z2_state_node",
        "node_id": f"unknown:firstword:cell:{source_cell}:Z2:{state_field}",
        "ids": [f"unknown:firstword:cell:{source_cell}:Z2:{state_field}"],
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "ordinary_Z_order": 2,
        "field": state_field,
        "missing": True,
        "error_classes": ["unknown:firstword_Z2_not_saved"],
    }
    return {
        "unknown_kind": "missing_B1_W1_Z2",
        "reason": "B1[0]*W[2] is outside the saved firstword Z0/Z1 profile",
        "origin": "B1[0]*W[2] (unknown)",
        "ordinary_Z_order": 1,
        "matrix_name": "B1",
        "matrix_order": "0",
        "matrix_key": matrix_key,
        "row": 5,
        "column": int(matrix_key.split(",")[1]),
        "state_field": state_field,
        "state_order": "2",
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "source_t": source_t,
        "target_cell": "x=2 seam",
        "target_cell_index": -1,
        "target_x": target_x,
        "target_t": target_t,
        "matrix_term_index": matrix_term_index,
        "matrix_node_ref": matrix_ref,
        "source_node_ref": matrix_ref,
        "state_node_ref": unknown_ref,
        "matrix_source_cell": matrix_term.get("source_cell", source_cell),
        "matrix_source_cell_index": matrix_term.get("source_cell_index", source_cell),
        "state_source_cell": source_cell,
        "state_source_cell_index": source_cell,
        "source_ids": matrix_ref["source_ids"],
        "state_ids": unknown_ref["ids"],
        "source_paths": matrix_ref["source_paths"],
        "state_paths": [],
        "error_classes": classes,
        "source_error_classes": _term_classes(matrix_term) or ["raw_B1_coefficient_interval"],
        "state_error_classes": ["unknown:firstword_Z2_not_saved"],
        "uncertainty": True,
        "error_only": True,
        "symmetric_error": False,
        "included": False,
        "materialized": False,
        "coefficient": None,
        "formal_powers": None,
        "kernel_mass": None,
        "no_cancellation": True,
    }


def _unknown_column(
    *,
    source_cell: int,
    source_t: Any,
    target_x: Any,
    target_t: Any,
    matrix_cell: Mapping[str, Any],
    target_cell: Mapping[str, Any],
    column: int,
    field: str,
) -> dict[str, Any]:
    matrix_refs: list[dict[str, Any]] = []
    for matrix_name in MATRIX_NAMES:
        for order in ORDERS:
            for term_index, term in enumerate(_optional_matrix_terms(matrix_cell, matrix_name, order, f"5,{column}")):
                matrix_refs.append(_node_ref(term, role="excluded_matrix_source_node", cell_index=source_cell, term_index=term_index, order=order))
    state_refs: list[dict[str, Any]] = []
    for order in ORDERS:
        for term_index, term in enumerate(_state_terms(target_cell, order, field)):
            state_refs.append(_node_ref(term, role="excluded_firstword_state_node", cell_index=source_cell, term_index=term_index, order=order, field=field))
    classes = ["unknown:firstword_column_not_admitted", "raw_matrix_or_firstword_profile_interval"]
    for ref in matrix_refs + state_refs:
        for class_name in ref.get("error_classes", []):
            if class_name not in classes:
                classes.append(class_name)
    return {
        "unknown_kind": "excluded_firstword_column",
        "reason": f"row 5 column {column} ({field}) is outside the admitted F1/Uz1 source subset",
        "origin": f"column {column} ({field}) unknown",
        "ordinary_Z_orders": [0, 1],
        "matrix_names": ["C0", "B1"],
        "matrix_keys": [f"5,{column}"],
        "row": 5,
        "column": column,
        "state_field": field,
        "state_orders": [0, 1],
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "source_t": source_t,
        "target_cell": "x=2 seam",
        "target_cell_index": -1,
        "target_x": target_x,
        "target_t": target_t,
        "matrix_node_refs": matrix_refs,
        "state_node_refs": state_refs,
        "source_node_refs": matrix_refs,
        "state_node_ref": state_refs[0] if state_refs else {"missing": True, "field": field, "source_cell": source_cell},
        "source_ids": [item for ref in matrix_refs for item in ref.get("source_ids", [])],
        "state_ids": [item for ref in state_refs for item in ref.get("source_ids", [])],
        "error_classes": classes,
        "source_error_classes": [],
        "state_error_classes": [],
        "uncertainty": True,
        "error_only": True,
        "included": False,
        "materialized": False,
        "coefficient": None,
        "formal_powers": None,
        "no_cancellation": True,
    }


def _flatten(values: Mapping[str, Mapping[str, list[dict[str, Any]]]], field: str) -> list[dict[str, Any]]:
    return list(values["Z0"][field]) + list(values["Z1"][field])


def build(c: Any, profile: Any, matrix_cells: Any = None, *, target_x: Any = None) -> dict[str, Any]:
    """Build the finite row-5 raw-K0 second-feedback source subset.

    ``profile`` is either the result mapping from the saved firstword profile
    (with ``target_cells``) or the eight target cells themselves.  ``matrix_cells``
    is the corresponding eight-cell raw C0/B1 list; it may be omitted only when
    the profile envelope or each target cell carries it.  The optional
    ``target_x`` is a guard and must equal 2 when supplied.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide mpmath interval mpf and make_mpf")
    if target_x is not None:
        x_guard = _interval(c, target_x, "target_x")
        if x_guard.a != 2 or x_guard.b != 2:
            raise ValueError("the partial source leaf is fixed at x=2")
    target_cells, raw_matrices, envelope = _resolve_inputs(profile, matrix_cells)
    target_t = c.mpf(1)
    target_x_iv = c.mpf(2)

    source_cells: list[dict[str, Any]] = []
    output = {"Z0": {field: [] for field in OUTPUT_FIELDS}, "Z1": {field: [] for field in OUTPUT_FIELDS}}
    unknowns: list[dict[str, Any]] = []
    kernels: dict[str, dict[str, Any]] = {}
    for source_index, (target_cell, matrix_cell) in enumerate(zip(target_cells, raw_matrices)):
        if "index" in target_cell and target_cell["index"] != source_index:
            raise ValueError(f"target cell index {source_index} is not canonical")
        if "index" in matrix_cell and matrix_cell["index"] != source_index:
            raise ValueError(f"matrix cell index {source_index} is not canonical")
        source_t = _validate_profile_cell(c, target_cell, source_index)
        matrix_t = _validate_matrix_cell(c, matrix_cell, source_index)
        if not _same_interval(source_t, matrix_t):
            raise ValueError(f"target and matrix cell {source_index} have different t intervals")
        mass_pair = target_seam_kernel_masses(c, source_index)
        kernels[str(source_index)] = {
            "source_cell": source_index,
            "source_cell_index": source_index,
            "t": source_t,
            "a": mass_pair["a"],
            "b": mass_pair["b"],
            "F1": mass_pair["F1"],
            "d_x_F1": mass_pair["d_x_F1"],
            "dxF1": mass_pair["d_x_F1"],
            "formula": mass_pair["formula"],
            "positive": True,
            "exact": True,
            "independent_polynomial": True,
        }
        cell_term_ids = {"Z0": {field: [] for field in OUTPUT_FIELDS}, "Z1": {field: [] for field in OUTPUT_FIELDS}}
        cell_unknown_ids: list[str] = []
        for order, rules in PRODUCT_RULES.items():
            for matrix_name, matrix_order, state_order, origin in rules:
                matrix_keys = C0_MATRIX_KEYS if matrix_name == "C0" else B1_MATRIX_KEYS
                for matrix_key in matrix_keys:
                    column = int(matrix_key.split(",")[1])
                    state_field = STATE_COLUMN[str(column)]
                    matrix_terms = _matrix_terms(matrix_cell, matrix_name, matrix_order, matrix_key)
                    state_terms = _state_terms(target_cell, state_order, state_field)
                    for matrix_term_index, matrix_term in enumerate(matrix_terms):
                        for state_term_index, state_term in enumerate(state_terms):
                            for field in OUTPUT_FIELDS:
                                term = _combine(
                                    c,
                                    matrix_term,
                                    state_term,
                                    mass_pair[field],
                                    source_cell=source_index,
                                    source_t=source_t,
                                    order=order[-1],
                                    matrix_name=matrix_name,
                                    matrix_order=matrix_order,
                                    state_order=state_order,
                                    matrix_key=matrix_key,
                                    matrix_term_index=matrix_term_index,
                                    state_field=state_field,
                                    state_term_index=state_term_index,
                                    origin=origin,
                                    target_x=target_x_iv,
                                    target_t=target_t,
                                )
                                term["kernel_mass_field"] = field
                                term["term_id"] = f"secondword:cell:{source_index}:Z{order[-1]}:{field}:{len(output[order][field])}"
                                cell_term_ids[order][field].append(term["term_id"])
                                output[order][field].append(term)
            for matrix_key in B1_MATRIX_KEYS:
                column = int(matrix_key.split(",")[1])
                state_field = STATE_COLUMN[str(column)]
                if order == "Z1":
                    for matrix_term_index, matrix_term in enumerate(_matrix_terms(matrix_cell, "B1", "0", matrix_key)):
                        unknown = _unknown_b1_z2(
                            matrix_term,
                            source_cell=source_index,
                            matrix_term_index=matrix_term_index,
                            matrix_key=matrix_key,
                            state_field=state_field,
                            source_t=source_t,
                            target_x=target_x_iv,
                            target_t=target_t,
                        )
                        unknown["unknown_id"] = f"secondword:cell:{source_index}:unknown:{len(unknowns)}"
                        cell_unknown_ids.append(unknown["unknown_id"])
                        unknowns.append(unknown)
        for column in UNKNOWN_COLUMNS:
            field = STATE_COLUMN[str(column)]
            unknown = _unknown_column(
                source_cell=source_index,
                source_t=source_t,
                target_x=target_x_iv,
                target_t=target_t,
                matrix_cell=matrix_cell,
                target_cell=target_cell,
                column=column,
                field=field,
            )
            unknown["unknown_id"] = f"secondword:cell:{source_index}:unknown:{len(unknowns)}"
            cell_unknown_ids.append(unknown["unknown_id"])
            unknowns.append(unknown)
        source_cells.append(
            {
                "index": source_index,
                "source_cell": source_index,
                "source_cell_index": source_index,
                "t": source_t,
                "rho": _interval(c, target_cell["rho"], f"target cell {source_index}.rho"),
                "kernel_id": str(source_index),
                "kernel_ref": f"kernels[{source_index}]",
                "known_term_ids_by_order": cell_term_ids,
                "known_term_counts_by_order": {
                    order: {field: len(values) for field, values in fields.items()}
                    for order, fields in cell_term_ids.items()
                },
                "unknown_ids": cell_unknown_ids,
                "unknown_count": len(cell_unknown_ids),
                "firstword_state_source_cell": source_index,
                "raw_matrix_source_cell": source_index,
            }
        )

    for field in OUTPUT_FIELDS:
        for order in ORDERS:
            for term in output[order][field]:
                term["target_field"] = field

    by_order = output
    term_counts = {
        order: {field: len(by_order[order][field]) for field in OUTPUT_FIELDS}
        for order in ORDERS
    }
    term_ids_by_field = {
        field: {
            order: [term["term_id"] for term in by_order[order][field]]
            for order in ORDERS
        }
        for field in OUTPUT_FIELDS
    }
    unknown_ids_by_kind: dict[str, list[str]] = {}
    for item in unknowns:
        unknown_ids_by_kind.setdefault(item["unknown_kind"], []).append(item["unknown_id"])
    return {
        "schema": SCHEMA,
        "input_schema": _extract_input_schema(envelope),
        "source_cells": source_cells,
        "target_seam": {
            "x": target_x_iv,
            "t": target_t,
            "rho": c.mpf(4),
            "fixed": True,
            "target_cell_index": -1,
            "kernel_mass_scope": "source cell [i/8,(i+1)/8] at x=2",
        },
        "kernels": kernels,
        "terms_by_order": by_order,
        "term_counts_by_order": term_counts,
        "term_ids_by_field": term_ids_by_field,
        "unknowns": unknowns,
        "unknown_ids_by_kind": unknown_ids_by_kind,
        "unknown_residual": {
            "symbol": "R",
            "name": "R_secondword_partial",
            "included": False,
            "included_in_products": False,
            "added_to_bounds": False,
            "materialized": False,
            "categories": ["B1[0]*W[2] (columns 1,2,3)"],
            "unknown_count": len(unknowns),
        },
        "ordinary_Z_orders": [0, 1],
        "known_product_rules": {
            "Z0": [rule[3] for rule in PRODUCT_RULE_Z0],
            "Z1": [rule[3] for rule in PRODUCT_RULE_Z1],
        },
        "missing_product_rule": "B1[0]*W[2]",
        "admitted_columns": {"1": "F1", "2": "Uz1", "3": "K1", "5": "d_x_F1"},
        "raw_matrix_columns": {"C0": [1, 2, 3, 5], "B1": [1, 2, 3]},
        "row5_structurally_absent_columns": {"C0": [4, 6], "B1": [4, 5, 6]},
        "known_matrix_keys": {"C0": list(C0_MATRIX_KEYS), "B1": list(B1_MATRIX_KEYS)},
        "excluded_columns_by_name": {},
        "source_power_basis": ["F0_power", "Lambda_power"],
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "source_node_refs_retained": True,
        "state_node_refs_retained": True,
        "error_classes_retained": True,
        "kernel_mass_exact": True,
        "kernel_mass_positive": True,
        "kernel_mass_independent_polynomial": True,
        "no_cancellation": True,
        "cancellation_performed": False,
        "formal_atoms_materialized": False,
        "F0_materialized": False,
        "Lambda_materialized": False,
        "ancestor_called": False,
        "ancestor_reconstructed": False,
        "firstword_reconstructed": False,
        "complete_secondword": False,
        "full_secondword_claimed": False,
        "Z2_materialized": False,
        "m15_computed": False,
        "ledger_appended": False,
        "owner_proved": False,
        "physical_field_claimed": False,
    }


build_secondword_partial = build

