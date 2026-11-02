"""Structural source binding for the finite ``B1[0]*W1_Z2`` seam.

This module is deliberately a consumer, not another Green or matrix leaf.  It
joins the explicit ``missing_B1_W1_Z2`` records from
``core_secondword_partial_v2`` to the already emitted terms from
``core_firstword_z2_profile_v2.build_missing_seam``.  No receipt is opened,
no ancestor module is imported, and no coefficient or kernel value is
recomputed here.

The bridge is finite and conditional.  It does not turn an unresolved old
unknown into a physical owner, a complete second word, an all-kernel result,
or an ``m15`` account.  Every input term is copied as raw metadata under the
returned binding; in particular directed coefficient endpoints, IDs, and all
error-class tags remain visible.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "core_z2_seam_bridge_v3"
OLD_SCHEMA = "core_secondword_partial_v2"
PROFILE_SCHEMA = "core_firstword_z2_profile_v2"
SEAM_SCHEMA = PROFILE_SCHEMA + ":missing_seam_v1"
CELL_COUNT = 8
SEAM_FIELDS = ("F1", "d_x_F1")
SEAM_FIELD_ALIASES = {"F1": "F1", "Fx": "d_x_F1", "d_x_F1": "d_x_F1", "dxF1": "d_x_F1"}
STATE_FIELD_BY_COLUMN = {1: "F1", 2: "Uz1", 3: "K1"}

__all__ = [
    "SCHEMA",
    "OLD_SCHEMA",
    "PROFILE_SCHEMA",
    "SEAM_SCHEMA",
    "CELL_COUNT",
    "SEAM_FIELDS",
    "build",
    "build_z2_seam_bridge",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    """Copy mappings/lists without touching or interpreting numeric payloads."""

    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _lookup(mapping: Mapping[Any, Any], key: Any, label: str) -> Any:
    if key in mapping:
        return mapping[key]
    alternate = int(key) if isinstance(key, str) and key.isdigit() else str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(f"{label} is missing key {key!r}")


def _mapping(value: Any, label: str) -> Mapping[Any, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def _items(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        return [value]
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a term list or mapping")
    result = list(value)
    if any(not isinstance(item, Mapping) for item in result):
        raise TypeError(f"{label} contains a non-mapping term")
    return result


def _as_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    return int(value)


def _as_order(value: Any, expected: int, label: str) -> None:
    if value not in (expected, str(expected), f"Z{expected}"):
        raise ValueError(f"{label} must be {expected!r}/{str(expected)!r}/'Z{expected}'")


def _exact_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int] | None:
    """Return an exact MPF endpoint tuple, if ``value`` visibly has one."""

    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        value = raw
    if not _is_sequence(value) or len(value) != 4:
        return None
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} exact MPF endpoint must contain integers")
    return tuple(int(item) for item in value)


def _endpoint_signature(value: Any, label: str) -> tuple[int, int, int, int] | None:
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return _exact_mpf_tuple(raw, label)
    return _exact_mpf_tuple(value, label)


def _interval_signature(value: Any, label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]] | None:
    """Extract directed endpoints without constructing an interval.

    Both live ``mpmath.iv`` values and the JSON records emitted by the frozen
    runners are accepted.  A missing signature returns ``None`` so callers
    comparing ordinary path/metadata values can use the same structural
    equality helper; source radii and masses explicitly require signatures.
    """

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        if not _is_sequence(raw) or len(raw) != 2:
            raise TypeError(f"{label} has malformed directed interval")
        lower = _endpoint_signature(raw[0], label + ".lower")
        upper = _endpoint_signature(raw[1], label + ".upper")
        if lower is None or upper is None:
            raise TypeError(f"{label} does not expose exact MPF endpoints")
        return lower, upper

    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            if not has_lower or not has_upper:
                raise TypeError(f"{label} requires both exact directed endpoints")
            lower = _endpoint_signature(value["lower_exact_mpf_tuple"], label + ".lower")
            upper = _endpoint_signature(value["upper_exact_mpf_tuple"], label + ".upper")
            if lower is None or upper is None:
                raise TypeError(f"{label} has malformed exact endpoint records")
            return lower, upper
        # Some callers retain a live interval-like ``a``/``b`` wrapper.
        if "a" in value or "b" in value:
            if "a" not in value or "b" not in value:
                raise TypeError(f"{label} requires both interval endpoints")
            lower = _endpoint_signature(value["a"], label + ".lower")
            upper = _endpoint_signature(value["b"], label + ".upper")
            if lower is not None and upper is not None:
                return lower, upper

    if hasattr(value, "a") and hasattr(value, "b"):
        lower = _endpoint_signature(value.a, label + ".lower")
        upper = _endpoint_signature(value.b, label + ".upper")
        if lower is None or upper is None:
            raise TypeError(f"{label} does not expose exact MPF endpoints")
        return lower, upper

    if _is_sequence(value) and len(value) == 2:
        lower = _endpoint_signature(value[0], label + ".lower")
        upper = _endpoint_signature(value[1], label + ".upper")
        if lower is not None and upper is not None:
            return lower, upper
    return None


def _structural_equal(left: Any, right: Any, label: str = "value") -> bool:
    """Compare metadata and intervals without selecting a midpoint."""

    left_interval = _interval_signature(left, label + ".left")
    right_interval = _interval_signature(right, label + ".right")
    if left_interval is not None or right_interval is not None:
        return left_interval is not None and left_interval == right_interval
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping):
            return False
        if list(left.keys()) != list(right.keys()):
            # For metadata dictionaries order is not semantic, but retaining
            # exact key identity avoids comparing an interval wrapper as a
            # generic dictionary.
            if set(left.keys()) != set(right.keys()):
                return False
        return all(key in right and _structural_equal(value, right[key], f"{label}.{key}") for key, value in left.items())
    if _is_sequence(left) or _is_sequence(right):
        if not _is_sequence(left) or not _is_sequence(right) or len(left) != len(right):
            return False
        return all(_structural_equal(a, b, f"{label}[{index}]") for index, (a, b) in enumerate(zip(left, right)))
    if isinstance(left, float) or isinstance(right, float):
        # A float is never an exact MPF endpoint.  Returning False makes a
        # mixed float/interval or float/scale binding fail closed.
        return False
    return type(left) is type(right) and left == right


def _same_exact_interval(left: Any, right: Any, label: str) -> None:
    left_signature = _interval_signature(left, label + ".left")
    right_signature = _interval_signature(right, label + ".right")
    if left_signature is None or right_signature is None:
        raise TypeError(f"{label} requires exact directed MPF endpoints")
    if left_signature != right_signature:
        raise ValueError(f"{label} directed endpoints differ")


def _same_exact_value(left: Any, right: Any, label: str) -> None:
    if not _structural_equal(left, right, label):
        raise ValueError(f"{label} differs")


def _unwrap(value: Any, label: str, preferred: str | None = None) -> Mapping[str, Any]:
    """Accept a direct result or the runner's ``output``/``result`` envelope."""

    current: Any = value
    for _ in range(4):
        if not isinstance(current, Mapping):
            raise TypeError(f"{label} must be a result mapping")
        if preferred is not None and preferred in current and isinstance(current[preferred], Mapping):
            current = current[preferred]
            preferred = None
            continue
        if "output" in current and isinstance(current["output"], Mapping):
            current = current["output"]
            continue
        if "result" in current and isinstance(current["result"], Mapping) and "schema" not in current:
            current = current["result"]
            continue
        break
    if not isinstance(current, Mapping):
        raise TypeError(f"{label} must resolve to a result mapping")
    return current


def _validate_schema(result: Mapping[str, Any], expected: str, label: str) -> None:
    schema = result.get("schema")
    if schema != expected:
        raise ValueError(f"{label} schema must be {expected!r}, got {schema!r}")


def _target_cells(profile: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    values = profile.get("target_cells", profile.get("W1_Z2_target_cells"))
    if not _is_sequence(values):
        raise KeyError("new profile result requires target_cells")
    cells = list(values)
    if len(cells) != CELL_COUNT or any(not isinstance(cell, Mapping) for cell in cells):
        raise ValueError("new profile result requires exactly eight target_cells")
    for index, cell in enumerate(cells):
        if "index" in cell and cell["index"] != index:
            raise ValueError(f"profile target cell {index} is not canonical")
    return cells


def _profile_state_terms(cell: Mapping[str, Any], field: str, label: str) -> list[Mapping[str, Any]]:
    family: Any = cell.get("W1_Z2")
    if family is None:
        family = cell.get("W1_by_order")
    if isinstance(family, Mapping) and "Z2" in family:
        family = family["Z2"]
    elif isinstance(family, Mapping) and 2 in family:
        family = family[2]
    if family is None:
        family = cell.get("terms_Z2", cell.get("terms"))
        if isinstance(family, Mapping) and "Z2" in family:
            family = family["Z2"]
    if not isinstance(family, Mapping):
        raise KeyError(f"{label} is missing W1_Z2")
    if field not in family:
        raise KeyError(f"{label} is missing W1_Z2[{field!r}]")
    terms = _items(family[field], f"{label}.W1_Z2[{field}]")
    if not terms:
        raise ValueError(f"{label}.W1_Z2[{field}] is empty")
    for index, term in enumerate(terms):
        term_id = term.get("term_id")
        if not isinstance(term_id, str) or not term_id:
            raise KeyError(f"{label}.W1_Z2[{field}][{index}] requires term_id")
    return terms


def _seam_terms(seam: Mapping[str, Any], field: str) -> list[Mapping[str, Any]]:
    by_order = seam.get("terms_by_order", seam.get("terms"))
    if not isinstance(by_order, Mapping):
        raise KeyError("new seam result requires terms_by_order")
    ordered: Any
    if "Z1" in by_order:
        ordered = by_order["Z1"]
    elif 1 in by_order:
        ordered = by_order[1]
    elif "1" in by_order:
        ordered = by_order["1"]
    else:
        raise KeyError("new seam result is missing ordinary Z1 terms")
    if not isinstance(ordered, Mapping) or field not in ordered:
        raise KeyError(f"new seam result is missing Z1[{field}]")
    return _items(ordered[field], f"new seam Z1[{field}]")


def _parse_matrix_key(value: Any, label: str) -> tuple[str, int]:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a row,column string")
    pieces = value.split(",")
    if len(pieces) != 2 or pieces[0] != "5" or not pieces[1].isdigit():
        raise ValueError(f"{label} must be a row-5 matrix key")
    column = int(pieces[1])
    if column not in STATE_FIELD_BY_COLUMN:
        raise ValueError(f"{label} must use B1 columns 1, 2, or 3")
    return value, column


def _source_paths(unknown: Mapping[str, Any], label: str) -> list[Any]:
    node = unknown.get("matrix_node_ref")
    if not isinstance(node, Mapping):
        raise KeyError(f"{label} requires matrix_node_ref")
    paths = node.get("source_paths")
    if not _is_sequence(paths) or not paths:
        raise KeyError(f"{label}.matrix_node_ref.source_paths must be a nonempty list")
    return list(paths)


def _find_path(paths: Sequence[Any], candidate: Any, label: str) -> Any:
    for path in paths:
        if _structural_equal(path, candidate, label + ".path"):
            return path
    raise ValueError(f"{label} does not contain the LOW matrix receipt input path")


def _kernel(result: Mapping[str, Any], source_cell: int, label: str) -> Mapping[str, Any]:
    kernels = result.get("kernels")
    if not isinstance(kernels, Mapping):
        raise KeyError(f"{label} requires kernels")
    kernel = _lookup(kernels, str(source_cell), f"{label}.kernels")
    return _mapping(kernel, f"{label}.kernels[{source_cell}]")


def _target_seam(result: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    seam = result.get("target_seam")
    if not isinstance(seam, Mapping):
        raise KeyError(f"{label} requires target_seam")
    if "x" not in seam or "t" not in seam:
        raise KeyError(f"{label}.target_seam requires x and t")
    return seam


def _term_key(term: Mapping[str, Any], label: str) -> tuple[int, str, int, str, int, str]:
    source_cell = _as_int(term.get("source_cell"), label + ".source_cell")
    if not 0 <= source_cell < CELL_COUNT:
        raise ValueError(f"{label}.source_cell is outside [0, 7]")
    matrix_key, column = _parse_matrix_key(term.get("matrix_key"), label + ".matrix_key")
    matrix_term_index = _as_int(term.get("matrix_term_index"), label + ".matrix_term_index")
    if matrix_term_index < 0:
        raise ValueError(f"{label}.matrix_term_index must be nonnegative")
    state_field = term.get("state_field")
    if state_field != STATE_FIELD_BY_COLUMN[column]:
        raise ValueError(f"{label}.state_field does not match matrix column")
    state_term_index = _as_int(term.get("state_term_index"), label + ".state_term_index")
    if state_term_index < 0:
        raise ValueError(f"{label}.state_term_index must be nonnegative")
    target_field = term.get("target_field")
    target_field = SEAM_FIELD_ALIASES.get(target_field, target_field)
    if target_field not in SEAM_FIELDS:
        raise ValueError(f"{label}.target_field must be F1 or d_x_F1")
    return source_cell, matrix_key, matrix_term_index, state_field, state_term_index, target_field


def _validate_old_unknown(unknown: Mapping[str, Any], old: Mapping[str, Any], old_seam: Mapping[str, Any]) -> tuple[str, tuple[int, str, int, str], list[Any]]:
    label = f"old unknown {unknown.get('unknown_id', '<missing>')!r}"
    unknown_id = unknown.get("unknown_id")
    if not isinstance(unknown_id, str) or not unknown_id:
        raise KeyError(f"{label} requires a nonempty unknown_id")
    if unknown.get("unknown_kind") != "missing_B1_W1_Z2":
        raise ValueError(f"{label} is not a missing_B1_W1_Z2 record")
    source_cell = _as_int(unknown.get("source_cell"), label + ".source_cell")
    if not 0 <= source_cell < CELL_COUNT:
        raise ValueError(f"{label}.source_cell is outside [0, 7]")
    matrix_key, column = _parse_matrix_key(unknown.get("matrix_key"), label + ".matrix_key")
    matrix_term_index = _as_int(unknown.get("matrix_term_index"), label + ".matrix_term_index")
    if matrix_term_index < 0:
        raise ValueError(f"{label}.matrix_term_index must be nonnegative")
    if unknown.get("matrix_name") != "B1":
        raise ValueError(f"{label}.matrix_name must be B1")
    _as_order(unknown.get("matrix_order"), 0, label + ".matrix_order")
    if unknown.get("row") != 5:
        raise ValueError(f"{label}.row must be 5")
    if unknown.get("column") != column:
        raise ValueError(f"{label}.column does not match matrix_key")
    state_field = unknown.get("state_field")
    if state_field != STATE_FIELD_BY_COLUMN[column]:
        raise ValueError(f"{label}.state_field does not match matrix column")
    _as_order(unknown.get("state_order"), 2, label + ".state_order")
    _as_order(unknown.get("ordinary_Z_order"), 1, label + ".ordinary_Z_order")
    paths = _source_paths(unknown, label)
    if "source_t" not in unknown:
        raise KeyError(f"{label} requires source_t")
    if "target_x" not in unknown or "target_t" not in unknown:
        raise KeyError(f"{label} requires target_x and target_t")
    old_seam_x = old_seam["x"]
    old_seam_t = old_seam["t"]
    _same_exact_value(unknown["target_x"], old_seam_x, label + ".target_x")
    _same_exact_value(unknown["target_t"], old_seam_t, label + ".target_t")
    return unknown_id, (source_cell, matrix_key, matrix_term_index, state_field), paths


def _validate_product(product: Mapping[str, Any], seam: Mapping[str, Any], label: str) -> tuple[tuple[int, str, int, str, int, str], str, Any]:
    term_id = product.get("term_id")
    if not isinstance(term_id, str) or not term_id:
        raise KeyError(f"{label} requires a nonempty term_id")
    key = _term_key(product, label)
    if product.get("matrix_name") != "B1":
        raise ValueError(f"{label}.matrix_name must be B1")
    _as_order(product.get("matrix_order"), 0, label + ".matrix_order")
    if product.get("row") != 5:
        raise ValueError(f"{label}.row must be 5")
    _as_order(product.get("state_order"), 2, label + ".state_order")
    _as_order(product.get("ordinary_Z_order"), 1, label + ".ordinary_Z_order")
    source_cell, matrix_key, matrix_term_index, state_field, _state_index, target_field = key
    if product.get("source_cell_index", source_cell) != source_cell:
        raise ValueError(f"{label}.source_cell_index differs")
    if product.get("matrix_source_cell", source_cell) != source_cell or product.get("matrix_source_cell_index", source_cell) != source_cell:
        raise ValueError(f"{label} matrix source cell differs")
    if "matrix_receipt_input_path" not in product:
        raise KeyError(f"{label} requires matrix_receipt_input_path")
    if "state_profile_term_id" not in product or not isinstance(product["state_profile_term_id"], str) or not product["state_profile_term_id"]:
        raise KeyError(f"{label} requires state_profile_term_id")
    if "source_t" not in product or "target_x" not in product or "target_t" not in product:
        raise KeyError(f"{label} requires source_t and target scales")
    seam_target = _target_seam(seam, "new seam result")
    _same_exact_value(product["target_x"], seam_target["x"], label + ".target_x")
    _same_exact_value(product["target_t"], seam_target["t"], label + ".target_t")
    if product.get("kernel_mass_field") not in (None, target_field):
        # The producer labels the field with the canonical name.  ``Fx`` is
        # accepted only as an alias for the derivative output field.
        if SEAM_FIELD_ALIASES.get(product.get("kernel_mass_field")) != target_field:
            raise ValueError(f"{label}.kernel_mass_field differs from target_field")
    if "kernel_mass" not in product:
        raise KeyError(f"{label} requires kernel_mass")
    return key, term_id, product["matrix_receipt_input_path"]


def _expected_product_key(base: tuple[int, str, int, str], target_field: str, state_term_index: int) -> tuple[int, str, int, str, int, str]:
    return base + (state_term_index, target_field)


_COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "coefficient_enclosure",
    "directed_interval",
    "enclosure",
    "interval",
)


def _unwrap_low_output(value: Any) -> Mapping[str, Any]:
    """Resolve a saved LOW result without opening a file or decoding bytes."""

    current: Any = value
    for _ in range(5):
        if not isinstance(current, Mapping):
            raise TypeError("low_output must resolve to a mapping")
        if "matrix_cells" in current:
            return current
        if "output" in current and isinstance(current["output"], Mapping):
            current = current["output"]
            continue
        if "result" in current and isinstance(current["result"], Mapping):
            current = current["result"]
            continue
        break
    raise KeyError("low_output requires saved matrix_cells")


def _follow_path(root: Any, path: Sequence[Any]) -> Any:
    current = root
    for key in path:
        if isinstance(current, Mapping):
            if key in current:
                current = current[key]
            elif isinstance(key, int) and str(key) in current:
                current = current[str(key)]
            elif isinstance(key, str) and key.isdigit() and int(key) in current:
                current = current[int(key)]
            else:
                raise KeyError(key)
        elif _is_sequence(current) and isinstance(key, int):
            current = current[key]
        else:
            raise KeyError(key)
    return current


def _locate_low_term(low_output: Mapping[str, Any], path: Any, label: str) -> Mapping[str, Any]:
    if not _is_sequence(path) or not path:
        raise TypeError(f"{label}.matrix_receipt_input_path must be a nonempty path")
    # A receipt path starts with the LOW receipt filename.  Depending on the
    # caller, low_output may be the complete receipt, its result envelope, or
    # just result.output.  Trying suffixes handles all three without reading a
    # path or reconstructing a source packet.
    failures: list[str] = []
    for start in range(len(path)):
        try:
            candidate = _follow_path(low_output, list(path[start:]))
        except (KeyError, IndexError, TypeError) as exc:
            failures.append(type(exc).__name__)
            continue
        if isinstance(candidate, Mapping):
            return candidate
        failures.append("not-mapping")
    raise KeyError(f"{label} cannot locate LOW matrix term from receipt path ({','.join(failures)})")


def _coefficient_value(term: Mapping[str, Any], label: str) -> Any:
    for key in _COEFFICIENT_KEYS:
        if key in term:
            return term[key]
    raise KeyError(f"{label} requires a directed coefficient")


def _formal_power(term: Mapping[str, Any], name: str, label: str) -> int:
    value = term.get(name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label}.{name} must be an integer")
    return int(value)


def _validate_low_payload(
    *,
    low_output: Mapping[str, Any],
    product: Mapping[str, Any],
    profile_term: Mapping[str, Any],
    label: str,
) -> None:
    """Check formal powers and raw input intervals against saved sources."""

    low_term = _locate_low_term(low_output, product["matrix_receipt_input_path"], label)
    matrix_f0 = _formal_power(low_term, "F0_power", label + ".LOW matrix")
    matrix_lam = _formal_power(low_term, "Lambda_power", label + ".LOW matrix")
    state_f0 = _formal_power(profile_term, "F0_power", label + ".profile state")
    state_lam = _formal_power(profile_term, "Lambda_power", label + ".profile state")
    product_f0 = _formal_power(product, "F0_power", label + ".seam product")
    product_lam = _formal_power(product, "Lambda_power", label + ".seam product")
    if product_f0 != matrix_f0 + state_f0:
        raise ValueError(f"{label}.F0_power does not equal LOW matrix + profile state powers")
    if product_lam != matrix_lam + state_lam:
        raise ValueError(f"{label}.Lambda_power does not equal LOW matrix + profile state powers")
    if "input_matrix_coefficient" not in product:
        raise KeyError(f"{label} requires input_matrix_coefficient")
    if "input_state_coefficient" not in product:
        raise KeyError(f"{label} requires input_state_coefficient")
    _same_exact_interval(
        product["input_matrix_coefficient"],
        _coefficient_value(low_term, label + ".LOW matrix"),
        label + ".input_matrix_coefficient",
    )
    _same_exact_interval(
        product["input_state_coefficient"],
        _coefficient_value(profile_term, label + ".profile state"),
        label + ".input_state_coefficient",
    )


def build(old: Any, new_profile: Any, new_seam: Any, low_output: Any = None) -> dict[str, Any]:
    """Bind old missing records to new seam products by source identity.

    Parameters are complete result objects (or their runner ``output`` /
    ``result`` envelopes).  ``low_output`` is the authenticated saved LOW
    result (or its envelope) and is used only to locate matrix terms by each
    seam term's ``matrix_receipt_input_path``.  The function intentionally has
    no interval context argument: exact MPF comparisons use endpoint tuples
    already present in the supplied objects and never recompute a number.
    """

    old_result = _unwrap(old, "old result")
    profile_result = _unwrap(new_profile, "new profile result")
    seam_result = _unwrap(new_seam, "new seam result", preferred="missing_seam_result")
    if low_output is None:
        raise TypeError("build requires the saved low_output fourth argument")
    low_result = _unwrap_low_output(low_output)
    _validate_schema(old_result, OLD_SCHEMA, "old result")
    _validate_schema(profile_result, PROFILE_SCHEMA, "new profile result")
    _validate_schema(seam_result, SEAM_SCHEMA, "new seam result")

    old_target_seam = _target_seam(old_result, "old result")
    new_target_seam = _target_seam(seam_result, "new seam result")
    _same_exact_value(old_target_seam["x"], new_target_seam["x"], "target seam x")
    _same_exact_value(old_target_seam["t"], new_target_seam["t"], "target seam t")

    profile_cells = _target_cells(profile_result)
    old_unknowns_value = old_result.get("unknowns")
    if not _is_sequence(old_unknowns_value):
        raise KeyError("old result requires unknowns")
    old_unknowns = list(old_unknowns_value)
    if any(not isinstance(item, Mapping) for item in old_unknowns):
        raise TypeError("old result unknowns must be mappings")

    # Validate and index every profile state term before reading seam products.
    state_index: dict[tuple[int, str], list[Mapping[str, Any]]] = {}
    state_id_owner: dict[str, tuple[int, str, int]] = {}
    for source_cell, cell in enumerate(profile_cells):
        if "t" not in cell:
            raise KeyError(f"profile target cell {source_cell} requires t")
        for field in STATE_FIELD_BY_COLUMN.values():
            terms = _profile_state_terms(cell, field, f"profile target cell {source_cell}")
            state_index[(source_cell, field)] = terms
            for state_term_index, term in enumerate(terms):
                state_id = term["term_id"]
                if state_id in state_id_owner:
                    raise ValueError(f"duplicate profile state term_id {state_id!r}")
                state_id_owner[state_id] = (source_cell, field, state_term_index)

    old_seam_x = old_target_seam["x"]
    old_seam_t = old_target_seam["t"]
    old_records: dict[str, dict[str, Any]] = {}
    old_key_owner: dict[tuple[int, str, int, str], str] = {}
    for unknown in old_unknowns:
        unknown_id, base_key, paths = _validate_old_unknown(unknown, old_result, old_target_seam)
        if unknown_id in old_records:
            raise ValueError(f"duplicate old unknown_id {unknown_id!r}")
        if base_key in old_key_owner:
            raise ValueError(f"duplicate old unknown binding key {base_key!r}")
        source_cell, matrix_key, matrix_term_index, state_field = base_key
        profile_terms = state_index[(source_cell, state_field)]
        old_kernel = _kernel(old_result, source_cell, "old result")
        new_kernel = _kernel(seam_result, source_cell, "new seam result")
        profile_cell = profile_cells[source_cell]
        _same_exact_interval(unknown["source_t"], profile_cell["t"], f"{unknown_id}.source_t/profile t")
        _same_exact_interval(unknown["source_t"], old_kernel.get("t"), f"{unknown_id}.source_t/old kernel t")
        _same_exact_interval(unknown["source_t"], new_kernel.get("t"), f"{unknown_id}.source_t/new kernel t")
        _same_exact_interval(old_kernel.get("F1"), new_kernel.get("F1"), f"source cell {source_cell} F kernel mass")
        _same_exact_interval(old_kernel.get("d_x_F1"), new_kernel.get("d_x_F1", new_kernel.get("Fx")), f"source cell {source_cell} Fx kernel mass")
        if "Fx" in new_kernel:
            _same_exact_interval(old_kernel.get("d_x_F1"), new_kernel["Fx"], f"source cell {source_cell} Fx alias kernel mass")
        old_records[unknown_id] = {
            "unknown": unknown,
            "base_key": base_key,
            "matrix_source_paths": paths,
            "profile_terms": profile_terms,
            "profile_cell": profile_cell,
            "old_kernel": old_kernel,
            "new_kernel": new_kernel,
        }
        old_key_owner[base_key] = unknown_id

    # The product table is indexed before binding so omissions and duplicate
    # identities are rejected rather than silently folded.
    product_index: dict[tuple[int, str, int, str, int, str], Mapping[str, Any]] = {}
    product_id_index: dict[str, tuple[int, str, int, str, int, str]] = {}
    product_count_by_field = {field: 0 for field in SEAM_FIELDS}
    for field in SEAM_FIELDS:
        for product_index_in_field, product in enumerate(_seam_terms(seam_result, field)):
            label = f"new seam {field}[{product_index_in_field}]"
            key, term_id, path = _validate_product(product, seam_result, label)
            if key[-1] != field:
                raise ValueError(f"{label}.target_field disagrees with its terms_by_order bucket")
            if term_id in product_id_index:
                raise ValueError(f"duplicate seam product term_id {term_id!r}")
            if key in product_index:
                raise ValueError(f"duplicate seam product binding key {key!r}")
            product_index[key] = product
            product_id_index[term_id] = key
            product_count_by_field[field] += 1

            source_cell, matrix_key, _matrix_term_index, state_field, state_term_index, _target_field = key
            profile_terms = state_index.get((source_cell, state_field))
            if profile_terms is None or state_term_index >= len(profile_terms):
                raise ValueError(f"{label} has no corresponding profile state term index")
            expected_state_id = profile_terms[state_term_index].get("term_id")
            if product.get("state_profile_term_id") != expected_state_id:
                raise ValueError(f"{label}.state_profile_term_id does not match profile target_cells[{source_cell}].W1_Z2[{state_field}]")
            profile_term = profile_terms[state_term_index]
            if profile_term.get("target_cell", source_cell) != source_cell:
                raise ValueError(f"{label} profile target cell differs from outer matrix cell")
            for key_name, state_name in (("state_source_cell", "source_cell"), ("state_source_cell_index", "source_cell_index")):
                if product.get(key_name, profile_term.get(state_name)) != profile_term.get(state_name):
                    raise ValueError(f"{label}.{key_name} differs from internal profile source cell")
            if product.get("state_profile_term_id") not in state_id_owner:
                raise ValueError(f"{label}.state_profile_term_id is not present in the profile")
            _same_exact_interval(product["source_t"], profile_cells[source_cell]["t"], f"{label}.source_t/profile t")
            _same_exact_interval(product["source_t"], _kernel(seam_result, source_cell, "new seam result").get("t"), f"{label}.source_t/kernel t")
            kernel = _kernel(seam_result, source_cell, "new seam result")
            mass_key = "F1" if field == "F1" else "d_x_F1"
            _same_exact_interval(product["kernel_mass"], kernel.get(mass_key, kernel.get("Fx")), f"{label}.kernel_mass")
            # The state term itself is retained as a raw reference in each
            # binding.  Its profile identity is checked above; no coefficient
            # is read or multiplied.
            if not isinstance(profile_term, Mapping):  # defensive for typed callers
                raise TypeError(f"{label} profile state term is not a mapping")
            _validate_low_payload(
                low_output=low_result,
                product=product,
                profile_term=profile_term,
                label=label,
            )

    bindings: dict[str, dict[str, Any]] = {}
    reverse_product_to_unknown: dict[str, str] = {}
    expected_keys: set[tuple[int, str, int, str, int, str]] = set()
    for unknown_id, record in old_records.items():
        unknown = record["unknown"]
        source_cell, matrix_key, matrix_term_index, state_field = record["base_key"]
        profile_terms = record["profile_terms"]
        by_field: dict[str, list[dict[str, Any]]] = {field: [] for field in SEAM_FIELDS}
        ids_by_field: dict[str, list[str]] = {field: [] for field in SEAM_FIELDS}
        state_ids: list[str] = []
        for state_term_index, state_term in enumerate(profile_terms):
            state_id = state_term["term_id"]
            state_ids.append(state_id)
            for target_field in SEAM_FIELDS:
                expected_key = _expected_product_key(record["base_key"], target_field, state_term_index)
                expected_keys.add(expected_key)
                if expected_key not in product_index:
                    raise ValueError(f"missing seam product for old unknown {unknown_id!r}: {expected_key!r}")
                product = product_index[expected_key]
                product_id = product["term_id"]
                if product_id in reverse_product_to_unknown:
                    raise ValueError(f"seam product {product_id!r} would bind more than one old unknown")
                reverse_product_to_unknown[product_id] = unknown_id
                by_field[target_field].append(_clone(product))
                ids_by_field[target_field].append(product_id)

                # The old unknown's LOW source path must be present verbatim in
                # the new product's direct receipt path field.
                _find_path(record["matrix_source_paths"], product["matrix_receipt_input_path"], f"old unknown {unknown_id!r}")
                _same_exact_interval(product["source_t"], unknown["source_t"], f"{product_id}.source_t/old unknown")
                _same_exact_value(product["target_x"], old_seam_x, f"{product_id}.target_x/old seam")
                _same_exact_value(product["target_t"], old_seam_t, f"{product_id}.target_t/old seam")

        for target_field in SEAM_FIELDS:
            if not ids_by_field[target_field]:
                raise ValueError(f"old unknown {unknown_id!r} has no {target_field} products")
        matrix_source = {
            "matrix_name": "B1",
            "matrix_order": "0",
            "matrix_key": matrix_key,
            "row": 5,
            "column": int(matrix_key.split(",")[1]),
            "matrix_term_index": matrix_term_index,
            "source_cell": source_cell,
            "source_cell_index": source_cell,
            "low_matrix_receipt_input_paths": _clone(record["matrix_source_paths"]),
            "matrix_node_ref": _clone(unknown.get("matrix_node_ref")),
        }
        state_source = {
            "source_cell": source_cell,
            "source_cell_index": source_cell,
            "state_field": state_field,
            "state_order": "2",
            "ordinary_Z_order": 1,
            "profile_target_cell_index": source_cell,
            "state_profile_term_ids": list(state_ids),
            "profile_terms": [_clone(term) for term in profile_terms],
        }
        product_ids = ids_by_field["F1"] + ids_by_field["d_x_F1"]
        binding = {
            "old_unknown_id": unknown_id,
            "unknown_id": unknown_id,
            "old_unknown": _clone(unknown),
            "unknown_modified": False,
            "unknown_uncertainty_eliminated": False,
            "matrix_source": matrix_source,
            "state_source": state_source,
            "product_ids_by_field": {field: list(ids_by_field[field]) for field in SEAM_FIELDS},
            "products_by_field": by_field,
            "product_counts_by_field": {field: len(ids_by_field[field]) for field in SEAM_FIELDS},
            "product_count": len(product_ids),
            "product_ids": product_ids,
            "relation": "one_old_unknown_to_many_new_seam_products_by_state_term_and_target_field",
            "one_to_many": len(product_ids) > 1,
            "many_to_one": False,
            "finite_K0_conditional_binding": True,
            "ordinary_Z_order": 1,
            "state_Z_order": 2,
            "row": 5,
            "target_x": _clone(old_seam_x),
            "target_t": _clone(old_seam_t),
            "error_classes_retained": True,
            "raw_directed_coefficients_retained": True,
            "raw_ids_retained": True,
            "no_cancellation": True,
            "cancellation_performed": False,
        }
        bindings[unknown_id] = binding

    if set(product_index) != expected_keys:
        extras = set(product_index) - expected_keys
        missing = expected_keys - set(product_index)
        raise ValueError(f"seam product coverage differs; extras={sorted(extras)!r}, missing={sorted(missing)!r}")
    if set(reverse_product_to_unknown) != set(product_id_index):
        raise ValueError("not every seam product is bound exactly once")

    product_ids_by_unknown = {unknown_id: list(binding["product_ids"]) for unknown_id, binding in bindings.items()}
    total_products = sum(binding["product_count"] for binding in bindings.values())
    return {
        "schema": SCHEMA,
        "input_schemas": {"old": OLD_SCHEMA, "profile": PROFILE_SCHEMA, "seam": SEAM_SCHEMA},
        "low_output_consumed": True,
        "low_matrix_source_binding_validated": True,
        "bindings": list(bindings.values()),
        "product_ids_by_unknown": product_ids_by_unknown,
        "product_id_to_old_unknown_id": reverse_product_to_unknown,
        "old_unknown_count": len(bindings),
        "unknown_count": len(bindings),
        "product_count": total_products,
        "product_counts_by_field": {
            field: sum(binding["product_counts_by_field"][field] for binding in bindings.values())
            for field in SEAM_FIELDS
        },
        "seam_product_counts_by_field": product_count_by_field,
        "all_old_unknowns_bound": True,
        "all_seam_products_bound": True,
        "missing_products": [],
        "duplicate_products": [],
        "duplicate_unknowns": [],
        "relation": {
            "kind": "old_unknown_to_new_seam_term",
            "direction": "one_old_unknown_to_many_new_terms",
            "reverse_product_to_old_unknown": reverse_product_to_unknown,
            "many_to_one_preserved": True,
            "no_product_reused": True,
        },
        "source_t_exact_checked": True,
        "target_scales_exact_checked": True,
        "F_and_Fx_kernel_masses_exact_checked": True,
        "directed_mpf_endpoints_compared_structurally": True,
        "midpoint_used": False,
        "ancestor_recomputed": False,
        "ancestor_called": False,
        "numeric_recomputed": False,
        "products_merged": False,
        "products_eliminated": False,
        "no_cancellation": True,
        "cancellation_performed": False,
        "finite_K0_condition_package": True,
        "conditional_on_same_LOW_matrix_and_profile_sources": True,
        "complete_physical_owner": False,
        "physical_owner_claimed": False,
        "physical_field_claimed": False,
        "complete_secondword": False,
        "complete_firstword_Z2": False,
        "Kchi_retained_untouched": True,
        "tail_retained_untouched": True,
        "old_error_classes_retained": True,
        "old_212_account_retained_untouched": True,
        "P1_is_m15": False,
        "m15_computed": False,
        "retention_scope": {
            "Kchi": "not touched by this finite K0 bridge",
            "tail": "not touched by this finite K0 bridge",
            "old_errors": "retained in old_unknown and raw seam product copies",
            "ledger_212": "not touched by this bridge",
            "P1": "not identified with m15",
        },
    }


build_z2_seam_bridge = build

