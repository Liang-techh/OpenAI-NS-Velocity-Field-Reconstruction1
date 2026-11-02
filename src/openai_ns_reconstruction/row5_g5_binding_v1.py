"""Exact, source-only binding of the row-5 ``g5`` occurrences.

The affine expression is the authority for occurrence identity.  This leaf
checks each ``h5``/``h5_Z1`` source-only occurrence against one first-word
source-only product and one baseline ``skipped_g5`` record.  It deliberately
does not evaluate intervals, take midpoints, infer ownership, or walk an
arbitrary source graph.  Exact MPF endpoint tuples are compared byte-for-byte
after validating their shape.

``build`` accepts either ``build(context, affine, firstword, baseline)`` (the
usual runner shape) or ``build(affine, firstword, baseline)``.  The context is
accepted for API compatibility and is not used.  The three payloads may also
be direct result mappings or explicit wrappers containing ``expression``,
``result``/``products`` and optional receipt metadata.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "row5_g5_binding_v1"
COMPONENT = "row5_g5_binding_v1"
SOURCE_SCHEMA = "row5_affine_forcing_v1"
FIRSTWORD_SCHEMA = "row5_firstword_product_v4"
BASELINE_SCHEMA = "row5_baseline_product_v4"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_ROWS = ("h5", "h5_Z1")
ROW_ORDERS = {"h5": 0, "h5_Z1": 1}
EXPECTED_G5_PER_ROW = 16

_MISSING = object()
_SEQUENCE_TYPES = (str, bytes, bytearray)
_INTERVAL_ALIASES = (
    "coefficient_interval",
    "coefficient",
    "source_interval",
    "coefficient_enclosure",
    "coefficient_endpoints",
    "source_interval_endpoints",
    "directed_interval",
    "enclosure",
    "interval",
)
_IDENTITY_FIELDS = (
    "source_family",
    "source_key",
    "matrix_name",
    "role",
    "term_role",
    "source_ordinary_Z_order",
    "ordinary_Z_order",
)
_REFERENCE_FIELDS = ("state_reference", "state_ref", "state", "reference")
_RECEIPT_KEYS = (
    "receipt_path",
    "source_receipt_path",
    "receipt",
    "source_receipt",
    "input_receipt_path",
)

__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_SCHEMA",
    "FIRSTWORD_SCHEMA",
    "BASELINE_SCHEMA",
    "SOURCE_POWER_BASIS",
    "OUTPUT_ROWS",
    "EXPECTED_G5_PER_ROW",
    "build",
    "build_binding",
    "bind",
    "build_g5_binding",
    "row5_g5_binding",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, _SEQUENCE_TYPES)


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def _mapping(value: Any, label: str) -> Mapping[Any, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def _sequence(value: Any, label: str, *, nonempty: bool = False) -> list[Any]:
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a sequence")
    result = list(value)
    if nonempty and not result:
        raise ValueError(f"{label} must not be empty")
    return result


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if value not in result:
            result.append(_clone(value))
    return result


def _unwrap(value: Any, expected_schema: str, names: Sequence[str], label: str) -> tuple[Mapping[Any, Any], Mapping[Any, Any]]:
    """Return (payload, outer metadata) from direct runner envelope paths.

    These are intentionally explicit paths, not a recursive search: affine is
    normally ``result.output.expression``, first-word is
    ``result.output.result``, and baseline is ``output.products`` after the
    runner has decoded its saved ``result.json`` record.
    """

    outer = _mapping(value, label)
    if outer.get("schema") == expected_schema:
        return outer, outer
    paths: list[tuple[str, ...]] = []
    for name in names:
        paths.extend(
            (
                (name,),
                ("output", name),
                ("result", name),
                ("result", "output", name),
            )
        )
    # A full saved receipt uses result.output.result for first-word and
    # result.output.expression for affine; the named paths above cover both.
    # Keep the explicit products path for the baseline output envelope.
    if "products" in names:
        paths.extend((("result", "output", "products"), ("output", "products")))
    for path in paths:
        nested: Any = outer
        for part in path:
            if not isinstance(nested, Mapping):
                nested = _MISSING
                break
            nested = nested.get(part, _MISSING)
        if isinstance(nested, Mapping) and nested.get("schema") == expected_schema:
            return nested, outer
    raise ValueError(f"{label} does not contain schema {expected_schema}")


def _receipt(outer: Mapping[Any, Any], payload: Mapping[Any, Any]) -> Any:
    for source in (outer, payload):
        for key in _RECEIPT_KEYS:
            if key in source and source[key] is not None:
                # Preserve the complete caller-owned receipt descriptor
                # (path, raw SHA, JSON prefix, encoding) so a cell reference
                # remains uniquely resolvable.
                return _clone(source[key])
    return None


def _receipt_prefix(receipt: Any) -> Any:
    if isinstance(receipt, Mapping) and "json_path_prefix" in receipt:
        return _clone(receipt["json_path_prefix"])
    return None


def _rows(payload: Mapping[Any, Any], label: str) -> dict[str, list[Any]]:
    rows = _mapping(payload.get("rows", _MISSING), label + ".rows")
    result: dict[str, list[Any]] = {}
    for row in OUTPUT_ROWS:
        result[row] = _sequence(rows.get(row, _MISSING), label + f".rows.{row}")
    return result


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label} must be a four-entry exact MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} must contain integers only")
    sign, mantissa, exponent, bitcount = (int(item) for item in value)
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(f"{label} has invalid exact MPF tuple fields")
    if mantissa == 0 and bitcount != 0:
        raise ValueError(f"{label} is non-finite")
    if mantissa and (mantissa.bit_length() != bitcount or not mantissa & 1):
        raise ValueError(f"{label} is not a canonical finite MPF tuple")
    return sign, mantissa, exponent, bitcount


def _compare_exact_mpf(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    # Compare dyadic endpoints without converting the potentially enormous
    # exponent to a float or constructing integers of astronomical size.
    ls, lm, le, lb = left
    rs, rm, re, rb = right
    if not lm and not rm:
        return 0
    if not lm:
        return 1 if rs else -1
    if not rm:
        return -1 if ls else 1
    if ls != rs:
        return -1 if ls else 1
    lt, rt = le + lb, re + rb
    if lt != rt:
        magnitude = -1 if lt < rt else 1
    else:
        bottom = min(le, re)
        lv, rv = lm << (le - bottom), rm << (re - bottom)
        magnitude = (lv > rv) - (lv < rv)
    return -magnitude if ls else magnitude


def _raw_endpoint(value: Any, upper: bool, label: str) -> tuple[int, int, int, int]:
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            key = "upper_exact_mpf_tuple" if upper else "lower_exact_mpf_tuple"
            if key not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            return _exact_tuple(value[key], label + "." + key)
        # First-word node endpoint payloads use lower/upper; they are still
        # required to contain exact MPF tuples, never decimal bounds.
        key = "upper" if upper else "lower"
        if key in value:
            return _exact_tuple(value[key], label + "." + key)
        key = "hi" if upper else "lo"
        if key in value:
            return _exact_tuple(value[key], label + "." + key)
        key = "right" if upper else "left"
        if key in value:
            return _exact_tuple(value[key], label + "." + key)
        for alias in _INTERVAL_ALIASES:
            if alias in value and value[alias] is not value:
                return _raw_endpoint(value[alias], upper, label + "." + alias)
        raise TypeError(f"{label} is not an exact directed interval")
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        if not _is_sequence(raw) or len(raw) != 2:
            raise TypeError(f"{label} has malformed MP interval")
        return _exact_tuple(raw[1 if upper else 0], label + (".upper" if upper else ".lower"))
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return _exact_tuple(raw, label + ".point")
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return _exact_tuple(nested, label + (".upper" if upper else ".lower"))
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None:
            return _exact_tuple(nested[1 if upper else 0], label + (".upper" if upper else ".lower"))
    raise TypeError(f"{label} is not an exact directed interval")


def _endpoints(value: Any, label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    lower = _raw_endpoint(value, False, label)
    upper = _raw_endpoint(value, True, label)
    if _compare_exact_mpf(lower, upper) > 0:
        raise ValueError(f"{label} has reversed directed endpoints")
    return lower, upper


def _interval_from_record(record: Mapping[Any, Any], label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    found = False
    expected: tuple[tuple[int, int, int, int], tuple[int, int, int, int]] | None = None
    for key in _INTERVAL_ALIASES:
        if key not in record:
            continue
        current = _endpoints(record[key], label + "." + key)
        if expected is None:
            expected = current
        elif current != expected:
            raise ValueError(f"{label} has inconsistent interval aliases")
        found = True
    if not found or expected is None:
        raise KeyError(f"{label} has no exact directed interval")
    return expected


def _powers(record: Mapping[Any, Any], label: str) -> list[int]:
    if "formal_powers" not in record:
        raise KeyError(f"{label}.formal_powers is required")
    value = record["formal_powers"]
    if isinstance(value, Mapping):
        value = [value.get(name, _MISSING) for name in SOURCE_POWER_BASIS]
    values = _sequence(value, label + ".formal_powers")
    if len(values) != 4 or any(isinstance(item, bool) or type(item) is not int for item in values):
        raise TypeError(f"{label}.formal_powers must contain four integers")
    result = [int(item) for item in values]
    if "powers" in record:
        power_alias = record["powers"]
        if isinstance(power_alias, Mapping):
            power_alias = [power_alias.get(name, _MISSING) for name in SOURCE_POWER_BASIS]
        alias_values = _sequence(power_alias, label + ".powers")
        if len(alias_values) != 4 or any(isinstance(item, bool) or type(item) is not int for item in alias_values):
            raise TypeError(f"{label}.powers must contain four integers")
        if [int(item) for item in alias_values] != result:
            raise ValueError(f"{label}.powers disagrees with formal_powers")
    return result


def _scalar(record: Mapping[Any, Any], label: str) -> int:
    value = record.get("scalar", _MISSING)
    if value is _MISSING or isinstance(value, bool) or type(value) is not int:
        raise TypeError(f"{label}.scalar must be an integer")
    return int(value)


def _sequence_field(record: Mapping[Any, Any], names: Sequence[str], label: str) -> list[Any]:
    for name in names:
        if name in record:
            value = record[name]
            if name.endswith("_ids") or name.endswith("_paths"):
                return _sequence(value, label + "." + name)
            return [value]
    raise KeyError(f"{label} is missing {names[0]}")


def _identity(record: Mapping[Any, Any], label: str) -> dict[str, Any]:
    ids = _unique(_sequence_field(record, ("source_ids", "provenance_ids", "parent_source_ids", "source_id"), label))
    paths = _unique(_sequence_field(record, ("source_paths", "parent_input_paths", "source_path", "path"), label))
    return {"source_ids": ids, "source_paths": paths}


def _check_identity(expected: Mapping[Any, Any], record: Mapping[Any, Any], label: str) -> None:
    expected_identity = _identity(expected, "affine source")
    actual_identity = _identity(record, label)
    if actual_identity != expected_identity:
        raise ValueError(f"{label} source identity differs")
    # Do not let a normalized alias hide a changed node identity.  Repeated
    # IDs/paths introduced by provenance wrappers are harmless; a distinct
    # alias value is not.
    for key in ("source_ids", "provenance_ids", "parent_source_ids", "source_id"):
        if key in record:
            if key == "source_id":
                if record[key] != expected_identity["source_ids"][0]:
                    raise ValueError(f"{label}.{key} differs")
            elif _unique(_sequence(record[key], label + "." + key)) != expected_identity["source_ids"]:
                raise ValueError(f"{label}.{key} differs")
    for key in ("source_paths", "parent_input_paths", "source_path", "path"):
        if key in record:
            if key.endswith("_paths"):
                values = _sequence(record[key], label + "." + key)
                if _unique(values) != expected_identity["source_paths"]:
                    raise ValueError(f"{label}.{key} differs")
            elif record[key] not in expected_identity["source_paths"]:
                raise ValueError(f"{label}.{key} is not an affine source path")
    for key in _IDENTITY_FIELDS:
        if key in expected and record.get(key, _MISSING) != expected[key]:
            raise ValueError(f"{label}.{key} differs from affine occurrence")


def _check_raw_identity(expected: Mapping[Any, Any], original: Mapping[Any, Any], label: str) -> None:
    """Require the original source IDs/paths, including order and repeats."""

    for key in ("source_ids", "source_paths"):
        if key not in expected:
            raise KeyError(f"affine source is missing {key}")
        expected_value = _sequence(expected[key], "affine source." + key)
        actual_value = _sequence(original.get(key, _MISSING), label + "." + key)
        if actual_value != expected_value:
            raise ValueError(f"{label}.{key} is not preserved itemwise")


def _check_state_free(record: Mapping[Any, Any], label: str) -> None:
    for key in _REFERENCE_FIELDS:
        if key in record and record[key] is not None:
            raise ValueError(f"{label}.{key} is not null for source-only g5")


def _node_original(node: Mapping[Any, Any], label: str, *, firstword: bool) -> Mapping[Any, Any]:
    key = "original_node" if firstword else "original_term"
    value = node.get(key, _MISSING)
    if not isinstance(value, Mapping):
        raise KeyError(f"{label}.{key} is required for source identity")
    return value


def _validate_affine_occurrences(affine: Mapping[Any, Any]) -> dict[str, list[dict[str, Any]]]:
    rows = _rows(affine, "affine")
    occurrences: dict[str, list[dict[str, Any]]] = {row: [] for row in OUTPUT_ROWS}
    for row in OUTPUT_ROWS:
        expected_order = ROW_ORDERS[row]
        for index, raw in enumerate(rows[row]):
            term = _mapping(raw, f"affine.rows.{row}[{index}]")
            family = term.get("source_family", term.get("family", _MISSING))
            marker = family == "g5" or term.get("source_only") is True or term.get("source_only_g5") is True
            if not marker:
                continue
            if family != "g5" or term.get("source_only") is not True or term.get("source_only_g5") is not True:
                raise ValueError(f"affine.rows.{row}[{index}] is not an explicit source-only g5 term")
            _check_state_free(term, f"affine.rows.{row}[{index}]")
            if term.get("ordinary_Z_order", _MISSING) != expected_order:
                raise ValueError(f"affine.rows.{row}[{index}] has wrong ordinary Z order")
            if term.get("source_ordinary_Z_order", _MISSING) != expected_order:
                raise ValueError(f"affine.rows.{row}[{index}] has wrong source ordinary Z order")
            _identity(term, f"affine.rows.{row}[{index}]")
            interval = _interval_from_record(term, f"affine.rows.{row}[{index}]")
            powers = _powers(term, f"affine.rows.{row}[{index}]")
            scalar = _scalar(term, f"affine.rows.{row}[{index}]")
            if scalar != 1:
                raise ValueError(f"affine.rows.{row}[{index}] source-only scalar must be 1")
            occurrences[row].append(
                {
                    "row": row,
                    "ordinary_Z_order": expected_order,
                    "source_index": index,
                    "node_id": f"expression:{row}:{index}",
                    "term": term,
                    "interval": interval,
                    "formal_powers": powers,
                    "scalar": scalar,
                }
            )
        if len(occurrences[row]) != EXPECTED_G5_PER_ROW:
            raise ValueError(
                f"affine row {row} must contain exactly {EXPECTED_G5_PER_ROW} source-only g5 occurrences; "
                f"got {len(occurrences[row])}"
            )
    return occurrences


def _node_table(payload: Mapping[Any, Any], label: str) -> Mapping[Any, Any]:
    table = payload.get("source_nodes", _MISSING)
    if not isinstance(table, Mapping):
        raise KeyError(f"{label}.source_nodes is required")
    return table


def _firstword_products(payload: Mapping[Any, Any], row: str, label: str) -> list[Mapping[Any, Any]]:
    rows = _rows(payload, label)
    values: list[Mapping[Any, Any]] = []
    for index, product in enumerate(rows[row]):
        values.append(_mapping(product, f"{label}.rows.{row}[{index}]") )
    return values


def _baseline_rows(payload: Mapping[Any, Any], label: str) -> Mapping[Any, Any]:
    rows = _mapping(payload.get("rows", _MISSING), label + ".rows")
    for row in OUTPUT_ROWS:
        row_payload = rows.get(row, _MISSING)
        if not isinstance(row_payload, Mapping) or not row_payload:
            raise ValueError(f"{label}.rows.{row} must be a non-empty cell mapping")
        for cell, cell_payload in row_payload.items():
            sectors = _mapping(cell_payload, f"{label}.rows.{row}[{cell!r}]")
            for sector in ("core", "local"):
                _sequence(sectors.get(sector, _MISSING), f"{label}.rows.{row}[{cell!r}].{sector}")
    return rows


def _find_one(values: Sequence[Mapping[Any, Any]], predicate: Any, label: str) -> tuple[int, Mapping[Any, Any]]:
    matches = [(index, value) for index, value in enumerate(values) if predicate(value)]
    if len(matches) != 1:
        raise ValueError(f"{label} requires exactly one match; got {len(matches)}")
    return matches[0]


def _check_source_node(
    affine: Mapping[Any, Any],
    node: Mapping[Any, Any],
    occurrence: Mapping[str, Any],
    label: str,
    *,
    firstword: bool,
) -> None:
    expected_node_id = occurrence["node_id"]
    if node.get("node_id", _MISSING) != expected_node_id:
        raise ValueError(f"{label}.node_id differs")
    _check_identity(affine, node, label)
    original = _node_original(node, label, firstword=firstword)
    _check_identity(affine, original, label + ".original")
    _check_raw_identity(affine, original, label + ".original")
    for key in _IDENTITY_FIELDS:
        if key in affine and original.get(key, _MISSING) != affine[key]:
            raise ValueError(f"{label}.original.{key} differs")
    for key in ("uncertainty", "error_only", "symmetric_error", "signed_field", "error_classes", "source_error_classes", "power_basis", "source_only", "source_only_g5", "origin"):
        if key in affine and original.get(key, _MISSING) != affine[key]:
            raise ValueError(f"{label}.original.{key} differs")
    _check_state_free(node, label)
    _check_state_free(original, label + ".original")
    if _powers(original, label + ".original") != occurrence["formal_powers"]:
        raise ValueError(f"{label}.original.formal_powers differs")
    if _scalar(original, label + ".original") != occurrence["scalar"]:
        raise ValueError(f"{label}.original.scalar differs")
    if _interval_from_record(original, label + ".original") != occurrence["interval"]:
        raise ValueError(f"{label}.original.interval differs")
    if _powers(node, label) != occurrence["formal_powers"]:
        raise ValueError(f"{label}.formal_powers differs")
    if _scalar(node, label) != occurrence["scalar"]:
        raise ValueError(f"{label}.scalar differs")
    if _interval_from_record(node, label) != occurrence["interval"]:
        raise ValueError(f"{label}.interval differs")
    for key in ("source_family", "source_key", "role", "term_role", "ordinary_Z_order", "source_ordinary_Z_order"):
        if key in affine and node.get(key, _MISSING) != affine[key]:
            raise ValueError(f"{label}.{key} differs")
    if node.get("function_value", False) is not False:
        raise ValueError(f"{label}.function_value must be false")


def _check_firstword(
    affine: Mapping[Any, Any],
    payload: Mapping[Any, Any],
    occurrences: Mapping[str, Sequence[Mapping[str, Any]]],
) -> dict[tuple[str, int], dict[str, Any]]:
    table = _node_table(payload, "firstword")
    all_products = {row: _firstword_products(payload, row, "firstword") for row in OUTPUT_ROWS}
    expected_ids = {item["node_id"] for row in OUTPUT_ROWS for item in occurrences[row]}
    actual_g5_ids = {
        str(node_id)
        for node_id, node in table.items()
        if isinstance(node, Mapping) and node.get("source_family") == "g5"
    }
    if actual_g5_ids != expected_ids:
        raise ValueError("firstword g5 source-node table does not exactly match affine occurrences")
    expected_keys = {item["node_id"]: (item["row"], int(item["source_index"])) for row in OUTPUT_ROWS for item in occurrences[row]}
    for row in OUTPUT_ROWS:
        for product in all_products[row]:
            source_node_id = product.get("source_node_id", _MISSING)
            if product.get("source_only") is True and source_node_id not in expected_ids:
                raise ValueError("firstword has an unbound source-only product")
            if source_node_id in expected_ids:
                expected_row, expected_index = expected_keys[source_node_id]
                expected_product = f"product:{expected_row}:{expected_index}:source-only"
                if (product.get("source_only") is not True or row != expected_row
                        or product.get("row") != expected_row or product.get("source_index") != expected_index
                        or product.get("product_id") != expected_product):
                    raise ValueError("firstword g5 source has unexpected row, index, product id or feedback reuse")
    result: dict[tuple[str, int], dict[str, Any]] = {}
    for occurrence in (item for row in OUTPUT_ROWS for item in occurrences[row]):
        row = occurrence["row"]
        index = int(occurrence["source_index"])
        node_id = occurrence["node_id"]
        node = table.get(node_id, _MISSING)
        if not isinstance(node, Mapping):
            raise KeyError(f"firstword.source_nodes[{node_id!r}] is missing")
        _check_source_node(occurrence["term"], node, occurrence, f"firstword.source_nodes[{node_id}]", firstword=True)
        product_index, product = _find_one(
            all_products[row],
            lambda item, node_id=node_id, index=index, row=row: item.get("source_node_id") == node_id
            and item.get("source_index", _MISSING) == index
            and item.get("row", _MISSING) == row
            and item.get("source_only") is True,
            f"firstword.rows.{row} source-only occurrence {index}",
        )
        if product.get("state_node_id", _MISSING) is not None or product.get("state_index", _MISSING) is not None:
            raise ValueError(f"firstword.rows.{row}[{product_index}] unexpectedly has state provenance")
        _check_state_free(product, f"firstword.rows.{row}[{product_index}]")
        if product.get("source_node_id") != node_id or product.get("source_index") != index:
            raise ValueError("firstword product source reference differs")
        if product.get("scalar", _MISSING) != occurrence["scalar"]:
            raise ValueError("firstword source-only product scalar differs")
        if product.get("scalar_applied_to_coefficient") is not False:
            raise ValueError("firstword source-only product scalar was applied")
        if _powers(product, f"firstword.rows.{row}[{product_index}]") != occurrence["formal_powers"]:
            raise ValueError("firstword source-only product formal powers differ")
        if _endpoints(product.get("coefficient_enclosure", _MISSING), f"firstword.rows.{row}[{product_index}].coefficient_enclosure") != occurrence["interval"]:
            raise ValueError("firstword source-only product interval differs")
        enclosure = _mapping(product.get("coefficient_enclosure", _MISSING), f"firstword.rows.{row}[{product_index}].coefficient_enclosure")
        if enclosure.get("directed") is not True or enclosure.get("function_value") is not False:
            raise ValueError("firstword source-only product is not a directed enclosure")
        if product.get("coefficient_is_function_value") is not False:
            raise ValueError("firstword source-only product function-value flag changed")
        result[(row, index)] = {"node_id": node_id, "product_index": product_index, "product": product}
    return result


def _check_baseline(
    affine: Mapping[Any, Any],
    payload: Mapping[Any, Any],
    occurrences: Mapping[str, Sequence[Mapping[str, Any]]],
) -> dict[tuple[str, int], dict[str, Any]]:
    table = _node_table(payload, "baseline")
    rows = _baseline_rows(payload, "baseline")
    skipped = _sequence(payload.get("skipped_g5", _MISSING), "baseline.skipped_g5")
    expected_ids = {item["node_id"] for row in OUTPUT_ROWS for item in occurrences[row]}
    actual_g5_ids = {
        str(node_id)
        for node_id, node in table.items()
        if isinstance(node, Mapping) and node.get("source_family") == "g5"
    }
    if actual_g5_ids != expected_ids:
        raise ValueError("baseline g5 source-node table does not exactly match affine occurrences")
    all_products: list[Mapping[Any, Any]] = []
    for row in OUTPUT_ROWS:
        for cell, cell_payload in rows[row].items():
            for sector in ("core", "local"):
                products = _sequence(cell_payload[sector], f"baseline.rows.{row}[{cell!r}].{sector}")
                for product_index, raw_product in enumerate(products):
                    product = _mapping(raw_product, f"baseline.rows.{row}[{cell!r}].{sector}[{product_index}]")
                    all_products.append(product)
                    source_node_id = product.get("source_node_id", _MISSING)
                    if source_node_id is _MISSING:
                        raise KeyError(f"baseline.rows.{row}[{cell!r}].{sector}[{product_index}] has no source_node_id")
                    if source_node_id in expected_ids:
                        raise ValueError(f"baseline emitted a product for skipped g5 source {source_node_id}")
                    if source_node_id in actual_g5_ids:
                        raise ValueError(f"baseline emitted a product for a g5 source node {source_node_id}")
    result: dict[tuple[str, int], dict[str, Any]] = {}
    for occurrence in (item for row in OUTPUT_ROWS for item in occurrences[row]):
        row = occurrence["row"]
        index = int(occurrence["source_index"])
        node_id = occurrence["node_id"]
        node = table.get(node_id, _MISSING)
        if not isinstance(node, Mapping):
            raise KeyError(f"baseline.source_nodes[{node_id!r}] is missing")
        _check_source_node(occurrence["term"], node, occurrence, f"baseline.source_nodes[{node_id}]", firstword=False)
        skip_index, skip = _find_one(
            skipped,
            lambda item, node_id=node_id, index=index, row=row: item.get("source_node_id") == node_id
            and item.get("source_index", _MISSING) == index
            and item.get("row", _MISSING) == row,
            f"baseline.skipped_g5 {row} occurrence {index}",
        )
        if skip.get("source_family") != "g5":
            raise ValueError(f"baseline.skipped_g5[{skip_index}] source_family differs")
        if skip.get("product_emitted") is not False:
            raise ValueError(f"baseline.skipped_g5[{skip_index}] product_emitted is not false")
        if skip.get("retained_by") != "row5_firstword_product_v4":
            raise ValueError(f"baseline.skipped_g5[{skip_index}] retained_by differs")
        result[(row, index)] = {"node_id": node_id, "skip_index": skip_index, "skip": skip}
    # No extra skipped g5 record may be hidden behind a duplicate identity.
    if len(skipped) != len(expected_ids):
        raise ValueError(f"baseline.skipped_g5 must contain exactly {len(expected_ids)} records")
    return result


def _cell_combinations(
    firstword: Mapping[Any, Any],
    firstword_outer: Mapping[Any, Any],
    baseline: Mapping[Any, Any],
    baseline_outer: Mapping[Any, Any],
) -> list[dict[str, Any]]:
    rows = _baseline_rows(baseline, "baseline")
    h5_cells = list(rows["h5"].keys())
    z1_cells = list(rows["h5_Z1"].keys())
    if h5_cells != z1_cells and set(h5_cells) != set(z1_cells):
        raise ValueError("baseline h5 and h5_Z1 cell sets differ")
    cells = h5_cells
    firstword_receipt = _receipt(firstword_outer, firstword)
    baseline_receipt = _receipt(baseline_outer, baseline)
    combinations: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(cells):
        refs = [
            {
                "branch": "firstword",
                "row": "h5",
                "reference_count": 1,
                "receipt": _clone(firstword_receipt),
                "json_path_prefix": _receipt_prefix(firstword_receipt),
                "json_path": ["rows", "h5"],
            },
            {
                "branch": "firstword",
                "row": "h5_Z1",
                "reference_count": 1,
                "receipt": _clone(firstword_receipt),
                "json_path_prefix": _receipt_prefix(firstword_receipt),
                "json_path": ["rows", "h5_Z1"],
            },
        ]
        for row in OUTPUT_ROWS:
            for sector in ("core", "local"):
                refs.append(
                    {
                        "branch": "baseline",
                        "row": row,
                        "cell": _clone(cell),
                        "sector": sector,
                        "reference_count": 1,
                        "receipt": _clone(baseline_receipt),
                        "json_path_prefix": _receipt_prefix(baseline_receipt),
                        "json_path": ["rows", row, str(cell), sector],
                    }
                )
        combinations.append(
            {
                "cell": _clone(cell),
                "cell_key": str(cell),
                "cell_index": cell_index,
                "references": refs,
                "firstword_scope": "full-collar rows; shared W branch is retained once per row reference",
                "baseline_scope": "this cell only; core and local sectors remain separate",
                "R_hi": {
                    "symbol": "R_hi",
                    "known": False,
                    "preserved": True,
                    "included_in_products": False,
                    "added_to_ledger": False,
                    "receipt": _clone(baseline_receipt),
                    "json_path_prefix": _receipt_prefix(baseline_receipt),
                    "json_path": ["R_hi"],
                },
            }
        )
    return combinations


def _residuals(firstword: Mapping[Any, Any], baseline: Mapping[Any, Any]) -> dict[str, Any]:
    firstword_residual = firstword.get("unknown_residual", _MISSING)
    if not isinstance(firstword_residual, Mapping):
        raise KeyError("firstword.unknown_residual is required")
    if firstword_residual.get("symbol") != "R" or firstword_residual.get("preserved") is not True:
        raise ValueError("firstword existing residual is not retained")
    if firstword_residual.get("included_in_products") is not False or firstword_residual.get("added_to_ledger") is not False:
        raise ValueError("firstword residual entered products or ledger")
    baseline_residual = baseline.get("R_hi", baseline.get("unknown_residual", _MISSING))
    if not isinstance(baseline_residual, Mapping):
        raise KeyError("baseline.R_hi is required")
    if baseline_residual.get("symbol") != "R_hi" or baseline_residual.get("known") is not False:
        raise ValueError("baseline R_hi is not explicitly unknown")
    if baseline_residual.get("preserved") is not True or baseline_residual.get("included_in_products") is not False:
        raise ValueError("baseline R_hi was not preserved outside products")
    return {
        "firstword_existing_W": {
            "symbol": "W",
            "preserved": True,
            "known": False,
            "materialized": False,
            "definition": "existing W retained by firstword branch; no value supplied here",
            "source": "firstword",
        },
        "firstword_unknown_residual": _clone(firstword_residual),
        "R_hi": _clone(baseline_residual),
        "conditional_relation": "R_firstword = W0_baseline + R_hi when both branches use the same W and W_firstword",
        "firstword_residual_is_additional_combined_term": False,
        "same_owner_premise_proved": False,
    }


def _split_arguments(first: Any, second: Any, third: Any, fourth: Any) -> tuple[Any, Mapping[Any, Any], Mapping[Any, Any], Mapping[Any, Any]]:
    if fourth is _MISSING:
        # No context: build(affine, firstword, baseline).
        affine, firstword, baseline = first, second, third
        context = None
    else:
        # Runner-compatible: build(context, affine, firstword, baseline).
        context, affine, firstword, baseline = first, second, third, fourth
    return context, _mapping(affine, "affine"), _mapping(firstword, "firstword"), _mapping(baseline, "baseline")


def build(first: Any, second: Any, third: Any, fourth: Any = _MISSING) -> dict[str, Any]:
    """Bind every source-only g5 occurrence exactly once across three outputs."""

    _context, affine_outer_candidate, firstword_outer_candidate, baseline_outer_candidate = _split_arguments(first, second, third, fourth)
    affine, affine_outer = _unwrap(affine_outer_candidate, SOURCE_SCHEMA, ("expression", "affine", "payload"), "affine")
    firstword, firstword_outer = _unwrap(firstword_outer_candidate, FIRSTWORD_SCHEMA, ("result", "product", "payload"), "firstword")
    baseline, baseline_outer = _unwrap(baseline_outer_candidate, BASELINE_SCHEMA, ("products", "result", "product", "payload"), "baseline")
    if tuple(affine.get("source_power_basis", ())) != SOURCE_POWER_BASIS:
        raise ValueError("affine source power basis differs")
    if tuple(firstword.get("source_power_basis", ())) != SOURCE_POWER_BASIS:
        raise ValueError("firstword source power basis differs")
    if tuple(baseline.get("source_power_basis", ())) != SOURCE_POWER_BASIS:
        raise ValueError("baseline source power basis differs")
    occurrences = _validate_affine_occurrences(affine)
    firstword_matches = _check_firstword(affine, firstword, occurrences)
    baseline_matches = _check_baseline(affine, baseline, occurrences)
    affine_receipt = _receipt(affine_outer, affine)
    firstword_receipt = _receipt(firstword_outer, firstword)
    baseline_receipt = _receipt(baseline_outer, baseline)
    mappings: list[dict[str, Any]] = []
    for row in OUTPUT_ROWS:
        for occurrence in occurrences[row]:
            index = int(occurrence["source_index"])
            key = (row, index)
            first_match = firstword_matches[key]
            baseline_match = baseline_matches[key]
            mappings.append(
                {
                    "mapping_id": f"g5:{row}:{index}",
                    "row": row,
                    "ordinary_Z_order": int(occurrence["ordinary_Z_order"]),
                    "source_index": index,
                    "source_node_references": {
                        "affine": {"row": row, "source_index": index, "node_id": occurrence["node_id"], "receipt": _clone(affine_receipt), "json_path": ["rows", row, index]},
                        "firstword": {"node_id": first_match["node_id"], "product_index": first_match["product_index"], "receipt": _clone(firstword_receipt), "json_path": ["rows", row, first_match["product_index"]]},
                        "baseline": {"node_id": baseline_match["node_id"], "skipped_index": baseline_match["skip_index"], "receipt": _clone(baseline_receipt), "json_path": ["skipped_g5", baseline_match["skip_index"]]},
                    },
                    "affine_source_node_id": occurrence["node_id"],
                    "firstword_source_node_id": first_match["node_id"],
                    "firstword_product_index": first_match["product_index"],
                    "firstword_product_id": first_match["product"].get("product_id"),
                    "baseline_source_node_id": baseline_match["node_id"],
                    "baseline_skipped_index": baseline_match["skip_index"],
                    "baseline_skipped_path": ["skipped_g5", baseline_match["skip_index"]],
                    "source_identity_checked": True,
                    "formal_powers": list(occurrence["formal_powers"]),
                    "scalar": int(occurrence["scalar"]),
                    "directed_exact_interval_endpoints": {
                        "lower_exact_mpf_tuple": list(occurrence["interval"][0]),
                        "upper_exact_mpf_tuple": list(occurrence["interval"][1]),
                        "directed": True,
                    },
                    "ordinary_Z_checked": True,
                    "firstword_source_only_once": True,
                    "baseline_skipped_once": True,
                    "baseline_product_emitted": False,
                }
            )
    if len(mappings) != 2 * EXPECTED_G5_PER_ROW:
        raise ValueError("binding did not retain exactly 32 mappings")
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "exact_source_only_g5_bound",
        "source_schema": SOURCE_SCHEMA,
        "firstword_schema": FIRSTWORD_SCHEMA,
        "baseline_schema": BASELINE_SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": {row: [item for item in mappings if item["row"] == row] for row in OUTPUT_ROWS},
        "mappings": mappings,
        "mapping_count": len(mappings),
        "source_only_g5_per_row": {row: len(occurrences[row]) for row in OUTPUT_ROWS},
        "exact_once": {
            "affine_occurrence_to_firstword_source_only_output": True,
            "affine_occurrence_to_baseline_skipped_record": True,
            "baseline_skipped_has_no_product": True,
            "duplicate_or_missing_fails_closed": True,
        },
        "cell_combinations": _cell_combinations(firstword, firstword_outer, baseline, baseline_outer),
        "unknowns": _residuals(firstword, baseline),
        "conditions": {
            "physical_ownership_validated": False,
            "same_owner_validated": False,
            "same_frame_validated": False,
            "same_scale_validated": False,
            "same_gauge_validated": False,
            "interval_overlap_used": False,
            "midpoint_used": False,
            "source_graph_traversal": False,
            "physical_h5_materialized": False,
            "interval_sum_performed": False,
            "formal_atoms_materialized": False,
        },
        "scope": {
            "firstword_is_full_collar_shared_W_branch": True,
            "baseline_is_cell_indexed_W0_branch": True,
            "core_local_rows_referenced_separately": True,
            "R_hi_unknown_preserved": True,
            "no_product_lists_copied": True,
        },
    }


def bind(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(*args, **kwargs)


def build_binding(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(*args, **kwargs)


def build_g5_binding(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(*args, **kwargs)


def row5_g5_binding(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(*args, **kwargs)

