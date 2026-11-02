"""Standalone four-cell ordinary ``W0`` ``Z2`` local enclosure, v2.

The source rows are already ordinary ``Z2`` derivatives.  This leaf therefore
keeps each source term and adds only the formal ``hb`` power supplied by the
directed Green-kernel enclosure.  It does not call an earlier leaf, producer,
or helper, and it never materialises a formal atom or a physical measure.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


FIELDS = ("F1", "Uz1", "K1")
ROWS = {"F1": "5", "Uz1": "6", "K1": "6"}
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDER = 2
CELL_COUNT = 4


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _is_singleton(value: Any) -> bool:
    try:
        raw = value._mpi_
    except AttributeError as exc:
        raise ValueError("directed interval required") from exc
    return raw[0] == raw[1]


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((lower._mpi_[0], upper._mpi_[1]))


def _nonnegative(c: Any, upper: Any) -> Any:
    return _interval(c, c.mpf(0), upper)


def _source_info(term: dict[str, Any]) -> tuple[list[Any], list[Any], list[Any], dict[str, Any]]:
    raw_origin = term.get("origin", {})
    if raw_origin is None:
        raw_origin = {}
    if not isinstance(raw_origin, dict):
        raise ValueError("origin must be a mapping")
    origin = deepcopy(raw_origin)

    ids = _as_list(term.get("source_ids")) or _as_list(origin.get("source_ids"))
    if not ids and term.get("input_id") is not None:
        ids = [term["input_id"]]
    if not ids and origin.get("input_id") is not None:
        ids = [origin["input_id"]]

    paths = _as_list(term.get("source_paths")) or _as_list(term.get("paths"))
    if not paths and term.get("source_path") is not None:
        paths = [term["source_path"]]
    if not paths and term.get("path") is not None:
        paths = [term["path"]]
    if not paths:
        paths = _as_list(origin.get("source_paths"))
    if not paths and origin.get("source_path") is not None:
        paths = [origin["source_path"]]
    if not paths and origin.get("path") is not None:
        paths = [origin["path"]]

    refs = term.get("source_node_refs")
    if refs is None:
        refs = origin.get("source_node_refs")
    refs = _as_list(deepcopy(refs))
    if not refs:
        for index in range(max(len(ids), len(paths))):
            node: dict[str, Any] = {}
            if index < len(ids):
                node["source_id"] = deepcopy(ids[index])
            if index < len(paths):
                node["source_path"] = deepcopy(paths[index])
            refs.append(node)
    if not ids and not paths and not refs:
        raise ValueError("source ids, source paths, or source node refs required")
    return ids, paths, refs, origin


def _validate_term(term: Any, label: str) -> dict[str, Any]:
    if not isinstance(term, dict):
        raise ValueError(f"{label} must be a mapping")
    powers = term.get("powers")
    if not isinstance(powers, (list, tuple)) or len(powers) != len(ATOMS):
        raise ValueError(f"{label} must have four formal powers")
    if any(isinstance(power, bool) or not isinstance(power, int) for power in powers):
        raise ValueError(f"{label} formal powers must be integers")
    coefficient = term.get("coefficient")
    try:
        coefficient._mpi_
    except AttributeError as exc:
        raise ValueError(f"{label} coefficient must be a directed interval") from exc
    _source_info(term)
    return term


def _local_term(
    source: dict[str, Any],
    *,
    factor: Any,
    cell: int,
    field: str,
    row: str,
    sign: int,
) -> dict[str, Any]:
    ids, paths, refs, source_origin = _source_info(source)
    out = deepcopy(source)
    powers = list(source["powers"])
    powers[0] += 2
    coefficient = source["coefficient"] * factor
    if sign < 0:
        coefficient = -coefficient

    origin = deepcopy(source_origin)
    origin["source_origin"] = deepcopy(source_origin)
    origin["source_ids"] = deepcopy(ids)
    if paths:
        origin["source_paths"] = deepcopy(paths)
    origin["source_node_refs"] = deepcopy(refs)
    origin["source_operation"] = source_origin.get("operation")
    origin["operation"] = "w0_z2_four_cell_v2_local"
    origin["phase_cell"] = cell
    origin["source_row"] = row
    origin["field"] = field
    origin["kernel_sign"] = sign
    origin["ordinary_Z_order"] = ORDINARY_Z_ORDER
    origin["source_z_derivative_is_ordinary"] = True
    origin["ordinary_Z_factorial_divided"] = False
    origin["source_is_whole_collar_enclosure"] = True
    origin.setdefault("uncertainty", source.get("uncertainty"))
    origin.setdefault("error_classes", deepcopy(source.get("error_classes", [])))

    out["powers"] = powers
    out["coefficient"] = coefficient
    out["origin"] = origin
    out["source_ids"] = deepcopy(ids)
    if paths:
        out["source_paths"] = deepcopy(paths)
        if len(paths) == 1:
            out["source_path"] = deepcopy(paths[0])
    out["source_node_refs"] = deepcopy(refs)
    out["uncertainty"] = deepcopy(source.get("uncertainty", source_origin.get("uncertainty")))
    out["error_classes"] = deepcopy(source.get("error_classes", source_origin.get("error_classes", [])))
    out["local_enclosure"] = {
        "kind": "uniform_source_times_directed_kernel_bound",
        "kernel_bound_without_formal_hb2": factor,
        "formal_hb_power_increment": 2,
        "formal_hb2_materialized": False,
        "whole_collar_source": True,
    }
    return out


def _require_keys(mapping: Any, keys: tuple[str, ...], label: str) -> None:
    if not isinstance(mapping, dict) or set(mapping) != set(keys):
        raise ValueError(f"{label} must contain exactly {keys}")


def build(c: Any, *, sc: Any, core: Any, g: Any) -> dict[str, Any]:
    """Build independent local ``W0 Z2`` enclosures for four q cells."""

    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
        raise ValueError("an interval context with mpf, make_mpf, and exp is required")
    _require_keys(core, FIELDS, "core")
    _require_keys(g, ("5", "6"), "g")
    sc = c.mpf(sc)
    if not _is_singleton(sc) or sc.a <= 0:
        raise ValueError("sc must be a positive directed singleton")

    for field in FIELDS:
        terms = core[field]
        if not isinstance(terms, (list, tuple)) or not terms:
            raise ValueError(f"core[{field}] Z2 source row must be non-empty")
        for index, term in enumerate(terms):
            _validate_term(term, f"core[{field}][{index}]")
    for row in ("5", "6"):
        terms = g[row]
        if not isinstance(terms, (list, tuple)) or not terms:
            raise ValueError(f"g[{row}] Z2 source row must be non-empty")
        for index, term in enumerate(terms):
            _validate_term(term, f"g[{row}][{index}]")

    core_transport = {field: deepcopy(core[field]) for field in FIELDS}
    local: dict[int, dict[str, list[dict[str, Any]]]] = {}
    q_cells: list[dict[str, Any]] = []
    node_refs: dict[str, Any] = {
        "core": {
            field: [{"container": "core_transport", "field": field}]
            for field in FIELDS
        },
        "g": {row: [] for row in ("5", "6")},
    }

    for cell in range(CELL_COUNT):
        q_lower = c.mpf(cell) / 4
        q_upper = c.mpf(cell + 1) / 4
        s_lower = (sc / 4) * q_lower
        s_upper = (sc / 4) * q_upper
        q_hi = q_upper
        s_q = (sc / 4) * q_hi * c.exp(sc / 16)
        bound = _nonnegative(c, (s_q * s_q) / 2)
        cell_data = {
            "index": cell,
            "q_lower": q_lower,
            "q_upper": q_upper,
            "q_hi": q_hi,
            "s_lower": s_lower,
            "s_upper": s_upper,
            "S_q": s_q,
            "S_q_formula": "(sc/4)*q_hi*exp(sc/16)",
            "x_minus_2_bound": {
                "formal_hb_multiplier": 1,
                "interval_without_hb": _nonnegative(c, s_q),
                "expression": "|X-2| <= hb*S_q",
            },
            "kernel_bound": bound,
            "kernel_bound_without_formal_hb2": bound,
            "kernel_bound_is_without_formal_hb2": True,
            "formal_hb2_materialized": False,
            "physical_or_cartesian_measure_included": False,
            "kernel_bound_formula": "[0,S_q^2/2]",
            "kernel_formulas": {
                "F1": "u/2*(1-(u/X)^2)",
                "Uz1": "u*log(X/u)",
                "K1": "-u/2*(1-(u/X)^2)",
            },
            "core_ref": {
                field: {"container": "core_transport", "field": field}
                for field in FIELDS
            },
            "definition": "W0_Z2[field] = core_transport[field] + local[cell][field]",
            "local_source_rows": dict(ROWS),
        }
        q_cells.append(cell_data)
        local[cell] = {field: [] for field in FIELDS}
        for field in FIELDS:
            row = ROWS[field]
            sign = -1 if field == "K1" else 1
            for source_index, source in enumerate(g[row]):
                mapped = _local_term(
                    source,
                    factor=bound,
                    cell=cell,
                    field=field,
                    row=row,
                    sign=sign,
                )
                original_refs = deepcopy(mapped["source_node_refs"])
                mapped["source_node_refs"] = [{
                    "container": "g",
                    "row": row,
                    "term_index": source_index,
                    "source_refs": original_refs,
                }]
                mapped["origin"]["source_node_refs"] = deepcopy(mapped["source_node_refs"])
                local[cell][field].append(mapped)
                node_refs["g"][row].append({
                    "cell": cell,
                    "field": field,
                    "term_index": source_index,
                })

    return {
        "ordinary_Z_order": ORDINARY_Z_ORDER,
        "ordinary_Z_orders": [ORDINARY_Z_ORDER],
        "fields": list(FIELDS),
        "rows": dict(ROWS),
        "source_power_basis": list(ATOMS),
        "core_transport": core_transport,
        "local": local,
        "q_cells": q_cells,
        "source_node_refs": node_refs,
        "source_z2_is_ordinary_derivative": True,
        "ordinary_Z_factorial_divided": False,
        "geometry_Z_independent": True,
        "hb_assumption": "0<hb<1/2",
        "hb_materialized": False,
        "formal_hb2_materialized": False,
        "X_formula": "X=2*exp(hb*s/2)",
        "X_minus_2_bound": "|X-2|<=hb*S_q",
        "kernel_inequalities": {
            "F1": "0<=u/2*(1-(u/X)^2)<=X-u",
            "Uz1": "0<=u*log(X/u)<=X-u",
            "K1": "K1_kernel=-F1_kernel",
        },
        "g_uniform_on_entire_first_collar_required": True,
        "source_z2_whole_collar_enclosure": True,
        "source_terms_are_values": False,
        "source_endpoint_not_selected": True,
        "kernel_bound_includes_formal_hb2": False,
        "physical_or_cartesian_measure_included": False,
        "no_physical_cartesian_measure": True,
        "core_history_reset": False,
        "core_retained_independently_per_cell": True,
        "local_independent_by_origin": True,
        "finite_word_definition": "core_transport[field] + local[cell][field]",
        "local_added_once_in_definition": True,
        "allword_tail_added": False,
        "m15_computed": False,
        "old_W1_consumed": False,
        "same_source_cancellation_claimed": False,
        "uncertainty_preserved": True,
        "error_classes_preserved": True,
        "live_owner_claimed": False,
        "physical_cartesian_frame_claimed": False,
        "uncertainty_correlation_proved": False,
        "external_W0_Z0_Z1_recomputed": False,
    }


__all__ = ["ATOMS", "CELL_COUNT", "FIELDS", "ORDINARY_Z_ORDER", "ROWS", "build"]
