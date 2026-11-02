"""Standalone four-cell ordinary ``W0`` ``Z2`` local enclosure.

This module is deliberately independent of the earlier collar leaves.  The
input rows are already ordinary ``Z2`` source rows, so the local construction
only adds the two formal ``hb`` powers supplied by the Green-kernel bound.  It
does not evaluate a formal atom, perform quadrature, merge source terms, or
make a claim about a physical Cartesian frame.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterable


FIELDS = ("F1", "Uz1", "K1")
ROWS = {"F1": "5", "Uz1": "6", "K1": "6"}
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDER = 2
CELL_COUNT = 4


def _as_list(value: Any) -> list[Any]:
    """Return a provenance value as a list without changing its members."""

    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _same_point(value: Any) -> bool:
    """Check that an interval is a directed singleton without a midpoint."""

    try:
        raw = value._mpi_
    except AttributeError as exc:  # pragma: no cover - defensive input guard
        raise ValueError("directed interval required") from exc
    return raw[0] == raw[1]


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    """Make an outward-directed interval from two interval endpoints."""

    return c.make_mpf((lower._mpi_[0], upper._mpi_[1]))


def _nonnegative_bound(c: Any, upper: Any) -> Any:
    """Return the exact directed interval ``[0, upper]``."""

    return _interval(c, c.mpf(0), upper)


def _provenance(term: dict[str, Any]) -> tuple[list[Any], list[Any], list[Any], dict[str, Any]]:
    """Collect source ids, paths, node refs, and the original origin mapping."""

    raw_origin = term.get("origin", {})
    if raw_origin is None:
        raw_origin = {}
    if not isinstance(raw_origin, dict):
        raise ValueError("origin must be a mapping")
    origin = deepcopy(raw_origin)

    source_ids = _as_list(term.get("source_ids"))
    if not source_ids:
        source_ids = _as_list(origin.get("source_ids"))
    if not source_ids and term.get("input_id") is not None:
        source_ids = [term["input_id"]]
    if not source_ids and origin.get("input_id") is not None:
        source_ids = [origin["input_id"]]

    source_paths = _as_list(term.get("source_paths"))
    if not source_paths:
        source_paths = _as_list(term.get("paths"))
    if not source_paths and term.get("source_path") is not None:
        source_paths = [term["source_path"]]
    if not source_paths and term.get("path") is not None:
        source_paths = [term["path"]]
    if not source_paths:
        source_paths = _as_list(origin.get("source_paths"))
    if not source_paths and origin.get("source_path") is not None:
        source_paths = [origin["source_path"]]
    if not source_paths and origin.get("path") is not None:
        source_paths = [origin["path"]]

    refs = term.get("source_node_refs")
    if refs is None:
        refs = origin.get("source_node_refs")
    refs_list = _as_list(deepcopy(refs))
    if not refs_list:
        # Keep an explicit node for every available source id/path.  This is
        # intentionally a reference record, not an inferred source identity.
        count = max(len(source_ids), len(source_paths))
        for index in range(count):
            node: dict[str, Any] = {}
            if index < len(source_ids):
                node["source_id"] = deepcopy(source_ids[index])
            if index < len(source_paths):
                node["source_path"] = deepcopy(source_paths[index])
            refs_list.append(node)

    if not source_ids and not source_paths and not refs_list:
        raise ValueError("source ids, source paths, or source node refs required")

    return source_ids, source_paths, refs_list, origin


def _validate_term(c: Any, term: Any, label: str) -> dict[str, Any]:
    """Validate the fixed formal shape without evaluating any source atom."""

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
    _provenance(term)
    return term


def _mapped_local_term(
    source: dict[str, Any],
    *,
    c: Any,
    factor: Any,
    phase_cell: int,
    field: str,
    source_row: str,
    kernel_sign: int,
) -> dict[str, Any]:
    """Map one source term to one independent local cell term."""

    source_ids, source_paths, source_refs, source_origin = _provenance(source)
    out = deepcopy(source)
    powers = list(source["powers"])
    powers[0] += 2
    coefficient = source["coefficient"] * factor
    if kernel_sign < 0:
        coefficient = -coefficient

    local_origin = deepcopy(source_origin)
    local_origin["source_origin"] = deepcopy(source_origin)
    local_origin["source_ids"] = deepcopy(source_ids)
    if source_paths:
        local_origin["source_paths"] = deepcopy(source_paths)
    local_origin["source_node_refs"] = deepcopy(source_refs)
    local_origin["source_operation"] = source_origin.get("operation")
    local_origin["operation"] = "w0_z2_four_cell_local"
    local_origin["phase_cell"] = phase_cell
    local_origin["source_row"] = source_row
    local_origin["field"] = field
    local_origin["kernel_sign"] = kernel_sign
    local_origin["ordinary_Z_order"] = ORDINARY_Z_ORDER
    local_origin["source_z_derivative_is_ordinary"] = True
    local_origin["ordinary_Z_factorial_divided"] = False
    local_origin["source_is_whole_collar_enclosure"] = True
    local_origin.setdefault("uncertainty", source.get("uncertainty"))
    local_origin.setdefault("error_classes", deepcopy(source.get("error_classes", [])))

    out["powers"] = powers
    out["coefficient"] = coefficient
    out["origin"] = local_origin
    out["source_ids"] = deepcopy(source_ids)
    if source_paths:
        out["source_paths"] = deepcopy(source_paths)
        if len(source_paths) == 1:
            out["source_path"] = deepcopy(source_paths[0])
    out["source_node_refs"] = deepcopy(source_refs)
    out["uncertainty"] = deepcopy(source.get("uncertainty", source_origin.get("uncertainty")))
    out["error_classes"] = deepcopy(source.get("error_classes", source_origin.get("error_classes", [])))
    out["local_enclosure"] = {
        "kind": "uniform_source_times_directed_kernel_bound",
        "kernel_bound": factor,
        "formal_hb_power_increment": 2,
        "whole_collar_source": True,
    }
    return out


def _validate_mapping(mapping: Any, keys: Iterable[str], label: str) -> None:
    if not isinstance(mapping, dict) or set(mapping) != set(keys):
        raise ValueError(f"{label} must contain exactly {tuple(keys)}")


def build(c: Any, *, sc: Any, core: Any, g: Any) -> dict[str, Any]:
    """Build the four independent local ordinary ``W0 Z2`` enclosures.

    ``core`` and ``g`` are intentionally field/row mappings, not the shape of
    any prior leaf.  ``core`` is copied as transport, while each local cell is
    made from source ``g['5']`` and ``g['6']`` without summing its terms.
    """

    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
        raise ValueError("an interval context with mpf, make_mpf, and exp is required")
    _validate_mapping(core, FIELDS, "core")
    _validate_mapping(g, ("5", "6"), "g")

    sc = c.mpf(sc)
    if not _same_point(sc) or sc.a <= 0:
        raise ValueError("sc must be a positive directed singleton")

    for field in FIELDS:
        terms = core[field]
        if not isinstance(terms, (list, tuple)):
            raise ValueError(f"core[{field}] must be a term list")
        if not terms:
            raise ValueError(f"core[{field}] Z2 source row must be non-empty")
        for index, term in enumerate(terms):
            _validate_term(c, term, f"core[{field}][{index}]")
    for row in ("5", "6"):
        terms = g[row]
        if not isinstance(terms, (list, tuple)):
            raise ValueError(f"g[{row}] must be a term list")
        if not terms:
            raise ValueError(f"g[{row}] Z2 source row must be non-empty")
        for index, term in enumerate(terms):
            _validate_term(c, term, f"g[{row}][{index}]")

    core_transport = {field: deepcopy(core[field]) for field in FIELDS}
    q_cells: list[dict[str, Any]] = []
    local: dict[int, dict[str, list[dict[str, Any]]]] = {}
    source_node_refs: dict[str, Any] = {
        "core": {
            field: [{"container": "core_transport", "field": field}]
            for field in FIELDS
        },
        "g": {row: [] for row in ("5", "6")},
    }

    quarter = c.mpf(1) / 4
    sixteenth = c.mpf(1) / 16
    for index in range(CELL_COUNT):
        q_lower = c.mpf(index) / 4
        q_upper = c.mpf(index + 1) / 4
        s_lower = (sc / 4) * q_lower
        s_upper = (sc / 4) * q_upper
        # q_hi is deliberately the cell endpoint; source endpoint selection
        # never enters this bound.
        s_q = (sc / 4) * q_upper * c.exp(sc / 16)
        kernel_bound = _nonnegative_bound(c, (s_q * s_q) / 2)
        x_minus_2_without_hb = _nonnegative_bound(c, s_q)
        cell = {
            "index": index,
            "q_lower": q_lower,
            "q_upper": q_upper,
            "q_hi": q_upper,
            "s_lower": s_lower,
            "s_upper": s_upper,
            "S_q": s_q,
            "S_q_formula": "(sc/4)*q_hi*exp(sc/16)",
            "x_minus_2_bound": {
                "formal_hb_multiplier": 1,
                "interval_without_hb": x_minus_2_without_hb,
                "expression": "|X-2| <= hb*S_q",
            },
            "kernel_bound": kernel_bound,
            "kernel_bound_without_formal_hb2": kernel_bound,
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
            "local_source_rows": {field: ROWS[field] for field in FIELDS},
        }
        q_cells.append(cell)

        local[index] = {field: [] for field in FIELDS}
        for field in FIELDS:
            row = ROWS[field]
            sign = -1 if field == "K1" else 1
            for source_index, source in enumerate(g[row]):
                mapped = _mapped_local_term(
                    source,
                    c=c,
                    factor=kernel_bound,
                    phase_cell=index,
                    field=field,
                    source_row=row,
                    kernel_sign=sign,
                )
                mapped["source_node_refs"] = [
                    {
                        "container": "g",
                        "row": row,
                        "term_index": source_index,
                        "source_refs": deepcopy(mapped["source_node_refs"]),
                    }
                ]
                # Keep the original node records in origin/source_origin while
                # exposing the explicit g-row reference needed for tracing.
                mapped["origin"]["source_node_refs"] = deepcopy(mapped["source_node_refs"])
                local[index][field].append(mapped)
                source_node_refs["g"][row].append(
                    {"cell": index, "field": field, "term_index": source_index}
                )

    return {
        "ordinary_Z_order": ORDINARY_Z_ORDER,
        "ordinary_Z_orders": [ORDINARY_Z_ORDER],
        "fields": list(FIELDS),
        "rows": dict(ROWS),
        "source_power_basis": list(ATOMS),
        "core_transport": core_transport,
        "local": local,
        "q_cells": q_cells,
        "source_node_refs": source_node_refs,
        "source_z2_is_ordinary_derivative": True,
        "ordinary_Z_factorial_divided": False,
        "geometry_Z_independent": True,
        "hb_assumption": "0<hb<1/2",
        "hb_materialized": False,
        "formal_hb2_materialized": False,
        "X_formula": "X=2*exp(hb*s/2)",
        "X_minus_2_bound": "|X-2|<=hb*S_q",
        "g_uniform_on_entire_first_collar_required": True,
        "source_z2_whole_collar_enclosure": True,
        "source_terms_are_values": False,
        "source_endpoint_not_selected": True,
        "kernel_bound_includes_formal_hb2": False,
        "physical_or_cartesian_measure_included": False,
        "no_physical_cartesian_measure": True,
        "kernel_inequalities": {
            "F1": "0<=u/2*(1-(u/X)^2)<=X-u",
            "Uz1": "0<=u*log(X/u)<=X-u",
            "K1": "K1_kernel=-F1_kernel",
        },
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
