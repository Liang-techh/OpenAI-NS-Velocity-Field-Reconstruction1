"""Symbolic K0-only W1 local m15 density products.

This leaf joins an already decoded :mod:`original_saved_k0_w1_cells_v1`
result with caller-owned full leading ``U`` and ``Phi`` rows.  It is a
symbolic product adapter, not a range evaluator: saved binary endpoint
records and background interval records are copied into each product without
conversion, endpoint arithmetic, midpoint selection, or integration.

The six formal source atoms are
``(hb, F0_anchor, Pstar_squared, Lambda, F0_local, L_local)``.  The first
four powers on a saved W1 term are retained in
``saved_formal_source_powers``.  ``formal_source_powers`` is the corresponding
six-atom expression after adding only the external physical factors owned by
this leaf: ``Lambda^-1`` for axial products, ``Lambda^-2 F0_local`` for
swirl products, and one ``L_local`` factor on the local logarithmic-
derivative swirl branch.

No old producer is imported or called.  The caller remains responsible for
establishing that the two saved families are same-source objects and for
providing any omitted Kchi, pressure, or all-word terms.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_k0_symbolic_m15_cells_v1"
SAVED_SCHEMA = "original_saved_k0_w1_cells_v1"
ORDERS = ("Z0", "Z1")
FIELDS = ("F1", "Uz1")
CHANNELS = ("axial", "swirl")
SOURCE_POWER_BASIS = (
    "hb",
    "F0_anchor",
    "Pstar_squared",
    "Lambda",
    "F0_local",
    "L_local",
)
SAVED_SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BACKGROUND_FIELDS = ("U_lead_Z0", "U_lead_Z1", "Phi_Z0", "Phi_Z1")


class SymbolicInputError(ValueError):
    """Raised when a saved sector or background row is structurally unsafe."""


def _copy(value: Any) -> Any:
    """Copy caller data without asking it for a numeric representation."""

    return deepcopy(value)


def _is_int(value: Any) -> bool:
    return type(value) is int


def _endpoint_record(value: Any, label: str) -> Mapping[str, Any]:
    """Validate an authenticated binary endpoint record structurally only."""

    if not isinstance(value, Mapping):
        raise SymbolicInputError(label + " must be a binary interval record")
    for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        item = value.get(key)
        if not isinstance(item, list) or len(item) != 4 or any(type(x) is not int for x in item):
            raise SymbolicInputError(label + " must retain four-integer binary endpoints")
    if "directed" in value and value.get("directed") is not True:
        raise SymbolicInputError(label + " must retain its directed marker")
    return value


def _powers(value: Any, label: str) -> list[int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not _is_int(power) for power in value)
    ):
        raise SymbolicInputError(label + " must contain four integer source powers")
    # These are formal integer labels, not interval values.  The source list
    # itself is copied separately so its original representation is retained.
    return [int(power) for power in value]


def _nonempty(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (str, bytes, bytearray, Sequence, Mapping)):
        return len(value) > 0
    return True


_PROVENANCE_KEYS = (
    "source_id",
    "source_ids",
    "source",
    "term_id",
    "term_ids",
    "source_path",
    "source_paths",
    "identity",
    "id",
)


def _source_values(value: Any) -> list[Any]:
    """Read shallow source labels without traversing a source graph."""

    if not isinstance(value, Mapping):
        return []
    result: list[Any] = []
    for key in _PROVENANCE_KEYS:
        if key not in value or not _nonempty(value.get(key)):
            continue
        raw = value.get(key)
        if isinstance(raw, (str, bytes, bytearray)) or not isinstance(raw, Sequence):
            result.append(_copy(raw))
        else:
            result.extend(_copy(list(raw)))
    for key in ("origin", "metadata", "provenance"):
        child = value.get(key)
        if isinstance(child, Mapping):
            result.extend(_source_values(child))
    return result


def _require_background_interval(value: Any, label: str) -> None:
    """Require an interval-shaped background record without decoding it."""

    if not isinstance(value, Mapping):
        raise SymbolicInputError(label + " must be a background interval record")
    if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
        _endpoint_record(value, label)
        return
    for key in ("interval", "bounds", "range", "value", "coefficient"):
        nested = value.get(key)
        if isinstance(nested, Mapping) and (
            "lower_exact_mpf_tuple" in nested or "upper_exact_mpf_tuple" in nested
        ):
            _endpoint_record(nested, label + "." + key)
            return
    raise SymbolicInputError(label + " must retain a binary interval")


def _background_items(value: Any, label: str) -> list[Mapping[str, Any]]:
    """Return background terms while preserving each original record."""

    if isinstance(value, Mapping):
        nested_key = next(
            (key for key in ("terms", "components", "rows", "values") if key in value),
            None,
        )
        if nested_key is None:
            items: list[Any] = [value]
        else:
            nested = value.get(nested_key)
            if not isinstance(nested, Sequence) or isinstance(nested, (str, bytes, bytearray)):
                raise SymbolicInputError(label + "." + nested_key + " must be a sequence")
            items = list(nested)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        items = list(value)
    else:
        raise SymbolicInputError(label + " must be an interval record or sequence")
    if not items:
        raise SymbolicInputError(label + " must not be empty")
    result: list[Mapping[str, Any]] = []
    for index, item in enumerate(items):
        item_label = f"{label}[{index}]"
        _require_background_interval(item, item_label)
        if not _source_values(item):
            raise SymbolicInputError(item_label + " has no source provenance")
        result.append(item)
    return result


def _background_value(background: Mapping[str, Any], name: str) -> Any:
    if name in background:
        return background[name]
    # A compact caller packet may group the same four named rows under U_lead
    # and Phi.  This is input routing only; output remains in canonical names.
    if name.startswith("U_lead_") and isinstance(background.get("U_lead"), Mapping):
        order = name.rsplit("_", 1)[-1]
        grouped = background["U_lead"]
        if order in grouped:
            return grouped[order]
    if name.startswith("Phi_") and isinstance(background.get("Phi"), Mapping):
        order = name.rsplit("_", 1)[-1]
        grouped = background["Phi"]
        if order in grouped:
            return grouped[order]
    raise SymbolicInputError("saved background is missing " + name)


def _normalise_background(value: Any) -> dict[str, list[Mapping[str, Any]]]:
    if not isinstance(value, Mapping):
        raise TypeError("saved_background mapping required")
    return {
        name: _background_items(_background_value(value, name), "saved_background." + name)
        for name in BACKGROUND_FIELDS
    }


def _validate_original_term(term: Any, label: str, target: int, source_limit: int) -> dict[str, Any]:
    if not isinstance(term, Mapping):
        raise SymbolicInputError(label + ".original_term must be a mapping")
    source = term.get("source_cell")
    term_target = term.get("target_cell")
    if not _is_int(source) or source < 0 or source > source_limit:
        raise SymbolicInputError(label + ".original_term source_cell is invalid")
    if not _is_int(term_target) or term_target != target:
        raise SymbolicInputError(label + ".original_term target_cell is invalid")
    for key in ("matrix_sector_identity", "W0_sector_identity"):
        identity = term.get(key)
        if not isinstance(identity, str) or not identity:
            raise SymbolicInputError(label + ".original_term requires " + key)
    _endpoint_record(term.get("coefficient"), label + ".original_term.coefficient")
    if not _source_values(term) and not all(term.get(key) for key in ("matrix_sector_identity", "W0_sector_identity")):
        raise SymbolicInputError(label + ".original_term has no source provenance")
    return dict(term)


def _validate_saved_term(term: Any, label: str, target: int, source_limit: int) -> dict[str, Any]:
    if not isinstance(term, Mapping):
        raise SymbolicInputError(label + " must be a saved W1 term")
    _endpoint_record(term.get("coefficient"), label + ".coefficient")
    original = _validate_original_term(term.get("original_term"), label, target, source_limit)
    saved_powers = _powers(term.get("formal_source_powers"), label + ".formal_source_powers")
    # The binder's selected coefficient and original record must stay aligned.
    _endpoint_record(original.get("coefficient"), label + ".original_term.coefficient")
    if term.get("coefficient") != original.get("coefficient"):
        raise SymbolicInputError(label + " coefficient differs from original_term")
    if saved_powers != [0, original.get("F0_power"), 0, original.get("Lambda_power")]:
        raise SymbolicInputError(label + " powers differ from saved producer convention")
    if term.get("source_cell", original.get("source_cell")) != original.get("source_cell"):
        raise SymbolicInputError(label + " source_cell differs from original_term")
    if term.get("target_cell", original.get("target_cell")) != original.get("target_cell"):
        raise SymbolicInputError(label + " target_cell differs from original_term")
    return {
        "coefficient": term.get("coefficient"),
        "formal_source_powers": saved_powers,
        "original_formal_source_powers": _copy(term.get("formal_source_powers")),
        "original_term": original,
    }


def _saved_cells(saved_binding: Any) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not isinstance(saved_binding, Mapping):
        raise TypeError("saved_binding mapping required")
    if saved_binding.get("schema") != SAVED_SCHEMA:
        raise SymbolicInputError("saved_binding must be the saved K0 W1 binder result")
    if saved_binding.get("formal_source_basis") != list(SAVED_SOURCE_POWER_BASIS):
        raise SymbolicInputError("saved_binding must retain the four source-power basis")
    scope = saved_binding.get("source_kernel_scope")
    if not isinstance(scope, str) or "K0" not in scope:
        raise SymbolicInputError("saved_binding must retain K0 kernel scope")
    for key in (
        "Kchi_minus_K0_included",
        "allword_tail_added",
        "pressure_added",
        "appended_to_cumulative_m15",
        "actual_W1_complete",
    ):
        if saved_binding.get(key) is not False:
            raise SymbolicInputError("saved_binding scope flag must be false: " + key)
    family = saved_binding.get("source_family")
    logs = saved_binding.get("source_logs")
    if not isinstance(family, Mapping) or not isinstance(logs, Mapping):
        raise SymbolicInputError("saved_binding source family and logs are required")

    raw_cells = saved_binding.get("cells")
    if not isinstance(raw_cells, Sequence) or isinstance(raw_cells, (str, bytes, bytearray)) or len(raw_cells) != 8:
        raise SymbolicInputError("saved_binding must contain eight target cells")
    cells: list[dict[str, Any]] = []
    for index, cell in enumerate(raw_cells):
        label = f"saved_binding.cells[{index}]"
        if not isinstance(cell, Mapping) or cell.get("index") != index:
            raise SymbolicInputError("saved target cell order or index differs")
        expected_sources = list(range(index + 1))
        if cell.get("source_indices") != expected_sources or cell.get("triangular_diagonal_source") != index:
            raise SymbolicInputError(label + " source-cell structure differs")
        for coordinate in ("t", "x", "rho"):
            _endpoint_record(cell.get(coordinate), label + "." + coordinate)
        rows = cell.get("ordinary_Z")
        if not isinstance(rows, Mapping):
            raise SymbolicInputError(label + ".ordinary_Z is missing")
        validated_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in ORDERS:
            row = rows.get(order)
            if not isinstance(row, Mapping):
                raise SymbolicInputError(label + ".ordinary_Z." + order + " is missing")
            validated_rows[order] = {}
            for field in FIELDS:
                terms = row.get(field)
                if not isinstance(terms, Sequence) or isinstance(terms, (str, bytes, bytearray)) or not terms:
                    raise SymbolicInputError(label + ".ordinary_Z." + order + "." + field + " is missing")
                validated_rows[order][field] = [
                    _validate_saved_term(term, f"{label}.{order}.{field}[{term_index}]", index, index)
                    for term_index, term in enumerate(terms)
                ]
        cells.append(
            {
                "index": index,
                "source_indices": _copy(cell["source_indices"]),
                "triangular_diagonal_source": index,
                "axis_endpoint_continuation": _copy(cell.get("axis_endpoint_continuation")),
                "t": _copy(cell["t"]),
                "x": _copy(cell["x"]),
                "rho": _copy(cell["rho"]),
                "ordinary_Z": validated_rows,
            }
        )
    return {
        "source_family": _copy(family),
        "source_logs": _copy(logs),
        "source_kernel_scope": scope,
    }, cells


def _formal_product_powers(saved: Sequence[int], channel: str, *, local_log: bool = False) -> list[int]:
    powers = [int(power) for power in saved] + [0, 0]
    if channel == "axial":
        powers[3] -= 1
    elif channel == "swirl":
        powers[3] -= 2
        powers[4] += 1
    else:
        raise SymbolicInputError("unknown symbolic channel " + channel)
    if local_log:
        powers[5] += 1
    return powers


def _term_source_ids(original: Mapping[str, Any], background: Mapping[str, Any]) -> list[Any]:
    values = _source_values(original)
    values.extend(_source_values(background))
    # Matrix and W0 identities are required source identities even when the
    # historical record did not duplicate them under source_id.
    for key in ("matrix_sector_identity", "W0_sector_identity"):
        if key in original and original[key] not in values:
            values.append(_copy(original[key]))
    return values


def _product_term(
    *,
    channel: str,
    ordinary_order: str,
    piece: str,
    expression: str,
    signed_factor: int,
    target_cell: Mapping[str, Any],
    saved_field: str,
    saved_order: str,
    saved_term: Mapping[str, Any],
    background_name: str,
    background_factor: Mapping[str, Any],
    rho_factor: Any = None,
    local_log: bool = False,
) -> dict[str, Any]:
    original = saved_term["original_term"]
    source_powers = saved_term["formal_source_powers"]
    output = {
        "channel": channel,
        "ordinary_Z": ordinary_order,
        "ordinary_Z_order": ORDERS.index(ordinary_order),
        "piece": piece,
        "expression": expression,
        "signed_factor": signed_factor,
        "sign": 1 if signed_factor > 0 else -1,
        "coefficient": _copy(saved_term["coefficient"]),
        "coefficient_is_unmodified_saved_binary_record": True,
        "original_term": _copy(original),
        "saved_field": saved_field,
        "saved_Z": saved_order,
        "saved_formal_source_powers": _copy(saved_term["original_formal_source_powers"]),
        "formal_source_powers": _formal_product_powers(source_powers, channel, local_log=local_log),
        "background_factor_name": background_name,
        "background_factor": _copy(background_factor),
        "target_cell": target_cell["index"],
        "source_cell": original["source_cell"],
        "matrix_sector_identity": _copy(original["matrix_sector_identity"]),
        "W0_sector_identity": _copy(original["W0_sector_identity"]),
        "source_indices": _copy(target_cell["source_indices"]),
        "source_cell_provenance": {
            "target_cell": target_cell["index"],
            "source_cell": original["source_cell"],
            "source_indices": _copy(target_cell["source_indices"]),
            "triangular_diagonal_source": target_cell["triangular_diagonal_source"],
            "matrix_sector_identity": _copy(original["matrix_sector_identity"]),
            "W0_sector_identity": _copy(original["W0_sector_identity"]),
        },
        "source_ids": _term_source_ids(original, background_factor),
        "source_identity_status": "source identity unverified",
        "product_factors": [],
        "unknown": False,
    }
    output["product_factors"] = [background_name]
    if rho_factor is not None:
        output["rho_factor"] = _copy(rho_factor)
        output["product_factors"].insert(0, "rho")
    if local_log:
        output["local_factor"] = {"symbol": "L_local", "formal_power": 1}
        output["product_factors"].insert(1 if rho_factor is not None else 0, "L_local")
    output["product_factors"].append(saved_field + ("_0" if saved_order == "Z0" else "_Z"))
    return output


def _unknown_terms() -> list[dict[str, str]]:
    return [
        {"name": "N1p", "status": "unknown", "reason": "pressure channel is outside this symbolic velocity leaf"},
        {"name": "Kchi", "status": "unknown", "reason": "unsupported Kchi-minus-K0 correction"},
        {"name": "allword_tail", "status": "unknown", "reason": "no saved all-word tail is supplied"},
        {"name": "K1", "status": "unknown", "reason": "omitted W1 component"},
        {"name": "P1", "status": "unknown", "reason": "omitted W1 component"},
        {"name": "d_x_F1", "status": "unknown", "reason": "omitted radial derivative component"},
        {"name": "d_x_Uz1", "status": "unknown", "reason": "omitted radial derivative component"},
    ]


def build(saved_binding: Any, saved_background: Any) -> dict[str, Any]:
    """Build unintegrated symbolic axial/swirl products for eight K0 cells.

    ``saved_binding`` must be the already decoded result of the saved K0 W1
    binder.  ``saved_background`` supplies ``U_lead_Z0``, ``U_lead_Z1``,
    ``Phi_Z0`` and ``Phi_Z1`` interval records (or finite sequences of such
    records).  Neither argument is numerically evaluated.
    """

    saved_meta, cells = _saved_cells(saved_binding)
    background = _normalise_background(saved_background)

    output_cells: list[dict[str, Any]] = []
    counts = {channel: {order: 0 for order in ORDERS} for channel in CHANNELS}
    for cell in cells:
        index = cell["index"]
        terms = {channel: {order: [] for order in ORDERS} for channel in CHANNELS}

        for background_factor in background["U_lead_Z0"]:
            for saved_term in cell["ordinary_Z"]["Z0"]["Uz1"]:
                terms["axial"]["Z0"].append(
                    _product_term(
                        channel="axial",
                        ordinary_order="Z0",
                        piece="U0*Uz1_0",
                        expression="2*U0*Uz1_0",
                        signed_factor=2,
                        target_cell=cell,
                        saved_field="Uz1",
                        saved_order="Z0",
                        saved_term=saved_term,
                        background_name="U_lead_Z0",
                        background_factor=background_factor,
                    )
                )
        for background_factor in background["U_lead_Z1"]:
            for saved_term in cell["ordinary_Z"]["Z0"]["Uz1"]:
                terms["axial"]["Z1"].append(
                    _product_term(
                        channel="axial",
                        ordinary_order="Z1",
                        piece="U1*Uz1_0",
                        expression="2*U1*Uz1_0",
                        signed_factor=2,
                        target_cell=cell,
                        saved_field="Uz1",
                        saved_order="Z0",
                        saved_term=saved_term,
                        background_name="U_lead_Z1",
                        background_factor=background_factor,
                    )
                )
        for background_factor in background["U_lead_Z0"]:
            for saved_term in cell["ordinary_Z"]["Z1"]["Uz1"]:
                terms["axial"]["Z1"].append(
                    _product_term(
                        channel="axial",
                        ordinary_order="Z1",
                        piece="U0*Uz1_Z",
                        expression="2*U0*Uz1_Z",
                        signed_factor=2,
                        target_cell=cell,
                        saved_field="Uz1",
                        saved_order="Z1",
                        saved_term=saved_term,
                        background_name="U_lead_Z0",
                        background_factor=background_factor,
                    )
                )

        for background_factor in background["Phi_Z0"]:
            for saved_term in cell["ordinary_Z"]["Z0"]["F1"]:
                terms["swirl"]["Z0"].append(
                    _product_term(
                        channel="swirl",
                        ordinary_order="Z0",
                        piece="rho*Phi0*F1_0",
                        expression="-2*rho*Phi0*F1_0",
                        signed_factor=-2,
                        target_cell=cell,
                        saved_field="F1",
                        saved_order="Z0",
                        saved_term=saved_term,
                        background_name="Phi_Z0",
                        background_factor=background_factor,
                        rho_factor=cell["rho"],
                    )
                )
        for background_factor in background["Phi_Z1"]:
            for saved_term in cell["ordinary_Z"]["Z0"]["F1"]:
                terms["swirl"]["Z1"].append(
                    _product_term(
                        channel="swirl",
                        ordinary_order="Z1",
                        piece="rho*Phi1*F1_0",
                        expression="-2*rho*Phi1*F1_0",
                        signed_factor=-2,
                        target_cell=cell,
                        saved_field="F1",
                        saved_order="Z0",
                        saved_term=saved_term,
                        background_name="Phi_Z1",
                        background_factor=background_factor,
                        rho_factor=cell["rho"],
                    )
                )
        for background_factor in background["Phi_Z0"]:
            for saved_term in cell["ordinary_Z"]["Z0"]["F1"]:
                terms["swirl"]["Z1"].append(
                    _product_term(
                        channel="swirl",
                        ordinary_order="Z1",
                        piece="L_local*Phi0*F1_0",
                        expression="-2*rho*L_local*Phi0*F1_0",
                        signed_factor=-2,
                        target_cell=cell,
                        saved_field="F1",
                        saved_order="Z0",
                        saved_term=saved_term,
                        background_name="Phi_Z0",
                        background_factor=background_factor,
                        rho_factor=cell["rho"],
                        local_log=True,
                    )
                )
        for background_factor in background["Phi_Z0"]:
            for saved_term in cell["ordinary_Z"]["Z1"]["F1"]:
                terms["swirl"]["Z1"].append(
                    _product_term(
                        channel="swirl",
                        ordinary_order="Z1",
                        piece="Phi0*F1_Z",
                        expression="-2*rho*Phi0*F1_Z",
                        signed_factor=-2,
                        target_cell=cell,
                        saved_field="F1",
                        saved_order="Z1",
                        saved_term=saved_term,
                        background_name="Phi_Z0",
                        background_factor=background_factor,
                        rho_factor=cell["rho"],
                    )
                )

        for channel in CHANNELS:
            for order in ORDERS:
                counts[channel][order] += len(terms[channel][order])
        output_cells.append(
            {
                "index": index,
                "target_cell": index,
                "source_indices": _copy(cell["source_indices"]),
                "triangular_diagonal_source": cell["triangular_diagonal_source"],
                "axis_endpoint_continuation": _copy(cell["axis_endpoint_continuation"]),
                "t": _copy(cell["t"]),
                "x": _copy(cell["x"]),
                "rho": _copy(cell["rho"]),
                "terms": terms,
                "term_counts": {
                    channel: {order: len(terms[channel][order]) for order in ORDERS}
                    for channel in CHANNELS
                },
            }
        )

    return {
        "schema": SCHEMA,
        "status": "symbolic_unintegrated",
        "source_family": saved_meta["source_family"],
        "source_logs": saved_meta["source_logs"],
        "source_kernel_scope": saved_meta["source_kernel_scope"],
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDERS),
        "channels": list(CHANNELS),
        "formulas": {
            "axial_Z0": "2*U0*Uz1_0",
            "axial_Z1": "2*(U1*Uz1_0 + U0*Uz1_Z)",
            "swirl_Z0": "-2*rho*Phi0*F1_0",
            "swirl_Z1": "-2*rho*(Phi1*F1_0 + L_local*Phi0*F1_0 + Phi0*F1_Z)",
        },
        "cells": output_cells,
        "term_counts": counts,
        "unknown_terms": _unknown_terms(),
        "unknown_terms_are_not_zero": True,
        "source_identity_status": "source identity unverified",
        "source_identity_verified": False,
        "epsilon_identity": {
            "expression": "epsilon = Lambda^-1",
            "status": "caller_declared_premise",
            "verified_here": False,
        },
        "geometry_and_cutoff_premises": {
            "rho_Z": "zero at fixed rho", "epsilon_Z": "zero",
            "eta": "one on the core under the pinned cutoff definition",
            "same_live_cutoff_owner_verified_here": False,
        },
        "F0_scale_identity": {
            "F0_anchor": "independent formal atom",
            "F0_local": "independent formal atom",
            "identified": False,
        },
        "background_U_is_full_lead_not_Delta": True,
        "ordinary_Z1_saved_rows_reused": True,
        "saved_Z1_coefficients_or_powers_differentiated": False,
        "dyadic_endpoints_embedded_unchanged": True,
        "fraction_conversion_performed": False,
        "interval_arithmetic_performed": False,
        "integration_performed": False,
        "numeric_operator_called": False,
        "source_values_computed": False,
        "point_value_selected": False,
        "N1p_included": False,
        "Kchi_included": False,
        "allword_tail_included": False,
        "appended_to_212_ledger": False,
        "ledger_consulted": False,
        "complete_W1": False,
        "complete_m15_computed": False,
        "whole_grid_recertified": False,
        "independent_value": False,
        "scope": "K0-only W1 local m15 density symbolic products; no integration",
    }


__all__ = [
    "SCHEMA",
    "ORDERS",
    "FIELDS",
    "CHANNELS",
    "SOURCE_POWER_BASIS",
    "SymbolicInputError",
    "build",
]
