"""Symbolic direct swirl products for the finite W0 Kchi first word.

This leaf joins an already-bound finite ``(-2*chi) T K0`` W0 packet to the
background ``Phi_Z0``/``Phi_Z1`` records retained by the earlier K0 symbolic
output.  It is intentionally a structural adapter:

* source coefficient interval records are copied, never multiplied;
* the four source powers are retained and extended only by formal labels;
* the one Z0 and three Z1 product branches remain unmerged;
* the ``j>=2`` packet remains an error-only branch; and
* ``d_x_F1`` and all Z2 rows remain available for later feedback.

No producer, archive, interval operator, numerical evaluator, quadrature, or
ancestor helper is imported or called here.  Equality of source-family/log
records is only a conditional compatibility check; it is not an owner proof.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_kchi_first_swirl_symbolic_v1"
KCHI_SCHEMA = "original_kchi_first_saved_cells_v2"
K0_SYMBOLIC_SCHEMA = "original_k0_symbolic_m15_cells_v1"
ORDERS = ("Z0", "Z1")
ALL_ORDERS = ("Z0", "Z1", "Z2")
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_POWER_BASIS = (
    "hb",
    "F0_anchor",
    "Pstar_squared",
    "Lambda",
    "F0_local",
    "L_local",
)
BACKGROUND_NAMES = ("Phi_Z0", "Phi_Z1")
TAIL_ERROR_CLASS = "Kchi_j_ge_2_kernel_tail"


class SwirlSymbolicInputError(ValueError):
    """Raised when the two saved symbolic packets are structurally unsafe."""


def _copy(value: Any) -> Any:
    """Copy without requesting a numeric representation."""

    return deepcopy(value)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SwirlSymbolicInputError(label + " must be a mapping")
    return value


def _sequence(value: Any, label: str, *, nonempty: bool = False) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise SwirlSymbolicInputError(label + " must be a sequence")
    if nonempty and not value:
        raise SwirlSymbolicInputError(label + " must not be empty")
    return value


def _int_powers(value: Any, label: str, *, size: int = 4) -> list[int]:
    raw = _sequence(value, label)
    if len(raw) != size or any(type(item) is not int for item in raw):
        raise SwirlSymbolicInputError(
            label + " must contain " + str(size) + " integer formal powers"
        )
    return [int(item) for item in raw]


def _nonempty_ids(value: Any, label: str) -> list[Any]:
    raw = _sequence(value, label, nonempty=True)
    result = _copy(list(raw))
    if any(item is None or item == "" for item in result):
        raise SwirlSymbolicInputError(label + " must contain nonempty source IDs")
    return result


def _require_endpoint_shape(value: Any, label: str) -> None:
    """Check the saved interval shape without decoding its endpoint values."""

    record = _mapping(value, label)
    for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        item = record.get(key)
        if (
            not isinstance(item, list)
            or len(item) != 4
            or any(type(part) is not int for part in item)
        ):
            raise SwirlSymbolicInputError(label + " must retain binary endpoints")
    if "directed" in record and record.get("directed") is not True:
        raise SwirlSymbolicInputError(label + " must retain directed=True")


def _require_family_and_logs(kchi: Mapping[str, Any], old: Mapping[str, Any]) -> tuple[Any, Any]:
    family = kchi.get("source_family")
    old_family = old.get("source_family")
    logs = kchi.get("source_logs")
    old_logs = old.get("source_logs")
    if not isinstance(family, Mapping) or not isinstance(old_family, Mapping):
        raise SwirlSymbolicInputError("both packets require source_family mappings")
    if not isinstance(logs, Mapping) or not isinstance(old_logs, Mapping):
        raise SwirlSymbolicInputError("both packets require source_logs mappings")
    if family != old_family:
        raise SwirlSymbolicInputError("source_family labels do not match")
    if logs != old_logs:
        raise SwirlSymbolicInputError("source_logs labels do not match")
    return _copy(dict(family)), _copy(dict(logs))


def _validate_packets(kchi: Any, old: Any) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    kchi_map = _mapping(kchi, "kchi_binding")
    old_map = _mapping(old, "k0_symbolic_output")
    if kchi_map.get("schema") != KCHI_SCHEMA:
        raise SwirlSymbolicInputError("kchi_binding must be the saved Kchi binder output")
    if old_map.get("schema") != K0_SYMBOLIC_SCHEMA:
        raise SwirlSymbolicInputError("k0_symbolic_output must be the saved K0 symbolic output")
    if kchi_map.get("finite_j1_correction") is not True:
        raise SwirlSymbolicInputError("finite j=1 Kchi correction is required")
    if kchi_map.get("finite_correction_is_signed") is not True:
        raise SwirlSymbolicInputError("finite Kchi correction must remain signed")
    if kchi_map.get("tail_is_error_only") is not True:
        raise SwirlSymbolicInputError("Kchi tail must remain error-only")
    if kchi_map.get("formal_source_basis") not in (None, list(SOURCE_POWER_BASIS)):
        raise SwirlSymbolicInputError("Kchi four-source basis differs")
    old_basis = old_map.get("source_power_basis")
    if old_basis is not None and list(old_basis) != list(OUTPUT_POWER_BASIS):
        raise SwirlSymbolicInputError("K0 symbolic six-source basis differs")
    _require_family_and_logs(kchi_map, old_map)
    return kchi_map, old_map


def _cell_map(packet: Mapping[str, Any], label: str) -> dict[int, Mapping[str, Any]]:
    raw = _sequence(packet.get("cells"), label, nonempty=True)
    result: dict[int, Mapping[str, Any]] = {}
    for position, cell in enumerate(raw):
        row = _mapping(cell, f"{label}[{position}]")
        index = row.get("index", row.get("target_cell"))
        if type(index) is not int or index < 0 or index != position:
            raise SwirlSymbolicInputError(label + " has noncanonical cell order")
        if index in result:
            raise SwirlSymbolicInputError(label + " contains duplicate cell index")
        result[index] = row
    return result


def _rho_record(kchi_cell: Mapping[str, Any], old_cell: Mapping[str, Any], index: int) -> tuple[Any, str]:
    kchi_rho = kchi_cell.get("rho")
    old_rho = old_cell.get("rho")
    if kchi_rho is not None:
        _require_endpoint_shape(kchi_rho, f"kchi_binding.cells[{index}].rho")
    if old_rho is not None:
        _require_endpoint_shape(old_rho, f"k0_symbolic_output.cells[{index}].rho")
    if kchi_rho is not None and old_rho is not None and kchi_rho != old_rho:
        raise SwirlSymbolicInputError("same-index rho records differ")
    if kchi_rho is not None:
        return _copy(kchi_rho), "kchi_binding.cells[index].rho"
    if old_rho is not None:
        return _copy(old_rho), "k0_symbolic_output.cells[index].rho"
    raise SwirlSymbolicInputError("a same-index rho record is required")


def _background_from_record(value: Any, label: str) -> dict[str, Any]:
    record = _mapping(value, label)
    # The old symbolic term normally stores the source interval as
    # ``{"interval": <directed endpoint record>, ...provenance}``.  Accept
    # the direct form as well, but always return the complete original wrapper
    # unchanged; the inner endpoint record is checked without decoding it.
    if "lower_exact_mpf_tuple" in record or "upper_exact_mpf_tuple" in record:
        _require_endpoint_shape(record, label)
    else:
        nested = None
        nested_key = None
        for key in ("interval", "bounds", "range", "value", "coefficient"):
            candidate = record.get(key)
            if isinstance(candidate, Mapping) and (
                "lower_exact_mpf_tuple" in candidate
                or "upper_exact_mpf_tuple" in candidate
            ):
                nested = candidate
                nested_key = key
                break
        if nested is None or nested_key is None:
            raise SwirlSymbolicInputError(label + " must retain an interval record")
        _require_endpoint_shape(nested, label + "." + nested_key)
    return _copy(dict(record))


def _background_from_old_cell(old_cell: Mapping[str, Any], name: str, index: int) -> list[dict[str, Any]]:
    terms_root = _mapping(old_cell.get("terms"), f"k0_symbolic_output.cells[{index}].terms")
    found: list[dict[str, Any]] = []
    seen_sources: list[Any] = []
    for order in ORDERS:
        channel = _mapping(terms_root.get("swirl"), f"k0_symbolic_output.cells[{index}].terms.swirl")
        rows = _sequence(channel.get(order), f"k0_symbolic_output.cells[{index}].terms.swirl.{order}")
        for position, term in enumerate(rows):
            row = _mapping(term, f"k0_symbolic_output.cells[{index}].terms.swirl.{order}[{position}]")
            if row.get("background_factor_name") != name:
                continue
            background = row.get("background_factor")
            source_ids = row.get("source_ids")
            original_background = _background_from_record(
                background, f"old {name}[{index}][{position}]"
            )
            # Key on the background record itself, not on the old W1 term's
            # union of source IDs: each saved W1 term can carry a different
            # source ID while pointing to the same Phi source.
            marker = _background_source_ids(original_background)
            if not marker:
                marker = [original_background]
            if marker in seen_sources:
                if found[seen_sources.index(marker)] != original_background:
                    raise SwirlSymbolicInputError(
                        "one background source ID has inconsistent saved records"
                    )
                continue
            seen_sources.append(marker)
            found.append(original_background)
    if not found:
        raise SwirlSymbolicInputError(
            f"old symbolic cell {index} has no retained {name} background record"
        )
    return found


def _normalise_explicit_background(
    background: Any,
    name: str,
    index: int,
) -> list[dict[str, Any]]:
    root = _mapping(background, "background")
    value = root.get(name)
    if value is None:
        raise SwirlSymbolicInputError("explicit background is missing " + name)
    # Accept either one global record, a sequence of records, or a cell-indexed
    # mapping.  This is routing only; records are copied unchanged.
    if isinstance(value, Mapping) and index in value:
        value = value[index]
    elif isinstance(value, Mapping) and str(index) in value:
        value = value[str(index)]
    if isinstance(value, Mapping):
        return [_background_from_record(value, f"background.{name}")]
    rows = _sequence(value, f"background.{name}", nonempty=True)
    return [
        _background_from_record(item, f"background.{name}[{position}]")
        for position, item in enumerate(rows)
    ]


def _background_rows(
    old_cell: Mapping[str, Any],
    name: str,
    index: int,
    background: Any,
) -> list[dict[str, Any]]:
    if background is not None:
        return _normalise_explicit_background(background, name, index)
    return _background_from_old_cell(old_cell, name, index)


def _finite_source_rows(kchi_cell: Mapping[str, Any], index: int) -> dict[str, list[dict[str, Any]]]:
    source_cells = _sequence(
        kchi_cell.get("source_cells"),
        f"kchi_binding.cells[{index}].source_cells",
        nonempty=True,
    )
    output = {order: [] for order in ORDERS}
    for source_position, source in enumerate(source_cells):
        source_row = _mapping(
            source,
            f"kchi_binding.cells[{index}].source_cells[{source_position}]",
        )
        source_index = source_row.get("source_index")
        target_index = source_row.get("target_index", source_row.get("target_cell"))
        if type(source_index) is not int or type(target_index) is not int:
            raise SwirlSymbolicInputError("Kchi source/target indices are required")
        if source_index < 0 or source_index > index or target_index != index:
            raise SwirlSymbolicInputError("Kchi source/target indices are inconsistent")
        ordinary = _mapping(
            source_row.get("ordinary_Z"),
            f"kchi_binding.cells[{index}].source_cells[{source_position}].ordinary_Z",
        )
        for order in ORDERS:
            rows = _mapping(
                ordinary.get(order),
                f"kchi_binding.cells[{index}].source_cells[{source_position}].ordinary_Z.{order}",
            )
            terms = _sequence(
                rows.get("F1"),
                f"kchi_binding.cells[{index}].source_cells[{source_position}].ordinary_Z.{order}.F1",
                nonempty=True,
            )
            for position, term in enumerate(terms):
                row = _mapping(
                    term,
                    f"kchi_binding.cells[{index}].source_cells[{source_position}].{order}.F1[{position}]",
                )
                if row.get("error_only") is True or row.get("symmetric_error") is True:
                    raise SwirlSymbolicInputError("error-only term entered finite swirl branch")
                if row.get("term_role") not in (None, "finite_j1_signed"):
                    raise SwirlSymbolicInputError("finite F1 term has an unexpected role")
                coefficient = row.get("coefficient")
                _require_endpoint_shape(coefficient, "finite F1 coefficient")
                powers_value = row.get("formal_source_powers", row.get("powers"))
                powers = _int_powers(powers_value, "finite F1 formal_source_powers")
                source_ids = _nonempty_ids(row.get("source_ids"), "finite F1 source_ids")
                if row.get("source_cell") != source_index or row.get("target_cell") != index:
                    raise SwirlSymbolicInputError("finite F1 source/target ownership differs")
                selected = _copy(dict(row))
                selected["_source_index"] = source_index
                selected["_source_term_position"] = position
                selected["_source_powers"] = powers
                selected["_source_ids"] = source_ids
                output[order].append(selected)
    for order in ORDERS:
        if not output[order]:
            raise SwirlSymbolicInputError("finite F1 rows are incomplete")
    return output


def _formal_output_powers(source_powers: Sequence[int], *, local_log: bool) -> list[int]:
    # These are formal labels, not values or interval arithmetic.
    powers = [int(item) for item in source_powers] + [1, 1 if local_log else 0]
    powers[3] -= 2
    return powers


def _background_source_ids(background: Mapping[str, Any]) -> list[Any]:
    for key in ("source_ids", "source_id", "term_ids", "term_id", "source_path", "source_paths"):
        value = background.get(key)
        if value is None or value == "":
            continue
        if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
            return [_copy(value)]
        return _copy(list(value))
    return []


def _archive_reference(term: Mapping[str, Any], background: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "finite_term_original": _copy(term.get("original_term", term)),
        "finite_term_origin": _copy(term.get("origin")),
        "finite_source_ids": _copy(term.get("source_ids")),
        "background_source_ids": _background_source_ids(background),
        "background_record": _copy(background),
    }


def _swirl_term(
    source_term: Mapping[str, Any],
    background: Mapping[str, Any],
    rho: Any,
    *,
    target_cell: int,
    branch: str,
    ordinary_order: str,
    piece: str,
    background_name: str,
    saved_order: str,
    local_log: bool,
    rho_source: str,
) -> dict[str, Any]:
    source_powers = source_term["_source_powers"]
    source_ids = _copy(source_term["_source_ids"])
    output = {
        "channel": "swirl",
        "branch": branch,
        "ordinary_Z": ordinary_order,
        "ordinary_Z_order": 0 if ordinary_order == "Z0" else 1,
        "piece": piece,
        "expression": piece,
        "signed_factor": -2,
        "sign": -1,
        "coefficient": _copy(source_term["coefficient"]),
        "coefficient_is_unmodified_saved_binary_record": True,
        "coefficient_was_not_multiplied": True,
        "source_formal_source_powers": _copy(source_powers),
        "saved_formal_source_powers": _copy(source_powers),
        "formal_source_powers": _formal_output_powers(source_powers, local_log=local_log),
        "formal_source_basis": list(OUTPUT_POWER_BASIS),
        "physical_scale_factors": ["epsilon^2=Lambda^-2", "F0_local", "rho", "-2"],
        "F0_anchor_kept_independent_from_F0_local": True,
        "L_local": {"symbol": "L_local", "formal_power": 1} if local_log else None,
        "background_factor_name": background_name,
        "background_factor": _copy(background),
        "rho_factor": _copy(rho),
        "rho_source": rho_source,
        "target_cell": target_cell,
        "source_cell": source_term.get("source_cell"),
        "source_cell_id": _copy(source_term.get("source_cell_id")),
        "source_term_id": _copy(source_term.get("source_term_id")),
        "source_index": source_term.get("source_index", source_term.get("_source_index")),
        "source_term_index": source_term.get("source_term_index", source_term.get("_source_term_position")),
        "source_ids": source_ids,
        "original_source_ids": _copy(source_ids),
        "source_ids_preserved": True,
        "archive_reference": _archive_reference(source_term, background),
        "original_kchi_term": _copy(source_term),
        "finite_j1_signed": True,
        "error_only": False,
        "symmetric_error": False,
        "signed_field": True,
        "unknown": False,
    }
    # Internal routing keys must not leak as a second source representation.
    output.pop("original_kchi_term", None)
    output["original_kchi_term"] = _copy(
        {key: value for key, value in source_term.items() if not key.startswith("_")}
    )
    return output


def _tail_error(kchi: Mapping[str, Any]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    raw = _mapping(kchi.get("symmetric_kernel_tail"), "kchi_binding.symmetric_kernel_tail")
    output: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for field in ("F1", "d_x_F1"):
        field_rows = _mapping(raw.get(field), "symmetric_kernel_tail." + field)
        output[field] = {}
        for order in ALL_ORDERS:
            rows = field_rows.get(order, field_rows.get(str(int(order.removeprefix("Z")))))
            terms = _sequence(rows, "symmetric_kernel_tail." + field + "." + order, nonempty=True)
            selected: list[dict[str, Any]] = []
            for position, term in enumerate(terms):
                row = _mapping(term, "symmetric_kernel_tail term")
                if row.get("error_only") is False or row.get("signed_field") is True:
                    raise SwirlSymbolicInputError("j>=2 tail cannot claim signed use")
                if TAIL_ERROR_CLASS not in _sequence(
                    row.get("error_classes"), "tail error_classes", nonempty=True
                ):
                    raise SwirlSymbolicInputError("j>=2 tail error class is missing")
                coefficient = row.get("coefficient")
                _require_endpoint_shape(coefficient, "tail coefficient")
                selected_row = _copy(dict(row))
                selected_row.update(
                    {
                        "field": field,
                        "ordinary_Z": order,
                        "term_role": "j_ge_2_symmetric_error",
                        "error_only": True,
                        "symmetric_error": True,
                        "signed_field": False,
                        "coefficient_was_not_multiplied": True,
                        "source_term_position": position,
                    }
                )
                selected.append(selected_row)
            output[field][order] = selected
    return output


def _feedback_reserved(kchi_cells: Mapping[int, Mapping[str, Any]]) -> dict[str, Any]:
    return {
        "cells": _copy([kchi_cells[index] for index in sorted(kchi_cells)]),
        "fields": ["F1", "d_x_F1"],
        "ordinary_Z_orders": list(ALL_ORDERS),
        "d_x_F1_retained": True,
        "Z2_retained": True,
        "purpose": "reserved unchanged Kchi W0 input for later C0/B1 feedback",
    }


def build(
    kchi_binding: Any,
    k0_symbolic_output: Any,
    *,
    background: Any | None = None,
) -> dict[str, Any]:
    """Build unintegrated W0_chi swirl products through ordinary Z1.

    ``background`` is optional.  When omitted, the original ``Phi_Z0`` and
    ``Phi_Z1`` records are recovered from the old symbolic term records.  An
    explicit mapping is useful when a runner has already extracted those
    records; it is still copied symbolically and never numerically evaluated.
    """

    kchi, old = _validate_packets(kchi_binding, k0_symbolic_output)
    family, logs = _require_family_and_logs(kchi, old)
    kchi_cells = _cell_map(kchi, "kchi_binding.cells")
    old_cells = _cell_map(old, "k0_symbolic_output.cells")
    if set(kchi_cells) != set(old_cells):
        raise SwirlSymbolicInputError("Kchi and K0 symbolic cells do not align")

    output_cells: list[dict[str, Any]] = []
    finite_by_order = {order: [] for order in ORDERS}
    finite_counts = {order: 0 for order in ORDERS}
    for index in sorted(kchi_cells):
        kchi_cell = kchi_cells[index]
        old_cell = old_cells[index]
        rho, rho_source = _rho_record(kchi_cell, old_cell, index)
        source_rows = _finite_source_rows(kchi_cell, index)
        phi0 = _background_rows(old_cell, "Phi_Z0", index, background)
        phi1 = _background_rows(old_cell, "Phi_Z1", index, background)

        terms = {"Z0": [], "Z1": []}
        # The one Z0 branch: Phi0 * F1_0.
        for source_term in source_rows["Z0"]:
            for background_record in phi0:
                terms["Z0"].append(
                    _swirl_term(
                        source_term,
                        background_record,
                        rho,
                        target_cell=index,
                        branch="Z0_Phi0_F1_0",
                        ordinary_order="Z0",
                        piece="-2*rho*F0_local*Phi0*W0_chi_F1_0",
                        background_name="Phi_Z0",
                        saved_order="Z0",
                        local_log=False,
                        rho_source=rho_source,
                    )
                )

        # The three Z1 branches: Phi1*F1_0, L*Phi0*F1_0, Phi0*F1_Z.
        for source_term in source_rows["Z0"]:
            for background_record in phi1:
                terms["Z1"].append(
                    _swirl_term(
                        source_term,
                        background_record,
                        rho,
                        target_cell=index,
                        branch="Z1_Phi1_F1_0",
                        ordinary_order="Z1",
                        piece="-2*rho*F0_local*Phi1*W0_chi_F1_0",
                        background_name="Phi_Z1",
                        saved_order="Z0",
                        local_log=False,
                        rho_source=rho_source,
                    )
                )
        for source_term in source_rows["Z0"]:
            for background_record in phi0:
                terms["Z1"].append(
                    _swirl_term(
                        source_term,
                        background_record,
                        rho,
                        target_cell=index,
                        branch="Z1_L_local_Phi0_F1_0",
                        ordinary_order="Z1",
                        piece="-2*rho*F0_local*L_local*Phi0*W0_chi_F1_0",
                        background_name="Phi_Z0",
                        saved_order="Z0",
                        local_log=True,
                        rho_source=rho_source,
                    )
                )
        for source_term in source_rows["Z1"]:
            for background_record in phi0:
                terms["Z1"].append(
                    _swirl_term(
                        source_term,
                        background_record,
                        rho,
                        target_cell=index,
                        branch="Z1_Phi0_F1_Z",
                        ordinary_order="Z1",
                        piece="-2*rho*F0_local*Phi0*W0_chi_F1_Z",
                        background_name="Phi_Z0",
                        saved_order="Z1",
                        local_log=False,
                        rho_source=rho_source,
                    )
                )

        for order in ORDERS:
            finite_by_order[order].extend(terms[order])
            finite_counts[order] += len(terms[order])
        output_cells.append(
            {
                "index": index,
                "target_cell": index,
                "source_indices": _copy(kchi_cell.get("source_indices")),
                "triangular_diagonal_source": _copy(kchi_cell.get("triangular_diagonal_source")),
                "x": _copy(kchi_cell.get("x")),
                "target_x": _copy(kchi_cell.get("target_x")),
                "rho": rho,
                "rho_source": rho_source,
                "terms": terms,
                "term_counts": {order: len(terms[order]) for order in ORDERS},
                "d_x_F1_and_Z2_not_consumed": True,
            }
        )

    tail = _tail_error(kchi)
    return {
        "schema": SCHEMA,
        "status": "symbolic_unintegrated",
        "channel": "swirl",
        "branch": "W0_chi",
        "source_family": family,
        "source_logs": logs,
        "source_identity_status": "conditional source_family/source_logs match only; owner unverified",
        "source_identity_verified": False,
        "owner_verified": False,
        "source_power_basis": list(OUTPUT_POWER_BASIS),
        "ordinary_Z_orders": list(ORDERS),
        "branch_order": [
            "Z0_Phi0_F1_0",
            "Z1_Phi1_F1_0",
            "Z1_L_local_Phi0_F1_0",
            "Z1_Phi0_F1_Z",
        ],
        "formula": {
            "Z0": "-2*Lambda^-2*F0_local*rho*Phi0*W0_chi_F1_0",
            "Z1": "-2*Lambda^-2*F0_local*rho*(Phi1*W0_chi_F1_0 + L_local*Phi0*W0_chi_F1_0 + Phi0*W0_chi_F1_Z)",
        },
        "cells": output_cells,
        "finite_terms": finite_by_order,
        "finite_term_counts": finite_counts,
        "finite_is_signed": True,
        "finite_coefficient_records_unmodified": True,
        "no_coefficient_multiplication": True,
        "no_summation": True,
        "no_integration": True,
        "no_interval_arithmetic": True,
        "rho_same_index_record": True,
        "rho_same_cell_identity_proved": False,
        "coordinate_owner_assumption": {
            "rho_source": "same-index K0 symbolic record; cell identity unverified",
            "x_rho_alignment": "conditional and unverified",
            "fixed_rho_Z": "caller premise",
            "epsilon_Z": "caller premise",
            "cutoff_owner_verified": False,
        },
        "F0_scale_identity": {
            "F0_anchor": "independent formal atom",
            "F0_local": "independent formal atom",
            "identified": False,
        },
        "background_source": {
            "names": list(BACKGROUND_NAMES),
            "records_are_original_copies": True,
            "explicit_background_supplied": background is not None,
            "owner_verified": False,
        },
        "symmetric_kernel_tail": tail,
        "j_ge_2_tail_is_error_only": True,
        "tail_is_not_a_signed_field": True,
        "tail_is_not_merged_with_finite": True,
        "tail_scope": "Kchi potential-kernel truncation only",
        "feedback_reserved": _feedback_reserved(kchi_cells),
        "d_x_F1_retained_for_feedback": True,
        "Z2_retained_for_feedback": True,
        "N1p_included": False,
        "Kchi_feedback_included": False,
        "allword_tail_included": False,
        "pressure_included": False,
        "complete_W0": False,
        "complete_W1": False,
        "complete_m15": False,
        "appended_to_212_ledger": False,
        "ledger_extended": False,
        "point_value_selected": False,
        "source_values_computed": False,
        "numeric_operator_called": False,
    }


__all__ = [
    "SCHEMA",
    "ORDERS",
    "ALL_ORDERS",
    "SOURCE_POWER_BASIS",
    "OUTPUT_POWER_BASIS",
    "SwirlSymbolicInputError",
    "build",
    "bind",
]

# ``build`` is the canonical name; the alias keeps the leaf convenient for
# runners that uniformly call symbolic adapters through ``bind``.
bind = build
