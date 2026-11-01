"""Conditional, source-only bridge for the saved K0 local logarithmic factor.

This module consumes an already materialized result from
``original_k0_symbolic_m15_cells_v1`` and a caller-owned connection packet.
It does not import the producer, decode a receipt, evaluate an interval, or
perform interval arithmetic.  The only operation is a structural substitution
on the affected swirl ``Z1`` records::

    L_local = Lambda * ell_over_lambda

The identity above is deliberately conditional.  A same-owner identity and
the original analytic recipe remain caller premises.  Receipt hashes and
source-log equality are evidence for that premise; they do not prove a
callable identity or equality of endpoint functions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_k0_log_connection_v1"
SYMBOLIC_SCHEMA = "original_k0_symbolic_m15_cells_v1"
SOURCE_POWER_BASIS = (
    "hb",
    "F0_anchor",
    "Pstar_squared",
    "Lambda",
    "F0_local",
    "L_local",
)
ORDERS = ("Z0", "Z1")
CHANNELS = ("axial", "swirl")


class ConnectionInputError(ValueError):
    """Raised when a symbolic record or connection packet is unsafe to bind."""


def _copy(value: Any) -> Any:
    """Copy records without asking them for numeric values."""

    return deepcopy(value)


def _is_int(value: Any) -> bool:
    return type(value) is int


def _is_nonempty(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (str, bytes, bytearray, Mapping, Sequence)):
        return len(value) > 0
    return True


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConnectionInputError(label + " must be a mapping")
    return value


def _binary_interval(value: Any, label: str) -> Mapping[str, Any]:
    """Require exact binary endpoint tuples, without comparing endpoints."""

    record = _require_mapping(value, label)
    for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        endpoint = record.get(key)
        if (
            not isinstance(endpoint, Sequence)
            or isinstance(endpoint, (str, bytes, bytearray))
            or len(endpoint) != 4
            or any(not _is_int(item) for item in endpoint)
        ):
            raise ConnectionInputError(label + " is missing a four-integer binary tuple")
    if record.get("directed") is not True:
        raise ConnectionInputError(label + ".directed must remain true")
    return record


def _powers(value: Any, label: str) -> list[int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != len(SOURCE_POWER_BASIS)
        or any(not _is_int(item) for item in value)
    ):
        raise ConnectionInputError(label + " must contain six integer source powers")
    return [int(item) for item in value]


def _source_family(value: Any, label: str) -> Any:
    if not _is_nonempty(value):
        raise ConnectionInputError(label + " is required")
    if isinstance(value, (Mapping, Sequence)) and not isinstance(value, (str, bytes, bytearray)):
        return _copy(value)
    if isinstance(value, (str, bytes)):
        return _copy(value)
    raise ConnectionInputError(label + " must be a stable source-family record")


def _source_logs(value: Any, label: str) -> Mapping[str, Any]:
    logs = _require_mapping(value, label)
    # The connection contract names four independent source-log records.  Do
    # not infer or synthesize missing logs from a hash or from a receipt path.
    if set(logs) != {"hb", "F0_anchor", "Pstar_squared", "Lambda"} or any(not _is_nonempty(item) for item in logs.values()):
        raise ConnectionInputError(label + " must contain exactly four nonempty source logs")
    return logs


def _receipt_hash(connection: Mapping[str, Any]) -> str:
    candidates = (
        "authenticated_raw_receipt_hash",
        "authenticated_raw_receipt_sha256",
        "raw_receipt_sha256",
        "receipt_sha256",
    )
    value = next((connection.get(key) for key in candidates if key in connection), None)
    if not isinstance(value, str) or len(value) != 64:
        raise ConnectionInputError("authenticated raw receipt hash is required")
    try:
        int(value, 16)
    except (TypeError, ValueError):
        raise ConnectionInputError("authenticated raw receipt hash must be hexadecimal") from None
    return value


def _provenance_correction(connection: Mapping[str, Any]) -> Any:
    if "provenance_correction" not in connection:
        raise ConnectionInputError("runner-declared provenance correction is required")
    correction = connection.get("provenance_correction")
    if not _is_nonempty(correction):
        raise ConnectionInputError("provenance correction must be nonempty")
    if isinstance(correction, Mapping):
        declared = None
        for key in ("declared_by_runner", "runner_declared", "declared"):
            if key in correction:
                declared = correction.get(key)
                break
        if declared is not True and connection.get("provenance_correction_declared") is not True:
            raise ConnectionInputError("provenance correction must be declared by the runner")
    elif connection.get("provenance_correction_declared") is not True:
        raise ConnectionInputError("provenance correction declaration is missing")
    return _copy(correction)


def _connection(connection: Any) -> dict[str, Any]:
    packet = _require_mapping(connection, "connection")
    family = _source_family(packet.get("source_family"), "connection.source_family")
    logs = _source_logs(packet.get("source_logs"), "connection.source_logs")
    ratio_value = packet.get("ell_over_lambda", packet.get("saved_ell_over_lambda"))
    ratio = _binary_interval(ratio_value, "connection.ell_over_lambda")
    receipt_hash = _receipt_hash(packet)
    correction = _provenance_correction(packet)
    domain = packet.get("saved_core_Z_domain")
    if not isinstance(domain, Mapping) or any(
        not isinstance(domain.get(key), (list, tuple)) or len(domain[key]) != 4
        for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple")
    ):
        raise ConnectionInputError("saved core Z domain endpoint records are required")
    if isinstance(correction, Mapping):
        for key in ("receipt_sha256", "raw_receipt_sha256", "source_receipt_sha256"):
            if key in correction and correction.get(key) != receipt_hash:
                raise ConnectionInputError("provenance correction receipt hash differs")
    return {
        "source_family": family,
        "source_logs": _copy(logs),
        "saved_ell_over_lambda": _copy(ratio),
        "authenticated_raw_receipt_hash": receipt_hash,
        "provenance_correction": correction,
        "saved_core_Z_domain": _copy(domain),
        "same_owner_identity_declared": packet.get("same_owner_identity_declared"),
        "analytic_recipe_declared": packet.get("analytic_recipe_declared"),
    }


_IDENTITY_CLAIM_KEYS = {
    "normalized_point_equality",
    "normalized_point_box_equality",
    "normalized_box_equality",
    "box_equality_as_function_identity",
    "point_or_function_equality_claimed",
    "function_identity",
    "callable_identity_proved",
    "endpoint_identity_proved",
}


def _reject_identity_claims(value: Any, label: str) -> None:
    """Reject explicit point/box/function identity claims, not ordinary flags."""

    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key).lower()
            normalized_name = key_text.replace("-", "_").replace(" ", "_")
            if (
                normalized_name in _IDENTITY_CLAIM_KEYS
                or ("normalized" in normalized_name and ("point" in normalized_name or "box" in normalized_name))
                or ("function" in normalized_name and "identity" in normalized_name)
            ) and item is True:
                raise ConnectionInputError(label + " contains an unsupported identity claim")
            # Source records are shallow in the saved schema.  Walking nested
            # mappings keeps an explicit claim from being smuggled in metadata.
            if isinstance(item, Mapping):
                _reject_identity_claims(item, label + "." + str(key))
            elif isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray)):
                for index, nested in enumerate(item):
                    if isinstance(nested, Mapping):
                        _reject_identity_claims(nested, label + "." + str(key) + "[" + str(index) + "]")


def _check_scope_flags(symbolic: Mapping[str, Any]) -> None:
    false_flags = (
        "source_identity_verified",
        "source_values_computed",
        "point_value_selected",
        "interval_arithmetic_performed",
        "integration_performed",
        "numeric_operator_called",
        "N1p_included",
        "Kchi_included",
        "allword_tail_included",
        "appended_to_212_ledger",
        "appended_to_ledger",
        "ledger_consulted",
        "complete_W1",
        "complete_m15_computed",
    )
    for key in false_flags:
        if key in symbolic and symbolic.get(key) is not False:
            raise ConnectionInputError("symbolic scope flag must remain false: " + key)
    scale_identity = symbolic.get("F0_scale_identity")
    if scale_identity is not None:
        if not isinstance(scale_identity, Mapping) or scale_identity.get("identified") is not False:
            raise ConnectionInputError("F0_anchor and F0_local must remain distinct")
    _reject_identity_claims(symbolic, "symbolic")


def _check_unknowns(symbolic: Mapping[str, Any]) -> list[Any]:
    unknowns = symbolic.get("unknown_terms")
    if (
        not isinstance(unknowns, Sequence)
        or isinstance(unknowns, (str, bytes, bytearray))
        or not unknowns
    ):
        raise ConnectionInputError("symbolic unknown_terms are required")
    copied = []
    for index, item in enumerate(unknowns):
        if not isinstance(item, Mapping) or item.get("status") != "unknown":
            raise ConnectionInputError("symbolic unknown_terms[" + str(index) + "] must remain unknown")
        copied.append(_copy(item))
    if not {"N1p", "Kchi", "allword_tail"}.issubset({item.get("name") for item in copied}):
        raise ConnectionInputError("pressure, Kchi and all-word tail must remain explicitly unknown")
    if symbolic.get("unknown_terms_are_not_zero") is not True:
        raise ConnectionInputError("symbolic unknown terms must not be replaced by zero")
    return copied


def _term_identity(term: Mapping[str, Any], label: str) -> dict[str, Any]:
    original = term.get("original_term")
    if original is not None and not isinstance(original, Mapping):
        raise ConnectionInputError(label + ".original_term must be a mapping")
    original = original if isinstance(original, Mapping) else {}
    identity: dict[str, Any] = {}
    for key in (
        "source_cell",
        "target_cell",
        "matrix_sector_identity",
        "W0_sector_identity",
    ):
        value = term.get(key, original.get(key))
        if not _is_nonempty(value):
            raise ConnectionInputError(label + " is missing " + key)
        if key in term and key in original and term.get(key) != original.get(key):
            raise ConnectionInputError(label + " disagrees with original_term." + key)
        identity[key] = _copy(value)
    for key in (
        "source_ids",
        "source_indices",
        "source_cell_provenance",
        "origin",
        "provenance",
        "error_metadata",
        "error_classes",
        "uncertainty",
        "uncertainty_metadata",
        "source_identity_status",
    ):
        if key in term:
            identity[key] = _copy(term.get(key))
    return identity


def _validate_term(term: Any, label: str, cell_index: int, *, target: str) -> dict[str, Any]:
    record = _require_mapping(term, label)
    powers = _powers(record.get("formal_source_powers"), label + ".formal_source_powers")
    if "source_power_basis" in record and list(record.get("source_power_basis")) != list(SOURCE_POWER_BASIS):
        raise ConnectionInputError(label + " has an unknown source-power basis")
    coefficient = _binary_interval(record.get("coefficient"), label + ".coefficient")
    identity = _term_identity(record, label)
    source = identity["source_cell"]
    destination = identity["target_cell"]
    if _is_int(source) and source < 0:
        raise ConnectionInputError(label + ".source_cell is negative")
    if _is_int(destination) and destination < 0:
        raise ConnectionInputError(label + ".target_cell is negative")
    if _is_int(destination) and destination != cell_index:
        raise ConnectionInputError(label + ".target_cell differs from its cell")
    if target == "swirl_Z1" and powers[5] not in (0, 1):
        raise ConnectionInputError(label + " has an unsupported or double L_local power")
    if target != "swirl_Z1" and powers[5] != 0:
        raise ConnectionInputError(label + " carries an L_local factor outside swirl Z1")
    return {
        "record": record,
        "powers": powers,
        "coefficient": coefficient,
        "identity": identity,
    }


def _rows(cell: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    terms = cell.get("terms")
    if not isinstance(terms, Mapping):
        raise ConnectionInputError(label + ".terms is missing")
    return terms


def _bind_cells(symbolic: Mapping[str, Any], packet: Mapping[str, Any]) -> tuple[list[dict[str, Any]], int, list[dict[str, Any]]]:
    ratio = packet["saved_ell_over_lambda"]
    cells = symbolic.get("cells")
    if not isinstance(cells, Sequence) or isinstance(cells, (str, bytes, bytearray)) or len(cells) != 8:
        raise ConnectionInputError("symbolic must contain eight target cells")
    output_cells: list[dict[str, Any]] = []
    all_identities: list[dict[str, Any]] = []
    total = 0
    for cell_index, raw_cell in enumerate(cells):
        cell = _require_mapping(raw_cell, "symbolic.cells[" + str(cell_index) + "]")
        if cell.get("index") != cell_index:
            raise ConnectionInputError("symbolic cell order or index differs")
        rows = _rows(cell, "symbolic.cells[" + str(cell_index) + "]")
        # Validate all supplied sectors so a hidden scale or an unsafe source
        # record cannot pass through merely because it is not substituted.
        for channel_name, channel_rows in rows.items():
            if channel_name not in CHANNELS:
                raise ConnectionInputError("unknown symbolic channel " + str(channel_name))
            if not isinstance(channel_rows, Mapping):
                raise ConnectionInputError("symbolic channel rows must be mappings")
            for order_name, raw_terms in channel_rows.items():
                if order_name not in ORDERS:
                    raise ConnectionInputError("unknown ordinary-Z order " + str(order_name))
                if not isinstance(raw_terms, Sequence) or isinstance(raw_terms, (str, bytes, bytearray)):
                    raise ConnectionInputError("symbolic terms must be sequences")
                target = channel_name + "_" + order_name
                for term_index, raw_term in enumerate(raw_terms):
                    _validate_term(
                        raw_term,
                        "symbolic.cells[" + str(cell_index) + "].terms." + target + "[" + str(term_index) + "]",
                        cell_index,
                        target=target,
                    )
        swirl_rows = rows.get("swirl")
        if not isinstance(swirl_rows, Mapping) or "Z1" not in swirl_rows:
            raise ConnectionInputError("every cell must provide swirl Z1 terms")
        affected: list[dict[str, Any]] = []
        for term_index, raw_term in enumerate(swirl_rows["Z1"]):
            label = "symbolic.cells[" + str(cell_index) + "].terms.swirl_Z1[" + str(term_index) + "]"
            validated = _validate_term(raw_term, label, cell_index, target="swirl_Z1")
            powers = validated["powers"]
            if powers[5] != 1:
                continue
            original = validated["record"]
            if original.get("piece") != "L_local*Phi0*F1_0":
                raise ConnectionInputError(label + " L_local power is not the expected derivative branch")
            saved_powers = original.get("saved_formal_source_powers")
            if (not isinstance(saved_powers, (list, tuple)) or len(saved_powers) != 4
                    or any(not _is_int(x) for x in saved_powers)
                    or powers != [saved_powers[0], saved_powers[1], saved_powers[2], saved_powers[3] - 2, 1, 1]
                    or original.get("signed_factor") != -2):
                raise ConnectionInputError(label + " does not retain the expected saved swirl powers and sign")
            prior_original = original.get("original_term")
            if (not isinstance(prior_original, Mapping)
                    or prior_original.get("coefficient") != original.get("coefficient")
                    or prior_original.get("F0_power") != saved_powers[1]
                    or prior_original.get("Lambda_power") != saved_powers[3]):
                raise ConnectionInputError(label + " changed the saved coefficient")
            original_product_factors = original.get("product_factors")
            if original_product_factors is not None:
                if (
                    not isinstance(original_product_factors, Sequence)
                    or isinstance(original_product_factors, (str, bytes, bytearray))
                ):
                    raise ConnectionInputError(label + ".product_factors must be a sequence")
                if "saved_ell_over_lambda" in original_product_factors:
                    raise ConnectionInputError(label + " already contains saved_ell_over_lambda")
                if sum(item == "L_local" for item in original_product_factors) != 1:
                    raise ConnectionInputError(label + " must contain exactly one L_local factor")
            else:
                raise ConnectionInputError(label + " must retain its original factor list")
            branch = _copy(original)
            branch["channel"] = "swirl"
            branch["ordinary_Z"] = "Z1"
            branch["substitution_status"] = "conditional"
            branch["conditional_substitution"] = True
            branch["substitution_identity"] = "L_local = Lambda * ell_over_lambda"
            branch["original_six_powers"] = _copy(powers)
            branch["substituted_from_formal_source_powers"] = _copy(powers)
            transformed = _copy(powers)
            transformed[3] += 1
            transformed[5] -= 1
            branch["formal_source_powers"] = transformed
            branch["saved_ell_over_lambda"] = _copy(ratio)
            branch["saved_ratio_is_an_enclosure_not_a_selected_value"] = True
            branch["connection_factor_provenance"] = {
                "source_id": "saved_F0_connection:ell_over_lambda",
                "source_receipt_sha256": packet["authenticated_raw_receipt_hash"],
                "provenance_correction": _copy(packet["provenance_correction"]),
            }
            branch["source_ids"] = list(dict.fromkeys(list(branch.get("source_ids", [])) + ["saved_F0_connection:ell_over_lambda"]))
            branch["error_classes"] = list(dict.fromkeys(list(branch.get("error_classes", [])) + ["F0_connection_interval"]))
            branch["uncertainty"] = True
            branch["symbolic_interval_factors"] = ["saved_ell_over_lambda"]
            branch["original_term_identity"] = _copy(validated["identity"])
            if original_product_factors is not None:
                factors = list(original_product_factors)
                if "L_local" in factors:
                    factors.remove("L_local")
                factors.append("saved_ell_over_lambda")
                branch["original_product_factors"] = _copy(original_product_factors)
                branch["product_factors"] = factors
            else:
                branch["product_factors"] = ["saved_ell_over_lambda"]
            if isinstance(original.get("local_factor"), Mapping):
                branch["original_local_factor"] = _copy(original["local_factor"])
                local_factor = _copy(original["local_factor"])
                local_factor["formal_power"] = 0
                local_factor["replaced_by"] = "Lambda*saved_ell_over_lambda"
                branch["local_factor"] = local_factor
            affected.append(branch)
            all_identities.append(_copy(validated["identity"]))
        cell_output = {
            "index": cell_index,
            "target_cell": _copy(cell.get("target_cell", cell_index)),
            "substituted_terms": affected,
            "terms": affected,
            "substituted_term_count": len(affected),
        }
        # Keep compact cell provenance where it is useful, without copying
        # coordinates or the complete 25 MB symbolic input.
        for key in ("source_indices", "triangular_diagonal_source", "axis_endpoint_continuation"):
            if key in cell:
                cell_output[key] = _copy(cell[key])
        output_cells.append(cell_output)
        total += len(affected)
    return output_cells, total, all_identities


def bind(symbolic: Any, connection: Any) -> dict[str, Any]:
    """Bind a saved symbolic result to a caller-owned conditional ratio.

    ``symbolic`` is the previously saved ``real_source_expressions`` mapping
    from ``original_k0_symbolic_m15_cells_v1``.  ``connection`` must carry
    ``ell_over_lambda`` (or its ``saved_ell_over_lambda`` alias), an exact
    source family, exactly four source logs, a 64-character authenticated raw
    receipt hash, and a runner-declared ``provenance_correction`` record.
    """

    source = _require_mapping(symbolic, "symbolic")
    if source.get("schema") != SYMBOLIC_SCHEMA:
        raise ConnectionInputError("symbolic must be original_k0_symbolic_m15_cells_v1 output")
    basis = source.get("source_power_basis")
    if (
        not isinstance(basis, Sequence)
        or isinstance(basis, (str, bytes, bytearray))
        or list(basis) != list(SOURCE_POWER_BASIS)
    ):
        raise ConnectionInputError("symbolic has the wrong six-axis source basis")
    if "formal_source_basis" in source and list(source.get("formal_source_basis")) != list(SOURCE_POWER_BASIS):
        raise ConnectionInputError("symbolic has an unknown formal source basis")
    _check_scope_flags(source)
    unknowns = _check_unknowns(source)
    packet = _connection(connection)
    symbolic_family = _source_family(source.get("source_family"), "symbolic.source_family")
    symbolic_logs = _source_logs(source.get("source_logs"), "symbolic.source_logs")
    if symbolic_family != packet["source_family"]:
        raise ConnectionInputError("symbolic and connection source_family differ")
    if symbolic_logs != packet["source_logs"]:
        raise ConnectionInputError("symbolic and connection source_logs differ")
    cells, count, identities = _bind_cells(source, packet)

    return {
        "schema": SCHEMA,
        "status": "conditional_symbolic_substitution",
        "source_family": _copy(packet["source_family"]),
        "source_logs": _copy(packet["source_logs"]),
        "saved_core_Z_domain": _copy(packet["saved_core_Z_domain"]),
        "saved_ell_over_lambda": _copy(packet["saved_ell_over_lambda"]),
        "connection_evidence": {
            "authenticated_raw_receipt_hash": packet["authenticated_raw_receipt_hash"],
            "provenance_correction": _copy(packet["provenance_correction"]),
        },
        "cells": cells,
        "counts": {
            "cells": [{"index": cell["index"], "substituted": cell["substituted_term_count"]} for cell in cells],
            "affected_terms": count,
        },
        "term_counts": {"swirl": {"Z1": count}},
        "original_term_identities": identities,
        "conditional_identity": {
            "equation": "L_local = Lambda * ell_over_lambda",
            "same_owner_identity_required": True,
            "analytic_recipe_required": True,
            "same_owner_identity_declared_by_connection": packet["same_owner_identity_declared"] is True,
            "analytic_recipe_declared_by_connection": packet["analytic_recipe_declared"] is True,
            "status": "conditional",
            "hash_and_source_log_matches_are_evidence_only": True,
            "callable_identity_proved": False,
            "endpoint_identity_proved": False,
        },
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "F0_scale_identity": {
            "F0_anchor": "independent formal atom",
            "F0_local": "independent formal atom",
            "identified": False,
        },
        "unknown_terms": unknowns,
        "unknown_terms_are_not_zero": True,
        "Kchi_included": False,
        "N1p_included": False,
        "allword_tail_included": False,
        "appended_to_212_ledger": False,
        "source_identity_status": "conditional only; source family/log matches are not identity proof",
        "source_identity_verified": False,
        "normalized_point_or_box_function_identity_claimed": False,
        "source_only": True,
        "source_values_computed": False,
        "point_value_selected": False,
        "interval_arithmetic_performed": False,
        "integration_performed": False,
        "numeric_operator_called": False,
        "old_source_imported": False,
        "old_source_called": False,
        "no_input_repeated": True,
        "scope": "affected swirl Z1 terms with exactly one L_local power only",
    }


__all__ = ["SCHEMA", "SOURCE_POWER_BASIS", "ConnectionInputError", "bind"]
