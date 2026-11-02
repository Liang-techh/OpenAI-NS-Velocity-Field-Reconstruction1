"""Cellwise whole-``Z`` partial reserved/m15 pressure comparison.

The caller supplies runner-authenticated directed intervals.  This leaf does
not read receipts, call an old producer, choose a midpoint/endpoint, merge
cells, or append a historical ledger.  For every source cell it evaluates

``A=exp(-17/10-5*mu/2)/(1+Z**2)``,
``A_Z=-2*Z*exp(-17/10-5*mu/2)/(1+Z**2)**2`` and
``own_k=-r0*(A*theta0)_k-r1*(A*theta1)_k``.

With ``J=-nonpulse`` and independent positive ``P`` (so ``T=-P``), the
partial-minus-independent-negative-pulse result is retained as
``own + r2*J + (1-r2)*P`` in both ordinary ``Z0`` and ``Z1`` rows.  Source
packets and error metadata remain attached to their own cell; no cross-cell
hull is formed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

SCHEMA = "reserved_partial_wholez_v1"
X_DECIMAL = "1.27"
RESPONSE_NAMES = ("r0", "r1", "r2")


class ReservedPartialWholeZInputError(ValueError):
    """Bad saved response/source-packet shape."""


__all__ = ["SCHEMA", "X_DECIMAL", "RESPONSE_NAMES", "ReservedPartialWholeZInputError", "build"]


def _seq(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _is_interval(value: Any) -> bool:
    return getattr(value, "_mpi_", None) is not None or (hasattr(value, "a") and hasattr(value, "b"))


def _raw(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None and len(raw) == 2:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        return getattr(endpoint, "_mpf_", endpoint)
    return value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    if hasattr(c, "make_mpf"):
        return c.make_mpf((_raw(lower, False), _raw(upper, True)))
    return c.mpf([_raw(lower, False), _raw(upper, True)])


def _scalar(c: Any, value: Any, label: str) -> Any:
    if isinstance(value, (bool, float)):
        raise TypeError(f"{label} must not be bool/float")
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not an exact scalar") from exc


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode a directed interval without selecting a representative."""

    if isinstance(value, (bool, float)):
        raise TypeError(f"{label} must be an exact scalar or directed interval")
    if _is_interval(value):
        result = _make(c, value, value)
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact endpoints")
            result = _make(c, tuple(value["lower_exact_mpf_tuple"]), tuple(value["upper_exact_mpf_tuple"]))
        else:
            low_key = next((k for k in ("lower", "lo", "left") if k in value), None)
            high_key = next((k for k in ("upper", "hi", "right") if k in value), None)
            if low_key is not None or high_key is not None:
                if low_key is None or high_key is None:
                    raise TypeError(f"{label} requires lower and upper")
                result = _make(c, _scalar(c, value[low_key], label + ".lower"), _scalar(c, value[high_key], label + ".upper"))
            else:
                nested = next((k for k in ("interval", "directed_interval", "coefficient", "enclosure") if k in value), None)
                if nested is None:
                    raise TypeError(f"{label} is not interval-compatible")
                return _iv(c, value[nested], label + "." + str(nested))
    elif _seq(value) and len(value) == 2:
        result = _make(c, _scalar(c, value[0], label + ".lower"), _scalar(c, value[1], label + ".upper"))
    else:
        point = _scalar(c, value, label)
        result = _make(c, point, point)
    if result.a > result.b:
        raise ValueError(f"{label} has reversed endpoints")
    return result


def _pair(c: Any, value: Any, label: str) -> tuple[Any, Any]:
    if isinstance(value, Mapping):
        k0 = next((k for k in ("Z0", "z0", "C0", "c0", 0) if k in value), None)
        k1 = next((k for k in ("Z1", "z1", "C1", "c1", 1) if k in value), None)
        if k0 is None or k1 is None:
            nested = next((k for k in ("values", "pair", "Z0_Z1", "C0_C1") if k in value), None)
            if nested is None:
                raise TypeError(f"{label} must contain Z0 and Z1")
            return _pair(c, value[nested], label + "." + str(nested))
        return _iv(c, value[k0], label + ".Z0"), _iv(c, value[k1], label + ".Z1")
    if not _seq(value) or len(value) != 2:
        raise TypeError(f"{label} must be a length-2 Z0/Z1 pair")
    return _iv(c, value[0], label + ".Z0"), _iv(c, value[1], label + ".Z1")


def _pair_record(pair: tuple[Any, Any], role: str, formula: str) -> dict[str, Any]:
    return {
        "Z0": pair[0], "Z1": pair[1], "Z0_Z1": pair, "intervals": pair,
        "role": role, "formula": formula, "signed": True,
        "signed_source_term": True, "interval_enclosure": True, "whole_Z_enclosure": True, "function_value": False,
        "error_only": False, "error_term": False, "symmetric_error": False,
    }


def _pick(mapping: Mapping[Any, Any], names: Sequence[Any]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _response(c: Any, value: Any) -> tuple[tuple[Any, Any, Any], bool, Any]:
    marker = None
    if isinstance(value, Mapping):
        marker = _pick(value, ("Z_independent", "z_independent", "fixed_in_Z", "fixed_in_z", "ordinary_Z_independent"))
        if marker is None and isinstance(value.get("metadata"), Mapping):
            marker = _pick(value["metadata"], ("Z_independent", "z_independent", "fixed_in_Z", "fixed_in_z"))
        if all(name in value for name in RESPONSE_NAMES):
            entries = tuple(value[name] for name in RESPONSE_NAMES)
        else:
            nested = next((k for k in ("response", "values", "coefficients") if k in value), None)
            if nested is None or not _seq(value[nested]):
                raise TypeError("response must provide exactly r0, r1, and r2")
            entries = tuple(value[nested])
    else:
        if not _seq(value):
            raise TypeError("response must provide exactly three entries")
        entries = tuple(value)
    if len(entries) != 3:
        raise ValueError(f"response must contain exactly three entries, got {len(entries)}")
    if marker is not True:
        raise ValueError("response must carry an explicit true fixed-in-Z premise")
    checked = tuple(_iv(c, item, f"response[{name}]") for name, item in zip(RESPONSE_NAMES, entries))
    return (checked[0], checked[1], checked[2]), marker is True, value


def _theta(c: Any, cell: Mapping[Any, Any], index: int) -> tuple[tuple[Any, Any], tuple[Any, Any]]:
    label = f"cells[{index}].normalized_theta_inputs"
    if "normalized_theta_inputs" in cell:
        raw = cell["normalized_theta_inputs"]
    elif "normalized_theta_inputs_Z0" in cell and "normalized_theta_inputs_Z1" in cell:
        raw = {"Z0": cell["normalized_theta_inputs_Z0"], "Z1": cell["normalized_theta_inputs_Z1"]}
    else:
        raise KeyError(f"{label} is required")
    if isinstance(raw, Mapping) and any(k in raw for k in ("Z0", "z0", "C0", "c0")):
        row0 = _pair(c, _pick(raw, ("Z0", "z0", "C0", "c0")), label + ".Z0")
        row1 = _pair(c, _pick(raw, ("Z1", "z1", "C1", "c1")), label + ".Z1")
        return (row0[0], row1[0]), (row0[1], row1[1])
    if isinstance(raw, Mapping):
        k0 = next((k for k in ("theta0", "theta_0", 0) if k in raw), None)
        k1 = next((k for k in ("theta1", "theta_1", 1) if k in raw), None)
        if k0 is None or k1 is None:
            raise TypeError(f"{label} must provide Z0/Z1 or theta0/theta1")
        return _pair(c, raw[k0], label + ".theta0"), _pair(c, raw[k1], label + ".theta1")
    if not _seq(raw) or len(raw) != 2:
        raise TypeError(f"{label} must contain two normalized theta inputs")
    return _pair(c, raw[0], label + "[0]"), _pair(c, raw[1], label + "[1]")


def _source(cell: Mapping[Any, Any], index: int) -> tuple[Any, Any, Any, dict[str, Any]]:
    packet = cell.get("source_packet", cell.get("packet", cell))
    nested = packet if isinstance(packet, Mapping) else {}
    source_id = _pick(cell, ("source_id", "cell_id", "source_cell_id", "id", "packet_id"))
    source_id = source_id if source_id is not None else _pick(nested, ("source_id", "cell_id", "source_cell_id", "id", "packet_id"))
    source_path = _pick(cell, ("source_path", "source_packet_path", "receipt_input_path", "input_path", "path"))
    source_path = source_path if source_path is not None else _pick(nested, ("source_path", "source_packet_path", "receipt_input_path", "input_path", "path"))
    if source_id is None or source_path is None:
        raise KeyError(f"cells[{index}] requires original source path and ID")
    if (isinstance(source_id, str) and not source_id) or (isinstance(source_path, str) and not source_path) or (_seq(source_path) and not source_path):
        raise ValueError(f"cells[{index}] source path/ID may not be empty")
    error_keys = ("error", "errors", "error_packet", "error_classes", "source_error_classes", "uncertainty", "error_only", "symmetric_error")
    errors = {key: cell[key] for key in error_keys if key in cell}
    if isinstance(nested, Mapping):
        errors.update({key: nested[key] for key in error_keys if key in nested and key not in errors})
    return source_id, source_path, packet, errors


def _cells(value: Any) -> list[Mapping[Any, Any]]:
    if isinstance(value, Mapping):
        nested = next((k for k in ("cells", "whole_Z_cells", "whole_z_cells", "whole_Z", "source_cells") if k in value), None)
        if nested is None:
            raise TypeError("cells mapping must contain whole-Z source cells")
        value = value[nested]
    if not _seq(value) or not value or any(not isinstance(cell, Mapping) for cell in value):
        raise TypeError("cells must be a nonempty sequence of source mappings")
    return list(value)


def _joint_or_pulse(c: Any, value: Any, label: str) -> tuple[Any, Any, Any]:
    if label == "pulse" and isinstance(value, Mapping):
        if value.get("is_P1", value.get("P1", False)) is True or str(value.get("label", "")).lower() == "p1":
            raise ValueError("independent pulse P must remain distinct from P1")
    pair = _pair(c, value, label)
    return pair[0], pair[1], value


def build(c: Any, *, response: Any, cells: Any, mu: Any, joint: Any, pulse: Any) -> dict[str, Any]:
    """Return per-cell signed ``own``, ``r2*J``, ``(1-r2)*P`` and result."""

    if c is None or not all(hasattr(c, name) for name in ("mpf", "exp", "make_mpf")):
        raise TypeError("c must be the runner-authenticated mpmath interval context")
    response_values, response_explicit, raw_response = _response(c, response)
    if isinstance(response, Mapping) and response.get("x") is not None:
        x_supplied = _iv(c, response["x"], "response.x")
        if x_supplied._mpi_ != _iv(c, X_DECIMAL, "x")._mpi_:
            raise ValueError("response.x must be the fixed x=1.27 packet")
    source_cells = _cells(cells)
    mu_iv = _iv(c, mu, "mu")
    joint_z0, joint_z1, raw_joint = _joint_or_pulse(c, joint, "joint")
    pulse_z0, pulse_z1, raw_pulse = _joint_or_pulse(c, pulse, "pulse")
    r0, r1, r2 = response_values
    one_minus_r2 = 1 - r2
    e = c.exp(-c.mpf(17) / 10 - c.mpf(5) * mu_iv / 2)
    results: list[dict[str, Any]] = []
    for index, cell in enumerate(source_cells):
        if "clipped_Z" not in cell:
            raise KeyError(f"cells[{index}].clipped_Z is required")
        z = _iv(c, cell["clipped_Z"], f"cells[{index}].clipped_Z")
        theta0, theta1 = _theta(c, cell, index)
        source_id, source_path, source_packet, source_errors = _source(cell, index)
        denominator = 1 + z ** 2
        A = e / denominator
        A_Z = -2 * z * e / (denominator * denominator)
        products = (
            (A * theta0[0], A_Z * theta0[0] + A * theta0[1]),
            (A * theta1[0], A_Z * theta1[0] + A * theta1[1]),
        )
        own = (-(r0 * products[0][0] + r1 * products[1][0]), -(r0 * products[0][1] + r1 * products[1][1]))
        joint_term = (r2 * joint_z0, r2 * joint_z1)
        pulse_term = (one_minus_r2 * pulse_z0, one_minus_r2 * pulse_z1)
        result_pair = (own[0] + joint_term[0] + pulse_term[0], own[1] + joint_term[1] + pulse_term[1])
        own_rec = _pair_record(own, "own", "-r0*(A*theta0)_k-r1*(A*theta1)_k")
        joint_rec = _pair_record(joint_term, "r2_times_J", "r2*J")
        pulse_rec = _pair_record(pulse_term, "one_minus_r2_times_P", "(1-r2)*P")
        result_rec = _pair_record(result_pair, "partial_minus_independent_negative_pulse", "own+r2*J+(1-r2)*P")
        record = dict(cell)
        record.update({
            "cell_index": index, "source_id": source_id, "source_path": source_path,
            "source_packet": source_packet, "source_error_metadata": source_errors,
            "source_errors_retained": True, "clipped_Z": z, "Z_interval": z,
            "normalized_theta_inputs_Z0": (theta0[0], theta1[0]),
            "normalized_theta_inputs_Z1": (theta0[1], theta1[1]),
            "A": A, "A_Z": A_Z, "amplitude_A": A, "amplitude_A_Z": A_Z,
            "theta_products": {"theta0": _pair_record(products[0], "A_times_theta0", "A*theta0 product rule"), "theta1": _pair_record(products[1], "A_times_theta1", "A*theta1 product rule")},
            "own": own_rec, "own_Z0_Z1": own, "own_C0_C1": own,
            "joint_term": joint_rec, "r2_times_J": joint_rec, "r2_times_joint": joint_rec, "joint_term_Z0_Z1": joint_term,
            "pulse_term": pulse_rec, "one_minus_r2_times_P": pulse_rec, "complement_r2_times_pulse": pulse_rec, "pulse_term_Z0_Z1": pulse_term,
            "result": result_rec, "result_Z0_Z1": result_pair, "partial_minus_independent_negative_pulse": result_rec,
            "signed_terms": {"own": own_rec, "r2_times_J": joint_rec, "one_minus_r2_times_P": pulse_rec},
            "signed_terms_separate": True, "error_terms_separate": False, "formula_identity": "own+r2*J+(1-r2)*P",
            "r2_fixed_in_Z": True, "P_is_not_P1": True,
        })
        results.append(record)
    x = _iv(c, X_DECIMAL, "x")
    return {
        "schema": SCHEMA, "status": "conditional_partial_minus_independent_negative_pulse", "conditional": True,
        "runner_authenticated_interval_context_required": True, "runner_authenticated_interval_context_consumed": True,
        "conditional_on": ("runner-authenticated interval context", "saved fixed-in-Z response", "caller-supplied whole-Z packets"),
        "x": x, "x_decimal": X_DECIMAL, "x_fixed": True, "Z_domain": "open", "open_Z_domain": True,
        "ordinary_Z_orders": ["Z0", "Z1"], "ordinary_Z1_product_rule": "A_Z*theta_Z0+A*theta_Z1",
        "A_formula": "exp(-17/10-5*mu/2)/(1+Z^2)", "A_Z_formula": "-2*Z*exp(-17/10-5*mu/2)/(1+Z^2)^2",
        "response": raw_response, "response_intervals": response_values, "response_count": 3,
        "response_Z_independent": True, "response_Z_independent_explicit": response_explicit,
        "response_fixed_in_Z_proved": False, "mu": mu_iv, "mu_interval": mu_iv,
        "joint": raw_joint, "joint_Z0_Z1": (joint_z0, joint_z1), "joint_definition": "J=-nonpulse",
        "joint_is_whole_Z_enclosure": True, "joint_is_function_value": False, "joint_error_separate": False,
        "pulse": raw_pulse, "pulse_Z0_Z1": (pulse_z0, pulse_z1), "pulse_sign_assumed": False,
        "pulse_is_whole_Z_enclosure": True, "pulse_is_function_value": False, "pulse_error_separate": False,
        "pulse_channel": "m15 pressure pulse channel", "P_is_not_P1": True, "P1_distinct": True,
        "original_owner_proved": False, "original_owner_premise_unproved": True, "fixed_in_Z_premise_unproved": True,
        "finite_N_closure": False, "not_finite_N_closure": True, "physical_Rp": False, "not_physical_Rp": True,
        "appended_to_212_ledger": False, "ledger_appended": False, "old_errors_deleted": False,
        "source_packets_retained": True, "source_paths_retained": True, "source_ids_retained": True, "source_errors_retained": True,
        "signed_terms_separate_from_errors": False, "error_terms_separate": False, "separate_numerical_error_bounds_computed": False, "midpoint_values_used": False,
        "midpoint_used": False, "endpoint_values_selected": False, "data_midpoint_or_error_replacement": False,
        "error_replacement": False, "errors_removed": False,
        "cross_cell_hull": False, "cross_cell_hull_performed": False, "cells": results, "cell_results": results,
        "cell_count": len(results), "input_packets": {"response": raw_response, "cells": cells, "mu": mu, "joint": raw_joint, "pulse": raw_pulse},
        "old_212_ledger_untouched": True, "actual_function_owner": "caller-supplied; not independently proved",
    }

