"""Pressure/width component lift of the flat-shape defect kernel.

The scalar flat-shape kernel is analytic in its amplitude ``B``.  This
adapter evaluates its Taylor coefficients at the scalar constant atom
``B0 = B[0, 0]`` and composes those coefficients with a rectangular
``PressureWidthJet``.  The pressure-tail and exit-width atoms stay separate
throughout the composition; ``B`` is never first evaluated at ``(1, 1)``.

This is a finite local component calculation.  It preserves every Taylor
term and every saddle receipt, while making no enclosure or five-moment field
claim.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_flat_shape_defect import (
        _work_precision,
        evaluate_flat_shape_defect,
        evaluate_flat_shape_derivative,
    )
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_flat_shape_defect import (  # type: ignore
        _work_precision,
        evaluate_flat_shape_defect,
        evaluate_flat_shape_derivative,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet  # type: ignore


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _orders(jet: PressureWidthJet) -> tuple[int, int]:
    return int(jet.pressure_order), int(jet.width_order)


def _zero(jet: PressureWidthJet) -> PressureWidthJet:
    return PressureWidthJet(
        0,
        pressure_order=jet.pressure_order,
        width_order=jet.width_order,
    )


def _signed_log_value(record: Any, *, derivative: bool = False) -> mp.mpf:
    """Read a scalar value from either a numeric or signed-log record."""

    if isinstance(record, Mapping):
        if derivative:
            direct_names = ("derivative_value", "dI_dB", "derivative")
        else:
            direct_names = ("value", "integral", "I")
        for name in direct_names:
            if name in record and not isinstance(record[name], Mapping):
                return _mp(record[name])
        sign_name = "derivative_sign" if derivative else "sign"
        log_name = "derivative_log_abs" if derivative else "log_abs"
        sign = record.get(sign_name)
        log_abs = record.get(log_name)
        if sign is not None and log_abs is not None:
            sign = int(sign)
            if sign == 0:
                return mp.mpf(0)
            return sign * mp.exp(_mp(log_abs))
        for name in ("signed_log", "signed_value"):
            if name in record:
                return _signed_log_value(record[name], derivative=derivative)
    names = (
        ("derivative_value", "dI_dB", "derivative")
        if derivative
        else ("value", "integral", "I")
    )
    for name in names:
        if hasattr(record, name):
            return _mp(getattr(record, name))
    sign = getattr(record, "derivative_sign" if derivative else "sign", None)
    log_abs = getattr(record, "derivative_log_abs" if derivative else "log_abs", None)
    if sign is not None and log_abs is not None:
        sign = int(sign)
        return mp.mpf(0) if sign == 0 else sign * mp.exp(_mp(log_abs))
    return _mp(record)


def _derivative_value(record: Any, order: int) -> mp.mpf:
    """Extract the actual B derivative represented by one evaluator record."""

    if order == 0:
        return _signed_log_value(record)
    # The base n = 1 record names its sensitivity fields ``derivative_*``;
    # the dedicated n >= 2 API returns an ordinary signed-log value record.
    if isinstance(record, Mapping) and (
        "derivative_value" in record or "derivative_log_abs" in record
    ):
        return _signed_log_value(record, derivative=True)
    return _signed_log_value(record)


def _record_receipt(record: Any, *, derivative_order: int = 0) -> dict[str, Any]:
    """Keep the requested derivative's mode/sign receipt.

    The base n = 1 record contains both the value-at-B0 receipt and a
    ``derivative_*`` receipt.  Normalize the latter to ``mode``/``width``
    when that record is used for the first derivative term, so each retained
    Taylor term points to its own saddle.
    """

    if not isinstance(record, Mapping):
        return {"value": _mp(record)}
    if derivative_order == 1 and "derivative_mode" in record:
        return {
            name: record[source]
            for name, source in (
                ("sign", "derivative_sign"),
                ("log_abs", "derivative_log_abs"),
                ("value", "derivative_value"),
                ("mode", "derivative_mode"),
                ("width", "derivative_width"),
                ("mode_method", "derivative_mode_method"),
                ("quadrature_order", "quadrature_order"),
                ("saddle_window", "saddle_window"),
                ("work_precision", "work_precision"),
            )
            if source in record
        }
    names = (
        "sign",
        "log_abs",
        "value",
        "derivative_sign",
        "derivative_log_abs",
        "derivative_value",
        "mode",
        "width",
        "mode_method",
        "derivative_mode",
        "derivative_width",
        "derivative_mode_method",
        "quadrature_order",
        "saddle_window",
        "work_precision",
    )
    return {name: record[name] for name in names if name in record}


def _scalar_work_precision(T: Any, precision: int) -> int:
    try:
        return int(_work_precision(_mp(T), int(precision)))
    except (AttributeError, TypeError, ValueError):
        T_value = _mp(T)
        digits = 0
        if T_value > 1:
            digits = max(0, int(mp.floor(mp.log10(T_value))) + 1)
        return max(320 + digits, int(precision))


def compose_flat_shape_defect(
    k: Any,
    m: Any,
    B: PressureWidthJet,
    B_Z: PressureWidthJet,
    T: Any,
    *,
    B_ZZ: PressureWidthJet | None = None,
    precision: int = 260,
    order: int = 32,
    window: Any = 24,
) -> dict[str, Any]:
    """Compose the scalar defect and its Z derivatives over a finite ring.

    ``B`` is expanded around its scalar constant atom ``B0``.  If ``P`` and
    ``W`` are the retained pressure and width orders, ``N = P + W`` and the
    returned value uses derivatives through ``I^(N)(B0)``.  The returned
    tangent is the chain-rule product through ``I^(N+1)(B0) * B_Z``.  When
    ``B_ZZ`` is supplied, the returned ``second_terms`` apply the full second
    chain rule through ``I^(N+2)(B0)``.
    """

    if not isinstance(B, PressureWidthJet) or not isinstance(B_Z, PressureWidthJet):
        raise TypeError("B and B_Z must be PressureWidthJet instances")
    if _orders(B) != _orders(B_Z):
        raise ValueError("B and B_Z must have matching pressure/width orders")
    if B_ZZ is not None:
        if not isinstance(B_ZZ, PressureWidthJet):
            raise TypeError("B_ZZ must be a PressureWidthJet instance")
        if _orders(B) != _orders(B_ZZ):
            raise ValueError("B and B_ZZ must have matching pressure/width orders")
    if int(precision) < 80:
        raise ValueError("precision must be at least 80")
    if int(order) < 8:
        raise ValueError("order must be at least 8")

    pressure_order, width_order = _orders(B)
    truncation_order = pressure_order + width_order
    derivative_cache_order = truncation_order + (2 if B_ZZ is not None else 1)
    workdps = _scalar_work_precision(T, precision)
    with mp.workdps(workdps):
        B0 = B.component(0, 0)
        baseline = PressureWidthJet(
            B0,
            pressure_order=pressure_order,
            width_order=width_order,
        )
        delta = B - baseline
        one = PressureWidthJet(
            1,
            pressure_order=pressure_order,
            width_order=width_order,
        )

        # Keep the scalar derivative records separate.  The n = 1 record is
        # already returned by the base API, so it is reused rather than
        # recomputed through a second saddle search.
        scalar_records: dict[int, Any] = {}
        base_record = evaluate_flat_shape_defect(
            k,
            m,
            B0,
            T,
            precision=precision,
            order=order,
            window=window,
        )
        scalar_records[0] = base_record
        scalar_records[1] = base_record
        for derivative_order in range(2, derivative_cache_order + 1):
            scalar_records[derivative_order] = evaluate_flat_shape_derivative(
                k,
                m,
                B0,
                T,
                derivative_order,
                precision=precision,
                order=order,
                window=window,
            )

        derivative_values = {
            derivative_order: _derivative_value(
                scalar_records[derivative_order], derivative_order
            )
            for derivative_order in range(derivative_cache_order + 1)
        }
        value_terms: dict[int, PressureWidthJet] = {}
        tangent_terms: dict[int, PressureWidthJet] = {}
        second_terms: dict[int, PressureWidthJet] | None = None
        second_chain_terms: dict[int, PressureWidthJet] | None = None
        second_curvature_terms: dict[int, PressureWidthJet] | None = None
        if B_ZZ is not None:
            # Keep the two second-chain contributions independently inspectable:
            # I''(B) * B_Z**2 and I'(B) * B_ZZ.
            second_terms = {}
            second_chain_terms = {}
            second_curvature_terms = {}
            B_Z_squared = B_Z * B_Z
        delta_power = one
        for term_order in range(truncation_order + 1):
            factorial = mp.factorial(term_order)
            value_terms[term_order] = delta_power * (
                derivative_values[term_order] / factorial
            )
            tangent_terms[term_order] = (
                delta_power
                * (derivative_values[term_order + 1] / factorial)
                * B_Z
            )
            if B_ZZ is not None:
                coefficient = delta_power / factorial
                second_chain_terms[term_order] = (
                    coefficient * derivative_values[term_order + 2] * B_Z_squared
                )
                second_curvature_terms[term_order] = (
                    coefficient * derivative_values[term_order + 1] * B_ZZ
                )
                second_terms[term_order] = (
                    second_chain_terms[term_order]
                    + second_curvature_terms[term_order]
                )
            delta_power = delta_power * delta

        value = _zero(B)
        tangent = _zero(B)
        for term_order in range(truncation_order + 1):
            value = value + value_terms[term_order]
            tangent = tangent + tangent_terms[term_order]

        result: dict[str, Any] = {
            "value": value,
            "tangent": tangent,
            "value_terms": value_terms,
            "tangent_terms": tangent_terms,
            "derivative_cache": scalar_records,
            "derivative_values": derivative_values,
            "B0": B0,
            "delta": delta,
            "B": B,
            "B_Z": B_Z,
            "k": _mp(k),
            "m": _mp(m),
            "T": _mp(T),
            "precision": int(precision),
            "quadrature_order": int(order),
            "saddle_window": _mp(window),
            "pressure_order": pressure_order,
            "width_order": width_order,
            "derivative_truncation_order": derivative_cache_order,
            "value_truncation_order": truncation_order,
            "term_receipts": {
                term_order: _record_receipt(
                    scalar_records[term_order], derivative_order=term_order
                )
                for term_order in range(truncation_order + 1)
            },
            "tangent_term_receipts": {
                term_order: _record_receipt(
                    scalar_records[term_order + 1],
                    derivative_order=term_order + 1,
                )
                for term_order in range(truncation_order + 1)
            },
            "metadata": {
                "pressure_order": pressure_order,
                "width_order": width_order,
                "value_truncation_order": truncation_order,
                "derivative_cache_order": derivative_cache_order,
                "work_precision": workdps,
                "quadrature_order": int(order),
                "saddle_window": _mp(window),
                "rectangular_pressure_width_truncation": True,
                "pressure_parameter_distinct_from_temporal": True,
                "terms_preserved_before_sum": True,
                "quadrature_remainder_enclosed": False,
                "all_five_moment_claim": False,
                "global_field_installed": False,
            },
            "rectangular_pressure_width_truncation": True,
            "terms_preserved_before_sum": True,
            "pressure_parameter_distinct_from_temporal": True,
            "quadrature_remainder_enclosed": False,
            "all_five_moment_claim": False,
            "global_field_installed": False,
            "limitations": (
                "Finite signed-log saddle receipts composed over a rectangular "
                "pressure/width ring; tails and quadrature remainder are unenclosed."
            ),
        }

        if B_ZZ is not None:
            assert (
                second_terms is not None
                and second_chain_terms is not None
                and second_curvature_terms is not None
            )
            fullsecond = _zero(B)
            for term_order in range(truncation_order + 1):
                fullsecond = fullsecond + second_terms[term_order]
            result.update(
                {
                    "B_ZZ": B_ZZ,
                    "second_terms": second_terms,
                    "second_chain_terms": second_chain_terms,
                    "second_curvature_terms": second_curvature_terms,
                    "fullsecond": fullsecond,
                    "second_term_receipts": {
                        term_order: _record_receipt(
                            scalar_records[term_order + 2],
                            derivative_order=term_order + 2,
                        )
                        for term_order in range(truncation_order + 1)
                    },
                    "metadata": {
                        **result["metadata"],
                        "second_Z_derivative": True,
                        "second_truncation_order": truncation_order,
                        "second_derivative_cache_order": derivative_cache_order,
                        "second_terms_preserved_before_sum": True,
                        "second_remainder_enclosed": False,
                    },
                    "second_Z_derivative": True,
                    "second_truncation_order": truncation_order,
                    "second_derivative_cache_order": derivative_cache_order,
                    "second_terms_preserved_before_sum": True,
                    "second_remainder_enclosed": False,
                }
            )
        return result


compose = compose_flat_shape_defect


__all__ = ["compose_flat_shape_defect", "compose"]
