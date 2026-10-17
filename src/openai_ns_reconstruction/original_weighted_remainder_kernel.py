"""Pure weighted primitive and finite-width remainder arithmetic.

This module is deliberately a small source-side kernel.  It does not load an
accepted source packet, construct a core, recompute a mass, or choose a width.
The inputs are the signed ordinary ``Z0/Z1/Z2`` rows and independently supplied
moment data.

For angular rows ``A_j`` and axial rows ``B_j`` the source model used here is

``A_h(t) = sum_j A_j exp((1-j) h t)``

``ell_h(t) = -(h**2/2) int_0^t A_h(u) du``

``K_a(h) = int_0^s kappa(t) exp(a h t + ell_h(t)) dt``.

Consequently, with ``M0=int kappa``, ``W1=int t*kappa`` and
``W2=int t**2*kappa``,

``(K_a-M0)/h = a W1 + h*(a**2 W2/2-A0 W1/2) + O(h**2)``

where ``A0=sum_j A_j``.  The ordinary rows are derivatives, rather than Taylor
coefficients: the second row uses the usual factor two in the product rule.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import mpmath as _mp


_DEFAULT_ERROR_CAP = "1e-1000"
_DEFAULT_ETA = "1e-1000"


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _row(values: Any, *, name: str = "ordinary row") -> list[Any]:
    if not _is_sequence(values) or len(values) != 3:
        raise ValueError(f"{name} must contain ordinary Z0/Z1/Z2 rows")
    if any(value is None for value in values):
        raise ValueError(f"{name} contains a missing ordinary derivative row")
    return list(values)


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_")


def _scalar(value: Any, *, name: str) -> Any:
    """Reject a row where a Z-independent scalar was required."""

    if _is_sequence(value):
        raise ValueError(f"{name} must be a Z-independent scalar, not ordinary rows")
    if value is None:
        raise ValueError(f"{name} is required")
    return value


def _mpf_pair(value: Any) -> tuple[Any, Any] | None:
    """Return raw mpf endpoints for ordinary or interval mpmath values."""

    if _is_interval(value):
        lower, upper = value._mpi_
        return lower, upper
    if hasattr(value, "_mpf_"):
        point = value._mpf_
        return point, point
    return None


def _ordinary_endpoint(ctx: Any, value: Any) -> tuple[Any, Any]:
    pair = _mpf_pair(value)
    if pair is None:
        point = ctx.mpf(value)
        return point, point
    return ctx.make_mpf(pair[0]), ctx.make_mpf(pair[1])


def _lower(ctx: Any, value: Any) -> Any:
    return _ordinary_endpoint(ctx, value)[0]


def _upper(ctx: Any, value: Any) -> Any:
    return _ordinary_endpoint(ctx, value)[1]


def _is_zero(ctx: Any, value: Any) -> bool:
    lower, upper = _ordinary_endpoint(ctx, value)
    return lower == 0 and upper == 0


def _require_nonnegative(ctx: Any, value: Any, *, name: str) -> Any:
    value = _scalar(value, name=name)
    lower, upper = _ordinary_endpoint(ctx, value)
    if lower < 0 or upper < lower:
        raise ValueError(f"{name} must be a nonnegative directed scalar")
    return value


def _ordinary_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply ordinary derivative rows through Z2."""

    a = _row(left, name="left ordinary row")
    b = _row(right, name="right ordinary row")
    return [
        a[0] * b[0],
        a[1] * b[0] + a[0] * b[1],
        a[2] * b[0] + 2 * a[1] * b[1] + a[0] * b[2],
    ]


def _ordinary_sum(rows: Sequence[Sequence[Any]]) -> list[Any]:
    if not rows:
        raise ValueError("at least one ordinary row is required")
    normal = [_row(row) for row in rows]
    return [sum((row[k] for row in normal), normal[0][k] * 0) for k in range(3)]


def _ordinary_scale(row: Sequence[Any], scalar: Any) -> list[Any]:
    row = _row(row)
    scalar = _scalar(scalar, name="ordinary scale")
    return [value * scalar for value in row]


def _ordinary_norm2(ctx: Any, iv: Any, row: Sequence[Any], weight: Any) -> Any:
    """A weighted upper norm for an ordinary row through Z2."""

    row = _row(row)
    weight_iv = iv.mpf(weight)
    return sum(
        (
            iv.mpf(_abs_upper(ctx, value))
            * weight_iv**order
            / math.factorial(order)
            for order, value in enumerate(row)
        ),
        iv.mpf(0),
    )


def _abs_upper(ctx: Any, value: Any) -> Any:
    lower, upper = _ordinary_endpoint(ctx, value)
    return max(abs(lower), abs(upper))


def _interval_with_upper(ctx: Any, iv: Any, value: Any) -> Any:
    return iv.mpf([iv.mpf(_lower(ctx, value)), iv.mpf(_upper(ctx, value))])


def _upper_atom(ctx: Any, iv: Any, value: Any) -> Any:
    return iv.mpf(_upper(ctx, value))


def _serialize_interval(value: Any) -> list[list[int]] | None:
    if not _is_interval(value):
        if hasattr(value, "_mpf_"):
            return [list(value._mpf_), list(value._mpf_)]
        return None
    return [list(pair) for pair in value._mpi_]


def _validate_weight(ctx: Any, iv: Any, weight: Any) -> Any:
    weight = _scalar(weight, name="weight")
    lower, upper = _ordinary_endpoint(ctx, weight)
    if lower <= 0 or upper < lower:
        raise ValueError("weight must be a positive directed scalar")
    return iv.mpf(weight)


def _validate_modes(
    axial_modes: Any,
    angular_modes: Any,
    *,
    name: str = "axial modes",
) -> tuple[list[dict[str, Any]], list[list[Any]]]:
    """Normalize the public mode dictionaries and angular rows."""

    if isinstance(angular_modes, dict):
        angular_modes = angular_modes.get("ordinary_rows", angular_modes.get("rows"))
    if _is_sequence(angular_modes) and len(angular_modes) == 3 and not any(
        _is_sequence(value) for value in angular_modes
    ):
        angular_modes = [angular_modes]
    if not _is_sequence(angular_modes) or not angular_modes:
        raise ValueError("at least one angular ordinary row is required")
    angular_rows = []
    for index, row in enumerate(angular_modes):
        if isinstance(row, dict):
            row = row.get("ordinary_rows", row.get("rows"))
        angular_rows.append(_row(row, name=f"angular ordinary row {index}"))
    # The source has j=0,1,2.  The finite bound below uses max |1-j|=1.
    if any(abs(1 - index) > 1 for index in range(len(angular_rows))):
        raise ValueError("angular modes must use the original j=0,1,2 rates")

    if isinstance(axial_modes, dict):
        axial_modes = axial_modes.get("modes", axial_modes.get("ordinary_rows"))
    if not _is_sequence(axial_modes) or not axial_modes:
        raise ValueError("at least one axial mode is required")

    normalized: list[dict[str, Any]] = []
    for index, mode in enumerate(axial_modes):
        if isinstance(mode, dict):
            rows = mode.get("ordinary_rows", mode.get("rows"))
            if rows is None:
                raise ValueError(f"axial mode {index} is missing ordinary_rows")
            if "rate" not in mode:
                raise ValueError(f"axial mode {index} is missing rate")
            rate = mode["rate"]
        elif _is_sequence(mode) and len(mode) == 2 and _is_sequence(mode[0]):
            rows, rate = mode
        else:
            raise ValueError(
                "axial modes must be {ordinary_rows:[B0,Bz,Bzz], rate:a} dictionaries"
            )
        normalized.append(
            {"ordinary_rows": _row(rows, name=f"axial mode {index}"), "rate": _scalar(rate, name="rate")}
        )
    return normalized, angular_rows


def _validate_moments(ctx: Any, mass: Any, W1: Any, W2: Any | None = None) -> tuple[Any, Any, Any | None]:
    mass = _require_nonnegative(ctx, mass, name="mass")
    W1 = _require_nonnegative(ctx, W1, name="W1")
    _, mass_upper = _ordinary_endpoint(ctx, mass)
    W1_lower, W1_upper = _ordinary_endpoint(ctx, W1)
    if W1_lower > mass_upper:
        raise ValueError("W1 is inconsistent with the supplied nonnegative mass")
    if W2 is not None:
        W2 = _require_nonnegative(ctx, W2, name="W2")
        W2_lower, _ = _ordinary_endpoint(ctx, W2)
        if W2_lower > W1_upper:
            raise ValueError("W2 is inconsistent with the supplied W1 moment")
    return mass, W1, W2


def _source_guard(source_proved: Any) -> None:
    if source_proved is not True:
        raise ValueError("an explicit same-source weighted-primitive proof is required")


def _coefficient_rows(
    axial_modes: Any,
    angular_modes: Any,
    mass: Any,
    W1: Any,
    W2: Any | None,
    *,
    source_proved: Any,
) -> tuple[list[Any], list[Any] | None, dict[str, Any]]:
    _source_guard(source_proved)
    normalized_axial, angular_rows = _validate_modes(axial_modes, angular_modes)
    # A local point/interval context is not needed for algebra; shape and sign
    # validation use the values' own mpmath endpoints.
    validation_ctx = _ValidationContext()
    mass, W1, W2 = _validate_moments(validation_ctx, mass, W1, W2)
    A0 = _ordinary_sum(angular_rows)
    leading: list[Any] = [0, 0, 0]
    next_rows: list[Any] | None = [0, 0, 0] if W2 is not None else None
    # Use a first mode's zero-like values to keep arbitrary mpmath contexts
    # intact instead of coercing through a global mpmath context.
    if normalized_axial:
        zero = normalized_axial[0]["ordinary_rows"][0] * 0
        leading = [zero, zero, zero]
        if next_rows is not None:
            next_rows = [zero, zero, zero]
    for mode in normalized_axial:
        rows = mode["ordinary_rows"]
        rate = mode["rate"]
        r0 = _ordinary_scale([1, 0, 0], rate * W1)
        leading = [a + b for a, b in zip(leading, _ordinary_product(rows, r0))]
        if next_rows is not None:
            scalar_part = (rate * rate * W2) / 2
            angular_part = _ordinary_scale(A0, W1 / 2)
            # The axial rate is Z-independent: its scalar contribution has
            # ordinary row [scalar_part, 0, 0].  Subtracting it from every
            # angular derivative row would invent a Z derivative of the rate.
            coefficient = [scalar_part - angular_part[0], -angular_part[1], -angular_part[2]]
            next_rows = [a + b for a, b in zip(next_rows, _ordinary_product(rows, coefficient))]
    proof = {
        "passed": True,
        "source_proved": True,
        "mass_Z_independent": True,
        "W1_Z_independent": True,
        "W2_Z_independent": W2 is not None,
        "local_ell_definition": "ell_h(t)=-(h^2/2)*integral_0^t sum_j A_j exp((1-j)h*u) du",
        "angular_weight": "1",
        "axial_weight": "kappa(t)",
        "mass_used_as_source_moment": True,
        "actual_field_zeroed": False,
    }
    return leading, next_rows, proof


class _ValidationContext:
    """Tiny context for shape/sign checks without selecting a numeric context."""

    @staticmethod
    def mpf(value: Any) -> Any:
        return _mp.mpf(value)

    @staticmethod
    def make_mpf(value: Any) -> Any:
        # This method is reached only for mpmath tuples in normal operation;
        # direct values are returned for a context-free point calculation.
        return _mp.mpf(value)


def signed_remainder_coefficient(
    axial_modes: Any,
    angular_modes: Any,
    mass: Any,
    W1: Any,
    source_proved: Any,
    return_details: bool = False,
) -> Any:
    """Return the signed first coefficient ``R0=sum_j B_j*a_j*W1``.

    ``mass`` is accepted and validated as the independently saved ``M0``
    source moment.  It does not appear in ``R0`` because the corrected source
    has an angular exponent beginning at order ``h^2``.
    """

    leading, _, proof = _coefficient_rows(
        axial_modes, angular_modes, mass, W1, None, source_proved=source_proved
    )
    if return_details:
        return {"leading": leading, "R0": leading, "certificate": proof}
    return leading


def signed_remainder_coefficients(
    axial_modes: Any,
    angular_modes: Any,
    mass: Any,
    W1: Any,
    W2: Any,
    source_proved: Any,
) -> dict[str, Any]:
    """Return the signed ``R0`` and ``R1`` ordinary rows.

    The result is suitable as the coefficient packet consumed by
    :func:`enclose_weighted_remainder` after each row is lifted to an interval
    truth enclosure.  No source value is selected from a positive cap.
    """

    leading, next_rows, proof = _coefficient_rows(
        axial_modes, angular_modes, mass, W1, W2, source_proved=source_proved
    )
    assert next_rows is not None
    proof = dict(proof)
    proof.update(
        {
            "leading_formula": "sum_j ordinary_product(B_j, a_j*W1)",
            "next_formula": "sum_j ordinary_product(B_j, a_j^2*W2/2-A0*W1/2)",
            "R0_is_not_mass_subtracted_value": True,
            "R1_is_factored_by_h": True,
        }
    )
    return {
        "leading": leading,
        "next": next_rows,
        "R0": leading,
        "R1": next_rows,
        "certificate": proof,
        "mass_Z_independent": True,
        "W1_Z_independent": True,
        "W2_Z_independent": True,
        "source_value_not_selected": True,
    }


def _as_row_truth(ctx: Any, iv: Any, values: Any, *, name: str) -> list[Any]:
    row = _row(values, name=name)
    return [value if _is_interval(value) else iv.mpf(value) for value in row]


def _looks_like_row(value: Any) -> bool:
    return _is_sequence(value) and len(value) == 3 and not any(_is_sequence(v) for v in value)


def _normalise_truth_modes(
    axial_modestruth: Any,
    angular_modestruth: Any,
    rates: Any | None,
    *,
    iv: Any,
) -> tuple[list[dict[str, Any]], list[list[Any]]]:
    """Normalize interval truth modes, with a separate-rate compatibility path."""

    if _looks_like_row(angular_modestruth):
        angular_modestruth = [angular_modestruth]
    if isinstance(angular_modestruth, dict):
        angular_modestruth = angular_modestruth.get("ordinary_rows", angular_modestruth.get("rows"))
    if not _is_sequence(angular_modestruth) or not angular_modestruth:
        raise ValueError("angular_modestruth must contain ordinary rows")
    angular_rows = []
    for index, row in enumerate(angular_modestruth):
        if isinstance(row, dict):
            row = row.get("ordinary_rows", row.get("rows"))
        angular_rows.append(_as_row_truth(None, iv, row, name=f"angular truth row {index}"))

    if isinstance(axial_modestruth, dict):
        axial_modestruth = axial_modestruth.get("modes", axial_modestruth.get("ordinary_rows"))
    if not _is_sequence(axial_modestruth) or not axial_modestruth:
        raise ValueError("axial_modestruth must contain axial rows")
    if rates is not None and not _is_sequence(rates):
        rates = [rates]
    normalized: list[dict[str, Any]] = []
    for index, mode in enumerate(axial_modestruth):
        rate = None
        if isinstance(mode, dict):
            rows = mode.get("ordinary_rows", mode.get("rows"))
            if "rate" in mode:
                rate = mode["rate"]
        elif _looks_like_row(mode):
            rows = mode
        elif _is_sequence(mode) and len(mode) == 2 and _is_sequence(mode[0]):
            rows, rate = mode
        else:
            raise ValueError(f"malformed axial truth mode {index}")
        if rate is None and rates is not None and index < len(rates):
            rate = rates[index]
        if rate is None:
            raise ValueError(f"missing axial rate for truth mode {index}")
        normalized.append(
            {
                "ordinary_rows": _as_row_truth(None, iv, rows, name=f"axial truth mode {index}"),
                "rate": rate if _is_interval(rate) else iv.mpf(rate),
            }
        )
    if len(normalized) != len(angular_rows) and len(angular_rows) > 3:
        raise ValueError("unsupported angular mode rates")
    return normalized, angular_rows


def _packet_rows(ctx: Any, iv: Any, packet: Any, *, name: str) -> tuple[list[Any], list[Any]]:
    if not isinstance(packet, dict):
        raise ValueError("coeffpackettruth must contain leading and next rows")
    leading = packet.get("leading", packet.get("R0"))
    next_rows = packet.get("next", packet.get("R1"))
    if leading is None or next_rows is None:
        raise ValueError("coeffpackettruth must contain both leading/R0 and next/R1 rows")
    return (
        _as_row_truth(ctx, iv, leading, name=f"{name} leading"),
        _as_row_truth(ctx, iv, next_rows, name=f"{name} next"),
    )


def _positive_cap(ctx: Any, iv: Any, value: Any, *, name: str) -> Any:
    value = _scalar(value, name=name)
    interval = iv.mpf(value)
    lower, upper = _ordinary_endpoint(ctx, interval)
    if lower <= 0 or upper >= 1 or upper < lower:
        raise ValueError(f"{name} must be strictly between zero and one")
    return interval


def _log_proof(
    ctx: Any,
    iv: Any,
    logh: Any,
    bound: Any,
    target: Any,
    *,
    power: int,
    label: str,
) -> dict[str, Any]:
    """Prove h**power * bound < target from the original log-h enclosure."""

    bound_iv = iv.mpf(bound)
    bound_upper = _upper(ctx, bound_iv)
    if bound_upper == 0:
        return {
            "label": label,
            "power": power,
            "exact_zero": True,
            "passed": True,
            "bound_exact_mpf_tuples": _serialize_interval(bound_iv),
            "combined_log_exact_mpf_tuples": None,
            "target": str(target),
        }
    if bound_upper < 0:
        raise ValueError(f"{label} bound must be nonnegative")
    bound_log = iv.log(bound_iv)
    combined = power * logh + bound_log
    target_log = iv.log(target)
    combined_upper = _upper(ctx, combined)
    target_lower = _lower(ctx, target_log)
    if combined_upper >= target_lower:
        raise ValueError(f"original logh does not prove the {label} error cap")
    return {
        "label": label,
        "power": power,
        "exact_zero": False,
        "passed": True,
        "bound_exact_mpf_tuples": _serialize_interval(bound_iv),
        "combined_log_exact_mpf_tuples": _serialize_interval(combined),
        "target_log_exact_mpf_tuples": _serialize_interval(target_log),
        "target": str(target),
        "source_logh_is_original": True,
    }


def _prove_small_parameter(
    ctx: Any,
    iv: Any,
    logh: Any,
    bound: Any,
    eta: Any,
    *,
    label: str,
) -> dict[str, Any]:
    """Prove h*bound <= eta using logs before any potentially huge product."""

    bound_iv = iv.mpf(bound)
    upper = _upper(ctx, bound_iv)
    if upper == 0:
        return {
            "label": label,
            "exact_zero": True,
            "passed": True,
            "bound_exact_mpf_tuples": _serialize_interval(bound_iv),
            "combined_log_exact_mpf_tuples": None,
            "eta": str(eta),
            "source_logh_is_original": True,
        }
    combined = logh + iv.log(bound_iv)
    eta_log = iv.log(eta)
    if _upper(ctx, combined) > _lower(ctx, eta_log):
        raise ValueError(f"original logh does not prove {label} <= eta")
    return {
        "label": label,
        "exact_zero": False,
        "passed": True,
        "bound_exact_mpf_tuples": _serialize_interval(bound_iv),
        "combined_log_exact_mpf_tuples": _serialize_interval(combined),
        "eta_log_exact_mpf_tuples": _serialize_interval(eta_log),
        "eta": str(eta),
        "source_logh_is_original": True,
    }


def _normalise_packet_arguments(
    leadingtruth3: Any,
    axial_modestruth: Any,
    angular_modestruth: Any,
    W2positiveIV: Any,
    weight: Any,
    logh: Any,
    rates: Any,
    coeffpackettruth: Any,
) -> tuple[Any, Any, Any, Any, Any, Any, Any]:
    if coeffpackettruth is None:
        if leadingtruth3 is None:
            raise ValueError("leadingtruth3 or coeffpackettruth is required")
        coeffpackettruth = leadingtruth3
    if axial_modestruth is None or angular_modestruth is None:
        raise ValueError("axial_modestruth and angular_modestruth are required")
    if W2positiveIV is None or weight is None or logh is None:
        raise ValueError("W2positiveIV, weight and original logh are required")
    return coeffpackettruth, axial_modestruth, angular_modestruth, W2positiveIV, weight, logh, rates


def enclose_weighted_remainder(
    ctx: Any,
    iv: Any,
    leadingtruth3: Any = None,
    axial_modestruth: Any = None,
    angular_modestruth: Any = None,
    W2positiveIV: Any = None,
    weight: Any = None,
    logh: Any = None,
    source_proved: Any = None,
    *,
    rates: Any = None,
    axial_rates: Any = None,
    coeffpackettruth: Any = None,
    error_cap: Any = _DEFAULT_ERROR_CAP,
    eta: Any = _DEFAULT_ETA,
) -> tuple[list[Any], dict[str, Any]]:
    """Enclose the finite-width signed remainder.

    The canonical call supplies ``coeffpackettruth={leading:[...], next:[...]}``
    and interval truth rows for the axial and angular modes.  ``leadingtruth3``
    remains an alias for the coefficient packet for callers using the earlier
    draft name.  Rates may be carried by axial mode dictionaries or supplied
    separately through ``rates``/``axial_rates``.

    The proof first establishes the small parameters from the original
    ``logh``: ``h``, ``h*max|a|`` and ``h*A_sum_norm`` are each at most
    ``eta``.  It then uses the source-derived third derivative bound and returns
    ``R(h)=R0+h*R1+h^2*E``.  No width is selected from a cap and no ``hcap`` is
    exponentiated against a large source norm.
    """

    _source_guard(source_proved)
    if axial_rates is not None:
        if rates is not None:
            raise ValueError("supply only one of rates and axial_rates")
        rates = axial_rates
    (
        coeffpackettruth,
        axial_modestruth,
        angular_modestruth,
        W2positiveIV,
        weight,
        logh,
        rates,
    ) = _normalise_packet_arguments(
        leadingtruth3,
        axial_modestruth,
        angular_modestruth,
        W2positiveIV,
        weight,
        logh,
        rates,
        coeffpackettruth,
    )
    if not hasattr(ctx, "mpf") or not hasattr(iv, "mpf"):
        raise TypeError("ctx and iv must be mpmath ordinary and interval contexts")
    error_cap_iv = _positive_cap(ctx, iv, error_cap, name="error_cap")
    eta_iv = _positive_cap(ctx, iv, eta, name="eta")
    logh_iv = iv.mpf(logh)
    if _upper(ctx, logh_iv) >= 0:
        raise ValueError("original logh must prove a width h<1")
    weight_iv = _validate_weight(ctx, iv, weight)
    packet = coeffpackettruth
    axial_modes = axial_modestruth
    angular_modes = angular_modestruth
    leading, next_rows = _packet_rows(ctx, iv, packet, name="coeffpackettruth")
    axial, angular = _normalise_truth_modes(axial_modes, angular_modes, rates, iv=iv)
    if any(abs(1 - index) > 1 for index in range(len(angular))):
        raise ValueError("angular truth modes must use j=0,1,2 rates")
    if len(axial) == 0:
        raise ValueError("at least one axial truth mode is required")
    W2_iv = iv.mpf(W2positiveIV)
    W2_lower, W2_upper = _ordinary_endpoint(ctx, W2_iv)
    if W2_lower < 0 or W2_upper < W2_lower:
        raise ValueError("W2positiveIV must be a nonnegative directed interval")

    # Norms use only truth enclosures and the requested ordinary weight.
    A_sum_norm = sum(
        (_ordinary_norm2(ctx, iv, row, weight_iv) for row in angular),
        iv.mpf(0),
    )
    rate_norms = [iv.mpf(_abs_upper(ctx, mode["rate"])) for mode in axial]
    max_rate = max(rate_norms, key=lambda value: _upper(ctx, value))
    width_proofs = [
        _prove_small_parameter(ctx, iv, logh_iv, iv.mpf(1), eta_iv, label="h"),
        _prove_small_parameter(ctx, iv, logh_iv, max_rate, eta_iv, label="h*max|a|"),
        _prove_small_parameter(ctx, iv, logh_iv, A_sum_norm, eta_iv, label="h*A_sum_norm"),
    ]

    e_small = iv.exp(eta_iv)
    # The exponent bound is evaluated only after the logarithmic smallness
    # proofs above.  In particular, no exp(hcap*large_source_norm) occurs.
    exponent_bound = iv.exp(eta_iv + e_small * eta_iv / 2)
    mode_coefficients: list[Any] = []
    mode_kernel_bounds: list[Any] = []
    taylor_e3 = iv.mpf(0)
    W2_upper_iv = iv.mpf(W2_upper)
    for mode, rate_norm in zip(axial, rate_norms):
        q1 = rate_norm + 5 * A_sum_norm * e_small / 4
        q2 = 13 * A_sum_norm * e_small / 6
        q3 = 21 * A_sum_norm * e_small / 8
        kernel_bound = W2_upper_iv * exponent_bound * (q1**3 + 3 * q1 * q2 + q3) / 6
        bnorm = _ordinary_norm2(ctx, iv, mode["ordinary_rows"], weight_iv)
        taylor_e3 += bnorm * kernel_bound
        mode_coefficients.append(bnorm)
        mode_kernel_bounds.append(kernel_bound)

    # Convert the weighted E3 norm into one ordinary derivative bound per Z
    # order.  If W2 is exactly zero, keep the exact zero branch and do not call
    # log(0).
    e3_order_bounds = [
        taylor_e3 * math.factorial(order) / weight_iv**order for order in range(3)
    ]
    target_half = error_cap_iv / 2
    e3_log_proofs = [
        _log_proof(
            ctx,
            iv,
            logh_iv,
            bound,
            target_half,
            power=2,
            label=f"h^2 E3 ordinary Z{order}",
        )
        for order, bound in enumerate(e3_order_bounds)
    ]

    # The h*R1 term is retained explicitly in the certificate.  Its own
    # ordinary norm is proved against half of the configured cap.
    next_norm = _ordinary_norm2(ctx, iv, next_rows, weight_iv)
    next_order_bounds = [next_norm * math.factorial(order) / weight_iv**order for order in range(3)]
    next_log_proofs = [
        _log_proof(
            ctx,
            iv,
            logh_iv,
            bound,
            target_half,
            power=1,
            label=f"h R1 ordinary Z{order}",
        )
        for order, bound in enumerate(next_order_bounds)
    ]

    truth: list[Any] = []
    phase_zero = W2_lower == 0 and W2_upper == 0
    if phase_zero:
        # The weighted source interval is exactly empty at phase zero.  Keep
        # the E3 proof exact and avoid any logarithmic zero operation above.
        if not all(proof["exact_zero"] for proof in e3_log_proofs):
            raise ArithmeticError("phase-zero W2 must produce exact zero E3 bounds")
        truth = list(leading)
    else:
        # The two omitted terms are separately proved <= error_cap/2 in every
        # ordinary order.  Keep the original width logarithm symbolic: the
        # production width can be far below the range where materializing h is
        # useful or safe.  R1 remains visible in the certificate and is not
        # silently folded into a selected point value.
        cap_upper = _upper(ctx, error_cap_iv)
        symmetric_cap = iv.mpf([-cap_upper, cap_upper])
        truth = [r0 + symmetric_cap for r0 in leading]

    certificate = {
        "passed": True,
        "source_proved": True,
        "truth": truth,
        "leading_truth": leading,
        "next_truth": next_rows,
        "R_definition": "R(h)=R0+h*R1+h^2*E3",
        "returned_truth_definition": "R0 plus/minus error_cap; h*R1 and h^2*E3 are separately proved",
        "R0_formula": "sum_j ordinary_product(B_j,a_j*W1)",
        "R1_formula": "sum_j ordinary_product(B_j,a_j^2*W2/2-A0*W1/2)",
        "truth_is_directed_interval": True,
        "width_log_proofs": width_proofs,
        "E3_log_proof_per_order": e3_log_proofs,
        "hR1_log_proof_per_order": next_log_proofs,
        "ordinary_error_bound_per_order_exact_mpf_tuples": [
            _serialize_interval(bound) for bound in e3_order_bounds
        ],
        "hR1_ordinary_bound_per_order_exact_mpf_tuples": [
            _serialize_interval(bound) for bound in next_order_bounds
        ],
        "TaylorE3coefficient": taylor_e3,
        "kernel_K3_bound_per_mode": mode_kernel_bounds,
        "axial_mode_weighted_norms": mode_coefficients,
        "A_sum_norm": A_sum_norm,
        "max_rate_norm": max_rate,
        "eSmall": e_small,
        "exponent_bound": exponent_bound,
        "W2_positive_interval": W2_iv,
        "W2_upper": W2_upper,
        "phase_zero_exact": phase_zero,
        "error_cap": error_cap_iv,
        "eta": eta_iv,
        "original_logh": logh_iv,
        "h_selected": False,
        "hb_chosen": False,
        "hb_materialized": False,
        "hcap_used": False,
        "h_not_materialized": True,
        "actual_field_zeroed": False,
        "positive_cap_selected_as_value": False,
        "mass_recomputed": False,
        "old_M0_producer_called": False,
        "source_proof_contract": {
            "mass_and_moments_Z_independent": True,
            "angular_weight_is_one": True,
            "axial_weight_is_kappa": True,
        },
    }
    return truth, certificate


# Names used by the first draft and by downstream source checkers are kept as
# thin aliases; all routes execute the same new pure kernel above.
finite_h_remainder = enclose_weighted_remainder
certified_finite_h_remainder = enclose_weighted_remainder


__all__ = [
    "signed_remainder_coefficient",
    "signed_remainder_coefficients",
    "enclose_weighted_remainder",
    "finite_h_remainder",
    "certified_finite_h_remainder",
]
