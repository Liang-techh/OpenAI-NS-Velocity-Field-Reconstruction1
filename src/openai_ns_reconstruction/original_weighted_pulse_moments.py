"""Independent weighted moments of the original flat pulse.

For the original first-switch cutoff

    sigma(t) = gp(t) / (gp(t) + gp(1 - t)),
    gp(t) = exp(-1/t**2) for t > 0,

write ``kappa = 1 - sigma`` and evaluate the genuinely new moments

    W_k(s) = integral_0^s t**k kappa(t) dt,  k >= 1.

The displayed approximation of each row is produced by a fresh scalar
``mpmath`` quadrature.  Its directed source enclosure is a separate sum of
positive endpoint rectangles.  On a cell ``[a,b]`` the rectangle uses the
exact weighted length ``integral_a^b t**k dt`` and the endpoint interval
``[kappa(b), kappa(a)]``.  A rectangle midpoint is never used as a source
value.

This module deliberately has no dependency on the accepted M 0 mass rows or
on an old pulse producer.  The endpoint ``W_k(0)=0`` is an exact branch; no
closed value is asserted for ``W_1(1)`` or ``W_2(1)``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import Any

from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import encode_pair


_SOURCE_FORMULA = (
    "gp(t)=exp(-1/t^2) for t>0; "
    "sigma(t)=gp(t)/(gp(t)+gp(1-t)); kappa(t)=1-sigma(t)"
)
_WEIGHTED_FORMULA = "W_k(s)=integral_0^s t^k*(1-sigma(t)) dt, k>=1"


def _number(ctx: Any, value: Any, digits: int = 90) -> str:
    """Serialize an ordinary mpmath value without a binary float hop."""

    value = ctx.mpf(value)
    return "0" if value == 0 else ctx.nstr(value, digits)


def _interval_bounds(ctx: Any, value: Any) -> tuple[Any, Any]:
    """Return outward ordinary-context endpoints of an interval value."""

    try:
        lower, upper = value._mpi_
    except AttributeError as exc:  # pragma: no cover - defensive API guard
        raise TypeError("expected an mpmath interval") from exc
    return ctx.mpf(lower), ctx.mpf(upper)


def _exact_interval_from_bounds(ctx: Any, iv: MPIntervalContext, lower: Any, upper: Any) -> Any:
    """Build an interval from exact directed mpf tuples."""

    lower_tuple = tuple(lower._mpi_[0]) if hasattr(lower, "_mpi_") else tuple(lower._mpf_)
    upper_tuple = tuple(upper._mpi_[1]) if hasattr(upper, "_mpi_") else tuple(upper._mpf_)
    return iv.mpf([ctx.make_mpf(lower_tuple), ctx.make_mpf(upper_tuple)])


def _serialized_error_number(ctx: Any, encoded: dict[str, Any]) -> str:
    error_tuple = encoded["serialized_approximation_absolute_error_upper_exact_mpf_tuple"]
    return _number(ctx, ctx.make_mpf(tuple(error_tuple)))


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _phase_label(ctx: Any, value: Any) -> tuple[str, Any]:
    """Preserve string phase labels while validating their numeric value."""

    label = value if isinstance(value, str) else _number(ctx, value)
    point = ctx.mpf(value)
    if point < 0 or point > 1:
        raise ValueError("phases must lie in [0, 1]")
    return label, point


def _gp_scalar(ctx: Any, value: Any) -> Any:
    point = ctx.mpf(value)
    if point <= 0:
        return ctx.mpf(0)
    return ctx.exp(-1 / point**2)


def _sigma_scalar(ctx: Any, value: Any) -> Any:
    """Evaluate the original gp ratio in an ordinary scalar context."""

    point = ctx.mpf(value)
    if point <= 0:
        return ctx.mpf(0)
    if point >= 1:
        return ctx.mpf(1)
    left = _gp_scalar(ctx, point)
    right = _gp_scalar(ctx, 1 - point)
    return left / (left + right)


def _gp_interval(iv: MPIntervalContext, ctx: Any, value: Any) -> Any:
    lo, hi = _interval_bounds(ctx, value)
    if hi <= 0:
        return iv.mpf(0)
    if lo <= 0:
        raise ValueError("a source node may not straddle t=0")
    return iv.exp(-1 / value**2)


def _sigma_interval(iv: MPIntervalContext, ctx: Any, value: Any) -> Any:
    """Evaluate sigma at an endpoint interval, retaining flat branches."""

    lo, hi = _interval_bounds(ctx, value)
    if hi <= 0:
        return iv.mpf(0)
    if lo >= 1:
        return iv.mpf(1)
    left = _gp_interval(iv, ctx, value)
    right = _gp_interval(iv, ctx, iv.mpf(1) - value)
    raw = left / (left + right)
    # The independent numerator/denominator intervals can straddle 0 or 1
    # by a few ulps near a flat endpoint.  The source identity itself gives
    # the exact range 0 <= sigma <= 1, so intersect that proved range before
    # forming kappa.  This is an enclosure tightening, never a selected value.
    raw_lower, raw_upper = _interval_bounds(ctx, raw)
    return iv.mpf([max(ctx.mpf(0), raw_lower), min(ctx.mpf(1), raw_upper)])


def _kappa_interval(iv: MPIntervalContext, ctx: Any, value: Any) -> Any:
    lo, hi = _interval_bounds(ctx, value)
    if hi <= 0:
        return iv.mpf(1)
    if lo >= 1:
        return iv.mpf(0)
    # Evaluate the reflected gp ratio directly.  It is algebraically
    # 1-sigma, but avoids subtracting a rounded sigma from 1 in the flat
    # right tail, where that subtraction could erase a positive source.
    left = _gp_interval(iv, ctx, value)
    right = _gp_interval(iv, ctx, iv.mpf(1) - value)
    raw = right / (left + right)
    raw_lower, raw_upper = _interval_bounds(ctx, raw)
    return iv.mpf([max(ctx.mpf(0), raw_lower), min(ctx.mpf(1), raw_upper)])


def _as_interval(iv: MPIntervalContext, value: Any) -> Any:
    return value if hasattr(value, "_mpi_") else iv.mpf(value)


def _weighted_length(iv: MPIntervalContext, left: Any, right: Any, power: int) -> Any:
    """Outward interval form of the exact integral of t**power on [left,right]."""

    left_iv = iv.mpf(left)
    right_iv = iv.mpf(right)
    return (right_iv ** (power + 1) - left_iv ** (power + 1)) / (power + 1)


def positive_weighted_rectangle_bounds(
    ctx: Any,
    iv: MPIntervalContext,
    phase: Any,
    power: int,
    cells: int,
    kappa: Callable[[Any], Any],
    *,
    return_details: bool = False,
) -> Any:
    """Enclose a positive weighted integral using monotone endpoint rectangles.

    ``kappa`` receives an ordinary-context endpoint and returns either an
    ordinary value or an interval.  It must be nonnegative and nonincreasing
    over the requested interval.  The callback form keeps this helper useful
    for independent analytic fixtures such as ``kappa=1`` and
    ``kappa=1-t`` while the production pulse passes its interval gp-ratio
    evaluator.

    The default return is the directed interval.  With ``return_details`` a
    ``(interval, metadata)`` pair is returned; metadata is JSON serializable.
    """

    power = _positive_integer(power, "power")
    cells = _positive_integer(cells, "cells")
    phase = ctx.mpf(phase)
    if phase < 0 or phase > 1:
        raise ValueError("phase must lie in [0, 1]")

    lower_sum = iv.mpf(0)
    upper_sum = iv.mpf(0)
    if phase != 0:
        for index in range(cells):
            left = phase * index / cells
            right = phase * (index + 1) / cells
            weight = _weighted_length(iv, left, right, power)
            weight_lower, weight_upper = _interval_bounds(ctx, weight)
            if weight_lower < 0 or weight_upper < weight_lower:
                raise ArithmeticError("weighted cell length is not positive")

            left_kappa = _as_interval(iv, kappa(left))
            right_kappa = _as_interval(iv, kappa(right))
            left_kappa_lower, left_kappa_upper = _interval_bounds(ctx, left_kappa)
            right_kappa_lower, right_kappa_upper = _interval_bounds(ctx, right_kappa)
            if (
                left_kappa_lower < 0
                or right_kappa_lower < 0
                or left_kappa_upper < left_kappa_lower
                or right_kappa_upper < right_kappa_lower
            ):
                raise ValueError("kappa must be nonnegative")
            if right_kappa_lower > left_kappa_upper:
                raise ValueError("kappa must be nonincreasing on each cell")

            # The interval products are deliberately kept separate.  The
            # lower endpoint uses kappa(right), and the upper endpoint uses
            # kappa(left); no midpoint is formed or selected.
            lower_cell = weight * right_kappa
            upper_cell = weight * left_kappa
            lower_cell_lower, _ = _interval_bounds(ctx, lower_cell)
            _, upper_cell_upper = _interval_bounds(ctx, upper_cell)
            lower_sum += iv.mpf(lower_cell_lower)
            upper_sum += iv.mpf(upper_cell_upper)

    lower, _ = _interval_bounds(ctx, lower_sum)
    _, upper = _interval_bounds(ctx, upper_sum)
    truth = _exact_interval_from_bounds(ctx, iv, lower, upper)
    if not return_details:
        return truth

    width = upper - lower
    details = {
        "cells": cells,
        "method": "positive endpoint rectangles",
        "cell_weight_formula": "(b^(power+1)-a^(power+1))/(power+1)",
        "lower_endpoint_rule": "cell_weight*kappa(b)",
        "upper_endpoint_rule": "cell_weight*kappa(a)",
        "lower": _number(ctx, lower),
        "upper": _number(ctx, upper),
        "enclosure_width": _number(ctx, width),
        "source_enclosure_exact_mpf_tuples": [list(value) for value in truth._mpi_],
        "source_value_not_selected": True,
        "rectangle_midpoint_used_as_value": False,
    }
    return truth, details


# A short alias makes the generic helper convenient in independent fixtures.
weighted_rectangle_bounds = positive_weighted_rectangle_bounds


def _adaptive_weighted_scalar(ctx: Any, phase: Any, power: int) -> Any:
    endpoint = ctx.mpf(phase)
    if endpoint == 0:
        return ctx.mpf(0)

    def integrand(value: Any) -> Any:
        return value**power * (1 - _sigma_scalar(ctx, value))

    return ctx.quad(integrand, [0, endpoint])


def _refinement_check(ctx: Any, boxes: dict[int, tuple[Any, Any]]) -> dict[str, Any]:
    ordered = sorted(boxes)
    comparisons: list[dict[str, Any]] = []
    for coarse_cells, fine_cells in zip(ordered, ordered[1:]):
        coarse_lower, coarse_upper = boxes[coarse_cells]
        fine_lower, fine_upper = boxes[fine_cells]
        nested_partition = fine_cells % coarse_cells == 0
        nested = bool(
            nested_partition
            and coarse_lower <= fine_lower
            and fine_upper <= coarse_upper
        )
        comparisons.append(
            {
                "coarse_cells": coarse_cells,
                "fine_cells": fine_cells,
                "nested_partition": nested_partition,
                "coarse_lower": _number(ctx, coarse_lower),
                "coarse_upper": _number(ctx, coarse_upper),
                "fine_lower": _number(ctx, fine_lower),
                "fine_upper": _number(ctx, fine_upper),
                "passed": nested,
            }
        )
    return {
        "comparisons": comparisons,
        "passed": bool(comparisons) and all(item["passed"] for item in comparisons),
    }


def _weighted_row(
    ctx: Any,
    iv: MPIntervalContext,
    phase_label: str,
    phase: Any,
    power: int,
    cells: Sequence[int],
) -> dict[str, Any]:
    boxes: dict[str, dict[str, Any]] = {}
    box_values: dict[int, tuple[Any, Any]] = {}

    def source_kappa(value: Any) -> Any:
        return _kappa_interval(iv, ctx, iv.mpf(value))

    for cell_count in cells:
        truth, details = positive_weighted_rectangle_bounds(
            ctx,
            iv,
            phase,
            power,
            cell_count,
            source_kappa,
            return_details=True,
        )
        lower, upper = _interval_bounds(ctx, truth)
        box_values[cell_count] = (lower, upper)
        boxes[str(cell_count)] = details

    scalar_approximation = _adaptive_weighted_scalar(ctx, phase, power)
    fine_cells = max(cells)
    fine_lower, fine_upper = box_values[fine_cells]
    fine_truth = _exact_interval_from_bounds(ctx, iv, fine_lower, fine_upper)
    encoded = encode_pair(ctx, iv, scalar_approximation, fine_truth)
    exact_zero = phase == 0
    if exact_zero:
        encoded["approximation"] = "0"
        encoded["approximation_is_exact_source_value"] = True
        encoded["source_enclosure_is_not_a_selected_value"] = False
    error = _serialized_error_number(ctx, encoded)
    refinement = _refinement_check(ctx, box_values)
    return {
        "phase": phase_label,
        "power": power,
        "kind": "weighted_source_integral",
        "formula": f"W_{power}(s)=integral_0^s t^{power}*(1-sigma(t)) dt",
        "source_formula": _SOURCE_FORMULA,
        "quadrature": {
            "method": "independent scalar adaptive quadrature (ctx.quad)",
            "approximation": encoded["approximation"],
            "source_value_not_selected": not exact_zero,
            "approximation_is_exact_source_value": exact_zero,
            "rectangle_midpoint_used_as_value": False,
        },
        "boxes": boxes,
        "refinement": refinement,
        **encoded,
        "lower": _number(ctx, fine_lower),
        "upper": _number(ctx, fine_upper),
        "enclosure_width": _number(ctx, fine_upper - fine_lower),
        "enclosure_radius": _number(ctx, (fine_upper - fine_lower) / 2),
        "absolute_error_upper": error,
        "serialized_approximation_absolute_error_upper": error,
        "approximation_is_not_exact_source_value": not exact_zero,
        "source_value_not_selected": not exact_zero,
        "source_enclosure_is_not_a_selected_value": not exact_zero,
        "rectangle_midpoint_used_as_value": False,
        "finite_quadrature_nodes_are_machinery": True,
        "accepted_M0_mass_reused": False,
    }


def build_weighted_moments(
    ctx: Any,
    iv: MPIntervalContext,
    phases: Iterable[Any],
    *,
    powers: Sequence[int] = (1, 2),
    cells: Sequence[int] = (128, 256),
) -> dict[str, Any]:
    """Build independent weighted pulse rows for the requested phases.

    The API intentionally accepts caller-owned ordinary and interval mpmath
    contexts so a checker can choose its own precision.  Every row is JSON
    serializable, and the exact directed endpoints emitted by
    :func:`encode_pair` remain available under
    ``source_enclosure_exact_mpf_tuples``.
    """

    if not hasattr(ctx, "mpf") or not hasattr(ctx, "quad"):
        raise TypeError("ctx must be an ordinary mpmath context")
    if not hasattr(iv, "mpf"):
        raise TypeError("iv must be an mpmath interval context")

    powers = tuple(_positive_integer(power, "power") for power in powers)
    if not powers:
        raise ValueError("powers must not be empty")
    if any(power < 1 for power in powers):
        raise ValueError("weighted moments require powers >= 1; M 0 is not part of this API")
    cells = tuple(_positive_integer(count, "cells") for count in cells)
    if not cells:
        raise ValueError("cells must not be empty")
    if len(set(cells)) != len(cells):
        raise ValueError("cells must be distinct")

    phase_items = tuple(_phase_label(ctx, phase) for phase in phases)
    rows: list[dict[str, Any]] = []
    for phase_label, phase in phase_items:
        for power in powers:
            rows.append(_weighted_row(ctx, iv, phase_label, phase, power, cells))

    zero_rows = [row for row in rows if ctx.mpf(row["phase"]) == 0]
    required_refinement = all(row["refinement"]["passed"] for row in rows)
    controls = {
        "powers_are_genuinely_weighted": all(row["power"] >= 1 for row in rows),
        "no_M0_mass_row": all(row["power"] >= 1 for row in rows),
        "endpoint_s0_exact_zero": bool(zero_rows)
        and all(
            row["approximation"] == "0"
            and row["absolute_error_upper"] == "0"
            and all(
                box["lower"] == "0" and box["upper"] == "0"
                for box in row["boxes"].values()
            )
            for row in zero_rows
        ),
        "independent_quadrature_present": all(
            row["quadrature"]["method"] == "independent scalar adaptive quadrature (ctx.quad)"
            for row in rows
        ),
        "separate_positive_boxes_present": all(
            row["boxes"] and row["rectangle_midpoint_used_as_value"] is False for row in rows
        ),
        "independent_refinement_passed": required_refinement,
        "accepted_M0_mass_reused": False,
        "old_mass_producer_called": False,
    }
    required = (
        "powers_are_genuinely_weighted",
        "no_M0_mass_row",
        "endpoint_s0_exact_zero",
        "independent_quadrature_present",
        "separate_positive_boxes_present",
        "independent_refinement_passed",
    )
    controls["required_positive_controls_passed"] = all(controls[name] for name in required)
    return {
        "source": {
            "formula": _SOURCE_FORMULA,
            "weighted_formula": _WEIGHTED_FORMULA,
            "kappa_monotonicity": "kappa(t)=1-sigma(t) is nonnegative and decreasing on [0,1]",
            "endpoint_convention": {"sigma(0)": "0", "sigma(1)": "1", "W_k(0)": "0"},
            "closed_symmetry_values_claimed": False,
            "accepted_M0_mass_reused": False,
            "source_value_not_selected": True,
        },
        "phases": [label for label, _ in phase_items],
        "powers": list(powers),
        "cells": list(cells),
        "rows": rows,
        "controls": controls,
        "limits": {
            "M0_rows": False,
            "W1_or_W2_closed_symmetry_values": False,
            "rectangle_midpoints_selected_as_values": False,
            "global_five_moment_closure": False,
            "bounded_weighted_interface_complete": controls["required_positive_controls_passed"],
            "required_positive_control_names": list(required),
            "required_positive_controls_passed": controls["required_positive_controls_passed"],
        },
        "passed": controls["required_positive_controls_passed"],
        "complete": False,
    }


__all__ = [
    "build_weighted_moments",
    "positive_weighted_rectangle_bounds",
    "weighted_rectangle_bounds",
]
