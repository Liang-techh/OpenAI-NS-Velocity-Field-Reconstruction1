"""Directed affine-cell coefficient enclosure for the flatten pressure stage.

On the stored flatten interval write ``y = y_v + L u`` with ``0 <= u <= 1``.
The stored schedule has

``ell(y) = ell(y_v) + (-1/2 - mu) L u``

and ``beta(u) = 2(1-sigma(u))``.  For the positive coefficient magnitudes,

``b_k(beta) = 2**(beta-3) (beta)_k/k!``.

Both the affine exponential and ``b_k(beta(u))`` are decreasing.  On each
``u`` cell this module integrates the affine exponential exactly with
``expm1`` and multiplies it by directed endpoint bounds for ``b_k``.  This
removes the large radial origin from the enclosure and is substantially
tighter than a full endpoint Riemann bound on the density.

The adapter route uses the existing ``ScheduleEndpointEnclosures`` object for
directed stored inputs and requires the exact 100-unit flatten length.  The
returned error is only the stored finite-Gauss flatten-stage error; source
parameter derivations, atom-generation roundoff, core propagation, and five
defects remain outside its scope.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any, Callable, Iterable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_schedule_endpoint_enclosures import (  # noqa: E402
    ScheduleEndpointEnclosures,
)
from lei_ren_part1_paper_uniform_flatten_error import (  # noqa: E402
    stored_flatten_coefficients_directed,
    uniform_error_from_coefficients_directed,
)


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _interval_bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    if hasattr(value, "_mpi_"):
        return tuple(mp.make_mpf(item) for item in value._mpi_)
    value = _mp(value)
    return value, value


def _order(value: Any) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 0:
        raise ValueError("coefficient_order must be a nonnegative integer")
    return int(value)


def _panels(value: Any) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError("panels must be a positive integer")
    return int(value)


def beta_factor(beta: Any, index: Any, context: MPIntervalContext) -> Any:
    """Directed ``2**(beta-3) (beta)_index/index!`` for beta in [0,2]."""

    index = _order(index)
    beta = beta if hasattr(beta, "_mpi_") else context.mpf(str(beta))
    blo, bhi = _interval_bounds(beta)
    if blo < 0 or bhi > 2:
        raise ValueError("flatten beta must lie in [0,2]")
    result = context.exp((beta - 3) * context.ln(2))
    for factor_index in range(1, index + 1):
        result *= (beta + factor_index - 1) / factor_index
    return result


def affine_exponential_integral(
    base: Any,
    slope: Any,
    length: Any,
    left: Any,
    right: Any,
    context: MPIntervalContext,
) -> Any:
    """Integrate ``exp(2*base + 2*slope*length*u)`` over one u cell."""

    base = base if hasattr(base, "_mpi_") else context.mpf(str(base))
    slope = slope if hasattr(slope, "_mpi_") else context.mpf(str(slope))
    length = length if hasattr(length, "_mpi_") else context.mpf(str(length))
    left = left if hasattr(left, "_mpi_") else context.mpf(str(left))
    right = right if hasattr(right, "_mpi_") else context.mpf(str(right))
    width = right - left
    rate = 2 * slope * length
    if _interval_bounds(rate)[1] >= 0:
        raise ValueError("flatten affine slope must be strictly negative")
    return (
        context.exp(2 * base + rate * left)
        * context.expm1(rate * width)
        / rate
    )


def affine_flatten_coefficient_intervals(
    *,
    base: Any,
    length: Any,
    slope: Any,
    beta_source: Callable[[Any], Any],
    coefficient_order: Any = 96,
    panels: Any = 512,
    context: MPIntervalContext | None = None,
) -> dict[str, Any]:
    """Enclose true flatten coefficients using affine exponential cells.

    ``beta_source(u)`` must return a directed interval for beta at the cell
    endpoint.  The analytic input condition is decreasing beta on [0,1].
    ``base``, ``length`` and ``slope`` may be directed intervals.
    """

    if context is None:
        context = MPIntervalContext()
        context.dps = 100
    order = _order(coefficient_order)
    panels = _panels(panels)
    base = base if hasattr(base, "_mpi_") else context.mpf(str(base))
    length = length if hasattr(length, "_mpi_") else context.mpf(str(length))
    slope = slope if hasattr(slope, "_mpi_") else context.mpf(str(slope))
    if _interval_bounds(length)[0] <= 0:
        raise ValueError("flatten length must be positive")
    if _interval_bounds(slope)[1] >= 0:
        raise ValueError("flatten affine slope must be strictly negative")

    step = context.mpf(1) / panels
    points = [context.mpf(i) / panels for i in range(panels + 1)]
    beta_values = []
    for point in points:
        beta = beta_source(point)
        beta = beta if hasattr(beta, "_mpi_") else context.mpf(str(beta))
        beta_lo, beta_hi = _interval_bounds(beta)
        if beta_lo < 0 or beta_hi > 2:
            raise ValueError("beta source interval leaves [0,2]")
        beta_values.append(beta)

    # The affine exponential cell mass is independent of k.  Likewise the
    # rising-factorial recurrence at each beta endpoint is shared by every
    # retained coefficient.  Precompute both tables before accumulating the
    # coefficient intervals; this keeps the actual K=96, 512-cell route
    # practical while retaining directed arithmetic.
    exponential_masses = [
        length
        * affine_exponential_integral(
            base, slope, length, points[panel], points[panel + 1], context
        )
        for panel in range(panels)
    ]
    factor_rows = []
    for beta in beta_values:
        row = [beta_factor(beta, 0, context)]
        for index in range(1, order + 1):
            row.append(row[-1] * (beta + index - 1) / index)
        factor_rows.append(row)

    coefficients = []
    widths = []
    for index in range(order + 1):
        lower = context.mpf(0)
        upper = context.mpf(0)
        for panel in range(panels):
            exponential_mass = exponential_masses[panel]
            factor_left = factor_rows[panel][index]
            factor_right = factor_rows[panel + 1][index]
            factor_lower = _interval_bounds(factor_right)[0]
            factor_upper = _interval_bounds(factor_left)[1]
            lower += exponential_mass * factor_lower
            upper += exponential_mass * factor_upper
        lower_endpoint = _interval_bounds(lower)[0]
        upper_endpoint = _interval_bounds(upper)[1]
        if index % 2:
            coefficient = context.mpf([-upper_endpoint, -lower_endpoint])
        else:
            coefficient = context.mpf([lower_endpoint, upper_endpoint])
        coefficients.append(coefficient)
        lo, hi = _interval_bounds(coefficient)
        widths.append(hi - lo)

    return {
        "coefficients": coefficients,
        "mass_interval": coefficients[0],
        "coefficient_order": order,
        "panels": panels,
        "coefficient_interval_widths": widths,
        "affine_exponential_integrated_exactly": True,
        "beta_endpoint_bounds_directed": True,
        "endpoint_beta_zero_convention": True,
        "directed_interval_arithmetic": True,
        "monotone_beta_assumed": True,
    }


def _adapter_atoms(adapter: Any, quadrature_order: int) -> list[dict[str, Any]]:
    data = adapter._compute_components(int(quadrature_order))
    for component in data["post_components"]:
        if component.get("stage") == "z_flatten":
            return list(component["atoms"])
    raise ValueError("continuous preheat adapter has no z_flatten atom component")


def enclose_affine_flatten_error(
    adapter: Any,
    *,
    quadrature_order: int | None = None,
    coefficient_order: int = 96,
    panels: int = 512,
    z_radius: Any = ".8",
) -> dict[str, Any]:
    """Enclose the stored adapter flatten error via normalized affine cells."""

    if quadrature_order is None:
        quadrature_order = int(adapter.quadrature_order)
    if isinstance(quadrature_order, bool) or int(quadrature_order) < 16:
        raise ValueError("quadrature_order must be at least 16")
    quadrature_order = int(quadrature_order)
    precision = max(80, int(getattr(adapter, "precision", 80)))
    with mp.workdps(precision):
        enclosures = ScheduleEndpointEnclosures(adapter.schedule, precision=precision)
        schedule = adapter.schedule
        raw_left, raw_right = schedule._stage_bounds["z_flatten"]
        if raw_right - raw_left != 100:
            raise ValueError("stored z_flatten length must equal 100")
        enclosures.stage_pressure_upper("z_flatten")
        iv = enclosures.iv
        left_iv = enclosures.scalar(raw_left)
        right_iv = enclosures.scalar(raw_right)
        if (_interval_bounds(left_iv)[0] < 1
                or _interval_bounds(left_iv - enclosures.scalar(schedule.y_d))[0] < 1
                or _interval_bounds(enclosures.scalar(schedule.y_rel) - right_iv)[0] < 0):
            raise ValueError("Flatten affine slope switches are not saturated on the stored interval")
        length = enclosures.scalar(raw_right - raw_left)
        slope = -enclosures.iv.mpf(".5") - enclosures.scalar(schedule.mu)
        base = enclosures.log_amplitude_ratio(raw_left)["interval"]

        def beta_source(u: Any) -> Any:
            u = u if hasattr(u, "_mpi_") else enclosures.iv.mpf(str(u))
            sigma = enclosures.sigma_interval(u)
            return 2 * (1 - sigma)

        true_data = affine_flatten_coefficient_intervals(
            base=base,
            length=length,
            slope=slope,
            beta_source=beta_source,
            coefficient_order=coefficient_order,
            panels=panels,
            context=enclosures.iv,
        )
        atoms = _adapter_atoms(adapter, quadrature_order)
        stored_coefficients, stored_mass = stored_flatten_coefficients_directed(
            atoms, coefficient_order, enclosures.iv
        )
        result = uniform_error_from_coefficients_directed(
            true_data["coefficients"],
            stored_coefficients,
            stored_mass,
            z_radius=z_radius,
            context=enclosures.iv,
        )
        result.update(
            {
                "stage": "z_flatten",
                "stored_quadrature_order": quadrature_order,
                "stored_atom_count": len(atoms),
                "stored_flatten_length": length,
                "flatten_slope": slope,
                "true_coefficient_panels": true_data["panels"],
                "true_coefficient_interval_widths": true_data[
                    "coefficient_interval_widths"
                ],
                "true_mass_interval": true_data["mass_interval"],
                "affine_exponential_integrated_exactly": True,
                "beta_endpoint_bounds_directed": True,
                "endpoint_beta_zero_convention": True,
                "directed_interval_arithmetic": True,
                "stored_parameters_only": True,
                "original_parameter_errors_enclosed": False,
                "stored_gauss_atom_generation_roundoff_unenclosed": True,
                "quadrature_roundoff_enclosed": False,
                "core_error_enclosed": False,
                "five_defect_interval_closure": False,
            }
        )
        return result


uniform_affine_flatten_error = enclose_affine_flatten_error


__all__ = [
    "beta_factor",
    "affine_exponential_integral",
    "affine_flatten_coefficient_intervals",
    "enclose_affine_flatten_error",
    "uniform_affine_flatten_error",
]
