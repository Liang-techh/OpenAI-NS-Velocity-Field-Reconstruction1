"""Uniform axial error bounds for the stored finite flatten quadrature.

The flatten pressure integral on ``[y_v, y_f]`` has the form

``I(Z) = integral w(y) (1 + Z**2)**(-beta(y)) dy``

with ``0 <= beta <= 2`` and positive ``w``.  The continuous source is
represented downstream by a finite positive Gauss sum.  This module compares
those two representations uniformly on ``abs(Z) <= z_radius``.  It does not
change the stored schedule, pressure datum, core, or any paper parameter.

For ``K >= 0`` we use

``(1 + Z**2)**(-beta) = sum c_k(beta) Z**(2 k)`` and
``abs(c_k(beta)) <= k + 1``.

The true coefficients are enclosed by monotone endpoint Riemann sums on the
flatten interval.  On the stored schedule, ``beta`` decreases from 2 to 0
and the positive coefficient magnitudes decrease: on this stage the stored
log-amplitude has slope ``-1/2-mu``.  The adapter route uses
``ScheduleEndpointEnclosures`` and directed interval arithmetic for endpoint
weights, beta, coefficient differences, and tails.  The generic scalar core
is retained as a nominal resolved diagnostic and is explicitly marked
conditional on its MP endpoint arithmetic.

The omitted series terms are bounded by weighted geometric tails, including
the first and second ``Z`` derivatives.  A positive-mass fallback is retained
because it is often a useful absolute bound when the flatten mass is tiny.
The result is a flatten-stage finite-datum bound only.  It does not enclose
source-parameter derivation, quadrature roundoff, the core, or five moments.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any, Callable, Iterable, Sequence

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_continuous_axial_pulse import (  # noqa: E402
    ContinuousAxialPulse,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import (  # noqa: E402
    ScheduleEndpointEnclosures,
    endpoints,
)


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _validate_order(value: Any, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _validate_panels(value: Any) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError("panels must be a positive integer")
    return int(value)


def _interval(value: Any) -> tuple[mp.mpf, mp.mpf]:
    """Convert a scalar or a two-endpoint value to an ordered MP interval."""

    if hasattr(value, "_mpi_"):
        left, right = tuple(mp.make_mpf(item) for item in value._mpi_)
    elif isinstance(value, (tuple, list)) and len(value) == 2:
        left, right = _mp(value[0]), _mp(value[1])
    else:
        left = right = _mp(value)
    if left > right:
        raise ValueError("interval lower endpoint exceeds upper endpoint")
    return left, right


def _coefficient(beta: Any, index: int) -> mp.mpf:
    """Return ``(-1)^index (beta)_index/index!``."""

    beta = _mp(beta)
    if beta < 0 or beta > 2:
        raise ValueError("flatten beta must lie in [0, 2]")
    index = int(index)
    if index == 0:
        return mp.mpf(1)
    return (-1 if index % 2 else 1) * mp.rf(beta, index) / mp.factorial(index)


def _weighted_geometric_moments(
    first: int, x: Any, context: MPIntervalContext | None = None
) -> tuple[Any, ...]:
    """Return ``sum_{k>=first} k^j x^k`` for ``j=0,1,2,3``.

    The formulas are obtained by shifting ``k=first+n`` and using the four
    standard geometric moments of ``n``.  They avoid a finite sampled tail.
    """

    first = int(first)
    if first < 0 or not (0 <= _interval(x)[0] and _interval(x)[1] < 1):
        raise ValueError("geometric tail requires first>=0 and 0<=x<1")
    zero = context.mpf(0) if context is not None else mp.mpf(0)
    one = context.mpf(1) if context is not None else mp.mpf(1)
    if x == 0:
        return (zero,) * 4 if first else (one,) + (zero,) * 3
    one_minus = one - x
    xm = x**first
    n0 = 1 / one_minus
    n1 = x / one_minus**2
    n2 = x * (1 + x) / one_minus**3
    n3 = x * (1 + 4 * x + x * x) / one_minus**4
    s0 = xm * n0
    s1 = xm * (first * n0 + n1)
    s2 = xm * (first * first * n0 + 2 * first * n1 + n2)
    s3 = xm * (
        first**3 * n0 + 3 * first * first * n1 + 3 * first * n2 + n3
    )
    return s0, s1, s2, s3


def weighted_geometric_tail(
    z_radius: Any,
    coefficient_order: Any,
    *,
    context: MPIntervalContext | None = None,
) -> dict[str, Any]:
    """Bound omitted value/first/second-Z binomial tails uniformly.

    If ``m=coefficient_order+1`` and ``a=abs(z_radius)``, the returned
    factors bound the tails of ``q^-beta``, its first derivative, and its
    second derivative for every ``0<=beta<=2``.  The factors exclude the
    positive radial mass; callers multiply by the true and stored masses.
    """

    if context is None:
        a = abs(_mp(z_radius))
    else:
        lo, hi = _interval(z_radius)
        a = context.mpf(str(max(abs(lo), abs(hi))))
    if not (0 < _interval(a)[0] and _interval(a)[1] < 1):
        raise ValueError("z_radius must satisfy 0 < z_radius < 1")
    order = _validate_order(coefficient_order, "coefficient_order")
    s0, s1, s2, s3 = _weighted_geometric_moments(
        order + 1, a * a, context
    )
    two = context.mpf(2) if context is not None else mp.mpf(2)
    value = s1 + s0
    first = two * (s2 + s1) / a
    second = (4 * s3 + 2 * s2 - 2 * s1) / (a * a)
    return {"value": value, "first_Z": first, "second_Z": second}


def _derivative_weight(
    index: int,
    derivative: int,
    radius: Any,
    context: MPIntervalContext | None = None,
) -> Any:
    zero = context.mpf(0) if context is not None else mp.mpf(0)
    if derivative == 0:
        return radius ** (2 * index)
    if derivative == 1:
        return 2 * index * radius ** (2 * index - 1) if index else zero
    if derivative == 2:
        return (
            2 * index * (2 * index - 1) * radius ** (2 * index - 2)
            if index
            else zero
        )
    raise ValueError("only value, first_Z and second_Z are supported")


def _atom_pair(atom: Any) -> tuple[mp.mpf, mp.mpf]:
    """Read one positive finite Gauss atom as ``(mass_at_Z0, beta)``."""

    if isinstance(atom, dict):
        mass = atom.get("atom_at_Z0", atom.get("mass", atom.get("weight")))
        beta = atom.get("beta")
        if mass is None or beta is None:
            raise ValueError("atom dict needs mass/weight and beta")
    else:
        try:
            mass, beta = atom
        except (TypeError, ValueError) as exc:
            raise ValueError("atom must be a pair or a dict") from exc
    mass, beta = _mp(mass), _mp(beta)
    if mass < 0:
        raise ValueError("stored flatten atoms must be nonnegative")
    if beta < 0 or beta > 2:
        raise ValueError("stored flatten beta must lie in [0, 2]")
    return mass, beta


def stored_flatten_coefficients(
    atoms: Iterable[Any], coefficient_order: Any
) -> tuple[list[mp.mpf], mp.mpf]:
    """Return finite-atom binomial coefficients and their positive mass."""

    order = _validate_order(coefficient_order, "coefficient_order")
    rows = [_atom_pair(atom) for atom in atoms]
    coefficients = [mp.mpf(0)] * (order + 1)
    mass = mp.fsum(row[0] for row in rows)
    for mass_i, beta in rows:
        magnitude = mp.mpf(1)
        sign = mp.mpf(1)
        coefficients[0] += mass_i
        for index in range(1, order + 1):
            magnitude *= (beta + index - 1) / index
            sign = -sign
            coefficients[index] += mass_i * sign * magnitude
    return coefficients, mass


def stored_flatten_coefficients_directed(
    atoms: Iterable[Any],
    coefficient_order: Any,
    context: MPIntervalContext,
) -> tuple[list[Any], Any]:
    """Directed interval version of :func:`stored_flatten_coefficients`."""

    order = _validate_order(coefficient_order, "coefficient_order")
    rows = [_atom_pair(atom) for atom in atoms]
    coefficients = [context.mpf(0)] * (order + 1)
    mass = context.mpf(0)
    for mass_i, beta in rows:
        mass_interval = context.mpf(str(mass_i))
        beta_interval = context.mpf(str(beta))
        mass += mass_interval
        magnitude = context.mpf(1)
        sign = context.mpf(1)
        coefficients[0] += mass_interval
        for index in range(1, order + 1):
            magnitude *= (beta_interval + index - 1) / index
            sign = -sign
            coefficients[index] += mass_interval * sign * magnitude
    return coefficients, mass


def monotone_coefficient_intervals(
    source: Callable[[mp.mpf], tuple[Any, Any]],
    left: Any,
    right: Any,
    *,
    coefficient_order: Any = 96,
    panels: Any = 1024,
    monotone_decreasing: bool = True,
) -> dict[str, Any]:
    """Enclose true flatten coefficients by endpoint sums.

    ``source(y)`` returns ``(w(y), beta(y))``.  The required analytic
    condition is that ``w(y)*abs(c_k(beta(y)))`` is non-increasing for every
    retained ``k``.  This follows on the stored flatten stage from positive
    decreasing amplitude and decreasing beta.  A caller using a different
    source must supply the same proof and set ``monotone_decreasing=True``.
    The sampled monotonicity flags in the receipt are diagnostics; they are
    not the proof itself.
    """

    if not monotone_decreasing:
        raise ValueError(
            "coefficient intervals require an analytic monotone-decrease declaration"
        )
    left, right = _mp(left), _mp(right)
    if right <= left:
        raise ValueError("flatten interval must have positive length")
    order = _validate_order(coefficient_order, "coefficient_order")
    panels = _validate_panels(panels)
    step = (right - left) / panels
    points = [left + step * i for i in range(panels + 1)]
    values: list[tuple[mp.mpf, mp.mpf]] = []
    for point in points:
        weight, beta = source(point)
        weight, beta = _mp(weight), _mp(beta)
        if weight < 0:
            raise ValueError("true flatten weight must be nonnegative")
        if beta < 0 or beta > 2:
            raise ValueError("true flatten beta must lie in [0, 2]")
        values.append((weight, beta))

    # Reuse the rising-factorial recurrence at every endpoint.  This keeps
    # the cost linear in ``panels * coefficient_order`` instead of repeatedly
    # constructing each rising factorial from scratch.
    magnitudes: list[list[mp.mpf]] = []
    for weight, beta in values:
        row = [mp.mpf(1)]
        for index in range(1, order + 1):
            row.append(row[-1] * (beta + index - 1) / index)
        magnitudes.append([weight * item for item in row])

    coefficients: list[tuple[mp.mpf, mp.mpf]] = []
    widths: list[mp.mpf] = []
    monotone_weight = True
    monotone_beta = True
    for index in range(1, len(values)):
        monotone_weight &= values[index][0] <= values[index - 1][0]
        monotone_beta &= values[index][1] <= values[index - 1][1]
    for coefficient_index in range(order + 1):
        lower = mp.mpf(0)
        upper = mp.mpf(0)
        for panel in range(panels):
            magnitude_left = magnitudes[panel][coefficient_index]
            magnitude_right = magnitudes[panel + 1][coefficient_index]
            lower += step * magnitude_right
            upper += step * magnitude_left
        if coefficient_index % 2:
            coefficients.append((-upper, -lower))
        else:
            coefficients.append((lower, upper))
        widths.append(upper - lower)

    return {
        "coefficients": coefficients,
        "mass_interval": coefficients[0],
        "panels": panels,
        "coefficient_order": order,
        "monotone_weight_diagnostic": bool(monotone_weight),
        "monotone_beta_diagnostic": bool(monotone_beta),
        "endpoint_Riemann_enclosure": True,
        "monotone_decrease_assumed": True,
        "directed_interval_arithmetic": False,
        "nominal_mp_roundoff_conditional": True,
        "coefficient_interval_widths": widths,
    }


def monotone_coefficient_intervals_directed(
    source: Callable[[Any], tuple[Any, Any]],
    left: Any,
    right: Any,
    *,
    coefficient_order: Any = 96,
    panels: Any = 1024,
    monotone_decreasing: bool = True,
    context: MPIntervalContext | None = None,
) -> dict[str, Any]:
    """Directed interval endpoint enclosure for true flatten coefficients.

    ``source`` returns interval enclosures for ``(w(y), beta(y))``.  The
    monotone-decrease declaration is analytic; endpoint diagnostics alone do
    not promote a nominal MP computation to a certificate.
    """

    if context is None:
        context = MPIntervalContext()
        context.dps = 80
    if not monotone_decreasing:
        raise ValueError(
            "coefficient intervals require an analytic monotone-decrease declaration"
        )
    left, right = _mp(left), _mp(right)
    if right <= left:
        raise ValueError("flatten interval must have positive length")
    order = _validate_order(coefficient_order, "coefficient_order")
    panels = _validate_panels(panels)
    step = context.mpf(str((right - left) / panels))
    points = [left + (right - left) * i / panels for i in range(panels + 1)]
    values: list[tuple[Any, Any]] = []
    for point in points:
        weight, beta = source(point)
        weight = weight if hasattr(weight, "_mpi_") else context.mpf(str(weight))
        beta = beta if hasattr(beta, "_mpi_") else context.mpf(str(beta))
        wlo, whi = _interval(weight)
        blo, bhi = _interval(beta)
        if wlo < 0 or blo < 0 or bhi > 2:
            raise ValueError("directed source interval leaves the allowed flatten range")
        values.append((weight, beta))

    magnitudes: list[list[Any]] = []
    for weight, beta in values:
        row = [weight]
        magnitude = context.mpf(1)
        for index in range(1, order + 1):
            magnitude *= (beta + index - 1) / index
            row.append(weight * magnitude)
        magnitudes.append(row)

    coefficients: list[Any] = []
    widths: list[mp.mpf] = []
    monotone_weight = True
    monotone_beta = True
    for index in range(1, len(values)):
        monotone_weight &= _interval(values[index][0])[1] <= _interval(values[index - 1][0])[1]
        monotone_beta &= _interval(values[index][1])[1] <= _interval(values[index - 1][1])[1]
    for coefficient_index in range(order + 1):
        lower = context.mpf(0)
        upper = context.mpf(0)
        for panel in range(panels):
            right_lower = _interval(magnitudes[panel + 1][coefficient_index])[0]
            left_upper = _interval(magnitudes[panel][coefficient_index])[1]
            lower += step * right_lower
            upper += step * left_upper
        lower_endpoint = _interval(lower)[0]
        upper_endpoint = _interval(upper)[1]
        coefficient = context.mpf(
            [
                lower_endpoint if coefficient_index % 2 == 0 else -upper_endpoint,
                upper_endpoint if coefficient_index % 2 == 0 else -lower_endpoint,
            ]
        )
        coefficients.append(coefficient)
        lo, hi = _interval(coefficient)
        widths.append(hi - lo)
    return {
        "coefficients": coefficients,
        "mass_interval": coefficients[0],
        "panels": panels,
        "coefficient_order": order,
        "monotone_weight_diagnostic": bool(monotone_weight),
        "monotone_beta_diagnostic": bool(monotone_beta),
        "endpoint_Riemann_enclosure": True,
        "monotone_decrease_assumed": True,
        "directed_interval_arithmetic": True,
        "coefficient_interval_widths": widths,
    }


def uniform_error_from_coefficients(
    true_coefficient_intervals: Sequence[Any],
    stored_coefficients: Sequence[Any],
    stored_mass: Any,
    *,
    z_radius: Any = ".8",
) -> dict[str, Any]:
    """Combine coefficient differences and analytic tails into C2 bounds."""

    if len(true_coefficient_intervals) != len(stored_coefficients):
        raise ValueError("true and stored coefficient orders must agree")
    radius = abs(_mp(z_radius))
    if not (0 < radius < 1):
        raise ValueError("z_radius must satisfy 0 < z_radius < 1")
    stored_mass = _mp(stored_mass)
    if stored_mass < 0:
        raise ValueError("stored mass must be nonnegative")
    order = len(stored_coefficients) - 1
    difference_abs: list[mp.mpf] = []
    for true_value, stored_value in zip(true_coefficient_intervals, stored_coefficients):
        lower, upper = _interval(true_value)
        stored = _mp(stored_value)
        difference_abs.append(max(abs(lower - stored), abs(upper - stored)))

    finite = {}
    for derivative, name in enumerate(("value", "first_Z", "second_Z")):
        finite[name] = mp.fsum(
            difference_abs[index] * _derivative_weight(index, derivative, radius)
            for index in range(order + 1)
        )

    true_mass_lower, true_mass_upper = _interval(true_coefficient_intervals[0])
    if true_mass_lower < 0:
        raise ValueError("true mass interval must be nonnegative")
    total_mass = true_mass_upper + stored_mass
    tails = weighted_geometric_tail(radius, order)
    coefficient_bound = {
        name: finite[name] + total_mass * tails[name]
        for name in ("value", "first_Z", "second_Z")
    }
    positive_fallback = {
        "value": total_mass,
        "first_Z": 4 * radius * total_mass,
        "second_Z": (4 + 24 * radius * radius) * total_mass,
    }
    final = {
        name: min(coefficient_bound[name], positive_fallback[name])
        for name in coefficient_bound
    }
    return {
        "uniform_error_upper_bounds": final,
        "coefficient_series_upper_bounds": coefficient_bound,
        "positive_mass_fallback_upper_bounds": positive_fallback,
        "finite_coefficient_contributions": finite,
        "coefficient_difference_abs_upper": difference_abs,
        "weighted_geometric_tail_factors": tails,
        "coefficient_order": order,
        "z_radius": radius,
        "true_mass_interval": (true_mass_lower, true_mass_upper),
        "stored_mass": stored_mass,
        "total_positive_mass": total_mass,
        "uniform_value_first_second_enclosed": False,
        "uniform_analytic_difference_enclosed": False,
        "finite_precision_roundoff_unenclosed": True,
        "nominal_mp_bound_conditional": True,
        "pressure_sign_invariant": True,
    }


def uniform_error_from_coefficients_directed(
    true_coefficient_intervals: Sequence[Any],
    stored_coefficients: Sequence[Any],
    stored_mass: Any,
    *,
    z_radius: Any = ".8",
    context: MPIntervalContext | None = None,
) -> dict[str, Any]:
    """Directed interval version of :func:`uniform_error_from_coefficients`."""

    if len(true_coefficient_intervals) != len(stored_coefficients):
        raise ValueError("true and stored coefficient orders must agree")
    if context is None:
        context = MPIntervalContext()
        context.dps = 80

    def as_interval(value: Any) -> Any:
        return value if hasattr(value, "_mpi_") else context.mpf(str(value))

    radius = context.mpf(str(abs(_mp(z_radius))))
    radius_lo, radius_hi = _interval(radius)
    if not (0 < radius_lo and radius_hi < 1):
        raise ValueError("z_radius must satisfy 0 < z_radius < 1")
    order = len(stored_coefficients) - 1
    difference_abs: list[mp.mpf] = []
    finite = {
        "value": context.mpf(0),
        "first_Z": context.mpf(0),
        "second_Z": context.mpf(0),
    }
    for index, (true_value, stored_value) in enumerate(
        zip(true_coefficient_intervals, stored_coefficients)
    ):
        difference = as_interval(true_value) - as_interval(stored_value)
        absolute = abs(difference)
        absolute_upper = _interval(absolute)[1]
        difference_abs.append(absolute_upper)
        for derivative, name in enumerate(("value", "first_Z", "second_Z")):
            finite[name] += context.mpf(str(absolute_upper)) * _derivative_weight(
                index, derivative, radius, context
            )
    true_mass = as_interval(true_coefficient_intervals[0])
    true_mass_lower, true_mass_upper = _interval(true_mass)
    if true_mass_lower < 0:
        raise ValueError("true mass interval must be nonnegative")
    stored_mass_interval = as_interval(stored_mass)
    if _interval(stored_mass_interval)[0] < 0:
        raise ValueError("stored mass must be nonnegative")
    total_mass = true_mass + stored_mass_interval
    tails = weighted_geometric_tail(radius, order, context=context)
    coefficient_bounds = {
        name: finite[name] + total_mass * tails[name]
        for name in ("value", "first_Z", "second_Z")
    }
    positive_fallback = {
        "value": total_mass,
        "first_Z": 4 * radius * total_mass,
        "second_Z": (4 + 24 * radius * radius) * total_mass,
    }
    coefficient_upper = {
        name: _interval(value)[1] for name, value in coefficient_bounds.items()
    }
    fallback_upper = {
        name: _interval(value)[1] for name, value in positive_fallback.items()
    }
    final = {
        name: min(coefficient_upper[name], fallback_upper[name])
        for name in coefficient_bounds
    }
    return {
        "uniform_error_upper_bounds": final,
        "coefficient_series_upper_bounds": coefficient_upper,
        "positive_mass_fallback_upper_bounds": fallback_upper,
        "finite_coefficient_contributions": {
            name: _interval(value)[1] for name, value in finite.items()
        },
        "coefficient_difference_abs_upper": difference_abs,
        "weighted_geometric_tail_factors": {
            name: _interval(value)[1] for name, value in tails.items()
        },
        "weighted_geometric_tail_intervals": tails,
        "coefficient_order": order,
        "z_radius": _interval(radius)[1],
        "true_mass_interval": (true_mass_lower, true_mass_upper),
        "stored_mass": _interval(stored_mass_interval)[1],
        "total_positive_mass": _interval(total_mass)[1],
        "uniform_value_first_second_enclosed": True,
        "uniform_analytic_difference_enclosed": True,
        "directed_interval_arithmetic": True,
        "finite_precision_roundoff_unenclosed": False,
        "pressure_sign_invariant": True,
    }


def _adapter_source(adapter: Any) -> tuple[Callable[[mp.mpf], tuple[mp.mpf, mp.mpf]], dict[str, Any]]:
    """Build the stored-schedule flatten density and verify its monotone stage."""

    schedule = adapter.schedule
    left, right = adapter._bounds("z_flatten")
    left, right = _mp(left), _mp(right)
    yv = _mp(getattr(adapter, "y_v", schedule.y_v))
    yf = _mp(getattr(adapter, "y_f", schedule.y_f))
    if left != yv or right != yf:
        raise ValueError("adapter z_flatten bounds do not match y_v/y_f")
    if right <= left:
        raise ValueError("adapter flatten interval must be positive")

    # On the stored flatten stage, the first two active switches are already
    # at one and the later switches are still at zero.  Thus
    # d(log A)/dy=-1/2-mu exactly.  The
    # switch beta=2(1-sigma((y-yv)/Tf)) is decreasing, so every positive
    # coefficient magnitude also decreases.  Check the stage inequalities
    # before relying on that analytic monotonicity statement.
    y_rel = getattr(schedule, "y_rel", None)
    ts = getattr(schedule, "Ts", None)
    mu = _mp(getattr(schedule, "mu", adapter.mu))
    if y_rel is None or ts is None or not (0 < mu < 1):
        raise ValueError("stored adapter lacks the flatten monotonicity checkpoints")
    y_rel, ts = _mp(y_rel), _mp(ts)
    if not (left >= 1 and right <= y_rel and right <= y_rel + 1 + ts):
        raise ValueError("flatten monotonicity checkpoints are not satisfied")

    def source(y: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        y = _mp(y)
        local = (y - yv) / (yf - yv)
        sigma, _ = ContinuousAxialPulse.sigma_pair(local)
        beta = 2 * (1 - sigma)
        log_amplitude = _mp(adapter._log_amplitude_ratio(y))
        weight = mp.exp(2 * log_amplitude + (beta - 3) * mp.log(2))
        return weight, beta

    return source, {
        "analytic_monotonicity_verified": True,
        "monotonicity_basis": (
            "beta decreases on [y_v,y_f], d log(A)/dy=-1/2-mu there, "
            "and w*abs(c_k(beta)) decreases"
        ),
        "stored_parameters_only": True,
        "original_parameter_errors_enclosed": False,
    }


def _adapter_interval_source(
    adapter: Any, enclosures: ScheduleEndpointEnclosures
) -> tuple[Callable[[Any], tuple[Any, Any]], dict[str, Any]]:
    """Build directed flatten ``(weight,beta)`` intervals from stored data."""

    schedule = adapter.schedule
    iv = enclosures.iv
    raw_left, raw_right = schedule._stage_bounds["z_flatten"]
    if raw_right - raw_left != 100:
        raise ValueError("stored z_flatten length must equal 100")
    left, right = _mp(raw_left), _mp(raw_right)
    yv = _mp(schedule.y_v)
    yf = _mp(schedule.y_f)
    if left != yv or right != yf:
        raise ValueError("stored z_flatten bounds do not match y_v/y_f")
    y_rel = _mp(schedule.y_rel)
    ts = _mp(schedule.Ts)
    y_d = _mp(schedule.y_d)
    mu = _mp(schedule.mu)
    if not (left >= 1 and left >= y_d + 1 and right <= y_rel and right <= y_rel + 1 + ts):
        raise ValueError("stored flatten monotonicity checkpoints are not satisfied")

    # The first and second incoming switches are saturated on this stage;
    # the relative/steep switches are still zero.  This is the exact stored
    # schedule slope, retained as an interval input rather than inferred from
    # sampled endpoint values.
    slope = -iv.mpf(".5") - enclosures.scalar(schedule.mu)
    left_log_amplitude = enclosures.log_amplitude_ratio(raw_left)["interval"]
    length_iv = enclosures.scalar(raw_right - raw_left)

    def source(y: Any) -> tuple[Any, Any]:
        # Integrate normalized u on [0,1]; retain the exact stored length
        # without converting the huge radial origins to nominal MP values.
        local = y if hasattr(y, "_mpi_") else iv.mpf(y)
        sigma = enclosures.sigma_interval(local)
        beta = 2 * (1 - sigma)
        ell = left_log_amplitude + slope * length_iv * local
        weight = length_iv * iv.exp(2 * ell + (beta - 3) * iv.ln(2))
        return weight, beta

    return source, {
        "analytic_monotonicity_verified": True,
        "monotonicity_basis": (
            "directed sigma_interval beta decrease plus exact stored slope "
            "d log(A)/dy=-1/2-mu on the 100-unit flatten stage"
        ),
        "flatten_slope": slope,
        "stored_flatten_length": length_iv,
        "stored_parameters_only": True,
        "original_parameter_errors_enclosed": False,
        "directed_interval_arithmetic": True,
    }


def _adapter_atoms(adapter: Any, quadrature_order: int) -> list[dict[str, Any]]:
    data = adapter._compute_components(int(quadrature_order))
    for component in data["post_components"]:
        if component.get("stage") == "z_flatten":
            return list(component["atoms"])
    raise ValueError("continuous preheat adapter has no z_flatten atom component")


def _enclose_uniform_flatten_error(
    adapter: Any,
    *,
    quadrature_order: int | None = None,
    coefficient_order: int = 96,
    panels: int = 1024,
    z_radius: Any = ".8",
) -> dict[str, Any]:
    """Enclose the stored adapter's flatten-stage error uniformly in ``Z``."""

    if quadrature_order is None:
        quadrature_order = int(adapter.quadrature_order)
    if isinstance(quadrature_order, bool) or int(quadrature_order) < 16:
        raise ValueError("quadrature_order must be at least 16")
    quadrature_order = int(quadrature_order)
    enclosures = ScheduleEndpointEnclosures(
        adapter.schedule, precision=max(80, int(getattr(adapter, "precision", 80)))
    )
    source, metadata = _adapter_interval_source(adapter, enclosures)
    left, right = adapter._bounds("z_flatten")
    true_data = monotone_coefficient_intervals_directed(
        source,
        0,
        1,
        coefficient_order=coefficient_order,
        panels=panels,
        monotone_decreasing=True,
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
    result.update(metadata)
    result.update(
        {
            "stage": "z_flatten",
            "left": _mp(left),
            "right": _mp(right),
            "stored_quadrature_order": quadrature_order,
            "stored_atom_count": len(atoms),
            "true_coefficient_panels": true_data["panels"],
            "true_coefficient_interval_widths": true_data[
                "coefficient_interval_widths"
            ],
            "true_mass_interval": true_data["mass_interval"],
            "monotone_weight_diagnostic": true_data[
                "monotone_weight_diagnostic"
            ],
            "monotone_beta_diagnostic": true_data["monotone_beta_diagnostic"],
            "finite_gauss_atom_error_only": True,
            "quadrature_roundoff_enclosed": False,
            "source_endpoint_and_tail_roundoff_enclosed": True,
            "stored_gauss_atom_generation_roundoff_unenclosed": True,
            "core_error_enclosed": False,
            "five_defect_interval_closure": False,
        }
    )
    return result


def enclose_uniform_flatten_error(
    adapter: Any,
    *,
    quadrature_order: int | None = None,
    coefficient_order: int = 96,
    panels: int = 1024,
    z_radius: Any = ".8",
) -> dict[str, Any]:
    """Run the stored-adapter enclosure at its declared MP precision."""

    precision = max(40, int(getattr(adapter, "precision", 80)))
    with mp.workdps(precision):
        return _enclose_uniform_flatten_error(
            adapter,
            quadrature_order=quadrature_order,
            coefficient_order=coefficient_order,
            panels=panels,
            z_radius=z_radius,
        )


# Descriptive aliases for callers following the existing pressure adapters.
uniform_flatten_error = enclose_uniform_flatten_error
bound_uniform_flatten_error = enclose_uniform_flatten_error


__all__ = [
    "weighted_geometric_tail",
    "stored_flatten_coefficients",
    "stored_flatten_coefficients_directed",
    "monotone_coefficient_intervals",
    "monotone_coefficient_intervals_directed",
    "uniform_error_from_coefficients",
    "uniform_error_from_coefficients_directed",
    "enclose_uniform_flatten_error",
    "uniform_flatten_error",
    "bound_uniform_flatten_error",
]
