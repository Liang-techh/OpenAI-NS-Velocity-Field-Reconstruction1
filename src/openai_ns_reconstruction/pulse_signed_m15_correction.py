"""Signed finite-endpoint pulse correction from physical source jets.

This module is deliberately a small numerical kernel, not a source binding.
The caller supplies two *real* physical field callables.  A field callable is
queried as ``field(chart_name, x, Z)`` at a shared source query for both sides
of each chart and returns the seven physical jets
``U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ``.  The caller also supplies a
source-owned radius callable ``radius(chart_name, x)`` for each side and must
explicitly declare that it is independent of ``Z``.

For the ordered pair ``physical - frame`` the kernel evaluates, chart by
chart, the stable signed difference

``R_P**2 * (S_P - S_F) + (R_P - R_F) * (R_P + R_F) * S_F``

and the corresponding endpoint term

``R_P * (Q_yP - Q_yF) + (R_P - R_F) * Q_yF``.

The ordinary ``Z1`` derivative is formed analytically before integration.  In
particular, ``Q_yZ`` and ``Q_ZZ`` are consumed; they are never silently
dropped.  ``main`` and ``gap`` use ``dy = dxi / mu`` and ``end`` uses
``dy = ds``.  Every chart carries its own finite source endpoints, while the
field callables use each chart's own source-query variable ``x``.  Endpoint terms
are retained with the orientation ``upper - lower``.

``scipy.integrate.quad_vec`` supplies a quadrature error estimate.  That
estimate is diagnostic only and is not a certification or a pulse-equality
proof.  The kernel does not bind source recipes, repeat ``Pstar``, assign
physical units, identify physical entities, or prove that the signed pulse
correction vanishes.  Those are caller responsibilities.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Literal

import numpy as np
from scipy.integrate import quad_vec


ChartName = Literal["main", "gap", "end"]
FieldCallable = Callable[[str, float, float], "FieldJets | Mapping[str, float]"]
RadiusCallable = Callable[[str, float], float]


_FIELD_NAMES = ("U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ")
_CHART_NAMES = ("main", "gap", "end")


def _finite(value: float, label: str) -> float:
    """Convert one source value to a finite real scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


@dataclass(frozen=True)
class FieldJets:
    """Physical source jets required by the pulse formula."""

    U: float
    U_Z: float
    Q: float
    Q_Z: float
    Q_ZZ: float
    Q_y: float
    Q_yZ: float

    def finite(self, *, label: str = "field jets") -> "FieldJets":
        """Return a finite scalar copy, rejecting incomplete source data."""

        values = {
            name: _finite(getattr(self, name), f"{label}.{name}")
            for name in _FIELD_NAMES
        }
        return type(self)(**values)


def _coerce_field_jets(value: FieldJets | Mapping[str, float], *, label: str) -> FieldJets:
    if isinstance(value, FieldJets):
        return value.finite(label=label)
    if not isinstance(value, Mapping):
        raise ValueError(
            f"{label} must return FieldJets or a mapping with {_FIELD_NAMES}"
        )
    missing = [name for name in _FIELD_NAMES if name not in value]
    if missing:
        raise ValueError(f"{label} is missing field jets: {', '.join(missing)}")
    return FieldJets(
        **{name: _finite(value[name], f"{label}.{name}") for name in _FIELD_NAMES}
    )


@dataclass(frozen=True)
class PulseSide:
    """One ordered side of the signed pair.

    ``field`` and ``radius`` are source-owned callables.  The first side
    passed to :func:`evaluate_signed_pulse_correction` is ``P`` and the
    second is ``F`` in the stable pair formulas in the module docstring.
    ``radius_z_independent`` is an explicit caller assertion; a false value
    is rejected rather than silently omitting a Leibniz term.
    ``source_radius_exponential_law`` is a second explicit assertion required
    by the endpoint integration-by-parts identity ``dR/dy = R``.
    """

    name: str
    field: FieldCallable
    radius: RadiusCallable
    radius_z_independent: bool
    source_radius_exponential_law: bool

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("PulseSide.name must be nonempty")
        if not callable(self.field):
            raise TypeError("PulseSide.field must be callable")
        if not callable(self.radius):
            raise TypeError("PulseSide.radius must be callable")


@dataclass(frozen=True)
class SourceInterval:
    """Optional fallback source interval for charts without local endpoints.

    New callers should put finite endpoints directly on each ``ChartSpec``.
    The fallback is retained for a genuinely common interval whose endpoints
    are explicitly declared Z independent.
    """

    lower: float
    upper: float
    z_independent: bool

    def finite(self) -> "SourceInterval":
        lower = _finite(self.lower, "source lower endpoint")
        upper = _finite(self.upper, "source upper endpoint")
        return type(self)(lower, upper, self.z_independent)


@dataclass(frozen=True)
class ChartSpec:
    """One of the three physical chart measures.

    ``main`` and ``gap`` must use coordinate ``xi`` and the physical measure
    ``1 / mu``.  ``end`` must use coordinate ``s`` and unit physical measure.
    Each chart may carry its own ``lower`` and ``upper`` endpoints; for
    example, ``[0, 11]``, ``[11, 13 - 4*mu]`` and ``[-4, 0]``.  ``mu`` is a
    scalar rather than a callable on purpose: a Z-dependent Jacobian would add
    Leibniz terms that this kernel cannot infer.
    """

    name: ChartName
    coordinate: str
    mu: float | None
    mu_z_independent: bool
    lower: float | None = None
    upper: float | None = None
    endpoints_z_independent: bool = False

    def weight(self) -> float:
        if self.name not in _CHART_NAMES:
            raise ValueError(f"unknown chart {self.name!r}; expected {_CHART_NAMES}")
        if self.name in ("main", "gap"):
            if self.coordinate != "xi":
                raise ValueError(f"{self.name} chart must use source coordinate xi")
            if self.mu is None:
                raise ValueError(f"{self.name} chart requires mu")
            if not self.mu_z_independent:
                raise ValueError(
                    f"{self.name} mu must be declared Z independent; "
                    "otherwise its Leibniz term is not available"
                )
            mu = _finite(self.mu, f"{self.name} mu")
            if mu <= 0:
                raise ValueError(f"{self.name} mu must be positive")
            return 1.0 / mu
        if self.coordinate != "s":
            raise ValueError("end chart must use source coordinate s")
        if self.mu is not None:
            raise ValueError("end chart uses dy=ds and must not provide mu")
        return 1.0

    def interval(self, fallback: SourceInterval | None) -> SourceInterval:
        """Return this chart's local interval, or a validated fallback."""

        has_lower = self.lower is not None
        has_upper = self.upper is not None
        if has_lower != has_upper:
            raise ValueError(f"{self.name} chart must provide both endpoints")
        if has_lower:
            if not self.endpoints_z_independent:
                raise ValueError(
                    f"{self.name} endpoints must be declared Z independent; "
                    "otherwise endpoint Leibniz terms are unavailable"
                )
            return SourceInterval(self.lower, self.upper, True).finite()
        if fallback is None:
            raise ValueError(
                f"{self.name} chart requires own endpoints or a source_interval fallback"
            )
        interval = fallback.finite()
        if not interval.z_independent:
            raise ValueError(
                "finite source endpoints must be declared Z independent; "
                "otherwise endpoint Leibniz terms are unavailable"
            )
        return interval


@dataclass(frozen=True)
class ResolvedFormulaParameters:
    """Formula coefficients at one ordinary slow coordinate ``Z``."""

    Z: float
    delta: float
    h: float
    L: float
    L_Z: float


@dataclass(frozen=True)
class FormulaParameters:
    """Z-independent ``delta`` source parameter for the physical formula.

    The formula coefficients are resolved at each requested ``Z`` as
    ``h = 1 - delta``, ``L = 1 - delta * Z**2`` and
    ``L_Z = -2 * delta * Z``.  Keeping ``delta`` as the only caller input
    prevents the ordinary Z1 derivative from silently dropping the varying
    denominator term.
    """

    delta: float

    def resolve(self, Z: float) -> ResolvedFormulaParameters:
        delta = _finite(self.delta, "delta")
        z = _finite(Z, "Z")
        h = 1.0 - delta
        L = 1.0 - delta * z * z
        L_Z = -2.0 * delta * z
        if L == 0:
            raise ValueError("L = 1 - delta*Z**2 must be nonzero")
        return ResolvedFormulaParameters(z, delta, h, L, L_Z)


@dataclass(frozen=True)
class QuadratureSettings:
    """Diagnostic settings passed to ``scipy.integrate.quad_vec``."""

    epsabs: float = 1.0e-10
    epsrel: float = 1.0e-10
    limit: int = 200

    def finite(self) -> "QuadratureSettings":
        epsabs = _finite(self.epsabs, "quadrature epsabs")
        epsrel = _finite(self.epsrel, "quadrature epsrel")
        if epsabs <= 0 or epsrel <= 0:
            raise ValueError("quadrature epsabs and epsrel must be positive")
        if int(self.limit) != self.limit or self.limit < 1:
            raise ValueError("quadrature limit must be a positive integer")
        return type(self)(epsabs, epsrel, int(self.limit))


@dataclass(frozen=True)
class ChartSignedPulseResult:
    """Stable signed ``P-F`` result for one chart."""

    chart: str
    coordinate: str
    source_lower: float
    source_upper: float
    physical_weight: float
    physical_side: str
    frame_side: str
    delta_integral: float
    radius_integral: float
    delta_integral_Z1: float
    radius_integral_Z1: float
    integral: float
    integral_Z1: float
    boundary_delta_upper: float
    boundary_radius_upper: float
    boundary_delta_lower: float
    boundary_radius_lower: float
    boundary_delta_upper_Z1: float
    boundary_radius_upper_Z1: float
    boundary_delta_lower_Z1: float
    boundary_radius_lower_Z1: float
    boundary: float
    boundary_Z1: float
    m15_X: float
    m15_X_Z1: float
    quadrature_error_estimate: float
    quadrature_is_certification: bool = False

    @property
    def z1(self) -> float:
        """Short alias for the ordinary ``Z1`` derivative."""

        return self.m15_X_Z1


CALLER_RESPONSIBILITIES = (
    "map each side to the real physical source field and physical radius",
    "use one common source query x while declaring each chart's finite endpoints",
    "provide physical U, U_Z, Q, Q_Z, Q_ZZ, Q_y and Q_yZ jets",
    "declare the source radius Z independence and dR/dy = R exponential law",
    "declare main/gap mu and each finite endpoint Z independent",
    "supply Pstar, units and real-entity identification outside this kernel",
    "prove any pulse equality or cancellation outside this diagnostic kernel",
)


@dataclass(frozen=True)
class SignedPulseCorrection:
    """All three chart pair differences and their signed aggregate."""

    Z: float
    physical_side: str
    frame_side: str
    charts: tuple[ChartSignedPulseResult, ...]
    signed_m15: float
    signed_m15_Z1: float
    quadrature_error_estimate: float
    quadrature_is_certification: bool = False
    caller_responsibilities: tuple[str, ...] = CALLER_RESPONSIBILITIES
    pulse_equality_proved: bool = False

    @property
    def chart_map(self) -> dict[str, ChartSignedPulseResult]:
        return {chart.chart: chart for chart in self.charts}


def source_density_and_Z1(
    jets: FieldJets | Mapping[str, float],
    *,
    Z: float,
    parameters: FormulaParameters,
) -> tuple[float, float]:
    """Return ``(S, S_Z)`` from one physical field jet.

    This is the ordinary slow ``Z1`` derivative at fixed source query ``x``.
    ``Q_yZ`` is part of the required input and appears in both the linear and
    nonlinear derivative terms.
    """

    q = _coerce_field_jets(jets, label="source field")
    z = _finite(Z, "Z")
    resolved = parameters.resolve(z)
    h, L, L_Z = resolved.h, resolved.L, resolved.L_Z

    A = q.Q + 0.5 * h * z * q.Q_Z + q.Q_y
    A_Z = q.Q_Z + 0.5 * h * (q.Q_Z + z * q.Q_ZZ) + q.Q_yZ

    B = q.Q * (0.5 * q.Q + q.Q_y)
    B_Z = q.Q_Z * (0.5 * q.Q + q.Q_y) + q.Q * (0.5 * q.Q_Z + q.Q_yZ)

    C = -2.0 * z * q.Q + (1.0 - z * z) * q.Q_Z - 2.0 * z * q.Q_y
    C_Z = (
        -2.0 * q.Q
        - 4.0 * z * q.Q_Z
        + (1.0 - z * z) * q.Q_ZZ
        - 2.0 * q.Q_y
        - 2.0 * z * q.Q_yZ
    )

    S = A / L + B + q.U * C / L
    numerator_Z = A_Z + q.U_Z * C + q.U * C_Z
    S_Z = numerator_Z / L + B_Z - (A + q.U * C) * L_Z / (L * L)
    return _finite(S, "S"), _finite(S_Z, "S_Z")


def _radius(side: PulseSide, chart: str, x: float) -> float:
    value = _finite(side.radius(chart, x), f"{side.name} radius at {chart}, x={x}")
    if value <= 0:
        raise ValueError(f"{side.name} radius must be positive")
    return value


def _boundary_components(
    physical: PulseSide,
    frame: PulseSide,
    chart: str,
    x: float,
    Z: float,
    parameters: FormulaParameters,
) -> np.ndarray:
    p = _coerce_field_jets(
        physical.field(chart, x, Z), label=f"{physical.name} field at {chart}, x={x}"
    )
    f = _coerce_field_jets(
        frame.field(chart, x, Z), label=f"{frame.name} field at {chart}, x={x}"
    )
    rp = _radius(physical, chart, x)
    rf = _radius(frame, chart, x)
    dr = rp - rf
    return np.asarray(
        (
            rp * (p.Q_y - f.Q_y),
            dr * f.Q_y,
            rp * (p.Q_yZ - f.Q_yZ),
            dr * f.Q_yZ,
        ),
        dtype=float,
    )


def _pair_integrand(
    physical: PulseSide,
    frame: PulseSide,
    chart: str,
    x: float,
    Z: float,
    parameters: FormulaParameters,
    weight: float,
) -> np.ndarray:
    p = _coerce_field_jets(
        physical.field(chart, x, Z), label=f"{physical.name} field at {chart}, x={x}"
    )
    f = _coerce_field_jets(
        frame.field(chart, x, Z), label=f"{frame.name} field at {chart}, x={x}"
    )
    sp, sp_Z = source_density_and_Z1(p, Z=Z, parameters=parameters)
    sf, sf_Z = source_density_and_Z1(f, Z=Z, parameters=parameters)
    rp = _radius(physical, chart, x)
    rf = _radius(frame, chart, x)
    dr = rp - rf
    radius_square_difference = dr * (rp + rf)
    # P-F is evaluated in a difference basis before quadrature.  This avoids
    # integrating two large values and subtracting them only at the end.
    return 0.5 * weight * np.asarray(
        (
            rp * rp * (sp - sf),
            radius_square_difference * sf,
            rp * rp * (sp_Z - sf_Z),
            radius_square_difference * sf_Z,
        ),
        dtype=float,
    )


def _validate_inputs(
    physical: PulseSide,
    frame: PulseSide,
    charts: Sequence[ChartSpec],
    source_interval: SourceInterval | None,
    parameters: FormulaParameters,
    settings: QuadratureSettings,
    Z: float,
) -> tuple[
    tuple[ChartSpec, ...],
    dict[str, SourceInterval],
    FormulaParameters,
    QuadratureSettings,
    float,
]:
    if not physical.radius_z_independent:
        raise ValueError(
            f"{physical.name} radius must be declared Z independent; "
            "otherwise its Leibniz term is unavailable"
        )
    if not frame.radius_z_independent:
        raise ValueError(
            f"{frame.name} radius must be declared Z independent; "
            "otherwise its Leibniz term is unavailable"
        )
    if not physical.source_radius_exponential_law:
        raise ValueError(
            f"{physical.name} must declare source_radius_exponential_law "
            "(dR/dy = R) for the endpoint integration-by-parts identity"
        )
    if not frame.source_radius_exponential_law:
        raise ValueError(
            f"{frame.name} must declare source_radius_exponential_law "
            "(dR/dy = R) for the endpoint integration-by-parts identity"
        )
    z = _finite(Z, "Z")
    normalized_parameters = parameters
    normalized_parameters.resolve(z)
    normalized_settings = settings.finite()
    normalized_charts = tuple(charts)
    if len(normalized_charts) != len(_CHART_NAMES):
        raise ValueError("exactly main, gap and end charts are required")
    names = tuple(chart.name for chart in normalized_charts)
    if set(names) != set(_CHART_NAMES) or len(set(names)) != len(_CHART_NAMES):
        raise ValueError("charts must contain main, gap and end exactly once")
    intervals: dict[str, SourceInterval] = {}
    for chart in normalized_charts:
        chart.weight()
        intervals[chart.name] = chart.interval(source_interval)
    return normalized_charts, intervals, normalized_parameters, normalized_settings, z


def evaluate_signed_pulse_correction(
    *,
    physical_side: PulseSide,
    frame_side: PulseSide,
    charts: Sequence[ChartSpec],
    parameters: FormulaParameters,
    Z: float,
    source_interval: SourceInterval | None = None,
    settings: QuadratureSettings = QuadratureSettings(),
) -> SignedPulseCorrection:
    """Evaluate signed ``m15_P - m15_F`` on all three physical charts.

    The returned chart values are

    ``m15_X = 1/2 * integral R_X^2 S_X dy - [R_X Q_y]``

    after stable ``P-F`` decomposition.  ``m15_X_Z1`` is its ordinary
    analytic derivative.  The aggregate is the direct sum of the three chart
    differences; no internal seam or pulse correction is assumed to vanish.
    """

    normalized_charts, intervals, p, q, z = _validate_inputs(
        physical_side,
        frame_side,
        charts,
        source_interval,
        parameters,
        settings,
        Z,
    )

    results: list[ChartSignedPulseResult] = []
    for chart in normalized_charts:
        weight = chart.weight()
        interval = intervals[chart.name]
        raw_integral, raw_error, quadrature_info = quad_vec(
            lambda x, chart_name=chart.name, chart_weight=weight: _pair_integrand(
                physical_side,
                frame_side,
                chart_name,
                float(x),
                z,
                p,
                chart_weight,
            ),
            interval.lower,
            interval.upper,
            epsabs=q.epsabs,
            epsrel=q.epsrel,
            norm="max",
            limit=q.limit,
            full_output=True,
        )
        if not quadrature_info.success:
            raise ArithmeticError(
                f"{chart.name} quadrature did not converge: {quadrature_info.message}"
            )
        integral_components = tuple(_finite(value, f"{chart.name} integral component") for value in raw_integral)
        if len(integral_components) != 4:
            raise ArithmeticError("pair integrand must return four stable components")
        boundary_upper = tuple(
            _finite(value, f"{chart.name} upper boundary component")
            for value in _boundary_components(
                physical_side, frame_side, chart.name, interval.upper, z, p
            )
        )
        boundary_lower = tuple(
            _finite(value, f"{chart.name} lower boundary component")
            for value in _boundary_components(
                physical_side, frame_side, chart.name, interval.lower, z, p
            )
        )
        boundary_components = tuple(
            upper - lower for upper, lower in zip(boundary_upper, boundary_lower)
        )
        delta_integral, radius_integral, delta_integral_Z1, radius_integral_Z1 = integral_components
        boundary_delta_upper, boundary_radius_upper, boundary_delta_upper_Z1, boundary_radius_upper_Z1 = boundary_upper
        boundary_delta_lower, boundary_radius_lower, boundary_delta_lower_Z1, boundary_radius_lower_Z1 = boundary_lower
        boundary_delta, boundary_radius, boundary_delta_Z1, boundary_radius_Z1 = boundary_components
        integral = delta_integral + radius_integral
        integral_Z1 = delta_integral_Z1 + radius_integral_Z1
        boundary = boundary_delta + boundary_radius
        boundary_Z1 = boundary_delta_Z1 + boundary_radius_Z1
        results.append(
            ChartSignedPulseResult(
                chart=chart.name,
                coordinate=chart.coordinate,
                source_lower=interval.lower,
                source_upper=interval.upper,
                physical_weight=weight,
                physical_side=physical_side.name,
                frame_side=frame_side.name,
                delta_integral=delta_integral,
                radius_integral=radius_integral,
                delta_integral_Z1=delta_integral_Z1,
                radius_integral_Z1=radius_integral_Z1,
                integral=integral,
                integral_Z1=integral_Z1,
                boundary_delta_upper=boundary_delta_upper,
                boundary_radius_upper=boundary_radius_upper,
                boundary_delta_lower=boundary_delta_lower,
                boundary_radius_lower=boundary_radius_lower,
                boundary_delta_upper_Z1=boundary_delta_upper_Z1,
                boundary_radius_upper_Z1=boundary_radius_upper_Z1,
                boundary_delta_lower_Z1=boundary_delta_lower_Z1,
                boundary_radius_lower_Z1=boundary_radius_lower_Z1,
                boundary=boundary,
                boundary_Z1=boundary_Z1,
                m15_X=integral - boundary,
                m15_X_Z1=integral_Z1 - boundary_Z1,
                quadrature_error_estimate=_finite(
                    raw_error, f"{chart.name} quadrature error estimate"
                ),
            )
        )

    return SignedPulseCorrection(
        Z=z,
        physical_side=physical_side.name,
        frame_side=frame_side.name,
        charts=tuple(results),
        signed_m15=sum(result.m15_X for result in results),
        signed_m15_Z1=sum(result.m15_X_Z1 for result in results),
        quadrature_error_estimate=sum(
            result.quadrature_error_estimate for result in results
        ),
    )


# A descriptive alias for callers that prefer ``compute`` terminology.
compute_signed_pulse_correction = evaluate_signed_pulse_correction


__all__ = [
    "CALLER_RESPONSIBILITIES",
    "ChartSignedPulseResult",
    "ChartSpec",
    "FieldJets",
    "FormulaParameters",
    "PulseSide",
    "QuadratureSettings",
    "ResolvedFormulaParameters",
    "SignedPulseCorrection",
    "SourceInterval",
    "compute_signed_pulse_correction",
    "evaluate_signed_pulse_correction",
    "source_density_and_Z1",
]
