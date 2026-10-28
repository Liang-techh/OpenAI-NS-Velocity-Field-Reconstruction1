"""Exact source-scale coordinate maps for the three pulse charts.

This module owns only coordinate and scale bookkeeping.  It deliberately does
not evaluate a physical field, bind a source consumer, compare a graph, or
claim pulse equality.  The returned objects keep ``mu`` and ``log_Rp`` as
SymPy expressions.  In particular, the actual source value

``mu = exp(-4 * (exp(40) + 11)) / 1000``

is returned as an unevaluated symbolic replacement, never as a binary64
value.  No radius value is formed: the only radius quantity exposed here is
``log_Rp + y``.

The gap chart is also exposed in its exact end-compatible parameterisation
``xi = 13 + mu * s``.  Consequently ``s`` runs from ``-2 / mu`` to ``-4``
and ``d xi / mu = d s``.  The end chart uses ``s`` in ``[-4, 0]``.  The
ordinary-Z independence flag is an explicit caller assumption and is always
marked as not automatically certified by this small kernel.

The optional ``sympy`` dependency is intentional.  Expressions are accepted
only when they are exact SymPy expressions (or integers/rationals that SymPy
can represent exactly); Python ``float`` values are rejected so a source
underflow cannot be hidden by a later conversion to arbitrary precision.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Mapping

import sympy as sp


ChartName = Literal["main", "gap", "end"]


MU_SYMBOL = sp.Symbol("mu", positive=True)
XI_SYMBOL = sp.Symbol("xi", real=True)
S_SYMBOL = sp.Symbol("s", real=True)
LOG_RP_SYMBOL = sp.Symbol("log_Rp", real=True)
SOURCE_MU_RECIPE = "mu = exp(-4 * (exp(40) + 11)) / 1000"


def source_mu_expression() -> sp.Expr:
    """Return the actual source ``mu`` as an exact SymPy expression.

    The expression is intentionally not numerically evaluated.  Its returned
    object contains the exact nested exponential source recipe and can be used
    as a replacement for :data:`MU_SYMBOL` in a separate source-scale view.
    """

    return sp.Rational(1, 1000) * sp.exp(
        -sp.Integer(4) * (sp.exp(sp.Integer(40)) + sp.Integer(11))
    )


# A descriptive alias makes the source replacement easy to discover without
# introducing another source expression or any numerical conversion.
real_source_mu_expression = source_mu_expression


def _exact_expr(value: Any, *, label: str) -> sp.Expr:
    """Coerce exact scalar input while rejecting binary floating values."""

    if isinstance(value, float):
        raise TypeError(f"{label} must not be a Python float; use SymPy Rational")
    expr = sp.sympify(value)
    if not isinstance(expr, sp.Expr):
        raise TypeError(f"{label} must be a SymPy scalar expression")
    if expr.has(sp.Float):
        raise TypeError(f"{label} must not contain a SymPy Float")
    return expr


def _positive_mu(value: Any) -> sp.Expr:
    """Validate a positive exact ``mu`` without materialising its value."""

    expr = _exact_expr(value, label="mu")
    if expr.is_positive is False:
        raise ValueError("mu must be positive")
    if expr.is_positive is not True:
        raise ValueError(
            "mu must carry a SymPy positive assumption; "
            "use Symbol(..., positive=True) or an exact positive expression"
        )
    return expr


@dataclass(frozen=True)
class DirectedBounds:
    """Exact oriented bounds under the stated source assumptions.

    ``directed`` describes the lower-to-upper orientation encoded by the two
    expressions.  It is not a numerical interval certificate: the ordering
    assumption is recorded explicitly and ``automatic_certificate`` remains
    false.  This distinction matters for the actual source, which is never
    materialised here.
    """

    lower: sp.Expr
    upper: sp.Expr
    ordering_assumption: str
    directed: bool = True
    automatic_certificate: bool = False
    numeric_materialized: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "lower", _exact_expr(self.lower, label="bound lower"))
        object.__setattr__(self, "upper", _exact_expr(self.upper, label="bound upper"))
        if not self.ordering_assumption:
            raise ValueError("DirectedBounds.ordering_assumption must be nonempty")

    def as_dict(self) -> dict[str, Any]:
        return {
            "lower": _expr_text(self.lower),
            "upper": _expr_text(self.upper),
            "ordering_assumption": self.ordering_assumption,
            "directed": self.directed,
            "automatic_certificate": self.automatic_certificate,
            "numeric_materialized": self.numeric_materialized,
        }


@dataclass(frozen=True)
class OrdinaryZIndependence:
    """Explicit ordinary-Z independence declaration.

    The declaration is accepted as a caller assumption only.  This coordinate
    module cannot infer or certify that a source radius, ``mu`` or endpoint is
    independent of ``Z``.
    """

    declared: bool
    automatic_certificate: bool = False

    @property
    def status(self) -> str:
        if self.declared:
            return "declared_assumption_not_automatically_certified"
        return "not_declared"

    def as_dict(self) -> dict[str, Any]:
        return {
            "declared": self.declared,
            "automatic_certificate": self.automatic_certificate,
            "status": self.status,
            "statement": (
                "mu, log_Rp, endpoints and the source coordinate map are "
                "ordinary-Z independent by caller assumption"
            ),
        }


def _expr_text(value: sp.Expr) -> str:
    """Stable human/JSON representation of one exact expression."""

    return sp.sstr(value)


@dataclass(frozen=True)
class ExactChartMap:
    """One exact chart map and its source-scale metadata."""

    name: ChartName
    coordinate: sp.Symbol
    coordinate_to_y: sp.Expr
    source_bounds: DirectedBounds
    y_bounds: DirectedBounds
    log_radius: sp.Expr
    log_radius_bounds: DirectedBounds
    physical_measure_density: sp.Expr
    physical_measure_differential: sp.Expr
    physical_measure_identity: str
    beta: sp.Expr
    theta_exponent_original: sp.Expr
    theta_exponent_log_rewrite: sp.Expr
    ordinary_z_independence: OrdinaryZIndependence
    xi_from_s: sp.Expr | None = None
    s_from_xi: sp.Expr | None = None
    reparameterized_coordinate: sp.Symbol | None = None
    reparameterized_to_y: sp.Expr | None = None
    reparameterized_bounds: DirectedBounds | None = None
    reparameterized_log_radius: sp.Expr | None = None
    reparameterized_measure_density: sp.Expr | None = None
    reparameterized_measure_differential: sp.Expr | None = None
    reparameterized_measure_identity: str | None = None
    reparameterized_theta_exponent_original: sp.Expr | None = None
    reparameterized_theta_exponent_log_rewrite: sp.Expr | None = None
    radius_value_materialized: bool = False

    def __post_init__(self) -> None:
        if self.name not in ("main", "gap", "end"):
            raise ValueError(f"unknown chart {self.name!r}")
        if not isinstance(self.coordinate, sp.Symbol):
            raise TypeError("coordinate must be a SymPy Symbol")
        for label in (
            "coordinate_to_y",
            "log_radius",
            "physical_measure_density",
            "physical_measure_differential",
            "beta",
            "theta_exponent_original",
            "theta_exponent_log_rewrite",
        ):
            object.__setattr__(
                self,
                label,
                _exact_expr(getattr(self, label), label=label),
            )
        optional_exprs = (
            "xi_from_s",
            "s_from_xi",
            "reparameterized_to_y",
            "reparameterized_log_radius",
            "reparameterized_measure_density",
            "reparameterized_measure_differential",
            "reparameterized_theta_exponent_original",
            "reparameterized_theta_exponent_log_rewrite",
        )
        for label in optional_exprs:
            value = getattr(self, label)
            if value is not None:
                object.__setattr__(self, label, _exact_expr(value, label=label))
        if self.name == "gap":
            required = (
                self.xi_from_s,
                self.s_from_xi,
                self.reparameterized_coordinate,
                self.reparameterized_to_y,
                self.reparameterized_bounds,
                self.reparameterized_log_radius,
                self.reparameterized_measure_density,
                self.reparameterized_measure_differential,
                self.reparameterized_measure_identity,
                self.reparameterized_theta_exponent_original,
                self.reparameterized_theta_exponent_log_rewrite,
            )
            if any(value is None for value in required):
                raise ValueError("gap map must expose its exact s reparameterisation")

    @property
    def coordinate_name(self) -> str:
        return str(self.coordinate)

    @property
    def ordinary_z_independent_declared(self) -> bool:
        return self.ordinary_z_independence.declared

    @property
    def ordinary_z_independent_certificate(self) -> bool:
        return self.ordinary_z_independence.automatic_certificate

    def source_substituted(self, source_mu: sp.Expr | None = None) -> dict[str, sp.Expr]:
        """Return a symbolic source view without evaluating any exponentials."""

        if self.beta - sp.Rational(1, 2) != MU_SYMBOL:
            raise ValueError("source substitution requires a system built with the default mu symbol")
        source = source_mu_expression() if source_mu is None else _positive_mu(source_mu)
        replacements = {MU_SYMBOL: source}
        values: dict[str, sp.Expr] = {
            "coordinate_to_y": self.coordinate_to_y,
            "source_lower": self.source_bounds.lower,
            "source_upper": self.source_bounds.upper,
            "y_lower": self.y_bounds.lower,
            "y_upper": self.y_bounds.upper,
            "log_radius_lower": self.log_radius_bounds.lower,
            "log_radius_upper": self.log_radius_bounds.upper,
            "log_radius": self.log_radius,
            "physical_measure_density": self.physical_measure_density,
            "physical_measure_differential": self.physical_measure_differential,
            "theta_exponent_original": self.theta_exponent_original,
            "theta_exponent_log_rewrite": self.theta_exponent_log_rewrite,
        }
        if self.reparameterized_bounds is not None:
            values["reparameterized_lower"] = self.reparameterized_bounds.lower
            values["reparameterized_upper"] = self.reparameterized_bounds.upper
        for label in (
            "xi_from_s",
            "s_from_xi",
            "reparameterized_to_y",
            "reparameterized_log_radius",
            "reparameterized_measure_density",
            "reparameterized_measure_differential",
            "reparameterized_theta_exponent_original",
            "reparameterized_theta_exponent_log_rewrite",
        ):
            value = getattr(self, label)
            if value is not None:
                values[label] = value
        return {
            label: expression.subs(replacements, simultaneous=True)
            for label, expression in values.items()
        }

    def as_dict(self, *, include_source_view: bool = False) -> dict[str, Any]:
        result: dict[str, Any] = {
            "name": self.name,
            "coordinate": self.coordinate_name,
            "coordinate_to_y": _expr_text(self.coordinate_to_y),
            "source_bounds": self.source_bounds.as_dict(),
            "y_bounds": self.y_bounds.as_dict(),
            "log_radius": _expr_text(self.log_radius),
            "log_radius_bounds": self.log_radius_bounds.as_dict(),
            "physical_measure_density": _expr_text(self.physical_measure_density),
            "physical_measure_differential": _expr_text(
                self.physical_measure_differential
            ),
            "physical_measure_identity": self.physical_measure_identity,
            "beta": _expr_text(self.beta),
            "theta_exponent_original": _expr_text(self.theta_exponent_original),
            "theta_exponent_log_rewrite": _expr_text(
                self.theta_exponent_log_rewrite
            ),
            "ordinary_z_independence": self.ordinary_z_independence.as_dict(),
            "radius_value_materialized": self.radius_value_materialized,
        }
        if self.name == "gap":
            result.update(
                {
                    "xi_from_s": _expr_text(self.xi_from_s),
                    "s_from_xi": _expr_text(self.s_from_xi),
                    "reparameterized_coordinate": str(self.reparameterized_coordinate),
                    "reparameterized_to_y": _expr_text(self.reparameterized_to_y),
                    "reparameterized_bounds": self.reparameterized_bounds.as_dict(),
                    "reparameterized_log_radius": _expr_text(
                        self.reparameterized_log_radius
                    ),
                    "reparameterized_measure_density": _expr_text(
                        self.reparameterized_measure_density
                    ),
                    "reparameterized_measure_differential": _expr_text(
                        self.reparameterized_measure_differential
                    ),
                    "reparameterized_measure_identity": self.reparameterized_measure_identity,
                    "reparameterized_theta_exponent_original": _expr_text(
                        self.reparameterized_theta_exponent_original
                    ),
                    "reparameterized_theta_exponent_log_rewrite": _expr_text(
                        self.reparameterized_theta_exponent_log_rewrite
                    ),
                }
            )
        if self.name == "end":
            result.update(
                {
                    "xi_from_s": _expr_text(self.xi_from_s),
                    "reparameterized_to_y": _expr_text(self.reparameterized_to_y),
                    "reparameterized_bounds": self.reparameterized_bounds.as_dict(),
                    "reparameterized_log_radius": _expr_text(
                        self.reparameterized_log_radius
                    ),
                    "reparameterized_measure_density": _expr_text(
                        self.reparameterized_measure_density
                    ),
                    "reparameterized_measure_differential": _expr_text(
                        self.reparameterized_measure_differential
                    ),
                    "reparameterized_measure_identity": self.reparameterized_measure_identity,
                    "reparameterized_theta_exponent_original": _expr_text(
                        self.reparameterized_theta_exponent_original
                    ),
                    "reparameterized_theta_exponent_log_rewrite": _expr_text(
                        self.reparameterized_theta_exponent_log_rewrite
                    ),
                }
            )
        if include_source_view:
            result["source_substituted"] = {
                label: _expr_text(value)
                for label, value in self.source_substituted().items()
            }
        return result


@dataclass(frozen=True)
class ExactPulseChartSystem:
    """The source-native main, gap and end chart maps."""

    mu: sp.Expr
    log_rp: sp.Expr
    beta: sp.Expr
    source_mu: sp.Expr
    source_mu_recipe: str
    source_substitution: Mapping[sp.Symbol, sp.Expr]
    ordinary_z_independence: OrdinaryZIndependence
    charts: tuple[ExactChartMap, ...]
    radius_value_materialized: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "mu", _positive_mu(self.mu))
        object.__setattr__(self, "log_rp", _exact_expr(self.log_rp, label="log_Rp"))
        object.__setattr__(self, "beta", _exact_expr(self.beta, label="beta"))
        object.__setattr__(self, "source_mu", _positive_mu(self.source_mu))
        if tuple(chart.name for chart in self.charts) != ("main", "gap", "end"):
            raise ValueError("charts must be ordered main, gap, end")

    @property
    def chart_map(self) -> dict[str, ExactChartMap]:
        return {chart.name: chart for chart in self.charts}

    def chart(self, name: ChartName) -> ExactChartMap:
        try:
            return self.chart_map[name]
        except KeyError as exc:
            raise KeyError(f"unknown chart {name!r}") from exc

    def source_expression_view(self) -> dict[str, dict[str, sp.Expr]]:
        """Return all chart expressions after the exact source substitution."""

        return {
            chart.name: chart.source_substituted(self.source_mu)
            for chart in self.charts
        }

    def as_dict(self, *, include_source_view: bool = False) -> dict[str, Any]:
        result: dict[str, Any] = {
            "mu": _expr_text(self.mu),
            "log_Rp": _expr_text(self.log_rp),
            "beta": _expr_text(self.beta),
            "source_mu": _expr_text(self.source_mu),
            "source_mu_recipe": self.source_mu_recipe,
            "source_substitution": {
                str(symbol): _expr_text(expression)
                for symbol, expression in self.source_substitution.items()
            },
            "ordinary_z_independence": self.ordinary_z_independence.as_dict(),
            "radius_value_materialized": self.radius_value_materialized,
            "charts": [
                chart.as_dict(include_source_view=include_source_view)
                for chart in self.charts
            ],
        }
        if include_source_view:
            result["source_expression_view"] = {
                name: {
                    label: _expr_text(expression)
                    for label, expression in values.items()
                }
                for name, values in self.source_expression_view().items()
            }
        return result


def _bounds(lower: sp.Expr, upper: sp.Expr) -> DirectedBounds:
    return DirectedBounds(
        lower=lower,
        upper=upper,
        ordering_assumption="0 < mu <= 1 / 2; chart orientation is source-assumed",
    )


def build_exact_chart_system(
    *,
    mu: Any = None,
    log_rp: Any = None,
    ordinary_z_independent: bool = False,
) -> ExactPulseChartSystem:
    """Build exact main/gap/end chart maps.

    Parameters are symbolic by default.  ``mu`` may be a positive SymPy
    symbol or an exact positive expression.  Passing a Python ``float`` is
    rejected.  The actual source replacement is always available through
    ``source_mu`` and ``source_substitution``; it is not automatically
    substituted into the generic map so nested expressions remain inspectable.
    """

    if mu is None:
        mu_expr = MU_SYMBOL
    else:
        mu_expr = _positive_mu(mu)
    log_rp_expr = (
        LOG_RP_SYMBOL if log_rp is None else _exact_expr(log_rp, label="log_Rp")
    )
    if not isinstance(ordinary_z_independent, bool):
        raise TypeError("ordinary_z_independent must be an explicit bool")
    z_declaration = OrdinaryZIndependence(declared=ordinary_z_independent)
    beta = sp.Rational(1, 2) + mu_expr
    source_mu = source_mu_expression()
    source_substitution = {MU_SYMBOL: source_mu}
    xi, s = XI_SYMBOL, S_SYMBOL
    gap_xi_upper = sp.Integer(13) - sp.Integer(4) * mu_expr
    gap_xi_from_s = sp.Integer(13) + mu_expr * s
    gap_s_from_xi = (xi - sp.Integer(13)) / mu_expr
    end_y = sp.Integer(13) / mu_expr + s
    gap_s_bounds = _bounds(-sp.Integer(2) / mu_expr, -sp.Integer(4))
    theta_end_original = -beta * end_y
    theta_end_rewrite = (
        -(sp.Integer(13) / (sp.Integer(2) * mu_expr) + sp.Integer(13))
        - beta * s
    )

    main_y = xi / mu_expr
    main = ExactChartMap(
        name="main",
        coordinate=xi,
        coordinate_to_y=main_y,
        source_bounds=_bounds(sp.Integer(0), sp.Integer(11)),
        y_bounds=_bounds(sp.Integer(0), sp.Integer(11) / mu_expr),
        log_radius=log_rp_expr + main_y,
        log_radius_bounds=_bounds(
            log_rp_expr, log_rp_expr + sp.Integer(11) / mu_expr
        ),
        physical_measure_density=sp.Integer(1) / mu_expr,
        physical_measure_differential=sp.Symbol("d_xi") / mu_expr,
        physical_measure_identity="dy = dxi / mu",
        beta=beta,
        theta_exponent_original=-beta * main_y,
        theta_exponent_log_rewrite=-beta * main_y,
        ordinary_z_independence=z_declaration,
    )

    gap_y = xi / mu_expr
    gap = ExactChartMap(
        name="gap",
        coordinate=xi,
        coordinate_to_y=gap_y,
        source_bounds=_bounds(sp.Integer(11), gap_xi_upper),
        y_bounds=_bounds(
            sp.Integer(11) / mu_expr, sp.Integer(13) / mu_expr - sp.Integer(4)
        ),
        log_radius=log_rp_expr + gap_y,
        log_radius_bounds=_bounds(
            log_rp_expr + sp.Integer(11) / mu_expr,
            log_rp_expr + sp.Integer(13) / mu_expr - sp.Integer(4),
        ),
        physical_measure_density=sp.Integer(1) / mu_expr,
        physical_measure_differential=sp.Symbol("d_xi") / mu_expr,
        physical_measure_identity="dy = dxi / mu",
        beta=beta,
        theta_exponent_original=-beta * gap_y,
        theta_exponent_log_rewrite=-beta * gap_y,
        ordinary_z_independence=z_declaration,
        xi_from_s=gap_xi_from_s,
        s_from_xi=gap_s_from_xi,
        reparameterized_coordinate=s,
        reparameterized_to_y=end_y,
        reparameterized_bounds=gap_s_bounds,
        reparameterized_log_radius=log_rp_expr + end_y,
        reparameterized_measure_density=sp.Integer(1),
        reparameterized_measure_differential=sp.Symbol("d_s"),
        reparameterized_measure_identity="d_xi / mu = ds",
        reparameterized_theta_exponent_original=theta_end_original,
        reparameterized_theta_exponent_log_rewrite=theta_end_rewrite,
    )

    end = ExactChartMap(
        name="end",
        coordinate=s,
        coordinate_to_y=end_y,
        source_bounds=_bounds(-sp.Integer(4), sp.Integer(0)),
        y_bounds=_bounds(
            sp.Integer(13) / mu_expr - sp.Integer(4), sp.Integer(13) / mu_expr
        ),
        log_radius=log_rp_expr + end_y,
        log_radius_bounds=_bounds(
            log_rp_expr + sp.Integer(13) / mu_expr - sp.Integer(4),
            log_rp_expr + sp.Integer(13) / mu_expr,
        ),
        physical_measure_density=sp.Integer(1),
        physical_measure_differential=sp.Symbol("d_s"),
        physical_measure_identity="dy = ds",
        beta=beta,
        theta_exponent_original=theta_end_original,
        theta_exponent_log_rewrite=theta_end_rewrite,
        ordinary_z_independence=z_declaration,
        xi_from_s=gap_xi_from_s,
        reparameterized_coordinate=s,
        reparameterized_to_y=end_y,
        reparameterized_bounds=_bounds(-sp.Integer(4), sp.Integer(0)),
        reparameterized_log_radius=log_rp_expr + end_y,
        reparameterized_measure_density=sp.Integer(1),
        reparameterized_measure_differential=sp.Symbol("d_s"),
        reparameterized_measure_identity="dy = ds",
        reparameterized_theta_exponent_original=theta_end_original,
        reparameterized_theta_exponent_log_rewrite=theta_end_rewrite,
    )

    return ExactPulseChartSystem(
        mu=mu_expr,
        log_rp=log_rp_expr,
        beta=beta,
        source_mu=source_mu,
        source_mu_recipe=SOURCE_MU_RECIPE,
        source_substitution=source_substitution,
        ordinary_z_independence=z_declaration,
        charts=(main, gap, end),
        radius_value_materialized=False,
    )


# Short aliases for callers that prefer ``chart_maps`` terminology.
build_source_chart_maps = build_exact_chart_system
make_exact_chart_system = build_exact_chart_system


__all__ = [
    "ChartName",
    "DirectedBounds",
    "ExactChartMap",
    "ExactPulseChartSystem",
    "LOG_RP_SYMBOL",
    "MU_SYMBOL",
    "OrdinaryZIndependence",
    "S_SYMBOL",
    "SOURCE_MU_RECIPE",
    "XI_SYMBOL",
    "build_exact_chart_system",
    "build_source_chart_maps",
    "make_exact_chart_system",
    "real_source_mu_expression",
    "source_mu_expression",
]
