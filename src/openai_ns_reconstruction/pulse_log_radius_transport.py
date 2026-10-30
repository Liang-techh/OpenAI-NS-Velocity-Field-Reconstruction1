"""Caller-owned log-radius transport for the pulse charts.

This is a deliberately small, source-free mathematical interface.  It only
uses point-valued caller inputs and ordinary arithmetic.  In particular, it
does not import a pulse source, a coordinate fixture, a quadrature, or a graph
and it never evaluates an exponential.  The physical radius is therefore not
materialised: ``R_Z`` is retained as the formal factor
``exp(ell)*ell_Z``.

For a main/gap chart the caller supplies ``xi`` and its ordinary-Z derivative
``xi_Z`` and the interface evaluates

    ell = logRp + xi/mu
    ell_Z = logRp_Z + xi_Z/mu - xi*mu_Z/mu**2.

For the end chart the caller supplies ``s`` and ``s_Z`` and the interface
evaluates

    ell = logRp + 13/mu + s
    ell_Z = logRp_Z - 13*mu_Z/mu**2 + s_Z.

The gap upper seam has a dedicated structured path.  It retains ``13/mu``
and ``-4`` as separate terms and retains ``-13*mu_Z/mu**2`` as a separate
derivative term.  This matters for a caller's very small ``mu``: finite
precision evaluation of ``13 - 4*mu`` may round to ``13`` before division.

No result from this module is evidence of source identity, ordinary-Z
independence, physical equality, or a common ``Pstar``.  A caller declaration
of one of those facts is recorded as a declaration only, never as proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from decimal import Decimal
from numbers import Number, Rational
import math
from typing import Any, Literal, Mapping


ChartName = Literal["main", "gap", "end"]

SCHEMA = "pulse_log_radius_transport_v1"
FORMAL_RADIUS_Z = "exp(ell)*ell_Z"
_MISSING = object()


class PointInputError(TypeError):
    """Raised when an input is not a determinate point scalar."""


def _comparison(value: Any, other: Any, operator: str, label: str) -> bool:
    """Evaluate one comparison without coercing the caller's numeric type.

    Interval and array libraries commonly return a non-bool comparison or
    raise when a truth value is requested.  Both cases are rejected rather
    than silently selecting a point or an endpoint.
    """

    try:
        raw = (
            value > other
            if operator == ">"
            else value < other
            if operator == "<"
            else value == other
        )
    except Exception as exc:  # pragma: no cover - type-specific path
        raise PointInputError(f"{label} must be a point numeric value") from exc
    if isinstance(raw, bool):
        return raw
    # numpy.bool_ is a determinate scalar; accepting it does not turn arrays
    # or intervals into points because their comparisons are not scalar bools.
    if type(raw).__name__ == "bool_":
        return bool(raw)
    raise PointInputError(f"{label} comparison is indeterminate")


def _require_point(value: Any, label: str) -> Any:
    """Reject bools, containers, intervals, and non-point scalar values.

    The value is returned unchanged.  No float conversion or arbitrary source
    substitution is performed, so Fraction/Decimal/mpmath/custom caller
    arithmetic can retain its own precision and type.
    """

    if value is _MISSING:
        raise PointInputError(f"{label} is required")
    if isinstance(value, bool):
        raise PointInputError(f"{label} must not be bool")
    if isinstance(value, (str, bytes, bytearray, list, tuple, set, frozenset, dict)):
        raise PointInputError(f"{label} must be a point numeric value")
    if getattr(value, "_mpi_", None) is not None:
        raise PointInputError(f"{label} must not be an interval")
    if not isinstance(value, Number) or isinstance(value, complex):
        raise PointInputError(f"{label} must be a finite real number")
    if isinstance(value, Decimal):
        finite = value.is_finite()
    elif isinstance(value, (int, Rational)):
        finite = True
    elif getattr(value, "_mpf_", None) is not None:
        import mpmath
        finite = bool(mpmath.isfinite(value))
    else:
        try:
            finite = math.isfinite(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise PointInputError(f"{label} must be finite") from exc
    if not finite:
        raise PointInputError(f"{label} must be finite")
    type_name = type(value).__name__.lower()
    module_name = type(value).__module__.lower()
    if "interval" in type_name or "interval" in module_name:
        raise PointInputError(f"{label} must not be an interval")
    if "ndarray" in type_name or "array" in type_name:
        raise PointInputError(f"{label} must be a scalar point")
    # A self-equality check catches NaN-like values without converting the
    # input.  If equality itself is uncertain, _comparison rejects it.
    if not _comparison(value, value, "==", label):
        raise PointInputError(f"{label} must be determinate (not NaN/unknown)")
    try:
        value + 0
    except Exception as exc:  # pragma: no cover - custom scalar path
        raise PointInputError(f"{label} must support point arithmetic") from exc
    return value


def _require_positive(value: Any, label: str) -> Any:
    value = _require_point(value, label)
    if not _comparison(value, 0, ">", label):
        raise ValueError(f"{label} must be positive")
    return value


def _divide(numerator: Any, denominator: Any) -> Any:
    """Keep rational caller inputs exact rather than implicitly using float."""
    if isinstance(numerator, Rational) and isinstance(denominator, Rational):
        return Fraction(numerator) / Fraction(denominator)
    return numerator / denominator


def _require_bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise PointInputError(f"{label} must be an explicit bool")
    return value


def _display(value: Any) -> str:
    """Stable, non-coercing display text for JSON/result metadata."""

    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return repr(value)


def _wire(value: Any) -> Any:
    """Convert a point scalar to JSON-safe metadata without changing math."""

    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Fraction):
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    # Decimal, mpmath, and caller-defined numeric objects are intentionally
    # represented by text.  This avoids a hidden float conversion at the
    # artifact boundary.
    return {"type": type(value).__name__, "repr": repr(value)}


@dataclass(frozen=True)
class LogRadiusTransport:
    """One direct log-radius/Z1 evaluation.

    ``ell_terms`` and ``ell_Z_terms`` are retained independently of the
    aggregate values.  Consumers needing a very small ``mu`` can therefore
    carry the large inverse-mu term and finite endpoint term separately.
    """

    chart: ChartName
    endpoint: str
    mu: Any
    logRp: Any
    logRp_Z: Any
    mu_Z: Any
    coordinate_name: Literal["xi", "s"]
    coordinate: Any
    coordinate_Z: Any
    coordinate_expression: str
    ell: Any
    ell_Z: Any
    ell_terms: tuple[Any, ...]
    ell_Z_terms: tuple[Any, ...]
    ell_term_names: tuple[str, ...]
    ell_Z_term_names: tuple[str, ...]
    radius_Z_formula: str = FORMAL_RADIUS_Z
    radius_materialized: bool = False
    source_free: bool = True
    moving_endpoint: bool = False
    claims: Mapping[str, Any] | None = None

    @property
    def ell_z(self) -> Any:
        """Snake-case alias for callers that use ``ell_z``."""

        return self.ell_Z

    @property
    def R_Z(self) -> str:
        """The formal radius derivative, never an evaluated radius value."""

        return self.radius_Z_formula

    @property
    def radius_Z(self) -> str:
        return self.radius_Z_formula

    @property
    def Z1(self) -> Any:
        """Common short name for the ordinary-Z log-radius derivative."""

        return self.ell_Z

    @property
    def logscale(self) -> Any:
        return self.ell

    @property
    def moving_endpoint_proved(self) -> bool:
        return False

    @property
    def source_identity_proved(self) -> bool:
        return False

    @property
    def ordinary_Z_independence_proved(self) -> bool:
        return False

    @property
    def physical_equality_proved(self) -> bool:
        return False

    @property
    def Pstar_proved(self) -> bool:
        return False

    @property
    def log_Rp(self) -> Any:
        return self.logRp

    @property
    def log_Rp_Z(self) -> Any:
        return self.logRp_Z

    def __getitem__(self, key: str) -> Any:
        """Allow lightweight dict-style access without losing exact values."""

        direct = {
            "chart": self.chart,
            "endpoint": self.endpoint,
            "mu": self.mu,
            "logRp": self.logRp,
            "log_Rp": self.logRp,
            "logRp_Z": self.logRp_Z,
            "log_Rp_Z": self.logRp_Z,
            "mu_Z": self.mu_Z,
            "coordinate": self.coordinate,
            "coordinate_Z": self.coordinate_Z,
            "ell": self.ell,
            "ell_Z": self.ell_Z,
            "Z1": self.ell_Z,
            "ell_terms": self.ell_terms,
            "ell_Z_terms": self.ell_Z_terms,
            "R_Z": self.radius_Z_formula,
            "radius_Z_formula": self.radius_Z_formula,
            "radius_materialized": self.radius_materialized,
            "moving_endpoint_proved": self.moving_endpoint_proved,
            "source_identity_proved": self.source_identity_proved,
            "ordinary_Z_independence_proved": self.ordinary_Z_independence_proved,
            "physical_equality_proved": self.physical_equality_proved,
            "Pstar_proved": self.Pstar_proved,
            "claims": self.claims,
        }
        if key == "as_dict":
            return self.as_dict
        try:
            return direct[key]
        except KeyError as exc:
            raise KeyError(key) from exc

    def as_dict(self) -> dict[str, Any]:
        """Return an artifact-safe view while preserving every term."""

        return {
            "schema": SCHEMA,
            "chart": self.chart,
            "endpoint": self.endpoint,
            "mu": _wire(self.mu),
            "logRp": _wire(self.logRp),
            "logRp_Z": _wire(self.logRp_Z),
            "mu_Z": _wire(self.mu_Z),
            "coordinate": {
                "name": self.coordinate_name,
                "value": _wire(self.coordinate),
                "Z1": _wire(self.coordinate_Z),
                "expression": self.coordinate_expression,
                "moving": self.moving_endpoint,
            },
            "ell": _wire(self.ell),
            "ell_Z": _wire(self.ell_Z),
            "ell_terms": [
                {"name": name, "value": _wire(value)}
                for name, value in zip(self.ell_term_names, self.ell_terms)
            ],
            "ell_Z_terms": [
                {"name": name, "value": _wire(value)}
                for name, value in zip(self.ell_Z_term_names, self.ell_Z_terms)
            ],
            "R_Z": self.radius_Z_formula,
            "radius_Z_formula": self.radius_Z_formula,
            "radius_materialized": self.radius_materialized,
            "source_free": self.source_free,
            "moving_endpoint_proved": self.moving_endpoint_proved,
            "source_identity_proved": self.source_identity_proved,
            "ordinary_Z_independence_proved": self.ordinary_Z_independence_proved,
            "physical_equality_proved": self.physical_equality_proved,
            "Pstar_proved": self.Pstar_proved,
            "claims": dict(self.claims or {}),
        }


def _unproved_claims(
    *,
    moving_endpoint: bool,
    source_identity_declared: bool,
    ordinary_Z_independence_declared: bool,
    pstar_supplied: bool,
) -> dict[str, Any]:
    """Metadata gate: caller declarations are not source evidence."""

    return {
        "moving_endpoint": {
            "present": moving_endpoint,
            "proved": False,
            "status": "recorded_not_proved",
        },
        "source_identity": {
            "caller_declared": source_identity_declared,
            "accepted_as_evidence": False,
            "proved": False,
        },
        "ordinary_Z_independence": {
            "caller_declared": ordinary_Z_independence_declared,
            "proved": False,
        },
        "physical_equality": {"proved": False},
        "Pstar": {
            "supplied": pstar_supplied,
            "proved": False,
            "physical_lift": False,
        },
        "rounded_zero_is_source_identity": False,
    }


def _take_alias(
    current: Any, aliases: dict[str, Any], names: tuple[str, ...], label: str
) -> Any:
    """Consume one optional spelling and reject duplicate values."""

    present = [(name, aliases.pop(name)) for name in names if name in aliases]
    if not present:
        return current
    if current is not _MISSING:
        raise TypeError(f"{label} supplied more than once")
    if len(present) != 1:
        raise TypeError(f"{label} supplied with multiple aliases")
    return present[0][1]


def pulse_log_radius_transport(
    chart: str = _MISSING,  # type: ignore[assignment]
    mu: Any = _MISSING,
    logRp: Any = _MISSING,
    logRp_Z: Any = _MISSING,
    mu_Z: Any = _MISSING,
    xi: Any = _MISSING,
    xi_Z: Any = _MISSING,
    s: Any = _MISSING,
    s_Z: Any = _MISSING,
    *,
    endpoint: str | None = None,
    Pstar: Any = _MISSING,
    source_identity_declared: bool = False,
    ordinary_Z_independence_declared: bool = False,
    **aliases: Any,
) -> LogRadiusTransport:
    """Compute a pulse chart's ``ell`` and ordinary-Z ``ell_Z``.

    The function is non-reentrant and caller-owned: all scale and endpoint
    values are supplied on each call, no module state is changed, and no
    source object is looked up.  ``chart`` is ``main``, ``gap``, or ``end``.
    Main/gap use ``xi``/``xi_Z``; end uses ``s``/``s_Z``.  Set
    ``endpoint='upper'`` for the structured gap upper seam, which uses
    ``xi = 13 - 4*mu`` only as a recorded relation and computes directly from
    ``13/mu - 4``.
    """

    # Friendly spellings do not alter the mathematical interface and are
    # useful when a caller's existing record uses snake case.
    chart = _take_alias(chart, aliases, ("chart_name",), "chart")
    logRp = _take_alias(logRp, aliases, ("log_rp", "log_Rp"), "logRp")
    logRp_Z = _take_alias(
        logRp_Z, aliases, ("log_rp_Z", "log_rp_z", "log_Rp_Z"), "logRp_Z"
    )
    mu_Z = _take_alias(mu_Z, aliases, ("mu_z",), "mu_Z")
    xi_Z = _take_alias(xi_Z, aliases, ("xi_z",), "xi_Z")
    s_Z = _take_alias(s_Z, aliases, ("s_z",), "s_Z")
    if aliases:
        unknown = ", ".join(sorted(aliases))
        raise TypeError(f"unexpected input(s): {unknown}")

    if not isinstance(chart, str):
        raise TypeError("chart must be 'main', 'gap', or 'end'")
    chart_key = chart.lower()
    if chart_key in {"gap-upper", "gap_upper", "gapupper"}:
        chart_key = "gap"
        endpoint = "upper"
    if chart_key not in {"main", "gap", "end"}:
        raise ValueError(f"unknown chart {chart!r}")
    if endpoint is not None and endpoint not in {"generic", "upper"}:
        raise ValueError("endpoint must be 'upper' or 'generic'")
    if endpoint == "upper" and chart_key != "gap":
        raise ValueError("the upper endpoint form is available only for gap")

    mu = _require_positive(mu, "mu")
    logRp = _require_point(logRp, "logRp")
    logRp_Z = _require_point(logRp_Z, "logRp_Z")
    mu_Z = _require_point(mu_Z, "mu_Z")
    source_identity_declared = _require_bool(
        source_identity_declared, "source_identity_declared"
    )
    ordinary_Z_independence_declared = _require_bool(
        ordinary_Z_independence_declared,
        "ordinary_Z_independence_declared",
    )
    if Pstar is not _MISSING:
        _require_point(Pstar, "Pstar")

    mu_square = mu * mu
    is_upper = endpoint == "upper"

    if chart_key in {"main", "gap"} and not is_upper:
        if s is not _MISSING or s_Z is not _MISSING:
            raise TypeError("main/gap charts use xi and xi_Z, not s and s_Z")
        xi = _require_point(xi, "xi")
        xi_Z = _require_point(xi_Z, "xi_Z")
        xi_over_mu = _divide(xi, mu)
        xi_Z_over_mu = _divide(xi_Z, mu)
        xi_mu_Z_over_mu_square = _divide(xi * mu_Z, mu_square)
        ell_terms = (logRp, xi_over_mu)
        ell_Z_terms = (logRp_Z, xi_Z_over_mu, -xi_mu_Z_over_mu_square)
        ell_names = ("logRp", "xi_over_mu")
        ell_Z_names = ("logRp_Z", "xi_Z_over_mu", "minus_xi_mu_Z_over_mu2")
        ell = logRp + xi_over_mu
        ell_Z = logRp_Z + xi_Z_over_mu - xi_mu_Z_over_mu_square
        coordinate_name: Literal["xi", "s"] = "xi"
        coordinate = xi
        coordinate_Z = xi_Z
        coordinate_expression = "xi"
        endpoint_name = "generic"
    elif is_upper:
        if xi is not _MISSING or xi_Z is not _MISSING:
            raise TypeError("gap upper uses its structured endpoint; omit xi/xi_Z")
        if s is not _MISSING or s_Z is not _MISSING:
            raise TypeError("gap upper does not use s/s_Z")
        # Keep the finite seam term separate.  Do not first form 13 - 4*mu.
        inverse_mu_13 = _divide(13, mu)
        minus_four = -4
        minus_13_mu_Z_over_mu_square = -_divide(13 * mu_Z, mu_square)
        ell_terms = (logRp, inverse_mu_13, minus_four)
        ell_Z_terms = (logRp_Z, minus_13_mu_Z_over_mu_square)
        ell_names = ("logRp", "13_over_mu", "minus_four")
        ell_Z_names = ("logRp_Z", "minus_13_mu_Z_over_mu2")
        ell = logRp + inverse_mu_13 + minus_four
        ell_Z = logRp_Z + minus_13_mu_Z_over_mu_square
        coordinate_name = "xi"
        # This is a point record for ordinary-size callers.  The actual
        # calculation above intentionally does not depend on this subtraction.
        coordinate = 13 - 4 * mu
        coordinate_Z = -4 * mu_Z
        coordinate_expression = "13 - 4*mu (recorded; direct scale is 13/mu - 4)"
        endpoint_name = "upper"
    else:
        if xi is not _MISSING or xi_Z is not _MISSING:
            raise TypeError("end chart uses s and s_Z, not xi and xi_Z")
        s = _require_point(s, "s")
        s_Z = _require_point(s_Z, "s_Z")
        inverse_mu_13 = _divide(13, mu)
        minus_13_mu_Z_over_mu_square = -_divide(13 * mu_Z, mu_square)
        ell_terms = (logRp, inverse_mu_13, s)
        ell_Z_terms = (logRp_Z, minus_13_mu_Z_over_mu_square, s_Z)
        ell_names = ("logRp", "13_over_mu", "s")
        ell_Z_names = ("logRp_Z", "minus_13_mu_Z_over_mu2", "s_Z")
        ell = logRp + inverse_mu_13 + s
        ell_Z = logRp_Z + minus_13_mu_Z_over_mu_square + s_Z
        coordinate_name = "s"
        coordinate = s
        coordinate_Z = s_Z
        coordinate_expression = "s"
        endpoint_name = "generic"

    # _comparison intentionally has only order/equality operators; use a
    # scalar equality check directly so no interval truth is inferred here.
    moving_endpoint = None
    try:
        endpoint_zero = coordinate_Z == 0
        if not isinstance(endpoint_zero, bool) and type(endpoint_zero).__name__ != "bool_":
            raise PointInputError("endpoint_Z comparison is indeterminate")
        moving_endpoint = not bool(endpoint_zero)
    except Exception as exc:
        raise PointInputError("endpoint_Z must be a determinate point") from exc

    claims = _unproved_claims(
        moving_endpoint=moving_endpoint,
        source_identity_declared=source_identity_declared,
        ordinary_Z_independence_declared=ordinary_Z_independence_declared,
        pstar_supplied=Pstar is not _MISSING,
    )
    return LogRadiusTransport(
        chart=chart_key,  # type: ignore[arg-type]
        endpoint=endpoint_name,
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        coordinate_name=coordinate_name,
        coordinate=coordinate,
        coordinate_Z=coordinate_Z,
        coordinate_expression=coordinate_expression,
        ell=ell,
        ell_Z=ell_Z,
        ell_terms=ell_terms,
        ell_Z_terms=ell_Z_terms,
        ell_term_names=ell_names,
        ell_Z_term_names=ell_Z_names,
        claims=claims,
    )


def main_gap_log_radius_transport(
    *, mu: Any, logRp: Any, logRp_Z: Any, mu_Z: Any, xi: Any, xi_Z: Any, chart: ChartName = "main"
) -> LogRadiusTransport:
    """Explicit main/gap wrapper for the ``xi/mu`` chart formula."""

    if chart not in {"main", "gap"}:
        raise ValueError("chart must be 'main' or 'gap'")
    return pulse_log_radius_transport(
        chart=chart,
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        xi=xi,
        xi_Z=xi_Z,
    )


def gap_upper_log_radius_transport(
    *, mu: Any, logRp: Any, logRp_Z: Any, mu_Z: Any
) -> LogRadiusTransport:
    """Direct structured transport at ``xi=13-4*mu``."""

    return pulse_log_radius_transport(
        chart="gap",
        endpoint="upper",
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
    )


def end_log_radius_transport(
    *, mu: Any, logRp: Any, logRp_Z: Any, mu_Z: Any, s: Any, s_Z: Any
) -> LogRadiusTransport:
    """Explicit end-chart wrapper for ``13/mu+s``."""

    return pulse_log_radius_transport(
        chart="end",
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        s=s,
        s_Z=s_Z,
    )


# Descriptive aliases for callers that prefer a verb or a snake-case name.
compute_log_radius_transport = pulse_log_radius_transport
compute_pulse_log_radius_transport = pulse_log_radius_transport
log_radius_transport = pulse_log_radius_transport
transport_log_radius = pulse_log_radius_transport
evaluate_log_radius_transport = pulse_log_radius_transport


def _synthetic_fixture() -> dict[str, Any]:
    """Run exact, source-free Fraction checks for this independent interface."""

    # Every nonzero derivative is intentional: this catches accidental
    # ``mu_Z=0`` or ``logRp_Z=0`` shortcuts and keeps the fixture independent
    # of the archived coordinate fixture and all quadrature/graph code.
    mu = Fraction(5, 7)
    logRp = Fraction(7, 11)
    mu_Z = Fraction(2, 9)
    logRp_Z = Fraction(-3, 13)
    gap_upper = gap_upper_log_radius_transport(
        mu=mu, logRp=logRp, logRp_Z=logRp_Z, mu_Z=mu_Z
    )
    end_upper = end_log_radius_transport(
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        s=Fraction(-4),
        s_Z=Fraction(0),
    )
    main_11 = main_gap_log_radius_transport(
        chart="main",
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        xi=Fraction(11),
        xi_Z=Fraction(0),
    )
    gap_11 = main_gap_log_radius_transport(
        chart="gap",
        mu=mu,
        logRp=logRp,
        logRp_Z=logRp_Z,
        mu_Z=mu_Z,
        xi=Fraction(11),
        xi_Z=Fraction(0),
    )

    xi = Fraction(13) - 4 * mu
    xi_Z = -4 * mu_Z
    # A finite difference of the moving endpoint itself is exact for the
    # affine relation xi(mu): no approximation or tolerance is involved.
    step = Fraction(1, 17)
    mu_next = mu + step * mu_Z
    xi_next = Fraction(13) - 4 * mu_next
    endpoint_difference = (xi_next - xi) / step
    quotient_difference_form = (xi_Z * mu - xi * mu_Z) / (mu * mu)
    checks = {
        "gap_upper_end_logscale_exact": gap_upper.ell == end_upper.ell,
        "gap_upper_end_Z1_exact": gap_upper.ell_Z == end_upper.ell_Z,
        "gap_upper_closed_form_exact": gap_upper.ell
        == logRp + Fraction(13) / mu - Fraction(4),
        "gap_upper_Z1_closed_form_exact": gap_upper.ell_Z
        == logRp_Z - Fraction(13) * mu_Z / (mu * mu),
        "main_gap_xi11_logscale_exact": main_11.ell == gap_11.ell,
        "main_gap_xi11_Z1_exact": main_11.ell_Z == gap_11.ell_Z,
        "moving_endpoint_difference_exact": endpoint_difference == xi_Z,
        "moving_endpoint_quotient_difference_exact": quotient_difference_form
        == xi_Z / mu - xi * mu_Z / (mu * mu),
        "nonzero_mu_Z": mu_Z != 0,
        "nonzero_logRp_Z": logRp_Z != 0,
        "R_Z_formula_preserved": gap_upper.R_Z == FORMAL_RADIUS_Z
        and end_upper.R_Z == FORMAL_RADIUS_Z,
        "radius_not_materialized": not gap_upper.radius_materialized
        and not end_upper.radius_materialized,
        "source_free": gap_upper.source_free and end_upper.source_free,
        "unproved_source_claims": not gap_upper.claims["source_identity"]["proved"],
        "unproved_Z_independence": not gap_upper.claims["ordinary_Z_independence"]["proved"],
        "unproved_physical_equality": not gap_upper.claims["physical_equality"]["proved"],
        "unproved_Pstar": not gap_upper.claims["Pstar"]["proved"],
    }
    integer_case = main_gap_log_radius_transport(chart="gap", mu=3, logRp=1,
        logRp_Z=2, mu_Z=1, xi=11, xi_Z=0)
    checks["integer_division_remains_exact"] = isinstance(integer_case.ell, Fraction) and integer_case.ell == Fraction(14,3)
    try:
        pulse_log_radius_transport(chart="end", mu=mu, logRp=logRp,
            logRp_Z=logRp_Z, mu_Z=mu_Z, s=Fraction(-4))
    except (PointInputError, TypeError):
        checks["missing_end_coordinate_derivative_rejected"] = True
    else:
        checks["missing_end_coordinate_derivative_rejected"] = False
    try:
        gap_upper_log_radius_transport(mu=float("inf"), logRp=logRp,
            logRp_Z=logRp_Z, mu_Z=mu_Z)
    except PointInputError:
        checks["nonfinite_mu_rejected"] = True
    else:
        checks["nonfinite_mu_rejected"] = False
    import mpmath
    try:
        gap_upper_log_radius_transport(mu=mpmath.iv.mpf([1,2]), logRp=logRp,
            logRp_Z=logRp_Z, mu_Z=mu_Z)
    except PointInputError:
        checks["interval_input_rejected"] = True
    else:
        checks["interval_input_rejected"] = False
    tiny_context = mpmath.mp.clone()
    tiny_context.dps = 50
    tiny = tiny_context.mpf("1e-100")
    tiny_case = gap_upper_log_radius_transport(mu=tiny, logRp=tiny_context.mpf(1),
        logRp_Z=tiny_context.mpf(0), mu_Z=tiny_context.mpf(0))
    checks["tiny_mu_keeps_finite_endpoint_term_separately"] = tiny_case.coordinate == 13 and tiny_case.ell_terms[2] == -4

    try:
        gap_upper_log_radius_transport(
            mu=Fraction(0), logRp=logRp, logRp_Z=logRp_Z, mu_Z=mu_Z
        )
    except (PointInputError, ValueError, ZeroDivisionError):
        checks["reject_mu_zero"] = True
    else:  # pragma: no cover - a regression guard
        checks["reject_mu_zero"] = False
    return {
        "schema": "pulse_log_radius_transport_synthetic_fixture_v1",
        "status": "source_free_synthetic_exact_checks",
        "source_free_synthetic": True,
        "not_actual_source": True,
        "different_derivative_interface_from_coordinate_fixture": True,
        "old_coordinate_fixture_imported": False,
        "quadrature_or_graph_called": False,
        "inputs": {
            "mu": _wire(mu),
            "logRp": _wire(logRp),
            "mu_Z": _wire(mu_Z),
            "logRp_Z": _wire(logRp_Z),
            "moving_endpoint": {
                "xi": _wire(xi),
                "xi_Z": _wire(xi_Z),
                "finite_difference_step": _wire(step),
            },
        },
        "checks": checks,
        "all_passed": all(checks.values()),
        "gap_upper": gap_upper.as_dict(),
        "end_upper": end_upper.as_dict(),
        "main_xi11": main_11.as_dict(),
        "gap_xi11": gap_11.as_dict(),
        "unproved": {
            "moving_endpoint_source_identity": True,
            "source_identity": True,
            "ordinary_Z_independence": True,
            "physical_equality": True,
            "Pstar": True,
        },
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Public one-call fixture entry used by the one-shot runner."""

    return _synthetic_fixture()


def run_fixture() -> dict[str, Any]:
    """Compatibility spelling for fixture harnesses; still one local call."""

    return verify_synthetic_fixture()


__all__ = [
    "ChartName",
    "FORMAL_RADIUS_Z",
    "LogRadiusTransport",
    "PointInputError",
    "SCHEMA",
    "compute_log_radius_transport",
    "compute_pulse_log_radius_transport",
    "end_log_radius_transport",
    "gap_upper_log_radius_transport",
    "log_radius_transport",
    "main_gap_log_radius_transport",
    "pulse_log_radius_transport",
    "run_fixture",
    "transport_log_radius",
    "evaluate_log_radius_transport",
    "verify_synthetic_fixture",
]
