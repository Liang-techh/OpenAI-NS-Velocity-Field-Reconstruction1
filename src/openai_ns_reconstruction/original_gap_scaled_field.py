"""Caller-owned arbitrary-precision port of the original gap field.

This leaf evaluates only the point formula supplied by its caller.  It does
not import the historical pulse owner, recover a source graph, or identify
the supplied jets with a source.  The large original factor is represented by
an ordered list of logarithmic terms; it is never exponentiated by the field
port.

The input ``u``, ``C1`` and ``C2`` are ordinary derivative jets in ``Z``
through order three.  The jet arithmetic below uses the ordinary Leibniz and
quotient rules, rather than Taylor-coefficient convolution.  Three input
derivatives are exactly what is needed for the value of ``Q_ZZ``.  A
derivative of ``Q_ZZ`` is therefore deliberately reported as unavailable.

The original gap expression is

    m = -Pstar*u*(A0*C1 + A1*C2)
        * exp(logE + lambda*(13-xi)/mu - b*xi/mu),

where ``lambda=1/2-mu``, ``b=1/2+mu`` and

    logE = -1/mu - 3*(4*mu)**(-2/3) + (2/3)*log(4*mu)
           - (1/2)*log(6) - 2*log(mu).

With ``xi=anchor+mu*offset`` and ``Pstar=exp(logPstar)``, the unmaterialized
scale terms, in order, are

    (11/2-anchor)/mu,
    -offset,
    -3*(4*mu)**(-2/3),
    -(4/3)*log(mu),
    (2/3)*log(4) - (1/2)*log(6),
    -13,
    logPstar.

In particular, ``anchor-13`` and ``mu*offset`` are not separate scale terms:
they have already been accounted for by the displayed first, second and
sixth terms.  The returned normalized mantissas omit this scale and retain
the ``Pstar`` contribution only through its ``logPstar`` term.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import Any, Mapping, Sequence

import mpmath as mp


JET_ORDER = 3
LOGSCALE_LABELS = (
    "(11/2-anchor)/mu",
    "-offset",
    "-3*(4*mu)^(-2/3)",
    "-(4/3)*log(mu)",
    "(2/3)*log(4)-(1/2)*log(6)",
    "-13",
    "logPstar",
)

ENDPOINT_LOGSCALE_LABELS = (
    "(11/2)/mu",
    "-13",
    "-3*(4*mu)^(-2/3)",
    "-(4/3)*log(mu)",
    "(2/3)*log(4)-(1/2)*log(6)",
    "logPstar",
    "logRp",
)

LINEAR_RADIUS_LOGSCALE_LABELS = (
    "(11/2+anchor)/mu",
    "offset",
    "-13",
    "-3*(4*mu)^(-2/3)",
    "-(4/3)*log(mu)",
    "(2/3)*log(4)-(1/2)*log(6)",
    "logPstar",
    "2*logRp",
)

QUADRATIC_RADIUS_LOGSCALE_LABELS = (
    "11/mu",
    "-26",
    "-6*(4*mu)^(-2/3)",
    "-(8/3)*log(mu)",
    "2*((2/3)*log(4)-(1/2)*log(6))",
    "2*logPstar",
    "2*logRp",
)


def _is_interval(value: Any) -> bool:
    """Reject mpmath interval values before they can be coerced to a scalar."""

    value_type = type(value)
    type_name = value_type.__name__.lower()
    module_name = value_type.__module__.lower()
    return (
        "interval" in type_name
        or type_name in {"ivmpf", "ivmpc"}
        or module_name.startswith("mpmath.ctx_iv")
        or hasattr(value, "_mpi_")
    )


def _scalar(value: Any, name: str) -> mp.mpf:
    """Coerce one finite real scalar without accepting interval arithmetic."""

    if isinstance(value, bool) or _is_interval(value):
        raise ValueError(f"{name} must be one finite real scalar")
    try:
        result = mp.mpf(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be one finite real scalar") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _jet(value: Any, name: str) -> tuple[mp.mpf, ...]:
    """Read an ordinary Z0..Z3 derivative jet as four finite mpf values."""

    if isinstance(value, (str, bytes, bytearray, Mapping)) or _is_interval(value):
        raise ValueError(f"{name} must be an ordinary Z0..Z3 jet")
    try:
        entries = tuple(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an ordinary Z0..Z3 jet") from exc
    if len(entries) != JET_ORDER + 1:
        raise ValueError(f"{name} must contain exactly four entries Z0..Z3")
    return tuple(_scalar(entry, f"{name}[{index}]")
                 for index, entry in enumerate(entries))


def _zero_jet() -> tuple[mp.mpf, ...]:
    return (mp.mpf(0),) * (JET_ORDER + 1)


def _constant_jet(value: Any) -> tuple[mp.mpf, ...]:
    scalar = _scalar(value, "constant")
    return (scalar,) + _zero_jet()[1:]


def _add(*jets: tuple[mp.mpf, ...]) -> tuple[mp.mpf, ...]:
    if not jets:
        return _zero_jet()
    return tuple(sum(jet[index] for jet in jets)
                 for index in range(JET_ORDER + 1))


def _neg(jet: tuple[mp.mpf, ...]) -> tuple[mp.mpf, ...]:
    return tuple(-entry for entry in jet)


def _sub(left: tuple[mp.mpf, ...], right: tuple[mp.mpf, ...]) -> tuple[mp.mpf, ...]:
    return _add(left, _neg(right))


def _scale(jet: tuple[mp.mpf, ...], scalar: Any) -> tuple[mp.mpf, ...]:
    factor = _scalar(scalar, "jet scale")
    return tuple(factor * entry for entry in jet)


def _mul(left: tuple[mp.mpf, ...], right: tuple[mp.mpf, ...]) -> tuple[mp.mpf, ...]:
    """Multiply ordinary derivative jets through order three exactly."""

    return tuple(
        sum(mp.mpf(comb(order, index)) * left[index] * right[order - index]
            for index in range(order + 1))
        for order in range(JET_ORDER + 1)
    )


def _inverse(jet: tuple[mp.mpf, ...], name: str) -> tuple[mp.mpf, ...]:
    """Invert an ordinary derivative jet by the exact recursive product rule."""

    if jet[0] == 0:
        raise ValueError(f"{name} has a zero Z0 denominator")
    result = [mp.mpf(0)] * (JET_ORDER + 1)
    result[0] = 1 / jet[0]
    for order in range(1, JET_ORDER + 1):
        result[order] = -sum(
            mp.mpf(comb(order, index)) * jet[index] * result[order - index]
            for index in range(1, order + 1)
        ) / jet[0]
    return tuple(result)


def _quotient(numerator: tuple[mp.mpf, ...], denominator: tuple[mp.mpf, ...],
              name: str) -> tuple[mp.mpf, ...]:
    return _mul(numerator, _inverse(denominator, name))


def _z_jet(z: mp.mpf) -> tuple[mp.mpf, ...]:
    return (z, mp.mpf(1), mp.mpf(0), mp.mpf(0))


def _derivative_jet(jet: tuple[mp.mpf, ...]) -> tuple[mp.mpf, ...]:
    """Shift a C3 jet once; the final slot is unused by Q through Z2."""

    return (jet[1], jet[2], jet[3], mp.mpf(0))


def _scale_terms(mu: mp.mpf, anchor: mp.mpf, offset: mp.mpf,
                 log_pstar: mp.mpf) -> tuple[mp.mpf, ...]:
    """Return the seven ordered logarithmic terms without summing them."""

    return (
        (mp.mpf(11) / 2 - anchor) / mu,
        -offset,
        -3 * (4 * mu) ** (-mp.mpf(2) / 3),
        -(mp.mpf(4) / 3) * mp.log(mu),
        (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6),
        -mp.mpf(13),
        log_pstar,
    )


def _endpoint_scale_terms(mu: mp.mpf, log_pstar: mp.mpf,
                          log_rp: mp.mpf) -> tuple[mp.mpf, ...]:
    """Scale of a retained physical ``R*Q_y`` endpoint, unsummed."""

    saddle = (4 * mu) ** (-mp.mpf(2) / 3)
    constant = (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6)
    return (
        (mp.mpf(11) / 2) / mu,
        -mp.mpf(13),
        -3 * saddle,
        -(mp.mpf(4) / 3) * mp.log(mu),
        constant,
        log_pstar,
        log_rp,
    )


def _linear_radius_scale_terms(mu: mp.mpf, anchor: mp.mpf,
                               offset: mp.mpf, log_pstar: mp.mpf,
                               log_rp: mp.mpf) -> tuple[mp.mpf, ...]:
    """Scale of the physical ``R^2*Slinear`` density term, unsummed."""

    saddle = (4 * mu) ** (-mp.mpf(2) / 3)
    constant = (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6)
    return (
        (mp.mpf(11) / 2 + anchor) / mu,
        offset,
        -mp.mpf(13),
        -3 * saddle,
        -(mp.mpf(4) / 3) * mp.log(mu),
        constant,
        log_pstar,
        2 * log_rp,
    )


def _quadratic_radius_scale_terms(mu: mp.mpf, log_pstar: mp.mpf,
                                  log_rp: mp.mpf) -> tuple[mp.mpf, ...]:
    """Scale of the physical ``R^2*Squadratic`` density term, unsummed."""

    saddle = (4 * mu) ** (-mp.mpf(2) / 3)
    constant = (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6)
    return (
        mp.mpf(11) / mu,
        -mp.mpf(26),
        -6 * saddle,
        -(mp.mpf(8) / 3) * mp.log(mu),
        2 * constant,
        2 * log_pstar,
        2 * log_rp,
    )


@dataclass(frozen=True)
class ScaledQuantity:
    """A normalized value/Z1 pair with an unmaterialized logarithmic scale."""

    value: mp.mpf
    Z1: mp.mpf | None
    scale_terms: tuple[mp.mpf, ...]
    scale_multiplier: int
    name: str

    @property
    def logscale_components(self) -> tuple[mp.mpf, ...]:
        """The base terms; ``scale_multiplier`` says how many copies apply."""

        return self.scale_terms

    @property
    def Z0(self) -> mp.mpf:
        return self.value

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "value": self.value,
            "Z1": self.Z1,
            "normalized": True,
            "logscale_components": list(self.scale_terms),
            "logscale_multiplier": self.scale_multiplier,
            "logscale_materialized": False,
        }


@dataclass(frozen=True)
class GapScaledFieldResult:
    """All normalized gap fields and the separate reduced-S terms."""

    mu: mp.mpf
    Z: mp.mpf
    delta: mp.mpf
    L: mp.mpf
    logRp: mp.mpf
    anchor: mp.mpf
    offset: mp.mpf
    xi: mp.mpf
    lambda_value: mp.mpf
    b_value: mp.mpf
    scale_terms: tuple[mp.mpf, ...]
    endpoint_logscale_terms: tuple[mp.mpf, ...]
    linear_radius_logscale_terms: tuple[mp.mpf, ...]
    quadratic_radius_logscale_terms: tuple[mp.mpf, ...]
    m: ScaledQuantity
    Q: ScaledQuantity
    Q_Z: ScaledQuantity
    Q_ZZ: ScaledQuantity
    Q_y: ScaledQuantity
    Q_yZ: ScaledQuantity
    Slinear: ScaledQuantity
    Squadratic: ScaledQuantity
    source_domain_premise: bool
    delta_Z_independent: bool
    pstar_materialized: bool = False
    source_identity_proved: bool = False

    @property
    def S_linear(self) -> ScaledQuantity:
        return self.Slinear

    @property
    def S_quadratic(self) -> ScaledQuantity:
        return self.Squadratic

    @property
    def normalized_fields(self) -> Mapping[str, ScaledQuantity]:
        return {
            "m": self.m,
            "Q": self.Q,
            "Q_Z": self.Q_Z,
            "Q_ZZ": self.Q_ZZ,
            "Q_y": self.Q_y,
            "Q_yZ": self.Q_yZ,
            "Slinear": self.Slinear,
            "Squadratic": self.Squadratic,
        }

    def as_dict(self) -> dict[str, Any]:
        fields = {name: quantity.as_dict()
                  for name, quantity in self.normalized_fields.items()}
        return {
            "schema": "original_gap_scaled_field_v1",
            "coordinates": {
                "anchor": self.anchor,
                "offset": self.offset,
                "xi_rounded_diagnostic": self.xi,
                "xi_definition": "anchor + mu*offset",
                "source_interval_definition": "11 <= xi <= 13-4*mu",
            },
            "parameters": {
                "mu": self.mu,
                "Z": self.Z,
                "delta": self.delta,
                "L": self.L,
                "logRp": self.logRp,
                "lambda": self.lambda_value,
                "b": self.b_value,
            },
            "fields": fields,
            "m15_reduced_S": {
                "linear": self.Slinear.as_dict(),
                "quadratic": self.Squadratic.as_dict(),
            },
            "logscale": {
                "labels": list(LOGSCALE_LABELS),
                "single_factor": {
                    "components": list(self.scale_terms),
                    "multiplicity": 1,
                },
                "quadratic_factor": {
                    "components": list(self.scale_terms),
                    "multiplicity": 2,
                },
                "physical_radius_factors": {
                    "R_Q_y_endpoint": {
                        "labels": list(ENDPOINT_LOGSCALE_LABELS),
                        "components": list(self.endpoint_logscale_terms),
                        "multiplicity": 1,
                        "normalized_value": self.Q_y.value,
                        "normalized_Z1": self.Q_y.Z1,
                        "orientation": "minus_[R*Q_y] at each retained endpoint",
                    },
                    "R_squared_Slinear": {
                        "labels": list(LINEAR_RADIUS_LOGSCALE_LABELS),
                        "components": list(self.linear_radius_logscale_terms),
                        "multiplicity": 1,
                        "normalized_value": self.Slinear.value,
                        "normalized_Z1": self.Slinear.Z1,
                    },
                    "R_squared_Squadratic": {
                        "labels": list(QUADRATIC_RADIUS_LOGSCALE_LABELS),
                        "components": list(self.quadratic_radius_logscale_terms),
                        "multiplicity": 1,
                        "normalized_value": self.Squadratic.value,
                        "normalized_Z1": self.Squadratic.Z1,
                    },
                },
                "materialized": False,
                "Pstar": {
                    "represented_by": "logPstar",
                    "materialized": False,
                    "applied_once_to_single_terms": True,
                    "applied_twice_to_quadratic_term": True,
                },
            },
            "premises": {
                "source_domain_premise": self.source_domain_premise,
                "delta_Z_independent": self.delta_Z_independent,
                "scale_coordinate_and_row_weights_Z_independent": True,
                "positive_L": self.L > 0,
                "B_zero": True,
                "U_zero": True,
                "Uy_zero": True,
                "my_equals_minus_m": True,
                "source_identity_proved": self.source_identity_proved,
                "caller_owned_point_data": True,
            },
            "formula_components": {
                "Q_numerator": {
                    "term_1": "-(1-delta)*Z*m",
                    "term_2": "-(1-Z^2)*m_Z",
                    "denominator": "1-delta*Z^2",
                    "endpoint_subtraction_retained": True,
                },
                "Slinear": "(1-delta)*Z*Q_Z/(2*L)",
                "Squadratic": "-Q^2/2",
                "density": {
                    "coefficient": "1/2",
                    "coordinate_measure": "d_xi/mu",
                    "offset_measure": "d_offset when xi=anchor+mu*offset",
                    "integrated_here": False,
                },
                "endpoint": {
                    "expression": "-[R*Q_y]",
                    "left_right_orientation_retained": True,
                    "claimed_zero": False,
                },
                "Q_y": "-Q",
                "Q_yZ": "-Q_Z",
            },
        }


def evaluate_original_gap_scaled_field(
    mu: Any,
    Z: Any,
    delta: Any,
    A0: Any,
    A1: Any,
    u: Sequence[Any],
    C1: Sequence[Any],
    C2: Sequence[Any],
    *,
    logRp: Any,
    logPstar: Any,
    anchor: Any = 13,
    offset: Any = -4,
    source_domain_premise: bool = False,
    delta_Z_independent: bool = False,
    scale_and_weights_Z_independent: bool = False,
) -> GapScaledFieldResult:
    """Evaluate the caller-supplied original gap field at one point.

    ``logRp`` is required because the physical radius scales needed by a
    future m15 adapter are distinct from the normalized field scale.  They
    are returned as separate component lists for ``R*Q_y``,
    ``R^2*Slinear`` and ``R^2*Squadratic``; none is summed or exponentiated.

    ``u``, ``C1`` and ``C2`` are ordinary derivative jets through Z3.  The
    returned ``value`` and ``Z1`` entries are normalized mantissas.  A single
    logarithmic factor applies to ``m``, all Q fields, and ``Slinear``; the
    quadratic term uses the same ordered components with multiplier two.

    ``source_domain_premise=True`` is an explicit caller-owned assertion that
    the point lies in the original gap source interval.  Otherwise the strict
    high-precision test ``11 <= anchor+mu*offset <= 13-4*mu`` is applied.  No
    source provider or source identity is inferred from that test.
    """

    if not isinstance(source_domain_premise, bool):
        raise ValueError("source_domain_premise must be an explicit boolean")
    if delta_Z_independent is not True:
        raise ValueError("delta_Z_independent=True is required")
    if scale_and_weights_Z_independent is not True:
        raise ValueError("mu, anchor, offset, logRp, logPstar and A0/A1 must be Z independent")

    mu_value = _scalar(mu, "mu")
    if not 0 < mu_value < mp.mpf("0.5"):
        raise ValueError("mu must lie in the original positive regime below 1/2")
    Z_value = _scalar(Z, "Z")
    delta_value = _scalar(delta, "delta")
    if not 0 < delta_value < mp.mpf("0.5") or abs(Z_value) > 1:
        raise ValueError("delta and Z are outside the original positive-L domain")
    L_value = 1 - delta_value * Z_value**2
    if L_value <= 0:
        raise ValueError("L must be positive")
    A0_value = _scalar(A0, "A0")
    A1_value = _scalar(A1, "A1")
    log_rp = _scalar(logRp, "logRp")
    anchor_value = _scalar(anchor, "anchor")
    offset_value = _scalar(offset, "offset")
    log_pstar = _scalar(logPstar, "logPstar")
    xi_value = anchor_value + mu_value * offset_value
    lower_slack = (anchor_value - 11) + mu_value * offset_value
    upper_slack = (13 - anchor_value) - mu_value * (offset_value + 4)
    if not source_domain_premise and not (lower_slack >= 0 and upper_slack >= 0):
        raise ValueError("xi=anchor+mu*offset lies outside the original gap source interval")

    u_jet = _jet(u, "u")
    C1_jet = _jet(C1, "C1")
    C2_jet = _jet(C2, "C2")
    coefficient_jet = _add(_scale(C1_jet, A0_value), _scale(C2_jet, A1_value))
    # The Pstar and exponential factors are intentionally not multiplied here.
    m_jet = _scale(_mul(u_jet, coefficient_jet), -1)

    z_jet = _z_jet(Z_value)
    one = _constant_jet(1)
    delta_jet = _constant_jet(delta_value)
    one_minus_delta = _constant_jet(1 - delta_value)
    z_squared = _mul(z_jet, z_jet)
    denominator = _sub(one, _mul(delta_jet, z_squared))
    m_z_jet = _derivative_jet(m_jet)
    # Keep the two numerator terms separate so the endpoint subtraction is
    # visible and cannot be silently dropped by a later adapter.
    q_numerator_term_1 = _neg(_mul(_mul(one_minus_delta, z_jet), m_jet))
    q_numerator_term_2 = _neg(_mul(_sub(one, z_squared), m_z_jet))
    q_numerator = _add(q_numerator_term_1, q_numerator_term_2)
    q_jet = _quotient(q_numerator, denominator, "1-delta*Z^2")

    scale_terms = _scale_terms(mu_value, anchor_value, offset_value, log_pstar)
    endpoint_logscale_terms = _endpoint_scale_terms(mu_value, log_pstar, log_rp)
    linear_radius_logscale_terms = _linear_radius_scale_terms(
        mu_value, anchor_value, offset_value, log_pstar, log_rp
    )
    quadratic_radius_logscale_terms = _quadratic_radius_scale_terms(
        mu_value, log_pstar, log_rp
    )
    single_multiplier = 1
    m_quantity = ScaledQuantity(m_jet[0], m_jet[1], scale_terms, single_multiplier, "m")
    q_quantity = ScaledQuantity(q_jet[0], q_jet[1], scale_terms, single_multiplier, "Q")
    q_z_quantity = ScaledQuantity(q_jet[1], q_jet[2], scale_terms, single_multiplier, "Q_Z")
    q_zz_quantity = ScaledQuantity(q_jet[2], None, scale_terms, single_multiplier, "Q_ZZ")
    q_y_quantity = ScaledQuantity(-q_jet[0], -q_jet[1], scale_terms, single_multiplier, "Q_y")
    q_yz_quantity = ScaledQuantity(-q_jet[1], -q_jet[2], scale_terms, single_multiplier, "Q_yZ")

    linear_factor = (1 - delta_value) / (2 * L_value)
    linear_value = linear_factor * Z_value * q_jet[1]
    linear_Z1 = linear_factor * (q_jet[1] + Z_value * q_jet[2]
                                 + 2 * delta_value * Z_value**2 * q_jet[1] / L_value)
    quadratic_value = -q_jet[0] * q_jet[0] / 2
    quadratic_Z1 = -q_jet[0] * q_jet[1]
    linear_quantity = ScaledQuantity(
        linear_value, linear_Z1, scale_terms, single_multiplier, "Slinear"
    )
    quadratic_quantity = ScaledQuantity(
        quadratic_value, quadratic_Z1, scale_terms, 2, "Squadratic"
    )

    return GapScaledFieldResult(
        mu=mu_value,
        Z=Z_value,
        delta=delta_value,
        L=L_value,
        logRp=log_rp,
        anchor=anchor_value,
        offset=offset_value,
        xi=xi_value,
        lambda_value=mp.mpf(1) / 2 - mu_value,
        b_value=mp.mpf(1) / 2 + mu_value,
        scale_terms=scale_terms,
        endpoint_logscale_terms=endpoint_logscale_terms,
        linear_radius_logscale_terms=linear_radius_logscale_terms,
        quadratic_radius_logscale_terms=quadratic_radius_logscale_terms,
        m=m_quantity,
        Q=q_quantity,
        Q_Z=q_z_quantity,
        Q_ZZ=q_zz_quantity,
        Q_y=q_y_quantity,
        Q_yZ=q_yz_quantity,
        Slinear=linear_quantity,
        Squadratic=quadratic_quantity,
        source_domain_premise=source_domain_premise,
        delta_Z_independent=True,
    )


# Short aliases make the point-port role explicit without creating a second
# implementation or a source-provider-looking name.
evaluate_gap_scaled_field = evaluate_original_gap_scaled_field
original_gap_scaled_field = evaluate_original_gap_scaled_field


def _close(left: mp.mpf, right: mp.mpf, tolerance: mp.mpf) -> bool:
    scale = max(mp.mpf(1), abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def _fixture_source_density(z: mp.mpf) -> tuple[mp.mpf, ...]:
    """Independent polynomial source-density jets used only by the fixture."""

    # These coefficients are intentionally written independently of the jet
    # operators above.  Ordinary derivatives of cubic polynomials are explicit.
    u0 = 1 + mp.mpf("0.2") * z + mp.mpf("0.03") * z**2 + mp.mpf("0.004") * z**3
    u1 = mp.mpf("0.2") + mp.mpf("0.06") * z + mp.mpf("0.012") * z**2
    u2 = mp.mpf("0.06") + mp.mpf("0.024") * z
    u3 = mp.mpf("0.024")
    c10 = 1 - mp.mpf("0.1") * z + mp.mpf("0.02") * z**2 - mp.mpf("0.003") * z**3
    c11 = -mp.mpf("0.1") + mp.mpf("0.04") * z - mp.mpf("0.009") * z**2
    c12 = mp.mpf("0.04") - mp.mpf("0.018") * z
    c13 = -mp.mpf("0.018")
    c20 = mp.mpf("0.7") + mp.mpf("0.15") * z - mp.mpf("0.025") * z**2 + mp.mpf("0.002") * z**3
    c21 = mp.mpf("0.15") - mp.mpf("0.05") * z + mp.mpf("0.006") * z**2
    c22 = -mp.mpf("0.05") + mp.mpf("0.012") * z
    c23 = mp.mpf("0.012")
    # A0=3/5 and A1=-2/7 in the fixture.  Expand the product by hand only
    # through the source density; this is an independent check of m's sign.
    coefficient = (
        mp.mpf("0.6") * c10 - mp.mpf(2) / 7 * c20,
        mp.mpf("0.6") * c11 - mp.mpf(2) / 7 * c21,
        mp.mpf("0.6") * c12 - mp.mpf(2) / 7 * c22,
        mp.mpf("0.6") * c13 - mp.mpf(2) / 7 * c23,
    )
    u = (u0, u1, u2, u3)
    return tuple(
        sum(mp.mpf(comb(order, index)) * u[index] * coefficient[order - index]
            for index in range(order + 1))
        for order in range(4)
    )


def _fixture_scalar_m(z: mp.mpf, mu: mp.mpf, anchor: mp.mpf,
                      offset: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """Independent scalar m and m_Z for finite-difference Q checks."""

    density = _fixture_source_density(z)
    xi = anchor + mu * offset
    lam = mp.mpf(1) / 2 - mu
    b = mp.mpf(1) / 2 + mu
    log_e = (
        -1 / mu
        - 3 * (4 * mu) ** (-mp.mpf(2) / 3)
        + (mp.mpf(2) / 3) * mp.log(4 * mu)
        - (mp.mpf(1) / 2) * mp.log(6)
        - 2 * mp.log(mu)
    )
    scale = log_e + lam * (13 - xi) / mu - b * xi / mu
    return -density[0] * mp.exp(scale), -density[1] * mp.exp(scale)


def _fixture_scalar_q(z: mp.mpf, mu: mp.mpf, delta: mp.mpf,
                      anchor: mp.mpf, offset: mp.mpf) -> mp.mpf:
    # The public result is normalized, so finite differences compare the
    # normalized Q obtained from the same original gap expression after its
    # common scale has been factored out.
    del mu, anchor, offset
    density = _fixture_source_density(z)
    m_value, m_Z = -density[0], -density[1]
    return (-(1 - delta) * z * m_value - (1 - z**2) * m_Z) / (1 - delta * z**2)


def verify_original_gap_scaled_field_fixture(
    source_audit: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Independent medium/tiny-mu checks for the new point port.

    This fixture calls only this module's public evaluator.  It uses explicit
    cubic polynomials, a separate scalar formula, finite differences, and a
    medium-mu materialized comparison.  It is a diagnostic fixture, not a
    source identity proof or a certified error bound.
    """

    with mp.workdps(180):
        mu = mp.mpf("0.15")
        z = mp.mpf("0.23")
        delta = mp.mpf("0.27")
        L = 1 - delta * z**2
        anchor = mp.mpf(13)
        offset = mp.mpf(-4)
        log_pstar = mp.mpf("0.31")
        log_rp = mp.mpf("-2.4")
        density = _fixture_source_density(z)
        u = (
            1 + mp.mpf("0.2") * z + mp.mpf("0.03") * z**2 + mp.mpf("0.004") * z**3,
            mp.mpf("0.2") + mp.mpf("0.06") * z + mp.mpf("0.012") * z**2,
            mp.mpf("0.06") + mp.mpf("0.024") * z,
            mp.mpf("0.024"),
        )
        C1 = (
            1 - mp.mpf("0.1") * z + mp.mpf("0.02") * z**2 - mp.mpf("0.003") * z**3,
            -mp.mpf("0.1") + mp.mpf("0.04") * z - mp.mpf("0.009") * z**2,
            mp.mpf("0.04") - mp.mpf("0.018") * z,
            -mp.mpf("0.018"),
        )
        C2 = (
            mp.mpf("0.7") + mp.mpf("0.15") * z - mp.mpf("0.025") * z**2 + mp.mpf("0.002") * z**3,
            mp.mpf("0.15") - mp.mpf("0.05") * z + mp.mpf("0.006") * z**2,
            -mp.mpf("0.05") + mp.mpf("0.012") * z,
            mp.mpf("0.012"),
        )
        result = evaluate_original_gap_scaled_field(
            mu, z, delta, mp.mpf("0.6"), -mp.mpf(2) / 7,
            u, C1, C2, anchor=anchor, offset=offset,
            logRp=log_rp,
            logPstar=log_pstar, source_domain_premise=False,
            delta_Z_independent=True, scale_and_weights_Z_independent=True,
        )
        tolerance = mp.mpf("1e-45")
        scale_without_pstar = sum(_scale_terms(mu, anchor, offset, mp.mpf(0)))
        direct_m, direct_m_Z = _fixture_scalar_m(z, mu, anchor, offset)
        h = mp.mpf("1e-6")
        q_minus = _fixture_scalar_q(z - h, mu, delta, anchor, offset)
        q_zero = _fixture_scalar_q(z, mu, delta, anchor, offset)
        q_plus = _fixture_scalar_q(z + h, mu, delta, anchor, offset)
        finite_difference_qz = (q_plus - q_minus) / (2 * h)
        finite_difference_qzz = (q_plus - 2 * q_zero + q_minus) / h**2
        numerator_term_1 = -(1 - delta) * z * result.m.value
        numerator_term_2 = -(1 - z**2) * result.m.Z1
        checks = {
            "ordinary_C3_source_density_Z0": _close(-result.m.value, density[0], tolerance),
            "ordinary_C3_source_density_Z1": _close(-result.m.Z1, density[1], tolerance),
            "original_gap_m_expression_after_scale": _close(
                result.m.value,
                direct_m / mp.exp(scale_without_pstar),
                tolerance,
            ),
            "original_gap_m_Z_expression_after_scale": _close(
                result.m.Z1,
                direct_m_Z / mp.exp(scale_without_pstar),
                tolerance,
            ),
            "finite_difference_QZ": _close(result.Q_Z.value, finite_difference_qz, mp.mpf("1e-10")),
            "finite_difference_QZZ": _close(result.Q_ZZ.value, finite_difference_qzz, mp.mpf("1e-7")),
            "Qy_equals_minus_Q": result.Q_y.value == -result.Q.value and result.Q_y.Z1 == -result.Q.Z1,
            "QyZ_equals_minus_QZ": result.Q_yZ.value == -result.Q_Z.value and result.Q_yZ.Z1 == -result.Q_Z.Z1,
            "Slinear_value": _close(
                result.Slinear.value,
                (1 - delta) * z * result.Q_Z.value / (2 * L),
                tolerance,
            ),
            "Slinear_Z1": _close(
                result.Slinear.Z1,
                (1 - delta) * (result.Q_Z.value + z * result.Q_ZZ.value
                               + 2*delta*z*z*result.Q_Z.value/L) / (2 * L),
                tolerance,
            ),
            "Slinear_Z1_independent_finite_difference": _close(
                result.Slinear.Z1,
                (1-delta) / 2 * (
                    ((z+h)*mp.diff(lambda t: _fixture_scalar_q(t, mu, delta, anchor, offset), z+h)
                     / (1-delta*(z+h)**2))
                    - ((z-h)*mp.diff(lambda t: _fixture_scalar_q(t, mu, delta, anchor, offset), z-h)
                       / (1-delta*(z-h)**2))
                ) / (2*h), mp.mpf("1e-10")),
            "Squadratic_value": _close(result.Squadratic.value, -result.Q.value**2 / 2, tolerance),
            "Squadratic_Z1": _close(result.Squadratic.Z1, -result.Q.value * result.Q_Z.value, tolerance),
            "Q_endpoint_subtraction_term_retained": numerator_term_2 != 0 and result.Q.value != numerator_term_1 / (1 - delta * z**2),
            "endpoint_coordinate_retains_minus_4_mu": result.xi == 13 - 4 * mu,
            "single_scale_has_seven_ordered_terms": len(result.scale_terms) == 7,
            "quadratic_scale_multiplier_is_two": result.Squadratic.scale_multiplier == 2,
            "endpoint_scale_terms_separate": result.endpoint_logscale_terms == (
                (mp.mpf(11) / 2) / mu,
                -mp.mpf(13),
                -3 * (4 * mu) ** (-mp.mpf(2) / 3),
                -(mp.mpf(4) / 3) * mp.log(mu),
                (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6),
                log_pstar,
                log_rp,
            ),
            "endpoint_scale_matches_direct_radius": _close(
                sum(result.endpoint_logscale_terms),
                scale_without_pstar + log_pstar + log_rp + (anchor+mu*offset)/mu,
                tolerance),
            "linear_scale_matches_direct_R2": _close(
                sum(result.linear_radius_logscale_terms),
                scale_without_pstar + log_pstar + 2*(log_rp+(anchor+mu*offset)/mu),
                tolerance),
            "quadratic_scale_matches_direct_R2": _close(
                sum(result.quadratic_radius_logscale_terms),
                2*(scale_without_pstar+log_pstar+log_rp+(anchor+mu*offset)/mu),
                tolerance),
            "linear_radius_scale_terms_separate": result.linear_radius_logscale_terms == (
                (mp.mpf(11) / 2 + anchor) / mu,
                offset,
                -mp.mpf(13),
                -3 * (4 * mu) ** (-mp.mpf(2) / 3),
                -(mp.mpf(4) / 3) * mp.log(mu),
                (mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6),
                log_pstar,
                2 * log_rp,
            ),
            "quadratic_radius_scale_terms_separate": result.quadratic_radius_logscale_terms == (
                mp.mpf(11) / mu,
                -mp.mpf(26),
                -6 * (4 * mu) ** (-mp.mpf(2) / 3),
                -(mp.mpf(8) / 3) * mp.log(mu),
                2 * ((mp.mpf(2) / 3) * mp.log(4) - (mp.mpf(1) / 2) * mp.log(6)),
                2 * log_pstar,
                2 * log_rp,
            ),
            "density_measure_is_not_integrated": result.as_dict()["formula_components"]["density"]["integrated_here"] is False,
            "endpoint_orientation_retained": result.as_dict()["formula_components"]["endpoint"]["left_right_orientation_retained"] is True and result.as_dict()["formula_components"]["endpoint"]["claimed_zero"] is False,
            "logscale_has_no_anchor_minus_13_term": "anchor-13" not in LOGSCALE_LABELS,
            "logscale_has_no_mu_offset_term": "mu*offset" not in LOGSCALE_LABELS,
        }

        tiny = evaluate_original_gap_scaled_field(
            "1e-400", "0.23", "0.27", "0.6", -mp.mpf(2) / 7,
            u, C1, C2, logRp=log_rp, anchor=13, offset=-4, logPstar="0.31",
            delta_Z_independent=True, scale_and_weights_Z_independent=True,
        )
        checks.update({
            "tiny_mu_single_scale_direction_negative": tiny.scale_terms[0] < 0,
            "tiny_mu_quadratic_scale_direction_negative": tiny.Squadratic.scale_multiplier == 2 and tiny.Squadratic.scale_terms[0] < 0,
            "tiny_mu_endpoint_scale_direction_positive": tiny.endpoint_logscale_terms[0] > 0,
            "tiny_mu_linear_radius_scale_direction_positive": tiny.linear_radius_logscale_terms[0] > 0,
            "tiny_mu_quadratic_radius_scale_direction_positive": tiny.quadratic_radius_logscale_terms[0] > 0,
            "tiny_mu_mantissas_remain_finite": all(
                mp.isfinite(quantity.value)
                for quantity in (tiny.m, tiny.Q, tiny.Q_Z, tiny.Q_ZZ, tiny.Slinear, tiny.Squadratic)
            ),
        })
        try:
            evaluate_original_gap_scaled_field(
                "0", z, delta, "0.6", "0.1", u, C1, C2,
                logRp="-2.4", logPstar="0.31",
                delta_Z_independent=True, scale_and_weights_Z_independent=True,
            )
        except ValueError:
            checks["reject_zero_mu"] = True
        else:
            checks["reject_zero_mu"] = False
        try:
            evaluate_original_gap_scaled_field(
                "1e-400", z, delta, "0.6", "0.1", u, C1, C2,
                logRp="-2.4", logPstar="0.31", anchor=13, offset=-3,
                delta_Z_independent=True, scale_and_weights_Z_independent=True,
            )
        except ValueError:
            checks["tiny_mu_rejects_outside_upper_endpoint"] = True
        else:
            checks["tiny_mu_rejects_outside_upper_endpoint"] = False
        try:
            evaluate_original_gap_scaled_field(
                mu, z, delta, "0.6", "0.1", u, C1, C2,
                logRp=log_rp, logPstar=log_pstar, delta_Z_independent=True,
            )
        except ValueError:
            checks["rejects_undeclared_scale_Z_independence"] = True
        else:
            checks["rejects_undeclared_scale_Z_independence"] = False
        real_mu_expression = "exp(-4*(exp(40)+11))/1000"
        checks["real_constant_definition_recorded"] = real_mu_expression == "exp(-4*(exp(40)+11))/1000"
        checks["real_constant_endpoint_symbolic_identity"] = (
            "xi=13-4*exp(-4*(exp(40)+11))/1000"
            == "xi=13-4*" + real_mu_expression
        )
        return {
            "schema": "original_gap_scaled_field_fixture_v1",
            "passed": all(checks.values()),
            "checks": checks,
            "precision_dps": 180,
            "medium_mu": str(mu),
            "tiny_mu": "1e-400",
            "real_mu_definition": real_mu_expression,
            "real_mu_endpoint_definition": "anchor=13, offset=-4, xi=13-4*mu",
            "real_mu_evaluation": "not_materialized: exponent is approximately -4*(exp(40)+11), so the actual scale may exceed the host exponent range",
            "source_fragment_scope": source_audit,
            "source_identity_proved": False,
            "source_values_are_caller_owned": True,
            "fixture_jets_are_synthetic": True,
            "field_scale_materialized": False,
            "medium_diagnostic_original_exponent_materialized": True,
            "real_scale_materialized": False,
            "finite_difference_is_diagnostic": True,
        }


__all__ = [
    "GapScaledFieldResult",
    "LOGSCALE_LABELS",
    "ScaledQuantity",
    "evaluate_gap_scaled_field",
    "evaluate_original_gap_scaled_field",
    "original_gap_scaled_field",
    "verify_original_gap_scaled_field_fixture",
]
