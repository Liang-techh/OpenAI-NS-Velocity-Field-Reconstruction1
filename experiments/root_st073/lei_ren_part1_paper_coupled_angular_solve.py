"""Signed, scale preserving solve of the coupled angular bump rows.

This module is a small algebra adapter for the two compact bumps in (7.21).
It consumes the weights already owned by
``ContinuousAngularCorrectionProvider`` and solves for the coefficients
directly.  In particular, it does not divide by ``r``: an angular target may
be negative, zero, or much smaller than the working scale.

For fixed ``mu`` and bump weights, the rows solved here are

``A*(exp(-3*(1-mu))*d1 + exp(-(1-mu))*d2) = r``

``B*(exp(3*(1+2*mu))*d1 + exp(1+2*mu)*d2)``
`` + D/2*(exp(3*(1+2*mu))*d1**2 + exp(1+2*mu)*d2**2) = s``.

The scalar quadratic is evaluated with the cancellation avoiding root and its
small branch is the branch continuous with the linear (``D=0``) solve.  The
``rZ`` and ``sZ`` arguments are supplied analytic target derivatives; the
coefficient tangent is obtained by solving the implicit two by two Jacobian.

The returned receipt reports direct and normalized residuals, branch and
smallness diagnostics, and a sufficient bump multiplier positivity check when
``beta_inf`` is available.  It makes no claim about a cone, a global field,
or pressure/NS closure.  The caller assembles actual targets and decides when
or whether to install a provider separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Mapping
import json
import math
from pathlib import Path
import sys

import mpmath as mp


_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_EQUATION = "(7.21)"
ALGEBRA_VERSION = "v2"
DEFAULT_PRECISION = 120


def _mp(value: Any, name: str = "value") -> mp.mpf:
    """Parse a scalar without routing it through binary64 first."""

    if isinstance(value, mp.mpf):
        result = value
    else:
        try:
            # ``str`` preserves Decimal and decimal-string exponents while
            # avoiding an accidental float conversion for numpy scalars.
            result = mp.mpf(str(value))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be a finite real scalar") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _read(source: Any, name: str, default: Any = None) -> Any:
    """Read either an attribute or a mapping entry from a provider/weights."""

    if hasattr(source, name):
        return getattr(source, name)
    if isinstance(source, Mapping) and name in source:
        return source[name]
    return default


def _signed(value: mp.mpf | None, precision: int) -> dict[str, Any] | None:
    """Serialize an MP value while retaining its sign and exponent."""

    if value is None:
        return None
    with mp.workdps(max(int(precision), 32)):
        x = _mp(value)
        if x == 0:
            return {
                "sign": 0,
                "log_abs": None,
                "arbitrary_exponent_value": "0",
            }
        return {
            "sign": 1 if x > 0 else -1,
            "log_abs": mp.nstr(mp.log(abs(x)), int(precision)),
            "arbitrary_exponent_value": mp.nstr(x, int(precision)),
        }


def _nstr(value: Any, precision: int) -> str | None:
    if value is None:
        return None
    return mp.nstr(_mp(value), int(precision))


def _relative_residual(residual: mp.mpf, left: mp.mpf, right: mp.mpf) -> mp.mpf:
    """Scale a row by the larger represented side, including zero rows."""

    scale = max(abs(left), abs(right))
    return mp.mpf(0) if scale == 0 else residual / scale


def _target_residual(residual: mp.mpf, target: mp.mpf) -> mp.mpf:
    """Return a target-relative residual without dividing by a zero target."""

    if target == 0:
        return mp.mpf(0) if residual == 0 else mp.sign(residual) * mp.inf
    return residual / abs(target)


def _log10_ratio(numerator: mp.mpf, denominator: mp.mpf) -> mp.mpf | None:
    if numerator == 0 or denominator == 0:
        return None
    return mp.log10(abs(numerator) / abs(denominator))


def _small_root(
    c2: mp.mpf, c1: mp.mpf, c0: mp.mpf, discriminant: mp.mpf
) -> tuple[mp.mpf, mp.mpf | None, bool, str]:
    """Return the q-continuation root, its alternate, and branch status."""

    if c2 == 0:
        if c1 == 0:
            if c0 == 0:
                return mp.mpf(0), None, True, "degenerate_zero_row"
            raise ValueError("angular scalar equation has no finite linear root")
        return -c0 / c1, None, True, "linear_D_zero"

    if discriminant < 0:
        raise ValueError("angular scalar equation has no real root")
    root_discriminant = mp.sqrt(discriminant)
    if c1 == 0:
        # There is no q=0 continuation at this singular linear limit.  The
        # lower magnitude real root is a deterministic diagnostic choice.
        plus = root_discriminant / (2 * c2)
        minus = -plus
        selected, alternate = (plus, minus) if abs(plus) <= abs(minus) else (minus, plus)
        return selected, alternate, False, "quadratic_c1_zero_minimum_magnitude"

    # sign(c1) keeps the denominator away from cancellation.  This is the
    # root that tends to -c0/c1 as q tends to zero for either sign of c1.
    signed_sqrt = mp.sign(c1) * root_discriminant
    denominator = c1 + signed_sqrt
    if denominator == 0:
        # This can only occur for a repeated/degenerate root.  The direct
        # expression remains exact and is preferable to inventing a scale.
        selected = (-c1 + signed_sqrt) / (2 * c2)
    else:
        selected = -2 * c0 / denominator
    alternate = -c1 / c2 - selected
    return selected, alternate, True, "small_q_continuation"


def _provider_parameters(provider: Any, precision: int | None) -> tuple[int, mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf | None]:
    """Read exactly the existing continuous provider weights interface."""

    weights = _read(provider, "weights")
    if weights is None:
        # A weights object is accepted for small independent fixtures.  It
        # still has the same A_mu/B_mu/D_mu/beta_inf fields.
        weights = provider
    mu_raw = _read(provider, "mu", _read(weights, "mu"))
    if mu_raw is None:
        raise TypeError("provider must expose mu and weights.A_mu/B_mu/D_mu")
    provider_precision = _read(provider, "precision", _read(provider, "decimal_precision", DEFAULT_PRECISION))
    dps = max(int(precision if precision is not None else provider_precision), 64)
    with mp.workdps(dps):
        mu = _mp(mu_raw, "mu")
        A = _mp(_read(weights, "A_mu"), "A_mu")
        B = _mp(_read(weights, "B_mu"), "B_mu")
        D = _mp(_read(weights, "D_mu"), "D_mu")
        beta_inf_raw = _read(weights, "beta_inf")
        beta_inf = None if beta_inf_raw is None else _mp(beta_inf_raw, "beta_inf")
    if A <= 0 or B <= 0 or D < 0:
        raise ValueError("A_mu and B_mu must be positive and D_mu must be nonnegative")
    if beta_inf is not None and beta_inf < 0:
        raise ValueError("beta_inf must be nonnegative")
    return dps, mu, A, B, D, beta_inf


@dataclass(frozen=True)
class CoupledAngularSolution:
    """High precision coefficient and tangent receipt.

    MP values are retained in the fields used by numerical consumers.  The
    ``as_dict`` method adds decimal and signed-log encodings for JSON receipts.
    ``solution["d1"]`` is supported as a convenience for provider-like code.
    """

    precision: int
    mu: mp.mpf
    A_mu: mp.mpf
    B_mu: mp.mpf
    D_mu: mp.mpf
    beta_inf: mp.mpf | None
    r: mp.mpf
    s: mp.mpf
    rZ: mp.mpf
    sZ: mp.mpf
    d1: mp.mpf | None
    d2: mp.mpf | None
    d1_Z: mp.mpf | None
    d2_Z: mp.mpf | None
    alpha: mp.mpf
    beta: mp.mpf
    rho: mp.mpf
    sigma: mp.mpf
    q: mp.mpf
    quadratic_coefficients: tuple[mp.mpf, mp.mpf, mp.mpf]
    discriminant: mp.mpf
    alternate_d1: mp.mpf | None
    alternate_d2: mp.mpf | None
    linear_limit_d1: mp.mpf | None
    jacobian: tuple[tuple[mp.mpf, mp.mpf], tuple[mp.mpf, mp.mpf]] | None
    jacobian_determinant: mp.mpf | None
    row_residuals: tuple[mp.mpf, mp.mpf] | None
    normalized_row_residuals: tuple[mp.mpf, mp.mpf] | None
    target_normalized_row_residuals: tuple[mp.mpf, mp.mpf] | None
    normalized_equation_residuals: tuple[mp.mpf, mp.mpf] | None
    derivative_residuals: tuple[mp.mpf, mp.mpf] | None
    normalized_derivative_residuals: tuple[mp.mpf, mp.mpf] | None
    branch: str
    branch_checks: dict[str, Any]
    smallness_checks: dict[str, Any]
    positivity_checked: bool
    positive: bool | None
    positivity_margin: mp.mpf | None
    target_scale_checks: dict[str, Any]
    accepted: bool
    rejection_reason: str | None

    @property
    def coefficients(self) -> tuple[mp.mpf | None, mp.mpf | None]:
        return self.d1, self.d2

    @property
    def tangent(self) -> tuple[mp.mpf | None, mp.mpf | None]:
        return self.d1_Z, self.d2_Z

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        return self.as_dict()[key]

    def as_dict(self, precision: int | None = None) -> dict[str, Any]:
        """Return a JSON-ready receipt without losing tiny signed values."""

        dps = int(precision if precision is not None else self.precision)

        def encode(value: Any) -> Any:
            if isinstance(value, mp.mpf):
                return mp.nstr(value, dps)
            if isinstance(value, tuple):
                return [encode(item) for item in value]
            if isinstance(value, list):
                return [encode(item) for item in value]
            if isinstance(value, dict):
                return {str(key): encode(item) for key, item in value.items()}
            return value

        receipt: dict[str, Any] = {
            "source": SOURCE,
            "source_equation": SOURCE_EQUATION,
            "algebra_version": ALGEBRA_VERSION,
            "precision": self.precision,
            "mu": encode(self.mu),
            "A_mu": encode(self.A_mu),
            "B_mu": encode(self.B_mu),
            "D_mu": encode(self.D_mu),
            "beta_inf": encode(self.beta_inf),
            "r": encode(self.r),
            "s": encode(self.s),
            "rZ": encode(self.rZ),
            "sZ": encode(self.sZ),
            "d1": encode(self.d1),
            "d2": encode(self.d2),
            "d1_Z": encode(self.d1_Z),
            "d2_Z": encode(self.d2_Z),
            "d1_signed": _signed(self.d1, dps),
            "d2_signed": _signed(self.d2, dps),
            "d1_Z_signed": _signed(self.d1_Z, dps),
            "d2_Z_signed": _signed(self.d2_Z, dps),
            "alpha": encode(self.alpha),
            "beta": encode(self.beta),
            "rho": encode(self.rho),
            "sigma": encode(self.sigma),
            "q": encode(self.q),
            "quadratic_coefficients": encode(self.quadratic_coefficients),
            "discriminant": encode(self.discriminant),
            "alternate_d1": encode(self.alternate_d1),
            "alternate_d2": encode(self.alternate_d2),
            "linear_limit_d1": encode(self.linear_limit_d1),
            "jacobian": encode(self.jacobian),
            "jacobian_determinant": encode(self.jacobian_determinant),
            "row_residuals": encode(self.row_residuals),
            "normalized_row_residuals": encode(self.normalized_row_residuals),
            "target_normalized_row_residuals": encode(self.target_normalized_row_residuals),
            "normalized_equation_residuals": encode(self.normalized_equation_residuals),
            "derivative_residuals": encode(self.derivative_residuals),
            "normalized_derivative_residuals": encode(self.normalized_derivative_residuals),
            "branch": self.branch,
            "branch_checks": encode(self.branch_checks),
            "smallness_checks": encode(self.smallness_checks),
            "positivity_checked": self.positivity_checked,
            "positive": self.positive,
            "positivity_margin": encode(self.positivity_margin),
            "target_scale_checks": encode(self.target_scale_checks),
            "accepted": self.accepted,
            "rejection_reason": self.rejection_reason,
            "scope": (
                "Finite-dimensional signed angular coefficient solve and its "
                "analytic target tangent only; no cone, global closure, or "
                "Navier--Stokes claim."
            ),
        }
        return receipt


def solve_coupled_angular(
    provider: Any,
    r: Any,
    s: Any | None = None,
    rZ: Any = 0,
    sZ: Any = 0,
    *,
    s_H: Any | None = None,
    precision: int | None = None,
    require_positive: bool = False,
    smallness_limit: Any | None = None,
    return_rejected: bool = True,
) -> CoupledAngularSolution:
    """Solve the signed (7.21) angular rows using provider-owned weights.

    ``r`` and ``s`` are direct row targets, not positive logarithms and not
    ratios by ``r``.  ``rZ`` and ``sZ`` are the supplied analytic target
    derivatives.  ``s_H`` is accepted as a spelling compatible with the
    legacy heat-row caller; it is only used when ``s`` is omitted.
    If the supplied rows are unreachable on a real branch, the default is a
    diagnostic receipt with ``accepted=False``.  Set ``return_rejected=False``
    for strict exception behavior.

    For callers that keep a target receipt together, ``r`` may also be a
    mapping/object exposing ``r``, ``s``, ``rZ`` and ``sZ``; this is merely an
    input convenience and does not mutate that receipt.
    """

    if s is None:
        bundled_s = _read(r, "s")
        if bundled_s is not None and s_H is None:
            bundled = r
            r = _read(bundled, "r")
            s = bundled_s
            rZ = _read(bundled, "rZ", _read(bundled, "r_Z", rZ))
            sZ = _read(bundled, "sZ", _read(bundled, "s_Z", sZ))
        else:
            if s_H is None:
                raise TypeError("s (or s_H) is required")
            s = s_H
    elif s_H is not None:
        raise TypeError("pass only one of s and s_H")

    dps, mu, A, B, D, beta_inf = _provider_parameters(provider, precision)
    # A moderate guard keeps the direct subtraction rho-alpha*d1 from
    # erasing a prescribed small second coefficient.  MP exponents do not
    # underflow, but cancellation still consumes significant digits.
    work_dps = dps + 40
    with mp.workdps(work_dps):
        r_mp = _mp(r, "r")
        s_mp = _mp(s, "s")
        rZ_mp = _mp(rZ, "rZ")
        sZ_mp = _mp(sZ, "sZ")
        lam = 1 - mu
        kap = 1 + 2 * mu
        alpha = mp.exp(-2 * lam)
        beta = mp.exp(-2 * kap)
        rho = mp.exp(lam) * r_mp / A
        sigma = mp.exp(-3 * kap) * s_mp / B
        q = D / (2 * B)

        c2 = q * (1 + beta * alpha * alpha)
        c1 = 1 - beta * alpha - 2 * q * beta * alpha * rho
        c0 = beta * rho + q * beta * rho * rho - sigma
        discriminant = c1 * c1 - 4 * c2 * c0
        if discriminant < 0:
            if not return_rejected:
                raise ValueError(
                    "angular data are unreachable on a real small-branch solve: "
                    f"discriminant={mp.nstr(discriminant, 18)}"
                )
            scale_log10 = _log10_ratio(sigma, rho)
            target_scale_checks = {
                "target_scale_preserved": True,
                "r_is_zero": r_mp == 0,
                "s_is_zero": s_mp == 0,
                "normalized_rho": rho,
                "normalized_sigma": sigma,
                "log10_abs_sigma_over_abs_rho": scale_log10,
                "pressure_target_below_working_precision_relative_to_r": bool(
                    scale_log10 is not None and scale_log10 < -(dps - 8)
                ),
                "arbitrary_exponent_representation": True,
                "warning": (
                    "s is more than the declared working precision below r in the "
                    "normalized rows; rejected branch diagnostics retain the "
                    "supplied target scale"
                    if scale_log10 is not None and scale_log10 < -(dps - 8)
                    else None
                ),
            }
            branch_checks = {
                "discriminant_nonnegative": False,
                "quadratic_coefficient_nonnegative": c2 >= 0,
                "linear_coefficient_nonzero": c1 != 0,
                "source_small_branch_linear_sign": c1 > 0,
                "small_branch_continuation": False,
                "selected_root_is_small_branch": False,
                "tangent_jacobian_nonsingular": False,
                "real_solution": False,
            }
            smallness_checks = {
                "selected_root_is_q_continuation": False,
                "bound_applied": smallness_limit is not None,
                "within_smallness_limit": False,
                "rejected_before_coefficient_selection": True,
            }
            if smallness_limit is not None:
                smallness_checks["limit"] = _mp(smallness_limit, "smallness_limit")
            return CoupledAngularSolution(
                precision=dps,
                mu=mu,
                A_mu=A,
                B_mu=B,
                D_mu=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                d1=None,
                d2=None,
                d1_Z=None,
                d2_Z=None,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                quadratic_coefficients=(c2, c1, c0),
                discriminant=discriminant,
                alternate_d1=None,
                alternate_d2=None,
                linear_limit_d1=None if c1 == 0 else -c0 / c1,
                jacobian=None,
                jacobian_determinant=None,
                row_residuals=None,
                normalized_row_residuals=None,
                target_normalized_row_residuals=None,
                normalized_equation_residuals=None,
                derivative_residuals=None,
                normalized_derivative_residuals=None,
                branch="rejected_unreachable_discriminant",
                branch_checks=branch_checks,
                smallness_checks=smallness_checks,
                positivity_checked=beta_inf is not None,
                positive=False if beta_inf is not None else None,
                positivity_margin=None,
                target_scale_checks=target_scale_checks,
                accepted=False,
                rejection_reason="negative_discriminant",
            )
        d1, alternate_d1, continuation, branch = _small_root(c2, c1, c0, discriminant)
        d2 = rho - alpha * d1
        alternate_d2 = None if alternate_d1 is None else rho - alpha * alternate_d1
        linear_limit_d1 = None if c1 == 0 else -c0 / c1

        # The implicit tangent differentiates the normalized equations.  The
        # weights are held fixed here; their own Z dependence belongs to the
        # target assembly/provider layer.
        j11, j12 = alpha, mp.mpf(1)
        j21 = 1 + 2 * q * d1
        j22 = beta * (1 + 2 * q * d2)
        rhs1 = mp.exp(lam) * rZ_mp / A
        rhs2 = mp.exp(-3 * kap) * sZ_mp / B
        jacobian_determinant = j11 * j22 - j12 * j21
        if jacobian_determinant == 0:
            raise ValueError("angular small-branch tangent Jacobian is singular")
        d1_Z = (rhs1 * j22 - j12 * rhs2) / jacobian_determinant
        d2_Z = (j11 * rhs2 - rhs1 * j21) / jacobian_determinant

        # Direct physical rows.
        angular_d1 = A * mp.exp(-3 * lam)
        angular_d2 = A * mp.exp(-lam)
        pressure_d1 = B * mp.exp(3 * kap)
        pressure_d2 = B * mp.exp(kap)
        pressure_q1 = (D / 2) * mp.exp(3 * kap)
        pressure_q2 = (D / 2) * mp.exp(kap)
        angular_lhs = angular_d1 * d1 + angular_d2 * d2
        pressure_lhs = (
            pressure_d1 * d1
            + pressure_d2 * d2
            + pressure_q1 * d1 * d1
            + pressure_q2 * d2 * d2
        )
        row_residuals = (angular_lhs - r_mp, pressure_lhs - s_mp)
        normalized_equation_residuals = (
            alpha * d1 + d2 - rho,
            d1 + beta * d2 + q * (d1 * d1 + beta * d2 * d2) - sigma,
        )
        normalized_row_residuals = (
            _relative_residual(normalized_equation_residuals[0], alpha * d1 + d2, rho),
            _relative_residual(
                normalized_equation_residuals[1],
                d1 + beta * d2 + q * (d1 * d1 + beta * d2 * d2),
                sigma,
            ),
        )
        target_normalized_row_residuals = (
            _target_residual(row_residuals[0], r_mp),
            _target_residual(row_residuals[1], s_mp),
        )

        angular_lhs_Z = angular_d1 * d1_Z + angular_d2 * d2_Z
        pressure_lhs_Z = (
            (pressure_d1 + 2 * pressure_q1 * d1) * d1_Z
            + (pressure_d2 + 2 * pressure_q2 * d2) * d2_Z
        )
        derivative_residuals = (angular_lhs_Z - rZ_mp, pressure_lhs_Z - sZ_mp)
        normalized_derivative_residuals = (
            _relative_residual(derivative_residuals[0], angular_lhs_Z, rZ_mp),
            _relative_residual(derivative_residuals[1], pressure_lhs_Z, sZ_mp),
        )

        if beta_inf is None:
            positivity_checked = False
            positive = None
            positivity_margin = None
        else:
            positivity_checked = True
            positivity_margin = min(1 + beta_inf * d1, 1 + beta_inf * d2)
            positive = positivity_margin >= 0

        selected_distance = None
        alternate_distance = None
        if linear_limit_d1 is not None:
            selected_distance = abs(d1 - linear_limit_d1)
            if alternate_d1 is not None:
                alternate_distance = abs(alternate_d1 - linear_limit_d1)
        q_quad = q * (d1 * d1 + beta * d2 * d2)
        linear_pressure = d1 + beta * d2
        root_magnitude_ratio = (
            None
            if alternate_d1 is None or alternate_d1 == 0
            else abs(d1 / alternate_d1)
        )
        smallness_checks: dict[str, Any] = {
            "selected_root_is_q_continuation": continuation,
            "linear_limit_d1": linear_limit_d1,
            "distance_to_linear_limit": selected_distance,
            "alternate_distance_to_linear_limit": alternate_distance,
            "selected_closer_to_linear_limit": (
                None
                if alternate_distance is None or selected_distance is None
                else selected_distance <= alternate_distance
            ),
            "selected_to_alternate_root_magnitude_ratio": root_magnitude_ratio,
            "selected_magnitude_no_larger_than_alternate": (
                None
                if root_magnitude_ratio is None
                else root_magnitude_ratio <= 1
            ),
            "max_abs_coefficient": max(abs(d1), abs(d2)),
            "quadratic_pressure_term": q_quad,
            "linear_pressure_term": linear_pressure,
            "quadratic_to_total_pressure_ratio": _relative_residual(
                q_quad, linear_pressure + q_quad, mp.mpf(0)
            ),
            "bound_applied": smallness_limit is not None,
        }
        if smallness_limit is None:
            smallness_checks["within_smallness_limit"] = None
        else:
            limit = _mp(smallness_limit, "smallness_limit")
            if limit < 0:
                raise ValueError("smallness_limit must be nonnegative")
            smallness_checks["limit"] = limit
            smallness_checks["within_smallness_limit"] = max(abs(d1), abs(d2)) <= limit

        scale_log10 = _log10_ratio(sigma, rho)
        target_scale_checks: dict[str, Any] = {
            "target_scale_preserved": True,
            "r_is_zero": r_mp == 0,
            "s_is_zero": s_mp == 0,
            "normalized_rho": rho,
            "normalized_sigma": sigma,
            "log10_abs_sigma_over_abs_rho": scale_log10,
            "pressure_target_below_working_precision_relative_to_r": (
                bool(scale_log10 is not None and scale_log10 < -(dps - 8))
            ),
            "arbitrary_exponent_representation": True,
            "warning": (
                "s is more than the declared working precision below r in the "
                "normalized rows; residuals remain diagnostic at supplied input "
                "precision"
                if scale_log10 is not None and scale_log10 < -(dps - 8)
                else None
            ),
        }
        branch_checks: dict[str, Any] = {
            "discriminant_nonnegative": discriminant >= 0,
            "quadratic_coefficient_nonnegative": c2 >= 0,
            "linear_coefficient_nonzero": c1 != 0,
            "source_small_branch_linear_sign": c1 > 0,
            "small_branch_continuation": continuation,
            "selected_root_is_small_branch": (
                continuation
                and c1 > 0
                and (
                    alternate_distance is None
                    or selected_distance is None
                    or selected_distance <= alternate_distance
                )
            ),
            "selected_magnitude_no_larger_than_alternate": smallness_checks[
                "selected_magnitude_no_larger_than_alternate"
            ],
            "tangent_jacobian_nonsingular": jacobian_determinant != 0,
            "real_solution": True,
        }

        positivity_ok = (not require_positive) or (positive is True)
        smallness_ok = smallness_checks["within_smallness_limit"]
        if smallness_ok is None:
            smallness_ok = True
        accepted = bool(
            branch_checks["selected_root_is_small_branch"]
            and branch_checks["tangent_jacobian_nonsingular"]
            and (positive is not False)
            and positivity_ok
            and smallness_ok
        )
        rejection_reason = None
        if not accepted:
            reasons: list[str] = []
            if not branch_checks["selected_root_is_small_branch"]:
                reasons.append("small_branch_check_failed")
            if not branch_checks["tangent_jacobian_nonsingular"]:
                reasons.append("singular_tangent_jacobian")
            if positive is False:
                reasons.append("positivity_check_failed")
            if smallness_ok is False:
                reasons.append("smallness_limit_failed")
            rejection_reason = ",".join(reasons) or "branch_check_failed"

        # Values are returned while the expanded context is active.  MP values
        # retain their arbitrary exponent and are not converted to float.
        return CoupledAngularSolution(
            precision=dps,
            mu=mu,
            A_mu=A,
            B_mu=B,
            D_mu=D,
            beta_inf=beta_inf,
            r=r_mp,
            s=s_mp,
            rZ=rZ_mp,
            sZ=sZ_mp,
            d1=d1,
            d2=d2,
            d1_Z=d1_Z,
            d2_Z=d2_Z,
            alpha=alpha,
            beta=beta,
            rho=rho,
            sigma=sigma,
            q=q,
            quadratic_coefficients=(c2, c1, c0),
            discriminant=discriminant,
            alternate_d1=alternate_d1,
            alternate_d2=alternate_d2,
            linear_limit_d1=linear_limit_d1,
            jacobian=((j11, j12), (j21, j22)),
            jacobian_determinant=jacobian_determinant,
            row_residuals=row_residuals,
            normalized_row_residuals=normalized_row_residuals,
            target_normalized_row_residuals=target_normalized_row_residuals,
            normalized_equation_residuals=normalized_equation_residuals,
            derivative_residuals=derivative_residuals,
            normalized_derivative_residuals=normalized_derivative_residuals,
            branch=branch,
            branch_checks=branch_checks,
            smallness_checks=smallness_checks,
            positivity_checked=positivity_checked,
            positive=positive,
            positivity_margin=positivity_margin,
            target_scale_checks=target_scale_checks,
            accepted=accepted,
            rejection_reason=rejection_reason,
        )


def solve_angular_bumps(*args: Any, **kwargs: Any) -> CoupledAngularSolution:
    """Compatibility wrapper for provider or explicit-weight callers.

    The preferred form is ``solve_angular_bumps(provider, r, s, rZ, sZ)``.
    For the older finite-dimensional API, pass ``mu``, ``A_mu``, ``B_mu``,
    ``D_mu`` and ``r`` together with either ``s`` or ``s_H`` as keywords.
    Both paths use the same signed MP implementation above.
    """

    if args or "provider" in kwargs:
        if "provider" in kwargs:
            if args:
                raise TypeError("provider must be positional or keyword, not both")
            provider = kwargs.pop("provider")
            return solve_coupled_angular(provider, **kwargs)
        return solve_coupled_angular(*args, **kwargs)
    required = ("mu", "A_mu", "B_mu", "D_mu", "r")
    missing = [name for name in required if name not in kwargs]
    if missing:
        raise TypeError("missing angular solve arguments: " + ", ".join(missing))
    provider = SimpleNamespace(
        mu=kwargs.pop("mu"),
        precision=kwargs.pop("precision", DEFAULT_PRECISION),
        weights=SimpleNamespace(
            A_mu=kwargs.pop("A_mu"),
            B_mu=kwargs.pop("B_mu"),
            D_mu=kwargs.pop("D_mu"),
            beta_inf=kwargs.pop("beta_inf", None),
        ),
    )
    return solve_coupled_angular(provider, **kwargs)


solve_angular_coefficients = solve_angular_bumps
coupled_angular_coefficients = solve_angular_bumps
solve_coupled_angular_from_targets = solve_coupled_angular


def _fixture_provider() -> Any:
    """Return a provider-shaped fixture with independent, known weights."""

    return SimpleNamespace(
        mu="0.01",
        precision=220,
        weights=SimpleNamespace(
            A_mu="0.73",
            B_mu="1.17",
            D_mu="0.29",
            beta_inf="2.4",
        ),
    )


def _fixture_targets(provider: Any, d1: Any, d2: Any, d1Z: Any, d2Z: Any) -> dict[str, mp.mpf]:
    """Generate exact row targets from prescribed coefficients and tangents."""

    dps, mu, A, B, D, _ = _provider_parameters(provider, None)
    with mp.workdps(dps + 40):
        one_minus_mu = 1 - mu
        one_plus_2mu = 1 + 2 * mu
        x1 = _mp(d1, "d1")
        x2 = _mp(d2, "d2")
        x1Z = _mp(d1Z, "d1Z")
        x2Z = _mp(d2Z, "d2Z")
        eang1 = mp.exp(-3 * one_minus_mu)
        eang2 = mp.exp(-one_minus_mu)
        ep1 = mp.exp(3 * one_plus_2mu)
        ep2 = mp.exp(one_plus_2mu)
        r = A * (eang1 * x1 + eang2 * x2)
        rZ = A * (eang1 * x1Z + eang2 * x2Z)
        s = B * (ep1 * x1 + ep2 * x2) + D / 2 * (ep1 * x1 * x1 + ep2 * x2 * x2)
        sZ = (
            B * (ep1 * x1Z + ep2 * x2Z)
            + D * (ep1 * x1 * x1Z + ep2 * x2 * x2Z)
        )
        return {"r": r, "s": s, "rZ": rZ, "sZ": sZ}


def run_fixture() -> dict[str, Any]:
    """Run independent signed, zero, and tiny coefficient replay cases."""

    provider = _fixture_provider()
    cases = (
        ("signed_r_positive_s_negative", "-.12", ".08", ".03", "-.02"),
        ("signed_r_negative_s_positive", ".07", "-.13", "-.04", ".05"),
        ("both_zero_with_nonzero_tangent", "0", "0", ".11", "-.04"),
        ("tiny_signed_coefficients", "-1e-90", "2e-90", "3e-91", "-4e-91"),
    )
    reports: list[dict[str, Any]] = []
    for name, d1, d2, d1Z, d2Z in cases:
        targets = _fixture_targets(provider, d1, d2, d1Z, d2Z)
        solution = solve_coupled_angular(provider, **targets)
        with mp.workdps(provider.precision):
            expected_d1 = _mp(d1)
            expected_d2 = _mp(d2)
            expected_d1Z = _mp(d1Z)
            expected_d2Z = _mp(d2Z)
            value_errors = {
                "d1": abs(solution.d1 - expected_d1),
                "d2": abs(solution.d2 - expected_d2),
                "d1_Z": abs(solution.d1_Z - expected_d1Z),
                "d2_Z": abs(solution.d2_Z - expected_d2Z),
            }
            max_error = max(value_errors.values())
            if max_error > mp.mpf("1e-120"):
                raise AssertionError(f"fixture {name} coefficient/tangent replay error {max_error}")
            if any(abs(value) > mp.mpf("1e-120") for value in solution.row_residuals):
                raise AssertionError(f"fixture {name} direct row replay failed")
            if any(abs(value) > mp.mpf("1e-120") for value in solution.derivative_residuals):
                raise AssertionError(f"fixture {name} derivative row replay failed")
            reports.append(
                {
                    "name": name,
                    "targets": {key: _signed(value, 80) for key, value in targets.items()},
                    "solution": solution.as_dict(80),
                    "max_coefficient_or_tangent_error": mp.nstr(max_error, 20),
                    "zero_target_case": all(value == 0 for value in (targets["r"], targets["s"])),
                    "signed_target_case": any(value < 0 for value in (targets["r"], targets["s"])),
                    "tiny_target_case": max(abs(targets["r"]), abs(targets["s"])) < mp.mpf("1e-80"),
                }
            )
    return {
        "source": SOURCE,
        "source_equation": SOURCE_EQUATION,
        "algebra_version": ALGEBRA_VERSION,
        "cases": reports,
        "scope": "Independent algebra fixture only; no actual field installation, cone, or global closure.",
    }


_fixture = run_fixture
run = run_fixture


if __name__ == "__main__":
    report = run_fixture()
    print(json.dumps({"cases": [case["name"] for case in report["cases"]]}, indent=2), flush=True)
