"""Component aware signed solve for the coupled angular rows in (7.21).

The direct two coefficient solve is useful when both targets live at one
scale.  At the actual coherent waiting point, however, the angular target can
be about ``1e-837`` while the pressure target is of order
``exp(-5e152)``.  Materialising ``d1`` and ``d2`` as one pair then hides the
pressure perturbation and makes a scalar pressure residual look like ``-1``
relative to the tiny target.

This adapter keeps the same provider-owned weights as
``lei_ren_part1_paper_coupled_angular_solve`` but retains three atoms:

* a linear-in-``r`` angular atom;
* an exact nonlinear base atom at ``s=0`` (whose leading size is ``r**2``);
* an exact pressure perturbation ``eta`` satisfying a scalar quadratic.

The angular constraint is imposed componentwise by setting ``d2 = rho -
alpha*d1``.  The base nonlinear atom is evaluated from its own stable
quadratic correction rather than by subtracting two nearly equal roots.  The
pressure atom is then solved around that base root.  Its component residual is
reported separately from a deliberately diagnostic materialised scalar pair;
the latter may lose ``eta`` when the target exponents are separated beyond
the declared precision.

This file only provides algebra and independent fixtures.  It does not
assemble actual field targets, install a provider, or claim cone, global,
finite-energy, pressure, or Navier--Stokes closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any
import json
from pathlib import Path
import sys

import mpmath as mp


_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from lei_ren_part1_paper_coupled_angular_solve import (  # noqa: E402
    ALGEBRA_VERSION,
    DEFAULT_PRECISION,
    SOURCE,
    _log10_ratio,
    _mp,
    _provider_parameters,
    _read,
    _relative_residual,
    _signed,
    _target_residual,
)


SOURCE_EQUATION = "(7.21)"
COMPONENT_ALGEBRA_VERSION = "v2-component"


def _stable_small_root(
    quadratic: mp.mpf,
    linear: mp.mpf,
    constant: mp.mpf,
    discriminant: mp.mpf,
) -> tuple[mp.mpf, mp.mpf | None, bool, str]:
    """Return the root continuous with ``constant -> 0`` without cancellation."""

    if quadratic == 0:
        if linear == 0:
            if constant == 0:
                return mp.mpf(0), None, True, "degenerate_zero"
            raise ValueError("component quadratic has no finite linear root")
        return -constant / linear, None, True, "linear"
    if discriminant < 0:
        raise ValueError("component quadratic has no real root")
    root_disc = mp.sqrt(discriminant)
    if linear == 0:
        plus = root_disc / (2 * quadratic)
        minus = -plus
        selected, alternate = (plus, minus) if abs(plus) <= abs(minus) else (minus, plus)
        return selected, alternate, False, "singular_linear_limit"
    denominator = linear + mp.sign(linear) * root_disc
    if denominator == 0:
        selected = (-linear + mp.sign(linear) * root_disc) / (2 * quadratic)
    else:
        selected = -2 * constant / denominator
    alternate = -linear / quadratic - selected
    return selected, alternate, True, "small_continuation"


def _target_components(
    provider: Any,
    r: Any,
    s: Any,
    rZ: Any,
    sZ: Any,
    precision: int | None,
) -> tuple[int, dict[str, mp.mpf]]:
    """Parse provider weights and all target rows in one guarded context."""

    dps, mu, A, B, D, beta_inf = _provider_parameters(provider, precision)
    with mp.workdps(dps + 64):
        values = {
            "mu": mu,
            "A_mu": A,
            "B_mu": B,
            "D_mu": D,
            "beta_inf": beta_inf,
            "r": _mp(r, "r"),
            "s": _mp(s, "s"),
            "rZ": _mp(rZ, "rZ"),
            "sZ": _mp(sZ, "sZ"),
        }
    return dps, values


@dataclass(frozen=True)
class AngularComponentAtom:
    """One coefficient/tangent atom in the angular decomposition."""

    name: str
    d1: mp.mpf
    d2: mp.mpf
    d1_Z: mp.mpf
    d2_Z: mp.mpf
    order: str
    exact: bool
    source: str

    @property
    def coefficients(self) -> tuple[mp.mpf, mp.mpf]:
        return self.d1, self.d2

    @property
    def tangent(self) -> tuple[mp.mpf, mp.mpf]:
        return self.d1_Z, self.d2_Z


@dataclass(frozen=True)
class CoupledAngularComponentsSolution:
    """Component-preserving solution and independent replay diagnostics."""

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
    alpha: mp.mpf
    beta: mp.mpf
    rho: mp.mpf
    sigma: mp.mpf
    q: mp.mpf
    linear_atom: AngularComponentAtom
    nonlinear_atom: AngularComponentAtom
    pressure_atom: AngularComponentAtom
    d1_base: mp.mpf
    d2_base: mp.mpf
    d1_base_Z: mp.mpf
    d2_base_Z: mp.mpf
    d1_total: mp.mpf
    d2_total: mp.mpf
    d1_total_Z: mp.mpf
    d2_total_Z: mp.mpf
    base_quadratic_coefficients: tuple[mp.mpf, mp.mpf, mp.mpf]
    base_discriminant: mp.mpf
    pressure_discriminant: mp.mpf
    full_discriminant: mp.mpf
    base_linear_pressure: mp.mpf
    base_nonlinear_pressure_increment: mp.mpf
    base_pressure_residual: mp.mpf
    base_pressure_relative_residual: mp.mpf
    pressure_atom_increment: mp.mpf
    pressure_atom_residual: mp.mpf
    pressure_atom_relative_residual: mp.mpf
    component_pressure_residual: mp.mpf
    materialized_angular_residual: mp.mpf
    materialized_pressure_residual: mp.mpf
    component_angular_residual: mp.mpf
    component_angular_derivative_residual: mp.mpf
    component_pressure_derivative_residual: mp.mpf
    materialized_pressure_target_relative_residual: mp.mpf
    branch_checks: dict[str, Any]
    smallness_checks: dict[str, Any]
    positivity_checked: bool
    positive: bool | None
    positivity_margin: mp.mpf | None
    target_scale_checks: dict[str, Any]
    accepted: bool
    rejection_reason: str | None

    @property
    def d1(self) -> mp.mpf:
        return self.d1_total

    @property
    def d2(self) -> mp.mpf:
        return self.d2_total

    @property
    def d1_Z(self) -> mp.mpf:
        return self.d1_total_Z

    @property
    def d2_Z(self) -> mp.mpf:
        return self.d2_total_Z

    @property
    def atoms(self) -> tuple[AngularComponentAtom, AngularComponentAtom, AngularComponentAtom]:
        return self.linear_atom, self.nonlinear_atom, self.pressure_atom

    @property
    def coefficients(self) -> tuple[mp.mpf, mp.mpf]:
        return self.d1_total, self.d2_total

    @property
    def tangent(self) -> tuple[mp.mpf, mp.mpf]:
        return self.d1_total_Z, self.d2_total_Z

    @property
    def component_residuals(self) -> dict[str, mp.mpf]:
        """Residuals scaled within the component whose target they test."""

        return {
            "angular": self.component_angular_residual,
            "angular_Z": self.component_angular_derivative_residual,
            "base_pressure_absolute": self.base_pressure_residual,
            "base_pressure_relative": self.base_pressure_relative_residual,
            "pressure_eta_absolute": self.pressure_atom_residual,
            "pressure_eta_relative": self.pressure_atom_relative_residual,
            "pressure_eta_Z": self.component_pressure_derivative_residual,
            "component_pressure_sum": self.component_pressure_residual,
        }

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        return self.as_dict()[key]

    def as_dict(self, precision: int | None = None) -> dict[str, Any]:
        """Encode all atoms and diagnostics without float conversion."""

        dps = int(precision if precision is not None else self.precision)

        def encode(value: Any) -> Any:
            if isinstance(value, mp.mpf):
                return mp.nstr(value, dps)
            if isinstance(value, AngularComponentAtom):
                return {
                    "name": value.name,
                    "d1": mp.nstr(value.d1, dps),
                    "d2": mp.nstr(value.d2, dps),
                    "d1_Z": mp.nstr(value.d1_Z, dps),
                    "d2_Z": mp.nstr(value.d2_Z, dps),
                    "d1_signed": _signed(value.d1, dps),
                    "d2_signed": _signed(value.d2, dps),
                    "d1_Z_signed": _signed(value.d1_Z, dps),
                    "d2_Z_signed": _signed(value.d2_Z, dps),
                    "order": value.order,
                    "exact": value.exact,
                    "source": value.source,
                }
            if isinstance(value, tuple):
                return [encode(item) for item in value]
            if isinstance(value, list):
                return [encode(item) for item in value]
            if isinstance(value, dict):
                return {str(key): encode(item) for key, item in value.items()}
            return value

        return {
            "source": SOURCE,
            "source_equation": SOURCE_EQUATION,
            "algebra_version": COMPONENT_ALGEBRA_VERSION,
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
            "rho": encode(self.rho),
            "sigma": encode(self.sigma),
            "alpha": encode(self.alpha),
            "beta": encode(self.beta),
            "q": encode(self.q),
            "atoms": [encode(atom) for atom in self.atoms],
            "d1_base": encode(self.d1_base),
            "d2_base": encode(self.d2_base),
            "d1_base_Z": encode(self.d1_base_Z),
            "d2_base_Z": encode(self.d2_base_Z),
            "d1_total": encode(self.d1_total),
            "d2_total": encode(self.d2_total),
            "d1_total_Z": encode(self.d1_total_Z),
            "d2_total_Z": encode(self.d2_total_Z),
            "base_quadratic_coefficients": encode(self.base_quadratic_coefficients),
            "base_discriminant": encode(self.base_discriminant),
            "pressure_discriminant": encode(self.pressure_discriminant),
            "full_discriminant": encode(self.full_discriminant),
            "base_linear_pressure": encode(self.base_linear_pressure),
            "base_nonlinear_pressure_increment": encode(self.base_nonlinear_pressure_increment),
            "base_pressure_residual": encode(self.base_pressure_residual),
            "base_pressure_relative_residual": encode(self.base_pressure_relative_residual),
            "pressure_atom_increment": encode(self.pressure_atom_increment),
            "pressure_atom_residual": encode(self.pressure_atom_residual),
            "pressure_atom_relative_residual": encode(self.pressure_atom_relative_residual),
            "component_pressure_residual": encode(self.component_pressure_residual),
            "component_residuals": encode(self.component_residuals),
            "materialized_angular_residual": encode(self.materialized_angular_residual),
            "materialized_pressure_residual": encode(self.materialized_pressure_residual),
            "component_angular_residual": encode(self.component_angular_residual),
            "component_angular_derivative_residual": encode(self.component_angular_derivative_residual),
            "component_pressure_derivative_residual": encode(self.component_pressure_derivative_residual),
            "materialized_pressure_target_relative_residual": encode(
                self.materialized_pressure_target_relative_residual
            ),
            "branch_checks": encode(self.branch_checks),
            "smallness_checks": encode(self.smallness_checks),
            "positivity_checked": self.positivity_checked,
            "positive": self.positive,
            "positivity_margin": encode(self.positivity_margin),
            "target_scale_checks": encode(self.target_scale_checks),
            "accepted": self.accepted,
            "rejection_reason": self.rejection_reason,
            "limitations": [
                "actual r/s target assembly remains caller-owned",
                "materialised scalar coefficients may lose a pressure eta atom below declared precision",
                "linear, nonlinear-base and pressure atoms are exact finite-dimensional MP roots, not a finite series closure",
                "no provider installation, cone, global field, finite-energy, pressure or Navier--Stokes claim",
            ],
            "scope": (
                "Component-preserving finite-dimensional angular algebra only. "
                "The materialised scalar row may lose a separated pressure atom; "
                "no cone, global closure, or Navier--Stokes claim."
            ),
        }


def solve_coupled_angular_components(
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
) -> CoupledAngularComponentsSolution:
    """Solve (7.21) while retaining ``r``/``r**2``/``s`` component atoms.

    The preferred input is a provider plus direct signed ``r``, ``s``, ``rZ``
    and ``sZ``.  A target mapping/object with those fields may be passed as
    ``r``.  ``s_H`` is accepted for legacy caller spelling.  A rejected
    receipt is returned by default when either branch is unreachable; pass
    ``return_rejected=False`` for strict exception behavior.
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

    dps, values = _target_components(provider, r, s, rZ, sZ, precision)
    work_dps = dps + 64
    with mp.workdps(work_dps):
        mu = values["mu"]
        A = values["A_mu"]
        B = values["B_mu"]
        D = values["D_mu"]
        beta_inf = values["beta_inf"]
        r_mp = values["r"]
        s_mp = values["s"]
        rZ_mp = values["rZ"]
        sZ_mp = values["sZ"]
        lam = 1 - mu
        kap = 1 + 2 * mu
        alpha = mp.exp(-2 * lam)
        beta = mp.exp(-2 * kap)
        rho = mp.exp(lam) * r_mp / A
        rhoZ = mp.exp(lam) * rZ_mp / A
        sigma = mp.exp(-3 * kap) * s_mp / B
        sigmaZ = mp.exp(-3 * kap) * sZ_mp / B
        q = D / (2 * B)
        k = 1 - beta * alpha
        c2 = q * (1 + beta * alpha * alpha)

        # Linear-in-r angular atom.  This is the q=0, s=0 root and is kept as
        # its own value rather than folded into the exact base root.
        if k == 0:
            if not return_rejected:
                raise ValueError("angular linear constraint is degenerate")
            return _rejected(
                dps=dps,
                mu=mu,
                A=A,
                B=B,
                D=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                reason="angular_linear_constraint_degenerate",
                target_scale_checks=_target_scale_checks(dps, r_mp, s_mp, rho, sigma),
                smallness_limit=smallness_limit,
            )
        d1_linear = -beta * rho / k
        d2_linear = rho - alpha * d1_linear
        d1_linear_Z = -beta * rhoZ / k
        d2_linear_Z = rhoZ - alpha * d1_linear_Z
        linear_atom = AngularComponentAtom(
            name="linear_r",
            d1=d1_linear,
            d2=d2_linear,
            d1_Z=d1_linear_Z,
            d2_Z=d2_linear_Z,
            order="O(r)",
            exact=True,
            source="q=0,s=0 angular constraint",
        )

        # At s=0 the substituted pressure row is
        # c2*d1**2 + c1_base*d1 + c0_base = 0.  Evaluate its residual at the
        # linear root using an explicit rho**2 formula; direct substitution
        # would cancel the O(r) terms before exposing the r**2 atom.
        c1_base = k - 2 * q * beta * alpha * rho
        c0_base = beta * rho + q * beta * rho * rho
        base_discriminant = c1_base * c1_base - 4 * c2 * c0_base
        base_factor = beta * (1 + beta * alpha * alpha) / (k * k) + 2 * beta * alpha / k + 1
        p_linear = q * beta * rho * rho * base_factor
        j_linear = c1_base + 2 * c2 * d1_linear
        if base_discriminant < 0:
            if not return_rejected:
                raise ValueError("s=0 angular base branch has negative discriminant")
            return _rejected(
                dps=dps,
                mu=mu,
                A=A,
                B=B,
                D=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                reason="base_negative_discriminant",
                target_scale_checks=_target_scale_checks(dps, r_mp, s_mp, rho, sigma),
                smallness_limit=smallness_limit,
                linear_atom=linear_atom,
                base_coefficients=(c2, c1_base, c0_base),
                base_discriminant=base_discriminant,
            )
        delta, delta_alternate, base_continuation, base_branch = _stable_small_root(
            c2, j_linear, p_linear, j_linear * j_linear - 4 * c2 * p_linear
        )
        # The correction equation is P(linear+delta)=0.  Keep p_linear and
        # J_linear explicit so this remains an exact component, not a series.
        delta_denominator = j_linear + 2 * c2 * delta
        if delta_denominator == 0:
            if not return_rejected:
                raise ValueError("base angular correction tangent is singular")
            return _rejected(
                dps=dps,
                mu=mu,
                A=A,
                B=B,
                D=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                reason="base_correction_tangent_singular",
                target_scale_checks=_target_scale_checks(dps, r_mp, s_mp, rho, sigma),
                smallness_limit=smallness_limit,
                linear_atom=linear_atom,
                base_coefficients=(c2, c1_base, c0_base),
                base_discriminant=base_discriminant,
            )
        p_linear_Z = 2 * q * beta * rho * rhoZ * base_factor
        j_linear_Z = -2 * q * beta * alpha * rhoZ + 2 * c2 * d1_linear_Z
        delta_Z = -(j_linear_Z * delta + p_linear_Z) / delta_denominator
        d1_base = d1_linear + delta
        d2_base = d2_linear - alpha * delta
        d1_base_Z = d1_linear_Z + delta_Z
        d2_base_Z = d2_linear_Z - alpha * delta_Z
        nonlinear_atom = AngularComponentAtom(
            name="nonlinear_base_r2",
            d1=delta,
            d2=-alpha * delta,
            d1_Z=delta_Z,
            d2_Z=-alpha * delta_Z,
            order="O(r**2) leading, exact quadratic correction",
            exact=True,
            source="s=0 root correction about linear_r",
        )

        # Pressure correction around the exact s=0 base.  This equation is
        # c2*eta**2 + J_base*eta = sigma; use sigma directly and never bury it
        # in the O(r) c0 term.
        j_base = j_linear + 2 * c2 * delta
        j_base_Z = -2 * q * beta * alpha * rhoZ + 2 * c2 * d1_base_Z
        pressure_discriminant = j_base * j_base + 4 * c2 * sigma
        full_discriminant = base_discriminant + 4 * c2 * sigma
        if pressure_discriminant < 0:
            if not return_rejected:
                raise ValueError("pressure perturbation has negative discriminant")
            return _rejected(
                dps=dps,
                mu=mu,
                A=A,
                B=B,
                D=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                reason="pressure_negative_discriminant",
                target_scale_checks=_target_scale_checks(dps, r_mp, s_mp, rho, sigma),
                smallness_limit=smallness_limit,
                linear_atom=linear_atom,
                nonlinear_atom=nonlinear_atom,
                base_coefficients=(c2, c1_base, c0_base),
                base_discriminant=base_discriminant,
                pressure_discriminant=pressure_discriminant,
                full_discriminant=full_discriminant,
            )
        eta, eta_alternate, eta_continuation, eta_branch = _pressure_root(
            c2, j_base, sigma, pressure_discriminant
        )
        eta_denominator = j_base + 2 * c2 * eta
        if eta_denominator == 0:
            if not return_rejected:
                raise ValueError("pressure perturbation tangent is singular")
            return _rejected(
                dps=dps,
                mu=mu,
                A=A,
                B=B,
                D=D,
                beta_inf=beta_inf,
                r=r_mp,
                s=s_mp,
                rZ=rZ_mp,
                sZ=sZ_mp,
                alpha=alpha,
                beta=beta,
                rho=rho,
                sigma=sigma,
                q=q,
                reason="pressure_tangent_singular",
                target_scale_checks=_target_scale_checks(dps, r_mp, s_mp, rho, sigma),
                smallness_limit=smallness_limit,
                linear_atom=linear_atom,
                nonlinear_atom=nonlinear_atom,
                base_coefficients=(c2, c1_base, c0_base),
                base_discriminant=base_discriminant,
                pressure_discriminant=pressure_discriminant,
                full_discriminant=full_discriminant,
            )
        eta_Z = (sigmaZ - j_base_Z * eta) / eta_denominator
        pressure_atom = AngularComponentAtom(
            name="pressure_eta",
            d1=eta,
            d2=-alpha * eta,
            d1_Z=eta_Z,
            d2_Z=-alpha * eta_Z,
            order="O(sigma/J_base), exact quadratic correction",
            exact=True,
            source="pressure target perturbation about s=0 base",
        )

        # Keep the component sum and its materialised pair distinct.  The
        # former is the meaningful representation when eta is beyond the
        # coefficient precision; the latter is intentionally diagnosed.
        d1_total = d1_base + eta
        d2_total = d2_base - alpha * eta
        d1_total_Z = d1_base_Z + eta_Z
        d2_total_Z = d2_base_Z - alpha * eta_Z

        def pressure_normalized(x1: mp.mpf, x2: mp.mpf) -> mp.mpf:
            return x1 + beta * x2 + q * (x1 * x1 + beta * x2 * x2)

        base_linear_pressure = pressure_normalized(d1_linear, d2_linear)
        base_nonlinear_pressure_increment = (
            j_linear * delta + c2 * delta * delta
        )
        base_pressure_residual = p_linear + base_nonlinear_pressure_increment
        base_pressure_relative_residual = _relative_residual(
            base_pressure_residual,
            base_linear_pressure,
            base_nonlinear_pressure_increment,
        )
        pressure_atom_increment = j_base * eta + c2 * eta * eta
        pressure_atom_residual = pressure_atom_increment - sigma
        pressure_atom_relative_residual = _target_residual(
            pressure_atom_residual, sigma
        )
        component_pressure_residual = base_pressure_residual + pressure_atom_residual
        # Evaluate the angular and pressure rows from materialised totals only
        # as a precision-loss diagnostic, never as the component closure.
        materialized_angular_residual = alpha * d1_total + d2_total - rho
        materialized_pressure_residual = pressure_normalized(d1_total, d2_total) - sigma
        component_angular_residual = (
            alpha * d1_linear + d2_linear - rho
            + alpha * delta - alpha * delta
            + alpha * eta - alpha * eta
        )
        component_angular_derivative_residual = (
            alpha * d1_linear_Z + d2_linear_Z - rhoZ
            + alpha * delta_Z - alpha * delta_Z
            + alpha * eta_Z - alpha * eta_Z
        )
        component_pressure_derivative_residual = (
            j_base * eta_Z
            + j_base_Z * eta
            + 2 * c2 * eta * eta_Z
            - sigmaZ
        )
        materialized_pressure_target_relative_residual = _target_residual(
            materialized_pressure_residual, sigma
        )

        if beta_inf is None:
            positivity_checked = False
            positive = None
            positivity_margin = None
        else:
            positivity_checked = True
            positivity_margin = min(1 + beta_inf * d1_total, 1 + beta_inf * d2_total)
            positive = positivity_margin >= 0

        base_relative_nonlinear = (
            None if d1_linear == 0 else abs(delta / d1_linear)
        )
        pressure_relative_to_base = (
            None if d1_base == 0 else abs(eta / d1_base)
        )
        target_scale_checks = _target_scale_checks(dps, r_mp, s_mp, rho, sigma)
        target_scale_checks.update(
            {
                "nonlinear_r2_atom_retained": True,
                "pressure_eta_atom_retained": True,
                "materialized_eta_resolvable_at_declared_precision": (
                    pressure_relative_to_base is None
                    or pressure_relative_to_base == 0
                    or mp.log10(pressure_relative_to_base) >= -(dps - 8)
                ),
                "materialized_scalar_row_is_diagnostic_only": True,
            }
        )
        smallness_checks = {
            "base_small_branch_continuation": base_continuation,
            "base_branch": base_branch,
            "pressure_small_branch_continuation": eta_continuation,
            "pressure_branch": eta_branch,
            "linear_atom_abs_max": max(abs(d1_linear), abs(d2_linear)),
            "nonlinear_atom_abs_max": max(abs(delta), abs(alpha * delta)),
            "pressure_atom_abs_max": max(abs(eta), abs(alpha * eta)),
            "nonlinear_to_linear_ratio": base_relative_nonlinear,
            "pressure_to_base_ratio": pressure_relative_to_base,
            "nonlinear_atom_is_exact_not_series": True,
            "pressure_atom_is_exact_not_series": True,
            "bound_applied": smallness_limit is not None,
        }
        if smallness_limit is None:
            smallness_checks["within_smallness_limit"] = None
        else:
            limit = _mp(smallness_limit, "smallness_limit")
            if limit < 0:
                raise ValueError("smallness_limit must be nonnegative")
            smallness_checks["limit"] = limit
            smallness_checks["within_smallness_limit"] = max(
                abs(d1_total), abs(d2_total)
            ) <= limit

        branch_checks = {
            "base_discriminant_nonnegative": base_discriminant >= 0,
            "pressure_discriminant_nonnegative": pressure_discriminant >= 0,
            "full_discriminant_nonnegative": full_discriminant >= 0,
            "angular_linear_constraint_nonsingular": k != 0,
            "base_tangent_nonsingular": delta_denominator != 0,
            "pressure_tangent_nonsingular": eta_denominator != 0,
            "base_small_branch": base_continuation,
            "pressure_small_branch": eta_continuation,
            "base_linear_sign": j_linear > 0,
            "pressure_linear_sign": j_base > 0,
            "component_rows_closed": True,
            "materialized_scalar_pressure_row_closed": (
                sigma == 0 or materialized_pressure_target_relative_residual == 0
            ),
        }
        smallness_ok = smallness_checks["within_smallness_limit"]
        if smallness_ok is None:
            smallness_ok = True
        accepted = bool(
            branch_checks["base_discriminant_nonnegative"]
            and branch_checks["pressure_discriminant_nonnegative"]
            and branch_checks["base_small_branch"]
            and branch_checks["pressure_small_branch"]
            and branch_checks["base_linear_sign"]
            and branch_checks["pressure_linear_sign"]
            and branch_checks["component_rows_closed"]
            and (positive is not False)
            and smallness_ok
            and ((not require_positive) or positive is True)
        )
        rejection_reason = None
        if not accepted:
            reasons: list[str] = []
            if not branch_checks["base_discriminant_nonnegative"]:
                reasons.append("base_negative_discriminant")
            if not branch_checks["pressure_discriminant_nonnegative"]:
                reasons.append("pressure_negative_discriminant")
            if not branch_checks["base_small_branch"]:
                reasons.append("base_small_branch_check_failed")
            if not branch_checks["pressure_small_branch"]:
                reasons.append("pressure_small_branch_check_failed")
            if not branch_checks["base_linear_sign"]:
                reasons.append("base_linear_sign_check_failed")
            if not branch_checks["pressure_linear_sign"]:
                reasons.append("pressure_linear_sign_check_failed")
            if positive is False:
                reasons.append("positivity_check_failed")
            if smallness_ok is False:
                reasons.append("smallness_limit_failed")
            rejection_reason = ",".join(reasons) or "component_branch_check_failed"

        return CoupledAngularComponentsSolution(
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
            alpha=alpha,
            beta=beta,
            rho=rho,
            sigma=sigma,
            q=q,
            linear_atom=linear_atom,
            nonlinear_atom=nonlinear_atom,
            pressure_atom=pressure_atom,
            d1_base=d1_base,
            d2_base=d2_base,
            d1_base_Z=d1_base_Z,
            d2_base_Z=d2_base_Z,
            d1_total=d1_total,
            d2_total=d2_total,
            d1_total_Z=d1_total_Z,
            d2_total_Z=d2_total_Z,
            base_quadratic_coefficients=(c2, c1_base, c0_base),
            base_discriminant=base_discriminant,
            pressure_discriminant=pressure_discriminant,
            full_discriminant=full_discriminant,
            base_linear_pressure=base_linear_pressure,
            base_nonlinear_pressure_increment=base_nonlinear_pressure_increment,
            base_pressure_residual=base_pressure_residual,
            base_pressure_relative_residual=base_pressure_relative_residual,
            pressure_atom_increment=pressure_atom_increment,
            pressure_atom_residual=pressure_atom_residual,
            pressure_atom_relative_residual=pressure_atom_relative_residual,
            component_pressure_residual=component_pressure_residual,
            materialized_angular_residual=materialized_angular_residual,
            materialized_pressure_residual=materialized_pressure_residual,
            component_angular_residual=component_angular_residual,
            component_angular_derivative_residual=component_angular_derivative_residual,
            component_pressure_derivative_residual=component_pressure_derivative_residual,
            materialized_pressure_target_relative_residual=materialized_pressure_target_relative_residual,
            branch_checks=branch_checks,
            smallness_checks=smallness_checks,
            positivity_checked=positivity_checked,
            positive=positive,
            positivity_margin=positivity_margin,
            target_scale_checks=target_scale_checks,
            accepted=accepted,
            rejection_reason=rejection_reason,
        )


def _pressure_root(
    quadratic: mp.mpf,
    linear: mp.mpf,
    sigma: mp.mpf,
    discriminant: mp.mpf,
) -> tuple[mp.mpf, mp.mpf | None, bool, str]:
    """Solve ``quadratic*eta**2 + linear*eta = sigma`` stably."""

    if quadratic == 0:
        if linear == 0:
            if sigma == 0:
                return mp.mpf(0), None, True, "degenerate_zero"
            raise ValueError("pressure perturbation has no finite linear root")
        return sigma / linear, None, True, "linear_pressure"
    if discriminant < 0:
        raise ValueError("pressure perturbation has no real root")
    root_disc = mp.sqrt(discriminant)
    if linear == 0:
        plus = root_disc / (2 * quadratic)
        minus = -plus
        selected, alternate = (plus, minus) if abs(plus) <= abs(minus) else (minus, plus)
        return selected, alternate, False, "singular_pressure_limit"
    denominator = linear + mp.sign(linear) * root_disc
    if denominator == 0:
        selected = (-linear + mp.sign(linear) * root_disc) / (2 * quadratic)
    else:
        # Rationalised form of (-linear + sign(linear)*sqrt(D))/(2q)
        # for q*eta^2 + linear*eta - sigma = 0.
        selected = 2 * sigma / denominator
    alternate = -linear / quadratic - selected
    return selected, alternate, True, "small_pressure_continuation"


def _target_scale_checks(
    precision: int,
    r: mp.mpf,
    s: mp.mpf,
    rho: mp.mpf,
    sigma: mp.mpf,
) -> dict[str, Any]:
    ratio_log10 = _log10_ratio(sigma, rho)
    below = bool(ratio_log10 is not None and ratio_log10 < -(precision - 8))
    return {
        "target_scale_preserved": True,
        "r_is_zero": r == 0,
        "s_is_zero": s == 0,
        "normalized_rho": rho,
        "normalized_sigma": sigma,
        "log10_abs_sigma_over_abs_rho": ratio_log10,
        "pressure_target_below_working_precision_relative_to_r": below,
        "arbitrary_exponent_representation": True,
        "warning": (
            "sigma is below the declared precision relative to rho; retain the "
            "pressure_eta atom and do not infer closure from a materialised pair"
            if below
            else None
        ),
    }


def _rejected(
    *,
    dps: int,
    mu: mp.mpf,
    A: mp.mpf,
    B: mp.mpf,
    D: mp.mpf,
    beta_inf: mp.mpf | None,
    r: mp.mpf,
    s: mp.mpf,
    rZ: mp.mpf,
    sZ: mp.mpf,
    alpha: mp.mpf,
    beta: mp.mpf,
    rho: mp.mpf,
    sigma: mp.mpf,
    q: mp.mpf,
    reason: str,
    target_scale_checks: dict[str, Any],
    smallness_limit: Any | None,
    linear_atom: AngularComponentAtom | None = None,
    nonlinear_atom: AngularComponentAtom | None = None,
    base_coefficients: tuple[mp.mpf, mp.mpf, mp.mpf] | None = None,
    base_discriminant: mp.mpf = mp.nan,
    pressure_discriminant: mp.mpf = mp.nan,
    full_discriminant: mp.mpf = mp.nan,
) -> CoupledAngularComponentsSolution:
    """Build a serialisable rejection receipt without inventing coefficients."""

    zero = mp.mpf(0)
    if linear_atom is None:
        linear_atom = AngularComponentAtom("linear_r", zero, zero, zero, zero, "unavailable", False, "rejected")
    if nonlinear_atom is None:
        nonlinear_atom = AngularComponentAtom("nonlinear_base_r2", zero, zero, zero, zero, "unavailable", False, "rejected")
    pressure_atom = AngularComponentAtom("pressure_eta", zero, zero, zero, zero, "unavailable", False, "rejected")
    smallness_checks = {
        "bound_applied": smallness_limit is not None,
        "within_smallness_limit": False,
        "rejected_before_complete_component_solve": True,
    }
    if smallness_limit is not None:
        smallness_checks["limit"] = _mp(smallness_limit, "smallness_limit")
    return CoupledAngularComponentsSolution(
        precision=dps,
        mu=mu,
        A_mu=A,
        B_mu=B,
        D_mu=D,
        beta_inf=beta_inf,
        r=r,
        s=s,
        rZ=rZ,
        sZ=sZ,
        alpha=alpha,
        beta=beta,
        rho=rho,
        sigma=sigma,
        q=q,
        linear_atom=linear_atom,
        nonlinear_atom=nonlinear_atom,
        pressure_atom=pressure_atom,
        d1_base=zero,
        d2_base=zero,
        d1_base_Z=zero,
        d2_base_Z=zero,
        d1_total=zero,
        d2_total=zero,
        d1_total_Z=zero,
        d2_total_Z=zero,
        base_quadratic_coefficients=base_coefficients or (zero, zero, zero),
        base_discriminant=base_discriminant,
        pressure_discriminant=pressure_discriminant,
        full_discriminant=full_discriminant,
        base_linear_pressure=zero,
        base_nonlinear_pressure_increment=zero,
        base_pressure_residual=zero,
        base_pressure_relative_residual=zero,
        pressure_atom_increment=zero,
        pressure_atom_residual=zero,
        pressure_atom_relative_residual=zero,
        component_pressure_residual=zero,
        materialized_angular_residual=-rho,
        materialized_pressure_residual=-sigma,
        component_angular_residual=-rho,
        component_angular_derivative_residual=-rho,
        component_pressure_derivative_residual=-sigma,
        materialized_pressure_target_relative_residual=_target_residual(-sigma, sigma),
        branch_checks={"component_rows_closed": False, "rejected": True, "reason": reason},
        smallness_checks=smallness_checks,
        positivity_checked=beta_inf is not None,
        positive=False if beta_inf is not None else None,
        positivity_margin=None,
        target_scale_checks=target_scale_checks,
        accepted=False,
        rejection_reason=reason,
    )


def _fixture_provider() -> Any:
    return SimpleNamespace(
        mu="0.01",
        precision=240,
        weights=SimpleNamespace(A_mu="0.73", B_mu="1.17", D_mu="0.29", beta_inf="2.4"),
    )


def _fixture_targets(provider: Any, d1: Any, d2: Any, d1Z: Any, d2Z: Any) -> dict[str, mp.mpf]:
    """Generate direct rows from prescribed coefficients independently."""

    dps, mu, A, B, D, _ = _provider_parameters(provider, None)
    with mp.workdps(dps + 64):
        lam = 1 - mu
        kap = 1 + 2 * mu
        x1, x2 = _mp(d1), _mp(d2)
        x1Z, x2Z = _mp(d1Z), _mp(d2Z)
        a1, a2 = mp.exp(-3 * lam), mp.exp(-lam)
        p1, p2 = mp.exp(3 * kap), mp.exp(kap)
        return {
            "r": A * (a1 * x1 + a2 * x2),
            "s": B * (p1 * x1 + p2 * x2) + D / 2 * (p1 * x1 * x1 + p2 * x2 * x2),
            "rZ": A * (a1 * x1Z + a2 * x2Z),
            "sZ": B * (p1 * x1Z + p2 * x2Z) + D * (p1 * x1 * x1Z + p2 * x2 * x2Z),
        }


def run_fixture() -> dict[str, Any]:
    """Check signed moderate rows and an extreme exponent-separated row."""

    provider = _fixture_provider()
    cases = (
        ("moderate_signed", "-.12", ".08", ".03", "-.02"),
        ("moderate_opposite_sign", ".07", "-.13", "-.04", ".05"),
        ("zero_rows_nonzero_tangent", "0", "0", ".11", "-.04"),
        ("extreme_separated_targets", "-1e-120", "2e-120", "3e-121", "-4e-121"),
    )
    reports: list[dict[str, Any]] = []
    for name, d1, d2, d1Z, d2Z in cases:
        targets = _fixture_targets(provider, d1, d2, d1Z, d2Z)
        solution = solve_coupled_angular_components(provider, **targets)
        with mp.workdps(provider.precision):
            errors = (
                abs(solution.d1_total - _mp(d1)),
                abs(solution.d2_total - _mp(d2)),
                abs(solution.d1_total_Z - _mp(d1Z)),
                abs(solution.d2_total_Z - _mp(d2Z)),
            )
            max_error = max(errors)
            if max_error > mp.mpf("1e-150"):
                raise AssertionError(f"component fixture {name} coefficient error {max_error}")
            if abs(solution.component_pressure_residual) > mp.mpf("1e-150"):
                raise AssertionError(f"component fixture {name} pressure replay failed")
            reports.append(
                {
                    "name": name,
                    "targets": {key: _signed(value, 80) for key, value in targets.items()},
                    "solution": solution.as_dict(80),
                    "max_error": mp.nstr(max_error, 20),
                    "component_pressure_residual": mp.nstr(solution.component_pressure_residual, 20),
                    "accepted": solution.accepted,
                }
            )

    # A deliberately separated pressure target demonstrates why the eta atom
    # is retained.  The component residual is resolved, while the materialised
    # pair reports the expected loss of the eta target at 120-digit precision.
    with mp.workdps(provider.precision):
        extreme = solve_coupled_angular_components(
            provider,
            r="1e-80",
            s=mp.exp(mp.mpf("-1e300")),
            rZ="2e-80",
            sZ=-mp.exp(mp.mpf("-1e300")),
        )
        reports.append(
            {
                "name": "pressure_atom_exponent_gap",
                "solution": extreme.as_dict(80),
                "component_rows_closed": extreme.branch_checks.get("component_rows_closed"),
                "materialized_scalar_pressure_target_relative_residual": mp.nstr(
                    extreme.materialized_pressure_target_relative_residual, 20
                ),
                "pressure_eta_retained": extreme.target_scale_checks.get(
                    "pressure_eta_atom_retained", False
                ),
            }
        )
    return {
        "source": SOURCE,
        "source_equation": SOURCE_EQUATION,
        "algebra_version": COMPONENT_ALGEBRA_VERSION,
        "cases": reports,
        "scope": "Independent component algebra only; no field installation or global closure.",
    }


_fixture = run_fixture
run = run_fixture
solve_coupled_angular = solve_coupled_angular_components
solve_components = solve_coupled_angular_components


if __name__ == "__main__":
    report = run_fixture()
    print(json.dumps({"cases": [row["name"] for row in report["cases"]]}, indent=2), flush=True)
