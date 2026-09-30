"""Continuous arbitrary-precision incoming axial primitive.

This is a small, reusable source-coordinate provider for the incoming axial
factor in the Lei--Ren Part I source schedule.  Put

``y = log(R) - log(R_ref)``

and define ``c(y)`` by the source flat transition

``c(y) = 1`` for ``y <= 1``,
``c(y) = 1 - sigma(log(y) / Md)`` for ``1 < y < exp(Md)``,
``c(y) = 0`` for ``y >= exp(Md)``.

The same continuous ``sigma`` switch and derivative used by
``continuous_axial_pulse.py`` are reused directly when that sibling is
importable, with an exact standalone fallback.  Point values and cumulative
primitives therefore use one definition, while this module does not depend
on the float-backed outer schedule.

The numerical quadrature is arbitrary precision ``mpmath`` quadrature.  It
is a numerical diagnostic/provider layer: it does not provide a quadrature
enclosure, a global finite-energy result, or recursive source certification.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mpmath as mp


try:  # Reuse the registered pulse switch whenever this folder is on sys.path.
    from lei_ren_part1_paper_continuous_axial_pulse import (
        ContinuousAxialPulse as _SharedPulse,
    )
except ImportError:  # Package imports from the repository root use this path.
    try:
        from experiments.root_st073.lei_ren_part1_paper_continuous_axial_pulse import (
            ContinuousAxialPulse as _SharedPulse,
        )
    except ImportError:  # pragma: no cover - standalone fallback only.
        _SharedPulse = None


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Part I Section 7; axial source cutoff from the Section 6 schedule"


def _as_mpf(value: Any) -> mp.mpf:
    """Convert a scalar without routing it through binary64 arithmetic."""

    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(value)


def _exp_integral(power: mp.mpf, left: mp.mpf, right: mp.mpf) -> mp.mpf:
    """Return ``integral_left^right exp(power*y) dy`` stably."""

    if right == left:
        return mp.mpf("0")
    if power == 0:
        return right - left
    return mp.exp(power * left) * mp.expm1(power * (right - left)) / power


class ContinuousIncomingAxial:
    """Point and cumulative incoming axial data at selected MP precision.

    Parameters are retained as mpmath values.  In particular, ``Z`` is
    converted only with :func:`mpmath.mpf` in each call and is never coerced
    to a float or used as a float cache key.  ``Md`` may be any positive
    arbitrary-precision value; the shared source candidate uses ``Md = .5``.
    """

    def __init__(
        self,
        Md: Any = ".5",
        *,
        precision: int = 160,
        quadrature_order: int = 48,
    ) -> None:
        try:
            precision_index = int(precision)
        except (TypeError, ValueError) as exc:
            raise TypeError("precision must be an integer") from exc
        if precision_index < 50:
            raise ValueError("precision must be at least 50 decimal digits")
        try:
            order_index = int(quadrature_order)
        except (TypeError, ValueError) as exc:
            raise TypeError("quadrature_order must be an integer") from exc
        if order_index < 4:
            raise ValueError("quadrature_order must be at least 4")

        self.precision = precision_index
        self.quadrature_order = order_index
        with mp.workdps(self.precision + 20):
            self.Md = _as_mpf(Md)
            if not mp.isfinite(self.Md) or self.Md <= 0:
                raise ValueError("Md must be finite and positive")
            self.cutoff_y = mp.exp(self.Md)

    @staticmethod
    def sigma_pair(x: Any) -> tuple[mp.mpf, mp.mpf]:
        """Return the shared flat switch ``sigma(x)`` and ``sigma'(x)``.

        The phase-stable form is exactly the one used by
        ``ContinuousAxialPulse.sigma_pair``.  The endpoint values are flat by
        definition.
        """

        x = _as_mpf(x)
        if _SharedPulse is not None:
            value, derivative = _SharedPulse.sigma_pair(x)
            return _as_mpf(value), _as_mpf(derivative)
        if x <= 0:
            return mp.mpf(0), mp.mpf(0)
        if x >= 1:
            return mp.mpf(1), mp.mpf(0)
        phase = 1 / x**2 - 1 / (1 - x) ** 2
        if phase >= 0:
            tiny = mp.exp(-phase)
            value = tiny / (1 + tiny)
        else:
            tiny = mp.exp(phase)
            value = 1 / (1 + tiny)
        product = tiny / (1 + tiny) ** 2
        derivative = product * (2 / x**3 + 2 / (1 - x) ** 3)
        return value, derivative

    def _dps(self) -> int:
        # A small guard protects endpoint comparisons and transformed
        # quadrature while preserving the requested instance precision.
        return self.precision + 12

    def _cutoff(self) -> mp.mpf:
        return self.cutoff_y

    def c(self, y: Any) -> mp.mpf:
        """Return the continuous source cutoff ``c(y)``."""

        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            cutoff = self._cutoff()
            if y_mp <= 1:
                return mp.mpf(1)
            if y_mp >= cutoff:
                return mp.mpf(0)
            argument = mp.log(y_mp) / self.Md
            # Reflection avoids losing the positive tail when sigma(argument)
            # rounds to one at the active precision.
            reflected, _ = self.sigma_pair(1 - argument)
            return reflected

    cutoff = c

    def c_y(self, y: Any) -> mp.mpf:
        """Return the analytic derivative ``dc/dy`` (flat at both joins)."""

        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            cutoff = self._cutoff()
            if y_mp <= 1 or y_mp >= cutoff:
                return mp.mpf(0)
            argument = mp.log(y_mp) / self.Md
            _, sigma_derivative = self.sigma_pair(argument)
            return -sigma_derivative / (self.Md * y_mp)

    cutoff_y_derivative = c_y

    def Uz(self, y: Any, Z: Any) -> mp.mpf:
        """Return ``U^z(y,Z) = 4 Z c(y)``."""

        with mp.workdps(self._dps()):
            return 4 * _as_mpf(Z) * self.c(y)

    def Uz_y(self, y: Any, Z: Any) -> mp.mpf:
        """Return the analytic source derivative ``partial_y U^z``."""

        with mp.workdps(self._dps()):
            return 4 * _as_mpf(Z) * self.c_y(y)

    def Uz_Z(self, y: Any, Z: Any | None = None) -> mp.mpf:
        """Return ``partial_Z U^z = 4 c(y)``.

        ``Z`` is accepted as an unused optional argument so callers that pass
        the same ``(y,Z)`` point tuple to every jet method can reuse it.
        """

        del Z
        with mp.workdps(self._dps()):
            return 4 * self.c(y)

    def values(self, y: Any, Z: Any) -> dict[str, mp.mpf]:
        """Return a serializable point-jet mapping for profile installers."""

        with mp.workdps(self._dps()):
            return {
                "c": self.c(y),
                "c_y": self.c_y(y),
                "Uz": self.Uz(y, Z),
                "Uz_y": self.Uz_y(y, Z),
                "Uz_Z": self.Uz_Z(y, Z),
            }

    def _transition_integral(
        self,
        power: mp.mpf,
        upper: mp.mpf,
        cutoff_power: mp.mpf,
        order: int,
    ) -> mp.mpf:
        """Integrate the finite transition after the exact unit plateau.

        The change of variables ``y = exp(Md*t)`` puts the shared sigma
        switch on ``t in [0,1]`` and keeps partial transition endpoints
        explicit.  ``mp.quad`` remains the only numerical integration layer.
        """

        cutoff = self._cutoff()
        if upper <= 1 or cutoff <= 1:
            return mp.mpf(0)
        endpoint = min(upper, cutoff)
        if endpoint <= 1:
            return mp.mpf(0)
        t_endpoint = mp.log(endpoint) / self.Md
        t_endpoint = min(mp.mpf(1), max(mp.mpf(0), t_endpoint))
        if t_endpoint == 0:
            return mp.mpf(0)

        def integrand(t: mp.mpf) -> mp.mpf:
            y_value = mp.exp(self.Md * t)
            # The reflected call is algebraically 1-sigma(t), but preserves
            # the positive cutoff tail near t=1.
            cutoff_value, _ = self.sigma_pair(1 - t)
            if cutoff_power == 0:
                cutoff_factor = mp.mpf(1)
            else:
                cutoff_factor = cutoff_value**cutoff_power
            return (
                mp.exp(power * y_value)
                * cutoff_factor
                * self.Md
                * y_value
            )

        # Splitting at fixed t knots makes order refinements meaningful and
        # handles the flat endpoint layers without asking one rule to cross
        # every shape scale at once.
        knots = [mp.mpf(0)]
        for fraction in (mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75")):
            if fraction < t_endpoint:
                knots.append(fraction)
        knots.append(t_endpoint)
        return mp.quad(
            integrand,
            knots,
            method="tanh-sinh",
            maxdegree=order,
        )

    def weighted_integral(
        self,
        power: Any = 1,
        upper: Any | None = None,
        cutoff_power: Any = 1,
        *,
        order: int | None = None,
    ) -> mp.mpf:
        """Return ``integral_0^upper exp(power*y)c(y)^cutoff_power dy``.

        The unit plateau is integrated analytically.  Only the finite smooth
        transition is sent to mpmath quadrature.  ``upper=None`` means the
        terminal source endpoint ``exp(Md)``; an upper point beyond that
        endpoint therefore returns the same full-support atom.  Negative
        finite upper points are supported by the oriented analytic plateau
        integral.
        """

        with mp.workdps(self._dps()):
            power_mp = _as_mpf(power)
            cutoff_power_mp = _as_mpf(cutoff_power)
            if not mp.isfinite(power_mp):
                raise ValueError("power must be finite")
            if not mp.isfinite(cutoff_power_mp) or cutoff_power_mp < 0:
                raise ValueError("cutoff_power must be finite and nonnegative")
            chosen_order = self.quadrature_order if order is None else int(order)
            if chosen_order < 4:
                raise ValueError("order must be at least 4")
            endpoint = self._cutoff() if upper is None else _as_mpf(upper)
            if not mp.isfinite(endpoint):
                raise ValueError("upper must be finite or None")
            if endpoint <= 1:
                return _exp_integral(power_mp, mp.mpf(0), endpoint)

            plateau = _exp_integral(power_mp, mp.mpf(0), mp.mpf(1))
            transition = self._transition_integral(
                power_mp,
                endpoint,
                cutoff_power_mp,
                chosen_order,
            )
            return plateau + transition

    def _base_mass(self, y: Any, *, order: int | None = None) -> mp.mpf:
        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            if y_mp < 0:
                return mp.exp(y_mp)
            return 1 + self.weighted_integral(power=1, upper=y_mp, cutoff_power=1, order=order)

    def I_z(self, y: Any, Z: Any, *, order: int | None = None) -> mp.mpf:
        """Return the continuous dimensionless axial mass primitive."""

        with mp.workdps(self._dps()):
            return 4 * _as_mpf(Z) * self._base_mass(y, order=order)

    def I_uz2(self, y: Any, Z: Any, *, order: int | None = None) -> mp.mpf:
        """Return cumulative ``integral exp(s) Uz(s,Z)^2 ds`` through ``y``."""

        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            z_mp = _as_mpf(Z)
            if y_mp < 0:
                base = mp.exp(y_mp)
            else:
                base = 1 + self.weighted_integral(
                    power=1,
                    upper=y_mp,
                    cutoff_power=2,
                    order=order,
                )
            return 16 * z_mp**2 * base

    def mean(self, y: Any, Z: Any, *, order: int | None = None) -> mp.mpf:
        """Return ``mean(y,Z) = exp(-y) I_z(y,Z)``."""

        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            z_mp = _as_mpf(Z)
            if y_mp < 0:
                return 4 * z_mp
            return mp.exp(-y_mp) * self.I_z(y_mp, z_mp, order=order)

    def mean_y(self, y: Any, Z: Any, *, order: int | None = None) -> mp.mpf:
        """Return the analytic identity ``mean_y = Uz - mean``."""

        with mp.workdps(self._dps()):
            return self.Uz(y, Z) - self.mean(y, Z, order=order)

    def mean_Z(self, y: Any, Z: Any | None = None, *, order: int | None = None) -> mp.mpf:
        """Return the analytic ``Z`` derivative of the mean.

        The optional ``Z`` argument mirrors :meth:`Uz_Z`; the result is
        independent of its value because the provider is linear in ``Z``.
        """

        del Z
        with mp.workdps(self._dps()):
            y_mp = _as_mpf(y)
            if y_mp < 0:
                return mp.mpf(4)
            return 4 * mp.exp(-y_mp) * self._base_mass(y_mp, order=order)

    def mean_jet(
        self,
        y: Any,
        Z: Any,
        *,
        order: int | None = None,
    ) -> dict[str, mp.mpf]:
        """Return ``mean`` and its analytic ``Z`` jet as a mapping."""

        with mp.workdps(self._dps()):
            return {
                "mean": self.mean(y, Z, order=order),
                "mean_y": self.mean_y(y, Z, order=order),
                "mean_Z": self.mean_Z(y, Z, order=order),
            }

    def cutoff_jet(self, y: Any) -> dict[str, mp.mpf]:
        """Return the cutoff value and analytic derivative as a mapping."""

        with mp.workdps(self._dps()):
            return {"value": self.c(y), "derivative": self.c_y(y)}

    def full_incoming_rows(
        self,
        Z: Any,
        *,
        order: int | None = None,
    ) -> dict[str, mp.mpf]:
        """Return full-support incoming rows ``I_z`` and ``I_uz2``."""

        with mp.workdps(self._dps()):
            terminal = self._cutoff()
            return {
                "I_z": self.I_z(terminal, Z, order=order),
                "I_uz2": self.I_uz2(terminal, Z, order=order),
            }

    incoming_rows = full_incoming_rows
    full_rows = full_incoming_rows

    def moments(self, Z: Any, *, order: int | None = None) -> dict[str, mp.mpf]:
        """Compatibility alias for installers requesting full incoming rows."""

        return self.full_incoming_rows(Z, order=order)


def _n(value: Any, digits: int = 60) -> str:
    return mp.nstr(_as_mpf(value), digits)


def _relative_difference(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    scale = abs(right)
    return abs(left - right) / scale if scale else abs(left - right)


def run() -> dict[str, Any]:
    """Run inexpensive point, endpoint, and quadrature refinement receipts."""

    precision = 120
    coarse_order = 20
    fine_order = 36
    with mp.workdps(precision + 20):
        provider = ContinuousIncomingAxial(
            Md=".5",
            precision=precision,
            quadrature_order=coarse_order,
        )
        fine = ContinuousIncomingAxial(
            Md=".5",
            precision=precision,
            quadrature_order=fine_order,
        )
        z = mp.mpf(".375")
        cutoff = provider.cutoff_y

        # A five-point differential replay compares the cumulative primitive
        # derivative with the direct point source term exp(y)*Uz.
        differential = []
        for y in (mp.mpf(".5"), mp.mpf("1.2"), (1 + cutoff) / 2):
            h = mp.mpf("1e-5")
            values = {
                offset: provider.I_z(y + offset * h, z)
                for offset in (-2, -1, 1, 2)
            }
            derivative = (
                values[-2]
                - 8 * values[-1]
                + 8 * values[1]
                - values[2]
            ) / (12 * h)
            point = mp.exp(y) * provider.Uz(y, z)
            differential.append(
                {
                    "y": _n(y),
                    "step": _n(h),
                    "primitive_derivative": _n(derivative),
                    "point_weighted_Uz": _n(point),
                    "relative_error": _n(_relative_difference(derivative, point)),
                }
            )

        # Flat joins are checked from both sides.  The source endpoint is
        # sampled at a moderate h so the receipt is not just a precision-zero
        # statement about an exponentially tiny correction.
        continuity = []
        for name, edge, h in (
            ("plateau_join", mp.mpf(1), mp.mpf("1e-5")),
            ("cutoff_join", cutoff, mp.mpf(".05")),
        ):
            left = edge - h
            right = edge + h
            continuity.append(
                {
                    "name": name,
                    "edge": _n(edge),
                    "step": _n(h),
                    "c_jump": _n(provider.c(right) - provider.c(left)),
                    "c_y_jump": _n(provider.c_y(right) - provider.c_y(left)),
                    "Uz_jump": _n(provider.Uz(right, z) - provider.Uz(left, z)),
                    "I_z_jump": _n(provider.I_z(right, z) - provider.I_z(left, z)),
                    "mean_jump": _n(provider.mean(right, z) - provider.mean(left, z)),
                }
            )

        order_refinement = []
        for power, cutoff_power in ((1, 1), (1, 2), (mp.mpf(".5"), 1)):
            coarse = provider.weighted_integral(
                power=power,
                cutoff_power=cutoff_power,
            )
            refined = fine.weighted_integral(
                power=power,
                cutoff_power=cutoff_power,
            )
            order_refinement.append(
                {
                    "power": _n(power),
                    "cutoff_power": _n(cutoff_power),
                    "coarse_order": coarse_order,
                    "fine_order": fine_order,
                    "coarse": _n(coarse),
                    "fine": _n(refined),
                    "relative_change": _n(_relative_difference(coarse, refined)),
                }
            )

        partial_upper = [
            mp.mpf(1),
            mp.mpf("1.2"),
            (1 + cutoff) / 2,
            cutoff,
            cutoff + mp.mpf(".2"),
        ]
        partial_transition = [
            {
                "upper": _n(upper),
                "weighted_integral": _n(
                    provider.weighted_integral(upper=upper, cutoff_power=1)
                ),
                "c_at_upper": _n(provider.c(upper)),
                "I_z": _n(provider.I_z(upper, z)),
                "I_uz2": _n(provider.I_uz2(upper, z)),
            }
            for upper in partial_upper
        ]

        rows = provider.full_incoming_rows(z)
        terminal_rows = {
            name: _n(value)
            for name, value in rows.items()
        }
        report = {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": SOURCE_SECTION,
            "api": {
                "class": "ContinuousIncomingAxial",
                "point": ["c", "c_y", "Uz", "Uz_y", "Uz_Z"],
                "cumulative": ["I_z", "I_uz2", "mean", "mean_y", "mean_Z"],
                "weighted_integral": "weighted_integral(power=1, upper=None, cutoff_power=1, order=None)",
                "full_rows": "full_incoming_rows(Z, order=None) -> I_z, I_uz2",
            },
            "parameters": {
                "Md": _n(provider.Md),
                "cutoff_y": _n(cutoff),
                "precision": precision,
                "coarse_order": coarse_order,
                "fine_order": fine_order,
                "Z": _n(z),
            },
            "differential_primitive_vs_point": differential,
            "continuity": continuity,
            "order_refinement": order_refinement,
            "partial_transition": partial_transition,
            "full_incoming_rows": terminal_rows,
            "definitions": {
                "Uz": "4*Z*c(y)",
                "I_z": "4*Z*exp(y) for y<0; 4*Z*(1+integral_0^y exp(s)c(s) ds) for y>=0",
                "I_uz2": "16*Z^2*exp(y) for y<0; 16*Z^2*(1+integral_0^y exp(s)c(s)^2 ds) for y>=0",
                "mean": "exp(-y)*I_z",
                "mean_y": "Uz-mean",
                "mean_Z": "analytic Z derivative of mean",
            },
            "claims": {
                "point_and_primitive_share_sigma": True,
                "arbitrary_precision_Z_path": True,
                "supports_partial_transition_endpoints": True,
                "quadrature_enclosure_certified": False,
                "global_energy_certified": False,
                "recursion_certified": False,
            },
            "scope": (
                "Continuous incoming axial source point/cumulative provider and "
                "focused numerical diagnostics only; no quadrature enclosure, "
                "global finite-energy certification, recursive closure, or "
                "installation into the global profile is claimed."
            ),
        }

    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "differential_errors": [
                    row["relative_error"] for row in differential
                ],
                "order_changes": [
                    row["relative_change"] for row in order_refinement
                ],
                "full_rows": terminal_rows,
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
