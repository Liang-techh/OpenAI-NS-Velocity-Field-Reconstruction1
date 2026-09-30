"""Arbitrary precision axis data and exact first radial jets for Part I.

This module is the small regular-core boundary adapter for Lei--Ren Part I,
Section 8.  It implements (8.7), rather than fitting a radial profile.  The
axis pressure ``P_0`` is deliberately an injected callback: this file does
not manufacture a pressure datum and makes no claim about the completed core.

The real-axis datum is

    U0z = 4 Z + j,
    H0 = (1-delta) Z / 2 + (1-Z^2) U0z,
    G(Z) = integral_{Z0}^Z L H0 / (H0^2 + sigma0^2) dZ,
    F0 = exp(-log(Cstar) - Lambda G),

where ``Z0`` is the unique real zero of ``H0``.  The first radial slopes are
the two exact equations (8.7).  All arithmetic for the construction is done
with ``mpmath``; the pressure derivative is a fourth-order centered finite
difference of the same pressure callback and is returned with refinement data.

The optional first-order radial values and cumulative moments are local Taylor
replays only.  They are useful to attach a higher-order core recurrence, but
do not certify the nonlinear core, cone, pressure closure, or PDE recursion.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
from typing import Any, Callable

import mpmath as mp


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"


def _mp(value: Any) -> mp.mpf:
    """Convert through text so Decimal and string inputs keep their digits."""

    if isinstance(value, bool):
        raise TypeError("boolean is not a real parameter")
    if isinstance(value, mp.mpf):
        return mp.mpf(value)
    return mp.mpf(str(value))


def _finite_mp(value: Any, name: str) -> mp.mpf:
    out = _mp(value)
    if not mp.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _nstr(value: Any, digits: int = 40) -> str:
    if value is None:
        return "null"
    return mp.nstr(_mp(value), digits)


def _relative_error(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    scale = max(abs(reference), mp.mpf("1e-120"))
    return abs(value - reference) / scale


def _pressure_scalar(result: Any) -> mp.mpf:
    """Extract P_0 from a scalar, pair, or source-profile mapping."""

    if isinstance(result, Mapping):
        for key in ("P0", "P", "pressure", "axis_pressure", "value"):
            if key in result:
                return _finite_mp(result[key], "pressure callback value")
        raise KeyError("pressure mapping must contain P0, P, pressure, or value")
    if isinstance(result, (tuple, list)):
        if not result:
            raise ValueError("pressure callback returned an empty sequence")
        return _finite_mp(result[0], "pressure callback value")
    return _finite_mp(result, "pressure callback value")


def _call_pressure(provider: Any, z: mp.mpf) -> mp.mpf:
    """Call a shared pressure source without replacing it by a local model."""

    callback = provider
    if not callable(callback):
        for name in ("axis_pressure", "pressure_at_axis", "pressure"):
            candidate = getattr(provider, name, None)
            if callable(candidate):
                callback = candidate
                break
    if not callable(callback):
        raise TypeError("pressure must be callable or expose axis_pressure(Z)")
    return _pressure_scalar(callback(z))


@dataclass(frozen=True)
class AxisJet:
    """Raw axis values and first radial jets at one Z."""

    Z: mp.mpf
    d: mp.mpf
    L: mp.mpf
    H0: mp.mpf
    G: mp.mpf
    U0z: mp.mpf
    U0z_Z: mp.mpf
    F0: mp.mpf
    F0_Z: mp.mpf
    F_R: mp.mpf
    Uz_R: mp.mpf
    P0: mp.mpf
    P0_Z: mp.mpf

    def as_dict(self) -> dict[str, mp.mpf]:
        return {
            "Z": self.Z,
            "d": self.d,
            "L": self.L,
            "H0": self.H0,
            "G": self.G,
            "U0z": self.U0z,
            "U0z_Z": self.U0z_Z,
            "F0": self.F0,
            "F0_Z": self.F0_Z,
            "F_R": self.F_R,
            "Uz_R": self.Uz_R,
            "P0": self.P0,
            "P0_Z": self.P0_Z,
        }


class RegularCoreAxisJets:
    """Section 8 axis data with the exact first radial slopes.

    Parameters are accepted as strings or ``Decimal`` values.  ``logC`` is
    ``log(C_*)``; if omitted, the minimal algebraic guard
    ``log(C_*) = 2 log(Lambda)`` is used.  The theorem's additional
    ``Lambda * A_Omega`` margin is a separate complex-domain requirement and
    is recorded as uncertified here.

    ``pressure`` is an optional shared callback.  It must evaluate the actual
    axis datum ``P_0(Z)`` used by the source profile.  Calling ``axis_jets``
    without it raises instead of silently inventing pressure.
    """

    def __init__(
        self,
        *,
        j: Any = "0.02",
        Lambda: Any = "2500",
        logC: Any | None = None,
        C: Any | None = None,
        delta: Any = "1e-32",
        pressure: Any | None = None,
        precision: int = 160,
        pressure_step: Any = "1e-5",
        quadrature_dps: int | None = None,
    ) -> None:
        self.precision = max(80, int(precision))
        with mp.workdps(self.precision):
            self.j = _finite_mp(j, "j")
            self.Lambda = _finite_mp(Lambda, "Lambda")
            self.delta = _finite_mp(delta, "delta")
            if not (self.j > 0 and self.j <= mp.mpf(1) / 20):
                raise ValueError("require 0 < j <= 1/20")
            if self.Lambda <= 0:
                raise ValueError("Lambda must be positive")
            if self.delta < 0 or self.delta >= 1:
                raise ValueError("require 0 <= delta < 1")
            if logC is not None and C is not None:
                raise ValueError("pass either logC or C, not both")
            if C is not None:
                c_value = _finite_mp(C, "C")
                if c_value <= 0:
                    raise ValueError("C must be positive")
                self.logC = mp.log(c_value)
            elif logC is None:
                self.logC = 2 * mp.log(self.Lambda)
            else:
                self.logC = _finite_mp(logC, "logC")
            if self.logC < 2 * mp.log(self.Lambda):
                raise ValueError("require logC >= 2 log(Lambda)")
            self.pressure_step = _finite_mp(pressure_step, "pressure_step")
            if self.pressure_step <= 0:
                raise ValueError("pressure_step must be positive")
        self.pressure = pressure
        self.quadrature_dps = max(
            self.precision,
            int(quadrature_dps) if quadrature_dps is not None else self.precision,
        )
        with mp.workdps(self.precision):
            self.sigma0 = self.j / 500
            self.Z0 = self._find_unique_root()
            self._root_width = self.sigma0 / abs(self.H0_Z(self.Z0))
        self._g_cache: dict[str, mp.mpf] = {}

    # ------------------------------------------------------------------
    # Axis coefficients and the narrow root layer.
    # ------------------------------------------------------------------
    def _validate_z(self, Z: Any) -> mp.mpf:
        z = _finite_mp(Z, "Z")
        if abs(z) > 1:
            raise ValueError("Z must satisfy |Z| <= 1")
        return z

    def d(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        return 1 - z * z

    def L(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        return 1 - self.delta * z * z

    def U0z(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        return 4 * z + self.j

    def U0z_Z(self, Z: Any) -> mp.mpf:
        self._validate_z(Z)
        return mp.mpf(4)

    def H0(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        return (1 - self.delta) * z / 2 + (1 - z * z) * (4 * z + self.j)

    def H0_Z(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        return (1 - self.delta) / 2 + 4 * (1 - z * z) - 2 * z * (4 * z + self.j)

    def _find_unique_root(self) -> mp.mpf:
        """Bracket the source's unique real zero, retaining MP precision."""

        left = mp.mpf(-1)
        right = mp.mpf(1)
        f_left = self.H0(left)
        f_right = self.H0(right)
        if not (f_left < 0 < f_right):
            raise ValueError("H0 does not have the required bracketed root")
        # Bisection avoids the precision loss that a float Newton seed would
        # introduce.  The extra iterations leave a comfortable guard below
        # the requested mpmath precision.
        iterations = max(320, 4 * self.precision)
        for _ in range(iterations):
            mid = (left + right) / 2
            f_mid = self.H0(mid)
            if f_mid == 0:
                return mid
            if f_mid < 0:
                left = mid
            else:
                right = mid
        return (left + right) / 2

    def G_integrand(self, Z: Any) -> mp.mpf:
        z = self._validate_z(Z)
        h = self.H0(z)
        return self.L(z) * h / (h * h + self.sigma0 * self.sigma0)

    def _integrate_segment(self, left: mp.mpf, right: mp.mpf) -> mp.mpf:
        if left == right:
            return mp.mpf(0)
        with mp.workdps(self.quadrature_dps):
            return mp.quad(self.G_integrand, [left, right])

    def G(self, Z: Any) -> mp.mpf:
        """Evaluate G with explicit subdivision at Z0 +/- sigma/H0_Z."""

        z = self._validate_z(Z)
        key = mp.nstr(z, self.precision)
        if key in self._g_cache:
            return self._g_cache[key]
        with mp.workdps(self.quadrature_dps):
            if z == self.Z0:
                value = mp.mpf(0)
            else:
                lo = min(self.Z0, z)
                hi = max(self.Z0, z)
                width = self._root_width
                split = [lo, hi]
                for point in (self.Z0 - width, self.Z0, self.Z0 + width):
                    if lo < point < hi:
                        split.append(point)
                split = sorted(set(split))
                value = mp.mpf(0)
                for a, b in zip(split, split[1:]):
                    value += self._integrate_segment(a, b)
                if z < self.Z0:
                    value = -value
            value = +value
        self._g_cache[key] = value
        return value

    def axis_data(self, Z: Any) -> dict[str, mp.mpf]:
        """Return pressure-independent axis values and angular jet data."""

        with mp.workdps(self.precision):
            z = self._validate_z(Z)
            d = 1 - z * z
            L = 1 - self.delta * z * z
            h0 = (1 - self.delta) * z / 2 + d * (4 * z + self.j)
            g = self.G(z)
            f0 = mp.exp(-self.logC - self.Lambda * g)
            f0_z = -self.Lambda * self.G_integrand(z) * f0
            return {
                "Z": z,
                "d": d,
                "L": L,
                "H0": h0,
                "G": g,
                "U0z": 4 * z + self.j,
                "U0z_Z": mp.mpf(4),
                "F0": f0,
                "F0_Z": f0_z,
            }

    # ------------------------------------------------------------------
    # Same-pressure derivative and exact Eq. (8.7) jets.
    # ------------------------------------------------------------------
    def pressure_value(self, Z: Any, *, pressure: Any | None = None) -> mp.mpf:
        provider = self.pressure if pressure is None else pressure
        if provider is None:
            raise ValueError("P0 callback is required for the axial jet")
        with mp.workdps(self.precision):
            return _call_pressure(provider, self._validate_z(Z))

    def pressure_derivative(
        self,
        Z: Any,
        *,
        pressure: Any | None = None,
        step: Any | None = None,
    ) -> mp.mpf:
        """Fourth-order centered P0_Z from the same callback as P0."""

        with mp.workdps(self.precision):
            z = self._validate_z(Z)
            h = self.pressure_step if step is None else _finite_mp(step, "step")
            if h <= 0 or abs(z) + 2 * h > 1:
                raise ValueError("pressure stencil must remain inside |Z| <= 1")
            p = lambda offset: self.pressure_value(z + offset * h, pressure=pressure)
            return (p(-2) - 8 * p(-1) + 8 * p(1) - p(2)) / (12 * h)

    def pressure_derivative_refinement(
        self,
        Z: Any,
        *,
        pressure: Any | None = None,
        steps: tuple[Any, Any] = ("1e-5", "5e-6"),
    ) -> dict[str, Any]:
        with mp.workdps(self.precision):
            values = [self.pressure_derivative(Z, pressure=pressure, step=s) for s in steps]
            return {
                "steps": [str(s) for s in steps],
                "values": [_nstr(v, self.precision) for v in values],
                "relative_difference": _nstr(_relative_error(values[0], values[1]), 30),
                "same_pressure_callback": True,
            }

    def axis_jets(
        self,
        Z: Any,
        *,
        pressure: Any | None = None,
        pressure_step: Any | None = None,
    ) -> AxisJet:
        """Return the exact Eq. (8.7) first radial slopes at Z."""

        with mp.workdps(self.precision):
            data = self.axis_data(Z)
            z = data["Z"]
            p0 = self.pressure_value(z, pressure=pressure)
            p0_z = self.pressure_derivative(z, pressure=pressure, step=pressure_step)
            d = data["d"]
            L = data["L"]
            h0 = data["H0"]
            u0 = data["U0z"]
            u0_z = data["U0z_Z"]
            f0 = data["F0"]
            f0_z = data["F0_Z"]
            f_r = (
                (1 + self.delta / 2 - z * u0 - d * u0_z) * f0 + h0 * f0_z
            ) / (4 * L)
            uz_r = (
                (1 + self.delta) / 2 * (1 - 2 * z * u0) * u0
                + h0 * u0_z
                + d * p0_z
                - 2 * (1 + self.delta) * z * p0
            ) / (2 * L)
            return AxisJet(
                Z=z,
                d=d,
                L=L,
                H0=h0,
                G=data["G"],
                U0z=u0,
                U0z_Z=u0_z,
                F0=f0,
                F0_Z=f0_z,
                F_R=f_r,
                Uz_R=uz_r,
                P0=p0,
                P0_Z=p0_z,
            )

    # ------------------------------------------------------------------
    # Local first-order replay, for parent recurrence wiring only.
    # ------------------------------------------------------------------
    def first_order_values(
        self,
        R: Any,
        Z: Any,
        *,
        pressure: Any | None = None,
        pressure_step: Any | None = None,
    ) -> dict[str, mp.mpf]:
        with mp.workdps(self.precision):
            r = _finite_mp(R, "R")
            if r < 0:
                raise ValueError("R must be nonnegative")
            jet = self.axis_jets(Z, pressure=pressure, pressure_step=pressure_step)
            F = jet.F0 + r * jet.F_R
            uz = jet.U0z + r * jet.Uz_R
            p = jet.P0 + r * jet.F0 * jet.F0
            return {
                "R": r,
                "Z": jet.Z,
                "F": F,
                "Uz": uz,
                "Utheta": mp.sqrt(2 * r) * F if r else mp.mpf(0),
                "P": p,
                "F_R": jet.F_R,
                "Uz_R": jet.Uz_R,
            }

    def first_order_moments(
        self,
        R: Any,
        Z: Any,
        *,
        pressure: Any | None = None,
        pressure_step: Any | None = None,
    ) -> dict[str, mp.mpf]:
        """Cumulative source moments of the first-order axis Taylor replay."""

        with mp.workdps(self.precision):
            r = _finite_mp(R, "R")
            if r < 0:
                raise ValueError("R must be nonnegative")
            jet = self.axis_jets(Z, pressure=pressure, pressure_step=pressure_step)
            f0, fr = jet.F0, jet.F_R
            u0, ur = jet.U0z, jet.Uz_R
            return {
                "theta": f0 * r * r + mp.mpf(2) / 3 * fr * r**3,
                "z": u0 * r + ur * r * r / 2,
                "theta_z": f0 * u0 * r * r
                + mp.mpf(2) / 3 * (fr * u0 + f0 * ur) * r**3,
                "z_theta": u0 * u0 * r
                + (u0 * ur - f0 * f0 / 2) * r * r,
                "p": f0 * f0 * r + f0 * fr * r * r,
            }

    def axis_balance_residual(
        self,
        Z: Any,
        *,
        pressure: Any | None = None,
        pressure_step: Any | None = None,
    ) -> dict[str, mp.mpf]:
        """Replay Eq. (8.2) at R = 0 using the exact axis jets."""

        with mp.workdps(self.precision):
            jet = self.axis_jets(Z, pressure=pressure, pressure_step=pressure_step)
            z, d, L, h0 = jet.Z, jet.d, jet.L, jet.H0
            u0, u0_z = jet.U0z, jet.U0z_Z
            # Mz/R -> U0z and (Mz)_Z/R -> U0z_Z at the axis.
            w0 = 1 - (1 - self.delta) * z * u0 - d * u0_z
            angular_rhs = (
                w0 * jet.F0
                + self.delta / 2 * (1 - 2 * z * u0) * jet.F0
                + h0 * jet.F0_Z
            )
            axial_rhs = (
                (1 + self.delta) / 2 * (1 - 2 * z * u0) * u0
                + h0 * u0_z
                + d * jet.P0_Z
                - 2 * (1 + self.delta) * z * jet.P0
            )
            return {
                "angular": 4 * L * jet.F_R - angular_rhs,
                "axial": 2 * L * jet.Uz_R - axial_rhs,
                "angular_relative": _relative_error(4 * L * jet.F_R, angular_rhs),
                "axial_relative": _relative_error(2 * L * jet.Uz_R, axial_rhs),
            }

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "equations": ["(8.1)", "(8.7)", "(8.9)", "(8.44)"],
            "j": _nstr(self.j, self.precision),
            "Lambda": _nstr(self.Lambda, self.precision),
            "delta": _nstr(self.delta, self.precision),
            "logC": _nstr(self.logC, self.precision),
            "sigma0": _nstr(self.sigma0, self.precision),
            "Z0": _nstr(self.Z0, self.precision),
            "root_width_sigma_over_abs_H0_Z": _nstr(self._root_width, self.precision),
            "precision": self.precision,
            "pressure_callback_supplied": self.pressure is not None,
            "pressure_derivative": "fourth-order centered FD of the same P0 callback",
            "theorem_complex_A_Omega_margin_supplied": False,
            "global_core_certified": False,
            "cone_certified": False,
            "pressure_closure_certified": False,
        }


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value, 50)
    if isinstance(value, Mapping):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(v) for v in value]
    return value


def run_focused_checks() -> dict[str, Any]:
    """Run algebraic checks with a clearly labelled zero-pressure fixture."""

    # This callback is only an algebraic fixture.  It is intentionally not
    # presented as the source pressure; parent integration must inject the
    # shared completed-pressure callback.
    def zero_pressure(_Z: Any) -> mp.mpf:
        return mp.mpf(0)

    jets = RegularCoreAxisJets(
        j="0.02",
        Lambda="2500",
        # Let the constructor evaluate the equality at its own precision;
        # computing this mpmath value outside that context would round it.
        logC=None,
        delta="1e-32",
        pressure=zero_pressure,
        precision=100,
        pressure_step="1e-6",
    )
    points: list[dict[str, Any]] = []
    with mp.workdps(jets.precision):
        for z in ("-0.75", "-0.25", "0", "0.25", "0.75"):
            axis = jets.axis_jets(z)
            h = mp.mpf("1e-6")
            if abs(_mp(z)) + 2 * h < 1:
                g_fd = (jets.G(_mp(z) + h) - jets.G(_mp(z) - h)) / (2 * h)
                f_fd = (jets.axis_data(_mp(z) + h)["F0"] - jets.axis_data(_mp(z) - h)["F0"]) / (2 * h)
                g_err = _relative_error(g_fd, jets.G_integrand(z))
                f_err = _relative_error(f_fd, axis.F0_Z)
            else:
                g_err = mp.nan
                f_err = mp.nan
            residual = jets.axis_balance_residual(z)
            points.append(
                {
                    "Z": z,
                    "H0": _nstr(axis.H0, 40),
                    "F0": _nstr(axis.F0, 40),
                    "F0_Z": _nstr(axis.F0_Z, 40),
                    "F_R": _nstr(axis.F_R, 40),
                    "Uz_R": _nstr(axis.Uz_R, 40),
                    "G_derivative_relative_error": _nstr(g_err, 30),
                    "F0_derivative_relative_error": _nstr(f_err, 30),
                    "axis_balance_relative": {
                        "angular": _nstr(residual["angular_relative"], 30),
                        "axial": _nstr(residual["axial_relative"], 30),
                    },
                }
            )
        max_g = max(mp.mpf(row["G_derivative_relative_error"]) for row in points)
        max_f = max(mp.mpf(row["F0_derivative_relative_error"]) for row in points)
        max_bal = max(
            mp.mpf(row["axis_balance_relative"][name])
            for row in points
            for name in ("angular", "axial")
        )
        report = {
            "metadata": jets.metadata(),
            "pressure_fixture": "P0(Z)=0 only for local Eq. (8.7) replay; not source pressure",
            "root_residual": _nstr(abs(jets.H0(jets.Z0)), 30),
            "points": points,
            "max_G_derivative_relative_error": _nstr(max_g, 30),
            "max_F0_derivative_relative_error": _nstr(max_f, 30),
            "max_axis_balance_relative_error": _nstr(max_bal, 30),
            "axis_jet_checks_passed": bool(max_bal < mp.mpf("1e-80")),
            "source_pressure_attached": False,
            "global_core_certified": False,
            "cone_certified": False,
            "pressure_closure_certified": False,
        }
    return report


if __name__ == "__main__":
    receipt = run_focused_checks()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(_json_value(receipt), indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "axis_jet_checks_passed": receipt["axis_jet_checks_passed"],
                "root_residual": receipt["root_residual"],
                "max_axis_balance_relative_error": receipt["max_axis_balance_relative_error"],
            }
        ),
        flush=True,
    )


__all__ = ["AxisJet", "RegularCoreAxisJets", "run_focused_checks"]
