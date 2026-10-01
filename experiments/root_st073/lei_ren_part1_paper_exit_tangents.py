"""Z tangents for the finite Section 9.25 exit bridge.

The bridge state is ``g = log(F/F_a)``, ``u = U^z`` and six normalized
cumulative moments.  Its tangent is integrated together with the state:

    g_y = A = -chi D/2,
    u_y = exp(g) B,
    g_{yZ} = A_Z,
    u_{yZ} = exp(g) (B_Z + B g_Z).

By default the comparison driver's ``A_Z`` and ``B_Z`` use an mpmath
fourth-order centered stencil on already normalized drivers.  An optional
``analytic_driver_provider`` can instead supply the four scalar fields
``A``, ``B``, ``A_Z`` and ``B_Z`` at the center ``(y, Z)``; that path performs
no off-center Z evaluations.  Neither path certifies the source bridge or
the relaxed cone.  The initial tangent comes from the actual core
coefficient and moment jets, including the derivative of the ``F_a``
normalization.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from typing import Any

import mpmath as mp


from lei_ren_part1_paper_exit_bridge import ExitBridge  # noqa: E402
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress  # noqa: E402


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _mp(value: Any) -> mp.mpf:
    """Coerce without routing an existing mpmath value through binary64."""

    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


class ExitTangents(ExitBridge):
    """Finite exit bridge carrying analytic core and numerical-driver Z jets."""

    def __init__(
        self,
        comparison: Any,
        *,
        epsilon: Any | None = None,
        steps: int = 16,
        derivative_step: Any = "1e-45",
        analytic_driver_provider: Any | None = None,
    ) -> None:
        super().__init__(comparison, epsilon=epsilon, steps=steps)
        with mp.workdps(self.precision):
            self.derivative_step = _mp(derivative_step)
            if self.derivative_step <= 0:
                raise ValueError("derivative_step must be positive")
            if 2 * self.derivative_step >= 1:
                raise ValueError("derivative_step must keep the centered stencil in |Z|<1")
        self.analytic_driver_provider = analytic_driver_provider
        self.driver_derivative_method = (
            "analytic driver provider"
            if analytic_driver_provider is not None
            else "mp fourth-order centered stencil"
        )
        self._driver_cache: dict[tuple[str, str], dict[str, Any]] = {}

    def _z_key(self, Z: Any) -> str:
        return mp.nstr(_mp(Z), self.precision)

    @lru_cache(maxsize=16)
    def initial(self, Z: Any) -> dict[str, Any]:
        """Return actual core bridge state and its Z tangent at ``y=0``."""

        with mp.workdps(self.precision):
            z = _mp(Z)
            if abs(z) >= 1:
                raise ValueError("ExitTangents requires |Z| < 1")
            base = dict(super().initial(z))
            base["state"] = list(base["state"])
            core = base["core"]
            Fa = _mp(base["Fa"])
            Fa_Z = _mp(core["F_Z"])
            r = _mp(self.r)
            r2 = r * r
            m = core["moments"]
            mz = core["moments_Z"]

            # The existing base state already contains axial and swirl
            # integrals.  Reconstruct their Z jets from the coefficient rows
            # so the normalization derivative is explicit and reproducible.
            coefficients = self.core.coefficients(z)
            f_rows = coefficients["F"]
            u_rows = coefficients["Uz"]
            axial = mp.mpf(0)
            axial_Z = mp.mpf(0)
            swirl = mp.mpf(0)
            swirl_Z = mp.mpf(0)
            for i, ui in enumerate(u_rows):
                for j, uj in enumerate(u_rows):
                    power = i + j + 1
                    factor = r ** power / power
                    axial += ui[0] * uj[0] * factor
                    axial_Z += (ui[1] * uj[0] + ui[0] * uj[1]) * factor
            for i, fi in enumerate(f_rows):
                for j, fj in enumerate(f_rows):
                    power = i + j + 2
                    factor = r ** power / power
                    swirl += fi[0] * fj[0] * factor
                    swirl_Z += (fi[1] * fj[0] + fi[0] * fj[1]) * factor

            Fa2r2 = Fa * Fa * r2
            base_p = _mp(m["p"])
            base_p_Z = _mp(mz["p"])
            theta_scale = Fa * r2
            theta = _mp(m["theta"])
            theta_Z = _mp(mz["theta"])
            theta_z = _mp(m["theta_z"])
            theta_z_Z = _mp(mz["theta_z"])
            m_z = _mp(m["z"])
            m_z_Z = _mp(mz["z"])
            state_Z = [
                mp.mpf(0),
                _mp(core["Uz_Z"]),
                (theta_Z - theta * Fa_Z / Fa) / theta_scale,
                m_z_Z / r,
                (theta_z_Z - theta_z * Fa_Z / Fa) / theta_scale,
                axial_Z / r,
                (swirl_Z - 2 * Fa_Z / Fa * swirl) / Fa2r2,
                (base_p_Z - 2 * Fa_Z / Fa * base_p) / (Fa * Fa * r),
            ]
            base["state_Z"] = state_Z
            base["Fa_Z"] = Fa_Z
            base["P0_Z"] = _mp(core["P_Z"]) - base_p_Z
            base["axial"] = axial
            base["axial_Z"] = axial_Z
            base["swirl"] = swirl
            base["swirl_Z"] = swirl_Z
            base["initial_tangent_source"] = (
                "core coefficient and moments_Z rows; Fa normalization differentiated analytically"
            )
            return base

    def _normalized_driver(self, y: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        """Evaluate normalized A and B before taking a Z difference."""

        with mp.workdps(self.precision):
            initial = self.initial(Z)
            comparison = self.comparison.evaluate(y, Z)
            R = _mp(self.r) * mp.exp(y)
            chi = _mp(self.multiplier(y))
            Fa = _mp(initial["Fa"])
            barF = _mp(comparison["F"])
            if barF == 0 or Fa == 0:
                raise ArithmeticError("normalized exit driver encountered zero F")
            Iz = _mp(comparison.get("I_z", comparison["E"] * comparison["F"]))
            A = -chi * _mp(comparison["D"]) / 2
            B = -chi * mp.sqrt(R / 2) * (Fa / barF) * Iz
            return {
                "A": A,
                "B": B,
                "R": R,
                "chi": chi,
                "Fa": Fa,
                "barF": barF,
                "I_z": Iz,
                "comparison": comparison,
                "normalization": "D and B differentiated after Fa/barF normalization",
            }

    def _analytic_driver(self, y: mp.mpf, Z: mp.mpf) -> dict[str, Any] | None:
        """Resolve the optional center-only analytic driver hook.

        The provider may be a callable ``(y, Z) -> Mapping`` or expose an
        ``exit_driver(y, Z)`` / ``analytic_exit_driver(y, Z)`` method.  The
        returned mapping must contain scalar ``A``, ``B``, ``A_Z`` and
        ``B_Z`` values.  No stencil points are constructed on this path.
        """

        provider = self.analytic_driver_provider
        if provider is None:
            return None
        method = getattr(provider, "exit_driver", None)
        if method is None:
            method = getattr(provider, "analytic_exit_driver", None)
        if method is None and callable(provider):
            method = provider
        if method is None or not callable(method):
            raise TypeError(
                "analytic_driver_provider must be callable or expose "
                "exit_driver(y, Z)"
            )
        raw = method(y, Z)
        if not isinstance(raw, Mapping):
            raise TypeError("analytic exit driver must return a mapping")
        required = ("A", "B", "A_Z", "B_Z")
        missing = [name for name in required if name not in raw]
        if missing:
            raise KeyError(
                "analytic exit driver is missing " + ", ".join(missing)
            )
        result = dict(raw)
        for name in required:
            result[name] = _mp(result[name])
        result["derivative_step"] = None
        result["derivative_method"] = "analytic driver provider"
        result["analytic_driver_used"] = True
        result["normalized_driver"] = True
        return result

    def _driver_tangent(self, y: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        """Return normalized driver values and their Z tangent."""

        with mp.workdps(self.precision):
            z = _mp(Z)
            key = (mp.nstr(y, self.precision), mp.nstr(z, self.precision))
            cached = self._driver_cache.get(key)
            if cached is not None:
                return cached
            analytic = self._analytic_driver(_mp(y), z)
            if analytic is not None:
                self._driver_cache[key] = analytic
                return analytic
            h = self.derivative_step
            values = {
                "minus_2": self._normalized_driver(y, z - 2 * h),
                "minus_1": self._normalized_driver(y, z - h),
                "center": self._normalized_driver(y, z),
                "plus_1": self._normalized_driver(y, z + h),
                "plus_2": self._normalized_driver(y, z + 2 * h),
            }
            derivative = lambda name: (
                -values["plus_2"][name]
                + 8 * values["plus_1"][name]
                - 8 * values["minus_1"][name]
                + values["minus_2"][name]
            ) / (12 * h)
            result = {
                "A": values["center"]["A"],
                "B": values["center"]["B"],
                "A_Z": derivative("A"),
                "B_Z": derivative("B"),
                "derivative_step": h,
                "derivative_method": self.driver_derivative_method,
                "analytic_driver_used": False,
                "normalized_driver": True,
            }
            self._driver_cache[key] = result
            return result

    def rhs(self, y: Any, state: list[mp.mpf], Z: Any) -> list[mp.mpf]:
        """Return the 16-component state and Z-tangent RHS."""

        with mp.workdps(self.precision):
            yy = _mp(y)
            z = _mp(Z)
            g, u = state[0], state[1]
            g_Z, u_Z = state[8], state[9]
            driver = self._driver_tangent(yy, z)
            A = driver["A"]
            B = driver["B"]
            A_Z = driver["A_Z"]
            B_Z = driver["B_Z"]
            eg = mp.exp(g)
            e2g = eg * eg
            e2y = mp.exp(2 * yy)
            ey = mp.exp(yy)
            return [
                A,
                eg * B,
                2 * e2y * eg,
                ey * u,
                2 * e2y * eg * u,
                ey * u * u,
                e2y * e2g,
                ey * e2g,
                A_Z,
                eg * (B_Z + B * g_Z),
                2 * e2y * eg * g_Z,
                ey * u_Z,
                2 * e2y * eg * (u_Z + u * g_Z),
                2 * ey * u * u_Z,
                2 * e2y * e2g * g_Z,
                2 * ey * e2g * g_Z,
            ]

    def _stress_with_shear_overrides(self, kwargs: dict[str, Any]) -> dict[str, Any]:
        """Call the shared stress API, preserving exact tiny bridge shears."""

        return evaluate_mp_stress(**kwargs)

    def evaluate(self, y: Any, Z: Any) -> dict[str, Any]:
        """Evaluate bridge fields, five moments, Z tangents and stress inputs."""

        with mp.workdps(self.precision):
            yy = _mp(y)
            z = _mp(Z)
            if not 0 <= yy <= 2 * self.hb:
                raise ValueError("Initial exit currently supports 0<=y<=2hb")
            if abs(z) >= 1:
                raise ValueError("ExitTangents requires |Z| < 1")
            init = self.initial(z)
            state = list(init["state"]) + list(init["state_Z"])
            if yy != 0:
                h = yy / self.steps

                def shifted(values: list[mp.mpf], direction: list[mp.mpf], factor: mp.mpf) -> list[mp.mpf]:
                    return [v + factor * h * w for v, w in zip(values, direction)]

                for n in range(self.steps):
                    t = n * h
                    k1 = self.rhs(t, state, z)
                    k2 = self.rhs(t + h / 2, shifted(state, k1, mp.mpf(".5")), z)
                    k3 = self.rhs(t + h / 2, shifted(state, k2, mp.mpf(".5")), z)
                    k4 = self.rhs(t + h, shifted(state, k3, mp.mpf(1)), z)
                    state = [
                        v + h * (a + 2 * b + 2 * c + d) / 6
                        for v, a, b, c, d in zip(state, k1, k2, k3, k4)
                    ]
            g, u, theta, z_moment, mixed, axial, swirl, pressure = state[:8]
            g_Z, u_Z, theta_Z, z_moment_Z, mixed_Z, axial_Z, swirl_Z, pressure_Z = state[8:]
            Fa = _mp(init["Fa"])
            Fa_Z = _mp(init["Fa_Z"])
            Fa2r2 = Fa * Fa * self.r * self.r
            R = _mp(self.r) * mp.exp(yy)
            root = mp.sqrt(2 * R)
            F = Fa * mp.exp(g)
            F_Z = mp.exp(g) * (Fa_Z + Fa * g_Z)
            driver = self._driver_tangent(yy, z)
            g_y = driver["A"]
            g_y_Z = driver["A_Z"]
            u_y = mp.exp(g) * driver["B"]
            u_y_Z = mp.exp(g) * (driver["B_Z"] + driver["B"] * g_Z)
            moments = {
                "theta": Fa * self.r * self.r * theta,
                "z": self.r * z_moment,
                "theta_z": Fa * self.r * self.r * mixed,
                "z_theta": self.r * axial - Fa2r2 * swirl,
                "p": Fa * Fa * self.r * pressure,
            }
            moments_Z = {
                "theta": self.r * self.r * (Fa_Z * theta + Fa * theta_Z),
                "z": self.r * z_moment_Z,
                "theta_z": self.r * self.r * (Fa_Z * mixed + Fa * mixed_Z),
                "z_theta": self.r * axial_Z
                - Fa2r2 * (swirl_Z + 2 * Fa_Z / Fa * swirl),
                "p": self.r * (2 * Fa * Fa_Z * pressure + Fa * Fa * pressure_Z),
            }
            P = _mp(init["P0"]) + moments["p"]
            P_Z = _mp(init["P0_Z"]) + moments_Z["p"]
            delta = _mp(self.core.delta)
            d = 1 - z * z
            L = 1 - delta * z * z
            Ur = (
                2 * z * R * u
                - (1 - delta) * z * moments["z"]
                - d * moments_Z["z"]
            ) / (L * root)
            F_R = F * g_y / R
            F_RZ = (F_Z * g_y + F * g_y_Z) / R
            Uz_R = u_y / R
            Uz_RZ = u_y_Z / R
            stress_kwargs = {
                "logR": mp.log(R),
                "Z": z,
                "delta": delta,
                "Utheta": root * F,
                "Uz": u,
                "Utheta_y": root * (F * g_y + F / 2),
                "Utheta_Z": root * F_Z,
                "Uz_y": u_y,
                "Uz_Z": u_Z,
                "moments": moments,
                "moments_Z": moments_Z,
                "P": P,
                "P_Z": P_Z,
                "shear_theta": 2 * F * g_y,
                "shear_z": root * u_y / R,
                "precision": self.precision,
            }
            stress = self._stress_with_shear_overrides(stress_kwargs)
            return {
                "R": R,
                "y": yy,
                "Z": z,
                "F": F,
                "F_Z": F_Z,
                "FZ": F_Z,
                "Uz": u,
                "Uz_Z": u_Z,
                "UZ": u_Z,
                "log_F_over_Fa": g,
                "log_F_over_Fa_Z": g_Z,
                "g_y": g_y,
                "g_y_Z": g_y_Z,
                "u_y": u_y,
                "u_y_Z": u_y_Z,
                "F_R": F_R,
                "F_RZ": F_RZ,
                "F_R_over_F": g_y / R,
                "Uz_R": Uz_R,
                "Uz_RZ": Uz_RZ,
                "moments": moments,
                "momentsZ": moments_Z,
                "moments_Z": moments_Z,
                "raw_quadratic_integrals": {
                    "axial": self.r * axial,
                    "axial_Z": self.r * axial_Z,
                    "swirl": Fa2r2 * swirl,
                    "swirl_Z": Fa2r2 * (swirl_Z + 2 * Fa_Z / Fa * swirl),
                },
                "P": P,
                "PZ": P_Z,
                "P_Z": P_Z,
                "P0_Z": init["P0_Z"],
                "Ur": Ur,
                "chi": self.multiplier(yy),
                "A": driver["A"],
                "B": driver["B"],
                "A_Z": driver["A_Z"],
                "B_Z": driver["B_Z"],
                "driver_derivative_step": driver.get("derivative_step"),
                "driver_derivative_method": driver.get(
                    "derivative_method", self.driver_derivative_method
                ),
                "analytic_driver_used": bool(driver.get("analytic_driver_used", False)),
                "normalized_driver": True,
                "stress": stress,
                "Z_moment_derivatives_available": True,
                "finite_Z_difference_used_for_moments": False,
                "cone_certified": False,
                "source_collar_constants_certified": False,
                "outer_matching_complete": False,
            }


ExitBridgeTangents = ExitTangents


__all__ = ["ExitTangents", "ExitBridgeTangents"]
