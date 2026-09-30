"""Actual Section 9.4 Step 3 switches from R=100 to R=110.

The provider supplies the already constructed epsilon exit through R=100.
This module integrates the two short source switches at their exact flat-step
boundaries and then uses analytic moment primitives on the final constant
``(a,b)=(4/5,0)`` segment.  The bridge state stores six scaled cumulative
moments so that the tiny swirl amplitude at R=100 is never used as a raw
normalization for an astronomical interval.

Source: Lei--Ren, arXiv:2609.35406v1, Sections 9.4--9.5, equations
(9.25)--(9.33).  This is a finite numerical connection candidate; it does
not certify the source parameter gate, the relaxed cone globally, moment
repair, pressure matching, PDE, or time-scale recursion.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_exit_continuation import ExitContinuation
from lei_ren_part1_paper_exit_tangents import ExitTangents
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _signed_log(value: Any, precision: int = 50) -> dict[str, Any]:
    with mp.workdps(max(80, int(precision) + 20)):
        x = _mp(value)
        if x == 0:
            return {"sign": 0, "log_abs": None, "value": "0"}
        return {
            "sign": int(mp.sign(x)),
            "log_abs": mp.nstr(mp.log(abs(x)), precision),
            "value": mp.nstr(x, precision),
        }


class ExitSwitches:
    """Piecewise actual-value extension of ``ExitContinuation`` to R=110."""

    def __init__(self, provider: ExitContinuation, *, steps: int = 32) -> None:
        if not isinstance(provider, ExitContinuation):
            raise TypeError("provider must be an ExitContinuation")
        if int(steps) != steps or steps < 8:
            raise ValueError("steps must be an integer at least 8")
        self.provider = provider
        self.tangent = provider.tangent
        self.comparison = provider.comparison
        self.core = provider.core
        self.precision = int(provider.precision)
        self.r = _mp(provider.r)
        self.hb = _mp(provider.hb)
        # ExitContinuation forwards the actual source epsilon through its
        # tangent object; it intentionally has no duplicate epsilon field.
        self.epsilon = _mp(provider.tangent.epsilon)
        self.R0 = mp.mpf("100")
        self.R110 = mp.mpf("110")
        self.steps = int(steps)
        self.derivative_step = _mp(getattr(self.tangent, "derivative_step", "1e-45"))
        self._start_cache: dict[str, dict[str, Any]] = {}
        self._driver_cache: dict[tuple[str, str, str], dict[str, Any]] = {}

    def _global_y(self, R: mp.mpf) -> mp.mpf:
        with mp.workdps(self.precision):
            return mp.log(R / self.r)

    def _initial_raw_integrals(self, Z: mp.mpf) -> dict[str, mp.mpf]:
        """Recover axial/swirl integrals at R=100 from the actual provider."""

        with mp.workdps(self.precision):
            collar = self.tangent.evaluate(2 * self.tangent.hb, Z)
            Re = _mp(collar["R"])
            R0 = self.R0
            u = _mp(collar["Uz"])
            u_Z = _mp(collar["UZ"])
            F = _mp(collar["F"])
            F_Z = _mp(collar["FZ"])
            endpoint_raw = collar.get("raw_quadratic_integrals")
            if endpoint_raw is None:
                raise KeyError(
                    "ExitTangents endpoint must expose raw_quadratic_integrals "
                    "for the actual axial/swirl moment continuation"
                )
            endpoint_axial = _mp(endpoint_raw["axial"])
            endpoint_axial_Z = _mp(endpoint_raw["axial_Z"])
            endpoint_swirl = _mp(endpoint_raw["swirl"])
            endpoint_swirl_Z = _mp(endpoint_raw["swirl_Z"])
            axial = endpoint_axial + u * u * (R0 - Re)
            axial_Z = endpoint_axial_Z + 2 * u * u_Z * (R0 - Re)
            swirl = endpoint_swirl + F * F * (R0 * R0 - Re * Re) / 2
            swirl_Z = endpoint_swirl_Z + F * F_Z * (R0 * R0 - Re * Re)
            return {
                "axial": axial,
                "axial_Z": axial_Z,
                "swirl": swirl,
                "swirl_Z": swirl_Z,
            }

    @lru_cache(maxsize=64)
    def _start(self, Z_key: str) -> dict[str, Any]:
        with mp.workdps(self.precision):
            Z = mp.mpf(Z_key)
            y0 = self._global_y(self.R0)
            value = self.provider.evaluate(y0, Z)
            raw = self._initial_raw_integrals(Z)
            F0 = _mp(value["F"])
            F0_Z = _mp(value["FZ"])
            if F0 <= 0:
                raise ArithmeticError("R=100 exit has nonpositive F")
            R02 = self.R0 * self.R0
            F0R02 = F0 * R02
            F02R02 = F0 * F0 * R02
            moments = value["moments"]
            moments_Z = value["momentsZ"]
            theta = _mp(moments["theta"])
            theta_Z = _mp(moments_Z["theta"])
            theta_z = _mp(moments["theta_z"])
            theta_z_Z = _mp(moments_Z["theta_z"])
            Mz = _mp(moments["z"])
            Mz_Z = _mp(moments_Z["z"])
            Mp = _mp(moments["p"])
            Mp_Z = _mp(moments_Z["p"])
            axial = raw["axial"]
            axial_Z = raw["axial_Z"]
            swirl = raw["swirl"]
            swirl_Z = raw["swirl_Z"]
            state = [
                mp.mpf(0),
                _mp(value["Uz"]),
                theta / F0R02,
                Mz / self.R0,
                theta_z / F0R02,
                axial / self.R0,
                swirl / F02R02,
                Mp / (F0 * F0 * self.R0),
            ]
            state_Z = [
                mp.mpf(0),
                _mp(value["UZ"]),
                (theta_Z - theta * F0_Z / F0) / F0R02,
                Mz_Z / self.R0,
                (theta_z_Z - theta_z * F0_Z / F0) / F0R02,
                axial_Z / self.R0,
                (swirl_Z - 2 * F0_Z / F0 * swirl) / F02R02,
                (Mp_Z - 2 * F0_Z / F0 * Mp) / (F0 * F0 * self.R0),
            ]
            return {
                "R": self.R0,
                "F0": F0,
                "F0_Z": F0_Z,
                "P0": _mp(value["P"]) - Mp,
                "P0_Z": _mp(value["PZ"]) - Mp_Z,
                "state": state,
                "state_Z": state_Z,
                "provider": value,
                "raw_integrals": raw,
                "raw_integrals_source": "provider endpoint raw_quadratic_integrals",
            }

    def _bar_driver(self, x: mp.mpf, Z: mp.mpf, blend: mp.mpf) -> dict[str, Any]:
        """Return normalized frozen comparison drivers at R=100 exp(x)."""

        with mp.workdps(self.precision):
            R = self.R0 * mp.exp(x)
            y_global = self._global_y(R)
            bar = self.comparison.evaluate(y_global, Z)
            start = self._start(mp.nstr(Z, self.precision))
            Fbar = _mp(bar["F"])
            if Fbar == 0:
                raise ArithmeticError("comparison frozen F vanished")
            D = _mp(bar["D"])
            Iz = _mp(bar["I_z"])
            B = -self.epsilon * mp.sqrt(R / 2) * (start["F0"] / Fbar) * Iz * blend
            return {
                "R": R,
                "bar": bar,
                "D": D,
                "B": B,
                "blend": blend,
                "F0": start["F0"],
                "F0_Z": start["F0_Z"],
            }

    def _driver_tangent(self, x: mp.mpf, Z: mp.mpf, blend: mp.mpf) -> dict[str, Any]:
        """MP fourth-order centered Z derivative of normalized D and B."""

        with mp.workdps(self.precision):
            key = (
                mp.nstr(x, self.precision),
                mp.nstr(Z, self.precision),
                mp.nstr(blend, self.precision),
            )
            if key in self._driver_cache:
                return self._driver_cache[key]
            h = self.derivative_step
            center = self._bar_driver(x, Z, blend)
            R = center["R"]
            y_global = self._global_y(R)
            # Only the inexpensive frozen comparison is sampled at the four
            # stencil points. Re-evaluating the actual collar there repeats
            # the expensive 32-step provider and is unnecessary: F_100 and
            # its analytic Z tangent are already available at the center.
            values = {
                "m2": self.comparison.evaluate(y_global, Z - 2 * h),
                "m1": self.comparison.evaluate(y_global, Z - h),
                "c": center["bar"],
                "p1": self.comparison.evaluate(y_global, Z + h),
                "p2": self.comparison.evaluate(y_global, Z + 2 * h),
            }

            def derivative(name: str) -> mp.mpf:
                return (
                    -values["p2"][name]
                    + 8 * values["p1"][name]
                    - 8 * values["m1"][name]
                    + values["m2"][name]
                ) / (12 * h)

            iz = _mp(center["bar"]["I_z"])
            iz_Z = derivative("I_z")
            barF = _mp(center["bar"]["F"])
            barF_Z = _mp(center["bar"].get("FZ", center["bar"].get("F_Z")))
            ratio_Z = center["F0_Z"] / barF - center["F0"] * barF_Z / (barF * barF)
            B_Z = -self.epsilon * mp.sqrt(R / 2) * blend * (
                ratio_Z * iz + center["F0"] * iz_Z / barF
            )
            result = {
                "R": R,
                "D": center["D"],
                "B": center["B"],
                "D_Z": derivative("D"),
                "B_Z": B_Z,
                "bar": center["bar"],
                "derivative_step": h,
                "derivative_method": "MP fourth-order centered driver stencil",
                "normalized_driver": True,
            }
            self._driver_cache[key] = result
            return result

    def _rhs(self, x: mp.mpf, state: list[mp.mpf], Z: mp.mpf, stage: int) -> list[mp.mpf]:
        with mp.workdps(self.precision):
            if stage == 1:
                blend = 1 - _sigma_mp(x / self.hb)
                driver = self._driver_tangent(x, Z, blend)
                a = self.epsilon * driver["D"]
                a_Z = self.epsilon * driver["D_Z"]
                B = driver["B"]
                B_Z = driver["B_Z"]
            elif stage == 2:
                s = (x - self.hb) / self.hb
                sigma = _sigma_mp(s)
                blend = mp.mpf(0)
                driver = self._driver_tangent(x, Z, blend)
                a = (1 - sigma) * self.epsilon * driver["D"] + mp.mpf(".8") * sigma
                a_Z = (1 - sigma) * self.epsilon * driver["D_Z"]
                B = mp.mpf(0)
                B_Z = mp.mpf(0)
            else:
                raise ValueError("stage must be 1 or 2")
            g, u = state[0], state[1]
            g_Z, u_Z = state[8], state[9]
            R = self.R0 * mp.exp(x)
            eg = mp.exp(g)
            e2g = eg * eg
            e2x = mp.exp(2 * x)
            ex = mp.exp(x)
            u_x = eg * B
            u_x_Z = eg * (B_Z + B * g_Z)
            return [
                -a / 2,
                u_x,
                2 * e2x * eg,
                ex * u,
                2 * e2x * eg * u,
                ex * u * u,
                e2x * e2g,
                ex * e2g,
                -a_Z / 2,
                u_x_Z,
                -a * e2x * eg * g_Z,
                ex * u_Z,
                2 * e2x * eg * (u_Z + u * g_Z),
                2 * ex * u * u_Z,
                2 * e2x * e2g * g_Z,
                2 * ex * e2g * g_Z,
            ]

    def _integrate_stage(
        self,
        left: mp.mpf,
        right: mp.mpf,
        state: list[mp.mpf],
        Z: mp.mpf,
        stage: int,
    ) -> list[mp.mpf]:
        with mp.workdps(self.precision):
            if right <= left:
                return state
            h = (right - left) / self.steps
            values = list(state)
            for n in range(self.steps):
                x = left + n * h
                k1 = self._rhs(x, values, Z, stage)
                k2_state = [v + h * q / 2 for v, q in zip(values, k1)]
                k2 = self._rhs(x + h / 2, k2_state, Z, stage)
                k3_state = [v + h * q / 2 for v, q in zip(values, k2)]
                k3 = self._rhs(x + h / 2, k3_state, Z, stage)
                k4_state = [v + h * q for v, q in zip(values, k3)]
                k4 = self._rhs(x + h, k4_state, Z, stage)
                values = [
                    v + h * (a + 2 * b + 2 * c + d) / 6
                    for v, a, b, c, d in zip(values, k1, k2, k3, k4)
                ]
            return values

    def _final_constant_power(self, left: mp.mpf, right: mp.mpf, state: list[mp.mpf]) -> list[mp.mpf]:
        """Analytic moments for a=.8, b=0 on [2hb, log(1.1)]."""

        with mp.workdps(self.precision):
            length = right - left
            g, u, theta, z_moment, mixed, axial, swirl, pressure = state[:8]
            g_Z, u_Z, theta_Z, z_moment_Z, mixed_Z, axial_Z, swirl_Z, pressure_Z = state[8:]
            eg = mp.exp(g)
            e2x = mp.exp(2 * left)
            ex = mp.exp(left)
            theta_inc = 2 * e2x * eg * mp.expm1(mp.mpf("1.6") * length) / mp.mpf("1.6")
            z_inc = ex * u * mp.expm1(length)
            axial_inc = ex * u * u * mp.expm1(length)
            swirl_inc = e2x * eg * eg * mp.expm1(mp.mpf("1.2") * length) / mp.mpf("1.2")
            pressure_inc = ex * eg * eg * mp.expm1(mp.mpf(".2") * length) / mp.mpf(".2")
            theta_inc_Z = theta_inc * g_Z
            z_inc_Z = ex * u_Z * mp.expm1(length)
            mixed_inc = u * theta_inc
            mixed_inc_Z = u_Z * theta_inc + u * theta_inc_Z
            axial_inc_Z = 2 * u * u_Z * ex * mp.expm1(length)
            swirl_inc_Z = 2 * g_Z * swirl_inc
            pressure_inc_Z = 2 * g_Z * pressure_inc
            return [
                g - mp.mpf(".4") * length,
                u,
                theta + theta_inc,
                z_moment + z_inc,
                mixed + mixed_inc,
                axial + axial_inc,
                swirl + swirl_inc,
                pressure + pressure_inc,
                g_Z,
                u_Z,
                theta_Z + theta_inc_Z,
                z_moment_Z + z_inc_Z,
                mixed_Z + mixed_inc_Z,
                axial_Z + axial_inc_Z,
                swirl_Z + swirl_inc_Z,
                pressure_Z + pressure_inc_Z,
            ]

    def _state_at(self, x: mp.mpf, Z: mp.mpf) -> tuple[list[mp.mpf], dict[str, Any]]:
        with mp.workdps(self.precision):
            start = self._start(mp.nstr(Z, self.precision))
            state = list(start["state"]) + list(start["state_Z"])
            if x <= 0:
                return state, {"region": "R100_endpoint"}
            x1 = min(x, self.hb)
            state = self._integrate_stage(mp.mpf(0), x1, state, Z, 1)
            if x <= self.hb:
                return state, {"region": "switch_1"}
            x2 = min(x, 2 * self.hb)
            state = self._integrate_stage(self.hb, x2, state, Z, 2)
            if x <= 2 * self.hb:
                return state, {"region": "switch_2"}
            state = self._final_constant_power(2 * self.hb, x, state)
            return state, {"region": "constant_power"}

    def _drivers_at(self, x: mp.mpf, Z: mp.mpf, region: str) -> dict[str, Any]:
        with mp.workdps(self.precision):
            if region in ("R100_endpoint", "switch_1"):
                blend = 1 - _sigma_mp(x / self.hb)
                d = self._driver_tangent(x, Z, blend)
                return {
                    "a": self.epsilon * d["D"],
                    "a_Z": self.epsilon * d["D_Z"],
                    "B": d["B"],
                    "B_Z": d["B_Z"],
                    "barD": d["D"],
                    "barE": d["bar"]["E"],
                    "bar": d["bar"],
                }
            if region == "switch_2":
                s = (x - self.hb) / self.hb
                sigma = _sigma_mp(s)
                d = self._driver_tangent(x, Z, mp.mpf(0))
                return {
                    "a": (1 - sigma) * self.epsilon * d["D"] + mp.mpf(".8") * sigma,
                    "a_Z": (1 - sigma) * self.epsilon * d["D_Z"],
                    "B": mp.mpf(0),
                    "B_Z": mp.mpf(0),
                    "barD": d["D"],
                    "barE": d["bar"]["E"],
                    "bar": d["bar"],
                }
            return {
                "a": mp.mpf(".8"),
                "a_Z": mp.mpf(0),
                "B": mp.mpf(0),
                "B_Z": mp.mpf(0),
                "barD": mp.mpf(0),
                "barE": mp.mpf(0),
                "bar": None,
            }

    def _output(self, x: mp.mpf, Z: mp.mpf, state: list[mp.mpf], metadata: dict[str, Any]) -> dict[str, Any]:
        with mp.workdps(self.precision):
            start = self._start(mp.nstr(Z, self.precision))
            F0 = start["F0"]
            F0_Z = start["F0_Z"]
            g, u, theta, z_moment, mixed, axial, swirl, pressure = state[:8]
            g_Z, u_Z, theta_Z, z_moment_Z, mixed_Z, axial_Z, swirl_Z, pressure_Z = state[8:]
            R = self.R0 * mp.exp(x)
            F = F0 * mp.exp(g)
            F_Z = mp.exp(g) * (F0_Z + F0 * g_Z)
            moments = {
                "theta": F0 * self.R0**2 * theta,
                "z": self.R0 * z_moment,
                "theta_z": F0 * self.R0**2 * mixed,
                "z_theta": self.R0 * axial - F0**2 * self.R0**2 * swirl,
                "p": F0**2 * self.R0 * pressure,
            }
            moments_Z = {
                "theta": self.R0**2 * (F0_Z * theta + F0 * theta_Z),
                "z": self.R0 * z_moment_Z,
                "theta_z": self.R0**2 * (F0_Z * mixed + F0 * mixed_Z),
                "z_theta": self.R0 * axial_Z
                - F0**2 * self.R0**2 * (swirl_Z + 2 * F0_Z / F0 * swirl),
                "p": self.R0 * (2 * F0 * F0_Z * pressure + F0**2 * pressure_Z),
            }
            P = start["P0"] + moments["p"]
            P_Z = start["P0_Z"] + moments_Z["p"]
            driver = self._drivers_at(x, Z, metadata["region"])
            g_y = -driver["a"] / 2
            g_y_Z = -driver["a_Z"] / 2
            u_y = mp.exp(g) * driver["B"]
            u_y_Z = mp.exp(g) * (driver["B_Z"] + driver["B"] * g_Z)
            delta = _mp(self.core.delta)
            d = 1 - Z * Z
            L = 1 - delta * Z * Z
            root = mp.sqrt(2 * R)
            Ur = (
                2 * Z * R * u
                - (1 - delta) * Z * moments["z"]
                - d * moments_Z["z"]
            ) / (L * root)
            F_R = F * g_y / R
            F_RZ = (F_Z * g_y + F * g_y_Z) / R
            Uz_R = u_y / R
            Uz_RZ = u_y_Z / R
            stress = evaluate_mp_stress(
                mp.log(R),
                Z,
                delta,
                Utheta=root * F,
                Uz=u,
                Utheta_y=root * (F * g_y + F / 2),
                Utheta_Z=root * F_Z,
                Uz_y=u_y,
                Uz_Z=u_Z,
                moments=moments,
                moments_Z=moments_Z,
                P=P,
                P_Z=P_Z,
                shear_theta=2 * F * g_y,
                shear_z=root * u_y / R,
                precision=self.precision,
            )
            return {
                "R": R,
                "y": self._global_y(R),
                "x": x,
                "Z": Z,
                "F": F,
                "FZ": F_Z,
                "F_Z": F_Z,
                "Uz": u,
                "UZ": u_Z,
                "Uz_Z": u_Z,
                "moments": moments,
                "momentsZ": moments_Z,
                "moments_Z": moments_Z,
                "P": P,
                "PZ": P_Z,
                "P_Z": P_Z,
                "P0_Z": start["P0_Z"],
                "Ur": Ur,
                "g_y": g_y,
                "g_y_Z": g_y_Z,
                "u_y": u_y,
                "u_y_Z": u_y_Z,
                "F_R": F_R,
                "F_RZ": F_RZ,
                "Uz_R": Uz_R,
                "Uz_RZ": Uz_RZ,
                "F_R_over_F": g_y / R,
                "a": driver["a"],
                "b": 2 * driver["B"] * mp.exp(g) / (root * F),
                "barD": driver["barD"],
                "barE": driver["barE"],
                "chi": None,
                "region": metadata["region"],
                "driver_derivative_step": self.derivative_step,
                "driver_derivative_method": "MP fourth-order centered normalized D/B stencil",
                "normalized_moment_state": state,
                "raw_quadratic_integrals": {
                    "axial": self.R0 * axial,
                    "axial_Z": self.R0 * axial_Z,
                    "swirl": F0**2 * self.R0**2 * swirl,
                    "swirl_Z": F0**2 * self.R0**2 * (swirl_Z + 2*F0_Z/F0*swirl),
                },
                "initial_raw_integrals_source": start["raw_integrals_source"],
                "stress": stress,
                "initial_moment_normalization": "actual provider values; axial/swirl split from core collar primitive",
                "finite_Z_difference_used_for_moments": False,
                "source_constants_certified": False,
                "cone_certified": False,
                "moment_repair_complete": False,
                "outer_matching_complete": False,
            }

    def evaluate(self, y: Any, Z: Any) -> dict[str, Any]:
        """Evaluate at global ``y=log(R/r)`` through R=110."""

        with mp.workdps(self.precision):
            yy = _mp(y)
            z = _mp(Z)
            if abs(z) >= 1:
                raise ValueError("ExitSwitches requires |Z| < 1")
            R = self.r * mp.exp(yy)
            if R < self.R0:
                return self.provider.evaluate(yy, z)
            if R > self.R110 * (1 + mp.power(10, -self.precision + 20)):
                raise ValueError("ExitSwitches supports R <= 110")
            x = mp.log(R / self.R0)
            state, metadata = self._state_at(x, z)
            return self._output(x, z, state, metadata)

    def evaluate_R(self, R: Any, Z: Any) -> dict[str, Any]:
        """Evaluate directly at a finite physical radius in [R_a,110]."""

        with mp.workdps(self.precision):
            radius = _mp(R)
            z = _mp(Z)
            if radius <= 0:
                raise ValueError("R must be positive")
            if radius < self.R0:
                return self.provider.evaluate(self._global_y(radius), z)
            if radius > self.R110:
                raise ValueError("ExitSwitches supports R <= 110")
            if abs(z) >= 1:
                raise ValueError("ExitSwitches requires |Z| < 1")
            # Preserve exact switch radii rather than round-tripping through
            # the much larger global log(R/Ra) coordinate.
            x = mp.log(radius / self.R0)
            state, metadata = self._state_at(x, z)
            return self._output(x, z, state, metadata)

    def evaluate_switch_phase(self, stage: int, phase: Any, Z: Any) -> dict[str, Any]:
        """Resolve a subprecision collar by its explicit logarithmic phase.

        Physical radii can round to R=100 while the shear transition remains
        distinct. This API does not claim a physical grid resolves that layer.
        """
        with mp.workdps(self.precision):
            p=_mp(phase); z=_mp(Z)
            if stage not in (1,2) or not 0<=p<=1 or not abs(z)<1:
                raise ValueError('Require stage 1/2, phase in [0,1], |Z|<1')
            x=self.hb*((stage-1)+p)
            state,metadata=self._state_at(x,z)
            result=self._output(x,z,state,metadata)
            result.update({'switch_stage':stage,'switch_phase':p,
                           'log_radius_offset':x,'explicit_phase_coordinate':True,
                           'physical_radius_offset_resolved':bool(result['R']!=self.R0)})
            return result


def build_default_switches(*, precision: int = 160, steps: int = 32) -> ExitSwitches:
    """Build the current finite source candidate through the Step 3 endpoint."""

    from lei_ren_part1_paper_exit_comparison import build_comparison

    comparison = build_comparison(
        precision=precision, degree=18, Lambda="1e36", h_b=".005"
    )
    return ExitSwitches(
        ExitContinuation(ExitTangents(comparison, steps=steps)),
        steps=steps,
    )


def _relative_error(left: Any, right: Any) -> mp.mpf:
    """Return a scale-free error without converting huge exponents to float."""

    a = _mp(left)
    b = _mp(right)
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else mp.mpf(0)


def _signed_log_map(values: dict[str, Any]) -> dict[str, Any]:
    return {key: _signed_log(value) for key, value in values.items()}


def endpoint_D_test(rows):
    """For b=0 the recorded relaxed margin is D-2 (source 9.19)."""
    with mp.workdps(80):
        row=rows[-1]
        cone=row['cone']
        if _mp(row['b'])!=0 or cone['branch']!='kappa<=2':
            raise ArithmeticError('Endpoint does not have the prescribed zero axial shear')
        margin=cone['branch_margin']
        D=2+margin['sign']*mp.exp(_mp(margin['log_abs']))
        return {'method':'D recovered from recorded b=0 relaxed margin D-2',
                'D':mp.nstr(D,35),'D_gt_3':bool(D>3)}


def run(*, precision: int = 160) -> dict[str, Any]:
    """Run focused endpoint, moment, switch, and same-provider refinement checks.

    The provider tangent is held at 32 steps while the two switch integrators
    use 16 and 32 steps, isolating switch quadrature from the inherited collar
    error.  The receipt deliberately reports numerical candidate evidence only.
    """

    from lei_ren_part1_paper_exit_comparison import build_comparison

    comparison = build_comparison(
        precision=precision, degree=18, Lambda="1e36", h_b=".005"
    )
    provider = ExitContinuation(
        ExitTangents(comparison, steps=32, derivative_step="1e-45")
    )
    coarse = ExitSwitches(provider, steps=16)
    fine = ExitSwitches(provider, steps=32)
    keys = ("F", "FZ", "Uz", "UZ", "P", "PZ", "Ur", "g_y", "u_y")
    moment_keys = MOMENT_KEYS

    def compare_output(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
        errors: dict[str, Any] = {
            key: mp.nstr(_relative_error(left[key], right[key]), 35) for key in keys
        }
        errors["moments"] = {
            key: mp.nstr(
                _relative_error(left["moments"][key], right["moments"][key]), 35
            )
            for key in moment_keys
        }
        errors["momentsZ"] = {
            key: mp.nstr(
                _relative_error(left["momentsZ"][key], right["momentsZ"][key]), 35
            )
            for key in moment_keys
        }
        return errors

    def max_error(tree: Any) -> mp.mpf:
        if isinstance(tree, dict):
            return max((max_error(value) for value in tree.values()), default=mp.mpf(0))
        return _mp(tree)

    with mp.workdps(precision):
        # Construct Z only after entering the high-precision context.  A
        # fifteen-digit ``.3`` changes the Lambda-scaled source amplitude.
        Z = mp.mpf(".3")
        y100 = fine._global_y(fine.R0)
        provider100 = provider.evaluate(y100, Z)
        fine100 = fine.evaluate_R(fine.R0, Z)
        initial_errors = compare_output(fine100, provider100)
        initial_match = mp.nstr(max_error(initial_errors), 35)
        if max_error(initial_errors) > mp.mpf("1e-100"):
            raise ArithmeticError("R=100 field/moment/jet matching failed")
        initial_slope = {
            "provider_g_y": mp.nstr(provider100["g_y"], 35),
            "switch_g_y": mp.nstr(fine100["g_y"], 35),
            "provider_u_y": mp.nstr(provider100["u_y"], 35),
            "switch_u_y": mp.nstr(fine100["u_y"], 35),
            "provider_a": mp.nstr(-2 * provider100["g_y"], 35),
            "switch_a": mp.nstr(fine100["a"], 35),
            "switch_b": mp.nstr(fine100["b"], 35),
            "raw_integrals_source": fine100["initial_raw_integrals_source"],
        }

        radii = [
            fine.R0 * mp.exp(fine.hb / 2),
            fine.R0 * mp.exp(fine.hb),
            fine.R0 * mp.exp(3 * fine.hb / 2),
            fine.R0 * mp.exp(2 * fine.hb),
            mp.mpf("105"),
            mp.mpf("110"),
        ]
        rows: list[dict[str, Any]] = []
        for radius in radii:
            coarse_value = coarse.evaluate_R(radius, Z)
            fine_value = fine.evaluate_R(radius, Z)
            row: dict[str, Any] = {
                "R": mp.nstr(radius, 35),
                "region": fine_value["region"],
                "switch_refinement": compare_output(coarse_value, fine_value),
                "switch_refinement_max": mp.nstr(
                    max_error(compare_output(coarse_value, fine_value)), 35
                ),
                "F_log_abs": mp.nstr(mp.log(abs(fine_value["F"])), 35),
                "FZ_over_F": mp.nstr(fine_value["FZ"] / fine_value["F"], 35),
                "Uz": mp.nstr(fine_value["Uz"], 35),
                "UZ": mp.nstr(fine_value["UZ"], 35),
                "P": mp.nstr(fine_value["P"], 35),
                "PZ": mp.nstr(fine_value["PZ"], 35),
                "a": mp.nstr(fine_value["a"], 35),
                "b": mp.nstr(fine_value["b"], 35),
                "stress": _signed_log_map(fine_value["stress"]),
            }
            try:
                from lei_ren_part1_paper_exit_field import cone_receipt

                row["cone"] = cone_receipt(fine_value)
            except Exception as exc:  # pragma: no cover - diagnostic only
                row["cone_error"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)

        receipt: dict[str, Any] = {
            "source": "Lei--Ren arXiv:2609.35406v1 Sections 9.4--9.5",
            "source_equations": ["(9.25)--(9.33)"],
            "domain": "100 <= R <= 110, |Z| < 1",
            "Z": mp.nstr(Z, 35),
            "provider_tangent_steps": 32,
            "switch_steps": [16, 32],
            "precision": precision,
            "derivative_step": mp.nstr(fine.derivative_step, 35),
            "hb": mp.nstr(fine.hb, 35),
            "epsilon": _signed_log(fine.epsilon),
            "initial_match_max_relative_error": initial_match,
            "initial_errors": initial_errors,
            "initial_slope": initial_slope,
            "rows": rows,
            "endpoint_D_test": endpoint_D_test(rows),
            "source_constants_certified": False,
            "relaxed_cone_certified": False,
            "moment_repair_complete": False,
            "pressure_matching_complete": False,
            "pde_or_global_solution_certified": False,
            "scope_note": (
                "Finite numerical Step 3 candidate using inherited actual provider "
                "values; no source parameter gate, global cone, PDE, or closure theorem."
            ),
        }
        path = Path(__file__).with_suffix(".json")
        path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt


__all__ = ["ExitSwitches", "build_default_switches", "run"]


if __name__ == "__main__":
    result = run()
    print(
        json.dumps(
            {
                "initial_match_max_relative_error": result[
                    "initial_match_max_relative_error"
                ],
                "rows": len(result["rows"]),
            },
            indent=2,
        )
    )
