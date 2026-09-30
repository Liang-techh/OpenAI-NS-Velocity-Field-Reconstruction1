"""Finite pressure/width-ring exit switches after the Section 9.4 collar.

The provider supplies the resolved collar and its analytic exponential
continuation to ``R = 100``.  This adapter keeps the two formal parameters
and the first axial derivative in ``AxialDual`` while the two short source
switches are integrated in the dimensionless coordinate ``s = x / W``.
After ``s = 2`` the ``(a, b) = (4/5, 0)`` branch is propagated by its exact
exponential moment primitives.  The construction is a finite numerical
candidate; its pressure/width and RK remainders are deliberately open.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Mapping

import mpmath as mp

from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_pressure_width_continuation import (
    PressureWidthExitContinuation,
)
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


class PressureWidthExitSwitches:
    """Propagate the finite Section 9.4 switches in the pressure/width ring."""

    branch = "pressure_width_exit_switches"

    def __init__(self, provider: PressureWidthExitContinuation, *, steps: int = 32):
        if not isinstance(provider, PressureWidthExitContinuation):
            raise TypeError("provider must be a PressureWidthExitContinuation")
        if int(steps) != steps or steps < 8:
            raise ValueError("steps must be an integer at least 8")
        comparison = provider.comparison
        if not getattr(comparison, "automatic_Z_tangent", False):
            raise TypeError("PressureWidthExitSwitches requires an axial comparison")
        self.provider = provider
        self.comparison = comparison
        self.precision = int(provider.precision)
        self.pressure_order = int(provider.pressure_order)
        self.width_order = int(provider.width_order)
        self.steps = int(steps)
        self.R0 = mp.mpf("100")
        self.R110 = mp.mpf("110")
        self.W = comparison.W
        self.epsilon = self.W
        # ``width`` is the resolved physical coefficient of the formal W
        # atom.  It is used only to select which short phase contains a
        # requested scalar radius; all state arithmetic remains dual-valued.
        with mp.workdps(self.precision):
            self.width = mp.mpf(str(comparison.width))
            if self.width <= 0:
                raise ValueError("comparison width must be positive")
            self.R0_ring = comparison._dual(self.R0, 0)
        self._start_cache: dict[str, dict[str, Any]] = {}
        self._state_cache: dict[tuple[str, str], tuple[list[AxialDual], str]] = {}

    def _dual(self, value: Any = 0, tangent: Any = 0) -> AxialDual:
        return self.comparison._dual(value, tangent)

    def _z_key(self, Z: Any) -> str:
        return mp.nstr(mp.mpf(str(Z)), self.precision)

    @staticmethod
    def _field_tangent(source: Mapping[str, Any], name: str) -> Any:
        """Read either the new ``name_Z`` or legacy ``nameZ`` spelling."""

        for candidate in (f"{name}_Z", f"{name}Z"):
            if candidate in source:
                return source[candidate]
        return 0

    def _as_dual(self, value: Any, tangent: Any = 0) -> AxialDual:
        if isinstance(value, AxialDual):
            return value
        return self._dual(value, tangent)

    def _raw_integral_dual(
        self,
        values: Mapping[str, Any],
        tangents: Mapping[str, Any] | None,
        name: str,
    ) -> AxialDual:
        """Recover one raw quadratic integral without total-minus-Mp."""

        if name not in values:
            raise KeyError(f"raw quadratic integrals missing {name!r}")
        value = values[name]
        if isinstance(value, AxialDual):
            return value
        tangent: Any = None
        for candidate in (f"{name}_Z", f"{name}Z"):
            if candidate in values:
                tangent = values[candidate]
                break
        if tangent is None and tangents is not None:
            for candidate in (name, f"{name}_Z", f"{name}Z"):
                if candidate in tangents:
                    tangent = tangents[candidate]
                    break
        if tangent is None:
            raise KeyError(f"raw quadratic integrals missing {name}_Z")
        return self._dual(value, tangent)

    def _start(self, Z: Any) -> dict[str, Any]:
        """Build normalized dual states at the exact physical radius 100."""

        with mp.workdps(self.precision):
            z = mp.mpf(str(Z))
            key = self._z_key(z)
            cached = self._start_cache.get(key)
            if cached is not None:
                return cached
            value = self.provider.evaluate_R(self.R0, z)
            bridge_initial = self.provider.bridge.initial(z)

            def field(name: str) -> AxialDual:
                return self._as_dual(value[name], self._field_tangent(value, name))

            F0 = field("F")
            U0 = field("Uz")
            P_value = value.get("P")
            P_tangent = self._field_tangent(value, "P")
            moments = {
                name: self._as_dual(value["moments"][name], value["moments_Z"][name])
                for name in MOMENT_KEYS
            }

            raw = value.get("raw_quadratic_integrals")
            raw_Z = value.get("raw_quadratic_integrals_Z")
            if raw is None or raw_Z is None:
                raise KeyError(
                    "PressureWidthExitContinuation must expose raw_quadratic_integrals "
                    "and raw_quadratic_integrals_Z at R=100"
                )
            raw_source = "provider raw_quadratic_integrals"
            raw_axial = self._raw_integral_dual(raw, raw_Z, "axial")
            raw_swirl = self._raw_integral_dual(raw, raw_Z, "swirl")

            F0_sq = F0 * F0
            state = [
                self._dual(0),
                U0,
                moments["theta"] / (F0 * self.R0_ring * self.R0_ring),
                moments["z"] / self.R0_ring,
                moments["theta_z"] / (F0 * self.R0_ring * self.R0_ring),
                raw_axial / self.R0_ring,
                raw_swirl / (F0_sq * self.R0_ring * self.R0_ring),
                moments["p"] / (F0_sq * self.R0_ring),
            ]
            Fbar = self.provider.auxiliary.endpoint(z)["F"]
            P0 = self._as_dual(bridge_initial["P0"])
            # Retain the provider's exact start values for an endpoint check;
            # the P0 field is intentionally the raw axis pressure coefficient.
            result = {
                "z": z,
                "provider": value,
                "F0": F0,
                "U0": U0,
                "P0": P0,
                "P0_provider": self._as_dual(P_value, P_tangent),
                "moments": moments,
                "raw_axial": raw_axial,
                "raw_swirl": raw_swirl,
                "state": state,
                "Fbar": self._as_dual(Fbar),
                "raw_integrals_source": raw_source,
            }
            self._start_cache[key] = result
            return result

    def _drivers(self, s: mp.mpf, Z: mp.mpf, stage: int) -> dict[str, Any]:
        """Return a and B from the continuation EPs at ``x = s W``."""

        with mp.workdps(self.precision):
            if stage not in (1, 2):
                raise ValueError("stage must be 1 or 2")
            start = self._start(Z)
            functions = self.provider.functions(Z)
            x = self.W * s
            radius = self.R0_ring * x.exp()
            y = (radius / functions["Rb"]).log()
            D = functions["D"].evaluate(y)
            J = functions["J"].evaluate(y)
            if stage == 1:
                blend = 1 - _sigma_mp(s)
                a = self.epsilon * D
                B = (
                    -self.epsilon
                    * (start["F0"] / start["Fbar"])
                    * J
                    * blend
                )
            else:
                sigma = _sigma_mp(s - 1)
                blend = mp.mpf(0)
                a = (1 - sigma) * self.epsilon * D + mp.mpf(".8") * sigma
                B = self._dual(0)
            return {
                "x": x,
                "R": radius,
                "D": D,
                "J": J,
                "a": a,
                "B": B,
                "blend": blend,
            }

    def _rhs_s(
        self, s: mp.mpf, state: list[AxialDual], Z: mp.mpf, stage: int
    ) -> list[AxialDual]:
        """Dimensionless ``d/ds`` RHS, with the leading W kept explicit."""

        driver = self._drivers(s, Z, stage)
        W = self.W
        x = driver["x"]
        g, u = state[0], state[1]
        eg = g.exp()
        e2g = eg * eg
        ex = x.exp()
        e2x = (2 * x).exp()
        expressions = [
            -driver["a"] / 2,
            eg * driver["B"],
            2 * e2x * eg,
            ex * u,
            2 * e2x * eg * u,
            ex * u * u,
            e2x * e2g,
            ex * e2g,
        ]
        return [W * expression for expression in expressions]

    def _integrate_stage(
        self,
        left: mp.mpf,
        right: mp.mpf,
        state: list[AxialDual],
        Z: mp.mpf,
        stage: int,
    ) -> list[AxialDual]:
        with mp.workdps(self.precision):
            if right <= left:
                return list(state)
            h = (right - left) / self.steps
            values = list(state)
            for n in range(self.steps):
                s = left + n * h
                k1 = self._rhs_s(s, values, Z, stage)
                k2_state = [v + h * q / 2 for v, q in zip(values, k1)]
                k2 = self._rhs_s(s + h / 2, k2_state, Z, stage)
                k3_state = [v + h * q / 2 for v, q in zip(values, k2)]
                k3 = self._rhs_s(s + h / 2, k3_state, Z, stage)
                k4_state = [v + h * q for v, q in zip(values, k3)]
                k4 = self._rhs_s(s + h, k4_state, Z, stage)
                values = [
                    v + h * (a + 2 * b + 2 * c + d) / 6
                    for v, a, b, c, d in zip(values, k1, k2, k3, k4)
                ]
            return values

    def _constant_power(
        self, left: AxialDual, right: AxialDual, state: list[AxialDual]
    ) -> list[AxialDual]:
        """Exact normalized moments on the terminal ``a = .8, B = 0`` branch."""

        with mp.workdps(self.precision):
            length = right - left
            g, u, theta, z_moment, mixed, axial, swirl, pressure = state
            eg = g.exp()
            ex = left.exp()
            e2x = (2 * left).exp()
            em16 = (mp.mpf("1.6") * length).exp() - 1
            em1 = length.exp() - 1
            em12 = (mp.mpf("1.2") * length).exp() - 1
            em02 = (mp.mpf(".2") * length).exp() - 1
            theta_inc = 2 * e2x * eg * em16 / mp.mpf("1.6")
            z_inc = ex * u * em1
            axial_inc = ex * u * u * em1
            swirl_inc = e2x * eg * eg * em12 / mp.mpf("1.2")
            pressure_inc = ex * eg * eg * em02 / mp.mpf(".2")
            return [
                g - mp.mpf(".4") * length,
                u,
                theta + theta_inc,
                z_moment + z_inc,
                mixed + u * theta_inc,
                axial + axial_inc,
                swirl + swirl_inc,
                pressure + pressure_inc,
            ]

    def _state_at(self, radius: mp.mpf, Z: mp.mpf) -> tuple[list[AxialDual], str, mp.mpf]:
        """Integrate only the switch phases reached by a scalar radius."""

        with mp.workdps(self.precision):
            x = mp.log(radius / self.R0)
            key = (mp.nstr(x, self.precision), self._z_key(Z))
            cached = self._state_cache.get(key)
            if cached is not None:
                state, region = cached
                return list(state), region, x
            start = self._start(Z)
            state = list(start["state"])
            width = self.width
            if x <= 0:
                result = (state, "R100_endpoint")
            elif x <= width:
                state = self._integrate_stage(0, x / width, state, Z, 1)
                result = (state, "switch_1")
            elif x <= 2 * width:
                state = self._integrate_stage(0, 1, state, Z, 1)
                state = self._integrate_stage(1, x / width, state, Z, 2)
                result = (state, "switch_2")
            else:
                state = self._integrate_stage(0, 1, state, Z, 1)
                state = self._integrate_stage(1, 2, state, Z, 2)
                state = self._constant_power(
                    2 * self.W,
                    self._dual(x, 0),
                    state,
                )
                result = (state, "constant_power")
            self._state_cache[key] = result
            return list(result[0]), result[1], x

    def _drivers_at(self, x: mp.mpf, Z: mp.mpf, region: str) -> dict[str, Any]:
        with mp.workdps(self.precision):
            if region in ("R100_endpoint", "switch_1"):
                return self._drivers(x / self.width, Z, 1)
            if region == "switch_2":
                return self._drivers(x / self.width, Z, 2)
            return {
                "a": self._dual(".8"),
                "B": self._dual(0),
                "D": self._dual(0),
                "J": self._dual(0),
                "blend": mp.mpf(0),
                "x": self._dual(x, 0),
                "R": self._dual(self.R0 * mp.exp(x), 0),
            }

    @staticmethod
    def _split_tree(value: Any, tangent: bool = False) -> Any:
        if isinstance(value, AxialDual):
            return value.tangent if tangent else value.value
        if isinstance(value, dict):
            return {key: PressureWidthExitSwitches._split_tree(item, tangent) for key, item in value.items()}
        if isinstance(value, tuple):
            return tuple(PressureWidthExitSwitches._split_tree(item, tangent) for item in value)
        if isinstance(value, list):
            return [PressureWidthExitSwitches._split_tree(item, tangent) for item in value]
        return value

    def _output(
        self,
        radius: mp.mpf,
        Z: mp.mpf,
        state: list[AxialDual],
        region: str,
        x: mp.mpf,
    ) -> dict[str, Any]:
        with mp.workdps(self.precision):
            start = self._start(Z)
            F0 = start["F0"]
            g, u, theta, z_moment, mixed, axial, swirl, pressure = state
            R = self._dual(radius, 0)
            F = F0 * g.exp()
            moments = {
                "theta": F0 * self.R0_ring * self.R0_ring * theta,
                "z": self.R0_ring * z_moment,
                "theta_z": F0 * self.R0_ring * self.R0_ring * mixed,
                "z_theta": self.R0_ring * axial
                - F0 * F0 * self.R0_ring * self.R0_ring * swirl,
                "p": F0 * F0 * self.R0_ring * pressure,
            }
            P = start["P0"] + moments["p"]
            driver = self._drivers_at(x, Z, region)
            a = driver["a"]
            B = driver["B"]
            g_y = -a / 2
            u_y = g.exp() * B
            F_R = F * g_y / R
            Uz_R = u_y / R
            delta = mp.mpf(str(self.comparison.delta))
            # The five moment slots and the primary fields carry first axial
            # tangents.  Stress itself consumes those slots as value and
            # first-Z inputs; its own Z derivative would require second-Z
            # data, so the public stress receipt is the value ring.
            moment_values = {key: value.value for key, value in moments.items()}
            moment_tangents = {key: value.tangent for key, value in moments.items()}
            R_value = R.value
            z_value = self.comparison._base_jet(Z)
            root = (2 * R_value).sqrt()
            Ur_direct = (
                2 * z_value * R_value * u.value
                - (1 - delta) * z_value * moment_values["z"]
                - (1 - z_value * z_value) * moment_tangents["z"]
            ) / ((1 - delta * z_value * z_value) * root)
            stress = evaluate_mp_stress(
                mp.log(radius),
                Z,
                delta,
                Utheta=root * F.value,
                Uz=u.value,
                Utheta_y=root * (F.value * g_y.value + F.value / 2),
                Utheta_Z=root * F.tangent,
                Uz_y=u_y.value,
                Uz_Z=u.tangent,
                moments=moment_values,
                moments_Z=moment_tangents,
                P=P.value,
                P_Z=P.tangent,
                precision=self.precision,
                radius_override=R_value,
                scalar_converter=self.comparison._base_jet,
                axial_override=z_value,
                shear_theta=2 * F.value * g_y.value,
                shear_z=root * u_y.value / R_value,
                include_components=True,
            )
            raw_quadratic = {
                "axial": self.R0_ring * axial,
                "swirl": F0 * F0 * self.R0_ring * self.R0_ring * swirl,
            }
            result = {
                "R": R.value,
                "R_Z": R.tangent,
                "x": x,
                "Z": Z,
                "F": F.value,
                "FZ": F.tangent,
                "F_Z": F.tangent,
                "Uz": u.value,
                "UZ": u.tangent,
                "Uz_Z": u.tangent,
                "moments": {key: value.value for key, value in moments.items()},
                "momentsZ": {key: value.tangent for key, value in moments.items()},
                "moments_Z": {key: value.tangent for key, value in moments.items()},
                "P": P.value,
                "PZ": P.tangent,
                "P_Z": P.tangent,
                "P0": start["P0"].value,
                "P0_Z": start["P0"].tangent,
                "Ur": stress["U_r"],
                "Utheta": root * F.value,
                "Utheta_Z": root * F.tangent,
                "g_y": g_y.value,
                "g_y_Z": g_y.tangent,
                "u_y": u_y.value,
                "u_y_Z": u_y.tangent,
                "Uz_y": u_y.value,
                "Uz_y_Z": u_y.tangent,
                "F_R": F_R.value,
                "F_RZ": F_R.tangent,
                "F_R_Z": F_R.tangent,
                "Uz_R": Uz_R.value,
                "Uz_RZ": Uz_R.tangent,
                "Uz_R_Z": Uz_R.tangent,
                "I_theta": stress["I_theta"],
                "I_z": stress["I_z"],
                "D": stress["I_theta"] / F.value,
                "E": stress["I_z"] / F.value,
                "a": a.value,
                "a_Z": a.tangent,
                "B": B.value,
                "B_Z": B.tangent,
                "raw_quadratic_integrals": {
                    key: value.value for key, value in raw_quadratic.items()
                },
                "raw_quadratic_integrals_Z": {
                    key: value.tangent for key, value in raw_quadratic.items()
                },
                "stress": stress,
                "direct_Ur_dual": Ur_direct,
                "normalized_state": tuple(value.value for value in state),
                "normalized_state_Z": tuple(value.tangent for value in state),
                "region": region,
                "metadata": self.metadata(),
                "initial_raw_integrals_source": start["raw_integrals_source"],
                "finite_Z_difference_used": False,
            }
            return result

    def evaluate_R(self, R: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the finite switch branch at ``100 <= R <= 110``."""

        with mp.workdps(self.precision):
            radius = mp.mpf(str(R))
            z = mp.mpf(str(Z))
            if not self.R0 <= radius <= self.R110:
                raise ValueError("PressureWidthExitSwitches supports 100 <= R <= 110")
            if abs(z) >= 1:
                raise ValueError("PressureWidthExitSwitches requires |Z| < 1")
            state, region, x = self._state_at(radius, z)
            return self._output(radius, z, state, region, x)

    __call__ = evaluate_R

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "independent_variable": "R; short phases use s = log(R/100) / W",
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "steps": self.steps,
            "positive_epsilon_changes_retained": True,
            "finite_pressure_width_ring": True,
            "finite_Z_difference_used": False,
            "switch_RK_error_enclosed": False,
            "ode_error_enclosed": False,
            "pressure_parameter_remainder_enclosed": False,
            "width_parameter_remainder_enclosed": False,
            "errors_unenclosed": True,
            "cone_certified": False,
            "global_field_installed": False,
            "outer_matching_complete": False,
            "finite_energy_certified": False,
            "temporal_recursion_claim": False,
        }


__all__ = ["PressureWidthExitSwitches"]
