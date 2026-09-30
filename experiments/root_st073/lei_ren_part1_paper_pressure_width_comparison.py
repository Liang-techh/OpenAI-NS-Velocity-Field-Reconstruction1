"""Section 9.23 comparison with separate pressure-tail and width jets.

The component radial core is polynomial in a formal pressure-tail parameter.
This module adds a second, independent Taylor parameter for the actual exit
width ``h_b``.  The comparison equations are the same finite auxiliary
equations as :mod:`lei_ren_part1_paper_exit_comparison`, with ``s=y/h_b`` as
the scalar control coordinate and ``R(s)=R_a exp(s W)`` retained as a jet.

This is a finite rectangular Taylor computation.  It does not install an
outer field or enclose the pressure or width remainders.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_component_pressure_core import (
    PressurePolynomial,
    evaluate_component_core_coefficients,
)
from lei_ren_part1_paper_exit_comparison import MOMENT_KEYS
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


class PressureWidthComparison:
    """Dimensionless Section 9.23 comparison in bivariate jet arithmetic."""

    def __init__(
        self,
        bundle: Any,
        pressure_order: int = 3,
        width_order: int = 2,
        h_b: Any = ".005",
        transition_steps: int = 4,
    ) -> None:
        component = (
            bundle.get("component_pressure_core")
            if isinstance(bundle, Mapping)
            else getattr(bundle, "component_pressure_core", bundle)
        )
        if component is None:
            raise ValueError("component_pressure_core is required")
        if int(transition_steps) != transition_steps or transition_steps < 1:
            raise ValueError("transition_steps must be a positive integer")
        self.bundle = bundle
        self.component_pressure_core = component
        self.pressure_order = int(pressure_order)
        self.width_order = int(width_order)
        self.transition_steps = int(transition_steps)
        bundle_precision = (
            bundle.get("precision", 0)
            if isinstance(bundle, Mapping)
            else getattr(bundle, "precision", 0)
        ) or 0
        self.precision = int(
            max(
                getattr(component, "precision", 80),
                getattr(getattr(component, "axis", None), "precision", 80),
                bundle_precision,
            )
        )
        self.axis = component.axis
        with mp.workdps(self.precision):
            self.Lambda = _mp(self.axis.Lambda)
            self.delta = _mp(self.axis.delta)
            self.Ra = mp.mpf(4) / self.Lambda
            self.width = _mp(h_b)
            if self.width <= 0:
                raise ValueError("h_b must be positive")
            self.h_b = self.width
            self.W = self.jet({(0, 1): self.width})
        self._coeff_cache: dict[str, dict[str, Any]] = {}
        self._snapshot_cache: dict[tuple[str, str], dict[str, Any]] = {}
        self._state_cache: dict[tuple[str, str], dict[str, Any]] = {}

    def jet(self, value: Any = 0) -> PressureWidthJet:
        """Construct a jet with this comparison's fixed rectangular orders."""

        if isinstance(value, PressureWidthJet):
            if value.orders != (self.pressure_order, self.width_order):
                return PressureWidthJet(
                    value.atoms,
                    pressure_order=self.pressure_order,
                    width_order=self.width_order,
                )
            return value
        return PressureWidthJet(
            value,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def _convert_coefficients(self, value: Any) -> Any:
        """Convert pressure-only atoms to ``(pressure, width=0)`` atoms."""

        if isinstance(value, PressureWidthJet):
            return self.jet(value)
        if isinstance(value, PressurePolynomial):
            return self.jet({(power, 0): coefficient for power, coefficient in value.atoms.items()})
        if isinstance(value, Mapping):
            return {key: self._convert_coefficients(item) for key, item in value.items()}
        if isinstance(value, list):
            return [self._convert_coefficients(item) for item in value]
        if isinstance(value, tuple):
            return tuple(self._convert_coefficients(item) for item in value)
        return value

    def core_coefficients(self, Z: Any) -> dict[str, Any]:
        """Return component core coefficients in bivariate jet arithmetic."""

        with mp.workdps(self.precision):
            key = mp.nstr(_mp(Z), self.precision)
            cached = self._coeff_cache.get(key)
            if cached is None:
                raw = self.component_pressure_core.coefficients(key)
                cached = self._convert_coefficients(raw)
                self._coeff_cache[key] = cached
            return cached

    def _radius(self, s: mp.mpf) -> PressureWidthJet:
        # Do not evaluate W at width = 1: retaining exp(s W) is the point of
        # this comparison.
        return self.jet(self.Ra) * (self.W * s).exp()

    @staticmethod
    def _mixed_rows(coefficients: Mapping[str, Any], radius: PressureWidthJet) -> tuple[PressureWidthJet, PressureWidthJet]:
        f_rz = PressureWidthJet(
            0,
            pressure_order=radius.pressure_order,
            width_order=radius.width_order,
        )
        uz_rz = PressureWidthJet(
            0,
            pressure_order=radius.pressure_order,
            width_order=radius.width_order,
        )
        for n, row in enumerate(coefficients["F"]):
            if n:
                f_rz += n * row[1] * radius ** (n - 1)
        for n, row in enumerate(coefficients["Uz"]):
            if n:
                uz_rz += n * row[1] * radius ** (n - 1)
        return f_rz, uz_rz

    def _evaluate_core(self, coefficients: dict[str, Any], radius: PressureWidthJet, Z: mp.mpf) -> dict[str, Any]:
        """Evaluate the shared component equations in bivariate jet form."""

        value = evaluate_component_core_coefficients(
            coefficients, radius, Z, self.delta
        )
        f_rz, uz_rz = self._mixed_rows(coefficients, radius)
        value["F_RZ"] = f_rz
        value["Uz_RZ"] = uz_rz
        value["R"] = radius
        value["Z"] = Z
        return value

    def core_snapshot(self, s: Any, Z: Any) -> dict[str, Any]:
        """Return the component core at ``R(s)`` and analytic mixed jets."""

        with mp.workdps(self.precision):
            s_value = _mp(s)
            z_value = _mp(Z)
            if s_value < 0 or abs(z_value) >= 1:
                raise ValueError("Require s >= 0 and |Z| < 1")
            key = (mp.nstr(s_value, self.precision), mp.nstr(z_value, self.precision))
            cached = self._snapshot_cache.get(key)
            if cached is None:
                radius = self._radius(s_value)
                cached = self._evaluate_core(self.core_coefficients(z_value), radius, z_value)
                self._snapshot_cache[key] = cached
            return cached

    def _state_from_core(self, core: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "F": core["F"],
            "F_Z": core["F_Z"],
            "Uz": core["Uz"],
            "Uz_Z": core["Uz_Z"],
            "P": core["P"],
            "P_Z": core["P_Z"],
            "moments": dict(core["moments"]),
            "moments_Z": dict(core["moments_Z"]),
        }

    @staticmethod
    def _copy_state(state: Mapping[str, Any]) -> dict[str, Any]:
        result = {key: state[key] for key in ("F", "F_Z", "Uz", "Uz_Z", "P", "P_Z")}
        result["moments"] = dict(state["moments"])
        result["moments_Z"] = dict(state["moments_Z"])
        return result

    def _state_add(self, *states: Mapping[str, Any]) -> dict[str, Any]:
        result = {
            key: sum((state[key] for state in states), self.jet(0))
            for key in ("F", "F_Z", "Uz", "Uz_Z", "P", "P_Z")
        }
        result["moments"] = {
            key: sum(
                (state["moments"][key] for state in states), self.jet(0)
            )
            for key in MOMENT_KEYS
        }
        result["moments_Z"] = {
            key: sum(
                (state["moments_Z"][key] for state in states), self.jet(0)
            )
            for key in MOMENT_KEYS
        }
        return result

    def _state_scale(self, state: Mapping[str, Any], scalar: Any) -> dict[str, Any]:
        factor = scalar if isinstance(scalar, PressureWidthJet) else self.jet(scalar)
        result = {
            key: factor * state[key]
            for key in ("F", "F_Z", "Uz", "Uz_Z", "P", "P_Z")
        }
        result["moments"] = {
            key: factor * state["moments"][key] for key in MOMENT_KEYS
        }
        result["moments_Z"] = {
            key: factor * state["moments_Z"][key] for key in MOMENT_KEYS
        }
        return result

    def _alpha_s(self, s: mp.mpf) -> mp.mpf:
        return 1 - _sigma_mp(s - 1)

    def _state_derivative_y(
        self, s: mp.mpf, state: Mapping[str, Any], Z: mp.mpf
    ) -> dict[str, Any]:
        """The original Section 9.23 d/dy state RHS, evaluated at s."""

        radius = self._radius(s)
        core = self.core_snapshot(s, Z)
        alpha = self._alpha_s(s)
        F_c = core["F"]
        F_c_R = core["F_R"]
        F_c_Z = core["F_Z"]
        U_c_R = core["Uz_R"]
        log_slope = radius * F_c_R / F_c
        log_slope_Z = radius * (
            core["F_RZ"] * F_c - F_c_R * F_c_Z
        ) / (F_c * F_c)
        result = {
            "F": alpha * log_slope * state["F"],
            "F_Z": alpha * (log_slope_Z * state["F"] + log_slope * state["F_Z"]),
            "Uz": alpha * radius * U_c_R,
            "Uz_Z": alpha * radius * core["Uz_RZ"],
            "P": radius * state["F"] ** 2,
            "P_Z": 2 * radius * state["F"] * state["F_Z"],
        }
        F = state["F"]
        F_Z = state["F_Z"]
        U = state["Uz"]
        U_Z = state["Uz_Z"]
        result["moments"] = {
            "theta": 2 * radius * radius * F,
            "z": radius * U,
            "theta_z": 2 * radius * radius * F * U,
            "z_theta": radius * (U * U - radius * F * F),
            "p": radius * F * F,
        }
        result["moments_Z"] = {
            "theta": 2 * radius * radius * F_Z,
            "z": radius * U_Z,
            "theta_z": 2 * radius * radius * (F_Z * U + F * U_Z),
            "z_theta": radius * (2 * U * U_Z - 2 * radius * F * F_Z),
            "p": 2 * radius * F * F_Z,
        }
        return result

    def _state_derivative_s(
        self, s: mp.mpf, state: Mapping[str, Any], Z: mp.mpf
    ) -> dict[str, Any]:
        return self._state_scale(self._state_derivative_y(s, state, Z), self.W)

    def _integrate_transition(self, s: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        start = self._state_from_core(self.core_snapshot(1, Z))
        # Keep the count in MP arithmetic.  Converting a tiny or very large
        # physical width through binary float would silently change the
        # dimensionless transition partition.
        count = max(1, int(mp.ceil((s - 1) * self.transition_steps)))
        step = (s - 1) / count
        current = start
        left = mp.mpf(1)
        for _ in range(count):
            k1 = self._state_derivative_s(left, current, Z)
            k2 = self._state_derivative_s(
                left + step / 2, self._state_add(current, self._state_scale(k1, step / 2)), Z
            )
            k3 = self._state_derivative_s(
                left + step / 2, self._state_add(current, self._state_scale(k2, step / 2)), Z
            )
            k4 = self._state_derivative_s(
                left + step, self._state_add(current, self._state_scale(k3, step)), Z
            )
            current = self._state_add(
                current,
                self._state_scale(
                    self._state_add(k1, self._state_scale(k2, 2), self._state_scale(k3, 2), k4),
                    step / 6,
                ),
            )
            left += step
        return current

    def _freeze_after_transition(
        self, s: mp.mpf, endpoint: Mapping[str, Any], Z: mp.mpf
    ) -> dict[str, Any]:
        radius = self._radius(s)
        endpoint_radius = self._radius(2)
        delta_radius = radius - endpoint_radius
        delta_radius_sq = radius * radius - endpoint_radius * endpoint_radius
        F = endpoint["F"]
        F_Z = endpoint["F_Z"]
        U = endpoint["Uz"]
        U_Z = endpoint["Uz_Z"]
        result = self._copy_state(endpoint)
        result["P"] += delta_radius * F * F
        result["P_Z"] += delta_radius * 2 * F * F_Z
        moments = result["moments"]
        moments_Z = result["moments_Z"]
        moments["theta"] += F * delta_radius_sq
        moments_Z["theta"] += F_Z * delta_radius_sq
        moments["z"] += U * delta_radius
        moments_Z["z"] += U_Z * delta_radius
        moments["theta_z"] += F * U * delta_radius_sq
        moments_Z["theta_z"] += (F_Z * U + F * U_Z) * delta_radius_sq
        moments["z_theta"] += U * U * delta_radius - F * F * delta_radius_sq / 2
        moments_Z["z_theta"] += (
            2 * U * U_Z * delta_radius - F * F_Z * delta_radius_sq
        )
        moments["p"] += F * F * delta_radius
        moments_Z["p"] += 2 * F * F_Z * delta_radius
        return result

    def _state_at(self, s: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        key = (mp.nstr(s, self.precision), mp.nstr(Z, self.precision))
        cached = self._state_cache.get(key)
        if cached is not None:
            return self._copy_state(cached)
        if s <= 1:
            result = self._state_from_core(self.core_snapshot(s, Z))
        elif s <= 2:
            result = self._integrate_transition(s, Z)
        else:
            endpoint = self._integrate_transition(mp.mpf(2), Z)
            result = self._freeze_after_transition(s, endpoint, Z)
        self._state_cache[key] = self._copy_state(result)
        return result

    def _derivatives(
        self, s: mp.mpf, state: Mapping[str, Any], Z: mp.mpf
    ) -> tuple[PressureWidthJet, PressureWidthJet]:
        if s <= 1:
            core = self.core_snapshot(s, Z)
            alpha = self._alpha_s(s)
            F_R = alpha * core["F_R"] * state["F"] / core["F"]
            Uz_R = alpha * core["Uz_R"]
            return F_R, Uz_R
        if s <= 2:
            derivative = self._state_derivative_y(s, state, Z)
            radius = self._radius(s)
            return derivative["F"] / radius, derivative["Uz"] / radius
        return self.jet(0), self.jet(0)

    def evaluate(self, s: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the comparison profile at scalar ``s=y/h_b``."""

        with mp.workdps(self.precision):
            s_value = _mp(s)
            z = _mp(Z)
            if s_value < 0 or abs(z) >= 1:
                raise ValueError("Require s >= 0 and |Z| < 1")
            radius = self._radius(s_value)
            state = self._state_at(s_value, z)
            F_R, Uz_R = self._derivatives(s_value, state, z)
            F = state["F"]
            F_Z = state["F_Z"]
            U = state["Uz"]
            U_Z = state["Uz_Z"]
            root = (2 * radius).sqrt()
            stress = evaluate_mp_stress(
                0,
                z,
                self.delta,
                Utheta=root * F,
                Uz=U,
                Utheta_y=root * (radius * F_R + F / 2),
                Utheta_Z=root * F_Z,
                Uz_y=radius * Uz_R,
                Uz_Z=U_Z,
                moments=state["moments"],
                moments_Z=state["moments_Z"],
                P=state["P"],
                P_Z=state["P_Z"],
                precision=self.precision,
                scalar_converter=self.jet,
                radius_override=radius,
            )
            I_theta = stress["I_theta"]
            I_z = stress["I_z"]
            return {
                "s": s_value,
                "R": radius,
                "F": F,
                "FZ": F_Z,
                "F_Z": F_Z,
                "FR": F_R,
                "F_R": F_R,
                "Uz": U,
                "UZ": U_Z,
                "Uz_Z": U_Z,
                "UR": Uz_R,
                "Uz_R": Uz_R,
                "P": state["P"],
                "PZ": state["P_Z"],
                "P_Z": state["P_Z"],
                "moments": state["moments"],
                "momentsZ": state["moments_Z"],
                "moments_Z": state["moments_Z"],
                "I_theta": I_theta,
                "I_z": I_z,
                "D": I_theta / F,
                "E": I_z / F,
                "stress": stress,
                "core_scope": "component_pressure_core_with_width_jet",
                "pressure_order": self.pressure_order,
                "width_order": self.width_order,
                "pressure_parameter_order_truncated": True,
                "exit_width_order_truncated": True,
                "ode_remainder_enclosed": False,
                "quadrature_error_enclosed": False,
                "global_field_installed": False,
            }

    __call__ = evaluate

    def metadata(self) -> dict[str, Any]:
        return {
            "comparison": "Section 9.23 auxiliary equations in s=y/h_b",
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "h_b": mp.nstr(self.width, self.precision),
            "R_a": mp.nstr(self.Ra, self.precision),
            "transition_steps": self.transition_steps,
            "pressure_parameter_order_truncated": True,
            "exit_width_order_truncated": True,
            "rectangular_jet_remainder_enclosed": False,
            "ode_remainder_enclosed": False,
            "quadrature_error_enclosed": False,
            "outer_connection_complete": False,
        }


def fixture() -> dict[str, Any]:
    """Resolved-family smoke fixture using the component-core test fixture."""

    return {
        "available": True,
        "scope": "API fixture requires a caller-supplied resolved component core",
        "pressure_parameter_order_truncated": True,
        "exit_width_order_truncated": True,
        "ode_remainder_enclosed": False,
    }
