"""Frozen auxiliary continuation of the pressure/width comparison.

This adapter starts at the existing comparison endpoint ``s = 2`` and
continues the exact frozen-branch powers in an arbitrary component radius
``R``.  It is an auxiliary profile evaluator: the radius is supplied by the
caller and no prescribed continuation or global field is installed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from math import comb
from types import SimpleNamespace
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_pressure_width_axial_comparison import (
    AxialPressureWidthComparison,
)
from lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


class PressureWidthFrozenComparison:
    """Evaluate the frozen Section 9.23 auxiliary branch at arbitrary ``R``."""

    branch = "frozen_auxiliary"

    def __init__(self, comparison: Any) -> None:
        if not hasattr(comparison, "evaluate") or not hasattr(comparison, "jet"):
            raise TypeError("comparison must provide evaluate() and jet()")
        self.comparison = comparison
        for name in (
            "precision",
            "pressure_order",
            "width_order",
            "delta",
            "width",
            "h_b",
            "Ra",
        ):
            if hasattr(comparison, name):
                setattr(self, name, getattr(comparison, name))
        self._endpoint_cache: dict[str, dict[str, Any]] = {}

    def jet(self, value: Any = 0) -> Any:
        """Coerce a radius or continuation value into the comparison ring."""

        return self.comparison.jet(value)

    def _key(self, Z: Any) -> str:
        precision = int(getattr(self.comparison, "precision", 80))
        return mp.nstr(mp.mpf(str(Z)), precision)

    def endpoint(self, Z: Any) -> dict[str, Any]:
        """Return and cache the inherited comparison endpoint at ``s = 2``."""

        with mp.workdps(int(getattr(self.comparison, "precision", 80))):
            key = self._key(Z)
            cached = self._endpoint_cache.get(key)
            if cached is None:
                cached = self.comparison.evaluate(2, mp.mpf(str(Z)))
                self._endpoint_cache[key] = cached
            return cached

    def _radius(self, R: Any) -> Any:
        radius = self.jet(R)
        constant = radius.value.constant if isinstance(radius, AxialDual) else radius.constant
        if constant <= 0:
            raise ValueError("frozen comparison radius must have a positive constant")
        return radius

    def _frozen_state(self, R: Any, Z: Any) -> tuple[Any, dict[str, Any], dict[str, Any], dict[str, Any]]:
        endpoint = self.endpoint(Z)
        radius = self._radius(R)
        endpoint_radius = endpoint["R"]
        delta_radius = radius - endpoint_radius
        delta_radius_sq = radius * radius - endpoint_radius * endpoint_radius

        F = endpoint["F"]
        F_Z = endpoint.get("F_Z", endpoint["FZ"])
        U = endpoint["Uz"]
        U_Z = endpoint.get("Uz_Z", endpoint["UZ"])
        result = {
            "F": F,
            "F_Z": F_Z,
            "Uz": U,
            "Uz_Z": U_Z,
            "P": endpoint.get("P", endpoint["P"]),
            "P_Z": endpoint.get("P_Z", endpoint["PZ"]),
            "moments": dict(endpoint["moments"]),
            "moments_Z": dict(endpoint.get("moments_Z", endpoint["momentsZ"])),
        }

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
        return radius, result, endpoint, {"R_b": endpoint_radius}

    def evaluate(self, R: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the frozen branch at component/dual radius ``R``."""

        precision = int(getattr(self.comparison, "precision", 80))
        with mp.workdps(precision):
            z = mp.mpf(str(Z))
            if abs(z) >= 1:
                raise ValueError("Require |Z| < 1")
            radius, state, endpoint, branch_data = self._frozen_state(R, z)
            F = state["F"]
            F_Z = state["F_Z"]
            U = state["Uz"]
            U_Z = state["Uz_Z"]
            zero = self.jet(0)
            root = (2 * radius).sqrt()
            stress_coordinate = (
                self.comparison._stress_axial_coordinate(z)
                if hasattr(self.comparison, "_stress_axial_coordinate")
                else z
            )
            stress = evaluate_mp_stress(
                0,
                z,
                self.comparison.delta,
                Utheta=root * F,
                Uz=U,
                Utheta_y=root * (F / 2),
                Utheta_Z=root * F_Z,
                Uz_y=zero,
                Uz_Z=U_Z,
                moments=state["moments"],
                moments_Z=state["moments_Z"],
                P=state["P"],
                P_Z=state["P_Z"],
                precision=precision,
                scalar_converter=self.comparison.jet,
                radius_override=radius,
                axial_override=stress_coordinate,
            )
            I_theta = stress["I_theta"]
            I_z = stress["I_z"]
            return {
                "R": radius,
                "Rb": branch_data["R_b"],
                "R_b": branch_data["R_b"],
                "F": F,
                "FZ": F_Z,
                "F_Z": F_Z,
                "FR": zero,
                "F_R": zero,
                "Uz": U,
                "UZ": U_Z,
                "Uz_Z": U_Z,
                "UR": zero,
                "Uz_R": zero,
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
                "branch": self.branch,
                "frozen_fields": True,
                "prescribed_continuation": False,
                "global_field_installed": False,
                "endpoint": endpoint,
            }

    __call__ = evaluate

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "description": "Frozen Section 9.23 auxiliary branch, valid by definition",
            "independent_variable": "component/dual radius R supplied by caller",
            "endpoint": "comparison.evaluate(s=2, Z)",
            "frozen_prescribed_field": False,
            "prescribed_continuation_installed": False,
            "outer_connection_complete": False,
            "global_field_installed": False,
            "finite_energy_certified": False,
            "pressure_parameter_remainder_enclosed": False,
            "exit_width_remainder_enclosed": False,
            "axial_tangent_remainder_enclosed": False,
        }


def fixture() -> dict[str, Any]:
    """Resolved plain/dual replay of the frozen formulas at ``s = 3``."""

    def shifted(formal: Mapping[str, Any], offset: mp.mpf) -> dict[str, Any]:
        result = dict(formal)
        for name in ("F", "Uz", "P"):
            result[name] = [
                [
                    sum(
                        (
                            mp.mpf(comb(k, j))
                            * offset ** (k - j)
                            * row[k]
                            for k in range(j, len(row))
                        ),
                        PressurePolynomial(0),
                    )
                    for j in range(len(row))
                ]
                for row in formal[name]
            ]
        return result

    with mp.workdps(100):
        f0 = [mp.mpf(v) for v in ("2", ".3", ".02", ".001", ".0004", ".0001")]
        u0 = [mp.mpf(v) for v in (".4", ".2", ".01", ".001", ".0004", ".0001")]
        p0 = [mp.mpf(v) for v in ("-1", ".2", ".03", ".002", ".0003", ".0001")]
        z0 = mp.mpf(".3")
        formal = core_coefficients(
            z0,
            ".01",
            F0_Z_taylor=f0,
            U0_Z_taylor=u0,
            P0_Z_taylor=[PressurePolynomial({0: value}) for value in p0],
            radial_degree=3,
            precision=100,
            scalar_converter=PressurePolynomial,
        )
        axis = SimpleNamespace(precision=100, Lambda=mp.mpf(10), delta=mp.mpf(".01"))
        component = SimpleNamespace(axis=axis, precision=100, coefficients=lambda _Z: formal)
        bundle = {
            "axis": axis,
            "Lambda": axis.Lambda,
            "precision": 100,
            "component_pressure_core": component,
        }

        plain = PressureWidthComparison(bundle, h_b=mp.mpf(".001"), pressure_order=2, width_order=1)
        plain_frozen = PressureWidthFrozenComparison(plain)
        plain_direct = plain.evaluate(3, z0)
        plain_replay = plain_frozen.evaluate(plain_direct["R"], z0)

        axial = AxialPressureWidthComparison(bundle, h_b=mp.mpf(".001"), pressure_order=2, width_order=1)
        axial_frozen = PressureWidthFrozenComparison(axial)
        axial_direct = axial.evaluate(3, z0)
        axial_replay = axial_frozen.evaluate(axial_direct["R"], z0)

        errors: list[mp.mpf] = []
        for direct, replay, dual in (
            (plain_direct, plain_replay, False),
            (axial_direct, axial_replay, True),
        ):
            for name in ("F", "Uz", "P", "I_theta", "I_z", "D", "E"):
                left = direct[name]
                right = replay[name]
                if dual:
                    for slot in ("value", "tangent"):
                        a = getattr(left, slot).evaluate(0, 1)
                        b = getattr(right, slot).evaluate(0, 1)
                        errors.append(abs(a - b) / max(mp.mpf(1), abs(a), abs(b)))
                else:
                    a = left.evaluate(0, 1)
                    b = right.evaluate(0, 1)
                    errors.append(abs(a - b) / max(mp.mpf(1), abs(a), abs(b)))
        maximum = max(errors)
        if maximum >= mp.mpf("1e-60"):
            raise AssertionError(maximum)
        return {
            "available": True,
            "maximum_scaled_difference": mp.nstr(maximum, 30),
            "plain_and_axial_checked": True,
            "branch": PressureWidthFrozenComparison.branch,
            "frozen_prescribed_field": False,
            "global_field_installed": False,
        }


__all__ = ["PressureWidthFrozenComparison", "fixture"]


if __name__ == "__main__":
    print(fixture())
