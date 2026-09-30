"""Finite pressure-width reference extension after the long reshape.

The long reshape ends at ``R_sh = 110 exp(T)``.  This adapter continues the
same component and first-``Z`` dual data through the exact reference power
law up to ``R_z = exp(log(R_ref) - 8)``.  Endpoint seeds, reshape increments,
and reference increments remain separate receipts; their public sums are
provided for downstream stress evaluation but do not replace the parts.

This is a finite analytic extension.  It does not establish global matching,
finite energy, a stress cone, or temporal recursion.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")
RAW_KEYS = ("axial", "swirl")


def _mp(value: Any) -> mp.mpf:
    """Keep an existing MP value without a decimal-string round trip."""

    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


class PressureWidthReferenceExtension:
    """Continue a pressure-width long reshape to the reference endpoint."""

    branch = "pressure_width_reference_extension"

    def __init__(self, reshape: Any, logRref: Any | None = None):
        if not hasattr(reshape, "evaluate_phase") or not hasattr(reshape, "comparison"):
            raise TypeError("reshape must expose evaluate_phase and comparison")
        comparison = reshape.comparison
        if not getattr(comparison, "automatic_Z_tangent", False):
            raise TypeError("PressureWidthReferenceExtension requires an axial comparison")
        self.reshape = reshape
        self.comparison = comparison
        self.precision = int(getattr(reshape, "precision", 80))
        self.work_precision = int(
            getattr(reshape, "work_precision", self.precision + 40)
        )
        with mp.workdps(self.work_precision):
            self.T = _mp(reshape.T)
            self.logRref = self._resolve_logRref(logRref)
            self.log_Rsh = mp.log(110) + self.T
            self.reference_end = self.logRref - mp.log(110) - 8
            if self.reference_end <= self.T:
                raise ValueError("reference endpoint must lie beyond R_sh")
            self.extension_length = self.reference_end - self.T
            self.Rsh = mp.exp(self.log_Rsh)
            self.Rz = mp.exp(self.logRref - 8)
        self._endpoint_cache: dict[str, dict[str, Any]] = {}

    def _resolve_logRref(self, explicit: Any | None) -> mp.mpf:
        bundle = getattr(self.comparison, "bundle", None)
        shared = bundle.get("shared_parameters", {}) if isinstance(bundle, Mapping) else {}
        shared_value = shared.get("logRref")
        if explicit is not None:
            value = _mp(explicit)
            if shared_value is not None and value != _mp(shared_value):
                raise ValueError("explicit logRref disagrees with shared logRref")
            return value
        if shared_value is None:
            raise ValueError("logRref is required when comparison has no shared logRref")
        return _mp(shared_value)

    def _base_jet(self, value: Any = 0) -> PressureWidthJet:
        if hasattr(self.comparison, "_base_jet"):
            return self.comparison._base_jet(value)
        return PressureWidthJet(
            value,
            pressure_order=self.comparison.pressure_order,
            width_order=self.comparison.width_order,
        )

    def _dual(self, value: Any = 0, tangent: Any = 0) -> AxialDual:
        if hasattr(self.comparison, "_dual"):
            return self.comparison._dual(value, tangent)
        return AxialDual(
            value,
            tangent,
            pressure_order=self.comparison.pressure_order,
            width_order=self.comparison.width_order,
        )

    def _z_key(self, Z: Any) -> str:
        return mp.nstr(_mp(Z), self.precision)

    @staticmethod
    def _tangent(source: Mapping[str, Any], name: str) -> Any:
        for key in (f"{name}_Z", f"{name}Z"):
            if key in source:
                return source[key]
        return 0

    def _field_dual(self, source: Mapping[str, Any], name: str) -> AxialDual:
        value = source[name]
        if isinstance(value, AxialDual):
            return value
        return self._dual(value, self._tangent(source, name))

    def _parts_dual(
        self,
        source: Mapping[str, Any],
        name: str,
        tangent_name: str,
        keys: tuple[str, ...],
    ) -> dict[str, AxialDual]:
        values = source.get(name, {})
        tangents = source.get(tangent_name, {})
        return {
            key: self._dual(values[key], tangents.get(key, 0)) for key in keys
        }

    def _endpoint_parts(self, endpoint: Mapping[str, Any]) -> dict[str, Any]:
        moment_inner = self._parts_dual(
            endpoint,
            "moment_seeds",
            "moment_seeds_Z",
            MOMENT_KEYS,
        )
        moment_reshape = self._parts_dual(
            endpoint,
            "moment_increments",
            "moment_increments_Z",
            MOMENT_KEYS,
        )
        raw_inner = self._parts_dual(
            endpoint,
            "raw_quadratic_seeds",
            "raw_quadratic_seeds_Z",
            RAW_KEYS,
        )
        raw_reshape = self._parts_dual(
            endpoint,
            "raw_quadratic_increments",
            "raw_quadratic_increments_Z",
            RAW_KEYS,
        )
        return {
            "moment_inner": moment_inner,
            "moment_reshape": moment_reshape,
            "raw_inner": raw_inner,
            "raw_reshape": raw_reshape,
        }

    @lru_cache(maxsize=64)
    def _endpoint(self, Z_key: str) -> dict[str, Any]:
        with mp.workdps(self.work_precision):
            return self.reshape.evaluate_phase(1, mp.mpf(Z_key))

    def _stress(
        self,
        *,
        log_radius: mp.mpf,
        radius: PressureWidthJet,
        z: mp.mpf,
        F: AxialDual,
        Utheta: AxialDual,
        Uz: AxialDual,
        moments: Mapping[str, AxialDual],
        pressure: AxialDual,
    ) -> dict[str, Any]:
        values = {name: moments[name].value for name in MOMENT_KEYS}
        tangents = {name: moments[name].tangent for name in MOMENT_KEYS}
        return evaluate_mp_stress(
            log_radius,
            z,
            self.comparison.delta,
            Utheta=Utheta.value,
            Uz=Uz.value,
            Utheta_y=Utheta.value * mp.mpf(".1"),
            Utheta_Z=Utheta.tangent,
            Uz_y=self._base_jet(0),
            Uz_Z=Uz.tangent,
            moments=values,
            moments_Z=tangents,
            P=pressure.value,
            P_Z=pressure.tangent,
            precision=self.work_precision,
            radius_override=radius,
            scalar_converter=self._base_jet,
            axial_override=self._base_jet(z),
            shear_theta=-mp.mpf(".8") * F.value,
            shear_z=self._base_jet(0),
            include_components=True,
        )

    @staticmethod
    def _public(duals: Mapping[str, AxialDual]) -> tuple[dict[str, Any], dict[str, Any]]:
        return (
            {key: value.value for key, value in duals.items()},
            {key: value.tangent for key, value in duals.items()},
        )

    def _zero_parts(self) -> tuple[dict[str, AxialDual], dict[str, AxialDual]]:
        moments = {key: self._dual(0) for key in MOMENT_KEYS}
        raw = {key: self._dual(0) for key in RAW_KEYS}
        return moments, raw

    def evaluate_log_offset(self, x: Any, Z: Any) -> dict[str, Any]:
        """Evaluate at ``x = log(R/R_sh)`` on the reference extension."""

        with mp.workdps(self.work_precision):
            x_value = _mp(x)
            z = _mp(Z)
            if not 0 <= x_value <= self.extension_length or abs(z) >= 1:
                raise ValueError("Require 0 <= x <= reference extension length and |Z| < 1")
            endpoint = self._endpoint(self._z_key(z))
            endpoint_parts = self._endpoint_parts(endpoint)

            Fsh = self._field_dual(endpoint, "F")
            Utheta_sh = self._field_dual(endpoint, "Utheta")
            Uz = self._field_dual(endpoint, "Uz")
            Psh = self._field_dual(endpoint, "P")
            axis_P0 = self._dual(
                endpoint.get("axis_P0", endpoint["P0"]),
                endpoint.get("axis_P0_Z", endpoint["P0_Z"]),
            )

            radius_sh = self._base_jet(self.Rsh)
            radius = radius_sh * mp.exp(x_value)
            radius_dual = self._dual(radius)
            F = Fsh * mp.exp(-mp.mpf(".4") * x_value)
            Utheta = Utheta_sh * mp.exp(mp.mpf(".1") * x_value)
            d_one = mp.expm1(x_value)
            dtheta_factor = -mp.expm1(-mp.mpf("1.6") * x_value) / mp.mpf("1.6")
            dpressure_factor = -mp.expm1(-mp.mpf(".2") * x_value) / mp.mpf(".2")
            dswirl_factor = -mp.expm1(-mp.mpf("1.2") * x_value) / mp.mpf("1.2")
            dtheta = (
                Utheta
                * radius_dual.sqrt()
                * radius_dual
                * mp.sqrt(2)
                * dtheta_factor
            )
            dpressure = Utheta * Utheta / 2 * dpressure_factor
            dswirl = radius_dual * Utheta * Utheta / 2 * dswirl_factor
            dz = Uz * self.Rsh * d_one
            dtheta_z = Uz * dtheta
            daxial = Uz * Uz * self.Rsh * d_one
            dz_theta = daxial - dswirl
            extension_moments = {
                "theta": dtheta,
                "z": dz,
                "theta_z": dtheta_z,
                "z_theta": dz_theta,
                "p": dpressure,
            }
            extension_raw = {"axial": daxial, "swirl": dswirl}

            if x_value == 0:
                # Preserve the inherited endpoint object exactly at x=0.
                total_moments = {
                    key: self._dual(endpoint["moments"][key], endpoint["moments_Z"][key])
                    for key in MOMENT_KEYS
                }
                total_raw = {
                    key: self._dual(
                        endpoint["raw_quadratic_integrals"][key],
                        endpoint["raw_quadratic_integrals_Z"][key],
                    )
                    for key in RAW_KEYS
                }
                F = self._field_dual(endpoint, "F")
                Utheta = self._field_dual(endpoint, "Utheta")
                P = self._field_dual(endpoint, "P")
                radius = self._base_jet(endpoint["R"])
                stress = endpoint["stress"]
                F_R = self._dual(endpoint.get("F_R", 0), endpoint.get("F_RZ", 0))
                Uz_R = self._dual(endpoint.get("Uz_R", 0), endpoint.get("Uz_RZ", 0))
            else:
                total_moments = {
                    key: endpoint_parts["moment_inner"][key]
                    + endpoint_parts["moment_reshape"][key]
                    + extension_moments[key]
                    for key in MOMENT_KEYS
                }
                total_raw = {
                    key: endpoint_parts["raw_inner"][key]
                    + endpoint_parts["raw_reshape"][key]
                    + extension_raw[key]
                    for key in RAW_KEYS
                }
                P = Psh + dpressure
                F_R = F * mp.mpf("-.4") / radius
                Uz_R = self._dual(0)
                stress = self._stress(
                    log_radius=self.log_Rsh + x_value,
                    radius=radius,
                    z=z,
                    F=F,
                    Utheta=Utheta,
                    Uz=Uz,
                    moments=total_moments,
                    pressure=P,
                )

            moment_public, moment_Z_public = self._public(total_moments)
            extension_public, extension_Z_public = self._public(extension_moments)
            raw_public, raw_Z_public = self._public(total_raw)
            extension_raw_public, extension_raw_Z_public = self._public(extension_raw)
            inner_public, inner_Z_public = self._public(endpoint_parts["moment_inner"])
            reshape_public, reshape_Z_public = self._public(endpoint_parts["moment_reshape"])
            raw_inner_public, raw_inner_Z_public = self._public(endpoint_parts["raw_inner"])
            raw_reshape_public, raw_reshape_Z_public = self._public(endpoint_parts["raw_reshape"])

            zero = self._base_jet(0)
            g_y = self._base_jet("-.4")
            result = {
                "R": radius,
                "R_Z": zero,
                "logR": self.log_Rsh + x_value,
                "log_radius_offset": x_value,
                "Z": z,
                "F": F.value,
                "FZ": F.tangent,
                "F_Z": F.tangent,
                "Utheta": Utheta.value,
                "Utheta_Z": Utheta.tangent,
                "Uz": Uz.value,
                "UZ": Uz.tangent,
                "Uz_Z": Uz.tangent,
                "P": P.value,
                "PZ": P.tangent,
                "P_Z": P.tangent,
                "P0": axis_P0.value,
                "P0_Z": axis_P0.tangent,
                "axis_P0": axis_P0.value,
                "axis_P0_Z": axis_P0.tangent,
                "Ur": stress["U_r"],
                "moments": moment_public,
                "momentsZ": moment_Z_public,
                "moments_Z": moment_Z_public,
                "stress": stress,
                "a": self._base_jet(".8"),
                "a_Z": zero,
                "b": zero,
                "g_y": g_y,
                "g_y_Z": zero,
                "u_y": self._base_jet(0),
                "Uz_y": zero,
                "F_R": F_R.value,
                "F_RZ": F_R.tangent,
                "F_R_Z": F_R.tangent,
                "Uz_R": Uz_R.value,
                "Uz_RZ": Uz_R.tangent,
                "Uz_R_Z": Uz_R.tangent,
                # Stable nested component receipts.
                "moment_parts": {
                    "inner_seed": inner_public,
                    "reshape": reshape_public,
                    "reference": extension_public,
                },
                "moment_parts_Z": {
                    "inner_seed": inner_Z_public,
                    "reshape": reshape_Z_public,
                    "reference": extension_Z_public,
                },
                "raw_quadratic_parts": {
                    "inner_seed": raw_inner_public,
                    "reshape": raw_reshape_public,
                    "reference": extension_raw_public,
                },
                "raw_quadratic_parts_Z": {
                    "inner_seed": raw_inner_Z_public,
                    "reshape": raw_reshape_Z_public,
                    "reference": extension_raw_Z_public,
                },
                "moment_seeds": inner_public,
                "moment_seeds_Z": inner_Z_public,
                "moment_increments": extension_public,
                "moment_increments_Z": extension_Z_public,
                "reference_moment_increments": extension_public,
                "reference_moment_increments_Z": extension_Z_public,
                "raw_quadratic_seeds": raw_inner_public,
                "raw_quadratic_seeds_Z": raw_inner_Z_public,
                "raw_quadratic_increments": extension_raw_public,
                "raw_quadratic_increments_Z": extension_raw_Z_public,
                "reference_raw_quadratic_increments": extension_raw_public,
                "reference_raw_quadratic_increments_Z": extension_raw_Z_public,
                "P_seed": endpoint["P_seed"],
                "P_seed_Z": endpoint["P_seed_Z"],
                "P_reshape_increment": endpoint["P_increment"],
                "P_reshape_increment_Z": endpoint["P_increment_Z"],
                "P_reference_increment": extension_public["p"],
                "P_reference_increment_Z": extension_Z_public["p"],
                "P_increment": extension_public["p"],
                "P_increment_Z": extension_Z_public["p"],
                "inherited_endpoint": endpoint,
                "extension_formula": {
                    "F": "F_sh exp(-.4 x)",
                    "Utheta": "Utheta_sh exp(.1 x)",
                    "Uz": "Uz_sh",
                    "P0_preserved": True,
                },
                "reference_end": self.reference_end,
                "Rsh": self.Rsh,
                "Rz": self.Rz,
                "region": "pressure_width_reference_extension",
                "metadata": self.metadata(),
            }
            return result

    def evaluate_total_log_offset(self, y: Any, Z: Any) -> dict[str, Any]:
        """Evaluate with the total coordinate ``y = log(R/110)``."""

        with mp.workdps(self.work_precision):
            y_value = _mp(y)
            if y_value < 0 or y_value > self.reference_end:
                raise ValueError("total log offset lies outside the reference domain")
            if y_value <= self.T:
                return self.reshape.evaluate_log_offset(y_value, Z)
            return self.evaluate_log_offset(y_value - self.T, Z)

    def evaluate_R(self, R: Any, Z: Any) -> dict[str, Any]:
        """Evaluate directly at a physical radius through ``R_z``."""

        with mp.workdps(self.work_precision):
            radius = _mp(R)
            if radius <= 0:
                raise ValueError("R must be positive")
            return self.evaluate_total_log_offset(mp.log(radius) - mp.log(110), Z)

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "reference_radius_log": mp.nstr(self.logRref, 40),
            "Rsh_log": mp.nstr(self.log_Rsh, 40),
            "Rz_log": mp.nstr(self.logRref - 8, 40),
            "reference_end": mp.nstr(self.reference_end, 40),
            "extension_length": mp.nstr(self.extension_length, 40),
            "pressure_order": int(self.comparison.pressure_order),
            "width_order": int(self.comparison.width_order),
            "inherited_F_Z_and_Utheta_Z": True,
            "Uz_constant": True,
            "axis_P0_preserved": True,
            "nested_seed_parts_preserved": True,
            "analytic_reference_primitives": True,
            "global_field_installed": False,
            "outer_matching_complete": False,
            "finite_energy_certified": False,
            "cone_certified": False,
            "temporal_recursion_claim": False,
            "pressure_width_remainder_enclosed": False,
        }


__all__ = ["PressureWidthReferenceExtension"]
