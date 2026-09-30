"""Finite Section 9.38 axial restoration after the reference extension.

The input is the pressure/width component-valued reference continuation.  The
restoration carries one automatic derivative in ``Z`` with :class:`AxialDual`
and keeps the endpoint, axial-restoration, and post-restoration moment parts
as separate receipts.  The implementation is a finite analytic construction;
it does not install a field or enclose the pressure-width remainder.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
    from .lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_pressure_width_long_reshape import _nodes
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_axial_primitive import _sigma_mp
    from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_pressure_width_long_reshape import _nodes


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")
RAW_KEYS = ("axial", "swirl")
PART_KEYS = ("inner_seed", "reshape", "reference")


def _mp(value: Any) -> mp.mpf:
    """Preserve an existing arbitrary-precision value."""

    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


class PressureWidthAxialRestore:
    """Restore the axial profile on ``0 <= t = log(R/Rz) <= 3``."""

    branch = "pressure_width_axial_restore"

    def __init__(self, reference: Any, *, order: int = 48):
        if not hasattr(reference, "evaluate_log_offset"):
            raise TypeError("reference must expose evaluate_log_offset")
        comparison = getattr(reference, "comparison", None)
        if comparison is None or not getattr(comparison, "automatic_Z_tangent", False):
            raise TypeError("reference must provide automatic first-Z tangents")
        if int(order) != order or order < 8:
            raise ValueError("order must be an integer at least 8")
        self.reference = reference
        self.comparison = comparison
        self.precision = int(getattr(reference, "precision", 80))
        self.work_precision = int(
            getattr(reference, "work_precision", self.precision + 40)
        )
        self.order = int(order)
        self.Rz = _mp(reference.Rz)
        self._endpoint_cache: dict[str, dict[str, Any]] = {}
        self._kernel_cache: dict[str, tuple[mp.mpf, mp.mpf, mp.mpf]] = {}
        self._restored_cache: dict[str, dict[str, Any]] = {}

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

    def _moment_duals(self, source: Mapping[str, Any]) -> dict[str, AxialDual]:
        values = source["moments"]
        tangents = source.get("moments_Z", source.get("momentsZ", {}))
        return {
            name: self._dual(values[name], tangents.get(name, 0))
            for name in MOMENT_KEYS
        }

    def _raw_duals(self, source: Mapping[str, Any]) -> dict[str, AxialDual]:
        values = source["raw_quadratic_integrals"]
        tangents = source.get("raw_quadratic_integrals_Z", {})
        result: dict[str, AxialDual] = {}
        for name in RAW_KEYS:
            tangent = tangents.get(name, values.get(f"{name}_Z", 0))
            result[name] = self._dual(values[name], tangent)
        return result

    def _part_duals(
        self,
        source: Mapping[str, Any],
        name: str,
        tangent_name: str,
        keys: tuple[str, ...],
    ) -> dict[str, dict[str, AxialDual]]:
        values = source.get(name, {})
        tangents = source.get(tangent_name, {})
        return {
            label: {
                key: self._dual(value, tangents.get(label, {}).get(key, 0))
                for key, value in entries.items()
                if key in keys
            }
            for label, entries in values.items()
        }

    @lru_cache(maxsize=64)
    def _endpoint(self, Z_key: str) -> dict[str, Any]:
        with mp.workdps(self.work_precision):
            # Pass the extension length itself so its exact endpoint survives
            # the same arbitrary-precision arithmetic as the reference.
            return self.reference.evaluate_log_offset(
                self.reference.extension_length, mp.mpf(Z_key)
            )

    @staticmethod
    def _sigma_prime(value: mp.mpf, sigma: mp.mpf) -> mp.mpf:
        if 0 < value < 1:
            return sigma * (1 - sigma) * (
                2 / value**3 + 2 / (1 - value) ** 3
            )
        return mp.mpf(0)

    def _kernels(self, t: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        """Return normalized Gauss kernels ``j1``, ``j16``, and ``j2``."""

        if t == 0:
            return mp.mpf(0), mp.mpf(0), mp.mpf(0)
        key = mp.nstr(t, self.precision)
        cached = self._kernel_cache.get(key)
        if cached is not None:
            return cached
        nodes, weights = _nodes(self.order, self.work_precision)
        half = t / 2
        j1 = mp.mpf(0)
        j16 = mp.mpf(0)
        j2 = mp.mpf(0)
        for node, weight in zip(nodes, weights):
            s = half * (node + 1)
            sigma = _sigma_mp(s)
            lag = t - s
            j1 += half * weight * mp.exp(-lag) * sigma
            j16 += half * weight * mp.exp(-mp.mpf("1.6") * lag) * sigma
            j2 += half * weight * mp.exp(-lag) * sigma * sigma
        result = (j1, j16, j2)
        self._kernel_cache[key] = result
        return result

    def _zero_moments(self) -> dict[str, AxialDual]:
        return {key: self._dual(0) for key in MOMENT_KEYS}

    def _zero_raw(self) -> dict[str, AxialDual]:
        return {key: self._dual(0) for key in RAW_KEYS}

    @staticmethod
    def _public(
        values: Mapping[str, AxialDual],
    ) -> tuple[dict[str, PressureWidthJet], dict[str, PressureWidthJet]]:
        return (
            {key: value.value for key, value in values.items()},
            {key: value.tangent for key, value in values.items()},
        )

    @staticmethod
    def _copy_parts(parts: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
        return {label: dict(values) for label, values in parts.items()}

    def _endpoint_parts(
        self, endpoint: Mapping[str, Any]
    ) -> tuple[dict[str, dict[str, AxialDual]], dict[str, dict[str, AxialDual]]]:
        moment_parts = self._part_duals(
            endpoint, "moment_parts", "moment_parts_Z", MOMENT_KEYS
        )
        raw_parts = self._part_duals(
            endpoint, "raw_quadratic_parts", "raw_quadratic_parts_Z", RAW_KEYS
        )
        for label in PART_KEYS:
            moment_parts.setdefault(label, {key: self._dual(0) for key in MOMENT_KEYS})
            raw_parts.setdefault(label, {key: self._dual(0) for key in RAW_KEYS})
        return moment_parts, raw_parts

    def _stress(
        self,
        *,
        log_radius: mp.mpf,
        radius: AxialDual,
        z: mp.mpf,
        F: AxialDual,
        Utheta: AxialDual,
        Uz: AxialDual,
        Uz_y: AxialDual,
        moments: Mapping[str, AxialDual],
        pressure: AxialDual,
    ) -> dict[str, Any]:
        values = {name: moments[name].value for name in MOMENT_KEYS}
        tangents = {name: moments[name].tangent for name in MOMENT_KEYS}
        root = radius.sqrt() * mp.sqrt(2)
        shear_z = root * Uz_y / radius
        return evaluate_mp_stress(
            log_radius,
            z,
            self.comparison.delta,
            Utheta=Utheta.value,
            Uz=Uz.value,
            Utheta_y=(Utheta * mp.mpf(".1")).value,
            Utheta_Z=Utheta.tangent,
            Uz_y=Uz_y.value,
            Uz_Z=Uz.tangent,
            moments=values,
            moments_Z=tangents,
            P=pressure.value,
            P_Z=pressure.tangent,
            precision=self.work_precision,
            radius_override=radius.value,
            scalar_converter=self._base_jet,
            axial_override=self._base_jet(z),
            shear_theta=(-mp.mpf(".8") * F).value,
            shear_z=shear_z.value,
            include_components=True,
        )

    def _base_receipt(
        self,
        endpoint: Mapping[str, Any],
        base_result: Mapping[str, Any] | None = None,
    ) -> tuple[
        dict[str, dict[str, Any]],
        dict[str, dict[str, Any]],
        dict[str, dict[str, Any]],
        dict[str, dict[str, Any]],
    ]:
        """Copy existing parts and add empty stage labels."""

        source = endpoint if base_result is None else base_result
        moment_parts = self._copy_parts(source.get("moment_parts", {}))
        moment_parts_Z = self._copy_parts(source.get("moment_parts_Z", {}))
        raw_parts = self._copy_parts(source.get("raw_quadratic_parts", {}))
        raw_parts_Z = self._copy_parts(source.get("raw_quadratic_parts_Z", {}))
        for label in PART_KEYS:
            moment_parts.setdefault(label, {key: self._base_jet(0) for key in MOMENT_KEYS})
            moment_parts_Z.setdefault(label, {key: self._base_jet(0) for key in MOMENT_KEYS})
            raw_parts.setdefault(label, {key: self._base_jet(0) for key in RAW_KEYS})
            raw_parts_Z.setdefault(label, {key: self._base_jet(0) for key in RAW_KEYS})
        return moment_parts, moment_parts_Z, raw_parts, raw_parts_Z

    def _endpoint_result(self, endpoint: Mapping[str, Any]) -> dict[str, Any]:
        """Return the exact inherited endpoint with zero restore receipts."""

        result = dict(endpoint)
        zero_mom, zero_mom_Z = self._public(self._zero_moments())
        zero_raw, zero_raw_Z = self._public(self._zero_raw())
        parts, parts_Z, raw_parts, raw_parts_Z = self._base_receipt(endpoint)
        parts["axial_restore"] = zero_mom
        parts_Z["axial_restore"] = zero_mom_Z
        parts["restored_reference"] = zero_mom
        parts_Z["restored_reference"] = zero_mom_Z
        raw_parts["axial_restore"] = zero_raw
        raw_parts_Z["axial_restore"] = zero_raw_Z
        raw_parts["restored_reference"] = zero_raw
        raw_parts_Z["restored_reference"] = zero_raw_Z
        result.update(
            moment_parts=parts,
            moment_parts_Z=parts_Z,
            raw_quadratic_parts=raw_parts,
            raw_quadratic_parts_Z=raw_parts_Z,
            axial_restore_increments=zero_mom,
            axial_restore_increments_Z=zero_mom_Z,
            restored_reference_increments=zero_mom,
            restored_reference_increments_Z=zero_mom_Z,
            stage_moment_increments=zero_mom,
            stage_moment_increments_Z=zero_mom_Z,
            axial_restore_raw_increments=zero_raw,
            axial_restore_raw_increments_Z=zero_raw_Z,
            restored_reference_raw_increments=zero_raw,
            restored_reference_raw_increments_Z=zero_raw_Z,
            stage_raw_quadratic_increments=zero_raw,
            stage_raw_quadratic_increments_Z=zero_raw_Z,
            P_axial_restore_increment=self._base_jet(0),
            P_axial_restore_increment_Z=self._base_jet(0),
            P_restored_reference_increment=self._base_jet(0),
            P_restored_reference_increment_Z=self._base_jet(0),
            inherited_endpoint=endpoint,
            axial_restore_phase=mp.mpf(0),
            region="pressure_width_axial_restore",
            metadata=self.metadata(),
            first_Z_tangent=True,
            second_Z_tangent=False,
        )
        return result

    def _stage_primitives(
        self,
        *,
        t: mp.mpf,
        Rz: AxialDual,
        Fz: AxialDual,
        Uthetaz: AxialDual,
        V0: AxialDual,
        use_sigma: bool,
    ) -> dict[str, Any]:
        R = Rz * mp.exp(t)
        F = Fz * mp.exp(-mp.mpf(".4") * t)
        Utheta = Uthetaz * mp.exp(mp.mpf(".1") * t)
        radius = R
        deltaR = Rz * mp.expm1(t)
        dtheta_factor = -mp.expm1(-mp.mpf("1.6") * t) / mp.mpf("1.6")
        dpressure_factor = -mp.expm1(-mp.mpf(".2") * t) / mp.mpf(".2")
        dswirl_factor = -mp.expm1(-mp.mpf("1.2") * t) / mp.mpf("1.2")
        dtheta = (
            Utheta
            * radius.sqrt()
            * radius
            * mp.sqrt(2)
            * dtheta_factor
        )
        dpressure = Utheta * Utheta / 2 * dpressure_factor
        dswirl = radius * Utheta * Utheta / 2 * dswirl_factor

        if use_sigma:
            sigma = _sigma_mp(t)
            sigma_prime = self._sigma_prime(t, sigma)
            j1, j16, j2 = self._kernels(t)
        else:
            sigma = mp.mpf(1)
            sigma_prime = mp.mpf(0)
            j1 = j16 = j2 = mp.mpf(0)
        return {
            "R": R,
            "F": F,
            "Utheta": Utheta,
            "deltaR": deltaR,
            "dtheta": dtheta,
            "dpressure": dpressure,
            "dswirl": dswirl,
            "sigma": sigma,
            "sigma_prime": sigma_prime,
            "j1": j1,
            "j16": j16,
            "j2": j2,
            "V0": V0,
        }

    def _stage_values(
        self,
        *,
        t: mp.mpf,
        z: mp.mpf,
        Rz: AxialDual,
        Fz: AxialDual,
        Uthetaz: AxialDual,
        V0: AxialDual,
        axial_phase: bool,
    ) -> dict[str, Any]:
        stage = self._stage_primitives(
            t=t,
            Rz=Rz,
            Fz=Fz,
            Uthetaz=Uthetaz,
            V0=V0,
            use_sigma=axial_phase,
        )
        if axial_phase:
            zdual = self._dual(z, 1)
            deltaV = 4 * zdual - V0
            V = V0 + deltaV * stage["sigma"]
            Vy = deltaV * stage["sigma_prime"]
            mass = V0 * stage["deltaR"] + deltaV * stage["R"] * stage["j1"]
            theta_z = V0 * stage["dtheta"] + deltaV * (
                mp.sqrt(2)
                * stage["R"].sqrt()
                * stage["R"]
                * stage["Utheta"]
                * stage["j16"]
            )
            axial = (
                V0 * V0 * stage["deltaR"]
                + 2 * V0 * deltaV * stage["R"] * stage["j1"]
                + deltaV * deltaV * stage["R"] * stage["j2"]
            )
        else:
            V = 4 * self._dual(z, 1)
            Vy = self._dual(0)
            mass = V * stage["deltaR"]
            theta_z = V * stage["dtheta"]
            axial = V * V * stage["deltaR"]
        moments = {
            "theta": stage["dtheta"],
            "z": mass,
            "theta_z": theta_z,
            "z_theta": axial - stage["dswirl"],
            "p": stage["dpressure"],
        }
        raw = {"axial": axial, "swirl": stage["dswirl"]}
        stage.update(V=V, Vy=Vy, moments=moments, raw=raw)
        return stage

    def _format_result(
        self,
        *,
        endpoint: Mapping[str, Any],
        base_result: Mapping[str, Any],
        stage: Mapping[str, Any],
        stage_label: str,
        t: mp.mpf,
        coordinate_increment: mp.mpf | None = None,
        z: mp.mpf,
        inherited_endpoint: Mapping[str, Any],
    ) -> dict[str, Any]:
        base_moments = self._moment_duals(base_result)
        base_raw = self._raw_duals(base_result)
        total_moments = {
            key: base_moments[key] + stage["moments"][key] for key in MOMENT_KEYS
        }
        total_raw = {
            key: base_raw[key] + stage["raw"][key] for key in RAW_KEYS
        }
        F = stage["F"]
        Utheta = stage["Utheta"]
        V = stage["V"]
        P = self._field_dual(base_result, "P") + stage["dpressure"]
        P0 = self._field_dual(base_result, "P0")
        radius = stage["R"]
        increment = t if coordinate_increment is None else coordinate_increment
        stress = self._stress(
            log_radius=_mp(base_result["logR"]) + increment,
            radius=radius,
            z=z,
            F=F,
            Utheta=Utheta,
            Uz=V,
            Uz_y=stage["Vy"],
            moments=total_moments,
            pressure=P,
        )

        moment_public, moment_Z_public = self._public(total_moments)
        raw_public, raw_Z_public = self._public(total_raw)
        stage_public, stage_Z_public = self._public(stage["moments"])
        stage_raw_public, stage_raw_Z_public = self._public(stage["raw"])
        old_parts, old_parts_Z, old_raw_parts, old_raw_parts_Z = self._base_receipt(
            endpoint, base_result
        )
        zero_mom, zero_mom_Z = self._public(self._zero_moments())
        zero_raw, zero_raw_Z = self._public(self._zero_raw())
        if stage_label == "axial_restore":
            old_parts["axial_restore"] = stage_public
            old_parts_Z["axial_restore"] = stage_Z_public
            old_raw_parts["axial_restore"] = stage_raw_public
            old_raw_parts_Z["axial_restore"] = stage_raw_Z_public
            old_parts["restored_reference"] = zero_mom
            old_parts_Z["restored_reference"] = zero_mom_Z
            old_raw_parts["restored_reference"] = zero_raw
            old_raw_parts_Z["restored_reference"] = zero_raw_Z
        else:
            old_parts["restored_reference"] = stage_public
            old_parts_Z["restored_reference"] = stage_Z_public
            old_raw_parts["restored_reference"] = stage_raw_public
            old_raw_parts_Z["restored_reference"] = stage_raw_Z_public

        axial_parts = old_parts.get("axial_restore", zero_mom)
        axial_parts_Z = old_parts_Z.get("axial_restore", zero_mom_Z)
        axial_raw_parts = old_raw_parts.get("axial_restore", zero_raw)
        axial_raw_parts_Z = old_raw_parts_Z.get("axial_restore", zero_raw_Z)
        restored_parts = old_parts.get("restored_reference", zero_mom)
        restored_parts_Z = old_parts_Z.get("restored_reference", zero_mom_Z)
        restored_raw_parts = old_raw_parts.get("restored_reference", zero_raw)
        restored_raw_parts_Z = old_raw_parts_Z.get("restored_reference", zero_raw_Z)

        result = {
            "R": radius.value,
            "R_Z": radius.tangent,
            "logR": _mp(base_result["logR"]) + increment,
            "log_radius_offset": t,
            "Z": z,
            "F": F.value,
            "FZ": F.tangent,
            "F_Z": F.tangent,
            "Utheta": Utheta.value,
            "Utheta_Z": Utheta.tangent,
            "Uz": V.value,
            "UZ": V.tangent,
            "Uz_Z": V.tangent,
            "P": P.value,
            "PZ": P.tangent,
            "P_Z": P.tangent,
            "P0": P0.value,
            "P0_Z": P0.tangent,
            "axis_P0": P0.value,
            "axis_P0_Z": P0.tangent,
            "Ur": stress["U_r"],
            "moments": moment_public,
            "momentsZ": moment_Z_public,
            "moments_Z": moment_Z_public,
            "raw_quadratic_integrals": raw_public,
            "raw_quadratic_integrals_Z": raw_Z_public,
            "stress": stress,
            "a": self._base_jet(".8"),
            "a_Z": self._base_jet(0),
            "b": self._base_jet(0),
            "g_y": self._base_jet("-.4"),
            "g_y_Z": self._base_jet(0),
            "u_y": stage["Vy"].value,
            "u_y_Z": stage["Vy"].tangent,
            "Uz_y": stage["Vy"].value,
            "Uz_y_Z": stage["Vy"].tangent,
            "Utheta_y": (Utheta * mp.mpf(".1")).value,
            "Utheta_y_Z": (Utheta * mp.mpf(".1")).tangent,
            "F_R": (F * mp.mpf("-.4") / radius).value,
            "F_RZ": (F * mp.mpf("-.4") / radius).tangent,
            "F_R_Z": (F * mp.mpf("-.4") / radius).tangent,
            "Uz_R": (stage["Vy"] / radius).value,
            "Uz_RZ": (stage["Vy"] / radius).tangent,
            "Uz_R_Z": (stage["Vy"] / radius).tangent,
            "shear_z": (radius.sqrt() * mp.sqrt(2) * stage["Vy"] / radius).value,
            "shear_z_Z": (radius.sqrt() * mp.sqrt(2) * stage["Vy"] / radius).tangent,
            "moment_parts": old_parts,
            "moment_parts_Z": old_parts_Z,
            "raw_quadratic_parts": old_raw_parts,
            "raw_quadratic_parts_Z": old_raw_parts_Z,
            "moment_seeds": endpoint.get("moment_seeds", {}).copy(),
            "moment_seeds_Z": endpoint.get("moment_seeds_Z", {}).copy(),
            "raw_quadratic_seeds": endpoint.get("raw_quadratic_seeds", {}).copy(),
            "raw_quadratic_seeds_Z": endpoint.get("raw_quadratic_seeds_Z", {}).copy(),
            "moment_increments": stage_public,
            "moment_increments_Z": stage_Z_public,
            "raw_quadratic_increments": stage_raw_public,
            "raw_quadratic_increments_Z": stage_raw_Z_public,
            "axial_restore_increments": axial_parts,
            "axial_restore_increments_Z": axial_parts_Z,
            "restored_reference_increments": restored_parts,
            "restored_reference_increments_Z": restored_parts_Z,
            "stage_moment_increments": stage_public,
            "stage_moment_increments_Z": stage_Z_public,
            "axial_restore_raw_increments": axial_raw_parts,
            "axial_restore_raw_increments_Z": axial_raw_parts_Z,
            "restored_reference_raw_increments": restored_raw_parts,
            "restored_reference_raw_increments_Z": restored_raw_parts_Z,
            "stage_raw_quadratic_increments": stage_raw_public,
            "stage_raw_quadratic_increments_Z": stage_raw_Z_public,
            "P_seed": endpoint["P_seed"],
            "P_seed_Z": endpoint["P_seed_Z"],
            "P_reshape_increment": endpoint["P_reshape_increment"],
            "P_reshape_increment_Z": endpoint["P_reshape_increment_Z"],
            "P_reference_increment": endpoint["P_reference_increment"],
            "P_reference_increment_Z": endpoint["P_reference_increment_Z"],
            "P_axial_restore_increment": axial_parts["p"],
            "P_axial_restore_increment_Z": axial_parts_Z["p"],
            "P_restored_reference_increment": restored_parts["p"],
            "P_restored_reference_increment_Z": restored_parts_Z["p"],
            "P_increment": stage_public["p"],
            "P_increment_Z": stage_Z_public["p"],
            "inherited_endpoint": inherited_endpoint,
            "axial_restore_phase": t,
            "region": stage_label,
            "normalized_kernels": {
                "j1": stage["j1"],
                "j16": stage["j16"],
                "j2": stage["j2"],
            },
            "metadata": self.metadata(),
            "first_Z_tangent": True,
            "second_Z_tangent": False,
        }
        return result

    def evaluate_phase(self, t: Any, Z: Any) -> dict[str, Any]:
        """Evaluate restoration at ``t = log(R/Rz)`` for ``0 <= t <= 3``."""

        with mp.workdps(self.work_precision):
            t_value = _mp(t)
            z = _mp(Z)
            if not 0 <= t_value <= 3 or abs(z) >= 1:
                raise ValueError("Require 0 <= t <= 3 and |Z| < 1")
            z_key = self._z_key(z)
            endpoint = self._endpoint(z_key)
            if t_value == 0:
                return self._endpoint_result(endpoint)

            if t_value <= 1:
                base_result = endpoint
                Rz = self._dual(endpoint["R"], endpoint.get("R_Z", 0))
                Fz = self._field_dual(endpoint, "F")
                Uthetaz = self._field_dual(endpoint, "Utheta")
                V0 = self._field_dual(endpoint, "Uz")
                stage = self._stage_values(
                    t=t_value,
                    z=z,
                    Rz=Rz,
                    Fz=Fz,
                    Uthetaz=Uthetaz,
                    V0=V0,
                    axial_phase=True,
                )
                return self._format_result(
                    endpoint=endpoint,
                    base_result=base_result,
                    stage=stage,
                    stage_label="axial_restore",
                    t=t_value,
                    coordinate_increment=t_value,
                    z=z,
                    inherited_endpoint=endpoint,
                )

            restored = self._restored_cache.get(z_key)
            if restored is None:
                restored = self.evaluate_phase(1, z)
                self._restored_cache[z_key] = restored
            base_result = restored
            Rz = self._dual(restored["R"], restored.get("R_Z", 0))
            Fz = self._field_dual(restored, "F")
            Uthetaz = self._field_dual(restored, "Utheta")
            V0 = 4 * self._dual(z, 1)
            x = t_value - 1
            stage = self._stage_values(
                t=x,
                z=z,
                Rz=Rz,
                Fz=Fz,
                Uthetaz=Uthetaz,
                V0=V0,
                axial_phase=False,
            )
            return self._format_result(
                endpoint=endpoint,
                base_result=base_result,
                stage=stage,
                stage_label="restored_reference",
                t=t_value,
                coordinate_increment=x,
                z=z,
                inherited_endpoint=endpoint,
            )

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "phase_domain": "0 <= log(R/Rz) <= 3",
            "restore_order": self.order,
            "pressure_order": int(self.comparison.pressure_order),
            "width_order": int(self.comparison.width_order),
            "first_Z_tangent": True,
            "second_Z_tangent": False,
            "normalized_gauss_kernels": True,
            "endpoint_parts_preserved": True,
            "pressure_axis_offset_preserved": True,
            "pressure_width_remainder_enclosed": False,
            "quadrature_error_certified": False,
            "global_field_installed": False,
            "finite_energy_certified": False,
            "cone_certified": False,
            "temporal_recursion_claim": False,
            "outer_matching_complete": False,
        }


__all__ = ["PressureWidthAxialRestore"]
