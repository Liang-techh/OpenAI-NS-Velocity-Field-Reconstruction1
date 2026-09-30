"""Reference background reconstructed from the centered source defects.

``ReferenceDefectBackground`` is the source-integral continuation on
``2 <= t = log(R/Rm) + 2 <= 3`` or, equivalently, ``1 <= x=R/Rm <= e``.
It uses the explicit five-row centered defect receipts supplied by
``CenteredComponentDefects``.  Every source label is converted to a physical
moment part and retained separately; the summed public moments are only a
convenience view and may round away atoms that are tiny relative to the
dominant reference profile.

The adapter does not solve the five-bump map, set targets, invent raw
quadratic integrals, or certify a cone, closure, finite energy, or global
field.  The axis pressure datum always comes from the centered source's
explicit ``P0`` receipt.
"""

from __future__ import annotations

from collections.abc import Mapping
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
SOURCE_ROWS = (1, 2, 3, 4, 5)


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _sqrt(value: Any) -> Any:
    method = getattr(value, "sqrt", None)
    return method() if method is not None else mp.sqrt(value)


def _pow_three_halves(value: Any) -> Any:
    return value * _sqrt(value)


def _orders(*values: Any) -> tuple[int, int] | None:
    for value in values:
        if isinstance(value, AxialDual):
            return value.orders
        if isinstance(value, PressureWidthJet):
            return value.orders
        if isinstance(value, Mapping):
            found = _orders(*value.values())
            if found is not None:
                return found
    return None


def _jet(value: Any, orders: tuple[int, int] | None) -> PressureWidthJet:
    if isinstance(value, PressureWidthJet):
        if orders is not None and value.orders != orders:
            raise ValueError(f"component jet orders {value.orders} do not match {orders}")
        return value
    kwargs = {}
    if orders is not None:
        kwargs = {"pressure_order": orders[0], "width_order": orders[1]}
    return PressureWidthJet(value, **kwargs)


def _dual(value: Any, tangent: Any = 0, orders: tuple[int, int] | None = None) -> AxialDual:
    if isinstance(value, AxialDual):
        if tangent == 0:
            return value
        return AxialDual(
            value.value,
            tangent,
            pressure_order=value.pressure_order,
            width_order=value.width_order,
        )
    return AxialDual(
        _jet(value, orders),
        _jet(tangent, orders),
        pressure_order=None if orders is None else orders[0],
        width_order=None if orders is None else orders[1],
    )


def _public(value: Any) -> tuple[Any, Any]:
    if isinstance(value, AxialDual):
        return value.value, value.tangent
    return value, _jet(0, _orders(value))


def _public_dict(values: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    result = {key: _public(value) for key, value in values.items()}
    return (
        {key: pair[0] for key, pair in result.items()},
        {key: pair[1] for key, pair in result.items()},
    )


def _field_tangent(source: Mapping[str, Any], name: str) -> Any:
    for candidate in (f"{name}_Z", f"{name}Z"):
        if candidate in source:
            return source[candidate]
    return 0


def _zero(orders: tuple[int, int] | None) -> AxialDual:
    return _dual(0, 0, orders)


def _sum_dual(values: list[AxialDual], orders: tuple[int, int] | None) -> AxialDual:
    result = _zero(orders)
    for value in values:
        result = result + value
    return result


class ReferenceDefectBackground:
    """Materialize the exact reference profile plus centered source defects."""

    branch = "five_moment_reference_defect_background"

    def __init__(self, centered_defects: Any, delta: Any = None) -> None:
        if not hasattr(centered_defects, "evaluate"):
            raise TypeError("centered_defects must expose evaluate(Z)")
        if not hasattr(centered_defects, "Rm"):
            raise AttributeError("centered_defects must expose fixed Rm")
        if not hasattr(centered_defects, "logPstar"):
            raise AttributeError("centered_defects must expose logPstar")
        self.centered_defects = centered_defects
        self.precision = int(
            getattr(centered_defects, "work_precision", getattr(centered_defects, "precision", 160))
        )
        self.work_precision = max(self.precision, int(getattr(centered_defects, "precision", self.precision)))
        with mp.workdps(self.work_precision):
            self.Rm = _mp(centered_defects.Rm)
            self.logPstar = _mp(centered_defects.logPstar)
            if delta is None:
                switches = getattr(centered_defects, 'switches', None)
                candidates = (centered_defects, switches, getattr(switches, 'comparison', None),
                              getattr(getattr(switches, 'provider', None), 'comparison', None))
                delta = next((getattr(c, 'delta') for c in candidates if c is not None and getattr(c, 'delta', None) is not None), None)
                if delta is None:
                    raise ValueError('Provide the canonical delta or a source exposing delta')
            self.delta = _mp(delta)
            if self.Rm <= 0 or not mp.isfinite(self.Rm):
                raise ValueError("centered defect Rm must be finite and positive")
            if not (0 <= self.delta < 1):
                raise ValueError("delta must satisfy 0 <= delta < 1")
            self._check_shared_delta()
        self._cache: dict[str, dict[str, Any]] = {}
        self._defect_cache: dict[str, dict[str, Any]] = {}

    def defect_data(self, Z: Any) -> dict[str, Any]:
        """Return one cached centered source evaluation at ``Z``."""

        z = _mp(Z)
        key = mp.nstr(z, self.work_precision)
        cached = self._defect_cache.get(key)
        if cached is None:
            cached = self.centered_defects.evaluate(z)
            self._defect_cache[key] = cached
        return cached

    def _check_shared_delta(self) -> None:
        switches = getattr(self.centered_defects, "switches", None)
        candidates = [
            self.centered_defects,
            switches,
            getattr(switches, "comparison", None),
            getattr(switches, "provider", None),
            getattr(getattr(switches, "provider", None), "comparison", None),
        ]
        for candidate in candidates:
            if candidate is None:
                continue
            value = candidate.get("delta") if isinstance(candidate, Mapping) else getattr(candidate, "delta", None)
            if value is not None and _mp(value) != self.delta:
                raise ValueError("reference background delta disagrees with centered source")

    def _source_parts(self, source: Mapping[str, Any], orders: tuple[int, int] | None) -> dict[str, dict[int, AxialDual]]:
        raw = source.get("parts_dual")
        if not isinstance(raw, Mapping):
            values = source.get("parts", {})
            tangents = source.get("parts_Z", {})
            raw = {
                label: {
                    row: _dual(value, tangents.get(label, {}).get(row, 0), orders)
                    for row, value in entries.items()
                }
                for label, entries in values.items()
            }
        result: dict[str, dict[int, AxialDual]] = {}
        for label, entries in raw.items():
            if not isinstance(entries, Mapping):
                continue
            result[str(label)] = {
                int(row): _dual(value, 0, orders) for row, value in entries.items()
            }
        return result

    def _physical_parts(
        self,
        source_parts: Mapping[str, Mapping[int, AxialDual]],
        *,
        Rm: Any,
        Am: AxialDual,
        Z: AxialDual,
        orders: tuple[int, int] | None,
    ) -> tuple[dict[str, dict[str, AxialDual]], dict[str, AxialDual]]:
        scale_theta = mp.sqrt(2) * _pow_three_halves(Rm) * Am
        scale_ztheta = Rm * Am**2
        result: dict[str, dict[str, AxialDual]] = {}
        total = {key: _zero(orders) for key in MOMENT_KEYS}
        for label, rows in source_parts.items():
            physical = {key: _zero(orders) for key in MOMENT_KEYS}
            for row, value in rows.items():
                if row == 1:
                    physical["z"] = physical["z"] + Rm * value
                    physical["z_theta"] = physical["z_theta"] + 8 * Z * Rm * value
                elif row == 2:
                    physical["theta_z"] = physical["theta_z"] + scale_theta * value
                elif row == 3:
                    physical["theta"] = physical["theta"] + scale_theta * value
                    physical["theta_z"] = physical["theta_z"] + 4 * Z * scale_theta * value
                elif row == 4:
                    physical["z_theta"] = physical["z_theta"] + scale_ztheta * value
                elif row == 5:
                    physical["p"] = physical["p"] + Am**2 * value
                else:
                    raise ValueError(f"unexpected centered defect row {row}")
            result[label] = physical
            for key in MOMENT_KEYS:
                total[key] = total[key] + physical[key]
        return result, total

    def _reference_moments(self, R: AxialDual, Am: AxialDual, Z: AxialDual, x: mp.mpf) -> tuple[dict[str, AxialDual], dict[str, AxialDual]]:
        u = Am * x ** mp.mpf(".1")
        theta = mp.mpf(5) / 8 * mp.sqrt(2) * _pow_three_halves(R) * u
        moments = {
            "z": 4 * Z * R,
            "theta": theta,
            "theta_z": 4 * Z * theta,
            "z_theta": 16 * Z * Z * R - mp.mpf(5) / 12 * R * u * u,
            "p": mp.mpf(5) / 2 * u * u,
        }
        return moments, {"u": u}

    def _stress(
        self,
        *,
        R: AxialDual,
        Z: AxialDual,
        u: AxialDual,
        uz: AxialDual,
        moments: Mapping[str, AxialDual],
        pressure: AxialDual,
        logR: mp.mpf,
        orders: tuple[int, int] | None,
    ) -> tuple[Any, str | None]:
        values = {key: value.value for key, value in moments.items()}
        tangents = {key: value.tangent for key, value in moments.items()}
        try:
            def convert(value: Any) -> PressureWidthJet:
                return _jet(value, orders)

            stress = evaluate_mp_stress(
                logR,
                _mp(Z.value.constant if isinstance(Z.value, PressureWidthJet) else Z.value),
                self.delta,
                Utheta=u.value,
                Uz=uz.value,
                Utheta_y=(u * mp.mpf(".1")).value,
                Utheta_Z=u.tangent,
                Uz_y=_jet(0, orders),
                Uz_Z=uz.tangent,
                moments=values,
                moments_Z=tangents,
                P=pressure.value,
                P_Z=pressure.tangent,
                precision=self.work_precision,
                radius_override=R.value,
                scalar_converter=convert,
                axial_override=convert(_mp(Z.value.constant if isinstance(Z.value, PressureWidthJet) else Z.value)),
                shear_theta=(-mp.mpf(".8") * (u / (mp.sqrt(2) * _sqrt(R)))).value,
                shear_z=_jet(0, orders),
                include_components=True,
            )
            return stress, None
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            return None, str(exc)

    def evaluate_x(self, x: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the reconstructed reference profile for ``1 <= x <= e``."""

        with mp.workdps(self.work_precision):
            x_value = _mp(x)
            if not 1 <= x_value <= mp.exp(1):
                raise ValueError("require 1 <= x <= exp(1)")
            z = _mp(Z)
            key = mp.nstr(z, self.work_precision) + ":" + mp.nstr(x_value, self.work_precision)
            cached = self._cache.get(key)
            if cached is not None:
                return cached
            source = self.defect_data(z)
            source_parts = self._source_parts(source, _orders(source.get("P0"), source.get("parts_dual", {})))
            orders = _orders(source_parts)
            zdual = _dual(z, 1, orders)
            rmdual = _dual(self.Rm, 0, orders)
            Am = mp.exp(self.logPstar - mp.mpf(".6")) / (1 + zdual * zdual)
            R = rmdual * x_value
            reference, auxiliaries = self._reference_moments(R, Am, zdual, x_value)
            parts, defect_total = self._physical_parts(
                source_parts, Rm=rmdual, Am=Am, Z=zdual, orders=orders
            )
            total_moments = {
                name: reference[name] + defect_total[name] for name in MOMENT_KEYS
            }
            # Keep the dominant analytic reference profile in the same
            # receipt map as the source defects.  The defect labels remain
            # individually addressable, while summing all physical parts now
            # reconstructs ``total_moments`` without silently dropping the
            # reference contribution.
            if "reference_power" in parts:
                raise ValueError("centered source label collides with reference_power")
            all_parts: dict[str, dict[str, AxialDual]] = {
                "reference_power": dict(reference),
                **parts,
            }
            P0 = _dual(source["P0"], source.get("P0_Z", 0), orders)
            pressure = P0 + total_moments["p"]
            u = auxiliaries["u"]
            uz = 4 * zdual
            F = u / (mp.sqrt(2) * _sqrt(R))
            stress, stress_error = self._stress(
                R=R,
                Z=zdual,
                u=u,
                uz=uz,
                moments=total_moments,
                pressure=pressure,
                logR=mp.log(self.Rm * x_value),
                orders=orders,
            )
            moments_public, moments_Z = _public_dict(total_moments)
            parts_public: dict[str, dict[str, PressureWidthJet]] = {}
            parts_Z_public: dict[str, dict[str, PressureWidthJet]] = {}
            parts_dual: dict[str, dict[str, AxialDual]] = {}
            for label, entries in all_parts.items():
                parts_dual[label] = entries
                value_map, tangent_map = _public_dict(entries)
                parts_public[label] = value_map
                parts_Z_public[label] = tangent_map
            d = source.get("d", {})
            d_Z = source.get("d_Z", {})
            result = {
                "x": x_value,
                "phase": mp.mpf(2) + mp.log(x_value),
                "R": R.value,
                "R_Z": R.tangent,
                "logR": mp.log(self.Rm * x_value),
                "Z": z,
                "F": F.value,
                "FZ": F.tangent,
                "F_Z": F.tangent,
                "Utheta": u.value,
                "Utheta_Z": u.tangent,
                "Utheta_y": (u * mp.mpf(".1")).value,
                "Utheta_y_Z": (u * mp.mpf(".1")).tangent,
                "Uz": uz.value,
                "UZ": uz.tangent,
                "Uz_Z": uz.tangent,
                "Uz_y": _jet(0, orders),
                "Uz_y_Z": _jet(0, orders),
                "F_R": (-mp.mpf(".4") * F / R).value,
                "F_R_Z": (-mp.mpf(".4") * F / R).tangent,
                "F_y": (-mp.mpf(".4") * F).value,
                "F_y_Z": (-mp.mpf(".4") * F).tangent,
                "P": pressure.value,
                "PZ": pressure.tangent,
                "P_Z": pressure.tangent,
                "P0": P0.value,
                "P0_Z": P0.tangent,
                "axis_P0": P0.value,
                "axis_P0_Z": P0.tangent,
                "moments": moments_public,
                "momentsZ": moments_Z,
                "moments_Z": moments_Z,
                "moment_parts": parts_public,
                "moment_parts_Z": parts_Z_public,
                "moment_parts_dual": parts_dual,
                "reference_moments": _public_dict(reference)[0],
                "reference_moments_Z": _public_dict(reference)[1],
                "defect_moments": _public_dict(defect_total)[0],
                "defect_moments_Z": _public_dict(defect_total)[1],
                "reference_part_label": "reference_power",
                "source_part_count": len(parts),
                "physical_part_count": len(all_parts),
                "moment_increments": {
                    key: value.value for key, value in defect_total.items()
                },
                "moment_increments_Z": {
                    key: value.tangent for key, value in defect_total.items()
                },
                "source_defects": d,
                "source_defects_Z": d_Z,
                "centered_defect_data": source,
                "source_parts": source.get("parts", {}),
                "source_parts_Z": source.get("parts_Z", {}),
                "source_parts_dual": source.get("parts_dual", {}),
                "Am": Am.value,
                "Am_Z": Am.tangent,
                "Rm": self.Rm,
                "Rm_Z": _jet(0, orders),
                "P0_receipt": source.get("P0_receipt", {}),
                "raw_quadratic_integrals": None,
                "raw_quadratic_integrals_Z": None,
                "raw_quadratic_parts": None,
                "raw_quadratic_parts_Z": None,
                "Ur": None if stress is None else stress["U_r"],
                "stress": stress,
                "stress_error": stress_error,
                "metadata": self.metadata(),
                "source_metadata": source.get("metadata", {}),
                "source_moments_may_round_tiny_atoms": True,
                "moment_parts_authoritative": True,
                "first_Z_from_same_dual": True,
            }
            if stress is None:
                root = mp.sqrt(2) * _sqrt(R)
                z_value = _mp(z)
                L = 1 - self.delta * z_value * z_value
                ur = (
                    2 * z_value * R.value * uz.value
                    - (1 - self.delta) * z_value * moments_public["z"]
                    - (1 - z_value * z_value) * moments_Z["z"]
                ) / (L * root)
                result["Ur"] = ur
            self._cache[key] = result
            return result

    def evaluate_phase(self, t: Any, Z: Any) -> dict[str, Any]:
        """Evaluate with ``t=2+log(x)`` on the reference interval ``[2,3]``."""

        t_value = _mp(t)
        if not 2 <= t_value <= 3:
            raise ValueError("require 2 <= phase t <= 3")
        return self.evaluate_x(mp.exp(t_value - 2), Z)

    __call__ = evaluate_phase

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "coordinate": "x=R/Rm; t=2+log(x)",
            "phase_domain": "2 <= t <= 3",
            "x_domain": "1 <= x <= exp(1)",
            "centered_source_integrals": True,
            "explicit_P0_preserved": True,
            "moment_parts_authoritative": True,
            "source_moments_may_round_tiny_atoms": True,
            "raw_quadratic_available": False,
            "finite_jet_quadrature_enclosed": False,
            "source_quadrature_enclosed": False,
            "five_bump_closure": False,
            "functional_closure": False,
            "cone_certified": False,
            "global_field_installed": False,
            "fixed_Rm_required": True,
            "delta": mp.nstr(self.delta, 40),
        }


__all__ = ["ReferenceDefectBackground", "MOMENT_KEYS", "SOURCE_ROWS"]
