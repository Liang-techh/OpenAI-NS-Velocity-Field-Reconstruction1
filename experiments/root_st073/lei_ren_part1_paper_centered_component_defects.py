"""Centered five-row component defects after the finite exit switches.

``CenteredComponentDefects`` is a local algebraic assembly at the physical
switch endpoint ``R = 110``.  It keeps the pressure-tail and exit-width atoms
in ``PressureWidthJet`` and carries the first axial derivative with
``AxialDual``.  The five rows are returned both as dual rows and as separate
value / Z-tangent jet maps.  Every subtractive constant and every flat-shape
Taylor term has its own receipt in ``parts`` and ``parts_Z``.

The assembly is intentionally local.  It preserves the explicit axis
pressure datum from the switch source and records that closure, smallness,
and error-enclosure claims are outside this adapter's scope.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
    from .lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual  # type: ignore
    from lei_ren_part1_paper_axial_primitive import _sigma_mp  # type: ignore
    from lei_ren_part1_paper_flat_shape_component import (  # type: ignore
        compose_flat_shape_defect,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet  # type: ignore


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")
KERNEL_SPECS = {
    "theta": ("1.6", "1"),
    "pressure": (".2", "2"),
    "energy": ("1.2", "2"),
}


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _digits(value: mp.mpf) -> int:
    if value <= 1:
        return 1
    return max(1, int(mp.floor(mp.log10(value))) + 1)


def _field_tangent(source: Mapping[str, Any], name: str) -> Any:
    for candidate in (f"{name}_Z", f"{name}Z"):
        if candidate in source:
            return source[candidate]
    raise KeyError(f"switch output is missing {name}_Z")


class CenteredComponentDefects:
    """Assemble the centered five-row defect over the component ring."""

    branch = "centered_component_defects"
    independent_variable = "Z; rows are centered component receipts at R=110"

    def __init__(
        self,
        switches: Any,
        A: Any = "1e150",
        logC: Any = "5e151",
        logPstar: Any = "14",
        *,
        order: int = 32,
        window: Any = 24,
        restore_order: int = 64,
    ) -> None:
        if not hasattr(switches, "evaluate_R"):
            raise TypeError("switches must expose evaluate_R(R, Z)")
        if int(order) != order or int(order) < 8:
            raise ValueError("order must be an integer at least 8")
        if int(restore_order) != restore_order or int(restore_order) < 8:
            raise ValueError("restore_order must be an integer at least 8")
        self.switches = switches
        self.order = int(order)
        self.window = _mp(window)
        self.restore_order = int(restore_order)
        self.precision = int(getattr(switches, "precision", 260))
        self.pressure_order = int(getattr(switches, "pressure_order", 9))
        self.width_order = int(getattr(switches, "width_order", 2))
        # Parse decimal controls inside a high precision context.  This
        # avoids default-context rounding of the very large production logs.
        with mp.workdps(max(320, self.precision)):
            self.A = _mp(A)
            self.logC = _mp(logC)
            self.logPstar = _mp(logPstar)
            if not (self.A > 0 and mp.isfinite(self.A)):
                raise ValueError("A must be a finite positive scalar")
            if not (mp.isfinite(self.logC) and mp.isfinite(self.logPstar)):
                raise ValueError("logC and logPstar must be finite")
            self.T = 400 * self.A
            rm_log = 10 * (self.logC + self.logPstar) - 6
            if self.T >= rm_log - 2:
                raise ValueError(
                    "require T < 10*(logC + logPstar) - 8 so Rsh < Rz"
                )
            self._check_shared_parameters()
        self.work_precision = max(
            320 + _digits(self.T),
            self.precision + _digits(self.T) + 30,
        )
        with mp.workdps(self.work_precision):
            rm_log = 10 * (self.logC + self.logPstar) - 6
            self.Rm = 110 * mp.exp(rm_log)
            # Use exponent differences directly; forming two enormous
            # exponentials and dividing them would spend precision on a
            # cancellation that is algebraically exact.
            self.rho0 = mp.exp(-rm_log)
            self.alpha = mp.exp(self.T - rm_log)
            self.Rsh = 110 * mp.exp(self.T)
            self.Rz = 110 * mp.exp(rm_log - 2)
            self.K1, self.K16, self.K2 = self._restore_constants()
        self._z_cache: dict[str, dict[str, Any]] = {}

    def _check_shared_parameters(self) -> None:
        """Reject a switch bundle whose shared log controls disagree."""

        candidates = [self.switches, getattr(self.switches, "provider", None)]
        comparison = getattr(self.switches, "comparison", None)
        provider = getattr(self.switches, "provider", None)
        provider_comparison = getattr(provider, "comparison", None)
        candidates.extend(
            (
                comparison,
                getattr(comparison, "bundle", None),
                provider_comparison,
                getattr(provider_comparison, "bundle", None),
            )
        )
        for candidate in candidates:
            if candidate is None:
                continue
            shared = (
                candidate.get("shared_parameters")
                if isinstance(candidate, Mapping)
                else getattr(candidate, "shared_parameters", None)
            )
            if shared is None and isinstance(candidate, Mapping):
                shared = candidate
            if not isinstance(shared, Mapping):
                continue
            for names, expected in (
                (("logCstar", "logC"), self.logC),
                (("logPstar", "logP"), self.logPstar),
            ):
                for name in names:
                    if name in shared:
                        actual = _mp(shared[name])
                        if actual != expected:
                            raise ValueError(
                                f"shared {name} disagrees with centered defect input"
                            )
                        break

    def _jet(self, value: Any = 0) -> PressureWidthJet:
        if isinstance(value, PressureWidthJet):
            if value.orders != (self.pressure_order, self.width_order):
                raise ValueError(
                    "switch component jet orders do not match centered defect orders"
                )
            return value
        return PressureWidthJet(
            value,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def _dual(self, value: Any = 0, tangent: Any = 0) -> AxialDual:
        return AxialDual(
            self._jet(value),
            self._jet(tangent),
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def _restore_constants(self) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        """Resolve the three scalar restoration integrals by composite Gauss."""

        nodes, weights = mp.gauss_quadrature(self.restore_order, "legendre")
        pieces = ((mp.mpf(0), mp.mpf(".5")), (mp.mpf(".5"), mp.mpf(1)))

        def integrate(function: Any) -> mp.mpf:
            total = mp.mpf(0)
            for left, right in pieces:
                midpoint = (left + right) / 2
                half = (right - left) / 2
                for node, weight in zip(nodes, weights):
                    t = midpoint + half * node
                    total += half * weight * function(t)
            return total

        K1 = integrate(lambda t: mp.exp(t) * (1 - _sigma_mp(t)))
        K16 = integrate(lambda t: mp.exp(mp.mpf("1.6") * t) * (1 - _sigma_mp(t)))
        K2 = integrate(lambda t: mp.exp(t) * (1 - _sigma_mp(t)) ** 2)
        return K1, K16, K2

    def _source(self, Z: Any) -> dict[str, Any]:
        """Read and dualize the exact switch endpoint, including P0."""

        z = _mp(Z)
        key = mp.nstr(z, self.work_precision)
        cached = self._z_cache.get(key)
        if cached is not None:
            return cached
        source = self.switches.evaluate_R(mp.mpf(110), z)
        for name in ("P0", "P0_Z"):
            if name not in source:
                raise KeyError(
                    "switch endpoint must expose explicit P0 and P0_Z; "
                    "total pressure cannot substitute for the axis datum"
                )
        F110 = self._dual(source["F"], _field_tangent(source, "F"))
        Uz110 = self._dual(source["Uz"], _field_tangent(source, "Uz"))
        moments: dict[str, AxialDual] = {}
        source_moments = source.get("moments")
        source_moments_Z = source.get("moments_Z", source.get("momentsZ"))
        if not isinstance(source_moments, Mapping) or not isinstance(
            source_moments_Z, Mapping
        ):
            raise KeyError("switch endpoint must expose moments and moments_Z")
        for name in MOMENT_KEYS:
            if name not in source_moments or name not in source_moments_Z:
                raise KeyError(f"switch endpoint is missing moment {name!r}")
            moments[name] = self._dual(source_moments[name], source_moments_Z[name])
        P0 = self._dual(source["P0"], source["P0_Z"])
        result = {
            "Z": z,
            "switch": source,
            "F110": F110,
            "Uz110": Uz110,
            "moments": moments,
            "P0": P0,
            "P110": self._dual(source.get("P", 0), _field_tangent(source, "P")),
        }
        self._z_cache[key] = result
        return result

    @staticmethod
    def _add_row(rows: dict[int, AxialDual], row: int, term: AxialDual) -> None:
        if row in rows:
            rows[row] = rows[row] + term
        else:
            rows[row] = term

    @staticmethod
    def _part_value(parts: Mapping[str, Mapping[int, AxialDual]]) -> dict[str, dict[int, PressureWidthJet]]:
        return {
            label: {row: value.value for row, value in rows.items()}
            for label, rows in parts.items()
        }

    @staticmethod
    def _part_tangent(parts: Mapping[str, Mapping[int, AxialDual]]) -> dict[str, dict[int, PressureWidthJet]]:
        return {
            label: {row: value.tangent for row, value in rows.items()}
            for label, rows in parts.items()
        }

    def _kernel_records(
        self,
        B: AxialDual,
    ) -> dict[str, dict[str, Any]]:
        records: dict[str, dict[str, Any]] = {}
        for name, (k_text, m_text) in KERNEL_SPECS.items():
            k, m = _mp(k_text), _mp(m_text)
            records[name] = compose_flat_shape_defect(
                k,
                m,
                B.value,
                B.tangent,
                self.T,
                precision=self.precision,
                order=self.order,
                window=self.window,
            )
        return records

    @staticmethod
    def _flat_term(record: Mapping[str, Any], n: int) -> AxialDual:
        return AxialDual(record["value_terms"][n], record["tangent_terms"][n])

    def evaluate(self, Z: Any) -> dict[str, Any]:
        """Evaluate all five centered defect rows at one axial coordinate."""

        with mp.workdps(self.work_precision):
            source = self._source(Z)
            zdual = self._dual(source["Z"], 1)
            F110 = source["F110"]
            Uz110 = source["Uz110"]
            moments = source["moments"]
            P0 = source["P0"]
            Am = mp.exp(self.logPstar - mp.mpf(".6")) / (1 + zdual * zdual)
            uref110 = mp.exp(-self.logC) / (1 + zdual * zdual)
            g = Uz110 - 4 * zdual
            B = (
                (mp.sqrt(220) * F110).log()
                + self.logC
                + (1 + zdual * zdual).log()
            )
            flat_records = self._kernel_records(B)
            S = mp.sqrt(2) * self.Rm ** mp.mpf("1.5") * Am

            parts: dict[str, dict[int, AxialDual]] = {}

            def add(label: str, row: int, value: AxialDual) -> None:
                parts.setdefault(label, {})[row] = value

            # Centered seed rows at the switch endpoint.  The labels retain
            # every subtractive constant as an independently inspectable term.
            add("centered_d1_Mz_over_Rm", 1, moments["z"] / self.Rm)
            add("centered_d1_sub_440Z_over_Rm", 1, -440 * zdual / self.Rm)
            add("centered_d2_Mtheta_z_over_S", 2, moments["theta_z"] / S)
            add(
                "centered_d2_sub_4Z_Mtheta_over_S",
                2,
                -4 * zdual * moments["theta"] / S,
            )
            add("centered_d3_Mtheta_over_S", 3, moments["theta"] / S)
            add(
                "centered_d3_sub_5_8_sqrt2_11032_uref_over_S",
                3,
                -mp.mpf(5) / 8 * mp.sqrt(2) * mp.mpf(110) ** mp.mpf("1.5") * uref110 / S,
            )
            add(
                "centered_d4_Mztheta_over_RmAm2",
                4,
                moments["z_theta"] / (self.Rm * Am * Am),
            )
            add(
                "centered_d4_sub_8Z_Mz_over_RmAm2",
                4,
                -8 * zdual * moments["z"] / (self.Rm * Am * Am),
            )
            add(
                "centered_d4_16Z2_110_over_RmAm2",
                4,
                16 * zdual * zdual * 110 / (self.Rm * Am * Am),
            )
            add(
                "centered_d4_5_12_110_uref2_over_RmAm2",
                4,
                mp.mpf(5) / 12 * 110 * uref110 * uref110 / (self.Rm * Am * Am),
            )
            add("centered_d5_Mp_over_Am2", 5, moments["p"] / (Am * Am))
            add(
                "centered_d5_sub_2p5_uref2_over_Am2",
                5,
                -mp.mpf("2.5") * uref110 * uref110 / (Am * Am),
            )

            # Direct source terms from 110 to the reference scale.  Constants
            # are split from the flat kernels so cancellation remains visible.
            add("source_d1_g_exp_minus2", 1, g * mp.exp(-2))
            add("source_d1_sub_g_rho0", 1, -g * self.rho0)
            add("source_d2_g_exp_minus3p2_over_1p6", 2, g * mp.exp(-mp.mpf("3.2")) / mp.mpf("1.6"))
            add("source_d2_sub_g_rho0p1p6_over_1p6", 2, -g * self.rho0 ** mp.mpf("1.6") / mp.mpf("1.6"))
            add("source_d4_g2_exp_minus2_over_Am2", 4, g * g / (Am * Am) * mp.exp(-2))
            add("source_d4_sub_g2_rho0_over_Am2", 4, -g * g / (Am * Am) * self.rho0)

            # Keep each composition Taylor term as a separate component row.
            theta_record = flat_records["theta"]
            energy_record = flat_records["energy"]
            pressure_record = flat_records["pressure"]
            nmax = self.pressure_order + self.width_order
            for n in range(nmax + 1):
                theta_term = self._flat_term(theta_record, n)
                energy_term = self._flat_term(energy_record, n)
                pressure_term = self._flat_term(pressure_record, n)
                add(
                    f"source_d2_g_alpha1p6_theta_n{n}",
                    2,
                    g * (self.alpha ** mp.mpf("1.6")) * theta_term,
                )
                add(
                    f"source_d3_alpha1p6_theta_n{n}",
                    3,
                    (self.alpha ** mp.mpf("1.6")) * theta_term,
                )
                add(
                    f"source_d4_sub_half_alpha1p2_energy_n{n}",
                    4,
                    -mp.mpf(".5") * (self.alpha ** mp.mpf("1.2")) * energy_term,
                )
                add(
                    f"source_d5_half_alpha0p2_pressure_n{n}",
                    5,
                    mp.mpf(".5") * (self.alpha ** mp.mpf(".2")) * pressure_term,
                )

            # Scalar restoration over the first unit phase.  No increments
            # are added after t = 1; the labels state that support explicitly.
            add("restore_d1_g_exp_minus2_K1", 1, g * mp.exp(-2) * self.K1)
            add("restore_d2_g_exp_minus3p2_K16", 2, g * mp.exp(-mp.mpf("3.2")) * self.K16)
            add("restore_d4_g2_exp_minus2_K2_over_Am2", 4, g * g / (Am * Am) * mp.exp(-2) * self.K2)

            rows: dict[int, AxialDual] = {
                row: self._dual(0) for row in range(1, 6)
            }
            for row_parts in parts.values():
                for row, value in row_parts.items():
                    self._add_row(rows, row, value)

            d = {row: value.value for row, value in rows.items()}
            d_Z = {row: value.tangent for row, value in rows.items()}
            metadata = self.metadata()
            metadata.update(
                {
                    "Z": mp.nstr(source["Z"], 40),
                    "work_precision": self.work_precision,
                    "flat_kernel_order": self.order,
                    "flat_kernel_window": mp.nstr(self.window, 40),
                    "restore_order": self.restore_order,
                    "Rm": mp.nstr(self.Rm, 40),
                    "Rsh": mp.nstr(self.Rsh, 40),
                    "Rz": mp.nstr(self.Rz, 40),
                    "rho0": mp.nstr(self.rho0, 40),
                    "alpha": mp.nstr(self.alpha, 40),
                    "K1": mp.nstr(self.K1, 40),
                    "K16": mp.nstr(self.K16, 40),
                    "K2": mp.nstr(self.K2, 40),
                }
            )
            return {
                "Z": source["Z"],
                "d": d,
                "d_Z": d_Z,
                "d_dual": rows,
                "parts": self._part_value(parts),
                "parts_dual": parts,
                "parts_Z": self._part_tangent(parts),
                "parts_value": self._part_value(parts),
                "P0": P0.value,
                "P0_Z": P0.tangent,
                "P0_dual": P0,
                "P0_receipt": {
                    "P0": P0.value,
                    "P0_Z": P0.tangent,
                    "source": "explicit switch axis datum",
                    "total_minus_moment_substitution": False,
                },
                "P110": source["P110"].value,
                "P110_Z": source["P110"].tangent,
                "F110": F110.value,
                "F110_Z": F110.tangent,
                "Uz110": Uz110.value,
                "Uz110_Z": Uz110.tangent,
                "B": B.value,
                "B_Z": B.tangent,
                "Am": Am.value,
                "Am_Z": Am.tangent,
                "uref110": uref110.value,
                "uref110_Z": uref110.tangent,
                "g": g.value,
                "g_Z": g.tangent,
                "flat_records": flat_records,
                "metadata": metadata,
            }

    __call__ = evaluate

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "independent_variable": self.independent_variable,
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "automatic_first_Z_derivative": True,
            "explicit_P0_preserved": True,
            "complete_formal_five_defect_assembly": True,
            "functional_closure": False,
            "smallness_certified": False,
            "error_enclosures": False,
            "quadrature_remainder_enclosed": False,
            "pressure_width_remainder_enclosed": False,
            "finite_Z_difference_used": False,
            "global_field_installed": False,
            "no_increments_after_restoration_t1": True,
            "pressure_parameter_distinct_from_temporal": True,
            "Rsh_less_than_Rz_guard": "T < 10*(logC + logPstar) - 8",
            "shared_log_parameter_consistency_checked": True,
        }


__all__ = ["CenteredComponentDefects"]
