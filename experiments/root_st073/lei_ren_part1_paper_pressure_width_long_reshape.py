"""Finite pressure/width-ring long reshape after the ``R = 110`` switches.

This is the Section 9.30 shape candidate on the same finite pressure/width
and first-axial-derivative ring as the preceding adapters.  The endpoint
moments and raw quadratic integrals are retained as separate seeds; the
backward normalized quadrature contributes a separate increment.  The public
sum is useful for the stress formulas, but does not erase the two receipts.
No source, quadrature, jet, cone, or global-field error is enclosed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


@lru_cache(maxsize=16)
def _nodes(order: int, precision: int):
    """Gauss nodes cached at the precision used by the reshape."""

    with mp.workdps(int(precision)):
        return mp.gauss_quadrature(int(order), "legendre")


def _split_dual(value: Any, tangent: bool = False) -> Any:
    if isinstance(value, AxialDual):
        return value.tangent if tangent else value.value
    return value


class PressureWidthLongReshape:
    """Section 9.30 long reshape over the pressure/width AxialDual ring."""

    branch = "pressure_width_long_reshape"

    def __init__(
        self,
        switches: Any,
        *,
        A: Any = "1e150",
        logC: Any = "5e151",
        order: int = 24,
        tail_digits: int = 60,
    ) -> None:
        if int(order) != order or order < 8:
            raise ValueError("order must be an integer at least 8")
        if int(tail_digits) != tail_digits or tail_digits < 1:
            raise ValueError("tail_digits must be a positive integer")
        if not hasattr(switches, "evaluate_R") or not hasattr(switches, "comparison"):
            raise TypeError("switches must expose evaluate_R and comparison")
        comparison = switches.comparison
        if not getattr(comparison, "automatic_Z_tangent", False):
            raise TypeError("PressureWidthLongReshape requires an axial comparison")
        self.switches = switches
        self.comparison = comparison
        self.precision = int(switches.precision)
        self.order = int(order)
        self.tail_digits = int(tail_digits)
        self.R0 = mp.mpf("110")
        with mp.workdps(self.precision):
            self.A = mp.mpf(str(A))
            self.logC = mp.mpf(str(logC))
            if self.A <= 0:
                raise ValueError("A must be positive")
            self.T = mp.mpf(400) * self.A
            self.work_precision = self.precision + max(
                0, int(mp.floor(mp.log10(self.T))) + 1
            ) + 30
            self.cutoff = (
                (self.tail_digits + 20)
                * mp.log(10)
                / mp.mpf(".1")
            )
            if self.T <= 0:
                raise ValueError("T must be positive")

            # Match the shared source input when the bundle carries it.  A
            # compact resolved fixture may omit that metadata, in which case
            # the caller's explicit logC is the complete input.
            bundle = getattr(comparison, "bundle", None)
            shared = bundle.get("shared_parameters", {}) if isinstance(bundle, Mapping) else {}
            expected_logC = shared.get("logCstar", shared.get("logC"))
            if expected_logC is not None and self.logC != mp.mpf(str(expected_logC)):
                raise ValueError("reshape logC must match the inherited shared candidate")

        self._start_cache: dict[str, dict[str, Any]] = {}

    def _dual(self, value: Any = 0, tangent: Any = 0) -> AxialDual:
        return self.comparison._dual(value, tangent)

    def _base_jet(self, value: Any = 0) -> PressureWidthJet:
        if hasattr(self.comparison, "_base_jet"):
            return self.comparison._base_jet(value)
        return PressureWidthJet(
            value,
            pressure_order=self.comparison.pressure_order,
            width_order=self.comparison.width_order,
        )

    def _z_key(self, Z: Any) -> str:
        return mp.nstr(mp.mpf(str(Z)), self.precision)

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
        tangent_source = source.get("moments_Z", source.get("momentsZ", {}))
        return {
            name: self._dual(source["moments"][name], tangent_source[name])
            for name in MOMENT_KEYS
        }

    def _raw_duals(self, source: Mapping[str, Any]) -> dict[str, AxialDual]:
        values = source["raw_quadratic_integrals"]
        tangents = source.get("raw_quadratic_integrals_Z", {})
        result: dict[str, AxialDual] = {}
        for name in ("axial", "swirl"):
            value = values[name]
            if isinstance(value, AxialDual):
                result[name] = value
                continue
            tangent = values.get(f"{name}_Z", tangents.get(name, tangents.get(f"{name}_Z", 0)))
            result[name] = self._dual(value, tangent)
        return result

    def _start(self, Z: Any) -> dict[str, Any]:
        """Cache the exact switch endpoint and its separate seed receipts."""

        with mp.workdps(self.work_precision):
            z = mp.mpf(str(Z))
            key = self._z_key(z)
            cached = self._start_cache.get(key)
            if cached is not None:
                return cached
            endpoint = self.switches.evaluate_R(self.R0, z)
            F = self._field_dual(endpoint, "F")
            V = self._field_dual(endpoint, "Uz")
            P = self._field_dual(endpoint, "P")
            axis_P0 = self._field_dual(endpoint, "P0")
            moments = self._moment_duals(endpoint)
            raw = self._raw_duals(endpoint)
            root220 = mp.sqrt(220)
            u1 = root220 * F
            logu1 = u1.log()
            zdual = self._dual(z, 1)
            B = logu1 + self.logC + (1 + zdual * zdual).log()

            # Check the nominal source input only.  The finite ring and the
            # first axial tangent remain separate; no C2 or uniform bound is
            # inferred from this pointwise budget check.
            B_nominal = B.value.evaluate(pressure=1, width=1)
            BZ_nominal = B.tangent.evaluate(pressure=1, width=1)
            budget_ok = max(abs(B_nominal), abs(BZ_nominal)) <= 2 * self.A
            if not budget_ok:
                raise ArithmeticError("nominal reshape B/BZ exceeds the 2A source budget")

            result = dict(
                endpoint=endpoint,
                F=F,
                V=V,
                P=P,
                axis_P0=axis_P0,
                moments=moments,
                raw=raw,
                u1=u1,
                logu1=logu1,
                B=B,
                zdual=zdual,
                B_nominal=B_nominal,
                BZ_nominal=BZ_nominal,
                budget_ok=budget_ok,
                seed_radius=self._base_jet(self.R0),
            )
            self._start_cache[key] = result
            return result

    @staticmethod
    def _sigma_prime(s: mp.mpf, sigma: mp.mpf) -> mp.mpf:
        if 0 < s < 1:
            return sigma * (1 - sigma) * (
                2 / s**3 + 2 / (1 - s) ** 3
            )
        return mp.mpf(0)

    def _shape(self, y: mp.mpf, start: Mapping[str, Any]) -> dict[str, Any]:
        with mp.workdps(self.work_precision):
            s = y / self.T
            sigma = _sigma_mp(s)
            sigma_prime = self._sigma_prime(s, sigma)
            logu = start["logu1"] + y / 10 - sigma * start["B"]
            u = logu.exp()
            q = mp.mpf(".1") - sigma_prime * start["B"] / self.T
            a = mp.mpf(".8") + 2 * sigma_prime * start["B"] / self.T
            R = mp.mpf(110) * mp.exp(y)
            F = u / mp.sqrt(2 * R)
            g_y = -a / 2
            return dict(
                s=s,
                sigma=sigma,
                sigma_prime=sigma_prime,
                logu=logu,
                u=u,
                q=q,
                a=a,
                R=R,
                F=F,
                g_y=g_y,
            )

    def _quadrature(
        self,
        y: mp.mpf,
        shape: Mapping[str, Any],
        start: Mapping[str, Any],
    ) -> tuple[dict[str, AxialDual], dict[str, mp.mpf]]:
        """Integrate endpoint-normalized theta/pressure/swirl increments."""

        with mp.workdps(self.work_precision):
            end = min(y, self.cutoff)
            totals = {
                "theta": self._dual(0),
                "pressure": self._dual(0),
                "swirl": self._dual(0),
            }
            nodes, weights = _nodes(self.order, self.work_precision)
            if end > 0:
                cuts = [mp.mpf(0)]
                next_cut = mp.mpf(1)
                while next_cut < end:
                    cuts.append(next_cut)
                    next_cut *= 2
                cuts.append(end)
                for left, right in zip(cuts, cuts[1:]):
                    midpoint = (left + right) / 2
                    half_width = (right - left) / 2
                    for node, weight in zip(nodes, weights):
                        ell = midpoint + half_width * node
                        sigma_prev = _sigma_mp((y - ell) / self.T)
                        ratio = (
                            -ell / 10
                            - (sigma_prev - shape["sigma"]) * start["B"]
                        ).exp()
                        angular = mp.exp(-mp.mpf("1.5") * ell) * ratio
                        pressure = ratio * ratio
                        swirl = mp.exp(-ell) * pressure
                        factor = half_width * weight
                        totals["theta"] += factor * angular
                        totals["pressure"] += factor * pressure
                        totals["swirl"] += factor * swirl

            omitted = y > end
            bounds = {
                "theta": mp.exp(-mp.mpf("1.55") * end) / mp.mpf("1.55")
                if omitted
                else mp.mpf(0),
                "pressure": mp.exp(-mp.mpf(".1") * end) / mp.mpf(".1")
                if omitted
                else mp.mpf(0),
                "swirl": mp.exp(-mp.mpf("1.1") * end) / mp.mpf("1.1")
                if omitted
                else mp.mpf(0),
            }
            return totals, bounds

    def _stress(
        self,
        z: mp.mpf,
        shape: Mapping[str, Any],
        V: AxialDual,
        moments: Mapping[str, AxialDual],
        P: AxialDual,
    ) -> dict[str, Any]:
        """Evaluate stress from the same public fields and moment receipts."""

        R = self._base_jet(shape["R"])
        zbase = self._base_jet(z)
        u = shape["u"]
        theta_y = u * shape["q"]
        values = {name: item.value for name, item in moments.items()}
        tangents = {name: item.tangent for name, item in moments.items()}
        return evaluate_mp_stress(
            mp.log(shape["R"]),
            z,
            self.comparison.delta,
            Utheta=u.value,
            Uz=V.value,
            Utheta_y=theta_y.value,
            Utheta_Z=u.tangent,
            Uz_y=self._base_jet(0),
            Uz_Z=V.tangent,
            moments=values,
            moments_Z=tangents,
            P=P.value,
            P_Z=P.tangent,
            precision=self.work_precision,
            radius_override=R,
            scalar_converter=self._base_jet,
            axial_override=zbase,
            shear_theta=2 * shape["F"].value * shape["g_y"].value,
            shear_z=self._base_jet(0),
            include_components=True,
        )

    @staticmethod
    def _public_dict(duals: Mapping[str, AxialDual]) -> tuple[dict[str, Any], dict[str, Any]]:
        return (
            {key: value.value for key, value in duals.items()},
            {key: value.tangent for key, value in duals.items()},
        )

    def evaluate_log_offset(self, y: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the reshape for ``0 <= y = log(R/110) <= T``."""

        with mp.workdps(self.work_precision):
            y_value = mp.mpf(str(y))
            z = mp.mpf(str(Z))
            if not 0 <= y_value <= self.T or abs(z) >= 1:
                raise ValueError("Require 0 <= y <= T and |Z| < 1")
            start = self._start(z)
            shape = self._shape(y_value, start)
            increments, tail_bounds = self._quadrature(y_value, shape, start)
            R_delta = self._dual(shape["R"] - self.R0, 0)

            # Keep the seed and every increment as separate duals before
            # forming the public totals.
            seed_moments = dict(start["moments"])
            increment_moments = {
                "theta": mp.sqrt(2) * shape["R"] ** mp.mpf("1.5") * shape["u"] * increments["theta"],
                "z": start["V"] * R_delta,
                "theta_z": start["V"] * (
                    mp.sqrt(2) * shape["R"] ** mp.mpf("1.5") * shape["u"] * increments["theta"]
                ),
                "z_theta": self._dual(0),
                "p": shape["u"] * shape["u"] / 2 * increments["pressure"],
            }
            axial_increment = start["V"] * start["V"] * R_delta
            swirl_increment = shape["R"] * shape["u"] * shape["u"] / 2 * increments["swirl"]
            increment_moments["z_theta"] = axial_increment - swirl_increment
            moments = {
                key: seed_moments[key] + increment_moments[key]
                for key in MOMENT_KEYS
            }
            P_increment = increment_moments["p"]
            P = start["P"] + P_increment
            raw_seed = dict(start["raw"])
            raw_increment = {"axial": axial_increment, "swirl": swirl_increment}
            raw = {
                key: raw_seed[key] + raw_increment[key]
                for key in ("axial", "swirl")
            }
            # At the inherited endpoint return the exact switch receipt.  A
            # finite reciprocal followed by its inverse can otherwise leave
            # a discarded high-order ring atom in F or U^r even though the
            # analytic increment is identically zero at y = 0.
            if y_value == 0:
                endpoint = start["endpoint"]
                F = start["F"]
                moments = dict(seed_moments)
                P = start["P"]
                raw = dict(raw_seed)
                stress = endpoint["stress"]
                F_R = self._field_dual(endpoint, "F_R") if "F_R" in endpoint else F * shape["g_y"] / shape["R"]
                Uz_R = self._field_dual(endpoint, "Uz_R") if "Uz_R" in endpoint else self._dual(0)
                u_public = start["u1"]
            else:
                stress = self._stress(z, shape, start["V"], moments, P)
                F = shape["F"]
                F_R = F * shape["g_y"] / shape["R"]
                Uz_R = self._dual(0)
                u_public = shape["u"]
            moments_public, moments_Z_public = self._public_dict(moments)
            seed_public, seed_Z_public = self._public_dict(seed_moments)
            increment_public, increment_Z_public = self._public_dict(increment_moments)
            raw_public, raw_Z_public = self._public_dict(raw)
            raw_seed_public, raw_seed_Z_public = self._public_dict(raw_seed)
            raw_increment_public, raw_increment_Z_public = self._public_dict(raw_increment)
            normalized_public, normalized_Z_public = self._public_dict(increments)
            # Keep the historical receipt spelling ``p`` while retaining the
            # descriptive ``pressure`` spelling used internally.
            normalized_public["p"] = normalized_public["pressure"]
            normalized_Z_public["p"] = normalized_Z_public["pressure"]
            tail_bounds["p"] = tail_bounds["pressure"]

            result = {
                "R": self._base_jet(shape["R"]),
                "R_Z": self._base_jet(0),
                "logR": mp.log(110) + y_value,
                "log_radius_offset": y_value,
                "Z": z,
                "F": F.value,
                "FZ": F.tangent,
                "F_Z": F.tangent,
                "Utheta": u_public.value,
                "Utheta_Z": u_public.tangent,
                "Uz": start["V"].value,
                "UZ": start["V"].tangent,
                "Uz_Z": start["V"].tangent,
                "P": P.value,
                "PZ": P.tangent,
                "P_Z": P.tangent,
                "P0": start["axis_P0"].value,
                "P0_Z": start["axis_P0"].tangent,
                "Ur": stress["U_r"],
                "moments": moments_public,
                "momentsZ": moments_Z_public,
                "moments_Z": moments_Z_public,
                "stress": stress,
                "a": shape["a"].value,
                "a_Z": shape["a"].tangent,
                "b": self._base_jet(0),
                "g_y": shape["g_y"].value,
                "g_y_Z": shape["g_y"].tangent,
                "u_y": self._base_jet(0),
                "Uz_y": self._base_jet(0),
                "F_R": F_R.value,
                "F_RZ": F_R.tangent,
                "F_R_Z": F_R.tangent,
                "Uz_R": Uz_R.value,
                "Uz_RZ": Uz_R.tangent,
                "Uz_R_Z": Uz_R.tangent,
                "moment_seeds": seed_public,
                "moment_seeds_Z": seed_Z_public,
                "moment_increments": increment_public,
                "moment_increments_Z": increment_Z_public,
                "seed_moments": seed_public,
                "seed_moments_Z": seed_Z_public,
                "increment_moments": increment_public,
                "increment_moments_Z": increment_Z_public,
                "raw_quadratic_integrals": raw_public,
                "raw_quadratic_integrals_Z": raw_Z_public,
                "raw_quadratic_seeds": raw_seed_public,
                "raw_quadratic_seeds_Z": raw_seed_Z_public,
                "raw_quadratic_increments": raw_increment_public,
                "raw_quadratic_increments_Z": raw_increment_Z_public,
                "P_seed": start["P"].value,
                "P_seed_Z": start["P"].tangent,
                "P_increment": P_increment.value,
                "P_increment_Z": P_increment.tangent,
                "axis_P0": start["axis_P0"].value,
                "axis_P0_Z": start["axis_P0"].tangent,
                "logu1": start["logu1"].value,
                "logu1_Z": start["logu1"].tangent,
                "B": start["B"].value,
                "B_Z": start["B"].tangent,
                "B_nominal": start["B_nominal"],
                "BZ_nominal": start["BZ_nominal"],
                "source_budget_ok": start["budget_ok"],
                "normalized_integrals": normalized_public,
                "normalized_integrals_Z": normalized_Z_public,
                "normalized_integrals_scope": "increments from R=110, endpoint-normalized backward ell quadrature",
                "normalized_tail_bounds": tail_bounds,
                "normalized_Z_tail_bounds": {
                    key: value
                    * max(mp.mpf(1), abs(start["BZ_nominal"]))
                    * (1 if key == "theta" else 2)
                    for key, value in tail_bounds.items()
                },
                "quadrature_order": self.order,
                "quadrature_error_certified": False,
                "tail_assumptions_locally_checked": True,
                "normalized_tail_envelopes_nominal_conditional": True,
                "tail_bounds_scope": "nominal conditional .05 growth bound; not atomwise enclosure",
                "source_constants_certified": False,
                "cone_certified": False,
                "outer_matching_complete": False,
                "functional_terminal_moments_closed": False,
                "region": "long_swirl_reshape",
                "metadata": self.metadata(),
            }
            return result

    def evaluate_phase(self, phase: Any, Z: Any) -> dict[str, Any]:
        """Evaluate at a normalized shape phase in ``[0, 1]``."""

        with mp.workdps(self.work_precision):
            phase_value = mp.mpf(str(phase))
            if not 0 <= phase_value <= 1:
                raise ValueError("phase must lie in [0, 1]")
            return self.evaluate_log_offset(phase_value * self.T, Z)

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "source_equation": "Section 9.30 long reshape",
            "reference_radius": "110",
            "shape_interval": "0 <= log(R/110) <= T",
            "A": mp.nstr(self.A, 40),
            "T": mp.nstr(self.T, 40),
            "logC": mp.nstr(self.logC, 40),
            "pressure_order": int(self.comparison.pressure_order),
            "width_order": int(self.comparison.width_order),
            "quadrature_order": self.order,
            "tail_digits": self.tail_digits,
            "automatic_Z_tangent": True,
            "endpoint_moment_seeds_preserved": True,
            "endpoint_raw_quadratic_seeds_preserved": True,
            "nominal_B_budget_checked": True,
            "nominal_C2_uniform_bounds_proved": False,
            "conditional_normalized_tail_envelopes": True,
            "tail_growth_bound": ".05",
            "normalized_tail_bounds_atomwise_enclosed": False,
            "quadrature_error_certified": False,
            "pressure_width_jet_remainder_enclosed": False,
            "source_error_enclosed": False,
            "cone_certified": False,
            "global_field_installed": False,
            "outer_matching_complete": False,
            "functional_terminal_moments_closed": False,
        }


__all__ = ["PressureWidthLongReshape"]
