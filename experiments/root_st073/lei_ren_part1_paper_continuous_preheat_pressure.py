"""Analytic pre-heat axis pressure for the continuous source schedule.

The source pressure datum in Lei--Ren (6.10) is

``P0_pre(Z) = -integral_0^infinity Utheta_pre(R,Z)^2/(2 R) dR``.

The pre-heat profile replaces the heat factor by ``H_delta = 1``.  Before
``Rv`` the angular factor is exactly ``(1 + Z**2)**-1``.  The interval
``[Rv, Rf]`` is different: the source uses the multiplicative flattening

``Utheta/A = 2**(theta - 1) * (1 + Z**2)**(-theta),
theta = 1 - sigma((y-yv)/Tf)``.

This module keeps that interval, every post-``Rv`` radial atom, and the
``H=1`` heat-collar/exterior atoms explicit.  It never folds tiny post-``Rv``
terms into the dominant pre-``Rv`` coefficient.  The returned pressure is
normalised by ``Pstar**2``; multiply by ``Pstar**2`` only at the API boundary
that requests physical units.

Variable radial stages are evaluated with mpmath Gauss--Legendre nodes.
Declared pure-power stages use exact exponential atoms.  Taylor coefficients
are coefficients in ``(Z-center)**n``.  For ``center=0`` the flatten factors
use the exact binomial coefficients ``(-1)**k (beta)_k/k!``, with
``beta=2*theta``.  At a general center the equivalent stable recurrence is
used, so no float-Z schedule evaluation is needed.

This is a pressure-datum adapter, not a field installer or a claim of full
five-moment, core, cone, PDE, or theorem compatibility.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_continuous_axial_pulse import (  # noqa: E402
    ContinuousAxialPulse,
)
from lei_ren_part1_paper_continuous_axis_pressure import (  # noqa: E402
    ContinuousAxisPressureJets,
    REFERENCE_EXTENSION_INTEGRAL,
)


SOURCE = "https://arxiv.org/html/2609.35406v2"
SOURCE_VERSION = "2609.35406v2"


def _mp(value: Any) -> mp.mpf:
    """Parse a value without routing it through binary64."""

    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _n(value: Any, digits: int) -> str:
    return mp.nstr(_mp(value), int(digits))


def _finite_z(value: Any, *, name: str = "Z") -> mp.mpf:
    try:
        result = _mp(value)
    except Exception as exc:  # pragma: no cover - defensive API boundary
        raise ValueError(f"{name} must be a finite real number in [-1, 1]") from exc
    if not mp.isfinite(result) or abs(result) > 1:
        raise ValueError(f"{name} must be a finite real number in [-1, 1]")
    return result


def _q_power_taylor(beta: mp.mpf, center: mp.mpf, degree: int) -> list[mp.mpf]:
    """Taylor coefficients of ``(1+Z**2)**(-beta)`` about ``center``.

    The centered recurrence follows from

    ``(1+Z**2) g'(Z) + 2 beta Z g(Z) = 0``.

    At zero the binomial form is preferable because all odd coefficients are
    exactly zero and the rising factorial is evaluated directly.
    """

    beta = _mp(beta)
    center = _mp(center)
    degree = int(degree)
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    Q = 1 + center * center
    c0 = mp.exp(-beta * mp.log(Q))
    if degree == 0:
        return [c0]
    coefficients = [mp.mpf(0)] * (degree + 1)
    coefficients[0] = c0
    if center == 0:
        for k in range(1, degree // 2 + 1):
            coefficients[2 * k] = (
                (-1) ** k * mp.rf(beta, k) / mp.factorial(k) * c0
            )
        return coefficients

    # n=0 with c[-1]=0.
    coefficients[1] = -2 * beta * center * coefficients[0] / Q
    for n in range(1, degree):
        lower = coefficients[n - 1]
        middle = coefficients[n]
        coefficients[n + 1] = -(
            2 * center * (n + beta) * middle
            + (n - 1 + 2 * beta) * lower
        ) / (Q * (n + 1))
    return coefficients


def _flat_edge(value: mp.mpf) -> mp.mpf:
    """The paper's flat edge ``f(t)=exp(-1/t**2)`` for ``t>0``."""

    value = _mp(value)
    if value <= 0:
        return mp.mpf(0)
    return mp.exp(-1 / (value * value))


class ContinuousPreheatPressure:
    """Standalone ``P0_pre`` datum from one shared continuous schedule.

    ``source`` can be an existing :class:`ContinuousAxisPressureJets` prefix
    adapter or a source profile accepted by that adapter.  Passing an
    existing prefix object is preferred because it preserves its shared
    schedule and its cached MP prefix receipt.
    """

    _VARIABLE_POST_RV = (
        "z_flatten",
        "steep_transition_in",
        "steep_transition_out",
    )
    _CONSTANT_POST_RV = (
        "power_buffer_rel",
        "steep_power",
        "waiting",
    )

    def __init__(
        self,
        source: Any,
        *,
        quadrature_order: int = 192,
    ) -> None:
        if int(quadrature_order) != quadrature_order or quadrature_order < 16:
            raise ValueError("quadrature_order must be an integer at least 16")
        if isinstance(source, ContinuousAxisPressureJets):
            prefix = source
        else:
            if not hasattr(source, "schedule"):
                raise TypeError(
                    "source must be ContinuousAxisPressureJets or a profile "
                    "with a shared schedule"
                )
            prefix = ContinuousAxisPressureJets(
                source, quadrature_order=int(quadrature_order)
            )
        schedule = prefix.schedule
        if getattr(schedule, "_continuous_angular_provider", None) is None:
            raise ValueError(
                "ContinuousPreheatPressure requires the shared continuous "
                "angular schedule provider"
            )
        self.prefix = prefix
        self.profile = prefix.profile
        self.schedule = schedule
        self.precision = int(prefix.precision)
        self.quadrature_order = int(quadrature_order)
        with mp.workdps(self.precision):
            self.mu = _mp(schedule.mu)
            self.delta = _mp(schedule.delta)
            self.epsilon = _mp(schedule.epsilon)
            self.log_pstar = _mp(schedule.logPstar)
            self.y_v = _mp(schedule.y_v)
            self.y_f = _mp(schedule.y_f)
            self.y_rel = _mp(schedule.y_rel)
            self.y_s = _mp(schedule.y_s)
            self.y_q = _mp(schedule.y_q)
            self.y_t = _mp(schedule.y_t)
            self.y_tail = _mp(schedule.y_tail)
            self.y_b = _mp(schedule.y_b)
            self.log_r_tail = _mp(schedule.logR_tail)
            self.log_c_inf = _mp(schedule._log_c_inf)
        self._component_cache: dict[int, dict[str, Any]] = {}

    def _mp_nodes(self, order: int) -> tuple[tuple[mp.mpf, mp.mpf], ...]:
        shared = getattr(self.prefix, "_mp_gauss_nodes", None)
        if shared is not None:
            return shared(int(order))
        with mp.workdps(self.precision):
            nodes, weights = mp.gauss_quadrature(int(order), "legendre")
            return tuple(
                ((mp.mpf(node) + 1) / 2, mp.mpf(weight) / 2)
                for node, weight in zip(nodes, weights)
            )

    def _bounds(self, name: str) -> tuple[mp.mpf, mp.mpf]:
        bounds = self.schedule._make_stage_bounds()
        left, right = bounds[name]
        if right is None:
            raise ValueError(f"stage {name} has no finite endpoint")
        return _mp(left), _mp(right)

    def _log_amplitude_ratio(self, y: mp.mpf) -> mp.mpf:
        """Return ``log(A(y)/Pstar)`` using the shared continuous primitive."""

        value = self.schedule._log_A(_mp(y))
        return _mp(value) - self.log_pstar

    @staticmethod
    def _atom(rate: mp.mpf, length: mp.mpf) -> mp.mpf:
        """Return ``integral_0^length exp(rate*s) ds`` stably."""

        rate = _mp(rate)
        length = _mp(length)
        if length <= 0:
            return mp.mpf(0)
        if rate == 0:
            return length
        return mp.expm1(rate * length) / rate

    @staticmethod
    def _heat_collar_factor(t: mp.mpf, epsilon: mp.mpf) -> mp.mpf:
        """The pre-heat terminal bracket ``K0`` with ``H_delta`` set to 1."""

        sigma, _ = ContinuousAxialPulse.sigma_pair(_mp(t))
        edge = (3 - _mp(t)) / 2
        flat = _flat_edge(edge)
        terminal = 1 - epsilon * flat
        return (1 - sigma) * (1 - epsilon) + sigma * terminal

    def _integrate_variable_scalar(
        self,
        left: mp.mpf,
        right: mp.mpf,
        function,
        order: int,
    ) -> mp.mpf:
        length = _mp(right) - _mp(left)
        if length <= 0:
            return mp.mpf(0)
        total = mp.fsum(
            weight * function(left + length * node)
            for node, weight in self._mp_nodes(order)
        )
        return length * total

    def _constant_power_value(
        self,
        name: str,
        left: mp.mpf,
        right: mp.mpf,
        slope: mp.mpf,
    ) -> tuple[mp.mpf, dict[str, Any]]:
        """Exact atom for a post-flatten stage with ``U=A/2``."""

        length = right - left
        log_e = self._log_amplitude_ratio(left)
        # U/Pstar = A/(2 Pstar), and pressure density is U^2/2.
        base = mp.exp(2 * log_e) / 8
        atom = self._atom(2 * slope, length)
        value = base * atom
        return value, {
            "stage": name,
            "kind": "exact_exponential_atom",
            "left": left,
            "right": right,
            "slope": slope,
            "base": base,
            "atom": atom,
            "value_at_Z0": value,
            "Z_independent": True,
        }

    def _pre_rv_components(self, order: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
        receipt = self.prefix.preflatten_integral(quadrature_order=order)
        components = []
        for row in receipt.get("intervals", ()):
            components.append(
                {
                    "stage": row["stage"],
                    "region": "pre_Rv",
                    "kind": "q_power",
                    "beta": mp.mpf(2),
                    "value_at_Z0": _mp(row["integral"]),
                    "source_receipt": row,
                }
            )
        reference = {
            "stage": "reference_extension",
            "region": "reference",
            "kind": "q_power",
            "beta": mp.mpf(2),
            "value_at_Z0": _mp(REFERENCE_EXTENSION_INTEGRAL),
            "source_receipt": {
                "formula": "5/2 from the reference branch on [0,Rref]"
            },
        }
        components.insert(0, reference)
        return components, receipt

    def _compute_components(self, order: int) -> dict[str, Any]:
        order = int(order)
        cached = self._component_cache.get(order)
        if cached is not None:
            return cached
        with mp.workdps(self.precision):
            prefix_components, prefix_receipt = self._pre_rv_components(order)
            # The flattening source factor is U/Pstar =
            # (A/Pstar) * 2**(theta-1) * q**(-theta), so the pressure
            # density has the weight (A/Pstar)^2 * 2**(2 theta-3).
            left, right = self._bounds("z_flatten")
            length = right - left
            flatten_atoms: list[dict[str, Any]] = []
            for node, weight in self._mp_nodes(order):
                y = left + length * node
                local = (y - self.y_v) / (self.y_f - self.y_v)
                sigma, _ = ContinuousAxialPulse.sigma_pair(local)
                theta = 1 - sigma
                beta = 2 * theta
                base = mp.exp(2 * self._log_amplitude_ratio(y)) * mp.power(
                    2, beta - 3
                )
                atom = length * weight * base
                flatten_atoms.append(
                    {
                        "node": node,
                        "weight": weight,
                        "y": y,
                        "theta": theta,
                        "beta": beta,
                        "base": base,
                        "atom_at_Z0": atom,
                    }
                )
            flatten_value = mp.fsum(row["atom_at_Z0"] for row in flatten_atoms)
            post_components: list[dict[str, Any]] = [
                {
                    "stage": "z_flatten",
                    "region": "post_Rv",
                    "kind": "q_power_atoms",
                    "left": left,
                    "right": right,
                    "value_at_Z0": flatten_value,
                    "atoms": flatten_atoms,
                    "Z_independent": False,
                    "formula": (
                        "integral (A/Pstar)^2 * 2^(2 theta-3) * "
                        "(1+Z^2)^(-2 theta) dy"
                    ),
                }
            ]

            # The two steep transitions have no Z factor after flattening;
            # retain MP quadrature because A(y) changes by sigma there.
            for name in ("steep_transition_in", "steep_transition_out"):
                stage_left, stage_right = self._bounds(name)
                value = self._integrate_variable_scalar(
                    stage_left,
                    stage_right,
                    lambda y: mp.exp(2 * self._log_amplitude_ratio(y)) / 8,
                    order,
                )
                post_components.append(
                    {
                        "stage": name,
                        "region": "post_Rv",
                        "kind": "MP_Gauss_constant_Z",
                        "left": stage_left,
                        "right": stage_right,
                        "value_at_Z0": value,
                        "Z_independent": True,
                    }
                )

            slopes = {
                "power_buffer_rel": -mp.mpf("0.5") - self.mu,
                "steep_power": -mp.mpf("1.5"),
                "waiting": -(1 + self.delta) / 2,
            }
            for name in self._CONSTANT_POST_RV:
                stage_left, stage_right = self._bounds(name)
                value, data = self._constant_power_value(
                    name, stage_left, stage_right, slopes[name]
                )
                data.update({"region": "post_Rv"})
                post_components.append(data)

            # In (6.4), replacing H_delta by 1 leaves the Z-independent
            # terminal bracket K0(t).  The collar is finite (t in [0,3]);
            # the power tail beyond Rb is an exact atom.
            heat_left, heat_right = self._bounds("heat_connection")
            heat_length = heat_right - heat_left
            lam = 1 + self.delta
            heat_scale = mp.exp(
                2 * (self.log_c_inf - self.log_pstar)
                - mp.log(2)
                - lam * self.log_r_tail
            )
            collar_integral = self._integrate_variable_scalar(
                mp.mpf(0), heat_length,
                lambda t: mp.exp(-lam * t)
                * self._heat_collar_factor(t, self.epsilon) ** 2,
                order,
            )
            collar_value = heat_scale * collar_integral
            post_components.append(
                {
                    "stage": "heat_collar",
                    "region": "post_Rv",
                    "kind": "MP_Gauss_K0_H1",
                    "left": heat_left,
                    "right": heat_right,
                    "value_at_Z0": collar_value,
                    "heat_scale": heat_scale,
                    "collar_integral": collar_integral,
                    "lambda": lam,
                    "Z_independent": True,
                    "H_replaced_by_one": True,
                    "formula": "heat_scale * integral_0^3 exp(-(1+delta)t) K0(t)^2 dt",
                }
            )
            exterior_value = heat_scale * mp.exp(-lam * heat_length) / lam
            post_components.append(
                {
                    "stage": "exterior_power_tail",
                    "region": "post_Rv",
                    "kind": "exact_exterior_power_atom",
                    "left": self.y_b,
                    "right": mp.inf,
                    "value_at_Z0": exterior_value,
                    "heat_scale": heat_scale,
                    "lambda": lam,
                    "tail_integral": mp.exp(-lam * heat_length) / lam,
                    "Z_independent": True,
                    "H_replaced_by_one": True,
                    "formula": "heat_scale * exp(-3*(1+delta))/(1+delta)",
                }
            )
            prefix_value = mp.fsum(row["value_at_Z0"] for row in prefix_components)
            post_value = mp.fsum(row["value_at_Z0"] for row in post_components)
            result = {
                "order": order,
                "prefix_components": prefix_components,
                "post_components": post_components,
                "prefix_value_at_Z0": prefix_value,
                "post_value_at_Z0": post_value,
                "total_value_at_Z0": prefix_value + post_value,
                "prefix_receipt": prefix_receipt,
                "normalization": "positive integral I=Pstar^-2 * (-P0_pre)",
                "dominant_K": -prefix_value,
                "post_Rv_preserved_separately": True,
                "flatten_theta": "theta=1-sigma((y-y_v)/T_f)",
                "heat_collar_K0": "(1-sigma(t))(1-epsilon)+sigma(t)(1-epsilon*f((3-t)/2))",
                "heat_tail_H": "H_delta replaced by 1",
                "quadrature_error_enclosed": False,
            }
        self._component_cache[order] = result
        return result

    @staticmethod
    def _zero_tail(degree: int) -> list[mp.mpf]:
        return [mp.mpf(0)] * (int(degree) + 1)

    def _component_integral_coefficients(
        self,
        component: dict[str, Any],
        *,
        center: mp.mpf,
        degree: int,
    ) -> list[mp.mpf]:
        kind = component["kind"]
        if kind == "q_power":
            return [
                component["value_at_Z0"] * value
                for value in _q_power_taylor(
                    component["beta"], center, degree
                )
            ]
        if kind == "q_power_atoms":
            result = self._zero_tail(degree)
            for atom in component["atoms"]:
                q_coefficients = _q_power_taylor(
                    atom["beta"], center, degree
                )
                for index, value in enumerate(q_coefficients):
                    result[index] += atom["atom_at_Z0"] * value
            return result
        return [component["value_at_Z0"]] + self._zero_tail(degree)[1:]

    def taylor_components(
        self,
        degree: int = 8,
        *,
        center: Any = 0,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return per-stage pressure Taylor components.

        Coefficients are for ``P0_pre(center+h)`` as powers of ``h``.  The
        positive integral coefficients are returned alongside them because
        the sign is easy to lose when comparing pressure and moment APIs.
        ``dominant`` contains the reference plus pre-``Rv`` prefix only;
        ``post_Rv`` contains every retained atom after ``Rv``.
        """

        if int(degree) != degree or degree < 0:
            raise ValueError("degree must be a nonnegative integer")
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        degree = int(degree)
        with mp.workdps(self.precision):
            center_mp = _finite_z(center, name="center")
            components_data = self._compute_components(order)
            all_components = (
                components_data["prefix_components"]
                + components_data["post_components"]
            )
            component_rows: dict[str, dict[str, Any]] = {}
            total = self._zero_tail(degree)
            dominant = self._zero_tail(degree)
            post = self._zero_tail(degree)
            prefix_names = {
                row["stage"] for row in components_data["prefix_components"]
            }
            post_names = {
                row["stage"] for row in components_data["post_components"]
            }
            for component in all_components:
                integral_coefficients = self._component_integral_coefficients(
                    component, center=center_mp, degree=degree
                )
                pressure_coefficients = [-value for value in integral_coefficients]
                for index, value in enumerate(integral_coefficients):
                    total[index] += value
                    if component["stage"] in prefix_names:
                        dominant[index] += value
                    elif component["stage"] in post_names:
                        post[index] += value
                row = dict(component)
                row.update(
                    {
                        "center": center_mp,
                        "degree": degree,
                        "integral_coefficients": integral_coefficients,
                        "pressure_coefficients": pressure_coefficients,
                        "integral_value": integral_coefficients[0],
                        "pressure_value": pressure_coefficients[0],
                        "derivative": (
                            pressure_coefficients[1]
                            if degree >= 1
                            else None
                        ),
                    }
                )
                component_rows[component["stage"]] = row
            pressure = [-value for value in total]
            dominant_pressure = [-value for value in dominant]
            post_pressure = [-value for value in post]
            return {
                "source": SOURCE,
                "source_version": SOURCE_VERSION,
                "center": center_mp,
                "degree": degree,
                "quadrature_order": order,
                "pressure_units": "P0_over_Pstar_squared",
                "Pstar_squared": mp.exp(2 * self.log_pstar),
                "coefficients": pressure,
                "pressure_coefficients": pressure,
                "integral_coefficients": total,
                "dominant_pressure_coefficients": dominant_pressure,
                "post_Rv_pressure_coefficients": post_pressure,
                "dominant_K": -components_data["prefix_value_at_Z0"],
                "post_Rv_value_at_center": post_pressure[0],
                "pressure_value": pressure[0],
                "pressure_Z": pressure[1] if degree >= 1 else None,
                "components": component_rows,
                "post_Rv_components": {
                    name: component_rows[name]
                    for name in sorted(post_names)
                },
                "dominant_components": {
                    name: component_rows[name]
                    for name in sorted(prefix_names)
                },
                "post_Rv_preserved_separately": True,
                "summed_coefficients_are_nominal": True,
                "consume_separate_atoms_for_tiny_tail": True,
                "Taylor_convention": "coefficient of (Z-center)^n; derivative=n!*coefficient",
                "q_power_jet_formula": (
                    "at center=0: c[2k]=(-1)^k*(beta)_k/k!, c[2k+1]=0; "
                    "general center uses Q(n+1)c[n+1]+2 center(n+beta)c[n]"
                    "+(n-1+2 beta)c[n-1]=0"
                ),
                "quadrature_error_enclosed": False,
            }

    def post_Rv_components_jet(
        self,
        degree: int = 8,
        *,
        center: Any = 0,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return only post-``Rv`` components and their pressure jet."""

        result = self.taylor_components(
            degree, center=center, quadrature_order=quadrature_order
        )
        return {
            "source": result["source"],
            "source_version": result["source_version"],
            "center": result["center"],
            "degree": result["degree"],
            "quadrature_order": result["quadrature_order"],
            "pressure_units": result["pressure_units"],
            "Pstar_squared": result["Pstar_squared"],
            "pressure_coefficients": result["post_Rv_pressure_coefficients"],
            "integral_coefficients": [
                -value for value in result["post_Rv_pressure_coefficients"]
            ],
            "pressure_value": result["post_Rv_pressure_coefficients"][0],
            "pressure_Z": (
                result["post_Rv_pressure_coefficients"][1]
                if degree >= 1
                else None
            ),
            "components": result["post_Rv_components"],
            "post_Rv_preserved_separately": True,
            "quadrature_error_enclosed": result["quadrature_error_enclosed"],
        }

    def pressure_jet(
        self,
        Z: Any,
        degree: int = 8,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return ``P0_pre`` and its Taylor jet centred at ``Z``."""

        return self.taylor_components(
            degree,
            center=Z,
            quadrature_order=quadrature_order,
        )

    def pressure_value(
        self, Z: Any, *, quadrature_order: int | None = None
    ) -> mp.mpf:
        return self.pressure_jet(Z, degree=0, quadrature_order=quadrature_order)[
            "pressure_value"
        ]

    def pressure_Z(
        self, Z: Any, *, quadrature_order: int | None = None
    ) -> mp.mpf:
        return self.pressure_jet(Z, degree=1, quadrature_order=quadrature_order)[
            "pressure_Z"
        ]

    def taylor(
        self,
        degree: int = 8,
        *,
        center: Any = 0,
        quadrature_order: int | None = None,
    ) -> list[mp.mpf]:
        """Short alias returning only pressure Taylor coefficients."""

        return self.taylor_components(
            degree, center=center, quadrature_order=quadrature_order
        )["pressure_coefficients"]

    # Names used by older pressure adapters and by small integration probes.
    value = pressure_value
    derivative = pressure_Z
    preheat_datum = pressure_value

    def receipt(
        self,
        Z: Any = 0,
        *,
        degree: int = 8,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return a compact receipt while retaining raw MP values in fields."""

        row = self.pressure_jet(
            Z, degree=degree, quadrature_order=quadrature_order
        )
        digits = self.precision
        return {
            "source": row["source"],
            "source_version": row["source_version"],
            "center": _n(row["center"], digits),
            "degree": row["degree"],
            "quadrature_order": row["quadrature_order"],
            "pressure_units": row["pressure_units"],
            "Pstar_squared": _n(row["Pstar_squared"], digits),
            "pressure_value": _n(row["pressure_value"], digits),
            "pressure_Z": (
                _n(row["pressure_Z"], digits)
                if row["pressure_Z"] is not None
                else None
            ),
            "pressure_coefficients": [
                _n(value, digits) for value in row["pressure_coefficients"]
            ],
            "dominant_pressure_coefficients": [
                _n(value, digits)
                for value in row["dominant_pressure_coefficients"]
            ],
            "post_Rv_pressure_coefficients": [
                _n(value, digits)
                for value in row["post_Rv_pressure_coefficients"]
            ],
            "dominant_K": _n(row["dominant_K"], digits),
            "post_Rv_value_at_center": _n(row["post_Rv_value_at_center"], digits),
            "component_values": {
                name: _n(component["pressure_value"], digits)
                for name, component in row["components"].items()
            },
            "post_Rv_preserved_separately": True,
            "quadrature_error_enclosed": False,
        }


# Descriptive aliases make the adapter easy to find without changing the
# canonical class name used by the parent integration lane.
ContinuousPreheatPressureDatum = ContinuousPreheatPressure
PreheatPressureDatum = ContinuousPreheatPressure


__all__ = [
    "ContinuousPreheatPressure",
    "ContinuousPreheatPressureDatum",
    "PreheatPressureDatum",
]
