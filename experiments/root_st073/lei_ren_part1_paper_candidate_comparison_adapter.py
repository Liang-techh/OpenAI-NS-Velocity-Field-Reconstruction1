"""Nominal MP adapter for the directed Lambda120 comparison jets.

``CandidateComparisonJets`` carries interval Taylor coefficients at the
single accepted axial center.  This adapter projects those coefficients to
nominal MP scalars for the existing transition classes while retaining the
candidate amplitude and pressure primitives algebraically.  It is a center-
only numerical adapter: it does not enlarge the axial domain, enclose RK4
discretization, or certify the transition.
"""

from __future__ import annotations

from typing import Any

import mpmath as mp

from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    """Return directed endpoints for an MP interval or a scalar."""

    raw = getattr(value, "_mpi_", None)
    if raw is None:
        scalar = mp.mpf(str(value))
        return scalar, scalar
    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _midpoint(value: Any) -> mp.mpf:
    left, right = _endpoints(value)
    return (left + right) / 2


class CandidateMPComparison:
    """Project one ``CandidateComparisonJets`` center to nominal MP fields.

    The constructor intentionally accepts the already-built ``calc`` and
    ``core`` objects.  It does not import or depend on the candidate exit
    module, avoiding a cycle with the transition builder.
    """

    branch = "candidate_lambda120_comparison_nominal_center"

    def __init__(self, calc: Any, core: Any) -> None:
        for name in ("evaluate", "rhs", "ctx", "precision", "center", "lam", "eps", "delta", "hb"):
            if not hasattr(calc, name):
                raise TypeError(f"calc must expose {name}")
        if not hasattr(core, "evaluate"):
            raise TypeError("core must expose evaluate(R, Z)")
        self.calc = calc
        self.core = core
        self.precision = int(calc.precision)
        with mp.workdps(self.precision + 40):
            self.center = _midpoint(calc.center)
            self.Lambda = _midpoint(calc.lam)
            self.epsilon = _midpoint(calc.eps)
            self.delta = _midpoint(calc.delta)
            self.h_b = mp.mpf(str(calc.hb))
            self.hb = self.h_b
            self.R_a = mp.mpf(4) / self.Lambda
            self.F0 = _midpoint(calc.F0)
        if self.Lambda <= 0:
            raise ValueError("candidate Lambda must be positive")
        if not 0 <= self.delta < 1:
            raise ValueError("candidate delta must satisfy 0 <= delta < 1")
        if self.h_b <= 0:
            raise ValueError("candidate h_b must be positive")
        self._cache: dict[tuple[str, str], dict[str, Any]] = {}
        self._amplitude = self._amplitude_jet()

    def _amplitude_jet(self) -> IntervalTaylor:
        """Recover F0(Z) through the retained axial order from ell=F0'/F0."""

        ctx = self.calc.ctx
        coefficients: list[Any] = [self.calc.F0]
        for order in range(1, 4):
            total = ctx.mpf(0)
            for index in range(order):
                total += self.calc.ell[index] * coefficients[order - 1 - index]
            coefficients.append(total / ctx.mpf(order))
        return IntervalTaylor(ctx, coefficients)

    def _alpha(self, y: mp.mpf) -> mp.mpf:
        if y <= self.h_b:
            return mp.mpf(1)
        if y >= 2 * self.h_b:
            return mp.mpf(0)
        x = (y - self.h_b) / self.h_b
        left = mp.exp(-1 / (x * x))
        right = mp.exp(-1 / ((1 - x) * (1 - x)))
        return right / (left + right)

    def _derivative_jets(
        self, y: mp.mpf, state: dict[str, IntervalTaylor]
    ) -> tuple[IntervalTaylor, IntervalTaylor]:
        """Return d(F)/dy and d(Uz)/dy as axial jets."""

        zero = state["phi"] * 0
        if y > 2 * self.h_b:
            return zero, state["U"] * 0
        derivative = self.calc.rhs(y, state)
        return self._amplitude * derivative["phi"], derivative["U"]

    def _physical_jets(
        self, state: dict[str, IntervalTaylor]
    ) -> tuple[IntervalTaylor, IntervalTaylor, dict[str, IntervalTaylor], IntervalTaylor, dict[str, IntervalTaylor]]:
        """Apply the physical epsilon/F0 scales without refitting any field."""

        eps_scalar = self.epsilon
        # S is the accepted F0^2 axial jet.  Keep it as supplied instead of
        # replacing it by a square of projected midpoint values.
        S = self.calc.S
        amplitude = self._amplitude
        F = amplitude * state["phi"]
        U = state["U"]
        moments = {
            "theta": amplitude * state["theta"] * (eps_scalar ** 2),
            "z": state["z"] * eps_scalar,
            "theta_z": amplitude * state["theta_z"] * (eps_scalar ** 2),
            "z_theta": state["u_squared"] * eps_scalar
            - S * state["weighted_phi_squared"] * (eps_scalar ** 2),
            "p": S * state["p"] * eps_scalar,
        }
        pressure = self.calc.p0 + moments["p"]
        return F, U, moments, pressure, {
            "raw_axial": state["u_squared"] * eps_scalar,
            "raw_swirl": S * state["weighted_phi_squared"] * (eps_scalar ** 2),
        }

    def evaluate(self, y: Any, Z: Any = ".3") -> dict[str, Any]:
        """Return nominal MP comparison fields at the accepted center only."""

        with mp.workdps(self.precision + 40):
            yy = mp.mpf(str(y))
            zz = mp.mpf(str(Z))
            if yy < 0:
                raise ValueError("comparison requires y >= 0")
            if abs(zz - self.center) > mp.mpf(10) ** (-self.precision + 5):
                raise ValueError("candidate comparison supplies only its accepted axial center")
            key = (mp.nstr(yy, self.precision), mp.nstr(zz, self.precision))
            cached = self._cache.get(key)
            if cached is not None:
                return dict(cached)

            packet = self.calc.evaluate(yy, zz)
            state = packet["state"]
            state_jets = {name: value for name, value in state.items()}
            F_jet, U_jet, moments_jet, pressure_jet, raw_jet = self._physical_jets(state_jets)
            F_y_jet, U_y_jet = self._derivative_jets(yy, state_jets)
            scaled_radius = _midpoint(packet["scaled_radius"])
            R = self.epsilon * scaled_radius
            root = mp.sqrt(2 * R)

            F = _midpoint(F_jet[0])
            F_Z = _midpoint(F_jet[1])
            Uz = _midpoint(U_jet[0])
            Uz_Z = _midpoint(U_jet[1])
            if F == 0:
                raise ArithmeticError("candidate comparison amplitude vanished")
            F_y = _midpoint(F_y_jet[0])
            F_y_Z = _midpoint(F_y_jet[1])
            Uz_y = _midpoint(U_y_jet[0])
            Uz_y_Z = _midpoint(U_y_jet[1])
            F_R = F_y / R
            F_RZ = F_y_Z / R
            Uz_R = Uz_y / R
            Uz_RZ = Uz_y_Z / R
            moments = {name: _midpoint(value[0]) for name, value in moments_jet.items()}
            moments_Z = {name: _midpoint(value[1]) for name, value in moments_jet.items()}
            P = _midpoint(pressure_jet[0])
            P_Z = _midpoint(pressure_jet[1])
            ur = (
                2 * zz * R * Uz
                - (1 - self.delta) * zz * moments["z"]
                - (1 - zz * zz) * moments_Z["z"]
            ) / ((1 - self.delta * zz * zz) * root)

            stress = evaluate_mp_stress(
                mp.log(R),
                zz,
                self.delta,
                Utheta=root * F,
                Uz=Uz,
                Utheta_y=root * (F_y + F / 2),
                Utheta_Z=root * F_Z,
                Uz_y=Uz_y,
                Uz_Z=Uz_Z,
                moments=moments,
                moments_Z=moments_Z,
                P=P,
                P_Z=P_Z,
                precision=self.precision + 20,
            )

            D = _midpoint(packet["D"][0])
            D_Z = _midpoint(packet["D"][1])
            I_z = _midpoint(packet["I_z"][0])
            I_z_Z = _midpoint(packet["I_z"][1])
            result = {
                "y": yy,
                "Z": zz,
                "R": R,
                "alpha": self._alpha(yy),
                "F": F,
                "F_Z": F_Z,
                "FZ": F_Z,
                "Uz": Uz,
                "Uz_Z": Uz_Z,
                "UZ": Uz_Z,
                "F_R": F_R,
                "FR": F_R,
                "F_RZ": F_RZ,
                "Uz_R": Uz_R,
                "UR": Uz_R,
                "Uz_RZ": Uz_RZ,
                "moments": moments,
                "moments_Z": moments_Z,
                "momentsZ": moments_Z,
                "P": P,
                "P_Z": P_Z,
                "PZ": P_Z,
                "P_R": F * F,
                "Ur": ur,
                "D": D,
                "D_Z": D_Z,
                "I_theta": D * F,
                "I_z": I_z,
                "I_z_Z": I_z_Z,
                "E": I_z / F,
                "F_R_over_F": F_R / F,
                "S_theta": stress["S_theta"],
                "S_z": stress["S_z"],
                "T_theta": stress["T_theta"],
                "T_z": stress["T_z"],
                "stress": stress,
                "raw_quadratic_integrals": {
                    "axial": _midpoint(raw_jet["raw_axial"][0]),
                    "axial_Z": _midpoint(raw_jet["raw_axial"][1]),
                    "swirl": _midpoint(raw_jet["raw_swirl"][0]),
                    "swirl_Z": _midpoint(raw_jet["raw_swirl"][1]),
                },
                "core_scope": "candidate comparison center nominal projection",
                "core_reference_R": R,
                "Z_derivative_source": "IntervalTaylor coefficient 1 projected to MP midpoint",
                "finite_Z_difference_used": False,
                "amplitude_factor_cancelled_algebraically": True,
                "pressure_source": "accepted p0 + S * epsilon * normalized p",
                "source_constants_certified": False,
                "cone_certified": False,
                "outer_matching_complete": False,
            }
            self._cache[key] = result
            return dict(result)

    def _transition_state(self, y: Any, Z: Any) -> dict[str, Any]:
        """Return the state keys consumed by Section923 frozen algebra."""

        value = self.evaluate(y, Z)
        return {
            "F": value["F"],
            "F_Z": value["F_Z"],
            "Uz": value["Uz"],
            "Uz_Z": value["Uz_Z"],
            "P": value["P"],
            "P_Z": value["P_Z"],
            "moments": dict(value["moments"]),
            "moments_Z": dict(value["moments_Z"]),
        }

    def frozen_coefficients(self, Z: Any) -> dict[str, Any]:
        """Reuse the established Section923 frozen power algebra."""

        from lei_ren_part1_paper_exit_comparison import Section923Comparison

        return Section923Comparison.frozen_coefficients(self, Z)

    def frozen_driver_coefficients(self, Z: Any) -> dict[str, Any]:
        """Return the unshifted frozen driver coefficients."""

        from lei_ren_part1_paper_exit_comparison import Section923Comparison

        return Section923Comparison.frozen_driver_coefficients(self, Z)


__all__ = ["CandidateMPComparison", "MOMENT_KEYS"]
