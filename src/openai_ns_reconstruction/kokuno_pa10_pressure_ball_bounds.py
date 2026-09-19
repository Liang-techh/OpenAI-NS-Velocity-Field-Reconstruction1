"""Conditional coefficient-space bounds for Kokuno's PA.10 pressure map.

Pinned public provenance is ``KokunoYumeto/yang-mills-interacting-workbench`` at
commit ``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction defines the coefficient norm

    a_{alpha,beta} = 20^{-alpha} rho^{-beta} beta!
                     binom(alpha+beta,beta)
                     / ((alpha+1)^2 (beta+1)^2)

and, after the leading rescaling,

    p = I(g^2 Phi^2),       I F(Y) = integral_0^Y F(v) dv,
    Y p_Y = Y g^2 Phi^2,
    partial_eta p = partial_eta I(g^2 Phi^2).

The source explicitly states that multiplication, ``I``, multiplication by
``Y`` and ``partial_eta I`` are bounded in this coefficient space. Its displayed
weight ratios give the conservative executable factors

    ||I F||_rho <= 80 ||F||_rho,
    ||Y F||_rho <= 80 ||F||_rho,
    ||partial_eta I F||_rho <= (80/rho) ||F||_rho.

Together with the source product algebra

    ||F G||_rho <= C_sq^2 ||F||_rho ||G||_rho,
    C_sq = 4 pi^2 / 3,

this module turns the pressure map into a fail-closed radius-ball calculator.
For the numerical envelope we use the elementary rational inequality
``pi < 22/7`` and exact ``Fraction.from_float`` arithmetic before directed
binary64 rounding. Thus every returned pressure bound is outward-rounded with
respect to the supplied binary64 inputs rather than depending on accidental
round-to-nearest direction.

The calculator is conditional. A caller must still supply a source-valid
``rho``, a coefficient-space norm for ``g`` and a radius-ball norm/Lipschitz
bound for ``Phi``. Merely passing diagnostic numbers does not make the source
pressure ball bound true. In particular this module does not promote source
R1/R2, M/K, B0/T_sh, PA.16, global matched pressure, global leading velocity,
or PDE validation.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field, replace
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_operator_primitives import KokunoPA10OperatorPrimitiveBounds
from .kokuno_pa10_remainder_ball_bounds import (
    BallFactorBound,
    KokunoPA10RemainderBallBounds,
    KokunoPA10RemainderPrimitiveBounds,
)


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-pressure-ball-bounds-v1"

_SOURCE_FORMULAS = {
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha} rho^{-beta} beta! "
        "binom(alpha+beta,beta)/((alpha+1)^2(beta+1)^2)"
    ),
    "pressure": "p=I(g^2 Phi^2)",
    "pressure_radial": "Y p_Y=Y g^2 Phi^2",
    "pressure_eta": "partial_eta p=partial_eta I(g^2 Phi^2)",
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "radial_integral": "||I F||_rho <= 80 ||F||_rho",
    "multiply_Y": "||Y F||_rho <= 80 ||F||_rho",
    "eta_radial_integral": "||partial_eta I F||_rho <= (80/rho)||F||_rho",
    "lipschitz": "replace moving Phi factors one at a time",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_pressure_formula_recorded": True,
    "source_pressure_operator_algebra_executable": True,
    "source_I_bound_derived_from_displayed_weight_ratio": True,
    "source_Y_bound_derived_from_displayed_weight_ratio": True,
    "source_deta_I_bound_derived_from_displayed_weight_ratio": True,
    "conditional_pressure_ball_norm_executable": True,
    "conditional_pressure_ball_lipschitz_executable": True,
    "conditional_pressure_to_R1_R2_bridge_executable": True,
    "binary64_bound_arithmetic_outward_rounded": True,
    "diagnostic_inputs_are_source_bounds": False,
    "source_rho_machine_bound": False,
    "source_g_coefficient_norm_machine_bound": False,
    "source_Phi_radius_one_ball_norm_machine_bound": False,
    "source_Phi_radius_one_ball_lipschitz_machine_bound": False,
    "source_pressure_radius_one_ball_norm_machine_bound": False,
    "source_pressure_radius_one_ball_lipschitz_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_R1_radius_one_ball_lipschitz_machine_bound": False,
    "source_R2_radius_one_ball_lipschitz_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _nonnegative_finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and >=0")
    return out


def _positive_finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and >0")
    return out


def _fraction_upper(value: Fraction) -> float:
    """Smallest convenient binary64 upper bound for a nonnegative Fraction."""

    if value < 0:
        raise ValueError("fraction upper-bound helper requires a nonnegative value")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper-bound arithmetic left binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return float(out)


def _as_fraction(value: float) -> Fraction:
    out = _nonnegative_finite(value, "bound")
    return Fraction.from_float(out)


@dataclass(frozen=True)
class KokunoPA10PressureBallBounds:
    """Source-derived conditional bounds for ``p=I(g^2 Phi^2)``."""

    operators: KokunoPA10OperatorPrimitiveBounds = field(
        default_factory=KokunoPA10OperatorPrimitiveBounds
    )

    def __post_init__(self) -> None:
        if not isinstance(self.operators, KokunoPA10OperatorPrimitiveBounds):
            raise TypeError("operators must be KokunoPA10OperatorPrimitiveBounds")

    @property
    def product_constant_upper(self) -> float:
        """Rationally certified upper bound for ``C_sq^2`` using ``pi<22/7``."""

        pi_upper = Fraction(22, 7)
        C_sq_upper = Fraction(4, 3) * pi_upper * pi_upper
        return _fraction_upper(C_sq_upper * C_sq_upper)

    @staticmethod
    def radial_integral_factor() -> float:
        """Conservative factor for ``I`` from the displayed weight ratio."""

        return 80.0

    @staticmethod
    def multiply_Y_factor() -> float:
        """Conservative factor for multiplication by ``Y``."""

        return 80.0

    @staticmethod
    def eta_radial_integral_factor(rho: float) -> float:
        """Conservative factor for ``partial_eta I``.

        The source's second weight ratio contributes ``(80/rho)(i+1)`` and
        the radial integration divisor ``i+1`` cancels that final factor.
        The supplied binary64 ``rho`` is treated as its exact rational value.
        """

        r = _positive_finite(rho, "rho")
        exact = Fraction(80, 1) / Fraction.from_float(r)
        return _fraction_upper(exact)

    def product(self, *factors: BallFactorBound) -> BallFactorBound:
        """Iterate product algebra with exact-rational outward arithmetic."""

        if not factors:
            return BallFactorBound(1.0, 0.0)
        if not all(isinstance(item, BallFactorBound) for item in factors):
            raise TypeError("all factors must be BallFactorBound values")

        P = Fraction.from_float(self.product_constant_upper)
        algebra = P ** max(len(factors) - 1, 0)
        norms = [_as_fraction(item.norm) for item in factors]
        lips = [_as_fraction(item.lipschitz) for item in factors]
        norm_exact = algebra * math.prod(norms, start=Fraction(1, 1))

        lip_exact = Fraction(0, 1)
        for index, lip in enumerate(lips):
            if lip == 0:
                continue
            other = math.prod(
                (norm for j, norm in enumerate(norms) if j != index),
                start=Fraction(1, 1),
            )
            lip_exact += lip * other
        lip_exact *= algebra
        return BallFactorBound(_fraction_upper(norm_exact), _fraction_upper(lip_exact))

    @staticmethod
    def scale(factor: BallFactorBound, scalar: float) -> BallFactorBound:
        if not isinstance(factor, BallFactorBound):
            raise TypeError("factor must be a BallFactorBound")
        s = _nonnegative_finite(abs(float(scalar)), "scalar")
        sf = Fraction.from_float(s)
        return BallFactorBound(
            _fraction_upper(sf * _as_fraction(factor.norm)),
            _fraction_upper(sf * _as_fraction(factor.lipschitz)),
        )

    def pressure_bounds(
        self,
        *,
        rho: float,
        g_norm: float,
        Phi: BallFactorBound,
    ) -> dict[str, BallFactorBound]:
        """Return conditional bounds for ``p``, ``p_eta`` and ``Y p_Y``.

        ``g`` is fixed by the selected/source rescaling, so its contraction-ball
        Lipschitz value is zero. ``Phi.lipschitz`` is supplied by the caller's
        ball metric. The two copies of ``Phi`` are replaced one at a time,
        exactly as in the source's local Lipschitz bookkeeping.
        """

        r = _positive_finite(rho, "rho")
        g = BallFactorBound(_nonnegative_finite(g_norm, "g_norm"), 0.0)
        if not isinstance(Phi, BallFactorBound):
            raise TypeError("Phi must be a BallFactorBound")

        integrand = self.product(g, g, Phi, Phi)
        p = self.scale(integrand, self.radial_integral_factor())
        p_eta = self.scale(integrand, self.eta_radial_integral_factor(r))
        Y_p_Y = self.scale(integrand, self.multiply_Y_factor())
        return {
            "integrand_g2_Phi2": integrand,
            "p": p,
            "p_eta": p_eta,
            "Y_p_Y": Y_p_Y,
        }

    def inject_pressure_bounds(
        self,
        primitives: KokunoPA10RemainderPrimitiveBounds,
        pressure: dict[str, BallFactorBound],
    ) -> KokunoPA10RemainderPrimitiveBounds:
        """Replace only the three pressure factors consumed by PA.10 ``R2``."""

        if not isinstance(primitives, KokunoPA10RemainderPrimitiveBounds):
            raise TypeError("primitives must be KokunoPA10RemainderPrimitiveBounds")
        required = {"p", "p_eta", "Y_p_Y"}
        if not required.issubset(pressure):
            raise ValueError("pressure mapping must contain p, p_eta and Y_p_Y")
        if not all(isinstance(pressure[name], BallFactorBound) for name in required):
            raise TypeError("pressure entries must be BallFactorBound values")
        return replace(
            primitives,
            p=pressure["p"],
            p_eta=pressure["p_eta"],
            Y_p_Y=pressure["Y_p_Y"],
        )

    def conditional_remainder_bridge(
        self,
        *,
        rho: float,
        g_norm: float,
        Phi: BallFactorBound,
        primitives: KokunoPA10RemainderPrimitiveBounds,
        remainder: KokunoPA10RemainderBallBounds,
    ) -> dict[str, Any]:
        """Feed conditional pressure bounds into the existing exact R1/R2 algebra."""

        if not isinstance(remainder, KokunoPA10RemainderBallBounds):
            raise TypeError("remainder must be KokunoPA10RemainderBallBounds")
        pressure = self.pressure_bounds(rho=rho, g_norm=g_norm, Phi=Phi)
        patched = self.inject_pressure_bounds(primitives, pressure)
        return {
            "pressure": pressure,
            "remainder": remainder.remainder_envelopes(patched),
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        """Deterministic arithmetic receipt using explicitly diagnostic inputs."""

        rho = 5.0e-4
        g_norm = 1.25
        Phi = BallFactorBound(norm=2.5, lipschitz=1.0)
        pressure = self.pressure_bounds(rho=rho, g_norm=g_norm, Phi=Phi)

        base = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
        calculator = KokunoPA10RemainderBallBounds(
            h=0.01,
            A=0.51,
            D=0.49,
            rescaling_lambda=2.0,
        )
        bridge = self.conditional_remainder_bridge(
            rho=rho,
            g_norm=g_norm,
            Phi=Phi,
            primitives=base,
            remainder=calculator,
        )

        pressure_json = {
            name: {"norm": bound.norm, "lipschitz": bound.lipschitz}
            for name, bound in pressure.items()
        }
        remainder_totals = {
            "R1": copy.deepcopy(bridge["remainder"]["R1"]),
            "R2": copy.deepcopy(bridge["remainder"]["R2"]),
        }
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "bound_arithmetic": {
                "pi_upper_rational": "22/7",
                "binary64_inputs_treated_as_exact_rationals": True,
                "final_conversion_directed_upward": True,
            },
            "derived_operator_factors": {
                "product_constant_upper": self.product_constant_upper,
                "source_product_constant_binary64_diagnostic": self.operators.product_constant,
                "radial_integral_I": self.radial_integral_factor(),
                "multiply_Y": self.multiply_Y_factor(),
                "partial_eta_radial_integral": self.eta_radial_integral_factor(rho),
            },
            "diagnostic_inputs_not_source_bounds": {
                "rho": rho,
                "g_norm": g_norm,
                "Phi": {"norm": Phi.norm, "lipschitz": Phi.lipschitz},
                "remainder_h": calculator.h,
                "remainder_A": calculator.A,
                "remainder_D": calculator.D,
                "remainder_lambda": calculator.rescaling_lambda,
            },
            "conditional_pressure_bounds": pressure_json,
            "conditional_remainder_totals_after_pressure_injection": remainder_totals,
            "truth_boundary": self.truth_boundary,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        }
        identity = copy.deepcopy(payload)
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(identity).encode("utf-8")
        ).hexdigest()
        return payload

    def save_report(self, path: str | Path) -> dict[str, Any]:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = self.report()
        destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = KokunoPA10PressureBallBounds().save_report(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
