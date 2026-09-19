"""Source-compatible PA.10 pressure ordinary slots after the radial inverse.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction contains the pressure part of the second PA.10
remainder

    L R2 = ... - 4 A eta p + d partial_eta p - 2 eta Y p_Y,
    p = I(g^2 Phi^2),   d = 1-eta^2.

Agent-1 #635 made the pressure-map coefficient algebra executable, #653 bound
the source-compatible rho/g/Phi inputs, and #668 introduced an exhaustive
post-J1 R2 ordinary-term ledger.  This module fills exactly the three pressure
slots in that ledger with source-compatible coefficient-space envelopes and
then applies #668's fixed L^-1 and J1 factors.  It does not fabricate any of
the remaining eight ordinary R2 slots and therefore does not promote full R2,
M/K, a fixed point, global pressure/velocity, or PDE validation.

The selected sigma/rho datum is a repository-autonomous source-compatible
existence choice, not recovered hidden Kokuno/OpenAI data.  Independent Agent-4
admission of the source-axis prerequisite is still required.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, field
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_post_j_remainder_bridge import KokunoPA10PostJRemainderBridge
from .kokuno_pa10_pressure_ball_bounds import KokunoPA10PressureBallBounds
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-pressure-ordinary-slots-v1"

R2_PRESSURE_ORDINARY_TERMS = (
    "4A_eta_p",
    "d_p_eta",
    "2eta_Y_p_Y",
)

_SOURCE_FORMULAS = {
    "pressure": "p=I(g^2 Phi^2)",
    "pressure_eta": "partial_eta p=partial_eta I(g^2 Phi^2)",
    "pressure_radial": "Y p_Y=Y g^2 Phi^2",
    "R2_pressure": "-4A eta p + d p_eta - 2 eta Y p_Y",
    "A": "A=1/2+h",
    "d": "d=1-eta^2",
    "post_J1": "J_1[L^-1 R2_numerator]",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_pressure_ball_consumed": True,
    "source_compatible_eta_coefficient_norm_machine_bound": True,
    "source_compatible_R2_pressure_ordinary_numerator_slots_machine_bound": True,
    "source_compatible_R2_pressure_ordinary_post_J1_slots_machine_bound": True,
    "source_axis_domain_independent_agent4_admission_required": True,
    "source_Phi_ball_independent_agent4_admission_required": False,
    "all_R2_ordinary_slots_source_bound": False,
    "all_R1_ordinary_slots_source_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "heldout_ns_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _q_nonnegative(value: float, name: str) -> Fraction:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return Fraction.from_float(out)


def _upper(value: Fraction) -> float:
    if value < 0:
        raise ValueError("upper-bound helper requires a nonnegative value")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper bound is outside binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return float(out)


def _scale_exact(bound: BallFactorBound, scalar: Fraction) -> BallFactorBound:
    if not isinstance(bound, BallFactorBound):
        raise TypeError("bound must be a BallFactorBound")
    if scalar < 0:
        scalar = -scalar
    return BallFactorBound(
        _upper(scalar * _q_nonnegative(bound.norm, "bound norm")),
        _upper(scalar * _q_nonnegative(bound.lipschitz, "bound lipschitz")),
    )


@dataclass(frozen=True)
class KokunoPA10SourcePressureOrdinarySlots:
    """Bind exactly the three pressure ordinary slots in the post-J1 R2 ledger."""

    bridge: KokunoPA10PostJRemainderBridge = field(
        default_factory=KokunoPA10PostJRemainderBridge
    )

    def __post_init__(self) -> None:
        if not isinstance(self.bridge, KokunoPA10PostJRemainderBridge):
            raise TypeError("bridge must be KokunoPA10PostJRemainderBridge")

    @property
    def domain(self):
        return self.bridge.domain

    @property
    def phi_source(self):
        return self.bridge.u_source.mixed.phi_source

    @property
    def pressure_calculator(self) -> KokunoPA10PressureBallBounds:
        return KokunoPA10PressureBallBounds()

    def eta_coefficient_ball(self) -> BallFactorBound:
        """Exact finite-derivative coefficient bound for fixed ``eta``.

        For alpha=0 the source weight gives ``a_{0,0}=1`` and
        ``a_{0,1}=1/(4 rho)``.  On the enlarged real interval
        ``|eta|<=E=1+margin`` the only nonzero derivatives are eta and 1, hence

            ||eta||_rho <= max(E, 4 rho).
        """
        rho = Fraction.from_float(float(self.domain.coefficient_rho))
        margin = Fraction.from_float(float(self.domain.enlarged_real_margin))
        if rho <= 0 or margin < 0:
            raise ValueError("source-compatible rho/margin must be valid")
        E = Fraction(1) + margin
        return BallFactorBound(_upper(max(E, 4 * rho)), 0.0)

    def source_pressure_ball(self) -> dict[str, BallFactorBound]:
        """Replay #635 with #653's source-compatible rho/g/Phi inputs."""
        inputs = self.phi_source.completed_pressure_operator_inputs()
        if not inputs["ready_for_source_compatible_pressure_operator"]:
            raise ValueError("source-compatible pressure inputs are not ready")
        return self.pressure_calculator.pressure_bounds(
            rho=float(inputs["rho"]),
            g_norm=float(inputs["g_norm_upper"]),
            Phi=BallFactorBound(
                norm=float(inputs["Phi_radius_one_ball_norm"]),
                lipschitz=float(inputs["Phi_radius_one_ball_lipschitz"]),
            ),
        )

    def r2_pressure_numerator_slots(self) -> dict[str, BallFactorBound]:
        """Return the three source-compatible pressure numerator envelopes."""
        pressure = self.source_pressure_ball()
        eta = self.eta_coefficient_ball()
        d = self.bridge.u_source.mixed.d_coefficient_ball()
        calc = self.pressure_calculator

        h = Fraction.from_float(float(self.domain.h))
        A = Fraction(1, 2) + h
        four_A = 4 * A

        eta_p = calc.product(eta, pressure["p"])
        d_p_eta = calc.product(d, pressure["p_eta"])
        eta_Y_p_Y = calc.product(eta, pressure["Y_p_Y"])
        return {
            "4A_eta_p": _scale_exact(eta_p, four_A),
            "d_p_eta": d_p_eta,
            "2eta_Y_p_Y": _scale_exact(eta_Y_p_Y, Fraction(2)),
        }

    def r2_pressure_after_j1(self) -> dict[str, BallFactorBound]:
        """Apply #668's L^-1 and J1 exactly once to the three pressure slots."""
        slots = self.r2_pressure_numerator_slots()
        return {
            name: self.bridge.ordinary_term_after_jnu(bound, nu=1)
            for name, bound in slots.items()
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta = self.eta_coefficient_ball()
        pressure = self.source_pressure_ball()
        numerator = self.r2_pressure_numerator_slots()
        after_j1 = self.r2_pressure_after_j1()
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
            "upstream_identity": {
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_Phi_ball_receipt_sha256": self.phi_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "h": float(self.domain.h),
                "A_exact_from_selected_h": str(
                    Fraction(1, 2) + Fraction.from_float(float(self.domain.h))
                ),
                "rho": float(self.domain.coefficient_rho),
                "enlarged_real_margin": float(self.domain.enlarged_real_margin),
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "fixed_multiplier_bounds": {
                "eta": asdict(eta),
                "d": asdict(self.bridge.u_source.mixed.d_coefficient_ball()),
                "L_inverse": asdict(self.bridge.l_inverse_coefficient_ball()),
            },
            "source_compatible_pressure_ball": {
                name: asdict(bound) for name, bound in pressure.items()
            },
            "R2_pressure_ordinary_numerator_slots": {
                name: asdict(bound) for name, bound in numerator.items()
            },
            "R2_pressure_ordinary_after_J1": {
                name: asdict(bound) for name, bound in after_j1.items()
            },
            "ledger_progress": {
                "R2_pressure_slots_filled": list(R2_PRESSURE_ORDINARY_TERMS),
                "R2_pressure_slots_filled_count": len(R2_PRESSURE_ORDINARY_TERMS),
                "R2_total_ordinary_slots": 11,
                "R2_other_ordinary_slots_still_required": 8,
                "R1_ordinary_slots_still_required": 11,
                "mixed_slots_handled_upstream_by_A1_668": True,
                "full_R1_R2_not_promoted": True,
            },
            "independent_audit_boundary": {
                "A4_source_Phi_ball_conditional_PASS_exists": True,
                "A4_source_axis_repair_PR": 671,
                "A4_source_axis_admission_required_before_downstream_promotion": True,
                "this_receipt_is_agent1_self_certificate_not_independent_validation": True,
            },
            "truth_boundary": self.truth_boundary,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.report()["receipt_sha256"])

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            KokunoPA10SourcePressureOrdinarySlots().save_report(args.output),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    _main()
