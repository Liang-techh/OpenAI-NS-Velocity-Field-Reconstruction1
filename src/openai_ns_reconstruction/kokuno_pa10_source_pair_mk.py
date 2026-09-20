"""Source-compatible PA.10 pair M/K self-certificate from full post-J bounds.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction defines the fixed-point pair

    ((1+T)^-1 J_2 R1 / 2, J_1 R2 / 2),

where ``T=J_2 chi/2``.  Agent-1 #751 now supplies a complete source-compatible
post-``J_2`` R1 self-envelope, while Agent-1 #693 supplies a complete
source-compatible post-``J_1`` R2 self-envelope.  The source-compatible Phi
ball already supplies an absolute-series upper bound for ``||(1+T)^-1||``.

This module performs exactly the next assembly step: it feeds those *already
post-radial-inverse* envelopes into the displayed pair once, obtaining
radius-one norm and Lipschitz envelopes M and K.  It deliberately does not call
an operator helper that would apply ``J_2`` or ``J_1`` again.  A replay against
#668's post-J diagnostic arithmetic is retained as a no-double-J shape guard.

The resulting M/K values are Agent-1 source-compatible self-certificates only.
They are not independent Agent-4 admission, do not establish a contraction or
fixed point, and do not promote B0/T_sh, PA.16, a global leading velocity, or
held-out Navier-Stokes validation.
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

from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_derivative_ordinary_slots import (
    KokunoPA10SourceDerivativeOrdinarySlots,
)
from .kokuno_pa10_source_r1_du_phi_eta_ordinary_slot import (
    KokunoPA10SourceR1DuPhiEtaOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-pair-mk-v1"

_SOURCE_FORMULAS = {
    "T": "T=J_2 chi/2",
    "inverse": (
        "||(1+T)^-1|| <= sum_{k>=0}(40 M_chi)^k/(k!(k+1)!)"
    ),
    "pair": "((1+T)^-1 J_2 R1/2, J_1 R2/2)",
    "M": "radius-one pair norm upper from the displayed fixed-point pair",
    "K": "radius-one pair Lipschitz upper from the same displayed pair",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_full_post_J2_R1_agent1_self_certificate_consumed": True,
    "source_compatible_full_post_J1_R2_agent1_self_certificate_consumed": True,
    "source_compatible_inverse_one_plus_T_machine_bound_consumed": True,
    "post_J2_R1_not_integrated_again": True,
    "post_J1_R2_not_integrated_again": True,
    "source_compatible_pair_M_agent1_self_certificate_machine_bound": True,
    "source_compatible_pair_K_agent1_self_certificate_machine_bound": True,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_operator_constant_M_independent_agent4_admission": False,
    "source_operator_constant_K_independent_agent4_admission": False,
    "source_full_post_J2_R1_independent_admission": False,
    "source_full_post_J1_R2_independent_admission": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_contraction_factor_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_fixed_point_solved": False,
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


@dataclass(frozen=True)
class KokunoPA10SourcePairMK:
    """Assemble Agent-1 full post-J R1/R2 self-envelopes into pair M/K."""

    r1_source: KokunoPA10SourceR1DuPhiEtaOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1DuPhiEtaOrdinarySlot
    )
    r2_source: KokunoPA10SourceDerivativeOrdinarySlots = field(
        default_factory=KokunoPA10SourceDerivativeOrdinarySlots
    )

    def __post_init__(self) -> None:
        if not isinstance(self.r1_source, KokunoPA10SourceR1DuPhiEtaOrdinarySlot):
            raise TypeError("r1_source must be KokunoPA10SourceR1DuPhiEtaOrdinarySlot")
        if not isinstance(self.r2_source, KokunoPA10SourceDerivativeOrdinarySlots):
            raise TypeError("r2_source must be KokunoPA10SourceDerivativeOrdinarySlots")
        if self.r1_source.domain.sha256 != self.r2_source.domain.sha256:
            raise ValueError("R1/R2 source-compatible axis domains do not match")
        if self.r1_source.bridge.sha256 != self.r2_source.bridge.sha256:
            raise ValueError("R1/R2 post-J bridge identities do not match")
        if self.phi_source.sha256 != self.r2_source.bridge.u_source.mixed.phi_source.sha256:
            raise ValueError("R1/R2 Phi source identities do not match")

    @property
    def bridge(self):
        return self.r1_source.bridge

    @property
    def phi_source(self):
        return self.bridge.u_source.mixed.phi_source

    def inverse_one_plus_t_upper(self) -> float:
        certificate = self.phi_source.phi_ball_certificate()
        value = float(certificate["inverse_one_plus_T_absolute_series_upper"])
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("inverse (1+T) bound must be finite and positive")
        return value

    def full_post_j_inputs(self) -> dict[str, BallFactorBound]:
        r1 = self.r1_source.full_r1_after_j2()
        r2 = self.r2_source.full_r2_after_j1()
        if not all(isinstance(item, BallFactorBound) for item in (r1, r2)):
            raise TypeError("full post-J inputs must be BallFactorBound values")
        return {"post_J2_R1": r1, "post_J1_R2": r2}

    def pair_mk_self_certificate(self) -> dict[str, float]:
        """Return M/K directly from already-post-J source-compatible envelopes."""
        values = self.full_post_j_inputs()
        r1 = values["post_J2_R1"]
        r2 = values["post_J1_R2"]
        inverse = _q_nonnegative(
            self.inverse_one_plus_t_upper(),
            "inverse_one_plus_T_absolute_series_upper",
        )
        half = Fraction(1, 2)
        m_phi = half * inverse * _q_nonnegative(r1.norm, "post J2 R1 norm")
        m_u = half * _q_nonnegative(r2.norm, "post J1 R2 norm")
        k_phi = half * inverse * _q_nonnegative(
            r1.lipschitz, "post J2 R1 lipschitz"
        )
        k_u = half * _q_nonnegative(r2.lipschitz, "post J1 R2 lipschitz")
        return {
            "inverse_one_plus_T_absolute_series_upper": _upper(inverse),
            "M_phi_component": _upper(m_phi),
            "M_u_component": _upper(m_u),
            "M_pair_max": _upper(max(m_phi, m_u)),
            "K_phi_component": _upper(k_phi),
            "K_u_component": _upper(k_u),
            "K_pair_max": _upper(max(k_phi, k_u)),
        }

    def bridge_post_j_replay(self) -> dict[str, float]:
        """Replay #668's post-J arithmetic as a no-double-J shape guard."""
        values = self.full_post_j_inputs()
        assembly = {
            "post_J2_R1": asdict(values["post_J2_R1"]),
            "post_J1_R2": asdict(values["post_J1_R2"]),
        }
        return self.bridge.diagnostic_pair_envelope_from_post_j(assembly)

    def replay_matches_self_certificate(self) -> bool:
        public = self.pair_mk_self_certificate()
        replay = self.bridge_post_j_replay()
        correspondence = {
            "inverse_one_plus_T_absolute_series_upper": (
                "inverse_one_plus_T_absolute_series_upper"
            ),
            "M_phi_component": "diagnostic_M_phi_component",
            "M_u_component": "diagnostic_M_u_component",
            "M_pair_max": "diagnostic_M_pair_max",
            "K_phi_component": "diagnostic_K_phi_component",
            "K_u_component": "diagnostic_K_u_component",
            "K_pair_max": "diagnostic_K_pair_max",
        }
        return all(public[left] == replay[right] for left, right in correspondence.items())

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        post_j = self.full_post_j_inputs()
        pair = self.pair_mk_self_certificate()
        replay = self.bridge_post_j_replay()
        if not self.replay_matches_self_certificate():
            raise RuntimeError("post-J pair replay drifted from source M/K self-certificate")
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
                "R1_full_post_J2_agent1_receipt_sha256": self.r1_source.sha256,
                "R2_full_post_J1_agent1_receipt_sha256": self.r2_source.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "Phi_ball_receipt_sha256": self.phi_source.sha256,
                "source_axis_receipt_sha256": self.r1_source.domain.sha256,
            },
            "operator_certificate": {
                "R1_input_is_already_post_J2": True,
                "R2_input_is_already_post_J1": True,
                "J2_reapplied_in_this_increment": False,
                "J1_reapplied_in_this_increment": False,
                "inverse_one_plus_T_source_ball_bound_reused": True,
                "pair_factor_one_half_applied_exactly_once": True,
                "bridge_post_J_replay_exact_match": True,
                "diagnostic_unit_ordinary_fixture_used": False,
            },
            "full_post_J_source_compatible_agent1_self_envelopes": {
                name: asdict(bound) for name, bound in post_j.items()
            },
            "source_compatible_pair_MK_agent1_self_certificate": pair,
            "post_J_bridge_arithmetic_replay": replay,
            "progress": {
                "R1_ordinary_post_J2_source_bound_count": 11,
                "R1_ordinary_total": 11,
                "R2_ordinary_post_J1_source_bound_count": 11,
                "R2_ordinary_total": 11,
                "full_R1_agent1_self_certificate_materialized": True,
                "full_R2_agent1_self_certificate_materialized": True,
                "pair_M_agent1_self_certificate_materialized": True,
                "pair_K_agent1_self_certificate_materialized": True,
                "full_R1_independent_admission_still_open": True,
                "full_R2_independent_admission_still_open": True,
                "M_K_independent_admission_still_open": True,
                "contraction_not_promoted": True,
                "fixed_point_not_promoted": True,
            },
            "independent_audit_boundary": {
                "R1_full_self_envelope_is_not_independent_admission": True,
                "R2_full_self_envelope_is_not_independent_admission": True,
                "M_K_self_certificate_is_not_independent_admission": True,
                "Agent4_independent_rebuild_required_before_promotion": True,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "The M/K values are conservative Agent-1 coefficient-space self-bounds, not an NS residual.",
                "No J2/J1 is applied again because the consumed R1/R2 envelopes are already post-J.",
                "Independent Agent-4 admission is still required before source operator constants can be promoted.",
                "No contraction, fixed point, B0/T_sh, PA.16, global velocity, or PDE gate is established here.",
                "No hidden Kokuno/OpenAI parameter is recovered or claimed paper-exact.",
            ],
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
    payload = KokunoPA10SourcePairMK().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print(
        "full_post_J_source_envelopes=",
        payload["full_post_J_source_compatible_agent1_self_envelopes"],
    )
    print(
        "source_compatible_pair_MK_agent1_self_certificate=",
        payload["source_compatible_pair_MK_agent1_self_certificate"],
    )
    print("progress=", payload["progress"])


if __name__ == "__main__":
    _main()
