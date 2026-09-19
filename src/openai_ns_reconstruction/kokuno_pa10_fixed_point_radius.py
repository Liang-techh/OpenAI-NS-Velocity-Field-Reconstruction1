"""Executable guard for Kokuno's rescaled-core fixed-point distance estimate.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

The corrected reconstruction places the rescaled-core map

    (Phi,u) -> (Phi_0 + (1+T)^-1 J_2 R_1/(2 Lambda),
                u_0 + J_1 R_2/(2 Lambda))

in a source-defined Banach ball, proves a contraction factor at most 1/2,
and states that its unique fixed point has norm distance at most ``M/Lambda``
from the center.  ``M`` is built from finite operator constants in that proof.

This module makes *that displayed radius relation* executable against the
repository's selected Lambda.  It deliberately does not invent a numerical
value for the source operator constant ``M`` and does not promote the selected
autonomous Lambda into the source's existential large-parameter threshold.
Therefore every returned radius is conditional on the source contraction
hypotheses, and all source PA.10 B_0/T_sh gates remain fail closed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-fixed-point-radius-v1"

_SOURCE_FORMULAS = {
    "fixed_point_map": (
        "(Phi,u)->(Phi_0+(1+T)^-1 J_2 R_1/(2 Lambda),"
        "u_0+J_1 R_2/(2 Lambda))"
    ),
    "contraction_factor": "q<=1/2",
    "fixed_point_distance": "||(Phi,u)-(Phi_0,u_0)||_rho <= M/Lambda",
    "M_provenance": "finite operator constants in the corrected reconstruction",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_fixed_point_radius_formula_recorded": True,
    "source_fixed_point_radius_formula_executable": True,
    "selected_lambda_executable": True,
    "selected_lambda_is_source_existential_threshold": False,
    "source_operator_constant_M_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_contraction_factor_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_fixed_point_solved": False,
    "source_M0_machine_bound": False,
    "combined_profile_C0_norm_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoPA10FixedPointRadiusGate:
    """Conditional executable form of the source ``M/Lambda`` radius.

    ``operator_constant_M`` is intentionally optional.  Supplying a finite
    positive value allows numerical sensitivity/replay of the displayed source
    formula, but it is *not* treated as a source-certified constant unless a
    future independent machine proof binds the finite operator constants.  The
    current schema therefore always keeps the source truth flags false.
    """

    seed: KokunoSourceRescaledCoreSeed = field(default_factory=KokunoSourceRescaledCoreSeed)
    operator_constant_M: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.seed, KokunoSourceRescaledCoreSeed):
            raise TypeError("seed must be a KokunoSourceRescaledCoreSeed")
        if self.operator_constant_M is not None:
            value = float(self.operator_constant_M)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError("operator_constant_M must be finite and >0")
            object.__setattr__(self, "operator_constant_M", value)

    @property
    def selected_lambda(self) -> float:
        return float(self.seed.rescaling_lambda)

    @property
    def conditional_radius(self) -> float | None:
        """Return ``M/Lambda`` when an explicit diagnostic ``M`` is supplied."""
        if self.operator_constant_M is None:
            return None
        radius = self.operator_constant_M / self.selected_lambda
        if not math.isfinite(radius) or radius < 0.0:
            raise OverflowError("conditional M/Lambda radius is not finite")
        return float(radius)

    def radius_for_M(self, operator_constant_M: float) -> float:
        """Evaluate the source-displayed radius for an explicit finite ``M``."""
        value = float(operator_constant_M)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("operator_constant_M must be finite and >0")
        radius = value / self.selected_lambda
        if not math.isfinite(radius) or radius < 0.0:
            raise OverflowError("conditional M/Lambda radius is not finite")
        return float(radius)

    def M_budget_for_radius(self, target_radius: float) -> float:
        """Return the largest ``M`` compatible with ``M/Lambda <= target``.

        This is a planning/sensitivity quantity only; it does not certify the
        source operator constants or the source large-Lambda threshold.
        """
        target = float(target_radius)
        if not math.isfinite(target) or target <= 0.0:
            raise ValueError("target_radius must be finite and >0")
        budget = target * self.selected_lambda
        if not math.isfinite(budget):
            raise OverflowError("M budget is outside float range")
        return float(budget)

    @cached_property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        radius = self.conditional_radius
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "selected_execution": {
                "Lambda": self.selected_lambda,
                "operator_constant_M_input": self.operator_constant_M,
                "conditional_fixed_point_radius_M_over_Lambda": radius,
                "contraction_factor_source_upper_bound": 0.5,
                "interpretation": (
                    "formula replay only; conditional on source invariant-ball and "
                    "contraction hypotheses"
                ),
            },
            "truth_boundary": self.truth_boundary,
        }

    def to_payload(self) -> dict[str, Any]:
        payload = self.report()
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode()).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.to_payload()["sha256"])

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n")

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10FixedPointRadiusGate":
        raw = copy.deepcopy(payload)
        supplied_sha = raw.pop("sha256", None)
        expected_sha = hashlib.sha256(_canonical_json(raw).encode()).hexdigest()
        if supplied_sha != expected_sha:
            raise ValueError("payload sha256 mismatch")
        if raw.get("schema") != SCHEMA:
            raise ValueError("schema mismatch")
        if raw.get("source") != cls().report()["source"]:
            raise ValueError("source provenance mismatch")
        if raw.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula metadata mismatch")
        if raw.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth boundary mismatch")
        selected = raw.get("selected_execution", {})
        value = selected.get("operator_constant_M_input")
        replay = cls(operator_constant_M=value)
        if replay.report()["selected_execution"] != selected:
            raise ValueError("selected execution replay mismatch")
        return replay

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPA10FixedPointRadiusGate":
        return cls.from_payload(json.loads(Path(path).read_text()))
