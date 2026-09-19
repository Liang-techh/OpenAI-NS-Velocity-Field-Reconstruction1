"""Propagate an oversized shared-C choice through the executable Appendix-B path.

KokunoYumeto's corrected 2026-09-09 reconstruction first fixes a constant
``C>1`` with ``C >= sup_Omega |phi_*|`` and later reuses the same ``C`` in
``X_R = 110 (C P_*)^10``.  The earlier executable path used only the minimum
real-axis normalization ``log C = Lambda max_I int zeta_*``.  That selected
choice is excluded by the PA.10 pointwise necessary screen, but the source
choice of C is a lower-bound choice, not an equality.

This module makes one larger, explicitly autonomous choice *upstream* and
propagates it through the rescaled core, reference continuation, and Appendix-B
activation before rerunning the same PA.10 screen.  It does not multiply C
after the fact.  The family parameter is

    log C_selected = log_C_factor * log C_real-axis-minimum,  factor >= 1.

The default factor 32 is a deliberately coarse dyadic oversize.  It is not a
recovered source value and is not a proof of the complex-neighborhood C bound
or of the analytic B_0/T_sh estimates.  Passing the pointwise necessary screen
therefore means only "not excluded"; PA.16 remains fail closed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_coupled_normalization_pa10_screen import (
    KokunoCoupledNormalizationPA10Screen,
)
from .kokuno_rescaled_appendix_b_boundary import (
    KokunoSourceRescaledAppendixBBoundary,
)
from .kokuno_rescaled_appendix_b_selected import (
    make_source_scale_aware_appendix_b_boundary,
)
from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed
from .kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-coupled-c-appendix-b-normalization-v1"

_SOURCE_FORMULAS = {
    "source_C_lower_bound": "choose fixed C>1 with C>=sup_{eta in Omega}|phi_*(eta)|",
    "normalized_swirl": "F=phi/C; E=sqrt(2X)F",
    "boundary_log": "ell_i=log(C E_i)",
    "outer_scale": "X_R=110(C P_*)^10",
    "PA10_T_sh": "T_sh>=20||sigma'||_inf(B_0+||log f||_inf)",
    "PA10_separation": "X_sep/X_R<exp(-8)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "autonomous_upstream_log_C_oversize_executable": True,
    "shared_C_propagated_through_core_reference_appendix_B": True,
    "pointwise_PA10_necessary_screen_executable": True,
    "posthoc_C_multiplier_applied": False,
    "passing_pointwise_screen_is_source_T_sh_certificate": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_hidden_C_recovered": False,
    "source_hidden_numeric_choices_recovered": False,
    "selected_pa16_handoff_allowed": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "amplitude_collapse_used_for_pde_claim": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoCoupledCRescaledCoreSeed(KokunoSourceRescaledCoreSeed):
    """Rescaled core seed with an upstream log-space C oversize.

    The base class uses the minimum real-axis normalization.  This subclass
    keeps every source/autonomous upstream datum fixed and increases C before
    any profile is evaluated.  Since ``F=phi/C``, this shifts ``log F`` by the
    exact constant ``-log_C_oversize`` and lets the downstream pressure
    primitive, reference stage, and Appendix-B ODE see the changed normalized
    angular amplitude.
    """

    log_C_factor: float = 32.0

    def __post_init__(self) -> None:
        super().__post_init__()
        factor = float(self.log_C_factor)
        if not math.isfinite(factor) or not 1.0 <= factor <= 1.0e4:
            raise ValueError("log_C_factor must lie in [1,1e4]")
        object.__setattr__(self, "log_C_factor", factor)

    @property
    def minimal_real_axis_log_C(self) -> float:
        return float(super().log_C_real_axis)

    @property
    def log_C_oversize(self) -> float:
        return (self.log_C_factor - 1.0) * self.minimal_real_axis_log_C

    @property
    def log_C_real_axis(self) -> float:
        value = self.log_C_factor * self.minimal_real_axis_log_C
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("coupled selected log C is outside the finite log range")
        return value

    def log_g(self, eta: Any) -> np.ndarray:
        phase = self.phase(eta)
        raw = (
            self.rescaling_lambda * (phase - self.real_axis_phase_max)
            - self.log_C_oversize
        )
        # Exact phase <= real-axis maximum.  Clip only positive quadrature
        # roundoff relative to the shifted maximum; preserve the C oversize.
        return np.minimum(raw, -self.log_C_oversize)

    def parent_seed_payload(self) -> dict[str, Any]:
        parent = KokunoSourceRescaledCoreSeed(
            binding=self.binding,
            rescaling_lambda_multiplier=self.rescaling_lambda_multiplier,
            pressure_quadrature_points=self.pressure_quadrature_points,
        )
        return parent.to_payload()

    def coupled_payload(self) -> dict[str, Any]:
        payload = {
            "schema": SCHEMA + "-seed",
            "parent_seed": self.parent_seed_payload(),
            "log_C_factor": self.log_C_factor,
            "minimal_real_axis_log_C": self.minimal_real_axis_log_C,
            "selected_log_C": self.log_C_real_axis,
            "log_C_oversize": self.log_C_oversize,
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload


@dataclass(frozen=True)
class KokunoCoupledCAppendixBNormalization:
    """One-dimensional autonomous C family propagated through Appendix B."""

    log_C_factor: float = 32.0
    kappa0_lambda_multiplier: float = 1.0

    _seed: KokunoCoupledCRescaledCoreSeed = field(
        init=False, repr=False, compare=False
    )
    _reference: KokunoSourceRescaledReferenceContinuation = field(
        init=False, repr=False, compare=False
    )
    _boundary: KokunoSourceRescaledAppendixBBoundary = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        factor = float(self.log_C_factor)
        multiplier = float(self.kappa0_lambda_multiplier)
        if not math.isfinite(factor) or not 1.0 <= factor <= 1.0e4:
            raise ValueError("log_C_factor must lie in [1,1e4]")
        if not math.isfinite(multiplier) or not 0.0 < multiplier <= 1.0e6:
            raise ValueError("kappa0_lambda_multiplier must lie in (0,1e6]")
        seed = KokunoCoupledCRescaledCoreSeed(log_C_factor=factor)
        reference = KokunoSourceRescaledReferenceContinuation(seed=seed)
        boundary = make_source_scale_aware_appendix_b_boundary(
            reference=reference,
            kappa0_lambda_multiplier=multiplier,
        )
        object.__setattr__(self, "log_C_factor", factor)
        object.__setattr__(self, "kappa0_lambda_multiplier", multiplier)
        object.__setattr__(self, "_seed", seed)
        object.__setattr__(self, "_reference", reference)
        object.__setattr__(self, "_boundary", boundary)

    @property
    def seed(self) -> KokunoCoupledCRescaledCoreSeed:
        return self._seed

    @property
    def reference(self) -> KokunoSourceRescaledReferenceContinuation:
        return self._reference

    @property
    def boundary(self) -> KokunoSourceRescaledAppendixBBoundary:
        return self._boundary

    @property
    def selected_log_C(self) -> float:
        return float(self.boundary.log_C)

    @property
    def log_P_star(self) -> float:
        return float(self.seed.binding.outer_schedule.log_P_star)

    @property
    def screen(self) -> KokunoCoupledNormalizationPA10Screen:
        return KokunoCoupledNormalizationPA10Screen(
            log_C=self.selected_log_C,
            log_P_star=self.log_P_star,
        )

    @property
    def witness_eta(self) -> float:
        return float(self.seed.phase_stationary_eta)

    def boundary_values(self, eta: Any) -> dict[str, np.ndarray]:
        return self.boundary.boundary_values(eta)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.boundary.velocity(x, y, z, t)

    __call__ = velocity

    def pointwise_report(self) -> dict[str, Any]:
        eta = self.witness_eta
        values = self.boundary_values(np.asarray(eta))
        ell_i = float(np.asarray(values["ell_i"]))
        G_i = float(np.asarray(values["G_i"]))
        log_E_i = ell_i - self.selected_log_C
        evaluation = self.screen.evaluation_report(log_E_i)
        return {
            "witness_eta": eta,
            "minimal_real_axis_log_C": self.seed.minimal_real_axis_log_C,
            "log_C_factor": self.log_C_factor,
            "selected_log_C": self.selected_log_C,
            "log_C_oversize": self.seed.log_C_oversize,
            "log_E_i": log_E_i,
            "ell_i": ell_i,
            "G_i": G_i,
            "materialized_E_i": float(np.asarray(values["E_i"])),
            "materialized_F_i": float(np.asarray(values["F_i"])),
            "materialized_angular_amplitude_underflowed": bool(
                float(np.asarray(values["E_i"])) == 0.0
            ),
            "ell_abs_strict_cap": self.screen.ell_abs_strict_cap,
            "pointwise_T_sh_lower_bound": evaluation[
                "pointwise_T_sh_lower_bound"
            ],
            "t_sh_strict_upper": evaluation["t_sh_strict_upper"],
            "strict_margin": evaluation["strict_margin"],
            "not_excluded_by_necessary_screen": evaluation[
                "not_excluded_by_necessary_screen"
            ],
            "source_B0_analytic_bound_proved": False,
            "source_T_sh_lower_bound_verified": False,
            "selected_pa16_handoff_allowed": False,
        }

    def report(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "configuration_sha256": self.sha256,
            "seed": self.seed.coupled_payload(),
            "pointwise_PA10": self.pointwise_report(),
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "parameters": {
                "log_C_factor": self.log_C_factor,
                "kappa0_lambda_multiplier": self.kappa0_lambda_multiplier,
            },
            "parent_seed_sha256": self.seed.parent_seed_payload()["sha256"],
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> "KokunoCoupledCAppendixBNormalization":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected coupled-C Appendix-B normalization schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("coupled-C Appendix-B source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("coupled-C Appendix-B truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("coupled-C Appendix-B payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing coupled-C Appendix-B parameters")
        obj = cls(**params)
        if obj.to_payload() != payload:
            raise ValueError("coupled-C Appendix-B replay changed payload")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "KokunoCoupledCAppendixBNormalization":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Propagate an autonomous oversized shared-C through Appendix B"
    )
    parser.add_argument("--log-c-factor", type=float, default=32.0)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    candidate = KokunoCoupledCAppendixBNormalization(
        log_C_factor=args.log_c_factor
    )
    text = json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
