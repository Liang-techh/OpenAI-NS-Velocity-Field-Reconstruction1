"""Executable selected-center envelope for the last PA.10 profile-norm dependency.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

The corrected reconstruction uses, before the later choice of ``C``,

    ||log phi_* + log Phi(4,.)||_{C^0_eta}

inside the displayed PA.10 bound.  For the repository's selected executable
*contraction center* the source rescaling is

    phi = phi_* Phi,             phi_* = exp(Lambda int_0^eta zeta_*),
    Phi_0(Y,eta) = f_0(Y chi),   Y = Lambda X,

so at ``X_0=4/Lambda`` the requested quantity is evaluated stably as

    b(eta) = Lambda int_0^eta zeta_*(s) ds + log f_0(4 chi(eta)).

This avoids forming ``log_C + log_F`` at the enormous selected scale.  The
eta derivative is also executable directly from the displayed center formulas,

    b'(eta) = Lambda zeta_*(eta)
              + 4 chi'(eta) f_0'(4 chi)/f_0(4 chi).

The code locates derivative sign changes on nested eta grids, refines roots
with Brent, and pads the observed maximum.  That is a reproducible engineering
envelope for the selected center, NOT an interval proof of the continuum
supremum and NOT a bound for Kokuno's source fixed point.  Consequently it does
not set the source PA.10 B_0/T_sh gates or open PA.16.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq

from .kokuno_pa10_l0_binding import KokunoPA10L0Binding
from .kokuno_pa10_selected_m0_envelope import KokunoPA10SelectedM0Envelope
from .kokuno_pa10_source_choice_order_guard import KokunoPA10SourceChoiceOrderGuard
from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed
from .kokuno_rescaled_reference_continuation import KokunoSourceRescaledReferenceContinuation


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-selected-profile-norm-v1"
X_I = 110.0

_SOURCE_FORMULAS = {
    "PA10_dependency": "||log phi_*+log Phi(4,.)||_{C^0_eta}",
    "source_rescaling": "phi=phi_* Phi; phi_*=exp(Lambda*int_0^eta zeta_*); Y=Lambda X",
    "selected_center": "Phi_0(Y,eta)=f_0(Y chi(eta))",
    "stable_boundary_profile": "b(eta)=Lambda*int_0^eta zeta_*+log(f_0(4 chi(eta)))",
    "selected_center_eta_derivative": (
        "b_eta=Lambda*zeta_*+4*chi_eta*f_0'(4chi)/f_0(4chi)"
    ),
    "B0_context": (
        "B0 = ||log phi_*+log Phi(4,.)||_C0 + "
        "0.5(M0+0.8)L0 + 0.5 log(2 Xi)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "combined_profile_C0_formula_recorded": True,
    "selected_center_combined_profile_executable": True,
    "selected_center_combined_profile_eta_derivative_executable": True,
    "selected_center_profile_norm_numerical_envelope_materialized": True,
    "later_C_cancels_from_selected_center_combined_profile": True,
    "selected_center_is_source_fixed_point": False,
    "selected_profile_norm_envelope_is_continuum_interval_proof": False,
    "combined_profile_C0_norm_machine_bound": False,
    "source_M0_machine_bound": False,
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


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoPA10SelectedProfileNormEnvelope:
    """Numerical C0 envelope for the selected contraction-center PA.10 datum."""

    seed: KokunoSourceRescaledCoreSeed = field(default_factory=KokunoSourceRescaledCoreSeed)
    coarse_eta_points: int = 65
    fine_eta_points: int = 129
    safety_fraction: float = 1.0e-10
    max_refinement_drift_fraction: float = 1.0e-8
    root_xtol: float = 1.0e-13

    def __post_init__(self) -> None:
        if not isinstance(self.seed, KokunoSourceRescaledCoreSeed):
            raise TypeError("seed must be a KokunoSourceRescaledCoreSeed")
        for name, value, minimum in (
            ("coarse_eta_points", self.coarse_eta_points, 17),
            ("fine_eta_points", self.fine_eta_points, 33),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer")
            if int(value) < minimum or int(value) % 2 == 0:
                raise ValueError(f"{name} must be an odd integer >= {minimum}")
            object.__setattr__(self, name, int(value))
        if self.fine_eta_points <= self.coarse_eta_points:
            raise ValueError("fine_eta_points must exceed coarse_eta_points")
        safety = float(self.safety_fraction)
        drift = float(self.max_refinement_drift_fraction)
        xtol = float(self.root_xtol)
        if not math.isfinite(safety) or not 1.0e-14 <= safety <= 1.0e-3:
            raise ValueError("safety_fraction must lie in [1e-14,1e-3]")
        if not math.isfinite(drift) or not 1.0e-12 <= drift <= 1.0e-3:
            raise ValueError("max_refinement_drift_fraction must lie in [1e-12,1e-3]")
        if not math.isfinite(xtol) or not 1.0e-15 <= xtol <= 1.0e-8:
            raise ValueError("root_xtol must lie in [1e-15,1e-8]")
        object.__setattr__(self, "safety_fraction", safety)
        object.__setattr__(self, "max_refinement_drift_fraction", drift)
        object.__setattr__(self, "root_xtol", xtol)

    def combined_log_profile(self, eta: Any) -> np.ndarray:
        """Return selected-center ``log phi_* + log Phi(4,eta)`` stably."""
        e = _finite_array(eta, "eta")
        if np.any(np.abs(e) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        state = self.seed.axis_state(e)
        phi0 = self.seed.f0(4.0 * state["chi"])
        if np.any(phi0 <= 0.0) or np.any(~np.isfinite(phi0)):
            raise RuntimeError("selected center f_0(4 chi) left the positive source range")
        return self.seed.rescaling_lambda * self.seed.phase(e) + np.log(phi0)

    def combined_log_profile_eta(self, eta: Any) -> np.ndarray:
        """Analytic eta derivative of the selected-center combined profile."""
        e = _finite_array(eta, "eta")
        if np.any(np.abs(e) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        state = self.seed.axis_state(e)
        z = 4.0 * state["chi"]
        phi0 = self.seed.f0(z)
        phi0_prime = self.seed.f0_prime(z)
        if np.any(phi0 <= 0.0) or np.any(~np.isfinite(phi0)):
            raise RuntimeError("selected center f_0(4 chi) left the positive source range")
        return (
            self.seed.rescaling_lambda * state["zeta_star"]
            + 4.0 * state["chi_eta"] * phi0_prime / phi0
        )

    def _eta_grid(self, points: int) -> np.ndarray:
        # Uniform spacing makes derivative sign-change coverage directly comparable
        # across the nested 65/129 grids; include known source-center landmarks.
        grid = np.linspace(-1.0, 1.0, int(points), dtype=float)
        landmarks = np.asarray(
            [-1.0, float(self.seed.phase_stationary_eta), 0.0, 1.0], dtype=float
        )
        return np.unique(np.concatenate((grid, landmarks)))

    def _stationary_points(self, points: int) -> np.ndarray:
        grid = self._eta_grid(points)
        deriv = np.asarray(self.combined_log_profile_eta(grid), dtype=float)
        roots: list[float] = []
        for left, right, dl, dr in zip(
            grid[:-1], grid[1:], deriv[:-1], deriv[1:], strict=True
        ):
            if dl == 0.0:
                roots.append(float(left))
            if dl * dr < 0.0:
                root = brentq(
                    lambda value: float(
                        np.asarray(self.combined_log_profile_eta(np.asarray(value)))
                    ),
                    float(left),
                    float(right),
                    xtol=self.root_xtol,
                    rtol=4.0 * np.finfo(float).eps,
                    maxiter=100,
                )
                roots.append(float(root))
        if deriv[-1] == 0.0:
            roots.append(float(grid[-1]))
        candidates = np.asarray([-1.0, *roots, 1.0], dtype=float)
        return np.unique(np.round(candidates, 15))

    def _sample(self, points: int) -> dict[str, Any]:
        grid = self._eta_grid(points)
        stationary = self._stationary_points(points)
        candidates = np.unique(np.concatenate((grid, stationary)))
        values = np.asarray(self.combined_log_profile(candidates), dtype=float)
        index = int(np.argmax(np.abs(values)))
        deriv_grid = np.asarray(self.combined_log_profile_eta(grid), dtype=float)
        return {
            "eta_grid_size": int(grid.size),
            "stationary_point_count": int(stationary.size),
            "stationary_points": stationary.tolist(),
            "candidate_count": int(candidates.size),
            "sample_max_abs_combined_log_profile": float(abs(values[index])),
            "sample_extremizer_eta": float(candidates[index]),
            "sample_extremizer_value": float(values[index]),
            "sample_min_combined_log_profile": float(np.min(values)),
            "sample_max_combined_log_profile": float(np.max(values)),
            "sample_max_abs_eta_derivative": float(np.max(np.abs(deriv_grid))),
        }

    @cached_property
    def envelope(self) -> dict[str, Any]:
        coarse = self._sample(self.coarse_eta_points)
        fine = self._sample(self.fine_eta_points)
        coarse_max = float(coarse["sample_max_abs_combined_log_profile"])
        fine_max = float(fine["sample_max_abs_combined_log_profile"])
        scale = max(1.0, coarse_max, fine_max)
        drift = abs(fine_max - coarse_max) / scale
        candidate = max(coarse_max, fine_max) * (1.0 + self.safety_fraction)
        guard = bool(
            math.isfinite(candidate)
            and candidate > fine_max
            and drift <= self.max_refinement_drift_fraction
        )
        return {
            "coarse": coarse,
            "fine": fine,
            "nested_grid_relative_drift": drift,
            "max_allowed_refinement_drift_fraction": self.max_refinement_drift_fraction,
            "safety_fraction": self.safety_fraction,
            "combined_profile_C0_candidate_upper_envelope": candidate,
            "numerical_envelope_guard_passed": guard,
            "continuum_supremum_proved": False,
        }

    def candidate_B0_diagnostic(self, M0: float, L0: float | None = None) -> float:
        if not self.envelope["numerical_envelope_guard_passed"]:
            raise RuntimeError("selected profile-norm envelope failed its refinement guard")
        if L0 is None:
            L0 = KokunoPA10L0Binding(seed=self.seed).L_0
        return KokunoPA10SourceChoiceOrderGuard(X_i=X_I).displayed_B0_upper_bound(
            combined_profile_C0_norm=float(
                self.envelope["combined_profile_C0_candidate_upper_envelope"]
            ),
            M0=float(M0),
            L0=float(L0),
        )

    def report(self, include_selected_B0_diagnostic: bool = True) -> dict[str, Any]:
        report: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "selected_seed_sha256": self.seed.sha256,
            "selected_rescaling_lambda": self.seed.rescaling_lambda,
            "selected_X0": 4.0 / self.seed.rescaling_lambda,
            "selected_profile_norm_envelope": copy.deepcopy(self.envelope),
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }
        if include_selected_B0_diagnostic:
            reference = KokunoSourceRescaledReferenceContinuation(seed=self.seed)
            m0 = KokunoPA10SelectedM0Envelope(reference=reference)
            M0 = float(m0.envelope["M0_candidate_upper_envelope"])
            L0 = float(KokunoPA10L0Binding(seed=self.seed).L_0)
            guard = KokunoPA10SourceChoiceOrderGuard(X_i=X_I)
            B0 = self.candidate_B0_diagnostic(M0=M0, L0=L0)
            report["selected_candidate_only_B0_diagnostic"] = {
                "M0_candidate_upper_envelope": M0,
                "L0": L0,
                "B0_from_three_selected_engineering_dependencies": B0,
                "T_sh_lower_from_selected_candidate_B0": (
                    guard.displayed_T_sh_lower_bound_from_B0(B0)
                ),
                "source_B0_analytic_bound_proved": False,
                "source_T_sh_lower_bound_verified": False,
            }
        return report

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": self.report(include_selected_B0_diagnostic=False)["source"],
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "seed": self.seed.to_payload(),
            "parameters": {
                "coarse_eta_points": self.coarse_eta_points,
                "fine_eta_points": self.fine_eta_points,
                "safety_fraction": self.safety_fraction,
                "max_refinement_drift_fraction": self.max_refinement_drift_fraction,
                "root_xtol": self.root_xtol,
            },
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10SelectedProfileNormEnvelope":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected PA.10 selected profile-norm schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("PA.10 selected profile-norm source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("PA.10 selected profile-norm truth boundary changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("PA.10 selected profile-norm payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing PA.10 selected profile-norm parameters")
        obj = cls(
            seed=KokunoSourceRescaledCoreSeed.from_payload(payload.get("seed")),
            coarse_eta_points=int(params["coarse_eta_points"]),
            fine_eta_points=int(params["fine_eta_points"]),
            safety_fraction=float(params["safety_fraction"]),
            max_refinement_drift_fraction=float(
                params["max_refinement_drift_fraction"]
            ),
            root_xtol=float(params["root_xtol"]),
        )
        if obj.to_payload() != payload:
            raise ValueError("PA.10 selected profile-norm replay changed payload")
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
    def load_json(cls, path: str | Path) -> "KokunoPA10SelectedProfileNormEnvelope":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit selected-center PA.10 combined-profile norm receipt"
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--skip-selected-b0-diagnostic",
        action="store_true",
        help="omit the slower candidate-only M0/L0/B0 diagnostic",
    )
    args = parser.parse_args(argv)

    envelope = KokunoPA10SelectedProfileNormEnvelope()
    report = envelope.report(
        include_selected_B0_diagnostic=not args.skip_selected_b0_diagnostic
    )
    report["profile_norm_sha256"] = envelope.sha256
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
