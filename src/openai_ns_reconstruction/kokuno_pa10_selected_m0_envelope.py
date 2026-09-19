"""Numerical selected-reference envelope for the PA.10 ``M_0`` dependency.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

The corrected reconstruction writes ``M_k,N_k`` for bounds on the reference
primitives ``p_{1,r},n_{s,r}`` before forming the PA.10 bound

    ||log phi(X_i)||_{C^k_eta}
      <= ||log phi_* + log Phi(4,.)||_{C^k_eta}
         + 1/2 (M_k + 0.8) L_0.

For ``k=0`` this module makes the repository's *selected executable reference*
``p_{1,r}`` directly inspectable.  It solves the displayed reference ODE

    D_X p_{1,r} = X S_{q,r}/L - l_r p_{1,r}

in ``y=log(X/X_0)`` on the full ``X_0 <= X <= X_i=110`` interval, samples a
nested Chebyshev/log-radial grid, and records a deliberately padded engineering
envelope.  The padding and nested-grid drift guard are numerical diagnostics;
they are NOT an interval proof of the continuum supremum.

This distinction matters because the current selected reference is built from
the executable contraction *center*, while the public reconstruction's fixed
point and analytic coefficient bounds have not been machine-certified here.
Accordingly this module does not set ``source_M0_machine_bound`` and cannot open
the PA.16 handoff.  Its purpose is to supply a reproducible candidate-side
scale for the next bound/certificate step without laundering a replay check
into independent PDE validation.
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
from scipy.integrate import solve_ivp

from .kokuno_pa10_l0_binding import KokunoPA10L0Binding
from .kokuno_pa10_source_choice_order_guard import KokunoPA10SourceChoiceOrderGuard
from .kokuno_rescaled_appendix_b_boundary import KokunoSourceRescaledAppendixBBoundary
from .kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-selected-m0-envelope-v1"
X_I = 110.0

_SOURCE_FORMULAS = {
    "reference_primitive": "D_X p1_r=X*S_q,r/L-l_r*p1_r",
    "source_bound_name": "M_0 is a reference C^0_eta bound on p1_r",
    "B0_context": (
        "||log phi(X_i)||_{C^0_eta} <= "
        "||log phi_*+log Phi(4,.)||_{C^0_eta} + 0.5(M_0+0.8)L_0"
    ),
    "radial_coordinate": "y=log(X/X_0), X_0=4/Lambda, X_i=110",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_M0_definition_recorded": True,
    "selected_reference_p1_ODE_executable": True,
    "selected_reference_full_radial_interval_sampled": True,
    "selected_reference_M0_numerical_envelope_materialized": True,
    "nested_grid_drift_guard_executable": True,
    "selected_reference_is_source_fixed_point": False,
    "selected_M0_envelope_is_continuum_interval_proof": False,
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


def _source_payload() -> dict[str, str]:
    return {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_COMMIT,
        "path": SOURCE_PATH,
        "corrected_release": CORRECTED_RELEASE,
        "corrected_release_date": CORRECTED_RELEASE_DATE,
    }


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoPA10SelectedM0Envelope:
    """Nested-grid engineering envelope for selected-reference ``p_{1,r}``."""

    reference: KokunoSourceRescaledReferenceContinuation = field(
        default_factory=KokunoSourceRescaledReferenceContinuation
    )
    coarse_eta_points: int = 9
    fine_eta_points: int = 17
    coarse_outer_y_points: int = 17
    fine_outer_y_points: int = 33
    safety_fraction: float = 0.10
    max_refinement_drift_fraction: float = 0.05
    rtol: float = 5.0e-10
    atol: float = 5.0e-12
    max_step: float = 0.50

    _helper: KokunoSourceRescaledAppendixBBoundary = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoSourceRescaledReferenceContinuation):
            raise TypeError("reference must be a KokunoSourceRescaledReferenceContinuation")
        integer_fields = (
            ("coarse_eta_points", self.coarse_eta_points, 9),
            ("fine_eta_points", self.fine_eta_points, 17),
            ("coarse_outer_y_points", self.coarse_outer_y_points, 17),
            ("fine_outer_y_points", self.fine_outer_y_points, 33),
        )
        for name, value, minimum in integer_fields:
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer")
            if int(value) < minimum or int(value) % 2 == 0:
                raise ValueError(f"{name} must be an odd integer >= {minimum}")
            object.__setattr__(self, name, int(value))
        if self.fine_eta_points <= self.coarse_eta_points:
            raise ValueError("fine_eta_points must exceed coarse_eta_points")
        if self.fine_outer_y_points <= self.coarse_outer_y_points:
            raise ValueError("fine_outer_y_points must exceed coarse_outer_y_points")

        safety = float(self.safety_fraction)
        drift = float(self.max_refinement_drift_fraction)
        rtol = float(self.rtol)
        atol = float(self.atol)
        max_step = float(self.max_step)
        if not math.isfinite(safety) or not 0.01 <= safety <= 0.50:
            raise ValueError("safety_fraction must lie in [0.01,0.50]")
        if not math.isfinite(drift) or not 0.005 <= drift <= 0.20:
            raise ValueError("max_refinement_drift_fraction must lie in [0.005,0.20]")
        if not math.isfinite(rtol) or not 1.0e-12 <= rtol <= 1.0e-7:
            raise ValueError("rtol must lie in [1e-12,1e-7]")
        if not math.isfinite(atol) or not 1.0e-14 <= atol <= 1.0e-9:
            raise ValueError("atol must lie in [1e-14,1e-9]")
        if not math.isfinite(max_step) or not 0.05 <= max_step <= 0.50:
            raise ValueError("max_step must lie in [0.05,0.50]")
        object.__setattr__(self, "safety_fraction", safety)
        object.__setattr__(self, "max_refinement_drift_fraction", drift)
        object.__setattr__(self, "rtol", rtol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(self, "max_step", max_step)
        object.__setattr__(
            self,
            "_helper",
            KokunoSourceRescaledAppendixBBoundary(reference=self.reference),
        )

    @property
    def X_0(self) -> float:
        return float(self.reference.X0)

    @property
    def y_i(self) -> float:
        return math.log(X_I / self.X_0)

    @property
    def transition_y_end(self) -> float:
        return 2.0 * float(self.reference.log_transition_width)

    def _eta_grid(self, points: int) -> np.ndarray:
        theta = np.linspace(0.0, math.pi, int(points), dtype=float)
        nodes = -np.cos(theta)
        stationary = float(self.reference.seed.phase_stationary_eta)
        return np.unique(np.concatenate((nodes, np.asarray([stationary], dtype=float))))

    def _y_grid(self, outer_points: int) -> np.ndarray:
        transition_points = max(17, (int(outer_points) + 1) // 2)
        first = np.linspace(0.0, self.transition_y_end, transition_points, dtype=float)
        second = np.linspace(
            self.transition_y_end, self.y_i, int(outer_points), dtype=float
        )
        return np.unique(np.concatenate((first, second)))

    def _solve_p1(self, eta: float, y_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        eta = float(eta)
        helper = self._helper
        initial = float(helper._initial_state(eta)[0])
        sl = helper._make_slice(eta)

        def rhs(y: float, state: np.ndarray) -> np.ndarray:
            X = helper._X_from_y(float(y))
            ref = helper._reference_state(X, eta, sl)
            p1 = float(state[0])
            value = X * ref["S_q"] / ref["L"] - ref["l"] * p1
            if not math.isfinite(value):
                raise OverflowError("selected reference p1 derivative became non-finite")
            return np.asarray([value], dtype=float)

        solved = solve_ivp(
            rhs,
            (0.0, self.y_i),
            np.asarray([initial], dtype=float),
            method="DOP853",
            t_eval=np.asarray(y_grid, dtype=float),
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
        )
        if not solved.success or solved.y.shape != (1, y_grid.size):
            raise RuntimeError(f"selected reference p1 solve failed: {solved.message}")
        p1 = np.asarray(solved.y[0], dtype=float)
        if np.any(~np.isfinite(p1)):
            raise OverflowError("selected reference p1 solve produced non-finite values")
        derivatives = np.empty_like(p1)
        for index, (y, value) in enumerate(zip(y_grid, p1, strict=True)):
            derivatives[index] = float(rhs(float(y), np.asarray([value]))[0])
        return p1, derivatives

    def _sample(self, eta_points: int, outer_y_points: int) -> dict[str, Any]:
        eta_grid = self._eta_grid(eta_points)
        y_grid = self._y_grid(outer_y_points)
        max_abs = -math.inf
        max_eta = math.nan
        max_y = math.nan
        max_p1 = math.nan
        min_p1 = math.inf
        min_radial_derivative = math.inf
        max_radial_derivative = -math.inf

        for eta in eta_grid:
            p1, derivatives = self._solve_p1(float(eta), y_grid)
            local_index = int(np.argmax(np.abs(p1)))
            local_abs = float(abs(p1[local_index]))
            if local_abs > max_abs:
                max_abs = local_abs
                max_eta = float(eta)
                max_y = float(y_grid[local_index])
                max_p1 = float(p1[local_index])
            min_p1 = min(min_p1, float(np.min(p1)))
            min_radial_derivative = min(
                min_radial_derivative, float(np.min(derivatives))
            )
            max_radial_derivative = max(
                max_radial_derivative, float(np.max(derivatives))
            )

        max_X = min(X_I, self.X_0 * math.exp(max_y))
        return {
            "eta_grid_size": int(eta_grid.size),
            "radial_grid_size": int(y_grid.size),
            "sample_count": int(eta_grid.size * y_grid.size),
            "sample_max_abs_p1_reference": max_abs,
            "sample_max_p1_reference": max_p1,
            "sample_max_eta": max_eta,
            "sample_max_y": max_y,
            "sample_max_X": max_X,
            "sample_min_p1_reference": min_p1,
            "sample_min_Dy_p1_reference": min_radial_derivative,
            "sample_max_Dy_p1_reference": max_radial_derivative,
        }

    @cached_property
    def envelope(self) -> dict[str, Any]:
        coarse = self._sample(self.coarse_eta_points, self.coarse_outer_y_points)
        fine = self._sample(self.fine_eta_points, self.fine_outer_y_points)
        coarse_max = float(coarse["sample_max_abs_p1_reference"])
        fine_max = float(fine["sample_max_abs_p1_reference"])
        scale = max(1.0, fine_max, coarse_max)
        drift = abs(fine_max - coarse_max) / scale
        candidate = max(fine_max, coarse_max) * (1.0 + self.safety_fraction)
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
            "M0_candidate_upper_envelope": candidate,
            "numerical_envelope_guard_passed": guard,
            "continuum_supremum_proved": False,
        }

    def candidate_B0_diagnostic(self, combined_profile_C0_norm: float) -> float:
        if not self.envelope["numerical_envelope_guard_passed"]:
            raise RuntimeError("selected M0 numerical envelope failed its refinement guard")
        binding = KokunoPA10L0Binding(seed=self.reference.seed)
        guard = KokunoPA10SourceChoiceOrderGuard(X_i=X_I)
        return guard.displayed_B0_upper_bound(
            combined_profile_C0_norm=float(combined_profile_C0_norm),
            M0=float(self.envelope["M0_candidate_upper_envelope"]),
            L0=binding.L_0,
        )

    def report(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": _source_payload(),
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "selected_reference_sha256": self.reference.sha256,
            "geometry": {
                "X_0": self.X_0,
                "X_i": X_I,
                "y_i": self.y_i,
                "reference_transition_y_end": self.transition_y_end,
            },
            "numerical_envelope": copy.deepcopy(self.envelope),
            "resolved_candidate_side_B0_inputs": ["L_0", "selected M_0 numerical envelope"],
            "unresolved_source_B0_dependencies": [
                "||log phi_*+log Phi(4,.)||_{C^0_eta}",
                "M_0 analytic/source-fixed-point bound",
            ],
            "interpretation": (
                "candidate-side scale only; selected contraction-center reference, "
                "not a continuum interval proof and not the source fixed point"
            ),
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": _source_payload(),
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "reference": self.reference.to_payload(),
            "parameters": {
                "coarse_eta_points": self.coarse_eta_points,
                "fine_eta_points": self.fine_eta_points,
                "coarse_outer_y_points": self.coarse_outer_y_points,
                "fine_outer_y_points": self.fine_outer_y_points,
                "safety_fraction": self.safety_fraction,
                "max_refinement_drift_fraction": self.max_refinement_drift_fraction,
                "rtol": self.rtol,
                "atol": self.atol,
                "max_step": self.max_step,
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10SelectedM0Envelope":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected selected M0 envelope schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("selected M0 source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("selected M0 truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("selected M0 payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing selected M0 parameters")
        obj = cls(
            reference=KokunoSourceRescaledReferenceContinuation.from_payload(
                payload.get("reference")
            ),
            **params,
        )
        if obj.to_payload() != payload:
            raise ValueError("selected M0 replay changed payload")
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
    def load_json(cls, path: str | Path) -> "KokunoPA10SelectedM0Envelope":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Materialize the selected-reference PA.10 M0 numerical envelope"
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    envelope = KokunoPA10SelectedM0Envelope()
    report = envelope.report()
    report["envelope_sha256"] = envelope.sha256
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
