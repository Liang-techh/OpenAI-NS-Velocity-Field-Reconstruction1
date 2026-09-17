"""Diagnose active-bound and local-basin behavior of the bounded Eq45 fit.

This module is deliberately downstream of ``CR005-EQ45-PROFILE-FORCE-OPT-015``.
It does not introduce a new optimizer target, pressure model, or forcing direction.
It reuses the exact registered two-profile objective and the same bounded
``RestrictedForce(a,c)`` variable projection to ask two narrow questions:

1. does the reported solution still improve when either profile coefficient is
   moved to its existing upper bound; and
2. does one deterministic alternate L-BFGS-B start return to the same basin?

Every force used on holdout samples is fitted on the training probes first and
then frozen.  The diagnostic is optimization/representation evidence only;
solver convergence, an active bound, or a lower objective is not Navier--Stokes
validation and is not visual-correspondence evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import minimize

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_force_optimization import (
    DEFAULT_DERIVATIVE_STEPS,
    DEFAULT_FIT_SEED,
    DEFAULT_FIT_STEP,
    DEFAULT_FTOL,
    DEFAULT_GTOL,
    DEFAULT_HOLDOUT_SEED,
    DEFAULT_MAX_FUNCTION_EVALUATIONS,
    DEFAULT_SAMPLE_COUNT,
    PROFILE_BOUNDS,
    _force_design,
    _project_force,
    candidate_with_profile_pair,
)
from .constrained_force import RestrictedForce
from .constrained_pressure_compatibility import curl_pressure_free_momentum
from .constrained_restricted_force_curl_capacity import (
    FORCE_BOUNDS,
    REGISTERED_NU,
    _curl_force,
    _vector_max,
    _vector_rms,
    deterministic_probe_cloud,
)

Array = np.ndarray
RECEIPT_SCHEMA = "eq45_profile_force_optimization_v1"
RECEIPT_TASK = "CR005-EQ45-PROFILE-FORCE-OPT-015"
TASK_ID = "CR005-EQ45-PROFILE-BOUNDARY-DIAGNOSTIC-016"
DEFAULT_ALTERNATE_START = (0.0, 0.0)
DEFAULT_INWARD_STEP = 0.02
DEFAULT_GRADIENT_STEP = 0.002


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_receipt_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimization.json"


@dataclass(frozen=True)
class BoundaryProbe:
    name: str
    phi_02: float
    F_02: float
    training_rms_after_force: float
    training_max_after_force: float
    force_a: float
    force_c: float


@dataclass(frozen=True)
class HoldoutLevel:
    spatial_step: float
    receipt_rms_after_force: float
    probe_rms_after_force: float
    probe_max_after_force: float
    relative_change_vs_receipt: float


@dataclass(frozen=True)
class AlternateRestart:
    start_phi_02: float
    start_F_02: float
    success: bool
    message: str
    function_evaluations: int
    iterations: int
    phi_02: float
    F_02: float
    training_rms_after_force: float
    force_a: float
    force_c: float
    parameter_distance_from_receipt: float
    objective_relative_change_vs_receipt: float


@dataclass(frozen=True)
class Eq45ProfileBoundaryDiagnosticReport:
    task_id: str
    receipt_schema: str
    fit_seed: int
    holdout_seed: int
    sample_count: int
    fit_step: float
    derivative_steps: tuple[float, ...]
    profile_bounds: tuple[float, float]
    force_bounds: tuple[float, float]
    max_function_evaluations: int
    receipt_phi_02: float
    receipt_F_02: float
    receipt_training_rms_after_force: float
    receipt_force_a: float
    receipt_force_c: float
    receipt_candidate_sha256: str
    reconstructed_candidate_sha256: str
    boundary_probes: tuple[BoundaryProbe, ...]
    best_boundary_probe: str
    best_boundary_training_relative_change: float
    local_gradient_phi: float
    local_gradient_F: float
    local_gradient_step: float
    alternate_restart: AlternateRestart
    holdout_levels_for_best_boundary: tuple[HoldoutLevel, ...]
    velocity_changed: bool = False
    pressure_fitted: bool = False
    forcing_family_changed: bool = False
    holdout_force_refit: bool = False
    pde_validated: bool = False
    visualization_ready_promoted: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["boundary_probes"] = [asdict(item) for item in self.boundary_probes]
        payload["holdout_levels_for_best_boundary"] = [
            asdict(item) for item in self.holdout_levels_for_best_boundary
        ]
        payload["interpretation"] = (
            "active-bound/local-basin diagnostic under the unchanged registered "
            "Phi(0,2)/F(0,2) objective and restricted-force variable projection; "
            "this is optimization-capacity evidence only, not PDE or visual validation"
        )
        return payload


def load_and_validate_receipt(path: str | Path | None = None) -> dict[str, Any]:
    receipt_path = default_receipt_path() if path is None else Path(path)
    payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    if payload.get("schema") != RECEIPT_SCHEMA:
        raise ValueError("unexpected Eq45 optimization receipt schema")
    if payload.get("task_id") != RECEIPT_TASK:
        raise ValueError("unexpected Eq45 optimization receipt task id")

    contract = payload.get("contract", {})
    if tuple(float(v) for v in contract.get("profile_bounds", ())) != PROFILE_BOUNDS:
        raise ValueError("profile-bound contract drift")
    if tuple(float(v) for v in contract.get("force_bounds", ())) != FORCE_BOUNDS:
        raise ValueError("force-bound contract drift")
    if int(contract.get("fit_seed", -1)) != DEFAULT_FIT_SEED:
        raise ValueError("fit-seed contract drift")
    if int(contract.get("holdout_seed", -1)) != DEFAULT_HOLDOUT_SEED:
        raise ValueError("holdout-seed contract drift")
    if int(contract.get("sample_count", -1)) != DEFAULT_SAMPLE_COUNT:
        raise ValueError("sample-count contract drift")
    if float(contract.get("fit_step", float("nan"))) != DEFAULT_FIT_STEP:
        raise ValueError("fit-step contract drift")
    if tuple(float(v) for v in contract.get("holdout_derivative_steps", ())) != tuple(
        DEFAULT_DERIVATIVE_STEPS
    ):
        raise ValueError("derivative-step contract drift")
    if int(contract.get("max_function_evaluations", -1)) != DEFAULT_MAX_FUNCTION_EVALUATIONS:
        raise ValueError("optimizer-budget contract drift")
    if contract.get("force_family") != "preregistered_restricted_two_parameter_family":
        raise ValueError("forcing-family contract drift")
    if not bool(contract.get("training_holdout_separate", False)):
        raise ValueError("training/holdout separation must remain explicit")

    fit = payload.get("fit", {})
    optimized = payload.get("optimized_candidate", {})
    for key in ("Phi_02", "F_02", "force_a", "force_c", "curl_rms_after_force"):
        value = float(fit.get(key, float("nan")))
        if not np.isfinite(value):
            raise ValueError(f"nonfinite receipt fit field: {key}")
    if not isinstance(optimized.get("sha256"), str) or len(optimized["sha256"]) != 64:
        raise ValueError("invalid optimized candidate identity")
    if bool(optimized.get("canonical_seed_promoted", True)):
        raise ValueError("diagnostic candidate must not be canonically promoted")
    status = payload.get("status", {})
    forbidden_true = (
        "pde_validated",
        "visualization_ready_promoted",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    )
    if any(bool(status.get(key, False)) for key in forbidden_true):
        raise ValueError("receipt unexpectedly promotes a scientific claim")
    return payload


def _training_evaluation(
    base: Eq45VelocityCandidate,
    phi_02: float,
    F_02: float,
    points: Array,
    times: Array,
    design: Array,
    step: float,
) -> tuple[Eq45VelocityCandidate, float, float, RestrictedForce]:
    trial = candidate_with_profile_pair(base, phi_02, F_02)
    _, target = curl_pressure_free_momentum(
        trial,
        points,
        times,
        nu=REGISTERED_NU,
        spatial_step=step,
        time_step=step,
    )
    force_solution, after = _project_force(target, design)
    force = RestrictedForce(a=float(force_solution.x[0]), c=float(force_solution.x[1]))
    return trial, _vector_rms(after), _vector_max(after), force


def _validate_steps(steps: Iterable[float]) -> tuple[float, ...]:
    result = tuple(float(step) for step in steps)
    if len(result) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in result):
        raise ValueError("derivative_steps must contain at least three positive finite levels")
    if any(b >= a for a, b in zip(result, result[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")
    return result


def audit_eq45_profile_boundary_pressure(
    *,
    receipt_path: str | Path | None = None,
    alternate_start: tuple[float, float] = DEFAULT_ALTERNATE_START,
    inward_step: float = DEFAULT_INWARD_STEP,
    gradient_step: float = DEFAULT_GRADIENT_STEP,
    derivative_steps: Iterable[float] = DEFAULT_DERIVATIVE_STEPS,
) -> Eq45ProfileBoundaryDiagnosticReport:
    """Replay the #104 optimum and diagnose bound pressure/local-basin sensitivity."""
    receipt = load_and_validate_receipt(receipt_path)
    contract = receipt["contract"]
    fit = receipt["fit"]
    expected_sha = receipt["optimized_candidate"]["sha256"]
    steps = _validate_steps(derivative_steps)

    if not np.isfinite(inward_step) or inward_step <= 0.0:
        raise ValueError("inward_step must be positive and finite")
    if not np.isfinite(gradient_step) or gradient_step <= 0.0:
        raise ValueError("gradient_step must be positive and finite")
    alt = np.asarray(alternate_start, dtype=float)
    if alt.shape != (2,) or not np.all(np.isfinite(alt)):
        raise ValueError("alternate_start must contain two finite values")
    if np.any(alt < PROFILE_BOUNDS[0]) or np.any(alt > PROFILE_BOUNDS[1]):
        raise ValueError("alternate_start lies outside registered profile bounds")

    base = Eq45VelocityCandidate.seed()
    phi_star = float(fit["Phi_02"])
    F_star = float(fit["F_02"])
    reconstructed = candidate_with_profile_pair(base, phi_star, F_star)
    if reconstructed.sha256 != expected_sha:
        raise ValueError("optimized candidate identity no longer matches #104 receipt")

    fit_points, fit_times = deterministic_probe_cloud(
        int(contract["fit_seed"]), int(contract["sample_count"])
    )
    holdout_points, holdout_times = deterministic_probe_cloud(
        int(contract["holdout_seed"]), int(contract["sample_count"])
    )
    fit_step = float(contract["fit_step"])
    design = _force_design(fit_points, fit_times, fit_step)

    cache: dict[tuple[float, float], tuple[Eq45VelocityCandidate, float, float, RestrictedForce]] = {}

    def evaluate(phi: float, swirl: float):
        key = (float(phi), float(swirl))
        if key not in cache:
            cache[key] = _training_evaluation(
                base, key[0], key[1], fit_points, fit_times, design, fit_step
            )
        return cache[key]

    receipt_trial, receipt_rms, _, receipt_force = evaluate(phi_star, F_star)
    recorded_rms = float(fit["curl_rms_after_force"])
    if not np.isclose(receipt_rms, recorded_rms, rtol=5e-10, atol=5e-10):
        raise ValueError("replayed training objective no longer matches #104 receipt")
    recorded_force = np.array([float(fit["force_a"]), float(fit["force_c"])])
    replay_force = np.array([receipt_force.a, receipt_force.c])
    if not np.allclose(replay_force, recorded_force, rtol=5e-9, atol=5e-9):
        raise ValueError("replayed restricted-force projection no longer matches #104 receipt")

    upper = PROFILE_BOUNDS[1]
    probes_raw = (
        ("receipt", phi_star, F_star),
        ("phi_upper", upper, F_star),
        ("F_upper", phi_star, upper),
        ("upper_corner", upper, upper),
        ("phi_inward", max(PROFILE_BOUNDS[0], phi_star - inward_step), F_star),
        ("F_inward", phi_star, max(PROFILE_BOUNDS[0], F_star - inward_step)),
    )
    probes: list[BoundaryProbe] = []
    for name, phi_value, F_value in probes_raw:
        _, rms, maximum, force = evaluate(phi_value, F_value)
        probes.append(
            BoundaryProbe(
                name=name,
                phi_02=float(phi_value),
                F_02=float(F_value),
                training_rms_after_force=float(rms),
                training_max_after_force=float(maximum),
                force_a=float(force.a),
                force_c=float(force.c),
            )
        )

    boundary_names = {"phi_upper", "F_upper", "upper_corner"}
    best_boundary = min(
        (probe for probe in probes if probe.name in boundary_names),
        key=lambda probe: probe.training_rms_after_force,
    )

    gstep = float(gradient_step)
    phi_step = min(gstep, 0.45 * (upper - phi_star), 0.45 * (phi_star - PROFILE_BOUNDS[0]))
    F_step = min(gstep, 0.45 * (upper - F_star), 0.45 * (F_star - PROFILE_BOUNDS[0]))
    if phi_step <= 0.0 or F_step <= 0.0:
        raise ValueError("receipt optimum does not admit a finite-difference gradient probe")
    phi_plus = evaluate(phi_star + phi_step, F_star)[1]
    phi_minus = evaluate(phi_star - phi_step, F_star)[1]
    F_plus = evaluate(phi_star, F_star + F_step)[1]
    F_minus = evaluate(phi_star, F_star - F_step)[1]
    grad_phi = float((phi_plus - phi_minus) / (2.0 * phi_step))
    grad_F = float((F_plus - F_minus) / (2.0 * F_step))

    def objective(parameters: Array) -> float:
        return float(evaluate(float(parameters[0]), float(parameters[1]))[1])

    restart = minimize(
        objective,
        alt,
        method="L-BFGS-B",
        bounds=(PROFILE_BOUNDS, PROFILE_BOUNDS),
        options={
            "maxfun": int(contract["max_function_evaluations"]),
            "maxiter": int(contract["max_function_evaluations"]),
            "ftol": DEFAULT_FTOL,
            "gtol": DEFAULT_GTOL,
            "maxls": 20,
        },
    )
    if not np.all(np.isfinite(restart.x)) or not np.isfinite(restart.fun):
        raise RuntimeError("alternate restart returned nonfinite state")
    _, restart_rms, _, restart_force = evaluate(float(restart.x[0]), float(restart.x[1]))
    restart_report = AlternateRestart(
        start_phi_02=float(alt[0]),
        start_F_02=float(alt[1]),
        success=bool(restart.success),
        message=str(restart.message),
        function_evaluations=int(restart.nfev),
        iterations=int(restart.nit),
        phi_02=float(restart.x[0]),
        F_02=float(restart.x[1]),
        training_rms_after_force=float(restart_rms),
        force_a=float(restart_force.a),
        force_c=float(restart_force.c),
        parameter_distance_from_receipt=float(
            np.linalg.norm(np.asarray(restart.x, dtype=float) - np.array([phi_star, F_star]))
        ),
        objective_relative_change_vs_receipt=float(
            (restart_rms - receipt_rms) / max(receipt_rms, 1e-15)
        ),
    )

    best_candidate, _, _, best_force = evaluate(
        best_boundary.phi_02, best_boundary.F_02
    )
    holdout_rows: list[HoldoutLevel] = []
    for step in steps:
        receipt_target = curl_pressure_free_momentum(
            receipt_trial,
            holdout_points,
            holdout_times,
            nu=REGISTERED_NU,
            spatial_step=step,
            time_step=step,
        )[1]
        best_target = curl_pressure_free_momentum(
            best_candidate,
            holdout_points,
            holdout_times,
            nu=REGISTERED_NU,
            spatial_step=step,
            time_step=step,
        )[1]
        receipt_after = receipt_target - _curl_force(
            receipt_force, holdout_points, holdout_times, step
        )
        best_after = best_target - _curl_force(
            best_force, holdout_points, holdout_times, step
        )
        receipt_holdout_rms = _vector_rms(receipt_after)
        best_holdout_rms = _vector_rms(best_after)
        holdout_rows.append(
            HoldoutLevel(
                spatial_step=float(step),
                receipt_rms_after_force=float(receipt_holdout_rms),
                probe_rms_after_force=float(best_holdout_rms),
                probe_max_after_force=float(_vector_max(best_after)),
                relative_change_vs_receipt=float(
                    (best_holdout_rms - receipt_holdout_rms)
                    / max(receipt_holdout_rms, 1e-15)
                ),
            )
        )

    return Eq45ProfileBoundaryDiagnosticReport(
        task_id=TASK_ID,
        receipt_schema=RECEIPT_SCHEMA,
        fit_seed=int(contract["fit_seed"]),
        holdout_seed=int(contract["holdout_seed"]),
        sample_count=int(contract["sample_count"]),
        fit_step=fit_step,
        derivative_steps=steps,
        profile_bounds=PROFILE_BOUNDS,
        force_bounds=FORCE_BOUNDS,
        max_function_evaluations=int(contract["max_function_evaluations"]),
        receipt_phi_02=phi_star,
        receipt_F_02=F_star,
        receipt_training_rms_after_force=float(receipt_rms),
        receipt_force_a=float(receipt_force.a),
        receipt_force_c=float(receipt_force.c),
        receipt_candidate_sha256=expected_sha,
        reconstructed_candidate_sha256=reconstructed.sha256,
        boundary_probes=tuple(probes),
        best_boundary_probe=best_boundary.name,
        best_boundary_training_relative_change=float(
            (best_boundary.training_rms_after_force - receipt_rms) / max(receipt_rms, 1e-15)
        ),
        local_gradient_phi=grad_phi,
        local_gradient_F=grad_F,
        local_gradient_step=gstep,
        alternate_restart=restart_report,
        holdout_levels_for_best_boundary=tuple(holdout_rows),
    )


def main() -> int:
    report = audit_eq45_profile_boundary_pressure()
    print(json.dumps(report.to_dict(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
