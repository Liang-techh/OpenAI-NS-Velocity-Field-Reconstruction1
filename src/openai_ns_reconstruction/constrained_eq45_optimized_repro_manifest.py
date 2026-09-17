"""Fail-closed reproducibility manifest for the optimized Eq45 research candidate.

The manifest binds one exact saved velocity artifact to the bounded CR005
optimization receipt and to the separately preregistered CR001 project
contract.  It deliberately keeps optimizer/holdout evidence separate from the
project validation seed and from visualization/PDE acceptance states.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_optimized_artifact import (
    CANONICAL_SEED_SHA256,
    OPTIMIZED_CANDIDATE_SHA256,
    load_checked_candidate,
)

MANIFEST_SCHEMA = "eq45_optimized_repro_manifest_v1"
TASK_ID = "CR011-EQ45-OPTIMIZED-REPRO-MANIFEST-009"
OPTIMIZATION_SCHEMA = "eq45_profile_force_optimization_v1"
OPTIMIZATION_TASK_ID = "CR005-EQ45-PROFILE-FORCE-OPT-015"
PROJECT_EXPERIMENT_ID = "compact_axisymmetric_window_v1"
MATERIALIZATION_COMMIT = "13f8ddb2a28b6edc80c06f9e1608148ff80b5ed3"
OPTIMIZER_RECEIPT_COMMIT = "40afef87a00490c79726a46b8a5bc3d1acf9db91"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_manifest_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimized_repro_manifest.json"


def default_constraints_path() -> Path:
    return _repo_root() / "configs" / "constraints.json"


def default_optimization_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimization.json"


def default_candidate_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimized_candidate.json"


def default_seed_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _load_mapping(path: Path, *, label: str) -> Mapping[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    return payload


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    return value


def _require_equal(actual: Any, expected: Any, *, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label} drifted: expected {expected!r}, got {actual!r}")


def _validate_manifest_identity(manifest: Mapping[str, Any]) -> None:
    _require_equal(manifest.get("schema"), MANIFEST_SCHEMA, label="manifest schema")
    _require_equal(manifest.get("task_id"), TASK_ID, label="manifest task id")

    candidate = _mapping(manifest.get("candidate"), label="candidate section")
    _require_equal(candidate.get("sha256"), OPTIMIZED_CANDIDATE_SHA256, label="optimized candidate SHA")
    _require_equal(candidate.get("source_seed_sha256"), CANONICAL_SEED_SHA256, label="canonical seed SHA")
    _require_equal(candidate.get("schema"), "eq45_velocity_candidate_v1", label="candidate schema")
    _require_equal(candidate.get("canonical_seed_promoted"), False, label="canonical seed promotion")

    replay = _mapping(manifest.get("replay"), label="replay section")
    _require_equal(
        replay.get("artifact_materialization_commit"),
        MATERIALIZATION_COMMIT,
        label="artifact materialization commit",
    )
    _require_equal(
        replay.get("optimizer_receipt_commit"),
        OPTIMIZER_RECEIPT_COMMIT,
        label="optimizer receipt commit",
    )
    _require_equal(
        replay.get("rebuild_entrypoint"),
        "openai_ns_reconstruction.constrained_eq45_optimized_artifact:candidate_from_optimization_receipt",
        label="rebuild entrypoint",
    )
    _require_equal(
        replay.get("load_entrypoint"),
        "openai_ns_reconstruction.constrained_eq45_optimized_artifact:load_checked_candidate",
        label="load entrypoint",
    )
    _require_equal(
        replay.get("public_velocity_entrypoints"),
        ["Eq45VelocityCandidate.at_points", "Eq45VelocityCandidate.velocity_xyz"],
        label="public velocity entrypoints",
    )


def _validate_optimization_manifest(manifest: Mapping[str, Any], receipt: Mapping[str, Any]) -> None:
    recorded = _mapping(manifest.get("optimization_receipt"), label="optimization manifest section")
    contract = _mapping(receipt.get("contract"), label="optimization contract")
    fit = _mapping(receipt.get("fit"), label="optimization fit")
    optimized = _mapping(receipt.get("optimized_candidate"), label="optimized candidate receipt")

    _require_equal(receipt.get("schema"), OPTIMIZATION_SCHEMA, label="optimization schema")
    _require_equal(receipt.get("task_id"), OPTIMIZATION_TASK_ID, label="optimization task id")
    _require_equal(recorded.get("schema"), receipt.get("schema"), label="manifest optimization schema")
    _require_equal(recorded.get("task_id"), receipt.get("task_id"), label="manifest optimization task id")

    field_map = {
        "fit_seed": "fit_seed",
        "diagnostic_holdout_seed": "holdout_seed",
        "sample_count": "sample_count",
        "probe_box": "probe_box",
        "probe_time_range": "probe_time_range",
        "fit_step": "fit_step",
        "holdout_derivative_steps": "holdout_derivative_steps",
        "profile_parameters": "profile_parameters",
        "profile_bounds": "profile_bounds",
        "force_family": "force_family",
        "force_bounds": "force_bounds",
        "maximum_function_evaluations": "max_function_evaluations",
    }
    for manifest_key, receipt_key in field_map.items():
        _require_equal(
            recorded.get(manifest_key),
            contract.get(receipt_key),
            label=f"optimization {manifest_key}",
        )

    fit_values = _mapping(recorded.get("fit_values"), label="optimization fit values")
    _require_equal(fit_values.get("Phi(0,2)"), fit.get("Phi_02"), label="Phi(0,2) fit value")
    _require_equal(fit_values.get("F(0,2)"), fit.get("F_02"), label="F(0,2) fit value")
    _require_equal(fit_values.get("force_a"), fit.get("force_a"), label="force a fit value")
    _require_equal(fit_values.get("force_c"), fit.get("force_c"), label="force c fit value")
    _require_equal(optimized.get("sha256"), OPTIMIZED_CANDIDATE_SHA256, label="receipt optimized SHA")
    _require_equal(optimized.get("canonical_seed_promoted"), False, label="receipt canonical promotion")

    results = _mapping(recorded.get("recorded_results"), label="recorded optimization results")
    _require_equal(
        results.get("training_curl_rms_after_force"),
        fit.get("curl_rms_after_force"),
        label="training curl RMS",
    )
    holdout = receipt.get("holdout")
    if not isinstance(holdout, list) or not holdout:
        raise ValueError("optimization receipt holdout must be a nonempty list")
    finest_step = results.get("finest_holdout_step")
    finest_rows = [row for row in holdout if isinstance(row, Mapping) and row.get("step") == finest_step]
    if len(finest_rows) != 1:
        raise ValueError("manifest finest holdout step must identify exactly one receipt row")
    finest = finest_rows[0]
    _require_equal(
        results.get("finest_holdout_curl_rms_after_force"),
        finest.get("optimized_curl_rms_after_force"),
        label="finest holdout curl RMS",
    )
    _require_equal(
        results.get("finest_holdout_curl_max_after_force"),
        finest.get("optimized_curl_max_after_force"),
        label="finest holdout curl max",
    )
    velocity_change = _mapping(receipt.get("velocity_change_on_holdout_probes"), label="velocity change receipt")
    _require_equal(
        results.get("relative_public_velocity_delta_rms"),
        velocity_change.get("relative_delta_rms"),
        label="relative public velocity delta RMS",
    )


def _validate_project_manifest(manifest: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    recorded = _mapping(manifest.get("project_contract"), label="project contract manifest section")
    domain = _mapping(constraints.get("domain"), label="project domain")
    validation = _mapping(constraints.get("validation"), label="project validation")

    _require_equal(recorded.get("schema_version"), constraints.get("schema_version"), label="constraints schema version")
    _require_equal(recorded.get("experiment_id"), constraints.get("experiment_id"), label="experiment id")
    _require_equal(recorded.get("status"), constraints.get("status"), label="project status")
    _require_equal(recorded.get("nu"), constraints.get("nu"), label="project viscosity")
    _require_equal(recorded.get("physical_domain"), domain.get("physical"), label="physical domain")
    _require_equal(recorded.get("evaluation_box"), domain.get("evaluation_box"), label="evaluation box")
    _require_equal(recorded.get("physical_support"), domain.get("support"), label="physical support contract")
    _require_equal(recorded.get("time_interval"), domain.get("time_interval"), label="project time interval")
    _require_equal(recorded.get("project_validation_seed"), validation.get("seed"), label="project validation seed")
    _require_equal(
        recorded.get("validation_derivative_steps"),
        validation.get("derivative_steps"),
        label="project validation derivative steps",
    )
    _require_equal(
        recorded.get("validation_thresholds"),
        validation.get("thresholds"),
        label="project validation thresholds",
    )
    _require_equal(recorded.get("seed_roles_distinct"), True, label="seed-role separation flag")

    optimization = _mapping(manifest.get("optimization_receipt"), label="optimization manifest section")
    project_seed = recorded.get("project_validation_seed")
    fit_seed = optimization.get("fit_seed")
    holdout_seed = optimization.get("diagnostic_holdout_seed")
    if len({project_seed, fit_seed, holdout_seed}) != 3:
        raise ValueError("project validation, optimization fit, and diagnostic holdout seeds must remain distinct")

    _require_equal(recorded.get("experiment_id"), PROJECT_EXPERIMENT_ID, label="governed project experiment")


def _validate_delivery_state(manifest: Mapping[str, Any], candidate_truth: Mapping[str, Any]) -> None:
    state = _mapping(manifest.get("delivery_state"), label="delivery state")
    _require_equal(state.get("velocity_export_ready"), True, label="velocity export readiness")
    _require_equal(state.get("callable_serializable"), True, label="callable/serializable readiness")
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(state.get(key), False, label=f"manifest {key}")
        _require_equal(candidate_truth.get(key), False, label=f"candidate {key}")
    _require_equal(candidate_truth.get("velocity_export_ready"), True, label="candidate velocity export readiness")
    _require_equal(candidate_truth.get("callable_serializable"), True, label="candidate callable readiness")


def audit_optimized_repro_manifest(
    manifest_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    optimization_path: str | Path | None = None,
    candidate_path: str | Path | None = None,
    seed_path: str | Path | None = None,
) -> dict[str, Any]:
    """Audit the exact optimized research artifact and its provenance contract."""
    manifest_file = Path(manifest_path) if manifest_path is not None else default_manifest_path()
    constraints_file = Path(constraints_path) if constraints_path is not None else default_constraints_path()
    optimization_file = Path(optimization_path) if optimization_path is not None else default_optimization_path()
    candidate_file = Path(candidate_path) if candidate_path is not None else default_candidate_path()
    seed_file = Path(seed_path) if seed_path is not None else default_seed_path()

    manifest = _load_mapping(manifest_file, label="reproducibility manifest")
    constraints = _load_mapping(constraints_file, label="project constraints")
    optimization = _load_mapping(optimization_file, label="optimization receipt")
    _validate_manifest_identity(manifest)
    _validate_optimization_manifest(manifest, optimization)
    _validate_project_manifest(manifest, constraints)

    candidate = load_checked_candidate(candidate_file, optimization_file, seed_file)
    _require_equal(candidate.sha256, OPTIMIZED_CANDIDATE_SHA256, label="loaded optimized candidate SHA")
    candidate_payload = candidate.to_dict()
    _require_equal(candidate_payload.get("schema"), "eq45_velocity_candidate_v1", label="loaded candidate schema")
    _validate_delivery_state(manifest, _mapping(candidate_payload.get("truth_boundary"), label="candidate truth boundary"))

    replay = _mapping(manifest.get("replay"), label="replay section")
    probe = _mapping(replay.get("public_probe"), label="public probe")
    points = np.asarray(probe.get("points"), dtype=float)
    times = np.asarray(probe.get("times"), dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or times.shape != (points.shape[0],):
        raise ValueError("public probe must provide matching (n,3) points and (n,) times")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(times)):
        raise ValueError("public probe coordinates/times must be finite")
    velocity = candidate.at_points(points, times)
    if velocity.shape != (points.shape[0], 3) or not np.all(np.isfinite(velocity)):
        raise ValueError("public velocity evaluator returned malformed/nonfinite output")
    velocity_rms = float(np.sqrt(np.mean(np.sum(velocity * velocity, axis=1))))
    if probe.get("require_finite") is not True or probe.get("require_nonzero_rms") is not True:
        raise ValueError("public probe fail-closed requirements drifted")
    if not np.isfinite(velocity_rms) or velocity_rms <= 0.0:
        raise ValueError("optimized candidate public velocity probe is inactive")

    reference = _mapping(manifest.get("reference_visualization"), label="reference visualization section")
    _require_equal(reference.get("times"), [0.25, 0.5, 0.75], label="reference visualization times")
    _require_equal(reference.get("classification"), "autonomous_design", label="reference time classification")
    _require_equal(
        reference.get("public_openai_frame_time_mapping"),
        "pending_unknown",
        label="public frame-time mapping",
    )

    scope = _mapping(manifest.get("scope"), label="scope section")
    _require_equal(scope.get("runtime_environment_receipt"), "separate_cr011_lane_not_duplicated_here", label="runtime scope")
    _require_equal(scope.get("canonical_promotion"), "pending_independent_visual_evidence", label="canonical promotion scope")
    _require_equal(scope.get("claim"), "reproducible_bounded_research_candidate_only", label="manifest claim scope")

    project = _mapping(manifest.get("project_contract"), label="project contract manifest section")
    opt = _mapping(manifest.get("optimization_receipt"), label="optimization manifest section")
    return {
        "schema": MANIFEST_SCHEMA,
        "candidate_sha256": candidate.sha256,
        "velocity_probe_rms": velocity_rms,
        "fit_seed": opt["fit_seed"],
        "diagnostic_holdout_seed": opt["diagnostic_holdout_seed"],
        "project_validation_seed": project["project_validation_seed"],
        "velocity_export_ready": True,
        "visualization_ready": False,
        "pde_validated": False,
    }
