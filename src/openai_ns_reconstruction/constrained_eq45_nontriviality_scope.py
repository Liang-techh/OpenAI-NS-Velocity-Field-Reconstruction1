"""Fail-closed nontriviality-scope governance for the Eq. (4.5) candidate.

The serializable Eq45 candidate has a small coefficient-norm guard and can be
sampled as a nonzero public velocity.  Those are candidate-construction facts,
not the CR001 physical-PDE kinetic-energy normalization/anti-collapse contract.
This module keeps the two scopes separate without making PDE acceptance a gate
for saving, loading, evaluating, or exporting a research/visualization field.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate, MIN_COEFFICIENT_NORM


_REQUIRED_FALSE_EQ45_CLAIMS = {
    "physical_support_validated",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _same_number(actual: Any, expected: Any) -> bool:
    try:
        return float(actual) == float(expected)
    except (TypeError, ValueError):
        return False


def audit_eq45_nontriviality_scope(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit candidate-local nonzero evidence versus CR001 energy acceptance.

    A pass means the repository does not confuse the Eq45 constructor/activity
    guard with the separately preregistered CR001 kinetic-energy requirements.
    It is not evidence that the Eq45 candidate satisfies those energy gates or
    the Navier--Stokes equations.
    """

    root = Path(repo_root)
    _require(
        contract.get("schema") == "eq45_nontriviality_scope_contract_v1",
        "unsupported Eq45 nontriviality scope schema",
    )

    eq45 = contract.get("eq45_candidate")
    cr001 = contract.get("cr001_nontriviality")
    semantics = contract.get("state_semantics")
    _require(isinstance(eq45, Mapping), "missing eq45_candidate")
    _require(isinstance(cr001, Mapping), "missing cr001_nontriviality")
    _require(isinstance(semantics, Mapping), "missing state_semantics")

    artifact_path = root / str(eq45.get("artifact", ""))
    constraints_path = root / str(cr001.get("constraints", ""))
    _require(artifact_path.is_file(), "Eq45 candidate artifact is missing")
    _require(constraints_path.is_file(), "CR001 constraints file is missing")

    artifact = _read_json(artifact_path)
    constraints = _read_json(constraints_path)
    _require(
        artifact.get("schema") == "eq45_velocity_candidate_v1",
        "unexpected Eq45 candidate schema",
    )

    _require(
        eq45.get("evaluator")
        == "openai_ns_reconstruction.constrained_eq45_candidate:Eq45VelocityCandidate.at_points",
        "unexpected Eq45 public evaluator",
    )
    _require(
        eq45.get("coefficient_guard_kind") == "profile_coefficient_l2_norm",
        "unexpected Eq45 coefficient guard kind",
    )
    _require(
        _same_number(eq45.get("coefficient_guard_min_exclusive"), MIN_COEFFICIENT_NORM),
        "Eq45 coefficient guard threshold drifted from implementation",
    )
    _require(
        eq45.get("coefficient_guard_classification") == "autonomous_design",
        "Eq45 coefficient guard must remain autonomous_design",
    )
    _require(
        eq45.get("activity_probe_classification") == "autonomous_design",
        "Eq45 activity probe must remain autonomous_design",
    )

    truth = artifact.get("truth_boundary")
    _require(isinstance(truth, Mapping), "Eq45 candidate lacks truth_boundary")
    _require(truth.get("callable_serializable") is True, "Eq45 artifact is not callable/serializable")
    _require(truth.get("velocity_export_ready") is True, "Eq45 artifact is not candidate-locally export ready")

    required_false = contract.get("required_false_eq45_claims")
    _require(isinstance(required_false, list), "required_false_eq45_claims must be a list")
    _require(
        set(required_false) == _REQUIRED_FALSE_EQ45_CLAIMS,
        "required false Eq45 claim set drifted",
    )
    for key in _REQUIRED_FALSE_EQ45_CLAIMS:
        _require(truth.get(key) is False, f"unsupported Eq45 claim promotion: {key}")

    candidate = Eq45VelocityCandidate.load_json(artifact_path)
    coefficients = np.asarray(
        candidate.profile_basis.phi_coefficients
        + candidate.profile_basis.swirl_coefficients,
        dtype=float,
    )
    coefficient_norm = float(np.linalg.norm(coefficients))
    _require(
        coefficient_norm > MIN_COEFFICIENT_NORM,
        "Eq45 candidate failed its constructor nontriviality guard",
    )

    probe_points = np.asarray(eq45.get("activity_probe_points"), dtype=float)
    _require(
        probe_points.ndim == 2 and probe_points.shape[1] == 3 and probe_points.shape[0] > 0,
        "activity_probe_points must have shape (n,3)",
    )
    _require(np.all(np.isfinite(probe_points)), "activity_probe_points must be finite")
    probe_time = eq45.get("activity_probe_time")
    _require(_same_number(probe_time, probe_time), "activity_probe_time must be numeric and finite")
    probe_time = float(probe_time)
    _require(np.isfinite(probe_time), "activity_probe_time must be numeric and finite")
    _require(
        candidate.time_start <= probe_time <= candidate.time_end,
        "activity_probe_time lies outside the candidate interval",
    )
    velocity = np.asarray(candidate.at_points(probe_points, probe_time), dtype=float)
    _require(velocity.shape == probe_points.shape, "Eq45 public velocity returned malformed shape")
    _require(np.all(np.isfinite(velocity)), "Eq45 public velocity returned nonfinite values")
    speeds = np.linalg.norm(velocity, axis=1)
    sampled_speed_max = float(np.max(speeds))
    sampled_speed_rms = float(np.sqrt(np.mean(speeds * speeds)))
    _require(
        sampled_speed_max > 0.0 and sampled_speed_rms > 0.0,
        "Eq45 activity probe is numerically zero",
    )

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, Mapping), "CR001 constraints lack nontriviality")
    expected_pairs = {
        "reference_time": "reference_time",
        "reference_energy": "reference_energy",
        "reference_energy_abs_tolerance": "reference_energy_abs_tolerance",
        "minimum_energy_each_validation_time": "minimum_energy_each_validation_time",
        "maximum_energy_each_validation_time": "maximum_energy_each_validation_time",
    }
    for contract_key, constraints_key in expected_pairs.items():
        _require(
            _same_number(cr001.get(contract_key), nontriviality.get(constraints_key)),
            f"CR001 nontriviality contract drifted: {contract_key}",
        )
    _require(
        cr001.get("classification") == "autonomous_design",
        "CR001 nontriviality values must remain autonomous_design",
    )
    _require(
        cr001.get("scope") == "physical_pde_experiment_acceptance",
        "CR001 nontriviality scope drifted",
    )

    _require(
        semantics.get("eq45_coefficient_guard_is_cr001_energy_evidence") is False,
        "Eq45 coefficient guard cannot become CR001 energy evidence",
    )
    _require(
        semantics.get("eq45_activity_probe_is_cr001_energy_evidence") is False,
        "Eq45 activity probe cannot become CR001 energy evidence",
    )
    _require(
        semantics.get("eq45_cr001_energy_status") == "pending_unknown",
        "Eq45 CR001 energy status must remain pending_unknown until independently measured",
    )
    _require(
        eq45.get("candidate_local_export_requires_cr001_energy_normalization") is False,
        "CR001 energy normalization must not become a candidate-local export gate",
    )
    _require(
        semantics.get("unassessed_cr001_energy_blocks_candidate_local_export") is False,
        "unassessed CR001 energy must not block candidate-local export",
    )
    _require(
        semantics.get("cr001_energy_pass_would_by_itself_promote_pde_validated") is False,
        "energy acceptance alone cannot promote PDE validity",
    )

    return {
        "contract_pass": True,
        "eq45_coefficient_norm": coefficient_norm,
        "eq45_coefficient_guard_min_exclusive": float(MIN_COEFFICIENT_NORM),
        "eq45_activity_probe_time": probe_time,
        "eq45_activity_speed_max": sampled_speed_max,
        "eq45_activity_speed_rms": sampled_speed_rms,
        "eq45_cr001_energy_status": "pending_unknown",
        "candidate_local_export_blocked_by_energy_status": False,
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_eq45_nontriviality_scope_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_eq45_nontriviality_scope(_read_json(Path(path)), repo_root=repo_root)
