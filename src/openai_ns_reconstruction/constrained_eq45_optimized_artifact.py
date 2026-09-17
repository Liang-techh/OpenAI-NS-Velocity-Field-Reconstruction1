"""Reproducible artifact bridge for the bounded optimized Eq45 diagnostic candidate.

This module does not rerun the optimizer.  It consumes the checked CR005
profile/force optimization receipt, replays only the recorded bounded
``Phi(0,2)`` and ``F(0,2)`` changes on the canonical frozen Eq45 seed, and
requires the reconstructed candidate identity to match the optimizer report.

The resulting object remains a research/visualization candidate.  Reproducible
serialization and residual reduction do not establish visual correspondence,
Navier--Stokes validity, paper exactness, or hidden OpenAI-field identity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_force_optimization import candidate_with_profile_pair

OPTIMIZATION_SCHEMA = "eq45_profile_force_optimization_v1"
OPTIMIZATION_TASK_ID = "CR005-EQ45-PROFILE-FORCE-OPT-015"
CANONICAL_SEED_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
OPTIMIZED_CANDIDATE_SHA256 = "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_optimization_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimization.json"


def default_seed_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def default_candidate_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_profile_force_optimized_candidate.json"


def _load_mapping(path: Path, *, label: str) -> Mapping[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    return payload


def _validated_profile_values(report: Mapping[str, Any], seed: Eq45VelocityCandidate) -> tuple[float, float]:
    if report.get("schema") != OPTIMIZATION_SCHEMA:
        raise ValueError("unexpected Eq45 optimization schema")
    if report.get("task_id") != OPTIMIZATION_TASK_ID:
        raise ValueError("unexpected Eq45 optimization task id")
    if report.get("candidate") != "Eq45VelocityCandidate.seed()":
        raise ValueError("optimization receipt is not bound to the canonical Eq45 seed")

    contract = report.get("contract")
    fit = report.get("fit")
    optimized_meta = report.get("optimized_candidate")
    status = report.get("status")
    if not all(isinstance(item, Mapping) for item in (contract, fit, optimized_meta, status)):
        raise ValueError("optimization receipt is missing required object sections")

    if contract.get("profile_parameters") != ["Phi(0,2)", "F(0,2)"]:
        raise ValueError("optimization receipt changed the registered profile parameters")
    if contract.get("training_holdout_separate") is not True:
        raise ValueError("optimization receipt must preserve train/holdout separation")
    if optimized_meta.get("canonical_seed_promoted") is not False:
        raise ValueError("diagnostic optimization must not promote the canonical seed")
    if optimized_meta.get("roundtrip_serializable") is not True:
        raise ValueError("optimizer did not record candidate round-trip serialization")
    if optimized_meta.get("sha256") != OPTIMIZED_CANDIDATE_SHA256:
        raise ValueError("optimization receipt candidate SHA drifted")

    forbidden_promotions = {
        "pde_validated": False,
        "visualization_ready_promoted": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, expected in forbidden_promotions.items():
        if status.get(key) is not expected:
            raise ValueError(f"unsupported scientific promotion in optimization receipt: {key}")

    try:
        phi_02 = float(fit["Phi_02"])
        F_02 = float(fit["F_02"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("optimization receipt is missing finite profile coefficients") from exc
    if not np.isfinite(phi_02) or not np.isfinite(F_02):
        raise ValueError("optimized profile coefficients must be finite")

    bounds = contract.get("profile_bounds")
    expected_limit = float(seed.profile_basis.coefficient_limit)
    if bounds != [-expected_limit, expected_limit]:
        raise ValueError("optimization profile bounds disagree with the candidate contract")
    if not (-expected_limit <= phi_02 <= expected_limit and -expected_limit <= F_02 <= expected_limit):
        raise ValueError("optimized profile coefficients exceed the governed bounds")
    return phi_02, F_02


def candidate_from_optimization_receipt(
    optimization_path: str | Path | None = None,
    seed_path: str | Path | None = None,
) -> Eq45VelocityCandidate:
    """Reconstruct the checked diagnostic candidate without rerunning optimization."""
    optimization = Path(optimization_path) if optimization_path is not None else default_optimization_path()
    seed_file = Path(seed_path) if seed_path is not None else default_seed_path()
    report = _load_mapping(optimization, label="optimization receipt")
    seed = Eq45VelocityCandidate.load_json(seed_file)
    if seed.sha256 != CANONICAL_SEED_SHA256:
        raise ValueError("canonical Eq45 seed SHA drifted")

    phi_02, F_02 = _validated_profile_values(report, seed)
    candidate = candidate_with_profile_pair(seed, phi_02, F_02)
    if candidate.sha256 != OPTIMIZED_CANDIDATE_SHA256:
        raise ValueError("reconstructed optimized candidate SHA does not match the checked receipt")

    truth = candidate.to_dict()["truth_boundary"]
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"optimized diagnostic candidate illegally promotes {key}")
    if truth.get("velocity_export_ready") is not True or truth.get("callable_serializable") is not True:
        raise ValueError("optimized diagnostic candidate lost callable/export readiness")
    return candidate


def load_checked_candidate(
    candidate_path: str | Path | None = None,
    optimization_path: str | Path | None = None,
    seed_path: str | Path | None = None,
) -> Eq45VelocityCandidate:
    """Load the committed artifact and require equality with receipt reconstruction."""
    path = Path(candidate_path) if candidate_path is not None else default_candidate_path()
    saved = Eq45VelocityCandidate.load_json(path)
    rebuilt = candidate_from_optimization_receipt(optimization_path, seed_path)
    if saved.sha256 != rebuilt.sha256 or saved.to_dict() != rebuilt.to_dict():
        raise ValueError("saved optimized candidate does not match receipt reconstruction")
    return saved


def write_checked_candidate(
    path: str | Path,
    optimization_path: str | Path | None = None,
    seed_path: str | Path | None = None,
) -> Eq45VelocityCandidate:
    """Materialize the receipt-reconstructed diagnostic candidate as governed JSON."""
    candidate = candidate_from_optimization_receipt(optimization_path, seed_path)
    candidate.save_json(path)
    roundtrip = Eq45VelocityCandidate.load_json(path)
    if roundtrip.sha256 != candidate.sha256:
        raise RuntimeError("optimized candidate JSON round trip changed candidate identity")
    return roundtrip
