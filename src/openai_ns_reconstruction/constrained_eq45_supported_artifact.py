"""Reproducible artifact bridge for the support-connected Eq45 child candidate.

This module does not modify the Eq. (4.5) profiles or rerun any optimizer.  It
reconstructs the production physical-support child from the canonical frozen
Eq45 parent and the default governed ``AxisymmetricPhysicalTaper``, then binds
that child to one committed JSON identity for downstream validation and
visualization consumers.

A loadable/exportable support-connected field is not, by itself, evidence of
physical-support validation, visual correspondence, Navier--Stokes validity,
paper exactness, hidden OpenAI-field identity, or blow-up.
"""
from __future__ import annotations

from pathlib import Path

from .constrained_axisymmetric_physical_taper import AxisymmetricPhysicalTaper
from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


CANONICAL_PARENT_SHA256 = (
    "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
)
SUPPORTED_CANDIDATE_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_parent_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def default_candidate_path() -> Path:
    return _repo_root() / "artifacts" / "constrained" / "eq45_supported_velocity_candidate.json"


def build_supported_candidate(
    parent_path: str | Path | None = None,
) -> Eq45SupportedVelocityCandidate:
    """Reconstruct the governed default support child from the canonical parent."""
    path = Path(parent_path) if parent_path is not None else default_parent_path()
    parent = Eq45VelocityCandidate.load_json(path)
    if parent.sha256 != CANONICAL_PARENT_SHA256:
        raise ValueError("canonical Eq45 parent SHA drifted")

    candidate = Eq45SupportedVelocityCandidate(
        parent=parent,
        taper=AxisymmetricPhysicalTaper(),
    )
    if candidate.sha256 != SUPPORTED_CANDIDATE_SHA256:
        raise ValueError("support-connected Eq45 candidate SHA drifted")

    truth = candidate.to_dict()["truth_boundary"]
    if truth.get("callable_serializable") is not True:
        raise ValueError("supported candidate lost callable/serialization readiness")
    if truth.get("velocity_export_ready") is not True:
        raise ValueError("supported candidate lost velocity export readiness")
    if truth.get("physical_support_connection_implemented") is not True:
        raise ValueError("supported candidate lost support-connection construction state")
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
            raise ValueError(f"supported candidate illegally promotes {key}")
    return candidate


def load_checked_supported_candidate(
    candidate_path: str | Path | None = None,
    parent_path: str | Path | None = None,
) -> Eq45SupportedVelocityCandidate:
    """Load the committed child and require exact equality with reconstruction."""
    path = Path(candidate_path) if candidate_path is not None else default_candidate_path()
    saved = Eq45SupportedVelocityCandidate.load_json(path)
    rebuilt = build_supported_candidate(parent_path)
    if saved.sha256 != SUPPORTED_CANDIDATE_SHA256:
        raise ValueError("saved supported candidate SHA drifted")
    if saved.to_dict() != rebuilt.to_dict():
        raise ValueError("saved supported candidate does not match governed reconstruction")
    return saved


def write_checked_supported_candidate(
    path: str | Path,
    parent_path: str | Path | None = None,
) -> Eq45SupportedVelocityCandidate:
    """Write the governed support child and verify its exact round-trip identity."""
    candidate = build_supported_candidate(parent_path)
    candidate.save_json(path)
    roundtrip = Eq45SupportedVelocityCandidate.load_json(path)
    if roundtrip.sha256 != SUPPORTED_CANDIDATE_SHA256:
        raise RuntimeError("supported candidate JSON round trip changed candidate identity")
    if roundtrip.to_dict() != candidate.to_dict():
        raise RuntimeError("supported candidate JSON round trip changed candidate payload")
    return roundtrip
