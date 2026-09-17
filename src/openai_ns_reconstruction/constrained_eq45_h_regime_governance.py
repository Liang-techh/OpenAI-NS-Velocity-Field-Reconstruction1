"""Fail-closed governance for the distinct Eq. (4.5) ``h`` scopes.

The public coordinate formula is meaningful on a broader ``0 < h < 1/2``
interval than the narrower theorem-specialized ``h < 0.01`` regime recorded in
the source contract.  A numerical candidate chooses its own finite ``h``.  This
module keeps those facts separate so callable/export readiness cannot be
mistaken for theorem identity, PDE validation, or visual correspondence.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


SCHEMA = "eq45_h_regime_contract_v1"
_SOURCE_CLASS = "public_source_fact"
_CANDIDATE_H_CLASS = "autonomous_design"
_REQUIRED_FALSE_CLAIMS = (
    "paper_exact",
    "pde_validated",
    "visual_correspondence_verified",
    "openai_field_identified",
    "blowup_proved",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def audit_eq45_h_regime(
    contract: Mapping[str, Any],
    *,
    source_contract: Mapping[str, Any],
    candidate: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit one serialized Eq. (4.5) candidate against the two public h scopes."""
    _require(contract.get("schema") == SCHEMA, "unsupported h-regime contract schema")

    source_binding = contract.get("source_binding")
    _require(isinstance(source_binding, Mapping), "missing source_binding")
    _require(
        source_binding.get("source_contract") == "configs/eq45_source_contract.json",
        "source contract path drift",
    )
    _require(
        source_binding.get("classification") == _SOURCE_CLASS,
        "h source range must remain public_source_fact",
    )

    coordinate = source_contract.get("coordinate_contract")
    _require(isinstance(coordinate, Mapping), "source contract lacks coordinate_contract")
    coordinate_interval = coordinate.get("h_open_interval")
    theorem_upper = coordinate.get("theorem_h_upper_bound_strict")
    _require(coordinate_interval == [0.0, 0.5], "source coordinate h interval drift")
    _require(theorem_upper == 0.01, "source theorem h bound drift")
    _require(
        source_binding.get("coordinate_h_open_interval") == coordinate_interval,
        "h-regime contract does not match source coordinate interval",
    )
    _require(
        source_binding.get("theorem_h_upper_bound_strict") == theorem_upper,
        "h-regime contract does not match source theorem bound",
    )

    candidate_binding = contract.get("candidate_binding")
    _require(isinstance(candidate_binding, Mapping), "missing candidate_binding")
    _require(
        candidate_binding.get("candidate_schema") == "eq45_velocity_candidate_v1",
        "candidate schema binding drift",
    )
    _require(
        candidate_binding.get("h_classification") == _CANDIDATE_H_CLASS,
        "candidate h must remain autonomous_design",
    )
    _require(
        candidate.get("schema") == candidate_binding.get("candidate_schema"),
        "candidate schema mismatch",
    )

    classification = candidate.get("classification")
    _require(isinstance(classification, Mapping), "candidate lacks classification metadata")
    _require(
        classification.get("h") == _CANDIDATE_H_CLASS,
        "chosen candidate h must be classified autonomous_design",
    )

    h = candidate.get("h")
    _require(not isinstance(h, bool) and np.isscalar(h), "candidate h must be a scalar")
    h_value = float(h)
    _require(np.isfinite(h_value), "candidate h must be finite")
    lower, upper = (float(value) for value in coordinate_interval)
    coordinate_regime_valid = lower < h_value < upper
    _require(coordinate_regime_valid, "candidate h lies outside the Eq. (4.1) coordinate regime")
    theorem_regime_matched = lower < h_value < float(theorem_upper)

    semantics = contract.get("delivery_semantics")
    _require(isinstance(semantics, Mapping), "missing delivery_semantics")
    _require(
        semantics.get("coordinate_regime_required_for_callable_candidate") is True,
        "callable candidates must remain inside the coordinate regime",
    )
    _require(
        semantics.get("theorem_regime_required_for_velocity_export") is False,
        "the theorem h regime must not become an export gate",
    )
    for key in (
        "theorem_regime_match_is_identity_evidence",
        "theorem_regime_match_is_pde_evidence",
        "theorem_regime_match_is_visual_correspondence_evidence",
    ):
        _require(semantics.get(key) is False, f"unsupported theorem-regime inference: {key}")

    forbidden = contract.get("forbidden_promotions")
    _require(
        isinstance(forbidden, list) and set(forbidden) == set(_REQUIRED_FALSE_CLAIMS),
        "forbidden claim promotions drifted",
    )

    truth = candidate.get("truth_boundary")
    _require(isinstance(truth, Mapping), "candidate lacks truth_boundary")
    _require(truth.get("velocity_export_ready") is True, "seed candidate must remain export-ready")
    for key in _REQUIRED_FALSE_CLAIMS:
        _require(truth.get(key) is False, f"unsupported candidate claim promotion: {key}")

    return {
        "contract_pass": True,
        "candidate_h": h_value,
        "candidate_h_classification": _CANDIDATE_H_CLASS,
        "source_h_regime_classification": _SOURCE_CLASS,
        "coordinate_h_open_interval": [lower, upper],
        "coordinate_regime_valid": True,
        "theorem_h_upper_bound_strict": float(theorem_upper),
        "theorem_regime_matched": theorem_regime_matched,
        "velocity_export_ready": True,
        "theorem_regime_required_for_velocity_export": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
    }


def audit_checked_eq45_h_regime(repo_root: str | Path) -> dict[str, Any]:
    """Audit the checked-in contract and frozen Eq. (4.5) seed artifact."""
    root = Path(repo_root)
    contract = _read_json(root / "configs" / "eq45_h_regime_contract.json")
    source_contract = _read_json(root / "configs" / "eq45_source_contract.json")
    candidate = _read_json(root / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json")
    return audit_eq45_h_regime(
        contract,
        source_contract=source_contract,
        candidate=candidate,
    )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    print(json.dumps(audit_checked_eq45_h_regime(root), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
