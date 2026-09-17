"""Fail-closed source governance for the paper Eq. (4.1)/(4.5) backbone.

This module validates provenance and representation metadata only. It does not
validate a numerical q solver, profile fit, Navier--Stokes residual, visual
correspondence, or paper-exact recovery.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_REQUIRED_FALSE_CLAIMS = (
    "final_profiles_numerically_identified",
    "full_constructed_field_recovered",
    "openai_correspondence_verified",
    "pde_validated",
    "paper_exact",
)
_EXPECTED_CLASSIFICATIONS = (
    ("Eq. (4.1) similarity coordinates", "public_source_fact"),
    ("callable velocity(x,y,z,t)", "user_requirement"),
    ("finite numerical choices", "autonomous_design"),
    ("final numerical profile data", "pending_unknown"),
)

_EXPECTED_COORDINATE = {
    "tau": "1-t",
    "A": "1/2+h",
    "D": "1/2-h",
    "z_relation": "z=q^D*eta",
    "tau_relation": "tau=q*(1-eta^2)",
    "implicit_q_relation": "q-z^2*q^(2*h)=1-t",
    "q_exponent_sign": "positive",
    "physical_root_condition": "q>|z|^(1/D)",
    "eta": "z/q^D",
    "X": "(x^2+y^2)/(2*q)",
}

_EXPECTED_VELOCITY = {
    "u": "x*v0(X,eta)/(2*q)-y*q^(-1-h)*F(X,eta)",
    "v": "y*v0(X,eta)/(2*q)+x*q^(-1-h)*F(X,eta)",
    "w": "q^(-1/2-h)*U(X,eta)",
    "profile_relation": "V0=X*v0",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _classification_for(rows: list[Mapping[str, Any]], fragment: str) -> str:
    matches = [row for row in rows if fragment in str(row.get("item", ""))]
    _require(len(matches) == 1, f"expected one source classification row matching {fragment!r}")
    return str(matches[0].get("classification"))


def audit_eq45_source_contract(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Validate the machine-readable public-source Eq. (4.1)/(4.5) contract."""
    root = Path(repo_root)
    _require(contract.get("schema_version") == 1, "unsupported schema_version")

    source = contract.get("source")
    _require(isinstance(source, Mapping), "missing source metadata")
    _require(source.get("coordinate_equation") == "4.1", "coordinate equation source drift")
    _require(source.get("coordinate_printed_page") == 24, "coordinate source page drift")
    _require(source.get("velocity_equation") == "4.5", "velocity equation source drift")
    _require(source.get("velocity_printed_page") == 25, "velocity source page drift")
    url = source.get("url")
    _require(isinstance(url, str) and url.startswith("https://cdn.openai.com/pdf/"), "unexpected source URL")

    coordinate = contract.get("coordinate_contract")
    _require(isinstance(coordinate, Mapping), "missing coordinate_contract")
    for key, expected in _EXPECTED_COORDINATE.items():
        _require(coordinate.get(key) == expected, f"Eq. (4.1) coordinate drift: {key}")
    _require(coordinate.get("h_open_interval") == [0.0, 0.5], "coordinate h-domain drift")
    _require(
        coordinate.get("theorem_h_upper_bound_strict") == 0.01,
        "theorem h upper bound drift",
    )

    velocity = contract.get("velocity_contract")
    _require(isinstance(velocity, Mapping), "missing velocity_contract")
    for key, expected in _EXPECTED_VELOCITY.items():
        _require(velocity.get(key) == expected, f"Eq. (4.5) velocity drift: {key}")
    _require(velocity.get("leading_field_only") is True, "Eq. (4.5) scope must remain leading-field only")

    rows_raw = contract.get("source_classification")
    _require(isinstance(rows_raw, list) and rows_raw, "source_classification must be nonempty")
    rows: list[Mapping[str, Any]] = []
    observed: set[str] = set()
    for row in rows_raw:
        _require(isinstance(row, Mapping), "source classification row must be an object")
        classification = row.get("classification")
        _require(classification in _ALLOWED_CLASSES, f"invalid source classification: {classification!r}")
        observed.add(str(classification))
        path = row.get("source")
        _require(isinstance(path, str) and path, "source classification row lacks source")
        _require((root / path).is_file(), f"source path does not exist: {path}")
        rows.append(row)
    _require(observed == _ALLOWED_CLASSES, "source classification coverage is incomplete")
    for fragment, expected_class in _EXPECTED_CLASSIFICATIONS:
        actual_class = _classification_for(rows, fragment)
        _require(
            actual_class == expected_class,
            f"source classification drift for {fragment!r}: expected {expected_class!r}",
        )

    status = contract.get("claim_status")
    _require(isinstance(status, Mapping), "missing claim_status")
    _require(status.get("coordinate_formula_publicly_sourced") is True, "coordinate source status drift")
    _require(status.get("leading_velocity_formula_publicly_sourced") is True, "velocity source status drift")
    for key in _REQUIRED_FALSE_CLAIMS:
        _require(status.get(key) is False, f"unsupported claim promotion: {key}")

    truth = contract.get("truth_boundary")
    _require(isinstance(truth, list) and truth, "truth_boundary must be nonempty")
    _require(
        any("q^(+2h)" in str(item) and "q^(-2h)" in str(item) for item in truth),
        "truth boundary must explicitly guard the q exponent sign",
    )

    return {
        "contract_pass": True,
        "q_relation": coordinate.get("implicit_q_relation"),
        "q_exponent_sign": coordinate.get("q_exponent_sign"),
        "coordinate_h_open_interval": coordinate.get("h_open_interval"),
        "leading_field_only": True,
        "source_classes": sorted(observed),
        "pde_validated": False,
        "paper_exact": False,
    }


def audit_eq45_source_contract_file(path: str | Path, *, repo_root: str | Path) -> dict[str, Any]:
    return audit_eq45_source_contract(_read_json(Path(path)), repo_root=repo_root)
