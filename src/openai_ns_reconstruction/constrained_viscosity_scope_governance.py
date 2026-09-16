"""Fail-closed governance for function-first versus CR001 viscosity scope."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping


_SCHEMA = "function_first_viscosity_scope_v1"
_ALLOWED_ORIGINS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _finite_positive(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite positive number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


def _require_origin(value: Any, name: str) -> str:
    if value not in _ALLOWED_ORIGINS:
        raise ValueError(
            f"{name} must use the canonical source-classification vocabulary"
        )
    return str(value)


def audit_viscosity_scope(
    problem_config: Mapping[str, Any],
    selected_profile_artifact: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit viscosity provenance without upgrading delivery or PDE claims."""

    contract = _mapping(contract, "contract")
    if contract.get("schema") != _SCHEMA:
        raise ValueError(f"contract schema must be {_SCHEMA}")

    active = _mapping(contract.get("active_problem"), "active_problem")
    profile = _mapping(
        contract.get("selected_function_first_profile"),
        "selected_function_first_profile",
    )
    compatibility = _mapping(contract.get("compatibility"), "compatibility")
    claims = _mapping(contract.get("claim_status"), "claim_status")

    if (
        _require_origin(active.get("classification"), "active_problem.classification")
        != "autonomous_design"
    ):
        raise ValueError(
            "the CR001 nu=0.01 choice must remain classified autonomous_design"
        )
    if (
        _require_origin(
            profile.get("classification"),
            "selected_function_first_profile.classification",
        )
        != "autonomous_design"
    ):
        raise ValueError(
            "the finite function-first viscosity normalization is an autonomous "
            "numerical choice unless separately sourced"
        )

    nu_key = active.get("nu_key")
    if not isinstance(nu_key, str) or not nu_key:
        raise ValueError("active_problem.nu_key must be a nonempty string")
    problem_nu = _finite_positive(
        problem_config.get(nu_key), f"problem_config.{nu_key}"
    )
    expected_problem_nu = _finite_positive(
        active.get("expected_nu"), "active_problem.expected_nu"
    )
    if not math.isclose(
        problem_nu, expected_problem_nu, rel_tol=0.0, abs_tol=1e-15
    ):
        raise ValueError("active problem viscosity drifted from the governed CR001 value")

    metadata = _mapping(
        selected_profile_artifact.get("metadata"),
        "selected_profile_artifact.metadata",
    )
    metadata_key = profile.get("metadata_key")
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError(
            "selected_function_first_profile.metadata_key must be a nonempty string"
        )
    profile_nu = _finite_positive(
        metadata.get(metadata_key),
        f"selected_profile_artifact.metadata.{metadata_key}",
    )
    expected_profile_nu = _finite_positive(
        profile.get("expected_profile_equation_nu"),
        "selected_function_first_profile.expected_profile_equation_nu",
    )
    if not math.isclose(
        profile_nu, expected_profile_nu, rel_tol=0.0, abs_tol=1e-15
    ):
        raise ValueError("selected function-first profile viscosity metadata drifted")

    same_viscosity = math.isclose(
        problem_nu, profile_nu, rel_tol=0.0, abs_tol=1e-15
    )
    if (
        compatibility.get("same_viscosity_required_for_direct_pde_evidence_transfer")
        is not True
    ):
        raise ValueError(
            "direct PDE evidence transfer must require viscosity compatibility"
        )
    if compatibility.get("explicit_mapping_required_when_values_differ") is not True:
        raise ValueError(
            "a viscosity mismatch must require an explicit mapping or new experiment"
        )
    if compatibility.get("mapping_status") != "pending_unknown":
        raise ValueError(
            "viscosity mapping remains pending_unknown until independent evidence is recorded"
        )
    if compatibility.get("direct_pde_evidence_transfer_allowed") is not False:
        raise ValueError(
            "direct PDE evidence transfer is forbidden for the current viscosity mismatch"
        )
    if (
        not same_viscosity
        and claims.get("function_first_profile_pde_contract_compatible") is not False
    ):
        raise ValueError(
            "mismatched viscosity scopes cannot be labeled PDE-contract-compatible"
        )

    for forbidden_true in (
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
    ):
        if claims.get(forbidden_true) is not False:
            raise ValueError(
                f"{forbidden_true} cannot be promoted by viscosity-scope governance"
            )
    if claims.get("function_first_profile_callable") is not True:
        raise ValueError(
            "the existing callable near-axis profile must not be disabled by a PDE-scope mismatch"
        )

    return {
        "schema": _SCHEMA,
        "problem_nu": problem_nu,
        "profile_equation_nu": profile_nu,
        "same_viscosity": same_viscosity,
        "direct_pde_evidence_transfer_allowed": False,
        "mapping_status": "pending_unknown",
        "function_first_profile_callable": True,
        "function_first_profile_pde_contract_compatible": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
    }


def audit_repository(root: str | Path) -> dict[str, Any]:
    """Audit the current repository files named by the machine-readable contract."""

    root = Path(root)
    contract_path = root / "configs" / "function_first_viscosity_scope.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    active = _mapping(contract.get("active_problem"), "active_problem")
    profile = _mapping(
        contract.get("selected_function_first_profile"),
        "selected_function_first_profile",
    )
    problem_config = json.loads(
        (root / str(active["config_path"])).read_text(encoding="utf-8")
    )
    selected = json.loads(
        (root / str(profile["artifact_path"])).read_text(encoding="utf-8")
    )
    return audit_viscosity_scope(problem_config, selected, contract)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    print(json.dumps(audit_repository(root), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
