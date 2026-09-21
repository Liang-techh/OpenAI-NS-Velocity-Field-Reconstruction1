"""Fail-closed CR002 audit for the current-I2 float64 observability boundary.

This module governs representation evidence only.  In particular, a nonzero
high-precision I2 correction or a new candidate semantic identity is not
evidence that the public float64 ``velocity(x,y,z,t)`` view changed.

The small Decimal witness below is autonomous mechanics evidence.  It is not a
Kokuno/OpenAI candidate value and must never be cited as source or PDE evidence.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT_RELATIVE_PATH = Path("configs/kokuno_i2_float64_observability_scope.json")
EXPECTED_CONTRACT_ID = "cr002-kokuno-i2-float64-observability-scope-v1"
EXPECTED_UPSTREAM_HEAD = "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
EXPECTED_PARENT_HEAD = "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

_REQUIRED_PROVENANCE_CLASSES = {
    "user_requirements",
    "public_source_facts",
    "autonomous_repository_design",
    "pending_or_unknown",
}

_EXPECTED_REPRESENTATION_TRUTH = {
    "high_precision_i2_correction_receipt_nonzero": True,
    "decimal_retained_until_float64_api_boundary": True,
    "public_float64_velocity_callable_through_i2": True,
    "configuration_save_load_semantic_identity_available": True,
    "public_float64_delta_from_i1_parent_verified": False,
    "public_float64_delta_probe_protocol_registered": False,
    "public_float64_visual_effect_verified": False,
    "high_precision_nonzero_implies_public_float64_delta": False,
    "semantic_identity_progress_implies_public_float64_delta": False,
    "callable_through_i2_implies_visual_correspondence": False,
}

_REQUIRED_FORBIDDEN_IMPLICATIONS = {
    "high_precision_nonzero_correction -> nonzero_public_float64_velocity_delta",
    "new_semantic_identity -> nonzero_public_float64_velocity_delta",
    "callable_through_I2 -> visual_correspondence_verified",
    "public_float64_velocity_delta -> visual_correspondence_verified",
    "visual_correspondence -> pde_validated",
    "scoped_or_local_audit -> complete_ns_validation",
    "ci_green -> pde_validated",
    "callable_or_save_load -> paper_exact_or_openai_field_identified",
}

_EXPECTED_CR001 = {
    "constraints_blob_sha": EXPECTED_CONSTRAINTS_BLOB,
    "nu": 0.01,
    "physical_domain": "R^3",
    "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
    "support": "r < 2 and abs(z) < 2",
    "time_interval": [0.25, 0.75],
    "forcing_mode": "restricted_two_parameter_family",
    "reference_energy": 1.0,
    "reference_energy_abs_tolerance": 0.001,
    "validation_seed": 914027,
    "held_out_points": 4096,
    "derivative_steps": [0.02, 0.01, 0.005],
    "quadrature_orders_per_axis": [24, 48, 96],
    "momentum_max": 0.001,
    "momentum_L2": 0.001,
    "divergence_max": 1e-5,
    "divergence_L2": 1e-5,
    "residual_defined_free_force_forbidden": True,
    "candidate_collapse_forbidden": True,
    "posthoc_threshold_relaxation_forbidden": True,
}

_EXPECTED_STATES = {
    "canonical_eq45_velocity_export_ready": True,
    "canonical_eq45_visualization_ready": False,
    "canonical_eq45_visual_correspondence_verified": False,
    "canonical_eq45_pde_validated": False,
    "kokuno_current_route_velocity_export_ready": False,
    "kokuno_current_route_visual_correspondence_verified": False,
    "kokuno_current_route_pde_validated": False,
    "kokuno_current_route_paper_exact": False,
    "kokuno_current_route_openai_field_identified": False,
}


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def contract_path(root: Path | None = None) -> Path:
    base = repository_root() if root is None else Path(root)
    return base / CONTRACT_RELATIVE_PATH


def load_contract(root: Path | None = None) -> dict[str, Any]:
    with contract_path(root).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError("CR002 observability contract must be a JSON object")
    return payload


def _mapping(value: Any, name: str, violations: list[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        violations.append(f"{name} must be an object")
        return {}
    return value


def _nonempty_string_list(value: Any, name: str, violations: list[str]) -> None:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item.strip() for item in value)
    ):
        violations.append(f"{name} must be a non-empty list of non-empty strings")


def audit_contract(payload: Mapping[str, Any]) -> list[str]:
    """Return all fail-closed contract violations."""

    violations: list[str] = []
    if payload.get("schema_version") != 1:
        violations.append("schema_version must remain 1")
    if payload.get("contract_id") != EXPECTED_CONTRACT_ID:
        violations.append("contract_id changed")
    if payload.get("scope") != "representation_and_delivery_evidence_only":
        violations.append("scope must remain representation/delivery evidence only")

    target = _mapping(payload.get("target"), "target", violations)
    if target.get("upstream_pr") != 1061:
        violations.append("target.upstream_pr must remain exact PR #1061")
    if target.get("upstream_head") != EXPECTED_UPSTREAM_HEAD:
        violations.append("target.upstream_head changed")
    if target.get("parent_i1_pr") != 1051:
        violations.append("target.parent_i1_pr must remain exact PR #1051")
    if target.get("parent_i1_head") != EXPECTED_PARENT_HEAD:
        violations.append("target.parent_i1_head changed")
    if target.get("module_path") != (
        "src/openai_ns_reconstruction/"
        "kokuno_pa16_current_cartesian_i2_heat_repair.py"
    ):
        violations.append("target.module_path changed")

    classes = _mapping(payload.get("provenance_classes"), "provenance_classes", violations)
    if set(classes) != _REQUIRED_PROVENANCE_CLASSES:
        violations.append("provenance_classes must contain exactly the four governed classes")
    for key in sorted(_REQUIRED_PROVENANCE_CLASSES):
        _nonempty_string_list(classes.get(key), f"provenance_classes.{key}", violations)

    representation = _mapping(
        payload.get("representation_truth"), "representation_truth", violations
    )
    for key, expected in _EXPECTED_REPRESENTATION_TRUTH.items():
        if representation.get(key) is not expected:
            violations.append(f"representation_truth.{key} must remain {expected!r}")

    forbidden = payload.get("forbidden_implications")
    if not isinstance(forbidden, list):
        violations.append("forbidden_implications must be a list")
    else:
        forbidden_set = set(forbidden)
        missing = sorted(_REQUIRED_FORBIDDEN_IMPLICATIONS - forbidden_set)
        if missing:
            violations.append(
                "forbidden_implications missing: " + ", ".join(missing)
            )

    promotion = _mapping(
        payload.get("promotion_requirements"), "promotion_requirements", violations
    )
    for key in (
        "public_float64_delta_from_i1_parent_verified",
        "different_public_precision_or_interface",
        "visual_correspondence_verified",
    ):
        _nonempty_string_list(promotion.get(key), f"promotion_requirements.{key}", violations)
    delta_requirements = " ".join(
        promotion.get("public_float64_delta_from_i1_parent_verified", [])
        if isinstance(
            promotion.get("public_float64_delta_from_i1_parent_verified"), list
        )
        else []
    ).lower()
    for token in ("save", "reload", "#1061", "#1051", "float64", "ulp", "nonzero", "negative"):
        if token not in delta_requirements:
            violations.append(
                "public float64 promotion requirements must mention " + token
            )
    if "do not change" not in delta_requirements:
        violations.append(
            "public float64 promotion requirements must prohibit retuning to manufacture observability"
        )

    cr001 = _mapping(payload.get("cr001_freeze"), "cr001_freeze", violations)
    for key, expected in _EXPECTED_CR001.items():
        if cr001.get(key) != expected:
            violations.append(f"cr001_freeze.{key} changed")

    states = _mapping(
        payload.get("independent_project_states"),
        "independent_project_states",
        violations,
    )
    for key, expected in _EXPECTED_STATES.items():
        if states.get(key) is not expected:
            violations.append(f"independent_project_states.{key} must remain {expected!r}")

    witness = _mapping(payload.get("mechanics_witness"), "mechanics_witness", violations)
    if witness.get("classification") != (
        "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data"
    ):
        violations.append("mechanics witness classification changed")

    boundary = payload.get("truth_boundary")
    if not isinstance(boundary, str) or "changes no velocity" not in boundary.lower():
        violations.append("truth_boundary must state that no velocity byte changes")

    return violations


def subulp_decimal_witness() -> dict[str, Any]:
    """Return an autonomous mechanics-only high-precision/float64 witness.

    This intentionally uses synthetic ``1 + 1e-40`` data.  It demonstrates a
    representation fact only: exact arithmetic and a semantic payload can change
    while the public float64 number remains unchanged after rounding.
    """

    with localcontext() as context:
        context.prec = 80
        base = Decimal("1")
        delta = Decimal("1e-40")
        shifted = base + delta

    base64 = float(base)
    shifted64 = float(shifted)
    payload_a = json.dumps(
        {"base": str(base), "delta": "0"}, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    payload_b = json.dumps(
        {"base": str(base), "delta": str(delta)},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    return {
        "classification": "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data",
        "high_precision_delta_nonzero": shifted != base,
        "float64_values_equal": shifted64 == base64,
        "float64_ulp_at_base": math.ulp(base64),
        "synthetic_delta": str(delta),
        "semantic_payload_sha_changed": (
            hashlib.sha256(payload_a).hexdigest()
            != hashlib.sha256(payload_b).hexdigest()
        ),
    }


def git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def repository_live_violations(root: Path | None = None) -> list[str]:
    """Replay the contract against repository state present in the checkout."""

    base = repository_root() if root is None else Path(root)
    violations: list[str] = []

    constraints_path = base / "configs/constraints.json"
    project_status_path = base / "project_status.json"
    upstream_path = base / (
        "src/openai_ns_reconstruction/"
        "kokuno_pa16_current_cartesian_i2_heat_repair.py"
    )
    for path in (constraints_path, project_status_path, upstream_path):
        if not path.is_file():
            violations.append(f"missing repository-live file: {path.relative_to(base)}")
    if violations:
        return violations

    constraints_bytes = constraints_path.read_bytes()
    if git_blob_sha(constraints_bytes) != EXPECTED_CONSTRAINTS_BLOB:
        violations.append("canonical configs/constraints.json blob changed")

    constraints = json.loads(constraints_bytes.decode("utf-8"))
    if constraints.get("nu") != 0.01:
        violations.append("repository-live nu changed")
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    if validation.get("seed") != 914027 or validation.get("held_out_points") != 4096:
        violations.append("repository-live held-out protocol changed")
    if validation.get("derivative_steps") != [0.02, 0.01, 0.005]:
        violations.append("repository-live derivative ladder changed")
    if validation.get("quadrature_orders_per_axis") != [24, 48, 96]:
        violations.append("repository-live quadrature ladder changed")
    if thresholds.get("pde_residual_max") != 0.001 or thresholds.get("pde_residual_L2") != 0.001:
        violations.append("repository-live momentum/PDE threshold changed")
    if thresholds.get("divergence_max") != 1e-5 or thresholds.get("divergence_L2") != 1e-5:
        violations.append("repository-live divergence threshold changed")
    restriction = constraints.get("forcing", {}).get("restriction", "")
    if "No residual-dependent basis or pointwise free force" not in restriction:
        violations.append("repository-live free-force prohibition changed")

    status = json.loads(project_status_path.read_text(encoding="utf-8"))
    states = status.get("states", {})
    if states.get("velocity_export_ready") is not True:
        violations.append("canonical Eq45 velocity_export_ready must remain true")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        if states.get(key) is not False:
            violations.append(f"canonical Eq45 project_status.states.{key} must remain false")

    upstream_text = upstream_path.read_text(encoding="utf-8")
    required_tokens = (
        "physical correction can be sub-ulp in float64",
        '"float_profile_view"',
        "materialized to float64 only at the candidate API boundary",
        "def i2_decimal_correction_report",
        '"current_I2_overlay_applied": True',
        '"cartesian_velocity_executable_through_I2": True',
        '"unified_global_cartesian_velocity_export_ready": False',
        '"pde_validated": False',
        '"paper_exact": False',
        '"openai_field_identified": False',
    )
    for token in required_tokens:
        if token not in upstream_text:
            violations.append(f"exact-I2 representation token missing: {token}")

    return violations


def _report(root: Path | None = None) -> dict[str, Any]:
    payload = load_contract(root)
    return {
        "contract_id": payload.get("contract_id"),
        "violations": audit_contract(payload),
        "repository_live_violations": repository_live_violations(root),
        "mechanics_witness": subulp_decimal_witness(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    report = _report(args.repo_root)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(json.dumps(report, sort_keys=True))
    witness = report["mechanics_witness"]
    ok = (
        not report["violations"]
        and not report["repository_live_violations"]
        and witness["high_precision_delta_nonzero"]
        and witness["float64_values_equal"]
        and witness["semantic_payload_sha_changed"]
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
