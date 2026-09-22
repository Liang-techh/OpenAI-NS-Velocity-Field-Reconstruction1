"""Fail-closed CR002 audit for downstream I2 observability evidence transfer.

A later I3/I4 descendant may have a different semantic identity, deterministic
save/load, or a scoped divergence audit while the original #1061 I2 repair is
still sub-ULP at the public float64 boundary. None of those descendant facts
retroactively proves that the I2 repair itself changed public ``[u,v,w]``.

The Decimal witness below is autonomous mechanics evidence only, not Kokuno or
OpenAI candidate data.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

CONTRACT_RELATIVE_PATH = Path("configs/kokuno_i2_downstream_observability_firewall.json")
EXPECTED_CONTRACT_ID = "cr002-kokuno-i2-downstream-observability-firewall-v1"
EXPECTED_BASE_A5_HEAD = "2886d888cafa80a30e9b21ffee5011e079d944a0"
EXPECTED_BASE_A5_BLOB = "505fc333463db23030a4eee72da30755317982fb"
EXPECTED_I2_HEAD = "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
EXPECTED_I1_HEAD = "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
EXPECTED_PRIOR_CR002_HEAD = "66438e9b4871d676ab8cc2eccb2ee959d1c18b93"
EXPECTED_PRIOR_CR002_CONTRACT_BLOB = "9c4acf91e3f0cc76979fa0c9119f00b1266a884e"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

_REQUIRED_PROVENANCE_CLASSES = {
    "user_requirements",
    "public_source_facts",
    "autonomous_repository_design",
    "pending_or_unknown",
}

_EXPECTED_DOWNSTREAM = {
    "agent1_i4": (1079, "b06742ca6e189499192ede3cce40f62cdc1e35ca", "6f04ce0a856b44430402576dad88438da90d1ebb"),
    "agent2_i4_composite": (1080, "c40d8ddecd2971544a6e07dab093436b423cf326", "2a0a5aa5966b02da856bdcf51940f3c186042802"),
    "agent4_i4_divergence_audit": (1082, "0ff663fcabaa52af984d273ba87356a72481c85c", "518208d1feac5b0f09d08e9989b2eba9275fbd82"),
    "agent3_i2_stress": (1081, "cb9ee25d0c466d42ec7ca30ee89cc10517d4b5c4", "876cdd2cbe4060b9d249cfc4960c26da0d577ca9"),
}

_EXPECTED_REPRESENTATION_TRUTH = {
    "i2_high_precision_correction_receipt_nonzero": True,
    "i2_decimal_retained_until_float64_api_boundary": True,
    "i2_public_float64_delta_from_i1_parent_verified": False,
    "i2_public_float64_delta_probe_protocol_executed": False,
    "i2_public_float64_visual_effect_verified": False,
    "downstream_i4_leading_materialized": True,
    "downstream_i4_composite_materialized": True,
    "downstream_i4_identity_preserving_save_load_available": True,
    "downstream_i4_scoped_divergence_audit_present": True,
    "downstream_i4_scoped_divergence_audit_admitted": False,
    "downstream_i2_compact_stress_materialized": True,
    "downstream_semantic_progress_implies_i2_public_delta": False,
    "downstream_save_load_implies_i2_public_delta": False,
    "downstream_scoped_divergence_implies_i2_public_delta": False,
    "downstream_i2_stress_implies_i2_public_delta": False,
    "descendant_public_difference_from_i1_parent_isolates_i2_repair": False,
    "i2_public_delta_implies_visual_correspondence": False,
    "i2_public_delta_implies_pde_validation": False,
}

_REQUIRED_FORBIDDEN = {
    "I4 descendant semantic identity -> I2 public float64 delta verified",
    "I4 descendant save/load replay -> I2 public float64 delta verified",
    "I4 scoped divergence audit -> I2 public float64 delta verified",
    "I2 nonlinear mean or compact stress -> I2 public float64 delta verified",
    "descendant public difference from I1 parent -> I2 repair isolated as the cause",
    "CI green -> I2 public float64 delta verified",
    "I2 public float64 delta -> visual_correspondence_verified",
    "visual resemblance -> pde_validated",
    "callable or save/load -> paper_exact or openai_field_identified",
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
        raise TypeError("CR002 downstream-observability contract must be a JSON object")
    return payload


def _mapping(value: Any, name: str, violations: list[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        violations.append(f"{name} must be an object")
        return {}
    return value


def _nonempty_string_list(value: Any, name: str, violations: list[str]) -> None:
    if not isinstance(value, list) or not value or any(not isinstance(v, str) or not v.strip() for v in value):
        violations.append(f"{name} must be a non-empty list of non-empty strings")


def audit_contract(payload: Mapping[str, Any]) -> list[str]:
    violations: list[str] = []
    if payload.get("schema_version") != 1:
        violations.append("schema_version must remain 1")
    if payload.get("contract_id") != EXPECTED_CONTRACT_ID:
        violations.append("contract_id changed")
    if payload.get("scope") != "representation_evidence_lineage_only":
        violations.append("scope must remain representation-evidence lineage only")

    target = _mapping(payload.get("target"), "target", violations)
    expected_target = {
        "base_a5_pr": 1083,
        "base_a5_head": EXPECTED_BASE_A5_HEAD,
        "base_a5_source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i4_routing_i2_stress_ingest.py",
        "base_a5_source_blob": EXPECTED_BASE_A5_BLOB,
        "origin_i2_pr": 1061,
        "origin_i2_head": EXPECTED_I2_HEAD,
        "origin_i1_parent_pr": 1051,
        "origin_i1_parent_head": EXPECTED_I1_HEAD,
        "prior_cr002_pr": 1066,
        "prior_cr002_head": EXPECTED_PRIOR_CR002_HEAD,
        "prior_cr002_contract_blob": EXPECTED_PRIOR_CR002_CONTRACT_BLOB,
    }
    for key, expected in expected_target.items():
        if target.get(key) != expected:
            violations.append(f"target.{key} changed")

    downstream = _mapping(payload.get("downstream_lineage"), "downstream_lineage", violations)
    if set(downstream) != set(_EXPECTED_DOWNSTREAM):
        violations.append("downstream_lineage must contain exactly the four pinned descendants")
    for name, (pr, head, blob) in _EXPECTED_DOWNSTREAM.items():
        item = _mapping(downstream.get(name), f"downstream_lineage.{name}", violations)
        if item.get("pr") != pr or item.get("head") != head or item.get("source_blob") != blob:
            violations.append(f"downstream_lineage.{name} identity changed")
        if not isinstance(item.get("role"), str) or not item.get("role", "").strip():
            violations.append(f"downstream_lineage.{name}.role must remain explicit")

    classes = _mapping(payload.get("provenance_classes"), "provenance_classes", violations)
    if set(classes) != _REQUIRED_PROVENANCE_CLASSES:
        violations.append("provenance_classes must contain exactly the four governed classes")
    for key in sorted(_REQUIRED_PROVENANCE_CLASSES):
        _nonempty_string_list(classes.get(key), f"provenance_classes.{key}", violations)

    truth = _mapping(payload.get("representation_truth"), "representation_truth", violations)
    for key, expected in _EXPECTED_REPRESENTATION_TRUTH.items():
        if truth.get(key) is not expected:
            violations.append(f"representation_truth.{key} must remain {expected!r}")

    forbidden = payload.get("forbidden_implications")
    if not isinstance(forbidden, list):
        violations.append("forbidden_implications must be a list")
    else:
        missing = sorted(_REQUIRED_FORBIDDEN - set(forbidden))
        if missing:
            violations.append("forbidden_implications missing: " + ", ".join(missing))

    promotion = _mapping(payload.get("promotion_requirements"), "promotion_requirements", violations)
    for key in ("i2_public_float64_delta_from_i1_parent_verified", "visual_correspondence_verified", "pde_validated"):
        _nonempty_string_list(promotion.get(key), f"promotion_requirements.{key}", violations)
    delta_req = " ".join(promotion.get("i2_public_float64_delta_from_i1_parent_verified", []) if isinstance(promotion.get("i2_public_float64_delta_from_i1_parent_verified"), list) else []).lower()
    required_tokens = ("#1061", "#1051", "do not substitute", "save", "reload", "strict-interior-i2", "float64", "ulp", "nonzero", "negative", "do not infer", "do not retune")
    for token in required_tokens:
        if token not in delta_req:
            violations.append(f"I2 public-delta promotion requirements must mention {token}")
    for forbidden_substitute in ("#1079", "#1080"):
        if forbidden_substitute not in delta_req:
            violations.append(f"I2 promotion firewall must name downstream non-substitute {forbidden_substitute}")

    cr001 = _mapping(payload.get("cr001_freeze"), "cr001_freeze", violations)
    for key, expected in _EXPECTED_CR001.items():
        if cr001.get(key) != expected:
            violations.append(f"cr001_freeze.{key} changed")

    states = _mapping(payload.get("independent_project_states"), "independent_project_states", violations)
    for key, expected in _EXPECTED_STATES.items():
        if states.get(key) is not expected:
            violations.append(f"independent_project_states.{key} must remain {expected!r}")

    witness = _mapping(payload.get("mechanics_witness"), "mechanics_witness", violations)
    if witness.get("classification") != "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data":
        violations.append("mechanics witness classification changed")

    boundary = payload.get("truth_boundary")
    if not isinstance(boundary, str) or "changes no velocity" not in boundary.lower():
        violations.append("truth_boundary must state that no velocity byte changes")
    return violations


def descendant_nonisolation_witness() -> dict[str, Any]:
    """Synthetic mechanics-only witness for downstream non-isolation."""
    with localcontext() as context:
        context.prec = 80
        parent = Decimal("1")
        i2_delta = Decimal("1e-40")
        later_delta = Decimal("1e-6")
        i2_exact = parent + i2_delta
        descendant_exact = i2_exact + later_delta
    parent64 = float(parent)
    i2_64 = float(i2_exact)
    descendant64 = float(descendant_exact)
    return {
        "classification": "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data",
        "i2_high_precision_delta_nonzero": i2_exact != parent,
        "i2_public_float64_unchanged": i2_64 == parent64,
        "descendant_public_float64_changed": descendant64 != parent64,
        "descendant_difference_does_not_isolate_i2": (i2_64 == parent64 and descendant64 != parent64),
        "synthetic_i2_delta": str(i2_delta),
        "synthetic_later_delta": str(later_delta),
    }


def git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def repository_live_violations(root: Path | None = None) -> list[str]:
    base = repository_root() if root is None else Path(root)
    violations: list[str] = []
    constraints_path = base / "configs/constraints.json"
    status_path = base / "project_status.json"
    a5_path = base / "src/openai_ns_reconstruction/kokuno_a5_current_i4_routing_i2_stress_ingest.py"
    for path in (constraints_path, status_path, a5_path):
        if not path.is_file():
            violations.append(f"missing repository-live file: {path.relative_to(base)}")
    if violations:
        return violations

    constraints_bytes = constraints_path.read_bytes()
    if git_blob_sha(constraints_bytes) != EXPECTED_CONSTRAINTS_BLOB:
        violations.append("canonical configs/constraints.json blob changed")
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    if constraints.get("nu") != 0.01:
        violations.append("repository-live nu changed")
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

    status = json.loads(status_path.read_text(encoding="utf-8"))
    states = status.get("states", {})
    if states.get("velocity_export_ready") is not True:
        violations.append("canonical Eq45 velocity_export_ready must remain true")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        if states.get(key) is not False:
            violations.append(f"canonical Eq45 project_status.states.{key} must remain false")

    a5_bytes = a5_path.read_bytes()
    if git_blob_sha(a5_bytes) != EXPECTED_BASE_A5_BLOB:
        violations.append("exact A5 #1083 source blob changed")
    a5_text = a5_bytes.decode("utf-8")
    required_a5_tokens = (
        '"velocity_frontier_stage": "current-I4"',
        '"correction_frontier_stage": "current-I2"',
        '"terminal_global_leading_completion_materialized": False',
        '"global_leading_plus_oscillatory_velocity_materialized": False',
        '"pde_validated": False',
        '"paper_exact": False',
        '"openai_field_identified": False',
    )
    for token in required_a5_tokens:
        if token not in a5_text:
            violations.append(f"exact A5 #1083 truth token missing: {token}")
    if '"public_float64_I2_delta_verified": True' in a5_text:
        violations.append("A5 #1083 unexpectedly promotes I2 public float64 observability")
    return violations


def report(root: Path | None = None) -> dict[str, Any]:
    payload = load_contract(root)
    return {
        "contract_id": payload.get("contract_id"),
        "violations": audit_contract(payload),
        "repository_live_violations": repository_live_violations(root),
        "mechanics_witness": descendant_nonisolation_witness(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    result = report(args.repo_root)
    print(json.dumps(result, indent=2 if args.json else None, sort_keys=True))
    witness = result["mechanics_witness"]
    ok = (
        not result["violations"]
        and not result["repository_live_violations"]
        and witness["i2_high_precision_delta_nonzero"]
        and witness["i2_public_float64_unchanged"]
        and witness["descendant_public_float64_changed"]
        and witness["descendant_difference_does_not_isolate_i2"]
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
