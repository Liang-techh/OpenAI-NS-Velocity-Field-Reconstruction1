from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/kokuno_x4_mean_audit_registration_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")

CONTRACT_ID = "CR002-KOKUNO-X4-MEAN-AUDIT-REGISTRATION-SCOPE-109"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_A5_HEAD = "2dab469f3d48268bfb7edac18ae454b8265a99c0"
EXPECTED_A3_HEAD = "4b7a4400dcf1ee0be6bb07f23c2f8af6f893de88"
EXPECTED_A4_HEAD = "8a34d11605d095334175df9618f8f5b18f5c210e"


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def audit_contract(contract: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    if contract.get("schema_version") != 1:
        errors.append("schema_version must remain 1")
    if contract.get("contract_id") != CONTRACT_ID:
        errors.append("unexpected contract_id")
    if contract.get("status") != "governance_only_no_scientific_promotion":
        errors.append("status must remain governance-only")

    base = contract.get("base_a5", {})
    if base.get("pr") != 1030:
        errors.append("base_a5.pr must remain 1030")
    if base.get("head") != EXPECTED_A5_HEAD:
        errors.append("base_a5.head drift")

    current = contract.get("current_repository_evidence", {})
    expected_current = {
        "agent3_x4_mean_pr": 1028,
        "agent3_x4_mean_head": EXPECTED_A3_HEAD,
        "agent4_x4_mean_audit_pr": 1029,
        "agent4_x4_mean_audit_head": EXPECTED_A4_HEAD,
        "agent4_x4_mean_audit_artifact_exists": True,
        "agent5_1030_consumes_agent3_1028": False,
        "agent5_1030_consumes_agent4_1029": False,
        "agent5_1030_registers_agent4_1029": False,
    }
    for key, value in expected_current.items():
        if current.get(key) != value:
            errors.append(f"current_repository_evidence.{key} must equal {value!r}")

    allowed = contract.get("allowed_true", {})
    for key in (
        "a3_1028_nonlinear_m0_mean_attribution_materialized_through_x4",
        "a4_1029_implementation_distinct_x4_mean_audit_artifact_exists",
        "a4_1029_uses_save_reloaded_public_cartesian_velocity_reference",
        "a5_1030_registers_prior_x2_radial_force_audit",
    ):
        if allowed.get(key) is not True:
            errors.append(f"allowed_true.{key} must remain true")

    required_false = contract.get("required_false", {})
    expected_false_keys = (
        "a5_1030_consumes_a3_1028_x4_mean",
        "a5_1030_consumes_or_registers_a4_1029_x4_mean_audit",
        "x4_compact_radial_stress_materialized",
        "x4_radial_force_partial_z_sigma1_materialized",
        "x4_mean_audit_means_x4_stress_or_force_audit",
        "x4_mean_audit_means_complete_ns_defect",
        "x4_mean_audit_means_authorized_correction_target",
        "cartesian_correction_velocity_materialized",
        "complete_candidate_api_materialized_on_kokuno_route",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    )
    for key in expected_false_keys:
        if required_false.get(key) is not False:
            errors.append(f"required_false.{key} must remain false")
    unexpected_false_keys = set(required_false) - set(expected_false_keys)
    if unexpected_false_keys:
        errors.append(
            "required_false contains unreviewed keys: "
            + ", ".join(sorted(unexpected_false_keys))
        )

    semantics = contract.get("evidence_semantics", {})
    if "without being consumed" not in semantics.get("artifact_vs_registration_rule", ""):
        errors.append("artifact-vs-registration rule must remain explicit")
    if "does not materialize or validate" not in semantics.get("representation_rule", ""):
        errors.append("mean-audit-vs-stress/force rule must remain explicit")
    if "not a complete NS defect" not in semantics.get("correction_authorization_rule", ""):
        errors.append("correction authorization firewall must remain explicit")
    if "queued/running" not in semantics.get("ci_rule", ""):
        errors.append("CI truth rule must preserve queued/running semantics")

    promotion = contract.get("promotion_requirements", {})
    registration_requirements = promotion.get("future_a5_x4_mean_audit_registration", [])
    if not any("A3 #1028" in item and "A4 #1029" in item for item in registration_requirements):
        errors.append("future A5 registration must bind exact A3 #1028 and A4 #1029 identities")
    force_requirements = promotion.get("future_x4_radial_force_materialization", [])
    if not any("compact stress" in item for item in force_requirements):
        errors.append("future X4 radial force must require compact-stress materialization")
    if not any("complete defect" in item for item in force_requirements):
        errors.append("future X4 radial force must remain unauthorized before complete defect")

    provenance = contract.get("provenance_classes", {})
    for category in ("user_requirement", "public_source_fact", "autonomous_governance", "pending_unknown"):
        if category not in provenance:
            errors.append(f"missing provenance class {category}")
    if provenance.get("public_source_fact") != []:
        errors.append("this increment must not add a new public-source fact")
    pending = " ".join(provenance.get("pending_unknown", []))
    for token in ("X4", "complete NS defect", "global"):
        if token not in pending:
            errors.append(f"pending_unknown must retain {token!r} blocker")

    ci = contract.get("ci_observation_at_creation", {})
    expected_runs = {
        "agent4_1029_tests_run": 35635278480,
        "agent4_1029_dedicated_run": 35635278504,
        "agent5_1030_tests_run": 35635405082,
        "agent5_1030_dedicated_run": 35635405080,
    }
    for key, value in expected_runs.items():
        if ci.get(key) != value:
            errors.append(f"ci_observation_at_creation.{key} must equal {value}")
    if ci.get("scientific_pass_claimed_from_these_observations") is not False:
        errors.append("queued creation-time CI may not be claimed as scientific PASS")

    cr001 = contract.get("cr001_snapshot", {})
    expected_cr001 = {
        "constraints_blob": EXPECTED_CONSTRAINTS_BLOB,
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
        "momentum_max_l2_threshold": 0.001,
        "divergence_max_l2_threshold": 0.00001,
        "residual_defined_free_force_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    for key, value in expected_cr001.items():
        if cr001.get(key) != value:
            errors.append(f"cr001_snapshot.{key} must equal {value!r}")

    canonical = contract.get("canonical_delivery_independence", {})
    expected_canonical = {
        "candidate": "eq45_supported_velocity_candidate_v1",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, value in expected_canonical.items():
        if canonical.get(key) != value:
            errors.append(f"canonical_delivery_independence.{key} must equal {value!r}")

    return errors


def audit_repository(repo_root: Path) -> list[str]:
    errors: list[str] = []
    contract_path = repo_root / CONTRACT_REL
    constraints_path = repo_root / CONSTRAINTS_REL

    if not contract_path.is_file():
        return [f"missing contract: {CONTRACT_REL}"]
    if not constraints_path.is_file():
        return [f"missing canonical constraints: {CONSTRAINTS_REL}"]

    try:
        contract = _load_json(contract_path)
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot parse contract: {exc}"]

    errors.extend(audit_contract(contract))

    constraints_bytes = constraints_path.read_bytes()
    observed_blob = _git_blob_sha(constraints_bytes)
    if observed_blob != EXPECTED_CONSTRAINTS_BLOB:
        errors.append(
            f"canonical constraints blob drift: {observed_blob} != {EXPECTED_CONSTRAINTS_BLOB}"
        )

    try:
        constraints = json.loads(constraints_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        errors.append(f"cannot parse canonical constraints: {exc}")
        return errors

    if constraints.get("nu") != 0.01:
        errors.append("canonical nu drift")
    if constraints.get("domain", {}).get("time_interval") != [0.25, 0.75]:
        errors.append("canonical time interval drift")
    if constraints.get("forcing", {}).get("mode") != "restricted_two_parameter_family":
        errors.append("canonical forcing mode drift")
    if constraints.get("nontriviality", {}).get("reference_energy") != 1.0:
        errors.append("canonical nontriviality reference energy drift")
    thresholds = constraints.get("validation", {}).get("thresholds", {})
    if thresholds.get("pde_residual_max") != 0.001 or thresholds.get("pde_residual_L2") != 0.001:
        errors.append("canonical momentum thresholds drift")
    if thresholds.get("divergence_max") != 1e-5 or thresholds.get("divergence_L2") != 1e-5:
        errors.append("canonical divergence thresholds drift")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fail-closed CR002 audit for RF40 X4 nonlinear-mean evidence registration scope."
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()

    errors = audit_repository(args.repo_root.resolve())
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS: CR002 X4 nonlinear-mean audit registration scope is intact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
