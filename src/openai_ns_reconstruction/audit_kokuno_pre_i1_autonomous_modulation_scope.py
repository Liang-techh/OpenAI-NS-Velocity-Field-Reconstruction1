from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/kokuno_pre_i1_autonomous_modulation_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")

CONTRACT_ID = "CR002-KOKUNO-PRE-I1-AUTONOMOUS-MODULATION-SCOPE-110"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_A1_1035_HEAD = "ddf71bab2a2038bab2fba5b854ad713c6b6c7e95"
EXPECTED_A1_1035_SOURCE_BLOB = "ef124c53163bb452096168648f0e22d9226750fe"
EXPECTED_A1_1005_HEAD = "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"
EXPECTED_A2_1010_HEAD = "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
EXPECTED_A3_1037_HEAD = "9a5cdcdf78b7862d5ebe171efb9ca56d53664612"
EXPECTED_A4_1012_HEAD = "9b2be0ffecc14479dbd941eebd663d63f89f2a2f"
EXPECTED_A4_1038_HEAD = "a4546cacfa224c82398ef7245ba41cca41159800"
EXPECTED_INTEGRATION_HEAD = "7f6658cac1a8ffea178bec5accf2fa5a784b8249"
EXPECTED_PUBLIC_SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
EXPECTED_PUBLIC_SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"


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

    base = contract.get("base_a1_modulated", {})
    if base.get("pr") != 1035:
        errors.append("base_a1_modulated.pr must remain 1035")
    if base.get("head") != EXPECTED_A1_1035_HEAD:
        errors.append("base_a1_modulated.head drift")
    if base.get("source_blob") != EXPECTED_A1_1035_SOURCE_BLOB:
        errors.append("base_a1_modulated.source_blob drift")

    source = contract.get("source_classification", {})
    public = source.get("public_corrected_reconstruction", {})
    expected_public = {
        "repository": "KokunoYumeto/yang-mills-interacting-workbench",
        "commit": EXPECTED_PUBLIC_SOURCE_COMMIT,
        "path": "navier-stokes/navier_stokes_workbench.tex",
        "blob": EXPECTED_PUBLIC_SOURCE_BLOB,
        "release_date": "2026-09-09",
    }
    for key, value in expected_public.items():
        if public.get(key) != value:
            errors.append(f"source_classification.public_corrected_reconstruction.{key} drift")
    if "not an official OpenAI field release" not in public.get("scope", ""):
        errors.append("public-source scope must reject official/OpenAI-field promotion")
    if "does not publish the hidden numerical loop" not in public.get("source_form_statement", ""):
        errors.append("public source-form statement must preserve hidden-loop unknown")

    autonomous = source.get("repository_autonomous_realization", {})
    expected_autonomous = {
        "frequency": 16,
        "alpha": 0.0005,
        "beta": 0.0005,
        "support_left_offset": -35.3,
        "support_right_offset": -29.1,
        "classification": "autonomous_design_not_public_source_parameter_recovery",
    }
    for key, value in expected_autonomous.items():
        if autonomous.get(key) != value:
            errors.append(f"source_classification.repository_autonomous_realization.{key} drift")
    if "Gauss-Legendre" not in autonomous.get("memory_quadrature", ""):
        errors.append("autonomous memory quadrature classification drift")

    lineage = contract.get("lineage_evidence", {})
    expected_lineage = {
        "a1_1005_pr": 1005,
        "a1_1005_head": EXPECTED_A1_1005_HEAD,
        "a1_1035_pr": 1035,
        "a1_1035_head": EXPECTED_A1_1035_HEAD,
        "a2_1010_pr": 1010,
        "a2_1010_head": EXPECTED_A2_1010_HEAD,
        "a2_1010_consumed_leading_head": EXPECTED_A1_1005_HEAD,
        "a3_1037_pr": 1037,
        "a3_1037_head": EXPECTED_A3_1037_HEAD,
        "a3_1037_consumes_unmodulated_a2_1010_lineage": True,
        "a4_1012_pr": 1012,
        "a4_1012_head": EXPECTED_A4_1012_HEAD,
        "a4_1012_audits_unmodulated_a2_1010_composite": True,
        "a4_1038_pr": 1038,
        "a4_1038_head": EXPECTED_A4_1038_HEAD,
        "a4_1038_audits_a3_1037_unmodulated_x4_stress": True,
    }
    for key, value in expected_lineage.items():
        if lineage.get(key) != value:
            errors.append(f"lineage_evidence.{key} must equal {value!r}")

    allowed = contract.get("allowed_true", {})
    expected_allowed = (
        "a1_1035_pre_i1_leading_cartesian_velocity_callable",
        "a1_1035_configuration_serializable",
        "a1_1035_semantic_identity_available",
        "a1_1035_repository_autonomous_modulation_bound",
        "a1_1035_prefix_incompressibility_memory_carried",
        "a1_1035_fails_closed_at_reserved_i1_boundary",
        "corrected_public_source_provenance_pinned",
        "existing_unmodulated_a2_a3_a4_evidence_remains_scoped_to_its_own_identity",
    )
    for key in expected_allowed:
        if allowed.get(key) is not True:
            errors.append(f"allowed_true.{key} must remain true")
    extra_allowed = set(allowed) - set(expected_allowed)
    if extra_allowed:
        errors.append("allowed_true contains unreviewed keys: " + ", ".join(sorted(extra_allowed)))

    required_false = contract.get("required_false", {})
    expected_false_keys = (
        "source_hidden_loop_parameters_recovered",
        "source_admissible_loop_reconstructed",
        "autonomous_modulation_equals_hidden_source_loop",
        "a2_1010_consumes_a1_1035_modulated_identity",
        "a4_1012_validates_a1_1035_modulated_velocity",
        "a3_1037_stress_applies_to_a1_1035_modulated_identity",
        "a4_1038_validates_a1_1035_modulated_identity",
        "leading_plus_oscillatory_modulated_pre_i1_materialized",
        "modulated_pre_i1_independent_cartesian_audit_available",
        "modulated_pre_i1_correction_mean_stress_force_materialized",
        "outer_global_leading_velocity_materialized",
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    )
    for key in expected_false_keys:
        if required_false.get(key) is not False:
            errors.append(f"required_false.{key} must remain false")
    extra_false = set(required_false) - set(expected_false_keys)
    if extra_false:
        errors.append("required_false contains unreviewed keys: " + ", ".join(sorted(extra_false)))

    semantics = contract.get("identity_semantics", {})
    if "does not transfer validation evidence" not in semantics.get("domain_overlap_rule", ""):
        errors.append("domain-overlap identity firewall must remain explicit")
    if "A1 #1005" not in semantics.get("lineage_rule", "") or "A1 #1035" not in semantics.get("lineage_rule", ""):
        errors.append("lineage rule must distinguish A1 #1005 from exact A1 #1035")
    if "does not recover hidden source parameters" not in semantics.get("source_rule", ""):
        errors.append("source/autonomous truth boundary must remain explicit")
    if "not a complete" not in semantics.get("delivery_rule", ""):
        errors.append("leading-only delivery scope must remain explicit")
    if "queued/running" not in semantics.get("ci_rule", ""):
        errors.append("CI truth rule must preserve queued/running semantics")

    promotion = contract.get("promotion_requirements", {})
    composite = promotion.get("future_modulated_composite", [])
    if not any("exact A1 #1035" in item for item in composite):
        errors.append("future modulated composite must consume exact A1 #1035")
    if not any("oscillatory runtime" in item for item in composite):
        errors.append("future modulated composite must bind oscillatory runtime")
    if not any("A4 audit" in item and "same modulated identity" in item for item in composite):
        errors.append("future modulated composite must require same-identity independent A4 audit")

    correction = promotion.get("future_modulated_correction", [])
    if not any("exact modulated composite identity" in item for item in correction):
        errors.append("future modulated correction must recompute from exact modulated identity")
    if not any("A3 #1037" in item and "A4 #1038" in item for item in correction):
        errors.append("future modulated correction must reject transplanted #1037/#1038 evidence")
    if not any("complete-NS correction target" in item for item in correction):
        errors.append("future correction must remain unauthorized before complete defect")

    source_exact = promotion.get("future_source_exactness", [])
    if not any("independently establish" in item for item in source_exact):
        errors.append("future source-exact promotion must require independent evidence")
    if not any("separate gates" in item for item in source_exact):
        errors.append("visual/PDE/exactness gates must remain independent")

    provenance = contract.get("provenance_classes", {})
    for category in ("user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"):
        if category not in provenance:
            errors.append(f"missing provenance class {category}")
    public_text = " ".join(provenance.get("public_source_fact", []))
    if "corrected Kokuno public reconstruction" not in public_text:
        errors.append("public-source provenance must remain specifically Kokuno corrected reconstruction")
    if "hidden OpenAI field" not in public_text:
        errors.append("public-source classification must reject hidden-field identification")
    autonomous_text = " ".join(provenance.get("autonomous_design", []))
    for token in ("N=16", "alpha=beta=5e-4", "-35.3", "-29.1", "Gauss-Legendre"):
        if token not in autonomous_text:
            errors.append(f"autonomous_design must retain {token!r}")
    pending_text = " ".join(provenance.get("pending_unknown", []))
    for token in ("hidden source loop", "A1 #1035", "global", "PDE"):
        if token not in pending_text:
            errors.append(f"pending_unknown must retain {token!r} blocker")

    ci = contract.get("ci_observation_at_creation", {})
    if ci.get("a1_1035_dedicated_run") != 35639014301:
        errors.append("unexpected A1 #1035 dedicated CI run id")
    if ci.get("a1_1035_repository_tests_run") != 35639014464:
        errors.append("unexpected A1 #1035 repository tests run id")
    if ci.get("a1_1035_observed_status") != "queued":
        errors.append("creation-time A1 #1035 CI status must remain queued")
    if ci.get("scientific_pass_claimed_from_these_observations") is not False:
        errors.append("queued creation-time CI may not be promoted to scientific PASS")

    integration = contract.get("active_integration_snapshot", {})
    if integration.get("branch") != "codex/cr001-constraints":
        errors.append("active integration branch drift")
    if integration.get("live_head") != EXPECTED_INTEGRATION_HEAD:
        errors.append("active integration live head drift")
    if "live git ref is authoritative" not in integration.get("note", ""):
        errors.append("snapshot-vs-live-ref rule must remain explicit")

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
    domain = constraints.get("domain", {})
    if domain.get("physical") != "R^3":
        errors.append("canonical physical domain drift")
    if domain.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]]:
        errors.append("canonical evaluation box drift")
    if domain.get("support") != "r < 2 and abs(z) < 2":
        errors.append("canonical support drift")
    if domain.get("time_interval") != [0.25, 0.75]:
        errors.append("canonical time interval drift")
    if constraints.get("forcing", {}).get("mode") != "restricted_two_parameter_family":
        errors.append("canonical forcing mode drift")
    restriction = constraints.get("forcing", {}).get("restriction", "")
    if "No residual-dependent basis or pointwise free force" not in restriction:
        errors.append("canonical free-force prohibition drift")
    nontrivial = constraints.get("nontriviality", {})
    if nontrivial.get("reference_energy") != 1.0:
        errors.append("canonical reference energy drift")
    if nontrivial.get("reference_energy_abs_tolerance") != 0.001:
        errors.append("canonical reference-energy tolerance drift")
    validation = constraints.get("validation", {})
    if validation.get("seed") != 914027 or validation.get("held_out_points") != 4096:
        errors.append("canonical validation split drift")
    if validation.get("derivative_steps") != [0.02, 0.01, 0.005]:
        errors.append("canonical derivative ladder drift")
    if validation.get("quadrature_orders_per_axis") != [24, 48, 96]:
        errors.append("canonical quadrature ladder drift")
    thresholds = validation.get("thresholds", {})
    if thresholds.get("pde_residual_max") != 0.001 or thresholds.get("pde_residual_L2") != 0.001:
        errors.append("canonical momentum thresholds drift")
    if thresholds.get("divergence_max") != 1e-5 or thresholds.get("divergence_L2") != 1e-5:
        errors.append("canonical divergence thresholds drift")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fail-closed CR002 audit for pre-I1 autonomous modulation provenance and candidate-identity scope."
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()

    errors = audit_repository(args.repo_root.resolve())
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS: CR002 pre-I1 autonomous modulation scope is intact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
