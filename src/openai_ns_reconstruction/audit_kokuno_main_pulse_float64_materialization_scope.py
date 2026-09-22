"""Fail-closed CR002 audit for the Kokuno main-pulse float64 materialization boundary."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping

CONTRACT_PATH = "configs/kokuno_main_pulse_float64_materialization_scope.json"
CANDIDATE_PATH = "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_main_pulse.py"
CONSTRAINTS_PATH = "configs/constraints.json"
PROJECT_STATUS_PATH = "project_status.json"

BASE_EXACT_HEAD = "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156"
CANDIDATE_BLOB_SHA = "efcce36612ffc1cb3a74f5eb9c9ae6c86748f60d"
CONSTRAINTS_BLOB_SHA = "6c559e42895a606e2ef025ade4cb448966d75814"


def _expect(condition: bool, message: str, violations: list[str]) -> None:
    if not condition:
        violations.append(message)


def _as_dict(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError("expected mapping")
    return dict(value)


def float64_coordinate_witness() -> dict[str, Any]:
    """Synthetic representation witness; not source/OpenAI/PDE numeric evidence."""
    X_p = 1.0e304
    lam = 1.0
    source_xi_end = 11.0
    log_float_max = math.log(sys.float_info.max)
    log_X_p = math.log(X_p)
    xi_float_max = lam * (log_float_max - log_X_p)
    source_endpoint_log_X = log_X_p + source_xi_end / lam
    return {
        "classification": "autonomous_mechanics_only",
        "X_p": X_p,
        "lambda": lam,
        "source_xi_end": source_xi_end,
        "xi_float_max": xi_float_max,
        "source_coordinate_endpoint_is_finite": math.isfinite(source_xi_end),
        "source_endpoint_requires_X_beyond_float64": source_endpoint_log_X > log_float_max,
        "float64_prefix_ends_before_source_xi_end": 0.0 < xi_float_max < source_xi_end,
        "not_source_numeric_evidence": True,
        "not_openai_numeric_evidence": True,
        "not_pde_evidence": True,
    }


def audit_contract(
    contract: Mapping[str, Any],
    candidate_text: str,
    constraints: Mapping[str, Any],
    project_status: Mapping[str, Any],
) -> list[str]:
    """Return violations; the governed baseline is admitted only when empty."""
    violations: list[str] = []
    cfg = _as_dict(contract)
    scope = _as_dict(cfg.get("scope", {}))
    classes = _as_dict(cfg.get("provenance_classes", {}))
    truth = _as_dict(cfg.get("current_truth", {}))
    rep = _as_dict(cfg.get("representation_contract", {}))
    promo = _as_dict(cfg.get("promotion_requirements", {}))
    witness_decl = _as_dict(cfg.get("mechanics_witness", {}))
    cr001 = _as_dict(cfg.get("cr001_lock", {}))
    eq45 = _as_dict(cfg.get("canonical_eq45_independent_state", {}))

    _expect(
        cfg.get("schema") == "cr002-kokuno-main-pulse-float64-materialization-scope-v1",
        "unexpected contract schema",
        violations,
    )
    _expect(
        cfg.get("task_id") == "CR002-KOKUNO-MAIN-PULSE-FLOAT64-MATERIALIZATION-SCOPE-130",
        "unexpected task id",
        violations,
    )
    _expect(scope.get("base_pr") == 1100, "base PR drifted", violations)
    _expect(scope.get("base_exact_head") == BASE_EXACT_HEAD, "base exact head drifted", violations)
    _expect(scope.get("audited_candidate_path") == CANDIDATE_PATH, "candidate path drifted", violations)
    _expect(scope.get("audited_candidate_blob_sha") == CANDIDATE_BLOB_SHA, "candidate blob drifted", violations)
    for name in (
        "candidate_bytes_changed",
        "source_formula_changed",
        "pressure_changed",
        "forcing_changed",
        "scientific_threshold_changed",
    ):
        _expect(scope.get(name) is False, f"governance scope unexpectedly changes {name}", violations)
    _expect(scope.get("change_type") == "governance_only", "not governance-only", violations)

    _expect(
        set(classes) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "four-way provenance classes changed",
        violations,
    )
    for key in classes:
        _expect(
            isinstance(classes.get(key), list) and bool(classes.get(key)),
            f"provenance class {key} is empty/non-list",
            violations,
        )
    source_facts = "\n".join(str(v).lower() for v in classes.get("public_source_fact", []))
    autonomous = "\n".join(str(v).lower() for v in classes.get("autonomous_design", []))
    pending = "\n".join(str(v).lower() for v in classes.get("pending_unknown", []))
    _expect("xi=11" in source_facts, "public source-coordinate endpoint classification missing", violations)
    _expect(
        "not evidence" in source_facts and "float64" in source_facts,
        "public-source class must deny source inference from float64 ceiling",
        violations,
    )
    _expect(
        "binary64" in autonomous or "float64" in autonomous,
        "float64 materialization choice not classified autonomous",
        violations,
    )
    _expect(
        "remaining source-coordinate interval" in pending,
        "missing prefix completion not kept pending",
        violations,
    )

    expected_true = (
        "source_coordinate_main_pulse_kernel_defined_through_xi_11",
        "current_cartesian_main_pulse_finite_X_prefix_materialized",
        "current_cartesian_main_pulse_velocity_callable_on_supported_prefix",
        "current_cartesian_main_pulse_configuration_serializable",
        "float64_materialization_ceiling_occurs_before_xi_11_for_frozen_default",
    )
    expected_false = (
        "float64_materialization_ceiling_is_public_source_support_endpoint",
        "float64_materialization_ceiling_is_physical_domain_boundary",
        "full_source_xi_11_current_cartesian_materialized",
        "terminal_global_velocity_materialized",
        "velocity_export_ready_for_kokuno_route",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in expected_true:
        _expect(truth.get(key) is True, f"{key} must remain true", violations)
    for key in expected_false:
        _expect(truth.get(key) is False, f"{key} must remain false", violations)

    _expect(rep.get("current_numeric_coordinate") == "finite float64 X", "current numeric coordinate boundary changed", violations)
    _expect(rep.get("source_coordinate") == "xi=lambda*log(X/X_p)", "source coordinate changed", violations)
    _expect(rep.get("frozen_default_relation") == "8 < xi_materializable_max < 11", "frozen default materialization relation changed", violations)
    _expect(rep.get("ceiling_classification") == "autonomous_numerical_representation_limit", "float64 ceiling misclassified", violations)
    _expect(rep.get("source_endpoint_classification") == "public_source_coordinate_domain", "source endpoint classification changed", violations)
    _expect(
        "new representation/candidate identity" in str(rep.get("future_representation_identity_rule", "")),
        "future representation identity firewall missing",
        violations,
    )
    shortcut = str(rep.get("forbidden_shortcut", "")).lower()
    _expect(
        all(term in shortcut for term in ("source support", "physical support", "global velocity")),
        "forbidden source/support/global shortcut incomplete",
        violations,
    )

    _expect(
        isinstance(promo.get("full_source_xi_11_current_cartesian_materialized"), list)
        and len(promo["full_source_xi_11_current_cartesian_materialized"]) >= 5,
        "xi=11 promotion requirements incomplete",
        violations,
    )
    _expect(
        isinstance(promo.get("velocity_export_ready_for_kokuno_route"), list)
        and len(promo["velocity_export_ready_for_kokuno_route"]) >= 4,
        "Kokuno export promotion requirements incomplete",
        violations,
    )

    implications = cfg.get("forbidden_implications")
    _expect(
        isinstance(implications, list) and len(implications) >= 10,
        "forbidden implication firewall incomplete",
        violations,
    )
    joined_implications = "\n".join(str(v).lower() for v in implications or [])
    for phrase in (
        "source-coordinate kernel defined through xi=11 -> float64 cartesian x can materialize xi=11",
        "float64 cartesian materialization ceiling -> public source pulse terminates there",
        "finite-x prefix callable -> terminal/global velocity materialized",
        "ci success -> scientific admission",
    ):
        _expect(phrase in joined_implications, f"missing forbidden implication: {phrase}", violations)

    _expect(witness_decl.get("classification") == "autonomous_mechanics_only", "mechanics witness classification changed", violations)
    for key in ("not_source_numeric_evidence", "not_openai_numeric_evidence", "not_pde_evidence"):
        _expect(witness_decl.get(key) is True, f"mechanics witness lost {key}", violations)
    witness = float64_coordinate_witness()
    _expect(witness["source_endpoint_requires_X_beyond_float64"] is True, "mechanics witness became vacuous", violations)
    _expect(
        witness["float64_prefix_ends_before_source_xi_end"] is True,
        "mechanics witness no longer separates coordinate/materialization domains",
        violations,
    )

    required_source_tokens = (
        'SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-v1"',
        "float_end = math.log(np.finfo(float).max / _FLOAT_X_SAFETY_DIVISOR)",
        "return min(source_end, float_end)",
        "return self.lambda_value * (self.log_X_materializable_end - self.log_X_p)",
        "self.xi_materializable_max < MAIN_XI_MAX",
        '"full_source_xi_11_current_cartesian_materialized": False',
        '"unified_global_cartesian_velocity_export_ready": False',
        'raise ValueError("Cartesian point lies beyond the finite-X current main-pulse domain")',
        "def velocity(self, x: Any, y: Any, z: Any, t: Any)",
    )
    for token in required_source_tokens:
        _expect(token in candidate_text, f"exact #1100 implementation token missing: {token}", violations)

    _expect(cr001.get("constraints_blob_sha") == CONSTRAINTS_BLOB_SHA, "canonical constraints blob pin changed", violations)
    expected_cr001 = {
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
    }
    for key, expected in expected_cr001.items():
        _expect(cr001.get(key) == expected, f"CR001 contract drift: {key}", violations)
    for key in (
        "residual_defined_free_forcing_forbidden",
        "candidate_collapse_forbidden",
        "post_hoc_threshold_relaxation_forbidden",
    ):
        _expect(cr001.get(key) is True, f"CR001 prohibition lost: {key}", violations)

    c = _as_dict(constraints)
    domain = _as_dict(c.get("domain", {}))
    forcing = _as_dict(c.get("forcing", {}))
    nontrivial = _as_dict(c.get("nontriviality", {}))
    validation = _as_dict(c.get("validation", {}))
    thresholds = _as_dict(validation.get("thresholds", {}))
    _expect(c.get("nu") == 0.01, "live constraints viscosity drifted", violations)
    _expect(domain.get("physical") == "R^3", "live physical domain drifted", violations)
    _expect(domain.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "live evaluation box drifted", violations)
    _expect(domain.get("support") == "r < 2 and abs(z) < 2", "live support drifted", violations)
    _expect(domain.get("time_interval") == [0.25, 0.75], "live time interval drifted", violations)
    _expect(forcing.get("mode") == "restricted_two_parameter_family", "live forcing mode drifted", violations)
    _expect(
        "No residual-dependent basis or pointwise free force" in str(forcing.get("restriction", "")),
        "free residual-defined forcing prohibition drifted",
        violations,
    )
    _expect(nontrivial.get("reference_energy") == 1.0, "live reference energy drifted", violations)
    _expect(nontrivial.get("reference_energy_abs_tolerance") == 0.001, "live energy tolerance drifted", violations)
    _expect(validation.get("seed") == 914027 and validation.get("held_out_points") == 4096, "live validation identity drifted", violations)
    _expect(validation.get("derivative_steps") == [0.02, 0.01, 0.005], "live FD ladder drifted", violations)
    _expect(validation.get("quadrature_orders_per_axis") == [24, 48, 96], "live quadrature ladder drifted", violations)
    _expect(
        thresholds.get("pde_residual_max") == 0.001 and thresholds.get("pde_residual_L2") == 0.001,
        "live momentum thresholds drifted",
        violations,
    )
    _expect(
        thresholds.get("divergence_max") == 1e-5 and thresholds.get("divergence_L2") == 1e-5,
        "live divergence thresholds drifted",
        violations,
    )

    _expect(eq45.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "Eq45 family pin changed", violations)
    _expect(eq45.get("velocity_export_ready") is True, "Eq45 export readiness was incorrectly downgraded", violations)
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _expect(eq45.get(key) is False, f"Eq45 independent state incorrectly promoted: {key}", violations)

    p = _as_dict(project_status)
    states = _as_dict(p.get("states", {}))
    _expect(p.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "live project canonical candidate changed", violations)
    _expect(
        p.get("velocity_api") == "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "live canonical velocity API changed",
        violations,
    )
    _expect(states.get("velocity_export_ready") is True, "live Eq45 velocity export readiness changed", violations)
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _expect(states.get(key) is False, f"live project state promoted unexpectedly: {key}", violations)

    return violations


def load_and_audit(repo_root: Path) -> list[str]:
    contract = json.loads((repo_root / CONTRACT_PATH).read_text(encoding="utf-8"))
    candidate_text = (repo_root / CANDIDATE_PATH).read_text(encoding="utf-8")
    constraints = json.loads((repo_root / CONSTRAINTS_PATH).read_text(encoding="utf-8"))
    project_status = json.loads((repo_root / PROJECT_STATUS_PATH).read_text(encoding="utf-8"))
    return audit_contract(contract, candidate_text, constraints, project_status)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    violations = load_and_audit(args.repo_root)
    if violations:
        for violation in violations:
            print(f"VIOLATION: {violation}")
        return 1
    print("CR002 Kokuno main-pulse float64 materialization scope: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
