"""Fail-closed CR002 audit for Kokuno pressure-increment gradient scope.

This module does not construct pressure or alter a velocity candidate.  It locks a
representation boundary exposed by the Agent-5 #911 ingest seam: a registered
pressure increment C_p=Pi-Pi_0(eta), even with independently checked C_p
 derivatives, is not yet a matched Cartesian pressure gradient while the
eta-dependent axis datum Pi_0(eta) and the eventual coordinate pullback are
missing.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


SCOPE_PATH = Path("configs/kokuno_pressure_increment_gradient_scope.json")
A5_SOURCE_PATH = Path(
    "src/openai_ns_reconstruction/kokuno_a5_reference_pressure_increment_ingest_contract.py"
)
CONSTRAINTS_PATH = Path("configs/constraints.json")

EXPECTED_A5_HEAD = "dd613c3636391f23a8fc8958af6f17d024482f5f"
EXPECTED_A5_SOURCE_BLOB = "c8c23490c62c981440e7e725f8345335df3e59f9"
EXPECTED_A1_SOURCE_BLOB = "4704fed4271da2c6070110a6935d2f33f7932f56"
EXPECTED_A4_SOURCE_BLOB = "027abf88fa1a916fcde3c6d818b19a9ff143c445"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
PROVENANCE_KEYS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}

EXPECTED_PUBLIC_FACTS = [
    "C_p_reference(X,eta)=Pi_r(X,eta)-Pi_0(eta)",
    "partial_X C_p_reference=F_r^2=E_r^2/(2X)",
    "Pi_0(eta) is an outer-construction datum and is not determined by the local radial pressure-increment identity alone",
]
EXPECTED_AUTONOMOUS = [
    "calculus identity grad p = C_p,X grad X + (C_p,eta + Pi_0'(eta)) grad eta for a two-coordinate pressure representation",
    "mechanics witness X=x, eta=z with p_A=X^2 and p_B=X^2+0.3 eta",
    "mutation regressions that reject pressure-increment-to-grad-p claim promotion",
]
EXPECTED_PENDING = [
    "numerical materialization of the source-consistent Pi_0(eta)",
    "Pi_0'(eta) on the eventual source-consistent pressure representation",
    "the full source-to-Cartesian coordinate pullback needed for matched global grad p",
    "matched/global pressure on a complete Kokuno candidate",
    "complete Kokuno velocity-pressure-forcing API and full-NS validation",
]
EXPECTED_RULES = {
    "pressure_increment_identity_is_absolute_pressure_identity": False,
    "pressure_increment_derivative_audit_authorizes_cartesian_grad_p": False,
    "eta_dependent_axis_datum_is_gradient_gauge": False,
    "spatially_constant_pressure_addition_leaves_spatial_gradient_unchanged": True,
    "canonical_pressure_gauge_and_boundary_still_must_be_satisfied": True,
    "Pi0_eta_required_for_full_eta_pressure_derivative": True,
    "coordinate_pullback_required_for_cartesian_grad_p": True,
    "matched_global_pressure_required_before_complete_ns_residual": True,
    "reference_increment_failure_blocks_existing_callable_eq45_delivery": False,
}
EXPECTED_TRUTH = {
    "kokuno_pressure_increment_registered": True,
    "kokuno_absolute_axis_pressure_materialized": False,
    "kokuno_matched_global_pressure_materialized": False,
    "kokuno_cartesian_pressure_gradient_ready": False,
    "kokuno_complete_candidate_api_ready": False,
    "kokuno_pde_validated": False,
    "visual_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "canonical_eq45_velocity_delivery_remains_independent": True,
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def mechanics_witness() -> dict[str, Any]:
    """Return a local calculus witness; this is not Kokuno candidate evidence."""
    alpha = 0.3
    points = [(0.0, -0.7), (0.2, 0.1), (0.9, 0.4)]
    max_increment_gap = 0.0
    max_increment_dx_gap = 0.0
    max_increment_deta_gap = 0.0
    max_gradient_gap = 0.0
    max_constant_gauge_gradient_gap = 0.0

    for x, eta in points:
        # Autonomous coordinate choice for this mechanics witness only: X=x, eta=z.
        p_a = x * x
        p_b = x * x + alpha * eta
        pi0_a = 0.0
        pi0_b = alpha * eta
        cp_a = p_a - pi0_a
        cp_b = p_b - pi0_b
        max_increment_gap = max(max_increment_gap, abs(cp_a - cp_b))

        cp_x_a = 2.0 * x
        cp_x_b = 2.0 * x
        cp_eta_a = 0.0
        cp_eta_b = 0.0
        max_increment_dx_gap = max(max_increment_dx_gap, abs(cp_x_a - cp_x_b))
        max_increment_deta_gap = max(max_increment_deta_gap, abs(cp_eta_a - cp_eta_b))

        grad_a = (2.0 * x, 0.0, 0.0)
        grad_b = (2.0 * x, 0.0, alpha)
        max_gradient_gap = max(
            max_gradient_gap,
            math.sqrt(sum((a - b) ** 2 for a, b in zip(grad_a, grad_b))),
        )

        # By contrast, adding a spatial constant changes p but not grad p.
        grad_c = (2.0 * x, 0.0, 0.0)
        max_constant_gauge_gradient_gap = max(
            max_constant_gauge_gradient_gap,
            math.sqrt(sum((a - c) ** 2 for a, c in zip(grad_a, grad_c))),
        )

    _require(max_increment_gap == 0.0, "mechanics witness increments must agree exactly")
    _require(max_increment_dx_gap == 0.0, "mechanics witness dX increments must agree")
    _require(max_increment_deta_gap == 0.0, "mechanics witness deta increments must agree")
    _require(abs(max_gradient_gap - alpha) <= 1.0e-15, "eta datum must change Cartesian gradient")
    _require(
        max_constant_gauge_gradient_gap == 0.0,
        "spatially constant pressure addition must leave gradient unchanged",
    )
    return {
        "alpha": alpha,
        "points": [[x, eta] for x, eta in points],
        "max_pressure_increment_gap": max_increment_gap,
        "max_pressure_increment_X_derivative_gap": max_increment_dx_gap,
        "max_pressure_increment_eta_derivative_gap": max_increment_deta_gap,
        "max_cartesian_gradient_gap": max_gradient_gap,
        "max_spatial_constant_gradient_gap": max_constant_gauge_gradient_gap,
        "candidate_evidence": False,
        "kokuno_or_openai_parameter_evidence": False,
    }


def _validate_canonical_constraints(scope: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    inv = scope.get("cr001_invariants")
    _require(isinstance(inv, Mapping), "cr001_invariants must be a mapping")
    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})

    _require(constraints.get("nu") == inv.get("nu") == 0.01, "nu drift")
    _require(domain.get("physical") == inv.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == inv.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drift")
    _require(domain.get("support") == inv.get("support") == "r < 2 and abs(z) < 2", "support drift")
    _require(domain.get("time_interval") == inv.get("time_interval") == [0.25, 0.75], "time interval drift")
    _require(domain.get("pressure_gauge") == inv.get("pressure_gauge") == "p=0 outside compact support", "pressure gauge drift")
    _require(forcing.get("mode") == inv.get("forcing_mode") == "restricted_two_parameter_family", "forcing family drift")
    _require(inv.get("residual_defined_pointwise_free_force_forbidden") is True, "free-force firewall disabled")
    _require("No residual-dependent basis or pointwise free force" in forcing.get("restriction", ""), "canonical free-force restriction drift")
    _require(nontriviality.get("reference_energy") == inv.get("reference_energy") == 1.0, "reference energy drift")
    _require(nontriviality.get("reference_energy_abs_tolerance") == inv.get("reference_energy_abs_tolerance") == 0.001, "energy tolerance drift")
    _require("reject collapsed candidates" in nontriviality.get("enforcement", ""), "amplitude-collapse firewall drift")
    _require(inv.get("amplitude_collapse_forbidden") is True, "amplitude-collapse governance disabled")
    _require(validation.get("seed") == inv.get("validation_seed") == 914027, "validation seed drift")
    _require(validation.get("held_out_points") == inv.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("derivative_steps") == inv.get("derivative_steps") == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(validation.get("quadrature_orders_per_axis") == inv.get("quadrature_orders_per_axis") == [24, 48, 96], "quadrature ladder drift")
    _require(thresholds.get("pde_residual_max") == inv.get("momentum_sampled_max_gate") == 0.001, "momentum max gate drift")
    _require(thresholds.get("pde_residual_L2") == inv.get("momentum_volume_l2_gate") == 0.001, "momentum L2 gate drift")
    _require(thresholds.get("divergence_max") == inv.get("divergence_sampled_max_gate") == 1.0e-5, "divergence max gate drift")
    _require(thresholds.get("divergence_L2") == inv.get("divergence_volume_l2_gate") == 1.0e-5, "divergence L2 gate drift")
    _require(thresholds.get("boundary_pressure_max") == inv.get("boundary_pressure_max_gate") == 1.0e-10, "boundary pressure gate drift")
    _require(inv.get("post_hoc_threshold_relaxation_forbidden") is True, "threshold-relaxation firewall disabled")
    _require("changing thresholds requires a new experiment version" in validation.get("failure_policy", ""), "canonical threshold-change policy drift")


def validate_scope(scope: Mapping[str, Any], *, repo_root: Path | None = None) -> dict[str, Any]:
    """Validate the governance payload and the exact local upstream identities."""
    root = _repo_root() if repo_root is None else Path(repo_root)
    c = copy.deepcopy(dict(scope))
    _require(c.get("schema_version") == 1, "schema version drift")
    _require(c.get("task_id") == "CR002-KOKUNO-PRESSURE-INCREMENT-GRADIENT-SCOPE-085", "task id drift")

    upstream = c.get("audited_upstream", {})
    _require(upstream.get("agent5_pr") == 911, "A5 PR drift")
    _require(upstream.get("agent5_head") == EXPECTED_A5_HEAD, "A5 head drift")
    _require(upstream.get("agent5_source_blob") == EXPECTED_A5_SOURCE_BLOB, "A5 source-blob registration drift")
    _require(upstream.get("agent1_source_blob") == EXPECTED_A1_SOURCE_BLOB, "A1 pressure source-blob drift")
    _require(upstream.get("agent4_source_blob") == EXPECTED_A4_SOURCE_BLOB, "A4 audit source-blob drift")
    _require(upstream.get("canonical_constraints_blob") == EXPECTED_CONSTRAINTS_BLOB, "canonical constraints registration drift")

    a5_path = root / A5_SOURCE_PATH
    constraints_path = root / CONSTRAINTS_PATH
    _require(a5_path.is_file(), f"missing audited A5 source: {a5_path}")
    _require(constraints_path.is_file(), f"missing canonical constraints: {constraints_path}")
    _require(_git_blob_sha(a5_path) == EXPECTED_A5_SOURCE_BLOB, "local A5 source no longer matches audited blob")
    _require(_git_blob_sha(constraints_path) == EXPECTED_CONSTRAINTS_BLOB, "local canonical constraints no longer match audited blob")

    a5_text = a5_path.read_text(encoding="utf-8")
    for snippet in (
        '"pressure_increment_authorized_as_cartesian_grad_p": False',
        '"Pi0_required_before_absolute_pressure": True',
        '"cartesian_pressure_gradient_materialized": False',
        '"absolute_axis_pressure_Pi0_materialized": False',
        '"source_blob": "4704fed4271da2c6070110a6935d2f33f7932f56"',
        '"source_blob": "027abf88fa1a916fcde3c6d818b19a9ff143c445"',
    ):
        _require(snippet in a5_text, f"A5 pressure truth-boundary snippet missing: {snippet}")

    provenance = c.get("provenance", {})
    _require(set(provenance) == PROVENANCE_KEYS, "provenance classes must remain exactly the canonical four")
    _require(provenance.get("public_source_fact") == EXPECTED_PUBLIC_FACTS, "public-source fact classification drift")
    _require(provenance.get("autonomous_design") == EXPECTED_AUTONOMOUS, "autonomous mechanics classification drift")
    _require(provenance.get("pending_unknown") == EXPECTED_PENDING, "pending/unknown classification drift")
    classified = [str(item) for key in PROVENANCE_KEYS for item in provenance.get(key, [])]
    _require(len(classified) == len(set(classified)), "provenance classes must not overlap")

    _require(c.get("representation_rules") == EXPECTED_RULES, "pressure representation rule drift/promotion")
    witness_cfg = c.get("mechanics_witness", {})
    _require(witness_cfg.get("classification") == "autonomous_design", "mechanics witness must remain autonomous")
    _require(witness_cfg.get("alpha") == 0.3, "mechanics witness alpha drift")
    _require(witness_cfg.get("candidate_evidence") is False, "mechanics witness cannot become candidate evidence")
    _require(witness_cfg.get("kokuno_or_openai_parameter_evidence") is False, "mechanics witness cannot become source parameter evidence")

    constraints = _load_json(constraints_path)
    _validate_canonical_constraints(c, constraints)
    _require(c.get("truth_state") == EXPECTED_TRUTH, "truth-state promotion/drift")

    witness = mechanics_witness()
    return {
        "ok": True,
        "agent5_source_blob": _git_blob_sha(a5_path),
        "canonical_constraints_blob": _git_blob_sha(constraints_path),
        "mechanics_witness": witness,
    }


def audit_registered_scope(*, repo_root: Path | None = None) -> dict[str, Any]:
    root = _repo_root() if repo_root is None else Path(repo_root)
    return validate_scope(_load_json(root / SCOPE_PATH), repo_root=root)


if __name__ == "__main__":
    print(json.dumps(audit_registered_scope(), indent=2, sort_keys=True))
