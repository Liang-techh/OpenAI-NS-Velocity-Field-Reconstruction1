"""Fail-closed CR002 audit for autonomous Kokuno axis-pressure identity scope.

This audit governs one narrow representation boundary: satisfying the displayed
axis-pressure admissibility conditions, and independently validating one chosen
repository realization, do not identify that realization with the
source-prepared Appendix-A ``Pi_0`` datum.

The numerical witness in this module is mechanics-only governance evidence.  It
is not a Kokuno/OpenAI parameter claim and it is not candidate evidence.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


CONFIG_REL = Path("configs/kokuno_autonomous_axis_pressure_identity_scope.json")
A1_REL = Path(
    "src/openai_ns_reconstruction/kokuno_pa10_autonomous_axis_pressure_seed_x100.py"
)
A4_REL = Path(
    "src/openai_ns_reconstruction/kokuno_a4_autonomous_reference_pressure_independent_audit.py"
)
CONSTRAINTS_REL = Path("configs/constraints.json")

EXPECTED_SCHEMA = "cr002-kokuno-autonomous-axis-pressure-identity-scope-v1"
EXPECTED_SCOPE_ID = "CR002-KOKUNO-AXIS-PRESSURE-ADMISSIBILITY-IDENTITY-086"
EXPECTED_A1_BLOB = "68dddc4f6ad1fb5bb6955498099b380240be6750"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_A1_HEAD = "a5903ce594bf1c1f775d6fbba8a5eeca226622fa"
EXPECTED_A4_HEAD = "a793bffd4739de834a81e2ed32fdbf9821dad04e"
EXPECTED_INTEGRATION_HEAD = "ca37d25d19b97d893030194ebd6364160ae4355e"

_PROVENANCE_KEYS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}

_REPRESENTATION_FLAGS = {
    "admissibility_conditions_satisfied_does_not_imply_source_identity",
    "envelope_saturation_is_autonomous_selection_not_public_source_equality",
    "implementation_distinct_audit_can_validate_realization_consistency_only",
    "implementation_distinct_audit_does_not_upgrade_provenance",
    "derived_quantities_using_autonomous_Pi0_remain_autonomous_unless_rebound_to_source_prepared_data",
    "source_prepared_Pi0_identity_requires_independent_source_binding",
    "pressure_scalar_identity_is_separate_from_matched_cartesian_gradient_identity",
}

_FALSE_TRUTH_FLAGS = {
    "source_prepared_appendixA_Pi0_materialized",
    "source_prepared_appendixA_Pi0_identified_with_autonomous_seed",
    "source_correspondence_verified",
    "matched_global_pressure_materialized",
    "cartesian_matched_pressure_gradient_materialized",
    "complete_velocity_pressure_forcing_candidate_ready",
    "kokuno_velocity_export_ready",
    "pde_validated",
    "visual_correspondence_verified",
    "paper_exact",
    "openai_field_identified",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha(raw: bytes) -> str:
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def mechanics_witness(*, epsilon: float = 0.1) -> dict[str, Any]:
    """Return a deterministic non-uniqueness witness for the admissibility data."""
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")

    etas = (-0.8, -0.4, -0.2, 0.2, 0.4, 0.8)
    rows: list[dict[str, float]] = []
    all_envelope = True
    all_seed_sign = True
    all_alt_sign = True
    max_profile_gap = 0.0
    max_derivative_gap = 0.0

    for eta in etas:
        q = 1.0 + eta * eta
        seed = -2.5 / (q * q)
        alt = seed - epsilon / (q * q * q)
        seed_eta = 10.0 * eta / (q * q * q)
        alt_eta = seed_eta + 6.0 * epsilon * eta / (q * q * q * q)
        all_envelope = all_envelope and alt <= seed + 1.0e-15
        all_seed_sign = all_seed_sign and eta * seed_eta > 0.0
        all_alt_sign = all_alt_sign and eta * alt_eta > 0.0
        max_profile_gap = max(max_profile_gap, abs(alt - seed))
        max_derivative_gap = max(max_derivative_gap, abs(alt_eta - seed_eta))
        rows.append(
            {
                "eta": eta,
                "seed": seed,
                "alternative": alt,
                "seed_eta": seed_eta,
                "alternative_eta": alt_eta,
            }
        )

    even_seed = all(
        math.isclose(row["seed"], rows[-idx - 1]["seed"], rel_tol=0.0, abs_tol=1.0e-15)
        for idx, row in enumerate(rows)
    )
    even_alt = all(
        math.isclose(
            row["alternative"],
            rows[-idx - 1]["alternative"],
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        for idx, row in enumerate(rows)
    )

    return {
        "epsilon": epsilon,
        "rows": rows,
        "alternative_obeys_seed_envelope": bool(all_envelope),
        "seed_derivative_sign_passes": bool(all_seed_sign),
        "alternative_derivative_sign_passes": bool(all_alt_sign),
        "seed_even_on_probe": bool(even_seed),
        "alternative_even_on_probe": bool(even_alt),
        "max_profile_gap": max_profile_gap,
        "max_derivative_gap": max_derivative_gap,
        "profiles_distinct": bool(max_profile_gap > 1.0e-6),
        "derivatives_distinct": bool(max_derivative_gap > 1.0e-6),
    }


def _require_text(text: str, token: str, label: str, errors: list[str]) -> None:
    if token not in text:
        errors.append(f"{label} is missing required token {token!r}")


def audit_scope(
    root: Path | None = None,
    *,
    config_override: Mapping[str, Any] | None = None,
    a1_text_override: str | None = None,
    a4_text_override: str | None = None,
) -> list[str]:
    """Return fail-closed errors for the registered scope; empty means audit pass."""
    base = Path(root) if root is not None else _repo_root()
    errors: list[str] = []

    config = (
        json.loads(json.dumps(config_override))
        if config_override is not None
        else _load_json(base / CONFIG_REL)
    )

    if config.get("schema") != EXPECTED_SCHEMA:
        errors.append("scope schema changed")
    if config.get("scope_id") != EXPECTED_SCOPE_ID:
        errors.append("scope_id changed")

    integration = config.get("integration_reference", {})
    if integration.get("branch") != "codex/cr001-constraints":
        errors.append("integration reference branch changed")
    if integration.get("head") != EXPECTED_INTEGRATION_HEAD:
        errors.append("integration reference head changed")

    stack = config.get("stack_base", {})
    if stack.get("agent1_head") != EXPECTED_A1_HEAD:
        errors.append("Agent-1 exact head binding changed")
    if stack.get("agent4_head") != EXPECTED_A4_HEAD:
        errors.append("Agent-4 exact head binding changed")
    if stack.get("agent1_axis_pressure_source_blob") != EXPECTED_A1_BLOB:
        errors.append("Agent-1 source-blob binding changed")
    if stack.get("canonical_constraints_blob") != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("canonical constraints-blob binding changed")

    provenance = config.get("provenance", {})
    if set(provenance) != _PROVENANCE_KEYS:
        errors.append("the four provenance classes are not exactly separated")
    public_text = "\n".join(str(x) for x in provenance.get("public_source_fact", [])).lower()
    for forbidden in ("pi0_seed", "p_*=1", "epsilon=0.1", "saturat"):
        if forbidden in public_text:
            errors.append(f"autonomous choice leaked into public_source_fact: {forbidden}")
    autonomous_text = "\n".join(str(x) for x in provenance.get("autonomous_design", [])).lower()
    if "p_*=1" not in autonomous_text:
        errors.append("autonomous P_*=1 selection is not classified as autonomous_design")
    if "pi0_seed" not in autonomous_text or "saturat" not in autonomous_text:
        errors.append("autonomous envelope-saturating Pi0 seed is not classified as autonomous_design")
    pending_text = "\n".join(str(x) for x in provenance.get("pending_unknown", [])).lower()
    if "source-prepared" not in pending_text or "unknown" not in pending_text:
        errors.append("source-prepared Pi_0 identity is not preserved as pending/unknown")

    representation = config.get("representation_contract", {})
    if set(representation) != _REPRESENTATION_FLAGS:
        errors.append("representation-contract keys changed")
    for flag in sorted(_REPRESENTATION_FLAGS):
        if representation.get(flag) is not True:
            errors.append(f"representation contract weakened: {flag}")

    witness_config = config.get("mechanics_witness", {})
    if witness_config.get("role") != "autonomous_nonuniqueness_witness_not_candidate_evidence":
        errors.append("mechanics witness role changed")
    if witness_config.get("P_star") != 1.0 or witness_config.get("epsilon") != 0.1:
        errors.append("mechanics witness constants changed")
    witness = mechanics_witness(epsilon=0.1)
    for flag in (
        "alternative_obeys_seed_envelope",
        "seed_derivative_sign_passes",
        "alternative_derivative_sign_passes",
        "seed_even_on_probe",
        "alternative_even_on_probe",
        "profiles_distinct",
        "derivatives_distinct",
    ):
        if witness.get(flag) is not True:
            errors.append(f"non-uniqueness witness failed: {flag}")

    cr001 = config.get("cr001_invariants", {})
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
        "momentum_volume_l2": 0.001,
        "divergence_max": 0.00001,
        "divergence_volume_l2": 0.00001,
        "residual_defined_pointwise_free_force_allowed": False,
        "amplitude_collapse_allowed": False,
        "post_hoc_threshold_relaxation_allowed": False,
    }
    if cr001 != expected_cr001:
        errors.append("CR001 frozen invariants changed")

    truth = config.get("truth_boundary", {})
    if truth.get("repository_autonomous_axis_pressure_seed_materialized") is not True:
        errors.append("repository-autonomous seed materialization was incorrectly removed")
    if truth.get("autonomous_reference_pressure_numerically_auditable") is not True:
        errors.append("autonomous pressure auditability was incorrectly removed")
    if truth.get("canonical_eq45_velocity_export_ready_remains_independent") is not True:
        errors.append("Eq45 callable velocity delivery was incorrectly coupled to this pressure audit")
    for flag in sorted(_FALSE_TRUTH_FLAGS):
        if truth.get(flag) is not False:
            errors.append(f"unsupported truth promotion: {flag}")

    a1_path = base / A1_REL
    a1_text = a1_text_override if a1_text_override is not None else a1_path.read_text(encoding="utf-8")
    if a1_text_override is None and _git_blob_sha(a1_path.read_bytes()) != EXPECTED_A1_BLOB:
        errors.append("Agent-1 autonomous pressure source blob drifted")
    _require_text(a1_text, "Pi0_seed(eta) = -(5/2) P_*^2 (1+eta^2)^-2", "Agent-1 source", errors)
    _require_text(a1_text, '"source_prepared_appendixA_Pi0_materialized": False', "Agent-1 source", errors)
    _require_text(a1_text, '"matched_global_pressure_materialized": False', "Agent-1 source", errors)
    _require_text(a1_text, '"pde_validated": False', "Agent-1 source", errors)
    if "does **not** pretend to recover it" not in a1_text:
        errors.append("Agent-1 source no longer preserves the autonomous-vs-source identity disclaimer")

    a4_path = base / A4_REL
    a4_text = a4_text_override if a4_text_override is not None else a4_path.read_text(encoding="utf-8")
    _require_text(a4_text, f'A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB = "{EXPECTED_A1_BLOB}"', "Agent-4 audit", errors)
    _require_text(a4_text, '"source_prepared_appendixA_Pi0_materialized": False', "Agent-4 audit", errors)
    _require_text(a4_text, '"matched_global_pressure_materialized": False', "Agent-4 audit", errors)
    _require_text(a4_text, '"cartesian_pressure_gradient_assessed": False', "Agent-4 audit", errors)
    _require_text(a4_text, '"pde_validated": False', "Agent-4 audit", errors)
    if "repository-autonomous rather than the source-prepared" not in a4_text:
        errors.append("Agent-4 audit no longer preserves the autonomous-vs-source identity disclaimer")

    constraints_path = base / CONSTRAINTS_REL
    if _git_blob_sha(constraints_path.read_bytes()) != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("canonical constraints blob drifted from the registered CR001 base")

    return errors


def report(root: Path | None = None) -> dict[str, Any]:
    errors = audit_scope(root)
    return {
        "schema": EXPECTED_SCHEMA,
        "scope_id": EXPECTED_SCOPE_ID,
        "ok": not errors,
        "errors": errors,
        "mechanics_witness": mechanics_witness(),
    }


def _main() -> None:
    payload = report()
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0 if payload["ok"] else 1)


if __name__ == "__main__":
    _main()
