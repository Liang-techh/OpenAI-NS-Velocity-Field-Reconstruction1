"""Fail-closed CR002 audit for the Kokuno main-pulse amplitude boundary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

CONTRACT_PATH = "configs/kokuno_main_pulse_amplitude_scope.json"
SOURCE_PATH = "src/openai_ns_reconstruction/kokuno_source_main_axial_pulse_kernel.py"
CONSTRAINTS_PATH = "configs/constraints.json"
PROJECT_STATUS_PATH = "project_status.json"

BASE_EXACT_HEAD = "dc33914d3dd60c1dc4c8dfaab0760ef071385d54"
CONSTRAINTS_BLOB_SHA = "6c559e42895a606e2ef025ade4cb448966d75814"


def amplitude_dependence_witness() -> dict[str, Any]:
    """Autonomous mechanics witness; not source/OpenAI/PDE numeric evidence."""
    E = 2.0
    R0 = 0.3
    amplitude_a = 1.0
    amplitude_b = 1.125
    U_a = E * amplitude_a * R0
    U_b = E * amplitude_b * R0
    return {
        "classification": "autonomous_mechanics_only",
        "E": E,
        "R0": R0,
        "amplitude_a": amplitude_a,
        "amplitude_b": amplitude_b,
        "U_a": U_a,
        "U_b": U_b,
        "same_public_shape_inputs": True,
        "different_executable_U": U_a != U_b,
        "not_source_numeric_evidence": True,
        "not_openai_numeric_evidence": True,
        "not_pde_evidence": True,
    }


def _expect(condition: bool, message: str, violations: list[str]) -> None:
    if not condition:
        violations.append(message)


def _as_dict(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError("expected mapping")
    return dict(value)


def audit_contract(
    contract: Mapping[str, Any],
    source_text: str,
    constraints: Mapping[str, Any],
    project_status: Mapping[str, Any],
) -> list[str]:
    """Return violations; an empty list is the only admitted baseline."""
    violations: list[str] = []
    cfg = _as_dict(contract)
    scope = _as_dict(cfg.get("scope", {}))
    classes = _as_dict(cfg.get("provenance_classes", {}))
    truth = _as_dict(cfg.get("current_truth", {}))
    promo = _as_dict(cfg.get("promotion_requirements", {}))
    cr001 = _as_dict(cfg.get("cr001_lock", {}))
    eq45 = _as_dict(cfg.get("canonical_eq45_independent_state", {}))
    witness = _as_dict(cfg.get("mechanics_witness", {}))

    _expect(
        cfg.get("schema") == "cr002-kokuno-main-pulse-amplitude-scope-v1",
        "unexpected contract schema",
        violations,
    )
    _expect(
        scope.get("base_exact_head") == BASE_EXACT_HEAD,
        "base exact head drifted",
        violations,
    )
    for key in (
        "candidate_bytes_changed",
        "source_formula_changed",
        "pressure_changed",
        "forcing_changed",
        "scientific_threshold_changed",
    ):
        _expect(scope.get(key) is False, f"governance-only scope mutated: {key}", violations)

    required_buckets = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    _expect(set(classes) == required_buckets, "four-way provenance buckets changed", violations)
    for name in required_buckets:
        _expect(
            isinstance(classes.get(name), list) and classes[name],
            f"empty provenance bucket: {name}",
            violations,
        )

    public_text = "\n".join(str(x) for x in classes.get("public_source_fact", []))
    autonomous_text = "\n".join(str(x) for x in classes.get("autonomous_design", []))
    pending_text = "\n".join(str(x) for x in classes.get("pending_unknown", []))
    _expect(
        "A_principal=sqrt" not in public_text,
        "autonomous A_principal laundered into public-source facts",
        violations,
    )
    _expect(
        "A_principal=sqrt" in autonomous_text,
        "autonomous principal amplitude classification missing",
        violations,
    )
    _expect(
        "full source Amp(eta) root" in pending_text,
        "source-exact amplitude root no longer pending",
        violations,
    )

    exact_truth = {
        "source_main_pulse_kernel_formula_materialized": True,
        "repository_autonomous_principal_amplitude_materialized": True,
        "source_exact_amplitude_root_materialized": False,
        "source_pulse_end_MJ_corrections_materialized": False,
        "source_hidden_parameters_recovered": False,
        "current_cartesian_pulse_velocity_composed": False,
        "terminal_global_velocity_materialized": False,
        "velocity_export_ready_for_kokuno_route": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    _expect(truth == exact_truth, "current truth boundary was promoted or altered", violations)

    forbidden = cfg.get("forbidden_implications")
    _expect(
        isinstance(forbidden, list) and len(forbidden) >= 9,
        "forbidden implication firewall weakened",
        violations,
    )
    if isinstance(forbidden, list):
        joined = "\n".join(str(x) for x in forbidden)
        for token in (
            "public main pulse shape/formulas -> source-exact executable amplitude",
            "repository-autonomous principal amplitude -> recovered hidden/source Amp(eta)",
            "nonzero low-dimensional U profile -> current Cartesian pulse velocity",
            "later source-exact amplitude implementation -> inheritance of proxy-candidate numerical evidence without a new identity",
            "CI success -> scientific admission",
        ):
            _expect(token in joined, f"missing forbidden implication: {token}", violations)

    for key in (
        "source_exact_amplitude_root_materialized",
        "current_cartesian_pulse_velocity_composed",
        "velocity_export_ready_for_kokuno_route",
        "evidence_transfer_after_amplitude_replacement",
    ):
        _expect(
            isinstance(promo.get(key), list) and promo[key],
            f"missing promotion gate: {key}",
            violations,
        )
    transfer = "\n".join(
        str(x) for x in promo.get("evidence_transfer_after_amplitude_replacement", [])
    )
    _expect(
        "new candidate identity" in transfer,
        "amplitude replacement no longer requires a new identity",
        violations,
    )
    _expect(
        "do not transfer proxy-specific numerical receipts" in transfer,
        "proxy evidence-transfer firewall weakened",
        violations,
    )

    _expect(
        witness.get("classification") == "autonomous_mechanics_only",
        "mechanics witness misclassified",
        violations,
    )
    for key in (
        "not_source_numeric_evidence",
        "not_openai_numeric_evidence",
        "not_pde_evidence",
    ):
        _expect(
            witness.get(key) is True,
            f"mechanics witness evidence boundary weakened: {key}",
            violations,
        )

    # Exact #1094 implementation tokens: these are repository facts, not a
    # re-derivation of the external public source.
    source_tokens = (
        'SCHEMA = "kokuno-source-main-axial-pulse-kernel-v1"',
        '"principal_amplitude": (',
        '"repository-autonomous A_principal=sqrt((1-exp(-26))/(4 K_b)); "',
        '"source_exact_amplitude_root_materialized": False',
        '"current_cartesian_pulse_velocity_composed": False',
        "return math.sqrt((1.0 - math.exp(-26.0)) / (4.0 * kb))",
        "U = E * self.A_principal * R0",
        "MAIN_XI_MAX = 11.0",
    )
    for token in source_tokens:
        _expect(
            token in source_text,
            f"exact #1094 source token missing/drifted: {token}",
            violations,
        )

    c = _as_dict(constraints)
    domain = _as_dict(c.get("domain", {}))
    forcing = _as_dict(c.get("forcing", {}))
    nontriviality = _as_dict(c.get("nontriviality", {}))
    validation = _as_dict(c.get("validation", {}))
    thresholds = _as_dict(validation.get("thresholds", {}))
    _expect(c.get("nu") == cr001.get("nu") == 0.01, "CR001 viscosity drift", violations)
    _expect(
        domain.get("physical") == cr001.get("physical_domain") == "R^3",
        "CR001 physical domain drift",
        violations,
    )
    _expect(
        domain.get("evaluation_box") == cr001.get("evaluation_box"),
        "CR001 evaluation box drift",
        violations,
    )
    _expect(domain.get("support") == cr001.get("support"), "CR001 support drift", violations)
    _expect(
        domain.get("time_interval") == cr001.get("time_interval"),
        "CR001 time interval drift",
        violations,
    )
    _expect(
        forcing.get("mode") == cr001.get("forcing_mode") == "restricted_two_parameter_family",
        "CR001 forcing family drift",
        violations,
    )
    _expect(
        "No residual-dependent basis or pointwise free force"
        in str(forcing.get("restriction", "")),
        "CR001 free-force prohibition weakened",
        violations,
    )
    _expect(
        nontriviality.get("reference_energy") == cr001.get("reference_energy") == 1.0,
        "CR001 reference energy drift",
        violations,
    )
    _expect(
        nontriviality.get("reference_energy_abs_tolerance")
        == cr001.get("reference_energy_abs_tolerance")
        == 0.001,
        "CR001 energy tolerance drift",
        violations,
    )
    _expect(
        validation.get("seed") == cr001.get("validation_seed") == 914027,
        "CR001 validation seed drift",
        violations,
    )
    _expect(
        validation.get("held_out_points") == cr001.get("held_out_points") == 4096,
        "CR001 held-out count drift",
        violations,
    )
    _expect(
        validation.get("derivative_steps") == cr001.get("derivative_steps"),
        "CR001 derivative ladder drift",
        violations,
    )
    _expect(
        validation.get("quadrature_orders_per_axis")
        == cr001.get("quadrature_orders_per_axis"),
        "CR001 quadrature ladder drift",
        violations,
    )
    _expect(
        thresholds.get("pde_residual_max") == cr001.get("momentum_max") == 0.001,
        "CR001 momentum max threshold drift",
        violations,
    )
    _expect(
        thresholds.get("pde_residual_L2") == cr001.get("momentum_L2") == 0.001,
        "CR001 momentum L2 threshold drift",
        violations,
    )
    _expect(
        thresholds.get("divergence_max") == cr001.get("divergence_max") == 1e-5,
        "CR001 divergence max threshold drift",
        violations,
    )
    _expect(
        thresholds.get("divergence_L2") == cr001.get("divergence_L2") == 1e-5,
        "CR001 divergence L2 threshold drift",
        violations,
    )
    _expect(
        cr001.get("constraints_blob_sha") == CONSTRAINTS_BLOB_SHA,
        "canonical constraints blob pin drifted",
        violations,
    )
    for key in (
        "residual_defined_free_forcing_forbidden",
        "candidate_collapse_forbidden",
        "post_hoc_threshold_relaxation_forbidden",
    ):
        _expect(cr001.get(key) is True, f"CR001 prohibition weakened: {key}", violations)

    status = _as_dict(project_status)
    states = _as_dict(status.get("states", {}))
    _expect(
        status.get("candidate_family")
        == eq45.get("candidate_family")
        == "eq45_supported_velocity_candidate_v1",
        "canonical Eq45 candidate identity drift",
        violations,
    )
    for key, expected in (
        ("velocity_export_ready", True),
        ("visualization_ready", False),
        ("visual_correspondence_verified", False),
        ("pde_validated", False),
        ("paper_exact", False),
        ("openai_field_identified", False),
    ):
        _expect(eq45.get(key) is expected, f"contract Eq45 state drift: {key}", violations)
        _expect(states.get(key) is expected, f"repository Eq45 state drift: {key}", violations)

    mechanical = amplitude_dependence_witness()
    _expect(
        mechanical["same_public_shape_inputs"] is True,
        "mechanics witness shape inputs changed",
        violations,
    )
    _expect(
        mechanical["different_executable_U"] is True,
        "mechanics witness became vacuous",
        violations,
    )
    return violations


def audit_repository(root: Path) -> dict[str, Any]:
    contract = json.loads((root / CONTRACT_PATH).read_text())
    constraints = json.loads((root / CONSTRAINTS_PATH).read_text())
    project_status = json.loads((root / PROJECT_STATUS_PATH).read_text())
    source_text = (root / SOURCE_PATH).read_text()
    violations = audit_contract(contract, source_text, constraints, project_status)
    return {
        "schema": "cr002-kokuno-main-pulse-amplitude-scope-audit-v1",
        "ok": not violations,
        "violations": violations,
        "mechanics_witness": amplitude_dependence_witness(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    result = audit_repository(args.root)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
