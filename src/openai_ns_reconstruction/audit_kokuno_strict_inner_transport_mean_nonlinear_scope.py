"""Fail-closed CR002 audit for the Kokuno strict-inner transport mean.

The target PR computes the nonlinear Cartesian transport first and only then
projects/averages around an azimuthal ring.  This audit prevents that quantity
from being promoted to the nonlinear transport of the azimuthally averaged
velocity without a separate closure argument.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


CONTRACT_REL = Path("configs/kokuno_strict_inner_transport_mean_nonlinear_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
TARGET_REL = Path("src/openai_ns_reconstruction/kokuno_strict_inner_transport_mean.py")
EXPECTED_TASK = "CR002-KOKUNO-STRICT-INNER-TRANSPORT-MEAN-NONLINEAR-SCOPE"
EXPECTED_TARGET_HEAD = "7c3e5c23df9125722f055df4696f593ea17dd0e8"
EXPECTED_TARGET_TASK = "KOKUNO-A3-STRICT-INNER-TRANSPORT-MEAN-072"
EXPECTED_TARGET_SCHEMA = "kokuno-a3-strict-inner-transport-mean-v1"
EXPECTED_TAXONOMY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class AuditError(RuntimeError):
    """Raised when a governance invariant fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read JSON {path}: {exc}") from exc
    _require(isinstance(value, dict), f"{path} must contain a JSON object")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _equal_number(actual: Any, expected: float, label: str) -> None:
    _require(isinstance(actual, (int, float)) and not isinstance(actual, bool), f"{label} must be numeric")
    _require(float(actual) == float(expected), f"{label} drifted: {actual!r} != {expected!r}")


def _check_classification(contract: dict[str, Any]) -> None:
    taxonomy = contract.get("classification_taxonomy")
    _require(isinstance(taxonomy, list), "classification_taxonomy must be a list")
    _require(set(taxonomy) == EXPECTED_TAXONOMY and len(taxonomy) == 4, "four-way classification taxonomy drifted")
    classification = contract.get("classification")
    _require(isinstance(classification, dict), "classification must be an object")
    _require(set(classification) == EXPECTED_TAXONOMY, "classification keys must exactly match taxonomy")
    seen: set[str] = set()
    for category in taxonomy:
        entries = classification.get(category)
        _require(isinstance(entries, list) and entries, f"classification {category} must be non-empty")
        for entry in entries:
            _require(isinstance(entry, str) and entry.strip(), f"classification {category} has invalid entry")
            _require(entry not in seen, f"classification entry appears in multiple categories: {entry}")
            seen.add(entry)
    autonomous = "\n".join(classification["autonomous_design"])
    public = "\n".join(classification["public_source_fact"])
    _require("angular ladder [32,64,128]" in autonomous, "angular ladder must remain autonomous_design")
    _require("mean-closure theorem" in public, "public-source boundary must deny adding a mean-closure theorem")


def _check_target_source(contract: dict[str, Any], source_text: str) -> None:
    target = contract.get("target")
    _require(isinstance(target, dict), "target must be an object")
    _require(target.get("pull_request") == 850, "target PR drifted")
    _require(target.get("head_sha") == EXPECTED_TARGET_HEAD, "target exact head drifted")
    _require(target.get("path") == str(TARGET_REL), "target path drifted")
    _require(target.get("task") == EXPECTED_TARGET_TASK, "target task drifted")
    _require(target.get("schema") == EXPECTED_TARGET_SCHEMA, "target schema drifted")

    anchors = [
        f'TASK = "{EXPECTED_TARGET_TASK}"',
        f'SCHEMA = "{EXPECTED_TARGET_SCHEMA}"',
        "raw = provider(xq, yq, zq, tq)",
        "advection = _validate_vector(raw[1]",
        "transport = _validate_vector(raw[3]",
        "_project_cartesian_ring_to_cylindrical_mean(advection, c, s)",
        "_project_cartesian_ring_to_cylindrical_mean(transport, c, s)",
        '"complete_ns_defect": False',
        '"pde_validated": False',
        '"caller_supplied_residual_allowed": False',
    ]
    for anchor in anchors:
        _require(anchor in source_text, f"target source anchor missing: {anchor}")

    observed = contract.get("observed_operator")
    _require(isinstance(observed, dict), "observed_operator must be an object")
    _require("after evaluating raw Cartesian" in str(observed.get("reported_mean", "")), "reported mean must remain post-operator/post-vector evaluation")
    not_reported = observed.get("not_reported")
    _require(isinstance(not_reported, list), "observed_operator.not_reported must be a list")
    _require("T(<u>_theta)" in not_reported, "transport-of-mean exclusion missing")
    _require("(<u>_theta dot grad)<u>_theta" in not_reported, "advection-of-mean exclusion missing")


def _check_nonlinear_boundary(contract: dict[str, Any]) -> None:
    boundary = contract.get("nonlinear_mean_boundary")
    _require(isinstance(boundary, dict), "nonlinear_mean_boundary must be an object")
    for key in (
        "mean_of_advection_equals_advection_of_mean",
        "mean_of_transport_equals_transport_of_mean",
        "closed_mean_velocity_pde_established",
        "fluctuation_or_closure_term_proved_zero",
    ):
        _require(boundary.get(key) is False, f"forbidden nonlinear promotion enabled: {key}")
    _require("Nonlinear averaging need not commute" in str(boundary.get("rule", "")), "noncommutation rule missing")

    truth = contract.get("truth_boundary")
    _require(isinstance(truth, dict), "truth_boundary must be an object")
    _require(truth.get("strict_inner_transport_mean_bookkeeping_scoped") is True, "scoped bookkeeping truth missing")
    for key in (
        "complete_ns_defect",
        "st006_same_protocol_comparable",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pde_pending_blocks_callable_velocity_delivery",
    ):
        _require(truth.get(key) is False, f"truth-boundary promotion enabled: {key}")

    promotions = contract.get("forbidden_promotions")
    _require(isinstance(promotions, dict) and promotions, "forbidden_promotions must be non-empty")
    for key, value in promotions.items():
        _require(value is False, f"forbidden promotion enabled: {key}")


def _check_cr001(contract: dict[str, Any], constraints: dict[str, Any]) -> None:
    snap = contract.get("cr001_snapshot")
    _require(isinstance(snap, dict), "cr001_snapshot must be an object")

    _equal_number(constraints.get("nu"), 0.01, "constraints.nu")
    _equal_number(snap.get("nu"), 0.01, "snapshot.nu")
    domain = constraints.get("domain", {})
    _require(domain.get("physical") == snap.get("physical_domain") == "R^3", "physical domain drifted")
    _require(domain.get("evaluation_box") == snap.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drifted")
    _require(domain.get("support") == snap.get("support") == "r < 2 and abs(z) < 2", "support drifted")
    _require(domain.get("time_interval") == snap.get("time_interval") == [0.25, 0.75], "time interval drifted")

    forcing = constraints.get("forcing", {})
    _require(forcing.get("mode") == snap.get("forcing_mode") == "restricted_two_parameter_family", "forcing mode drifted")
    params = forcing.get("parameters", {})
    _require(params.get("a") == snap.get("forcing_parameter_bounds", {}).get("a") == [0.0, 10.0], "forcing a bounds drifted")
    _require(params.get("c") == snap.get("forcing_parameter_bounds", {}).get("c") == [0.0, 10.0], "forcing c bounds drifted")
    restriction = str(forcing.get("restriction", ""))
    _require("No residual-dependent basis or pointwise free force" in restriction, "free-force prohibition drifted")
    _require(snap.get("residual_defined_pointwise_free_force_allowed") is False, "snapshot enabled free residual force")

    nontriviality = constraints.get("nontriviality", {})
    _equal_number(nontriviality.get("reference_energy"), 1.0, "reference energy")
    _equal_number(nontriviality.get("reference_energy_abs_tolerance"), 0.001, "energy tolerance")
    _equal_number(snap.get("reference_energy"), 1.0, "snapshot reference energy")
    _equal_number(snap.get("reference_energy_abs_tolerance"), 0.001, "snapshot energy tolerance")
    _require("reject collapsed candidates" in str(nontriviality.get("enforcement", "")), "amplitude-collapse rejection drifted")
    _require(snap.get("amplitude_collapse_allowed") is False, "snapshot enabled amplitude collapse")

    validation = constraints.get("validation", {})
    _require(validation.get("seed") == snap.get("validation_seed") == 914027, "validation seed drifted")
    _require(validation.get("held_out_points") == snap.get("held_out_points") == 4096, "held-out count drifted")
    _require(validation.get("times") == snap.get("validation_times"), "validation times drifted")
    _require(validation.get("derivative_steps") == snap.get("derivative_steps") == [0.02, 0.01, 0.005], "derivative ladder drifted")
    _require(validation.get("quadrature_orders_per_axis") == snap.get("quadrature_orders_per_axis") == [24, 48, 96], "quadrature ladder drifted")
    thresholds = validation.get("thresholds", {})
    for key, expected in (
        ("divergence_max", 1e-5),
        ("divergence_L2", 1e-5),
        ("pde_residual_max", 1e-3),
        ("pde_residual_L2", 1e-3),
    ):
        _equal_number(thresholds.get(key), expected, f"constraints threshold {key}")
        _equal_number(snap.get(key), expected, f"snapshot threshold {key}")
    _require("changing thresholds requires a new experiment version" in str(validation.get("failure_policy", "")), "threshold-change policy drifted")
    _require(snap.get("post_hoc_threshold_relaxation_allowed") is False, "snapshot permits threshold relaxation")


def _manufactured_witness(contract: dict[str, Any]) -> dict[str, Any]:
    witness = contract.get("manufactured_annular_witness")
    _require(isinstance(witness, dict), "manufactured_annular_witness must be an object")
    _require(witness.get("classification") == "autonomous_design", "manufactured witness provenance drifted")
    _require(witness.get("axis_or_global_smoothness_claimed") is False, "annular witness promoted to axis/global smoothness")
    _require(witness.get("candidate_or_source_fact_claimed") is False, "annular witness promoted to candidate/source fact")
    params = witness.get("parameters", {})
    a = float(params.get("a"))
    radius = float(params.get("r"))
    _require(radius > 0.0, "manufactured witness must remain off-axis")

    n = 4096
    sum_ur = 0.0
    sum_ut = 0.0
    sum_ar = 0.0
    sum_at = 0.0
    for k in range(n):
        theta = 2.0 * math.pi * k / n
        s = math.sin(theta)
        c = math.cos(theta)
        ur = a * c
        ut = a * s
        # No r dependence. Cylindrical convective terms retain basis geometry.
        ar = (ut / radius) * (-a * s) - (ut * ut) / radius
        at = (ut / radius) * (a * c) + (ur * ut) / radius
        sum_ur += ur
        sum_ut += ut
        sum_ar += ar
        sum_at += at
    mean_velocity = [sum_ur / n, sum_ut / n, 0.0]
    mean_advection = [sum_ar / n, sum_at / n, 0.0]
    analytic = [-a * a / radius, 0.0, 0.0]
    declared_velocity = [float(x) for x in witness.get("mean_velocity_cylindrical", [])]
    declared_advection = [float(x) for x in witness.get("mean_advection_cylindrical", [])]
    _require(len(declared_velocity) == 3 and len(declared_advection) == 3, "manufactured witness vectors malformed")
    for actual, expected in zip(mean_velocity, [0.0, 0.0, 0.0]):
        _require(abs(actual - expected) < 1e-12, "sampled manufactured mean velocity is not zero")
    for actual, expected in zip(mean_advection, analytic):
        _require(abs(actual - expected) < 1e-12, "sampled mean advection disagrees with analytic counterexample")
    for actual, expected in zip(declared_velocity, [0.0, 0.0, 0.0]):
        _require(abs(actual - expected) < 1e-15, "declared manufactured mean velocity drifted")
    for actual, expected in zip(declared_advection, analytic):
        _require(abs(actual - expected) < 1e-15, "declared manufactured mean advection drifted")
    _require(abs(analytic[0]) > 0.0, "manufactured witness lost nonlinear noncommutation")
    return {
        "samples": n,
        "a": a,
        "r": radius,
        "sampled_mean_velocity_cylindrical": mean_velocity,
        "sampled_mean_advection_cylindrical": mean_advection,
        "analytic_mean_advection_cylindrical": analytic,
        "advection_of_mean_velocity_cylindrical": [0.0, 0.0, 0.0],
        "noncommutation_radial_gap": analytic[0],
    }


def audit_repo(repo_root: Path) -> dict[str, Any]:
    repo_root = Path(repo_root)
    contract_path = repo_root / CONTRACT_REL
    constraints_path = repo_root / CONSTRAINTS_REL
    target_path = repo_root / TARGET_REL
    contract = _load_json(contract_path)
    constraints = _load_json(constraints_path)
    try:
        source_text = target_path.read_text()
    except OSError as exc:
        raise AuditError(f"cannot read target source {target_path}: {exc}") from exc

    _require(contract.get("schema_version") == 1, "contract schema_version drifted")
    _require(contract.get("task") == EXPECTED_TASK, "contract task drifted")
    _require(contract.get("status") == "governance_only_not_scientific_admission", "contract status drifted")
    _check_classification(contract)
    _check_target_source(contract, source_text)
    _check_nonlinear_boundary(contract)
    _check_cr001(contract, constraints)
    witness = _manufactured_witness(contract)
    return {
        "ok": True,
        "task": EXPECTED_TASK,
        "contract_sha256": _sha256(contract_path),
        "target_source_sha256": _sha256(target_path),
        "constraints_sha256": _sha256(constraints_path),
        "target_head": EXPECTED_TARGET_HEAD,
        "manufactured_witness": witness,
        "claim": "mean of nonlinear transport is governed separately from transport of mean velocity",
        "scientific_status_promoted": False,
    }


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = audit_repo(args.repo_root)
    except AuditError as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, sort_keys=True))
        return 1
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
