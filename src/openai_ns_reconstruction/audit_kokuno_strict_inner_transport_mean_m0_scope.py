"""Fail-closed CR002 audit for the Kokuno strict-inner transport m=0 scope.

PR #850 evaluates a raw strict-inner transport on azimuthal rings, projects it
into the rotating cylindrical frame, and reports both an m=0 mean RMS and a
raw full-ring RMS.  This audit governs the information lost by the m=0
projection: a small or zero mean does not control non-axisymmetric angular
modes and neither finite-ring RMS is the canonical CR001 whole-volume momentum
L2 norm.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any


CONTRACT_REL = Path("configs/kokuno_strict_inner_transport_mean_m0_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
TARGET_REL = Path("src/openai_ns_reconstruction/kokuno_strict_inner_transport_mean.py")

EXPECTED_TASK = "CR002-KOKUNO-STRICT-INNER-TRANSPORT-M0-NULLSPACE-SCOPE"
EXPECTED_TARGET_HEAD = "7c3e5c23df9125722f055df4696f593ea17dd0e8"
EXPECTED_TARGET_BLOB = "244006a59e4d8a68e9acd7c7cc1c885be279f7f7"
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
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read JSON {path}: {exc}") from exc
    _require(isinstance(value, dict), f"{path} must contain a JSON object")
    return value


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _equal_number(actual: Any, expected: float, label: str) -> None:
    _require(
        isinstance(actual, (int, float)) and not isinstance(actual, bool),
        f"{label} must be numeric",
    )
    _require(
        float(actual) == float(expected),
        f"{label} drifted: {actual!r} != {expected!r}",
    )


def _check_classification(contract: dict[str, Any]) -> None:
    taxonomy = contract.get("classification_taxonomy")
    _require(isinstance(taxonomy, list), "classification_taxonomy must be a list")
    _require(
        set(taxonomy) == EXPECTED_TAXONOMY and len(taxonomy) == 4,
        "four-way classification taxonomy drifted",
    )
    classification = contract.get("classification")
    _require(isinstance(classification, dict), "classification must be an object")
    _require(
        set(classification) == EXPECTED_TAXONOMY,
        "classification keys must exactly match taxonomy",
    )
    seen: set[str] = set()
    for category in taxonomy:
        entries = classification.get(category)
        _require(
            isinstance(entries, list) and entries,
            f"classification {category} must be non-empty",
        )
        for entry in entries:
            _require(
                isinstance(entry, str) and entry.strip(),
                f"classification {category} has invalid entry",
            )
            _require(
                entry not in seen,
                f"classification entry appears in multiple categories: {entry}",
            )
            seen.add(entry)

    autonomous = "\n".join(classification["autonomous_design"]).lower()
    public = "\n".join(classification["public_source_fact"]).lower()
    _require(
        "angular ladder [32,64,128]" in autonomous,
        "angular ladder must remain autonomous_design",
    )
    _require(
        "manufactured pure-m=2" in autonomous,
        "manufactured nullspace witness must remain autonomous_design",
    )
    _require(
        "no public source theorem" in public,
        "public-source bucket must preserve the no-control-theorem boundary",
    )
    _require(
        "angular ladder" not in public and "manufactured" not in public,
        "repository validation mechanics cannot be laundered into public_source_fact",
    )


def _check_target(contract: dict[str, Any], source_path: Path, source_text: str) -> str:
    target = contract.get("target")
    _require(isinstance(target, dict), "target must be an object")
    expected = {
        "pull_request": 850,
        "head_sha": EXPECTED_TARGET_HEAD,
        "path": str(TARGET_REL),
        "source_blob_sha": EXPECTED_TARGET_BLOB,
        "task": EXPECTED_TARGET_TASK,
        "schema": EXPECTED_TARGET_SCHEMA,
    }
    for key, value in expected.items():
        _require(target.get(key) == value, f"target {key} drifted")

    try:
        actual_blob = _git_blob_sha(source_path.read_bytes())
    except OSError as exc:
        raise AuditError(f"cannot read target source {source_path}: {exc}") from exc
    _require(actual_blob == EXPECTED_TARGET_BLOB, "target source blob drifted")

    anchors = [
        f'TASK = "{EXPECTED_TARGET_TASK}"',
        f'SCHEMA = "{EXPECTED_TARGET_SCHEMA}"',
        "mean_transport = _project_cartesian_ring_to_cylindrical_mean(transport, c, s)",
        "full_ring_rms = float(",
        "np.sqrt(np.mean(np.sum(transport * transport, axis=-1)))",
        "mean_rms = float(",
        "ratio = mean_rms / max(float(full_rms), np.finfo(float).tiny)",
        '"complete_ns_defect": False',
        '"pde_validated": False',
        '"caller_supplied_residual_allowed": False',
    ]
    for anchor in anchors:
        _require(anchor in source_text, f"target source anchor missing: {anchor}")
    return actual_blob


def _check_projection_boundary(contract: dict[str, Any]) -> None:
    metrics = contract.get("reported_metrics")
    _require(isinstance(metrics, dict), "reported_metrics must be an object")
    _require(
        metrics.get("metrics_are_interchangeable") is False,
        "reported metrics must remain non-interchangeable",
    )
    _require(
        "azimuthal m=0" in str(metrics.get("transport_mean_rms")),
        "transport_mean_rms scope drifted",
    )
    _require(
        "unweighted RMS" in str(metrics.get("full_ring_transport_rms")),
        "full_ring_transport_rms scope drifted",
    )
    _require(
        "volume-weighted spatial L2" in str(metrics.get("canonical_cr001_momentum_L2")),
        "canonical momentum-L2 description drifted",
    )

    boundary = contract.get("projection_boundary")
    _require(isinstance(boundary, dict), "projection_boundary must be an object")
    _require(boundary.get("m0_nullspace_nontrivial") is True, "m0 nullspace must remain nontrivial")
    _require(
        boundary.get("reverse_bound_from_mean_to_full_supported") is False,
        "reverse mean-to-full bound cannot be promoted",
    )
    _require(
        boundary.get("zero_m0_mean_implies_zero_full_ring_transport") is False,
        "zero m0 mean cannot imply zero full-ring transport",
    )
    _require(
        boundary.get("small_m0_mean_implies_small_full_ring_transport") is False,
        "small m0 mean cannot imply small full-ring transport",
    )
    _require(
        boundary.get("m0_mean_controls_nonaxisymmetric_modes") is False,
        "m0 mean cannot control non-axisymmetric modes",
    )
    rule = str(boundary.get("rule", ""))
    _require("m!=0" in rule and "cancel" in rule, "projection-nullspace rule drifted")


def _manufactured_witness(contract: dict[str, Any]) -> dict[str, float]:
    witness = contract.get("manufactured_nullspace_witness")
    _require(isinstance(witness, dict), "manufactured_nullspace_witness must be an object")
    _require(witness.get("classification") == "autonomous_design", "witness classification drifted")
    _require(witness.get("field") == "T(theta)=amplitude*cos(2*theta)*e_z", "witness field drifted")
    _require(witness.get("axis_or_global_smoothness_claimed") is False, "witness cannot claim global smoothness")
    _require(witness.get("candidate_or_source_fact_claimed") is False, "witness cannot become candidate/source data")

    order = witness.get("angular_order")
    _require(isinstance(order, int) and order == 128, "manufactured angular order drifted")
    amplitude = witness.get("amplitude")
    _equal_number(amplitude, 0.6, "manufactured amplitude")

    values = [
        float(amplitude) * math.cos(2.0 * (2.0 * math.pi * k / order))
        for k in range(order)
    ]
    mean_z = sum(values) / order
    full_rms = math.sqrt(sum(value * value for value in values) / order)
    expected_rms = float(amplitude) / math.sqrt(2.0)
    _require(abs(mean_z) <= 2e-15, f"manufactured m0 mean is not numerically zero: {mean_z}")
    _require(abs(full_rms - expected_rms) <= 2e-15, "manufactured full-ring RMS drifted")
    expected_vector = witness.get("expected_m0_mean_vector")
    _require(expected_vector == [0.0, 0.0, 0.0], "expected m0 mean vector drifted")
    _require(abs(float(witness.get("expected_full_ring_rms")) - expected_rms) <= 1e-15, "expected full-ring RMS drifted")
    _equal_number(witness.get("expected_mean_to_full_ratio"), 0.0, "expected mean/full ratio")
    _require(full_rms > 0.4, "manufactured full-ring transport must remain nonzero")
    return {
        "computed_m0_mean_z": mean_z,
        "computed_full_ring_rms": full_rms,
        "computed_mean_to_full_ratio": abs(mean_z) / full_rms,
    }


def _check_norm_scope(contract: dict[str, Any]) -> None:
    scope = contract.get("norm_scope")
    _require(isinstance(scope, dict), "norm_scope must be an object")
    false_keys = [
        "pr850_full_ring_rms_is_whole_domain_volume_weighted_L2",
        "pr850_mean_rms_is_whole_domain_volume_weighted_L2",
        "pr850_transport_contains_pressure_gradient",
        "pr850_transport_contains_restricted_forcing",
        "pr850_transport_is_complete_ns_momentum_residual",
        "st006_same_protocol_comparable",
    ]
    for key in false_keys:
        _require(scope.get(key) is False, f"norm-scope promotion forbidden: {key}")


def _check_cr001(contract: dict[str, Any], constraints: dict[str, Any]) -> None:
    snapshot = contract.get("cr001_snapshot")
    _require(isinstance(snapshot, dict), "cr001_snapshot must be an object")

    _equal_number(snapshot.get("nu"), constraints.get("nu"), "nu")
    domain = constraints.get("domain")
    _require(isinstance(domain, dict), "canonical domain missing")
    _require(snapshot.get("physical_domain") == domain.get("physical"), "physical domain drifted")
    _require(snapshot.get("evaluation_box") == domain.get("evaluation_box"), "evaluation box drifted")
    _require(snapshot.get("support") == domain.get("support"), "support drifted")
    _require(snapshot.get("time_interval") == domain.get("time_interval"), "time interval drifted")

    forcing = constraints.get("forcing")
    _require(isinstance(forcing, dict), "canonical forcing missing")
    _require(snapshot.get("forcing_mode") == forcing.get("mode"), "forcing mode drifted")
    _require(snapshot.get("forcing_parameter_bounds") == forcing.get("parameters"), "forcing bounds drifted")
    _require(
        "No residual-dependent basis or pointwise free force" in str(forcing.get("restriction")),
        "canonical free-force prohibition drifted",
    )

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, dict), "canonical nontriviality missing")
    _equal_number(snapshot.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy")
    _equal_number(
        snapshot.get("reference_energy_abs_tolerance"),
        nontriviality.get("reference_energy_abs_tolerance"),
        "reference energy tolerance",
    )
    _require(
        "reject collapsed candidates" in str(nontriviality.get("enforcement")),
        "canonical amplitude-collapse prohibition drifted",
    )

    validation = constraints.get("validation")
    _require(isinstance(validation, dict), "canonical validation missing")
    _require(snapshot.get("validation_seed") == validation.get("seed"), "validation seed drifted")
    _require(snapshot.get("held_out_points") == validation.get("held_out_points"), "held-out count drifted")
    _require(snapshot.get("validation_times") == validation.get("times"), "validation times drifted")
    _require(snapshot.get("derivative_steps") == validation.get("derivative_steps"), "derivative ladder drifted")
    _require(
        snapshot.get("quadrature_orders_per_axis") == validation.get("quadrature_orders_per_axis"),
        "quadrature ladder drifted",
    )
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, dict), "canonical thresholds missing")
    for snap_key, canonical_key in (
        ("divergence_max", "divergence_max"),
        ("divergence_L2", "divergence_L2"),
        ("pde_residual_max", "pde_residual_max"),
        ("pde_residual_L2", "pde_residual_L2"),
    ):
        _equal_number(snapshot.get(snap_key), thresholds.get(canonical_key), snap_key)

    for key in (
        "residual_defined_pointwise_free_force_allowed",
        "amplitude_collapse_allowed",
        "post_hoc_threshold_relaxation_allowed",
    ):
        _require(snapshot.get(key) is False, f"{key} must remain false")
    _require(
        "changing thresholds requires a new experiment version" in str(validation.get("failure_policy")),
        "canonical threshold-relaxation policy drifted",
    )


def _check_truth(contract: dict[str, Any]) -> None:
    truth = contract.get("truth_boundary")
    _require(isinstance(truth, dict), "truth_boundary must be an object")
    _require(truth.get("strict_inner_transport_mean_bookkeeping_scoped") is True, "scoped bookkeeping flag drifted")
    _require(truth.get("m0_projection_nullspace_governed") is True, "m0 governance flag drifted")
    for key in (
        "complete_ns_defect",
        "whole_domain_momentum_L2_assessed",
        "st006_same_protocol_comparable",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pde_pending_blocks_callable_velocity_delivery",
    ):
        _require(truth.get(key) is False, f"truth promotion forbidden: {key}")

    forbidden = contract.get("forbidden_promotions")
    _require(isinstance(forbidden, dict) and forbidden, "forbidden_promotions must be non-empty")
    for key, value in forbidden.items():
        _require(value is False, f"forbidden promotion enabled: {key}")

    sibling = contract.get("relationship_to_sibling_governance")
    _require(isinstance(sibling, dict), "relationship_to_sibling_governance missing")
    _require(sibling.get("duplicate") is False, "this increment must remain non-duplicate")
    _require("nonlinear" in str(sibling.get("pr853", "")).lower(), "PR #853 boundary drifted")
    _require("nullspace" in str(sibling.get("this_increment", "")).lower(), "this increment boundary drifted")


def audit_repository(
    root: Path,
    *,
    contract_override: dict[str, Any] | None = None,
    constraints_override: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract = copy.deepcopy(contract_override) if contract_override is not None else _load_json(root / CONTRACT_REL)
    constraints = copy.deepcopy(constraints_override) if constraints_override is not None else _load_json(root / CONSTRAINTS_REL)

    _require(contract.get("task") == EXPECTED_TASK, "contract task drifted")
    _require(contract.get("status") == "governance_only_not_scientific_admission", "contract status drifted")
    _check_classification(contract)

    source_path = root / TARGET_REL
    try:
        source_text = source_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise AuditError(f"cannot read target source {source_path}: {exc}") from exc
    actual_blob = _check_target(contract, source_path, source_text)
    _check_projection_boundary(contract)
    witness = _manufactured_witness(contract)
    _check_norm_scope(contract)
    _check_cr001(contract, constraints)
    _check_truth(contract)

    return {
        "task": EXPECTED_TASK,
        "status": "PASS",
        "target_pr": 850,
        "target_head": EXPECTED_TARGET_HEAD,
        "target_source_blob_sha": actual_blob,
        "manufactured_nullspace_witness": witness,
        "scope": "governance_only_not_scientific_admission",
        "pde_validated": False,
        "callable_velocity_delivery_blocked": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = audit_repository(args.root)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
