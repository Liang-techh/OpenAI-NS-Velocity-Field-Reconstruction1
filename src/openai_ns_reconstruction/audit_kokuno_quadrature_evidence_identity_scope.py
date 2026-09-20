"""CR002 audit: a quadrature boolean is not canonical quadrature evidence.

Agent-4 PR #896 deliberately keeps sampled/Monte-Carlo quantities separate
from the preregistered CR001 24/48/96 volume quadrature.  Its current
``validate_complete_candidate`` API nevertheless accepts caller-provided
``volume_weights`` and a bare ``quadrature_ladder_assessed`` boolean.  The
boolean can participate directly in the final ``pde_validated`` expression,
but it carries no candidate, physical-contract, operator, quadrature-level,
metric, or receipt identity.

This module changes no candidate, threshold, pressure, forcing or scientific
state.  It records that representation/provenance boundary and fail-closes any
future claim that a bare flag or held-out weighting itself is the canonical
CR001 volume-L2 assessment.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-quadrature-evidence-identity-scope-v1"
TASK = "CR002-KOKUNO-QUADRATURE-EVIDENCE-IDENTITY-083"
CONFIG_PATH = Path("configs/kokuno_quadrature_evidence_identity_scope.json")
CONSTRAINTS_PATH = Path("configs/constraints.json")
AGENT4_SOURCE_PATH = Path(
    "src/openai_ns_reconstruction/kokuno_a4_blackbox_full_ns_validator.py"
)

EXPECTED_AGENT4_PR = 896
EXPECTED_AGENT4_HEAD = "e8a7712515f75bb8a0a5d86b9a6177ab9447f2b1"
EXPECTED_AGENT4_SOURCE_BLOB = "a9e46daf60d909a63261dd10821e1a08101dd872"
EXPECTED_AGENT5_PR = 897
EXPECTED_AGENT5_HEAD = "21ff4a86638a1d8e179e25ed25886a68a34b4e6b"
EXPECTED_AGENT5_SOURCE_BLOB = "fcf2f8945592782b2dd4008007a970d5d4fdf83d"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

EXPECTED_OBSERVED_A4 = {
    "validate_complete_candidate_accepts_caller_volume_weights": True,
    "validate_complete_candidate_accepts_bare_quadrature_ladder_assessed_boolean": True,
    "sample_volume_l2_is_computed_from_caller_weights": True,
    "pde_validated_formula_depends_on_bare_quadrature_boolean": True,
    "typed_canonical_quadrature_receipt_parameter_present": False,
    "quadrature_receipt_sha256_present": False,
    "quadrature_receipt_candidate_sha256_binding_present": False,
    "quadrature_receipt_physical_contract_sha256_binding_present": False,
    "quadrature_receipt_operator_identity_binding_present": False,
    "quadrature_receipt_24_48_96_metrics_binding_present": False,
}

EXPECTED_PROMOTION_GUARDS = {
    "bare_boolean_is_canonical_quadrature_evidence": False,
    "caller_supplied_sample_weights_are_canonical_24_48_96_quadrature": False,
    "sampled_or_monte_carlo_l2_may_replace_canonical_volume_l2": False,
    "quadrature_ladder_assessed_true_alone_may_authorize_pde_validated": False,
    "quadrature_evidence_may_be_reused_across_candidate_sha": False,
    "quadrature_evidence_may_be_reused_across_physical_contract_sha": False,
    "quadrature_evidence_may_omit_operator_and_derivative_identity": False,
    "quadrature_evidence_may_omit_per_level_metrics": False,
    "pde_pending_blocks_callable_velocity_delivery": False,
    "validator_contract_registration_implies_pde_validated": False,
}

EXPECTED_TRUTH_BOUNDARY = {
    "candidate_bytes_changed": False,
    "velocity_changed": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "threshold_changed": False,
    "canonical_quadrature_receipt_exists": False,
    "real_kokuno_complete_ns_residual_assessed": False,
    "velocity_export_ready_promoted": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}

REQUIRED_FUTURE_RECEIPT_KEYS = (
    "receipt_schema",
    "receipt_sha256",
    "candidate_sha256",
    "physical_contract_sha256",
    "validator_operator_identity",
    "derivative_step_ladder",
    "validation_times",
    "quadrature_orders_per_axis",
    "momentum_volume_l2_by_time_and_order",
    "divergence_volume_l2_by_time_and_order",
    "energy_by_time_and_order",
    "convergence_assessment",
    "thresholds_bound",
    "independent_of_training_and_held_in_selection",
    "residual_defined_forcing_forbidden",
)


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label} drifted: {actual!r} != {expected!r}")


def _canonical_snapshot(constraints: Mapping[str, Any]) -> dict[str, Any]:
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    return {
        "nu": constraints["nu"],
        "physical_domain": domain["physical"],
        "evaluation_box": domain["evaluation_box"],
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "validation_times": validation["times"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "norms": validation["norms"],
        "divergence_max": thresholds["divergence_max"],
        "divergence_L2": thresholds["divergence_L2"],
        "pde_residual_max": thresholds["pde_residual_max"],
        "pde_residual_L2": thresholds["pde_residual_L2"],
        "reference_energy": nontriviality["reference_energy"],
        "reference_energy_abs_tolerance": nontriviality[
            "reference_energy_abs_tolerance"
        ],
        "forcing_mode": forcing["mode"],
        "residual_defined_pointwise_free_force_allowed": False,
        "amplitude_collapse_allowed": False,
        "post_hoc_threshold_relaxation_allowed": False,
    }


def _function(tree: ast.AST, name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise ValueError(f"required function missing: {name}")


def _class(tree: ast.AST, name: str) -> ast.ClassDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == name:
            return node
    raise ValueError(f"required class missing: {name}")


def observed_agent4_semantics(source: str) -> dict[str, bool]:
    """Extract only the representation facts needed by this CR002 audit."""

    tree = ast.parse(source)
    validate_fn = _function(tree, "validate_complete_candidate")
    arg_names = [arg.arg for arg in validate_fn.args.args + validate_fn.args.kwonlyargs]
    receipt_cls = _class(tree, "FullNSValidationReceipt")
    receipt_fields = {
        node.target.id
        for node in receipt_cls.body
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
    }
    compact = "".join(source.split())
    pde_formula = (
        "pde_validated=bool(heldout_assessedandquadrature_ladder_assessedandnotfailures)"
        in compact
    )
    sample_l2_from_weights = (
        "momentum_l2=tuple(_volume_l2(level.residual/RESIDUAL_SCALE,weights)"
        in compact
        and "divergence_l2=tuple(_volume_l2(level.divergence,weights)" in compact
    )
    return {
        "validate_complete_candidate_accepts_caller_volume_weights": (
            "volume_weights" in arg_names
        ),
        "validate_complete_candidate_accepts_bare_quadrature_ladder_assessed_boolean": (
            "quadrature_ladder_assessed" in arg_names
        ),
        "sample_volume_l2_is_computed_from_caller_weights": sample_l2_from_weights,
        "pde_validated_formula_depends_on_bare_quadrature_boolean": pde_formula,
        "typed_canonical_quadrature_receipt_parameter_present": (
            "quadrature_receipt" in arg_names
            or "canonical_quadrature_receipt" in arg_names
        ),
        "quadrature_receipt_sha256_present": "quadrature_receipt_sha256" in receipt_fields,
        "quadrature_receipt_candidate_sha256_binding_present": (
            "quadrature_candidate_sha256" in receipt_fields
        ),
        "quadrature_receipt_physical_contract_sha256_binding_present": (
            "quadrature_physical_contract_sha256" in receipt_fields
        ),
        "quadrature_receipt_operator_identity_binding_present": (
            "quadrature_operator_identity" in receipt_fields
        ),
        "quadrature_receipt_24_48_96_metrics_binding_present": (
            "quadrature_metrics_by_order" in receipt_fields
        ),
    }


def quadrature_receipt_identity_failures(
    receipt: Any,
    *,
    expected_candidate_sha256: str,
    expected_physical_contract_sha256: str,
) -> tuple[str, ...]:
    """Check provenance identity only; this does not certify scientific metrics.

    A bare bool is intentionally invalid.  A future receipt remains responsible
    for carrying independently computed scientific values; this helper merely
    prevents identity-free evidence laundering across candidates/contracts.
    """

    if not isinstance(receipt, Mapping):
        return ("typed_receipt_required",)
    value = dict(receipt)
    failures: list[str] = []
    for key in REQUIRED_FUTURE_RECEIPT_KEYS:
        if key not in value:
            failures.append(f"missing:{key}")
    if failures:
        return tuple(failures)
    if value["candidate_sha256"] != expected_candidate_sha256:
        failures.append("candidate_sha256_mismatch")
    if value["physical_contract_sha256"] != expected_physical_contract_sha256:
        failures.append("physical_contract_sha256_mismatch")
    if value["derivative_step_ladder"] != [0.02, 0.01, 0.005]:
        failures.append("derivative_step_ladder_drift")
    if value["validation_times"] != [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]:
        failures.append("validation_times_drift")
    if value["quadrature_orders_per_axis"] != [24, 48, 96]:
        failures.append("quadrature_ladder_drift")
    if value["thresholds_bound"] != {
        "momentum_volume_l2": 0.001,
        "divergence_volume_l2": 1e-05,
        "reference_energy_abs_tolerance": 0.001,
    }:
        failures.append("threshold_binding_drift")
    if value["independent_of_training_and_held_in_selection"] is not True:
        failures.append("independence_not_bound")
    if value["residual_defined_forcing_forbidden"] is not True:
        failures.append("residual_defined_forcing_not_forbidden")
    sha = value.get("receipt_sha256")
    if not isinstance(sha, str) or len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha):
        failures.append("receipt_sha256_invalid")
    return tuple(failures)


def audit(root: str | Path = ".") -> dict[str, Any]:
    root_path = Path(root)
    config_path = root_path / CONFIG_PATH
    constraints_path = root_path / CONSTRAINTS_PATH
    source_path = root_path / AGENT4_SOURCE_PATH

    config_bytes = config_path.read_bytes()
    constraints_bytes = constraints_path.read_bytes()
    source_bytes = source_path.read_bytes()
    config = json.loads(config_bytes)
    constraints = json.loads(constraints_bytes)
    source = source_bytes.decode("utf-8")

    _require_equal(config.get("schema"), SCHEMA, "schema")
    _require_equal(config.get("task"), TASK, "task")

    upstream = config["audited_upstream"]
    _require_equal(upstream["agent4_pr"], EXPECTED_AGENT4_PR, "Agent-4 PR")
    _require_equal(upstream["agent4_head"], EXPECTED_AGENT4_HEAD, "Agent-4 head")
    _require_equal(
        upstream["agent4_source_blob_sha"], EXPECTED_AGENT4_SOURCE_BLOB, "Agent-4 source blob"
    )
    _require_equal(upstream["agent5_registration_pr"], EXPECTED_AGENT5_PR, "Agent-5 PR")
    _require_equal(upstream["agent5_registration_head"], EXPECTED_AGENT5_HEAD, "Agent-5 head")
    _require_equal(
        upstream["agent5_registration_source_blob_sha"],
        EXPECTED_AGENT5_SOURCE_BLOB,
        "Agent-5 source blob",
    )
    _require_equal(
        upstream["canonical_constraints_blob_sha"],
        EXPECTED_CONSTRAINTS_BLOB,
        "constraints blob registration",
    )
    _require_equal(
        _git_blob_sha(source_bytes), EXPECTED_AGENT4_SOURCE_BLOB, "checked Agent-4 source blob"
    )
    _require_equal(
        _git_blob_sha(constraints_bytes), EXPECTED_CONSTRAINTS_BLOB, "checked constraints blob"
    )

    _require_equal(
        config["canonical_cr001_snapshot"], _canonical_snapshot(constraints), "CR001 snapshot"
    )
    observed = observed_agent4_semantics(source)
    _require_equal(observed, EXPECTED_OBSERVED_A4, "Agent-4 validator semantics")
    _require_equal(
        config["observed_agent4_validator_semantics"],
        EXPECTED_OBSERVED_A4,
        "registered Agent-4 semantics",
    )
    _require_equal(
        config["promotion_guards"], EXPECTED_PROMOTION_GUARDS, "promotion guards"
    )
    _require_equal(
        config["truth_boundary"], EXPECTED_TRUTH_BOUNDARY, "truth boundary"
    )

    a5 = config["observed_agent5_registration_semantics"]
    _require_equal(a5["canonical_quadrature_ladder_registered"], [24, 48, 96], "A5 quadrature ladder")
    _require(a5["quadrature_ladder_assessed_required_for_pde_validated"] is True, "A5 quadrature requirement weakened")
    _require(a5["monte_carlo_volume_l2_does_not_replace_canonical_quadrature_ladder"] is True, "A5 Monte-Carlo firewall weakened")
    _require(a5["scientific_admission"] is False, "registration cannot be scientific admission")

    required = config["required_future_quadrature_receipt"]
    _require_equal(required["derivative_step_ladder"], [0.02, 0.01, 0.005], "future receipt derivative ladder")
    _require_equal(required["validation_times"], [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75], "future receipt times")
    _require_equal(required["quadrature_orders_per_axis"], [24, 48, 96], "future receipt quadrature ladder")
    _require_equal(required["thresholds_bound"], {
        "momentum_volume_l2": 0.001,
        "divergence_volume_l2": 1e-05,
        "reference_energy_abs_tolerance": 0.001,
    }, "future receipt thresholds")
    _require(required["independent_of_training_and_held_in_selection"] is True, "future quadrature must stay independent")
    _require(required["residual_defined_forcing_forbidden"] is True, "free residual forcing cannot enter quadrature evidence")

    witness = config["mechanics_witness"]
    _require(witness["candidate_scientific_evidence"] is False, "mechanics witness cannot become scientific evidence")
    _require(witness["canonical_quadrature_was_run"] is False, "mechanics witness cannot claim quadrature was run")
    _require(witness["caller_sets_quadrature_ladder_assessed"] is True, "mechanics witness must preserve the boolean seam")

    provenance = config["provenance"]
    _require(set(provenance) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}, "four provenance classes required")
    _require(provenance["public_source_fact"] == [], "this increment adds no public-source fact")
    _require(len(provenance["pending_unknown"]) >= 5, "pending scientific evidence was dropped")

    public_boundary = config["public_source_boundary"]
    _require(public_boundary["this_increment_adds_new_public_source_fact"] is False, "audit cannot invent public-source facts")
    _require(public_boundary["repository_validator_mechanics_are_public_source_fact"] is False, "repository mechanics cannot be laundered into public-source facts")
    _require(public_boundary["mechanics_witness_is_public_source_data"] is False, "mechanics witness cannot be source data")
    _require(public_boundary["canonical_cr001_thresholds_are_openai_or_kokuno_paper_parameters"] is False, "repository thresholds are not paper parameters")

    # A bare True is exactly the representation seam under audit: it is not a
    # provenance-bearing receipt and must fail the future-receipt identity guard.
    bool_failures = quadrature_receipt_identity_failures(
        True,
        expected_candidate_sha256="0" * 64,
        expected_physical_contract_sha256="1" * 64,
    )
    _require_equal(bool_failures, ("typed_receipt_required",), "bare-boolean guard")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "config_blob_sha": _git_blob_sha(config_bytes),
        "agent4_source_blob_sha": _git_blob_sha(source_bytes),
        "constraints_blob_sha": _git_blob_sha(constraints_bytes),
        "bare_boolean_receipt_failures": list(bool_failures),
        "canonical_quadrature_receipt_exists": False,
        "pde_validated": False,
        "velocity_or_scientific_state_changed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--output")
    args = parser.parse_args()
    result = audit(args.root)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
