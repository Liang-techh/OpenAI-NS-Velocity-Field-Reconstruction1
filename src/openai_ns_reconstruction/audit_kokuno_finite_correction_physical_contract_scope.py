"""CR002 audit: finite-correction protocol identity is not physical-contract identity.

Agent-3 #888 and Agent-5 #890 correctly require a complete residual, disjoint
held-in/held-out partitions, one residual-protocol SHA, preregistered restricted
forcing flags, and no residual-as-forcing shortcut.  The current typed candidate
provenance does not, however, bind the *before* and *after* candidate to one
explicit physical-contract identity (nu/domain/time/support/normalization,
pressure identity, and restricted-forcing parameter identity).

This module changes no candidate or scientific gate.  It only fail-closes the
claim that an observed finite-stage contraction was obtained under one fixed
physical contract until that identity is explicitly bound, or until a future
stage is honestly typed as a joint physical-contract update.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_finite_correction_stage_ingest_contract import (
    AGENT3_HEAD,
    AGENT3_SOURCE_BLOB,
    AGENT3_STAGE_PROTOCOL,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
)

SCHEMA = "cr002-kokuno-finite-correction-physical-contract-scope-v1"
TASK = "CR002-KOKUNO-FINITE-CORRECTION-PHYSICAL-CONTRACT-082"
CONFIG_PATH = Path("configs/kokuno_finite_correction_physical_contract_scope.json")
CONSTRAINTS_PATH = Path("configs/constraints.json")
EXPECTED_AGENT5_HEAD = "e833441279ed1f1df28c83b91b6226257e478e3f"
EXPECTED_AGENT3_HEAD = "3fd33d3ca20d3ebed7dad872d3264036bdd9cb51"
EXPECTED_AGENT3_SOURCE_BLOB = "4fc931fe7703b8b8e05efb41af957b54c9f7a4f9"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def mechanics_protocol_identity_witness() -> dict[str, Any]:
    """Show that one evaluator protocol can compare two different contracts.

    This is deliberately a one-dimensional repository mechanics example, not a
    Kokuno/OpenAI datum and not candidate residual evidence.  The force remains
    inside one fixed scalar family and is not chosen as ``f=R``; only its
    parameter changes.  Thus the residual can contract while the evaluator
    protocol label remains identical even though the physical contract changed.
    """

    transport = 1.0
    before_parameter = 0.2
    after_parameter = 0.4
    before = abs(transport - before_parameter)
    after = abs(transport - after_parameter)
    return {
        "after_forcing_parameter": after_parameter,
        "after_residual_abs": after,
        "apparent_contraction_factor": after / before,
        "before_forcing_parameter": before_parameter,
        "before_residual_abs": before,
        "candidate_scientific_evidence": False,
        "free_residual_defined_forcing_used": False,
        "kind": "autonomous_protocol-vs-physical-contract counterexample",
        "physical_contract_same": False,
        "residual_protocol_sha_same": True,
        "restricted_forcing_family": "fixed_one_basis_scalar_parameter",
        "transport_scalar": transport,
    }


def _canonical_snapshot(constraints: Mapping[str, Any]) -> dict[str, Any]:
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    forcing = constraints["forcing"]
    domain = constraints["domain"]
    nontriviality = constraints["nontriviality"]
    return {
        "amplitude_collapse_allowed": False,
        "derivative_steps": validation["derivative_steps"],
        "divergence_L2": thresholds["divergence_L2"],
        "divergence_max": thresholds["divergence_max"],
        "evaluation_box": domain["evaluation_box"],
        "forcing_mode": forcing["mode"],
        "forcing_parameters": forcing["parameters"],
        "held_out_points": validation["held_out_points"],
        "nu": constraints["nu"],
        "pde_residual_L2": thresholds["pde_residual_L2"],
        "pde_residual_max": thresholds["pde_residual_max"],
        "physical_domain": domain["physical"],
        "post_hoc_threshold_relaxation_allowed": False,
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "reference_energy": nontriviality["reference_energy"],
        "reference_energy_abs_tolerance": nontriviality[
            "reference_energy_abs_tolerance"
        ],
        "residual_defined_pointwise_free_force_allowed": False,
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "validation_seed": validation["seed"],
    }


EXPECTED_STAGE_SEMANTICS = {
    "candidate_hash_must_change_after_correction": True,
    "complete_ns_residual_flag_required": True,
    "residual_as_forcing_shortcut_forbidden": True,
    "restricted_forcing_included_and_preregistered_flags_required": True,
    "same_physical_contract_sha_before_after_enforced": False,
    "same_pressure_identity_before_after_enforced": False,
    "same_residual_protocol_sha_before_after_enforced": True,
    "same_restricted_forcing_parameter_identity_before_after_enforced": False,
    "same_viscosity_domain_support_time_normalization_identity_before_after_enforced": False,
}

EXPECTED_PROMOTION_GUARDS = {
    "joint_physical_contract_update_may_be_labeled_velocity_only_correction": False,
    "pde_pending_blocks_callable_velocity_delivery": False,
    "residual_contraction_claimable_under_fixed_physical_contract_without_identity_binding": False,
    "same_forcing_family_flags_imply_same_forcing_parameters": False,
    "same_residual_protocol_sha_implies_same_physical_contract": False,
    "silent_pressure_or_forcing_drift_allowed_in_velocity_only_correction_claim": False,
    "stage_acceptance_implies_pde_validated": False,
    "stage_accepted_implies_velocity_only_residual_improvement": False,
}

EXPECTED_REQUIRED_BINDINGS = [
    "physical_contract_sha256_before_equals_after",
    "nu_domain_support_time_normalization_identity_before_equals_after",
    "matched_pressure_identity_before_equals_after_or_explicitly_typed_joint_update",
    "restricted_forcing_family_and_parameter_identity_before_equals_after_or_explicitly_typed_joint_update",
    "same_residual_protocol_sha256",
    "disjoint_held_in_held_out_partitions",
]

EXPECTED_TRUTH_BOUNDARY = {
    "candidate_bytes_changed": False,
    "fixed_physical_contract_residual_contraction_verified": False,
    "forcing_changed": False,
    "openai_field_identified": False,
    "paper_exact": False,
    "pde_validated": False,
    "pressure_changed": False,
    "real_candidate_finite_correction_cycle_run": False,
    "threshold_changed": False,
    "velocity_changed": False,
    "visual_correspondence_verified": False,
}

EXPECTED_PROVENANCE = {
    "autonomous_design": [
        "Agent-3 typed finite-correction API and same-residual-protocol-SHA check",
        "Agent-5 registration of that finite-correction interface",
        "this CR002 physical-contract invariance audit and mechanics witness",
    ],
    "pending_unknown": [
        "before/after candidate physical-contract SHA binding",
        "before/after matched-pressure identity binding",
        "before/after restricted-forcing family and parameter identity binding",
        "before/after viscosity, domain, support, time-window and normalization identity binding",
        "whether a future accepted stage is a velocity-only correction or an explicitly typed joint physical-contract update",
    ],
    "public_source_fact": [],
    "user_requirement": [
        "CR001 physical, forcing, nontriviality and threshold contracts remain preregistered and separately governed",
        "PDE status must not block delivery of an independently governed callable velocity artifact",
    ],
}


def load_contract(root: str | Path = ".") -> dict[str, Any]:
    path = Path(root) / CONFIG_PATH
    return json.loads(path.read_text(encoding="utf-8"))


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label} drifted")


def validate_contract(payload: Mapping[str, Any], root: str | Path = ".") -> None:
    if not isinstance(payload, Mapping):
        raise ValueError("contract must be a mapping")
    _require_equal(payload.get("schema"), SCHEMA, "schema")
    _require_equal(payload.get("task"), TASK, "task")

    upstream = payload.get("audited_upstream")
    _require_equal(
        upstream,
        {
            "agent3_head": EXPECTED_AGENT3_HEAD,
            "agent3_pr": 888,
            "agent3_source_blob_sha": EXPECTED_AGENT3_SOURCE_BLOB,
            "agent5_head": EXPECTED_AGENT5_HEAD,
            "agent5_pr": 890,
            "canonical_constraints_blob_sha": EXPECTED_CONSTRAINTS_BLOB,
        },
        "audited upstream identity",
    )

    _require_equal(payload.get("provenance"), EXPECTED_PROVENANCE, "four-way provenance")
    _require_equal(
        payload.get("public_source_boundary"),
        {
            "mechanics_witness_is_public_source_data": False,
            "repository_protocol_or_threshold_is_public_source_data": False,
            "this_increment_adds_new_public_source_fact": False,
        },
        "public-source boundary",
    )
    _require_equal(
        payload.get("existing_stage_semantics"),
        EXPECTED_STAGE_SEMANTICS,
        "existing stage semantics",
    )
    _require_equal(
        payload.get("promotion_guards"), EXPECTED_PROMOTION_GUARDS, "promotion guards"
    )
    _require_equal(
        payload.get("required_future_binding_for_fixed_contract_contraction_claim"),
        EXPECTED_REQUIRED_BINDINGS,
        "future physical-contract binding",
    )
    _require_equal(
        payload.get("mechanics_witness"),
        mechanics_protocol_identity_witness(),
        "mechanics witness",
    )
    _require_equal(
        payload.get("truth_boundary"), EXPECTED_TRUTH_BOUNDARY, "truth boundary"
    )

    root_path = Path(root)
    constraints_path = root_path / CONSTRAINTS_PATH
    constraints_bytes = constraints_path.read_bytes()
    _require_equal(
        _git_blob_sha(constraints_bytes),
        EXPECTED_CONSTRAINTS_BLOB,
        "canonical constraints blob",
    )
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    _require_equal(
        payload.get("canonical_cr001_snapshot"),
        _canonical_snapshot(constraints),
        "canonical CR001 snapshot",
    )

    # Bind the exact A3/A5 registration facts consumed by this audit.  The
    # current protocol guarantees evaluator-protocol continuity but contains no
    # explicit before/after physical-contract identity field.
    _require_equal(AGENT3_HEAD, EXPECTED_AGENT3_HEAD, "registered Agent-3 head")
    _require_equal(
        AGENT3_SOURCE_BLOB, EXPECTED_AGENT3_SOURCE_BLOB, "registered Agent-3 source blob"
    )
    _require_equal(
        AGENT3_STAGE_PROTOCOL.get("same_residual_protocol_sha_required_before_and_after"),
        True,
        "same residual-protocol requirement",
    )
    _require_equal(
        AGENT3_STAGE_PROTOCOL.get("caller_supplied_pressure_or_forcing_allowed"),
        False,
        "caller pressure/forcing injection guard",
    )
    _require_equal(
        AGENT3_STAGE_PROTOCOL.get("post_observation_retuning_allowed"),
        False,
        "post-observation retuning guard",
    )
    for absent_key in (
        "same_physical_contract_sha_required_before_and_after",
        "same_pressure_identity_required_before_and_after",
        "same_restricted_forcing_parameter_identity_required_before_and_after",
        "same_viscosity_domain_support_time_normalization_identity_required_before_and_after",
    ):
        if absent_key in AGENT3_STAGE_PROTOCOL:
            raise ValueError(
                "upstream protocol now exposes a physical-contract identity binding; "
                "refresh this scoped audit instead of preserving a stale negative claim"
            )

    _require_equal(FINAL_NORMALIZED_MOMENTUM_GATE, 1.0e-3, "momentum gate")
    _require_equal(FINAL_NORMALIZED_DIVERGENCE_GATE, 1.0e-5, "divergence gate")


def audit(root: str | Path = ".") -> dict[str, Any]:
    payload = load_contract(root)
    validate_contract(payload, root)
    witness = mechanics_protocol_identity_witness()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "status": "PASS",
        "scope": "governance-only; not candidate scientific evidence",
        "agent3_head": EXPECTED_AGENT3_HEAD,
        "agent5_head": EXPECTED_AGENT5_HEAD,
        "same_residual_protocol_sha_enforced": True,
        "same_physical_contract_identity_enforced": False,
        "mechanics_same_protocol_with_contract_drift_replayed": bool(
            witness["residual_protocol_sha_same"]
            and not witness["physical_contract_same"]
            and witness["after_residual_abs"] < witness["before_residual_abs"]
        ),
        "fixed_physical_contract_residual_contraction_verified": False,
        "pde_validated": False,
        "callable_velocity_delivery_blocked": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    args = parser.parse_args(argv)
    print(json.dumps(audit(args.root), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
