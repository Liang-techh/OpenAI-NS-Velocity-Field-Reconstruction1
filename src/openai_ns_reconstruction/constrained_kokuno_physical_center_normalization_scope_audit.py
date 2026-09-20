"""CR002 audit separating Kokuno PA.10 C-normalization from CR001 energy normalization.

The PA.10 physical-center map uses ``F=(phi_*/C)Phi`` and carries the source
condition ``C >= sup_Omega |phi_*|``.  CR001 separately preregisters a project
candidate nontriviality condition ``E(0.25)=1+-0.001``.  These are different
normalization contracts and neither may be used to certify or replace the
other.

This audit is intentionally claim/representation governance only.  It changes
no velocity, coefficient, pressure, forcing family, residual operator, sample,
threshold, source datum, or candidate selection.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_physical_center_ingest_contract as a5


SCHEMA = "cr002-kokuno-physical-center-normalization-scope-v1"
AUDITED_PARENT_PR = 791
AUDITED_PARENT_HEAD = "566a3fe7d53918c19016aa11d55c3e071dff779d"
ACTIVE_INTEGRATION_HEAD_AT_FREEZE = "43d4c0309cdd22763eaba966de439cc1b4e1d9d9"
CONTRACT_RELATIVE_PATH = Path("configs/kokuno_physical_center_normalization_scope_contract.json")
CONSTRAINTS_RELATIVE_PATH = Path("configs/constraints.json")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _canonical_digest(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode(
            "utf-8"
        )
    ).hexdigest()


def _expect_equal(observed: Any, expected: Any, label: str) -> None:
    if observed != expected:
        raise ValueError(f"{label} drift: observed={observed!r} expected={expected!r}")


def _expect_false(mapping: Mapping[str, Any], key: str) -> None:
    if mapping.get(key) is not False:
        raise ValueError(f"{key} must remain false")


def _constraints_projection(constraints: Mapping[str, Any]) -> dict[str, Any]:
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    return {
        "nu": constraints["nu"],
        "physical_domain": domain["physical"],
        "evaluation_box": domain["evaluation_box"],
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "forcing_mode": forcing["mode"],
        "forcing_parameter_bounds": forcing["parameters"],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "validation_times": validation["times"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "divergence_max": thresholds["divergence_max"],
        "divergence_L2": thresholds["divergence_L2"],
        "momentum_max": thresholds["pde_residual_max"],
        "momentum_L2": thresholds["pde_residual_L2"],
        "residual_defined_or_pointwise_free_force_allowed": False,
        "amplitude_collapse_success_route_allowed": False,
        "posthoc_threshold_relaxation_allowed": False,
        "_forcing_restriction": forcing["restriction"],
        "_nontriviality": nontriviality,
        "_failure_policy": validation["failure_policy"],
    }


def _default_a5_payload() -> dict[str, Any]:
    payload = a5.deterministic_physical_center_ingest_contract(
        exact_head=AUDITED_PARENT_HEAD
    )
    a5.validate_physical_center_ingest_contract(payload)
    return payload


def validate_normalization_scope_contract(
    contract: Mapping[str, Any],
    constraints: Mapping[str, Any],
    *,
    a5_payload: Mapping[str, Any] | None = None,
) -> None:
    """Fail closed on normalization-role, state, or CR001 contract laundering."""

    _expect_equal(contract.get("schema"), SCHEMA, "schema")
    _expect_equal(
        contract.get("audited_parent"),
        {
            "pr": AUDITED_PARENT_PR,
            "head": AUDITED_PARENT_HEAD,
            "task": a5.TASK,
        },
        "audited parent",
    )
    _expect_equal(
        contract.get("active_integration_at_freeze"),
        {
            "branch": "codex/cr001-constraints",
            "head": ACTIVE_INTEGRATION_HEAD_AT_FREEZE,
        },
        "active integration freeze",
    )

    classification = contract.get("classification")
    if not isinstance(classification, Mapping):
        raise ValueError("classification must be a mapping")
    _expect_equal(
        set(classification),
        {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "four-way classification vocabulary",
    )
    public_source = classification["public_source_fact"]
    autonomous = classification["autonomous_design"]
    pending = classification["pending_unknown"]
    if not any("F=(phi_*/C)Phi" in item for item in public_source):
        raise ValueError("public source F/C mapping missing")
    if not any("C >= sup" in item for item in public_source):
        raise ValueError("public source complex-C condition missing")
    if not any("configured numerical C" in item for item in autonomous):
        raise ValueError("current configured C must remain autonomous")
    if not any("configured C satisfies" in item for item in pending):
        raise ValueError("configured-C complex-domain status must remain pending")

    source_c = contract.get("source_C_role")
    if not isinstance(source_c, Mapping):
        raise ValueError("source_C_role must be a mapping")
    _expect_equal(source_c.get("formula"), "F=(phi_*/C)Phi", "source C formula")
    _expect_equal(
        source_c.get("source_condition"),
        "C >= sup_{eta in Omega}|phi_*(eta)|",
        "source C condition",
    )
    for key in (
        "current_complex_domain_certified",
        "current_physical_center_ingest_admitted",
        "is_CR001_energy_normalizer",
        "may_be_tuned_posthoc_to_reduce_NS_residual",
        "may_be_inflated_posthoc_to_shrink_F_or_velocity_amplitude",
        "may_be_selected_from_CR001_heldout_validation",
        "unresolved_source_normalization_authorizes_free_amplitude_tuning",
    ):
        _expect_false(source_c, key)

    energy = contract.get("cr001_energy_normalization_role")
    if not isinstance(energy, Mapping):
        raise ValueError("cr001_energy_normalization_role must be a mapping")
    _expect_equal(energy.get("reference_time"), 0.25, "energy reference time")
    _expect_equal(energy.get("reference_energy"), 1.0, "energy reference")
    _expect_equal(energy.get("absolute_tolerance"), 0.001, "energy tolerance")
    _expect_equal(energy.get("minimum_energy_each_validation_time"), 0.1, "minimum energy")
    _expect_equal(energy.get("maximum_energy_each_validation_time"), 10.0, "maximum energy")
    if energy.get("applies_to_complete_candidate_acceptance") is not True:
        raise ValueError("CR001 energy normalization must remain a complete-candidate acceptance rule")
    for key in (
        "physical_center_profile_energy_assessed_here",
        "certifies_source_complex_C_condition",
        "may_replace_source_C_condition",
    ):
        _expect_false(energy, key)

    projection = _constraints_projection(constraints)
    frozen = contract.get("cr001_frozen_contract")
    if not isinstance(frozen, Mapping):
        raise ValueError("cr001_frozen_contract must be a mapping")
    for key, expected in frozen.items():
        if key not in projection:
            raise ValueError(f"unexpected CR001 frozen key {key}")
        _expect_equal(projection[key], expected, f"CR001 {key}")

    restriction = projection["_forcing_restriction"]
    if "No residual-dependent basis or pointwise free force" not in restriction:
        raise ValueError("CR001 free-force prohibition drift")
    nontriviality = projection["_nontriviality"]
    _expect_equal(nontriviality["reference_time"], 0.25, "CR001 nontriviality time")
    _expect_equal(nontriviality["reference_energy"], 1.0, "CR001 nontriviality energy")
    _expect_equal(
        nontriviality["reference_energy_abs_tolerance"], 0.001, "CR001 nontriviality tolerance"
    )
    if "reject collapsed candidates" not in nontriviality["enforcement"]:
        raise ValueError("CR001 amplitude-collapse rejection drift")
    if "changing thresholds requires a new experiment version" not in projection["_failure_policy"]:
        raise ValueError("CR001 threshold failure policy drift")

    _expect_equal(
        a5.SOURCE_FORMULAS.get("physical_swirl_profile"),
        "F=(phi_*/C)Phi",
        "A5 registered source-C formula",
    )
    if a5.SOURCE_COMPLEX_C_NORMALIZATION_CERTIFIED is not False:
        raise ValueError("A5 source complex-C certification unexpectedly promoted")

    live = dict(_default_a5_payload() if a5_payload is None else a5_payload)
    status = live.get("ingest_status")
    truth = live.get("truth_boundary")
    if not isinstance(status, Mapping) or not isinstance(truth, Mapping):
        raise ValueError("A5 payload missing status/truth boundary")
    if status.get("typed_physical_center_profile_interface_registered") is not True:
        raise ValueError("A5 physical-center interface registration missing")
    for key in (
        "physical_center_profile_ingest_admitted",
        "source_complex_C_normalization_certified",
        "fixed_point_correction_materialized",
        "global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_composite_validated",
        "leading_ready",
    ):
        _expect_false(status, key)
    for key in (
        "source_complex_C_normalization_claimed",
        "source_native_center_promoted_to_cartesian_velocity",
        "free_residual_defined_forcing_allowed",
        "threshold_relaxed",
        "velocity_export_ready",
        "pde_validated",
    ):
        _expect_false(truth, key)

    states = contract.get("state_separation")
    if not isinstance(states, Mapping):
        raise ValueError("state_separation must be a mapping")
    if states.get("physical_center_formula_interface_registered") is not True:
        raise ValueError("formula registration must remain true")
    for key in (
        "physical_center_profile_ingest_admitted",
        "kokuno_global_cartesian_velocity_materialized",
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _expect_false(states, key)

    boundary = contract.get("claim_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("claim_boundary must be a mapping")
    for key in (
        "source_C_normalization_failure_or_pending_does_not_retroactively_invalidate_other_honestly_labeled_callable_repository_candidates",
        "a_future_source_C_certification_alone_does_not_establish_CR001_energy_normalization",
        "a_future_CR001_energy_normalization_pass_alone_does_not_establish_source_C_certification",
        "no_state_promotion_in_this_contract",
    ):
        if boundary.get(key) is not True:
            raise ValueError(f"claim boundary {key} drift")


def audit_normalization_scope(
    *,
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    a5_payload: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if contract is None:
        contract = _load_json(_repo_root() / CONTRACT_RELATIVE_PATH)
    if constraints is None:
        constraints = _load_json(_repo_root() / CONSTRAINTS_RELATIVE_PATH)
    validate_normalization_scope_contract(contract, constraints, a5_payload=a5_payload)
    return {
        "schema": SCHEMA,
        "audited_parent_pr": AUDITED_PARENT_PR,
        "audited_parent_head": AUDITED_PARENT_HEAD,
        "contract_sha256": _canonical_digest(contract),
        "cr001_constraints_sha256": _canonical_digest(constraints),
        "source_C_role_separated_from_CR001_energy_normalization": True,
        "source_complex_C_normalization_certified": False,
        "physical_center_profile_ingest_admitted": False,
        "kokuno_velocity_export_ready": False,
        "pde_validated": False,
        "audit_passed": True,
    }


def main() -> int:
    receipt = audit_normalization_scope()
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
