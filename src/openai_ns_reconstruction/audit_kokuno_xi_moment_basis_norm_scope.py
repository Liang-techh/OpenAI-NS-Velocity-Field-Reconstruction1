"""Fail-closed CR002 audit for the Kokuno Xi five-moment basis/norm seam.

This module governs representation identity only.  It does not evaluate or modify
any velocity, pressure, forcing, candidate coefficient, residual, or scientific
threshold.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-xi-moment-basis-norm-scope-v1"
TASK_ID = "CR002-KOKUNO-XI-MOMENT-BASIS-NORM-090"
DEFAULT_SCOPE = Path("configs/kokuno_xi_moment_basis_norm_scope.json")
DEFAULT_CONSTRAINTS = Path("configs/constraints.json")
PARENT_SOURCE = Path(
    "src/openai_ns_reconstruction/"
    "kokuno_a5_current_xi_moment_discrepancy_ingest_contract.py"
)
EXPECTED_PARENT_BLOB = "52562076ed2b54ea069a55d245ede4b244c5d2b9"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

EXPECTED_RAW_ROWS = ["M", "I", "J", "S", "C_p"]
EXPECTED_PA16_ROWS = ["M", "J-4*eta*I", "I", "S-8*eta*M", "C_p"]
REQUIRED_FUTURE_BINDINGS = {
    "one checksum-bound moment_basis_id shared by discrepancy_d, matrix_B and bilinear_Q",
    "exact row ordering and PA.16 transform identity",
    "exact PA.15/source-coordinate scaling and current outer-normalization identity",
    "one explicit vector/matrix/bilinear norm convention",
    "proof or machine check that any row/basis transformation is applied consistently to B, Q and d",
    "discrepancy_from_complete_ns_defect=true on the same complete candidate/physical contract",
    "no residual-as-forcing shortcut and no held-out leakage",
}
REQUIRED_FORBIDDEN_PROMOTIONS = {
    "small current Xi discrepancy norm => source smallness condition passed",
    "A4 Xi moment audit pass => gain-stage authorization",
    "rescale discrepancy_d alone => same five-moment correction system",
    "change PA.15 or PA.16 row normalization without transforming B and Q => same five-moment correction system",
    "source-coordinate moment discrepancy => complete NS defect",
    "local five-moment contraction => PDE validation",
}


class MomentBasisNormScopeError(RuntimeError):
    """Raised when the governance contract drifts or is promoted unsafely."""


def _git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise MomentBasisNormScopeError(f"{path} must contain a JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise MomentBasisNormScopeError(message)


def _smallness(B: float, Q: float, d: float) -> float:
    nu = abs(1.0 / B)
    q = abs(Q)
    return 8.0 * nu * nu * q * abs(d)


def _audit_mechanics_witness(scope: Mapping[str, Any]) -> None:
    witness = scope["mechanics_witness"]
    _require(witness["classification"] == "autonomous mechanics only",
             "mechanics witness must remain autonomous-only")
    _require(witness["not_source_or_candidate_data"] is True,
             "mechanics witness must not be promoted to source/candidate data")

    original = witness["scalar_original"]
    consistent = witness["consistent_common_row_scale"]
    inconsistent = witness["inconsistent_d_only_scale"]

    original_value = _smallness(original["B"], original["Q"], original["d"])
    consistent_value = _smallness(consistent["B"], consistent["Q"], consistent["d"])
    inconsistent_value = _smallness(inconsistent["B"], inconsistent["Q"], inconsistent["d"])

    _require(abs(original_value - 1.6) <= 1e-12, "original mechanics witness drift")
    _require(abs(consistent_value - original_value) <= 1e-12,
             "consistent row scaling must preserve the scalar gain parameter")
    _require(abs(inconsistent_value - 0.16) <= 1e-12,
             "d-only scaling witness drift")
    _require(inconsistent_value < 1.0 < original_value,
             "d-only scaling witness must demonstrate a false pass possibility")


def _audit_cr001(scope: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    frozen = scope["frozen_cr001"]
    _require(constraints["nu"] == frozen["nu"] == 0.01, "nu drift")
    domain = constraints["domain"]
    _require(domain["physical"] == frozen["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == frozen["evaluation_box"], "evaluation box drift")
    _require(domain["support"] == frozen["support"], "support drift")
    _require(domain["time_interval"] == frozen["time_interval"], "time interval drift")
    _require(constraints["forcing"]["mode"] == frozen["forcing_mode"], "forcing family drift")
    _require(
        "No residual-dependent basis or pointwise free force"
        in constraints["forcing"]["restriction"],
        "residual-defined free-force prohibition missing",
    )
    nontriviality = constraints["nontriviality"]
    _require(nontriviality["reference_energy"] == frozen["reference_energy"], "energy target drift")
    _require(
        nontriviality["reference_energy_abs_tolerance"]
        == frozen["reference_energy_abs_tolerance"],
        "energy tolerance drift",
    )
    validation = constraints["validation"]
    _require(validation["seed"] == frozen["validation_seed"], "validation seed drift")
    _require(validation["held_out_points"] == frozen["held_out_points"], "held-out count drift")
    _require(validation["derivative_steps"] == frozen["derivative_steps"], "derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == frozen["quadrature_orders_per_axis"],
        "quadrature ladder drift",
    )
    gates = validation["thresholds"]
    _require(gates["pde_residual_max"] == frozen["momentum_max_gate"] == 1e-3,
             "momentum max gate drift")
    _require(gates["pde_residual_L2"] == frozen["momentum_l2_gate"] == 1e-3,
             "momentum L2 gate drift")
    _require(gates["divergence_max"] == frozen["divergence_max_gate"] == 1e-5,
             "divergence max gate drift")
    _require(gates["divergence_L2"] == frozen["divergence_l2_gate"] == 1e-5,
             "divergence L2 gate drift")
    _require(frozen["residual_defined_pointwise_free_force_forbidden"] is True,
             "free-force firewall disabled")
    _require(frozen["collapsed_candidate_forbidden"] is True, "collapse firewall disabled")
    _require(frozen["post_hoc_threshold_relaxation_forbidden"] is True,
             "threshold-relaxation firewall disabled")


def audit_scope(
    scope_path: str | Path = DEFAULT_SCOPE,
    constraints_path: str | Path = DEFAULT_CONSTRAINTS,
    *,
    require_parent_blob: bool = True,
) -> dict[str, Any]:
    scope_path = Path(scope_path)
    constraints_path = Path(constraints_path)
    scope = _load_json(scope_path)
    constraints = _load_json(constraints_path)

    _require(scope["schema"] == SCHEMA, "schema drift")
    _require(scope["task_id"] == TASK_ID, "task id drift")

    pinned = scope["pinned_inputs"]
    _require(
        pinned["parent_a5"]["head"] == "b3fb36d91239a2ce2669971b7b07bfe892cfdf1c",
        "parent A5 head drift",
    )
    _require(pinned["parent_a5"]["source_blob"] == EXPECTED_PARENT_BLOB,
             "parent A5 source-blob declaration drift")
    _require(
        pinned["canonical_constraints"]["blob"] == EXPECTED_CONSTRAINTS_BLOB,
        "canonical constraints-blob declaration drift",
    )
    if require_parent_blob:
        _require(PARENT_SOURCE.exists(), "pinned parent A5 source is missing")
        _require(_git_blob_sha(PARENT_SOURCE) == EXPECTED_PARENT_BLOB,
                 "pinned parent A5 source blob mismatch")
        _require(_git_blob_sha(constraints_path) == EXPECTED_CONSTRAINTS_BLOB,
                 "canonical constraints blob mismatch")

    xi = pinned["current_xi_discrepancy"]
    _require(xi["raw_rows"] == EXPECTED_RAW_ROWS, "raw moment row order drift")
    _require(xi["pa16_rows"] == EXPECTED_PA16_ROWS, "PA.16 row transform drift")
    _require(xi["source_normalization_C"] == 1000.0,
             "current autonomous outer-normalization identity drift")

    bridge = pinned["current_xi_row_transform"]
    _require(bridge["discrepancy_from_complete_ns_defect"] is False,
             "source moment discrepancy laundered into complete-NS defect")
    _require(bridge["authorized_for_gain_gated_ns_stage"] is False,
             "current Xi discrepancy prematurely authorized for gain stage")

    gain = pinned["existing_gain_firewall"]
    _require(gain["source_blob"] == "6625db71ec8478db5d6ec09cabb7033a03e658a4",
             "gain-firewall source identity drift")
    _require(gain["basis_or_row_normalization_identity_field_present"] is False,
             "basis-binding gap silently rewritten")
    _require(
        gain["norm_convention"]
        == "Euclidean vectors / spectral matrix inverse / slice-spectral bilinear upper bound",
        "gain norm convention drift",
    )

    contract = scope["moment_coordinate_contract"]
    _require(contract["current_discrepancy_coordinate_surface_registered"] is True,
             "current discrepancy surface lost")
    _require(contract["current_gain_firewall_uses_euclidean_norms"] is True,
             "gain norm truth drift")
    for key in (
        "current_gain_firewall_schema_explicitly_binds_moment_basis_id",
        "current_gain_firewall_schema_explicitly_binds_row_scale_or_pa15_normalization",
        "current_gain_firewall_schema_explicitly_binds_pa16_transform_identity",
        "current_gain_firewall_schema_explicitly_binds_outer_normalization_C",
        "current_xi_discrepancy_authorized_for_gain_stage_now",
    ):
        _require(contract[key] is False, f"unsafe promotion: {key}")
    _require(set(contract["future_gain_authorization_requires"]) == REQUIRED_FUTURE_BINDINGS,
             "future basis-binding requirements drift")
    _require(set(contract["forbidden_promotions"]) == REQUIRED_FORBIDDEN_PROMOTIONS,
             "forbidden-promotion firewall drift")

    classes = scope["provenance_classes"]
    _require(
        "current outer normalization C=1000 is a repository realization, not a recovered hidden source value"
        in classes["autonomous_design"],
        "autonomous C=1000 provenance lost",
    )
    _require(
        "the exact moment-coordinate/basis/row-normalization identity tying d, B and Q together is not yet bound"
        in classes["pending_unknown"],
        "pending basis identity silently promoted",
    )

    delivery = scope["delivery_state"]
    _require(delivery["canonical_eq45_velocity_export_ready"] is True,
             "unresolved Kokuno basis work must not cancel Eq45 delivery")
    for key in (
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(delivery[key] is False, f"unsafe scientific promotion: {key}")

    _audit_mechanics_witness(scope)
    _audit_cr001(scope, constraints)

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "status": "PASS",
        "scope": "five-moment coordinate/basis/norm governance only",
        "current_xi_gain_authorized": False,
        "pde_validated": False,
    }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", default=str(DEFAULT_SCOPE))
    parser.add_argument("--constraints", default=str(DEFAULT_CONSTRAINTS))
    args = parser.parse_args()
    print(json.dumps(audit_scope(args.scope, args.constraints), sort_keys=True))


if __name__ == "__main__":
    _main()
