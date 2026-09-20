"""Fail-closed CR002 audit for Kokuno inner convective consistency metric scope.

The upstream Agent-4 quantity is a relative implementation-consistency error for
one strict-inner convective subterm.  It is intentionally not interchangeable
with the absolute complete-momentum residual norms preregistered by CR001.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "cr002-kokuno-inner-convective-consistency-norm-scope-v1"
CONTRACT_PATH = "configs/kokuno_inner_convective_consistency_norm_scope.json"
UPSTREAM_PATH = (
    "src/openai_ns_reconstruction/"
    "kokuno_a4_inner_leading_oscillatory_full_advection_independent_audit.py"
)
UPSTREAM_BLOB = "5de0114405f1907454d5f0da7616e15dba29bd4e"

_REQUIRED_CLASSIFICATIONS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_REQUIRED_UPSTREAM_CONSTANTS = {
    "FD6_STEPS": (1.6e-3, 8.0e-4, 4.0e-4),
    "FINE_RELATIVE_RMS_GATE": 2.0e-3,
    "FINE_RELATIVE_MAX_GATE": 5.0e-3,
    "AGENT2_FULL_ADVECTION_PR": 830,
    "AGENT2_FULL_ADVECTION_HEAD": "3479bed94ac6ef5c67547227f81a224ad43fb409",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object at {path}")
    return value


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()  # noqa: S324 - Git object identity


def _literal_module_constants(source: str) -> dict[str, Any]:
    tree = ast.parse(source)
    values: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name = node.targets[0].id
            value_node = node.value
        else:
            if not isinstance(node.target, ast.Name) or node.value is None:
                continue
            name = node.target.id
            value_node = node.value
        try:
            values[name] = ast.literal_eval(value_node)
        except (ValueError, TypeError):
            continue
    return values


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _validate_contract(
    contract: dict[str, Any],
    cr001: dict[str, Any],
    upstream_constants: dict[str, Any],
) -> None:
    _require(contract.get("schema") == SCHEMA, "scope contract schema drift")
    _require(
        set(contract.get("classification", {})) == _REQUIRED_CLASSIFICATIONS,
        "four-way provenance classification drift",
    )

    upstream = contract["upstream"]
    _require(upstream["pr"] == 832, "upstream PR drift")
    _require(
        upstream["head"] == "2cb98c5e7dc97d17af650c673ad70da0f9ac3c2c",
        "upstream head drift",
    )
    _require(upstream["module"] == UPSTREAM_PATH, "upstream module drift")
    _require(upstream["module_blob"] == UPSTREAM_BLOB, "upstream blob drift")

    for name, expected in _REQUIRED_UPSTREAM_CONSTANTS.items():
        _require(upstream_constants.get(name) == expected, f"upstream constant drift: {name}")

    metric = contract["upstream_consistency_metric"]
    _require(metric["fine_relative_rms_gate"] == 2.0e-3, "relative RMS gate drift")
    _require(metric["fine_relative_sampled_max_gate"] == 5.0e-3, "relative max gate drift")
    _require(metric["dimensionless_relative_error"] is True, "relative metric semantics drift")
    for false_key in (
        "complete_momentum_residual",
        "includes_velocity_dt",
        "includes_pressure_gradient",
        "includes_viscous_laplacian",
        "includes_restricted_forcing",
        "includes_correction_velocity",
        "canonical_cr001_momentum_norm",
        "directly_comparable_to_st006_full_residual",
    ):
        _require(metric[false_key] is False, f"metric-scope laundering: {false_key}")

    canonical = contract["canonical_cr001"]
    validation = cr001["validation"]
    forcing = cr001["forcing"]
    nontriviality = cr001["nontriviality"]
    domain = cr001["domain"]

    _require(cr001["nu"] == canonical["nu"] == 0.01, "nu drift")
    _require(domain["physical"] == canonical["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == canonical["evaluation_box"], "evaluation box drift")
    _require(domain["support"] == canonical["support"], "support drift")
    _require(domain["time_interval"] == canonical["time_interval"], "time interval drift")
    _require(
        forcing["mode"] == canonical["forcing_mode"] == "restricted_two_parameter_family",
        "forcing mode drift",
    )
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "free-force prohibition missing",
    )
    _require(canonical["free_residual_defined_force_allowed"] is False, "free residual force enabled")
    _require(nontriviality["reference_energy"] == canonical["reference_energy"] == 1.0, "energy target drift")
    _require(
        nontriviality["reference_energy_abs_tolerance"]
        == canonical["reference_energy_abs_tolerance"]
        == 1.0e-3,
        "energy tolerance drift",
    )
    _require("reject collapsed candidates" in nontriviality["enforcement"], "amplitude collapse guard missing")
    _require(canonical["reject_amplitude_collapse"] is True, "amplitude collapse allowed")

    _require(validation["seed"] == canonical["validation_seed"] == 914027, "validation seed drift")
    _require(validation["held_out_points"] == canonical["held_out_points"] == 4096, "held-out count drift")
    _require(validation["times"] == canonical["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == canonical["derivative_steps"], "derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == canonical["quadrature_orders_per_axis"],
        "quadrature ladder drift",
    )
    _require(
        validation["norms"]
        == ["max Euclidean vector norm", "volume-weighted L2 spatial norm at each time"],
        "canonical validation norms drift",
    )
    thresholds = validation["thresholds"]
    _require(thresholds["pde_residual_max"] == canonical["momentum_max_threshold"] == 1.0e-3, "momentum max threshold drift")
    _require(thresholds["pde_residual_L2"] == canonical["momentum_l2_threshold"] == 1.0e-3, "momentum L2 threshold drift")
    _require(thresholds["divergence_max"] == canonical["divergence_max_threshold"] == 1.0e-5, "divergence max threshold drift")
    _require(thresholds["divergence_L2"] == canonical["divergence_l2_threshold"] == 1.0e-5, "divergence L2 threshold drift")
    _require("changing thresholds requires a new experiment version" in validation["failure_policy"], "post-hoc threshold guard missing")
    _require(canonical["posthoc_threshold_relaxation_allowed"] is False, "post-hoc relaxation enabled")

    st006 = contract["st006_reference"]
    _require(st006["momentum_sampled_max"] == 0.1082289305112118, "ST006 max drift")
    _require(st006["momentum_volume_l2"] == 0.10758432876230622, "ST006 L2 drift")
    _require(st006["pde_validated"] is False, "ST006 promoted to PDE validated")
    _require(st006["comparison_from_upstream_consistency_metric_authorized"] is False, "ST006 comparison laundering")

    states = contract["states"]
    _require(states["strict_inner_convective_consistency_registered"] is True, "upstream consistency registration lost")
    _require(states["strict_inner_convective_consistency_ci_passed"] is None, "queued/unresolved CI pre-promoted")
    for false_key in (
        "canonical_cr001_momentum_assessed_by_this_contract",
        "st006_comparable",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(states[false_key] is False, f"truth-state promotion: {false_key}")

    promotion = contract["promotion_boundary"]
    required_nonclaims = {
        "canonical_cr001_momentum_max_acceptance",
        "canonical_cr001_momentum_l2_acceptance",
        "st006_same_protocol_improvement",
        "complete_global_kokuno_candidate",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    }
    _require(
        required_nonclaims.issubset(set(promotion["future_upstream_pass_does_not_establish"])),
        "promotion boundary weakened",
    )
    _require(
        promotion["upstream_failure_invalidates_existing_callable_velocity_delivery"] is False,
        "PDE/subterm failure incorrectly blocks existing callable delivery",
    )


def audit(
    contract_path: Path | None = None,
    cr001_path: Path | None = None,
    upstream_path: Path | None = None,
) -> dict[str, Any]:
    root = _repo_root()
    contract_path = contract_path or root / CONTRACT_PATH
    cr001_path = cr001_path or root / "configs/constraints.json"
    upstream_path = upstream_path or root / UPSTREAM_PATH

    contract = _load_json(contract_path)
    cr001 = _load_json(cr001_path)
    upstream_bytes = upstream_path.read_bytes()
    actual_blob = _git_blob_sha(upstream_bytes)
    _require(actual_blob == UPSTREAM_BLOB, f"upstream source blob drift: {actual_blob}")
    upstream_constants = _literal_module_constants(upstream_bytes.decode("utf-8"))
    _validate_contract(contract, cr001, upstream_constants)
    return {
        "schema": SCHEMA,
        "status": "pass",
        "upstream_blob": actual_blob,
        "scope": "strict_inner_convective_relative_consistency_not_complete_cr001_momentum_residual",
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
