"""Fail-closed CR002 audit for the CR001/Kokuno ``h=0.005`` namespace collision.

The current canonical CR001 contract and the current Kokuno A1 lineage happen to
use the same symbol family and the same numerical value. That coincidence is
not a semantic parameter identity and must not transfer source, delivery, or
scientific evidence between the lineages.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .kokuno_pa16_current_cartesian_postswirl_lminus1_hold import (
    KokunoPA16CurrentCartesianPostSwirlLMinus1Hold,
)

CONTRACT_PATH = "configs/kokuno_cr001_h_namespace_scope.json"
_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _fail(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parameter_namespace_key(record: Mapping[str, Any]) -> tuple[str, str]:
    """Return semantic namespace key, deliberately excluding the scalar value."""
    return str(record["namespace"]), str(record["symbol"])


def _load_contract(root: Path) -> dict[str, Any]:
    return json.loads((root / CONTRACT_PATH).read_text())


def audit_contract(
    repo_root: str | Path | None = None,
    *,
    contract_override: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit the cross-lineage namespace firewall against exact parent bytes."""

    root = Path(repo_root) if repo_root is not None else _repo_root()
    contract = (
        copy.deepcopy(dict(contract_override))
        if contract_override is not None
        else _load_contract(root)
    )

    _fail(contract.get("schema_version") == 1, "unexpected contract schema")
    _fail(
        set(contract.get("provenance_classes", {})) == _CANONICAL_CLASSES,
        "four-class provenance vocabulary drift",
    )
    _fail(
        contract.get("exact_parent", {}).get("head")
        == "e0cf691884fd371c477de5ae55a18cdabbd72e81",
        "A1 #1196 exact-parent drift",
    )

    for name, pinned in contract["pinned_files"].items():
        path = root / pinned["path"]
        _fail(path.is_file(), f"missing pinned file: {name}")
        actual = _git_blob_sha1(path)
        _fail(
            actual == pinned["git_blob_sha1"],
            f"pinned Git blob drift for {name}: {actual}",
        )

    constraints = json.loads((root / "configs/constraints.json").read_text())
    status = json.loads((root / "project_status.json").read_text())
    freeze = contract["canonical_cr001_freeze"]

    _fail(constraints["nu"] == freeze["nu"] == 0.01, "canonical nu drift")
    _fail(
        constraints["domain"]["physical"] == freeze["physical_domain"] == "R^3",
        "canonical physical domain drift",
    )
    _fail(
        constraints["domain"]["evaluation_box"] == freeze["evaluation_box"],
        "canonical evaluation box drift",
    )
    _fail(
        constraints["domain"]["support"] == freeze["support"],
        "canonical support drift",
    )
    _fail(
        constraints["domain"]["time_interval"] == freeze["time_interval"],
        "canonical time interval drift",
    )
    _fail(
        constraints["forcing"]["mode"] == freeze["forcing_mode"]
        == "restricted_two_parameter_family",
        "canonical forcing family drift",
    )
    _fail(
        freeze["residual_defined_free_force_allowed"] is False
        and "No residual-dependent basis or pointwise free force"
        in constraints["forcing"]["restriction"],
        "residual-defined forcing firewall drift",
    )
    _fail(
        constraints["nontriviality"]["reference_energy"] == freeze["reference_energy"]
        and constraints["nontriviality"]["reference_energy_abs_tolerance"]
        == freeze["reference_energy_abs_tolerance"],
        "canonical nontriviality drift",
    )
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    _fail(validation["seed"] == freeze["validation_seed"] == 914027, "validation seed drift")
    _fail(
        validation["held_out_points"] == freeze["held_out_points"] == 4096,
        "held-out sample count drift",
    )
    _fail(
        validation["derivative_steps"] == freeze["derivative_steps"] == [0.02, 0.01, 0.005],
        "canonical derivative ladder drift",
    )
    _fail(
        validation["quadrature_orders_per_axis"]
        == freeze["quadrature_orders_per_axis"]
        == [24, 48, 96],
        "canonical quadrature ladder drift",
    )
    for key, expected in (
        ("pde_residual_max", 0.001),
        ("pde_residual_L2", 0.001),
        ("divergence_max", 1e-5),
        ("divergence_L2", 1e-5),
    ):
        _fail(
            thresholds[key] == freeze[key] == expected,
            f"canonical threshold drift: {key}",
        )

    namespaces = contract["parameter_namespaces"]
    cr001_h = namespaces["canonical_cr001_h"]
    kokuno_h = namespaces["kokuno_a1_h_current"]

    _fail(cr001_h["classification"] == "autonomous_design", "CR001 h provenance laundering")
    _fail(kokuno_h["classification"] == "autonomous_design", "Kokuno numeric h provenance laundering")
    _fail(
        kokuno_h["public_relation_classification"] == "public_source_fact",
        "Kokuno public relation classification drift",
    )
    _fail(cr001_h["source_exact_numeric_value"] is False, "CR001 h cannot become source-exact here")
    _fail(kokuno_h["source_exact_numeric_value"] is False, "Kokuno h cannot become source-exact here")

    actual_cr001_h = float(constraints["structure"]["h"])
    _fail(
        math.isclose(actual_cr001_h, float(cr001_h["current_value"]), rel_tol=0.0, abs_tol=0.0),
        "canonical CR001 h value drift",
    )

    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    actual_kokuno_h = float(field.h_value)
    _fail(
        math.isclose(actual_kokuno_h, float(kokuno_h["current_value"]), rel_tol=0.0, abs_tol=2e-15),
        "Kokuno current h value drift",
    )
    _fail(field.truth_boundary["source_exact_h_recovered"] is False, "Kokuno source-exact h promotion")
    _fail(
        field.truth_boundary["unified_global_cartesian_velocity_export_ready"] is False,
        "Kokuno global-export promotion",
    )
    _fail(field.truth_boundary["pde_validated"] is False, "Kokuno PDE promotion")

    firewall = contract["identity_firewall"]
    _fail(
        firewall == {
            "same_numeric_value_observed": True,
            "shared_parameter_identity": False,
            "shared_source_provenance": False,
            "automatic_value_synchronization": False,
            "canonical_cr001_evidence_transfers_to_kokuno": False,
            "kokuno_evidence_transfers_to_canonical_cr001": False,
            "numeric_equality_can_promote_velocity_visual_pde_or_exactness_state": False,
        },
        "cross-lineage identity/evidence firewall drift",
    )
    _fail(
        parameter_namespace_key(cr001_h) != parameter_namespace_key(kokuno_h),
        "independent h namespaces collapsed",
    )
    _fail(
        math.isclose(actual_cr001_h, actual_kokuno_h, rel_tol=0.0, abs_tol=2e-15),
        "registered same-value witness no longer holds",
    )

    eq45 = contract["independent_state_boundary"]["canonical_eq45"]
    _fail(
        status["states"]["velocity_export_ready"] is eq45["velocity_export_ready"] is True,
        "canonical callable delivery demoted",
    )
    for key in (
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _fail(status["states"][key] is eq45[key] is False, f"canonical state promotion: {key}")

    kokuno_state = contract["independent_state_boundary"]["kokuno_a1_lminus1_hold"]
    _fail(kokuno_state["stage_callable_velocity_materialized"] is True, "stage callability demoted")
    for key in (
        "unified_global_cartesian_velocity_export_ready",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _fail(kokuno_state[key] is False, f"Kokuno state promotion: {key}")

    return {
        "contract_id": contract["contract_id"],
        "canonical_cr001_h": actual_cr001_h,
        "kokuno_a1_h_current": actual_kokuno_h,
        "same_numeric_value_observed": True,
        "shared_parameter_identity": False,
        "canonical_eq45_velocity_export_ready": True,
        "kokuno_global_velocity_export_ready": False,
        "pde_validated": False,
    }
