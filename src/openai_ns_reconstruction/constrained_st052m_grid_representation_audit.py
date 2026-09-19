"""CR002 fail-closed audit for the ST052-M sampled-grid representation scope.

The integrated CR-A9-071 grid exporter is a useful save/load and visualization
adapter for the authenticated whole-child velocity callable.  A finite set of
serialized grid samples is nevertheless not, by itself, an identity for the
continuous function ``velocity(x,y,z,t)``.  This audit keeps those two facts
separate without adding a scientific/PDE gate to callable velocity delivery.
"""
from __future__ import annotations

import inspect
import json
import math
from pathlib import Path
from typing import Any

from . import st052_grid_export as grid_export
from . import st052_linear_temporal_capsule as temporal_capsule


_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT_PATH = _REPO_ROOT / "configs" / "st052m_sampled_grid_representation_contract.json"


class GridRepresentationAuditError(RuntimeError):
    """Raised when sampled-grid scope is promoted beyond verified evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GridRepresentationAuditError(message)


def load_contract(path: str | Path = DEFAULT_CONTRACT_PATH) -> dict[str, Any]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    _require(isinstance(data, dict), "grid representation contract must be a JSON object")
    return data


def sampled_grid_nonuniqueness_witness() -> dict[str, float]:
    """Return a deterministic witness that frozen grid nodes do not fix off-grid values.

    For the default uniform x-grid, ``delta(x)=sin(pi*(x-x_min)/h)`` vanishes
    at every serialized x node but equals one at the first half-grid point.
    Extending this perturbation constantly in y/z/t and adding it to one
    component produces a distinct continuous field with identical serialized
    grid values.  The witness is representation-only and is never applied to a
    candidate.
    """
    lo, hi = map(float, grid_export.DEFAULT_DOMAIN)
    n = int(grid_export.DEFAULT_GRID_POINTS)
    h = (hi - lo) / float(n - 1)
    node_values = [math.sin(math.pi * ((lo + i * h) - lo) / h) for i in range(n)]
    midpoint = lo + 0.5 * h
    midpoint_value = math.sin(math.pi * (midpoint - lo) / h)
    return {
        "grid_spacing": h,
        "max_abs_at_serialized_nodes": max(abs(v) for v in node_values),
        "abs_at_first_half_grid_point": abs(midpoint_value),
    }


def audit_grid_representation_contract(
    path: str | Path = DEFAULT_CONTRACT_PATH,
) -> dict[str, Any]:
    """Audit the machine-readable scope against the integrated live adapter."""
    contract = load_contract(path)
    _require(contract.get("schema_version") == 1, "unexpected contract schema")
    _require(
        contract.get("contract_id") == "st052m_sampled_grid_representation_scope_v1",
        "unexpected contract id",
    )

    classification = contract.get("classification")
    _require(isinstance(classification, dict), "classification block is required")
    for key in ("user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"):
        value = classification.get(key)
        _require(isinstance(value, list) and len(value) > 0, f"classification {key} must be nonempty")
        _require(all(isinstance(item, str) and item.strip() for item in value), f"classification {key} contains an invalid entry")

    representation = contract.get("representation")
    _require(isinstance(representation, dict), "representation block is required")
    _require(
        representation.get("grid_shape_xyz") == [grid_export.DEFAULT_GRID_POINTS] * 3,
        "contract grid shape drifted from the live exporter",
    )
    _require(
        representation.get("domain") == list(grid_export.DEFAULT_DOMAIN),
        "contract domain drifted from the live exporter",
    )
    _require(
        representation.get("times") == list(grid_export.DEFAULT_TIMES),
        "contract times drifted from the live exporter",
    )
    _require(
        representation.get("callable_grid_node_replay_atol") == grid_export.CALLABLE_GRID_PARITY_ATOL,
        "contract callable-grid replay tolerance drifted from the live exporter",
    )
    _require(
        representation.get("callable_grid_node_replay_tolerance_role")
        == "engineering_serialization_replay_only",
        "the 5e-12 software replay tolerance must not be relabelled as a scientific threshold",
    )
    _require(
        grid_export.callable_grid_parity_passes(grid_export.CALLABLE_GRID_PARITY_ATOL),
        "live exporter must accept its frozen software replay tolerance",
    )
    _require(
        not grid_export.callable_grid_parity_passes(
            math.nextafter(grid_export.CALLABLE_GRID_PARITY_ATOL, math.inf)
        ),
        "live exporter must fail closed immediately above the frozen software replay tolerance",
    )

    truth = contract.get("truth_boundary")
    _require(isinstance(truth, dict), "truth_boundary block is required")
    required_true = (
        "whole_child_runtime_callable_exists_with_authenticated_source_runtime",
        "underlying_runtime_callable_remains_a_separate_delivery_representation",
        "sampled_grid_export_path_integrated",
        "grid_payload_bound_to_whole_candidate_identity",
        "npz_mat_payload_exact_equality_verified_on_exact_ci_receipt",
        "callable_values_replayed_at_serialized_grid_nodes",
    )
    required_false = (
        "sampled_grid_is_complete_continuous_velocity_representation",
        "grid_export_alone_satisfies_continuous_velocity_api",
        "off_grid_interpolant_specified_and_identity_bound",
        "off_grid_interpolation_error_bound_verified",
        "continuous_derivative_equivalence_verified",
        "pde_validation_transfers_to_grid_interpolant",
        "pde_validated_for_grid_interpolant",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    )
    for key in required_true:
        _require(truth.get(key) is True, f"verified delivery fact must remain true: {key}")
    for key in required_false:
        _require(truth.get(key) is False, f"unverified sampled-grid claim must remain false: {key}")

    _require(callable(getattr(temporal_capsule, "load_bundle_runtime", None)), "whole-child runtime loader is missing")
    export_source = inspect.getsource(grid_export)
    _require(
        "candidate.velocity(xx, yy, zz, float(time))" in export_source,
        "grid export no longer visibly samples the underlying candidate velocity callable",
    )
    _require(
        "This is a delivery adapter, not a scientific validator." in export_source,
        "live exporter lost its delivery-adapter truth boundary",
    )
    _require(
        "It is not a PDE, visualization, or source-correspondence acceptance gate." in export_source,
        "live exporter lost the non-scientific replay-tolerance truth boundary",
    )

    witness = sampled_grid_nonuniqueness_witness()
    _require(
        witness["max_abs_at_serialized_nodes"] < 1e-12,
        "nonuniqueness witness must vanish numerically at all frozen x-grid nodes",
    )
    _require(
        witness["abs_at_first_half_grid_point"] > 1.0 - 1e-12,
        "nonuniqueness witness must differ off grid",
    )

    return {
        "contract_id": contract["contract_id"],
        "grid_shape_xyz": representation["grid_shape_xyz"],
        "times": representation["times"],
        "callable_grid_node_replay_atol": representation["callable_grid_node_replay_atol"],
        "nonuniqueness_witness": witness,
        "sampled_grid_is_complete_continuous_velocity_representation": False,
        "underlying_runtime_callable_delivery_blocked_by_this_audit": False,
    }
