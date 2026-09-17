"""Fail-closed evidence snapshot for the derivative-balanced quartic Phi(1,0) field.

This binds completed sibling receipts to the exact replayable quartic candidate without
claiming those PRs are code ancestry. It changes no velocity value and does not promote
visualization, PDE, paper-exact, OpenAI-field, or blow-up status.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_supported_phi10_quartic_delivery_capsule import (
    EXPECTED_BASE_SHA256,
    EXPECTED_CANDIDATE_SHA256,
    load_recipe,
)

SCHEMA = "eq45_supported_phi10_quartic_evidence_snapshot_v1"
TASK_ID = "CR012-EQ45-PHI10-QUARTIC-EVIDENCE-SNAPSHOT-017"
BASE_DELIVERY_HEAD = "379ea2a72a96b8b396fdfa21f4657e5a91f1f850"
MORPHOLOGY_PR = 170
MORPHOLOGY_HEAD = "2dc27213359a2ae9687bf716f6fd484b30cc30d8"
MORPHOLOGY_ACTIONS_RUN = 35223571477
FRESH_SEED_PR = 171
FRESH_SEED_HEAD = "7b56055b5a6bd7be2e9dad23d060d5f63cb12880"
FRESH_SEED_ACTIONS_RUN = 35223536811
EXPECTED_SEEDS = (914531, 914547, 914563)
EXPECTED_QUARTIC_VS_STATIC_PCT = (-1.4678, -0.5585, -0.7939)
EXPECTED_QUARTIC_VS_CUBIC_PCT = (2.3828, 0.4938, -1.9274)

_SNAPSHOT_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_quartic_evidence_snapshot.json"
)

_EXPECTED_STATUS = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _require_close(name: str, actual: Any, expected: float, *, atol: float = 1.0e-12) -> None:
    value = float(actual)
    if abs(value - expected) > atol:
        raise ValueError(f"{name} drifted: {value!r} != {expected!r}")


def _require_float_tuple(
    name: str, actual: Any, expected: tuple[float, ...], *, atol: float = 1.0e-12
) -> None:
    try:
        values = tuple(float(x) for x in actual)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric sequence") from exc
    if len(values) != len(expected):
        raise ValueError(f"{name} length drifted")
    if any(abs(a - b) > atol for a, b in zip(values, expected)):
        raise ValueError(f"{name} drifted")


def _validate_snapshot(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("quartic evidence snapshot must be an object")
    data = dict(payload)

    if data.get("schema") != SCHEMA:
        raise ValueError("quartic evidence snapshot schema drifted")
    if data.get("task_id") != TASK_ID:
        raise ValueError("quartic evidence task identity drifted")
    if data.get("base_delivery_pr") != 172:
        raise ValueError("base delivery PR drifted")
    if data.get("base_delivery_head") != BASE_DELIVERY_HEAD:
        raise ValueError("base delivery head drifted")
    if data.get("base_supported_sha256") != EXPECTED_BASE_SHA256:
        raise ValueError("supported parent identity drifted")
    if data.get("candidate_sha256") != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("quartic candidate identity drifted")

    morphology = data.get("morphology_evidence")
    if not isinstance(morphology, Mapping):
        raise ValueError("morphology_evidence must be an object")
    if (
        morphology.get("pr") != MORPHOLOGY_PR
        or morphology.get("exact_head") != MORPHOLOGY_HEAD
        or morphology.get("actions_run") != MORPHOLOGY_ACTIONS_RUN
    ):
        raise ValueError("quartic morphology provenance drifted")
    if (
        morphology.get("evidence_class")
        != "target_free_three_resolution_public_velocity_vorticity_morphology"
    ):
        raise ValueError("quartic morphology evidence class drifted")
    _require_close("morphology time", morphology.get("time"), 0.6875)
    if morphology.get("grid_sizes") != [49, 65, 81]:
        raise ValueError("morphology resolution ladder drifted")
    finest = morphology.get("finest_81")
    if not isinstance(finest, Mapping):
        raise ValueError("finest_81 must be an object")
    if finest.get("field_order") != [
        "static_supported",
        "cubic_localized",
        "quartic_derivative_balanced",
    ]:
        raise ValueError("morphology field ordering drifted")
    _require_float_tuple(
        "radial_q99",
        finest.get("radial_q99"),
        (0.8514693182963201, 0.8514693182963199, 0.8514693182963201),
    )
    _require_float_tuple("axial_q99", finest.get("axial_q99"), (0.15, 0.15, 0.15))
    _require_float_tuple(
        "radial_rms",
        finest.get("vorticity2_weighted_radial_rms"),
        (0.5563493225057393, 0.5555062336698936, 0.5564014545337104),
    )
    _require_float_tuple(
        "axial_rms",
        finest.get("vorticity2_weighted_axial_rms"),
        (0.13754962418235045, 0.13749256401032742, 0.13755313769546357),
    )
    _require_float_tuple(
        "weighted_aspect",
        finest.get("weighted_aspect"),
        (0.24723607743933484, 0.24750858888116742, 0.24721922736657712),
    )
    _require_float_tuple(
        "collar_vorticity2_fraction",
        finest.get("whole_grid_support_collar_vorticity2_fraction"),
        (2.1450261710738248e-05, 2.0843928754785833e-05, 2.148799898240067e-05),
    )
    _require_float_tuple(
        "superlevel collar occupancy",
        finest.get("core_superlevel_support_collar_occupancy"),
        (0.0, 0.0, 0.0),
    )
    if morphology.get("public_openai_reference_used") is not False:
        raise ValueError("target-free morphology cannot be relabeled as public-image fitting")
    if morphology.get("visual_correspondence_verified") is not False:
        raise ValueError("target-free morphology cannot establish public visual correspondence")

    seed = data.get("fresh_seed_pde_evidence")
    if not isinstance(seed, Mapping):
        raise ValueError("fresh_seed_pde_evidence must be an object")
    if (
        seed.get("pr") != FRESH_SEED_PR
        or seed.get("exact_head") != FRESH_SEED_HEAD
        or seed.get("actions_run") != FRESH_SEED_ACTIONS_RUN
    ):
        raise ValueError("quartic fresh-seed provenance drifted")
    if seed.get("evidence_class") != "pressure_free_vorticity_equation_fresh_seed_diagnostic":
        raise ValueError("fresh-seed evidence class drifted")
    if tuple(seed.get("seeds", ())) != EXPECTED_SEEDS:
        raise ValueError("fresh-seed set drifted")
    if seed.get("spatial_derivative_steps") != [0.02, 0.01, 0.005]:
        raise ValueError("fresh-seed spatial derivative ladder drifted")
    _require_close("fresh-seed dt", seed.get("time_step"), 0.0025)
    _require_close("fresh-seed nu", seed.get("nu"), 0.01)
    _require_float_tuple(
        "quartic_vs_static_finest_pct",
        seed.get("quartic_vs_static_finest_pct"),
        EXPECTED_QUARTIC_VS_STATIC_PCT,
        atol=1.0e-10,
    )
    _require_float_tuple(
        "quartic_vs_cubic_finest_pct",
        seed.get("quartic_vs_cubic_finest_pct"),
        EXPECTED_QUARTIC_VS_CUBIC_PCT,
        atol=1.0e-10,
    )
    if seed.get("quartic_vs_static_wins") != 3:
        raise ValueError("quartic-vs-static fresh-seed win count drifted")
    if seed.get("quartic_vs_cubic_wins") != 1:
        raise ValueError("quartic-vs-cubic fresh-seed win count drifted")
    _require_close("quartic-vs-static mean pct", seed.get("quartic_vs_static_mean_pct"), -0.9401)
    _require_close("quartic-vs-cubic mean pct", seed.get("quartic_vs_cubic_mean_pct"), 0.3164)
    _require_close(
        "quartic-vs-cubic projected mean pct",
        seed.get("quartic_vs_cubic_projected_mean_pct"),
        0.3535,
    )
    if seed.get("ordering") != "unresolved_seed_sensitive":
        raise ValueError("quartic-vs-cubic PDE ordering must remain seed-sensitive/unresolved")
    if seed.get("formal_pde_gate_assessed") is not False or seed.get("pde_validated") is not False:
        raise ValueError("pressure-free fresh-seed diagnostic cannot promote formal PDE validity")

    routing = data.get("routing")
    if not isinstance(routing, Mapping):
        raise ValueError("routing must be an object")
    if routing.get("visual_baseline_for_next_public_comparison") != "quartic_derivative_balanced_phi10":
        raise ValueError("quartic public-visual baseline route drifted")
    if routing.get("candidate_selection_resolved") is not False:
        raise ValueError("competing visualization candidate selection remains unresolved")
    if routing.get("quartic_preferred_over_cubic_for_late_morphology_collateral") is not True:
        raise ValueError("late-morphology routing conclusion drifted")
    if routing.get("quartic_pde_superiority_over_cubic") is not False:
        raise ValueError("seed-sensitive PDE evidence cannot establish quartic superiority")
    if routing.get("quartic_vs_cubic_pde_ordering") != "unresolved_seed_sensitive":
        raise ValueError("routing must preserve unresolved seed-sensitive PDE ordering")
    if routing.get("public_visual_comparison_still_required") is not True:
        raise ValueError("public visual comparison must remain pending")
    if routing.get("formal_pde_gate_still_required") is not True:
        raise ValueError("formal PDE gate must remain pending")
    if routing.get("acceptance_ready") is not False:
        raise ValueError("quartic candidate is not final-acceptance ready")

    siblings = data.get("unconsumed_sibling_evidence")
    if siblings != {
        "quartic_visualization_smoke_pr": 177,
        "quartic_temporal_derivative_pr": 179,
        "c2_compact_family_prs": [175, 178, 180],
        "c2_compact_candidate_selection_consumed": False,
    }:
        raise ValueError("newer unconsumed sibling-evidence inventory drifted")

    if data.get("status") != _EXPECTED_STATUS:
        raise ValueError("quartic evidence snapshot truth-boundary status drifted")
    return data


def load_snapshot(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the committed evidence snapshot."""
    snapshot_path = Path(path) if path is not None else _SNAPSHOT_PATH
    return _validate_snapshot(json.loads(snapshot_path.read_text(encoding="utf-8")))


def audit_snapshot(payload: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Cross-bind the evidence snapshot to the exact #172 replayable candidate recipe."""
    snapshot = load_snapshot() if payload is None else _validate_snapshot(payload)
    recipe = load_recipe()
    if recipe["candidate_sha256"] != snapshot["candidate_sha256"]:
        raise ValueError("evidence snapshot candidate does not match the replay recipe")
    if recipe["base_supported_sha256"] != snapshot["base_supported_sha256"]:
        raise ValueError("evidence snapshot parent does not match the replay recipe")
    recipe_status = recipe["status"]
    if recipe_status.get("velocity_export_ready") is not True:
        raise ValueError("base delivery recipe lost velocity-export readiness")
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if recipe_status.get(key) is not False:
            raise ValueError(f"base delivery recipe over-promoted {key}")

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_sha256": snapshot["candidate_sha256"],
        "velocity_export_ready": True,
        "visual_baseline_for_next_public_comparison": snapshot["routing"][
            "visual_baseline_for_next_public_comparison"
        ],
        "candidate_selection_resolved": False,
        "quartic_vs_cubic_pde_ordering": "unresolved_seed_sensitive",
        "visualization_ready": False,
        "pde_validated": False,
        "acceptance_ready": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_snapshot(), indent=2, sort_keys=True))
