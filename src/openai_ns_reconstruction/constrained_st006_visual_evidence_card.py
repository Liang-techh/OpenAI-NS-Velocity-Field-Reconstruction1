"""Truth-bounded candidate-side visual evidence card for retained ST006.

This module binds already-produced, candidate-specific receipts into one fixed
public-observable checklist. It does not score resemblance, fit hidden OpenAI
quantities, or provide PDE acceptance evidence.
"""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from math import isfinite
from typing import Any, Mapping


SCHEMA = "st006_public_visual_evidence_card_v1"
TASK_ID = "CR-A9-041"
BASE_COMMIT = "f0193d66c9d92948b4820ebcb70263673995b324"
CANDIDATE = "ST006"
CANDIDATE_SHA256 = "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3"
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"
PUBLIC_SOURCE_DATE = "2026-09-08"

OBSERVABLE_IDS = (
    "inward_spiraling",
    "axial_stretching",
    "vortex_swirl",
    "increasing_elongation",
    "central_region_shrinks",
    "speed_increases",
    "relative_angular_rotation_color_encoding",
    "circulating_speed_depends_on_radius",
)

SOURCE_BINDINGS = {
    "radial_circulation": {
        "task_id": "CR-A9-038",
        "pr": 334,
        "head_sha": "0a9cfca629c0fc992fa1496dd30e5916608347d7",
        "workflow_run": 35313443818,
        "artifact_id": 10534586681,
        "artifact_digest": "sha256:fcf5a999663e05f1a5daa31d5762c1416907d29cb4d511d4d308f69f63f0a965",
        "report_sha256": "cb07626d9196e087a6728473f0b844d323f743c071bafa5743c28a6782833553",
        "schema": "st006_centerplane_radial_circulation_v1",
    },
    "vorticity_morphology": {
        "task_id": "CR-A9-039",
        "pr": 341,
        "head_sha": "93405a1cb73b6283e09d89795457a6f67362a8b9",
        "workflow_run": 35317812917,
        "artifact_id": 10535552821,
        "artifact_digest": "sha256:928c179151a283fecee886a18cb6b673952c9b7ad1734a25e47b4da378292e92",
        "report_sha256": "b2e15394dfd0fbc7b6d4cfc92a29225553c74f3d36dd06ca58b42dc42ebf3985",
        "schema": "st006_vorticity_morphology_timeseries_v1",
    },
    "material_paths": {
        "task_id": "CR-A9-040",
        "pr": 349,
        "head_sha": "3911718884f021f0c0e3ba41ed694b9688764084",
        "workflow_run": 35322722094,
        "artifact_id": 10536949231,
        "artifact_digest": "sha256:e73165cbc4ae75e2d9b437ed1b8606a4af9c3e7439ffb50d65c315e33378c88d",
        "report_sha256": "2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772",
        "schema": "st006_material_path_observables_v1",
    },
}


def _canonical_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _with_digest(payload: dict[str, Any]) -> dict[str, Any]:
    body = deepcopy(payload)
    body.pop("card_sha256", None)
    payload = deepcopy(body)
    payload["card_sha256"] = sha256(_canonical_bytes(body)).hexdigest()
    return payload


def make_st006_visual_evidence_card() -> dict[str, Any]:
    """Return the frozen ST006 candidate-side public-observable evidence card."""

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_commit": BASE_COMMIT,
        "candidate": CANDIDATE,
        "candidate_sha256": CANDIDATE_SHA256,
        "public_reference": {
            "publisher": "OpenAI",
            "url": PUBLIC_SOURCE_URL,
            "published": PUBLIC_SOURCE_DATE,
            "classification": "direct_public_observable_use_only",
            "numerical_velocity_samples_used": False,
            "hidden_frame_time_used": False,
            "camera_registration_used": False,
            "streamline_seed_recovery_used": False,
            "quantitative_visual_target_used": False,
        },
        "evidence_sources": deepcopy(SOURCE_BINDINGS),
        "observables": {
            "inward_spiraling": {
                "status": "measured_candidate_side_only",
                "source": "material_paths",
                "measurements": {
                    "path_count": 48,
                    "inward_path_count": 48,
                    "mean_radius_change": -0.08601065505093726,
                    "mean_absolute_turns": 0.024485286231482325,
                    "minimum_absolute_turns": 0.003880821959757115,
                    "maximum_absolute_turns": 0.0436185396886357,
                },
                "interpretation": "inward transport is ensemble-wide under the frozen seed contract; angular travel is nonzero but modest",
            },
            "axial_stretching": {
                "status": "measured_candidate_side_only",
                "source": "material_paths",
                "measurements": {
                    "paired_material_line_count": 24,
                    "pair_growth_count": 16,
                    "pair_shrink_count": 8,
                    "mean_pair_separation_change": 0.07130825328092755,
                    "mean_pair_separation_ratio": 1.1188470888015458,
                },
                "interpretation": "axial material-line separation grows on average but is spatially heterogeneous",
            },
            "vortex_swirl": {
                "status": "measured_candidate_side_only",
                "source": "radial_circulation",
                "measurements": {
                    "peak_abs_circulating_speed_t025": 0.11911048491205961,
                    "peak_abs_circulating_speed_t075": 0.19942107771770323,
                    "peak_abs_circulating_speed_delta": 0.08031059280564362,
                },
                "interpretation": "registered center-plane rings carry a nonzero coherent circulating component",
            },
            "increasing_elongation": {
                "status": "measured_candidate_side_only",
                "source": "vorticity_morphology",
                "measurements": {
                    "vorticity_aspect_ratio_final_over_initial": 1.1215380885092212,
                    "principal_extent_final_over_initial": 0.9911917886513485,
                    "transverse_extent_final_over_initial": 0.8837789806754289,
                },
                "interpretation": "vorticity aspect increases mainly through transverse contraction, not measured axial-extent growth",
            },
            "central_region_shrinks": {
                "status": "measured_candidate_side_only",
                "source": "vorticity_morphology+radial_circulation",
                "measurements": {
                    "transverse_extent_final_over_initial": 0.8837789806754289,
                    "covariance_volume_proxy_final_over_initial": 0.774185498561495,
                    "swirl_weighted_radius_delta": -0.036714143473095096,
                },
                "interpretation": "multiple registered-frame concentration measures contract over the finite window",
            },
            "speed_increases": {
                "status": "measured_candidate_side_only",
                "source": "vorticity_morphology+radial_circulation+material_paths",
                "measurements": {
                    "velocity_rms_final_over_initial": 1.0728912572743456,
                    "peak_centerplane_mean_speed_delta": 0.0733587158961791,
                    "mean_material_path_speed_change": 0.07706346903942406,
                },
                "interpretation": "three candidate-side speed summaries increase over their frozen contracts",
            },
            "relative_angular_rotation_color_encoding": {
                "status": "not_measured_renderer_semantics",
                "source": None,
                "measurements": {},
                "interpretation": "the OpenAI color legend is a rendering semantic; no candidate renderer/color mapping is inferred here",
            },
            "circulating_speed_depends_on_radius": {
                "status": "measured_candidate_side_only",
                "source": "radial_circulation",
                "measurements": {
                    "relative_radial_range_t025": 0.9879188220042697,
                    "relative_radial_range_t075": 0.9645746143371134,
                    "radial_sample_count": 31,
                },
                "interpretation": "ring-mean absolute circulating speed has strong radial variation at fixed registered times",
            },
        },
        "routing_summary": {
            "candidate_strengths": [
                "ensemble-wide inward radial material transport",
                "nonzero radius-dependent circulation",
                "transverse concentration and increasing vorticity aspect",
                "candidate-side speed increase",
            ],
            "candidate_limits": [
                "material paths rotate only a small fraction of a turn over t=0.25..0.75",
                "axial material-line stretching is positive on average but not uniform",
                "vorticity principal extent does not itself grow over the measured window",
                "renderer color semantics remain unmeasured",
            ],
            "next_candidate_contract": "rerun the same three frozen evidence contracts before comparing any promoted candidate; do not tune seeds, times, camera, or thresholds after seeing the result",
            "visual_score_defined": False,
            "candidate_selected_by_this_card": False,
        },
        "truth_boundary": {
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "visual_acceptance_threshold_defined": False,
            "public_observable_coverage_is_acceptance": False,
            "used_for_pde_acceptance": False,
            "free_residual_force_used": False,
            "collapsed_velocity_accepted": False,
        },
    }
    return _with_digest(payload)


def _require_finite_numbers(value: Any) -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not isfinite(float(value)):
            raise ValueError("nonfinite numeric value in evidence card")
        return
    if isinstance(value, Mapping):
        for item in value.values():
            _require_finite_numbers(item)
        return
    if isinstance(value, (list, tuple)):
        for item in value:
            _require_finite_numbers(item)
        return
    raise TypeError(f"unsupported evidence-card value type: {type(value).__name__}")


def audit_st006_visual_evidence_card(card: Mapping[str, Any]) -> None:
    """Fail closed if the ST006 evidence card drifts or launders truth claims."""

    _require_finite_numbers(card)
    if card.get("schema") != SCHEMA or card.get("task_id") != TASK_ID:
        raise ValueError("unexpected evidence-card schema/task")
    if card.get("base_commit") != BASE_COMMIT:
        raise ValueError("evidence card is not bound to the audited main base")
    if card.get("candidate") != CANDIDATE or card.get("candidate_sha256") != CANDIDATE_SHA256:
        raise ValueError("candidate identity drift")
    if card.get("evidence_sources") != SOURCE_BINDINGS:
        raise ValueError("upstream evidence provenance drift")

    reference = card.get("public_reference")
    if not isinstance(reference, Mapping) or reference.get("url") != PUBLIC_SOURCE_URL:
        raise ValueError("public-reference provenance drift")
    for key in (
        "numerical_velocity_samples_used",
        "hidden_frame_time_used",
        "camera_registration_used",
        "streamline_seed_recovery_used",
        "quantitative_visual_target_used",
    ):
        if reference.get(key) is not False:
            raise ValueError(f"forbidden public-reference promotion: {key}")

    observables = card.get("observables")
    if not isinstance(observables, Mapping) or tuple(observables.keys()) != OBSERVABLE_IDS:
        raise ValueError("public observable set/order drift")
    allowed_status = {"measured_candidate_side_only", "not_measured_renderer_semantics"}
    for observable_id, entry in observables.items():
        if not isinstance(entry, Mapping) or entry.get("status") not in allowed_status:
            raise ValueError(f"invalid status for {observable_id}")
        if "target" in entry or "threshold" in entry:
            raise ValueError("numerical target/threshold injection is forbidden")

    inward = observables["inward_spiraling"]["measurements"]
    if inward["path_count"] != 48 or inward["inward_path_count"] != 48:
        raise ValueError("material-path population drift")
    if not inward["mean_radius_change"] < 0 or not inward["mean_absolute_turns"] > 0:
        raise ValueError("material-path measurement sign drift")

    axial = observables["axial_stretching"]["measurements"]
    if axial["paired_material_line_count"] != axial["pair_growth_count"] + axial["pair_shrink_count"]:
        raise ValueError("paired material-line accounting drift")
    if axial["mean_pair_separation_ratio"] <= 0:
        raise ValueError("invalid material-line separation ratio")

    elongation = observables["increasing_elongation"]["measurements"]
    if elongation["vorticity_aspect_ratio_final_over_initial"] <= 1:
        raise ValueError("ST006 aspect-trend receipt drift")
    if elongation["transverse_extent_final_over_initial"] >= 1:
        raise ValueError("ST006 transverse-contraction receipt drift")

    shrink = observables["central_region_shrinks"]["measurements"]
    if shrink["covariance_volume_proxy_final_over_initial"] >= 1 or shrink["swirl_weighted_radius_delta"] >= 0:
        raise ValueError("ST006 contraction receipt drift")

    speed = observables["speed_increases"]["measurements"]
    if speed["velocity_rms_final_over_initial"] <= 1:
        raise ValueError("ST006 speed-trend receipt drift")

    color = observables["relative_angular_rotation_color_encoding"]
    if color["status"] != "not_measured_renderer_semantics" or color["measurements"]:
        raise ValueError("renderer color semantics must stay unmeasured")

    routing = card.get("routing_summary")
    if not isinstance(routing, Mapping) or routing.get("visual_score_defined") is not False:
        raise ValueError("visual score injection is forbidden")
    if routing.get("candidate_selected_by_this_card") is not False:
        raise ValueError("candidate selection cannot be promoted by this card")

    truth = card.get("truth_boundary")
    required_false = (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "visual_acceptance_threshold_defined",
        "public_observable_coverage_is_acceptance",
        "used_for_pde_acceptance",
        "free_residual_force_used",
        "collapsed_velocity_accepted",
    )
    if not isinstance(truth, Mapping) or any(truth.get(key) is not False for key in required_false):
        raise ValueError("truth-boundary laundering detected")

    digest = card.get("card_sha256")
    if not isinstance(digest, str) or len(digest) != 64:
        raise ValueError("missing/invalid card digest")
    body = deepcopy(dict(card))
    body.pop("card_sha256", None)
    if sha256(_canonical_bytes(body)).hexdigest() != digest:
        raise ValueError("evidence-card digest mismatch")
