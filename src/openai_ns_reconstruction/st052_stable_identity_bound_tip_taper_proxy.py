"""Stable-identity-bound tip-vs-core taper diagnostic for frozen ST052-M.

Preregistered in issue #1194.  This is a Constrained Agent-7 basis-routing
increment, not a candidate mutation.  It reuses the admitted stable ST052-M
candidate/callable identity and the already-frozen renderer-independent
vorticity morphology geometry to answer one bounded question: are both smooth
vorticity-tip regions radially narrower than the central vorticity band?

The comparison boundary ``ratio < 1`` is an autonomous candidate-side geometry
proxy only.  It is not an OpenAI numerical target and cannot establish visual
correspondence, PDE validity, paper exactness, or exact OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from . import st052_stable_identity_bound_axial_vorticity_aspect as axial

SCHEMA = "st052-stable-identity-bound-tip-taper-proxy/v1"
TASK_ID = "CR003-ST052M-STABLE-TIP-TAPER-PROXY-141"
PREREG_ISSUE = 1194
SOURCE_AXIAL_EVIDENCE_PR = 1163
SOURCE_AXIAL_EVIDENCE_HEAD = "8fec2517677adb37754df8c0164330d6f346db57"
SOURCE_AXIAL_EVIDENCE_MAIN_MERGE = "b0093e04039946923de10e421e001cd363213dcc"
SOURCE_MORPHOLOGY_PR = 1069
SOURCE_MORPHOLOGY_HEAD = "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7"
SOURCE_MORPHOLOGY_BLOB_SHA = "c072712ced553eb6aef09a5789f53d5cea1b7aa3"

CANDIDATE_ID = axial.CANDIDATE_ID
EXPECTED_STABLE_CANDIDATE_IDENTITY = axial.EXPECTED_STABLE_CANDIDATE_IDENTITY
EXPECTED_STABLE_VELOCITY_IDENTITY = axial.EXPECTED_STABLE_VELOCITY_IDENTITY
TIME = axial.TIME
GRID_RESOLUTION = axial.GRID_RESOLUTION
BOX = axial.BOX
SUPPORT_HALF_HEIGHT = 2.0
TIP_WEIGHT_WINDOW = (0.50, 0.80)
CENTRAL_BAND_MAX = 0.35

METRIC_PROVENANCE = {
    "stable_execution_source_pr": SOURCE_AXIAL_EVIDENCE_PR,
    "stable_execution_source_head": SOURCE_AXIAL_EVIDENCE_HEAD,
    "stable_execution_main_merge": SOURCE_AXIAL_EVIDENCE_MAIN_MERGE,
    "morphology_source_pr": SOURCE_MORPHOLOGY_PR,
    "morphology_source_head": SOURCE_MORPHOLOGY_HEAD,
    "morphology_metric_origin_path": "experiments/root_st052/agent7_st052m_threshold_free_morphology.py",
    "morphology_metric_origin_blob_sha": SOURCE_MORPHOLOGY_BLOB_SHA,
    "classification": "internal_target_free_protocol_reuse",
    "difference": (
        "the prior morphology asset pooled both tips; this increment splits upper/lower "
        "smooth-tip radial RMS and compares the worse tip with the frozen central band "
        "on the admitted stable ST052 callable identity"
    ),
}

PUBLIC_CONTEXT = {
    "qualitative_relation": "tip taper is an autonomous streamline/vorticity morphology diagnostic",
    "public_numerical_target": None,
    "source_correspondence_inferred": False,
}

TRUTH_BOUNDARY = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "renderer_or_camera_fit_used": False,
    "pixel_similarity_objective_used": False,
    "stable_semantic_identity_bound": True,
    "tip_taper_proxy_measured": True,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


class VelocityField(Protocol):
    def velocity(self, x: Any, y: Any, z: Any, t: float) -> np.ndarray:
        """Return Cartesian velocity with final component axis of length three."""


def _canonical_sha256(value: Any) -> str:
    return axial._canonical_sha256(value)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _assert_upstream_lock() -> None:
    if axial.TASK_ID != "CR-A9-109":
        raise RuntimeError("stable axial execution source task drifted")
    if axial.CANDIDATE_ID != CANDIDATE_ID:
        raise RuntimeError("stable candidate id drifted")
    if axial.EXPECTED_STABLE_CANDIDATE_IDENTITY != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise RuntimeError("stable candidate semantic identity drifted")
    if axial.EXPECTED_STABLE_VELOCITY_IDENTITY != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise RuntimeError("stable callable semantic identity drifted")
    if axial.TIME != 0.50 or axial.GRID_RESOLUTION != 25 or tuple(axial.BOX) != (-2.0, 2.0):
        raise RuntimeError("stable morphology sampling protocol drifted")


def smooth_tip_weight(normalized_abs_z: np.ndarray) -> np.ndarray:
    """Frozen smooth tip bump from the threshold-free morphology protocol."""
    s = np.asarray(normalized_abs_z, dtype=float)
    if not np.isfinite(s).all():
        raise ValueError("normalized_abs_z must be finite")
    a, b = TIP_WEIGHT_WINDOW
    out = np.zeros_like(s)
    mask = (s > a) & (s < b)
    xi = (s[mask] - a) / (b - a)
    out[mask] = np.exp(4.0 - 1.0 / (xi * (1.0 - xi)))
    return out


def _weighted_rms(values: np.ndarray, weight: np.ndarray, *, label: str) -> float:
    v = np.asarray(values, dtype=float)
    w = np.asarray(weight, dtype=float)
    if v.shape != w.shape or not np.isfinite(v).all() or not np.isfinite(w).all():
        raise ValueError(f"invalid {label} weighted-RMS arrays")
    if np.any(w < 0.0):
        raise ValueError(f"negative {label} weights")
    total = float(np.sum(w))
    if total <= np.finfo(float).tiny:
        raise ValueError(f"empty {label} weighted region")
    return float(np.sqrt(np.sum(w * np.square(v)) / total))


def tip_taper_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict[str, float | bool]:
    """Measure upper/lower smooth-tip radial RMS against the central radial RMS."""
    mag = np.asarray(magnitude, dtype=float)
    coords = np.asarray(axis, dtype=float)
    n = int(coords.size)
    if n < 3 or mag.shape != (n, n, n):
        raise ValueError("magnitude must match a one-dimensional axis")
    if not np.isfinite(mag).all() or np.any(mag < 0.0) or not np.isfinite(coords).all():
        raise ValueError("magnitude and axis must be finite, with nonnegative magnitude")
    if not np.all(np.diff(coords) > 0.0):
        raise ValueError("axis must be strictly increasing")

    xx, yy, zz = np.meshgrid(coords, coords, coords, indexing="ij")
    radius = np.hypot(xx, yy)
    normalized_abs_z = np.abs(zz) / SUPPORT_HALF_HEIGHT
    enstrophy = np.square(mag) * axial.tensor_trapezoid_weights(n)
    if float(np.sum(enstrophy)) <= np.finfo(float).tiny:
        raise ValueError("zero full-grid enstrophy")

    tip_shape = smooth_tip_weight(normalized_abs_z)
    upper_weight = enstrophy * tip_shape * (zz > 0.0)
    lower_weight = enstrophy * tip_shape * (zz < 0.0)
    central_weight = enstrophy * (normalized_abs_z <= CENTRAL_BAND_MAX)

    upper = _weighted_rms(radius, upper_weight, label="upper-tip")
    lower = _weighted_rms(radius, lower_weight, label="lower-tip")
    central = _weighted_rms(radius, central_weight, label="central-band")
    if central <= np.finfo(float).tiny:
        raise ValueError("zero central radial RMS")

    worst = float(max(upper, lower))
    ratio = float(worst / central)
    mismatch = float(abs(upper - lower) / max(upper, lower))
    return {
        "tip_radial_rms_upper": upper,
        "tip_radial_rms_lower": lower,
        "central_radial_rms": central,
        "worst_tip_radial_rms": worst,
        "worst_tip_to_central_ratio": ratio,
        "upper_lower_tip_relative_mismatch": mismatch,
        "both_tips_radially_tapered_proxy": bool(worst < central),
    }


def _protocol() -> dict[str, Any]:
    return {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "cartesian_curl_order": 2,
        "integration_rule": "tensor_product_trapezoid",
        "support_half_height": SUPPORT_HALF_HEIGHT,
        "tip_weight_window_abs_z_over_2": list(TIP_WEIGHT_WINDOW),
        "tip_weight_formula": "exp(4-1/(xi*(1-xi))) inside fixed window; zero outside",
        "central_band_abs_z_over_2_max": CENTRAL_BAND_MAX,
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
    }


def measure_field(field: VelocityField) -> dict[str, Any]:
    _assert_upstream_lock()
    axis, velocity = axial.sample_velocity_grid(
        field, time=TIME, resolution=GRID_RESOLUTION, box=BOX
    )
    spacing = float(axis[1] - axis[0])
    _omega, magnitude = axial.cartesian_vorticity(velocity, spacing)
    return {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "spacing": spacing,
        "metrics": tip_taper_metrics(magnitude, axis),
        "proxy_interpretation": (
            "ratio<1 means only that both frozen candidate smooth-tip radial RMS values "
            "are narrower than the central radial RMS under this autonomous diagnostic; "
            "it is not a public numerical target or visual/source correspondence verdict"
        ),
    }


def routing_from_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    tapered = bool(metrics["both_tips_radially_tapered_proxy"])
    if tapered:
        return {
            "classification": "no_candidate_side_blunt_tip_deficit_established",
            "tip_specific_basis_growth_authorized": False,
            "axial_turnover_change_authorized": False,
            "next_minimal_action": "require_distinct_stable_identity_bound_trajectory_or_return_flow_discrepancy",
        }
    return {
        "classification": "candidate_side_tip_taper_deficit_established",
        "tip_specific_basis_growth_authorized": False,
        "axial_turnover_change_authorized": False,
        "next_minimal_action": "route_to_existing_tip_thickness_control_sensitivity_before_basis_growth",
    }


def _measurement_binding_sha256(identity: dict[str, Any], measurement_sha256: str) -> str:
    return _canonical_sha256(
        {
            "schema": SCHEMA,
            "stable_candidate_identity": identity,
            "protocol": _protocol(),
            "measurement_sha256": measurement_sha256,
            "metric_provenance": METRIC_PROVENANCE,
        }
    )


def validate_receipt(receipt: dict[str, Any]) -> None:
    _assert_upstream_lock()
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected stable tip-taper receipt schema/task")
    if receipt.get("prereg_issue") != PREREG_ISSUE:
        raise ValueError("preregistration identity drift")

    identity = receipt.get("stable_candidate_identity")
    materialization = receipt.get("materialization_join")
    measurement = receipt.get("measurement")
    routing = receipt.get("routing")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, materialization, measurement, routing, truth)):
        raise ValueError("malformed stable tip-taper receipt")

    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    if identity != expected_identity:
        raise ValueError("stable candidate/callable identity drift")
    if receipt.get("protocol") != _protocol():
        raise ValueError("frozen tip-taper protocol drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("metric provenance drift")
    if receipt.get("public_context") != PUBLIC_CONTEXT:
        raise ValueError("public-context boundary drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")

    for key in (
        "stable_identity_receipt_sha256",
        "legacy_whole_candidate_identity_sha256",
        "materialization_evidence_sha256",
    ):
        if not _is_sha256(materialization.get(key)):
            raise ValueError(f"invalid materialization join {key}")
    if materialization.get("included_in_measurement_binding") is not False:
        raise ValueError("materialization evidence leaked into semantic measurement binding")

    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing tip-taper metrics")
    required = (
        "tip_radial_rms_upper",
        "tip_radial_rms_lower",
        "central_radial_rms",
        "worst_tip_radial_rms",
        "worst_tip_to_central_ratio",
        "upper_lower_tip_relative_mismatch",
    )
    values = {key: float(metrics.get(key)) for key in required}
    if not all(np.isfinite(v) and v >= 0.0 for v in values.values()):
        raise ValueError("nonfinite/negative tip-taper metric")
    upper = values["tip_radial_rms_upper"]
    lower = values["tip_radial_rms_lower"]
    central = values["central_radial_rms"]
    worst = values["worst_tip_radial_rms"]
    ratio = values["worst_tip_to_central_ratio"]
    mismatch = values["upper_lower_tip_relative_mismatch"]
    if min(upper, lower, central, worst) <= np.finfo(float).tiny:
        raise ValueError("degenerate tip-taper metric")
    if not np.isclose(worst, max(upper, lower), rtol=0.0, atol=2.0e-15):
        raise ValueError("worst-tip arithmetic drift")
    if not np.isclose(ratio, worst / central, rtol=0.0, atol=2.0e-15):
        raise ValueError("tip/central ratio arithmetic drift")
    expected_mismatch = abs(upper - lower) / max(upper, lower)
    if not np.isclose(mismatch, expected_mismatch, rtol=0.0, atol=2.0e-15):
        raise ValueError("tip-symmetry arithmetic drift")
    tapered = bool(worst < central)
    if metrics.get("both_tips_radially_tapered_proxy") is not tapered:
        raise ValueError("tip-taper proxy drift")
    if routing != routing_from_metrics(metrics):
        raise ValueError("tip-taper routing drift")
    if routing.get("tip_specific_basis_growth_authorized") is not False:
        raise ValueError("tip basis growth cannot be authorized by this receipt")

    measurement_sha = _canonical_sha256(measurement)
    if receipt.get("measurement_sha256") != measurement_sha:
        raise ValueError("measurement checksum mismatch")
    if receipt.get("candidate_measurement_binding_sha256") != _measurement_binding_sha256(identity, measurement_sha):
        raise ValueError("stable candidate/measurement binding mismatch")
    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("receipt checksum mismatch")


def execute(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    """Measure tip taper on the same stable semantic ST052 candidate identity."""
    _assert_upstream_lock()
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    stable_receipt = axial.stable_identity.execute(
        constrained_root=constrained_root, bundle_dir=bundle_dir
    )
    axial.stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    candidate_semantic_id = stable["candidate_semantic_identity_sha256"]
    velocity_semantic_id = stable["velocity_semantic_identity_sha256"]
    if candidate_semantic_id != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("fresh stable candidate semantic identity mismatch")
    if velocity_semantic_id != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("fresh stable callable semantic identity mismatch")

    axial.identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    legacy_id = str(candidate.identity_sha256)
    stable_legacy_id = stable_receipt["materialization_evidence"][
        "legacy_whole_candidate_identity_sha256"
    ]
    if legacy_id != stable_legacy_id:
        raise ValueError("loaded callable materialization disagrees with stable-identity bundle")

    identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": candidate_semantic_id,
        "velocity_semantic_identity_sha256": velocity_semantic_id,
    }
    measurement = measure_field(candidate)
    measurement_sha = _canonical_sha256(measurement)
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
            "legacy_whole_candidate_identity_sha256": legacy_id,
            "materialization_evidence_sha256": stable_receipt["materialization_evidence_sha256"],
            "included_in_measurement_binding": False,
        },
        "protocol": _protocol(),
        "metric_provenance": dict(METRIC_PROVENANCE),
        "public_context": dict(PUBLIC_CONTEXT),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "routing": routing_from_metrics(measurement["metrics"]),
        "direct_contribution": (
            "closes one stable candidate-side blunt-tip basis-growth trigger without changing velocity; "
            "a negative taper proxy routes first to existing tip-thickness controls, while a positive "
            "proxy requires a different trajectory/return-flow discrepancy"
        ),
        "remaining_limits": [
            "the public source supplies no numerical tip-to-core taper target",
            "this candidate-side proxy does not establish visual/source correspondence",
            "no coefficient or new basis is selected",
            "no pressure, restricted forcing, or complete held-out NS residual is evaluated",
        ],
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = _canonical_sha256(receipt)
    validate_receipt(receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constrained-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    receipt = execute(
        constrained_root=args.constrained_root,
        source_root=args.source_root,
        bundle_dir=args.bundle_dir,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    m = receipt["measurement"]["metrics"]
    print(
        json.dumps(
            {
                "candidate_semantic_identity_sha256": receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"],
                "velocity_semantic_identity_sha256": receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"],
                "tip_radial_rms_upper": m["tip_radial_rms_upper"],
                "tip_radial_rms_lower": m["tip_radial_rms_lower"],
                "central_radial_rms": m["central_radial_rms"],
                "worst_tip_to_central_ratio": m["worst_tip_to_central_ratio"],
                "both_tips_radially_tapered_proxy": m["both_tips_radially_tapered_proxy"],
                "routing": receipt["routing"]["classification"],
                "visualization_ready": receipt["truth_boundary"]["visualization_ready"],
                "pde_validated": receipt["truth_boundary"]["pde_validated"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
