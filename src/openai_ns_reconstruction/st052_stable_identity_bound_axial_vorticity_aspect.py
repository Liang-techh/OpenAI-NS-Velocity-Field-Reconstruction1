"""Fresh stable-identity-bound axial/vorticity aspect diagnostic for frozen ST052-M.

This module supersedes the historical CR-A9-104 receipt only at the evidence layer.
It evaluates the unchanged save/reloaded ST052-M callable, but binds the result to
the stable semantic candidate/callable identities admitted by CR-A9-105 instead of
the legacy rematerialization-sensitive whole-child hashes.  It does not modify the
candidate, fit public pixels, infer a public numerical target, or perform Navier--
Stokes acceptance.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from . import st052_identity_bound_morphology_execution as identity_bridge
from . import st052_stable_semantic_identity as stable_identity

SCHEMA = "st052-stable-identity-bound-axial-vorticity-aspect/v1"
TASK_ID = "CR-A9-109"
CANDIDATE_ID = stable_identity.CANDIDATE_ID
STABLE_IDENTITY_MAIN_MERGE = "0713513575178f11b60a59fc7bcd1d1ecd27b1da"
EXPECTED_STABLE_CANDIDATE_IDENTITY = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
EXPECTED_STABLE_VELOCITY_IDENTITY = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
TIME = 0.50
GRID_RESOLUTION = 25
BOX = (-2.0, 2.0)

PUBLIC_OBSERVABLE = {
    "id": "axial_stretching_trajectories",
    "kind": "qualitative_trajectory",
    "publisher": "OpenAI",
    "source_title": "On the Navier-Stokes Millennium Prize Problem",
    "source_url": "https://openai.com/index/navier-stokes-solution/",
    "source_contract_pr": 906,
    "source_contract_head": "cdfaab178cfd710f196250ffa4d73a38f3676644",
    "source_contract_blob_sha": "a222dbf88c4cf61e1dcc21e33cecefb7bb1aa4d1",
    "public_observation": "Displayed trajectories stretch or elongate along the vortex axis.",
    "numerical_target": None,
    "source_contract_admitted_by_this_increment": False,
}

METRIC_PROVENANCE = {
    "source_pr": 1069,
    "source_head": "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7",
    "source_path": "experiments/root_st052/agent7_st052m_observable_space_five_channel_jacobian.py",
    "source_blob_sha": "862f5b25554b1a6e7c038418092e0b20d9435bef",
    "metric_origin_path": "experiments/root_st052/agent7_st052m_threshold_free_morphology.py",
    "metric_origin_blob_sha": "c072712ced553eb6aef09a5789f53d5cea1b7aa3",
    "classification": "internal_protocol_reimplementation",
    "migration_scope": "full enstrophy-weighted axial RMS, radial RMS, and aspect ratio only",
    "difference": (
        "Agent-7 computes morphology-control sensitivities; this Agent-9 diagnostic "
        "evaluates only the frozen authenticated ST052 callable and binds the result "
        "to the admitted stable semantic candidate/callable identities."
    ),
}

EXTERNAL_METHOD_SCREEN = {
    "source_repo": "pyvista/pyvista",
    "source_commit": "f749a1b0a10a5a4c3c5ca3eedbc8f7860c6f9d9e",
    "license": "MIT",
    "classification": "screened_not_adopted",
    "migration_scope": "none",
    "difference": "direct NumPy callable/grid diagnostics are sufficient; no renderer dependency is needed",
    "new_external_method_migrated": False,
}

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "velocity_coefficients_changed": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "scientific_threshold_changed": False,
    "public_source_numeric_targets_used": False,
    "renderer_or_camera_used": False,
    "pixel_loss_used": False,
    "stable_semantic_identity_bound": True,
    "historical_cr_a9_104_receipt_reused": False,
    "axial_vorticity_aspect_measured": True,
    "axial_stretching_correspondence_verified": False,
    "velocity_export_ready": False,
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
    return stable_identity.canonical_sha256(value)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def tensor_trapezoid_weights(resolution: int) -> np.ndarray:
    n = int(resolution)
    if n < 3:
        raise ValueError("resolution must be at least 3")
    w = np.ones(n, dtype=float)
    w[[0, -1]] = 0.5
    return w[:, None, None] * w[None, :, None] * w[None, None, :]


def sample_velocity_grid(
    field: VelocityField,
    *,
    time: float = TIME,
    resolution: int = GRID_RESOLUTION,
    box: tuple[float, float] = BOX,
) -> tuple[np.ndarray, np.ndarray]:
    n = int(resolution)
    if n < 9 or n % 2 == 0:
        raise ValueError("resolution must be an odd integer >= 9")
    lo, hi = (float(box[0]), float(box[1]))
    if not np.isfinite([lo, hi]).all() or not lo < hi:
        raise ValueError("box must be finite and increasing")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    axis = np.linspace(lo, hi, n, dtype=float)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    values = np.asarray(field.velocity(xx.ravel(), yy.ravel(), zz.ravel(), float(time)), dtype=float)
    if values.shape != (n**3, 3) or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape-(n^3,3) data")
    return axis, values.reshape(n, n, n, 3)


def cartesian_vorticity(field: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(field, dtype=float)
    if values.ndim != 4 or values.shape[-1] != 3 or not np.isfinite(values).all():
        raise ValueError("field must be finite with shape (n,n,n,3)")
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("spacing must be finite and positive")
    u, v, w = (values[..., i] for i in range(3))
    _du_dx, du_dy, du_dz = np.gradient(u, spacing, spacing, spacing, edge_order=2)
    dv_dx, _dv_dy, dv_dz = np.gradient(v, spacing, spacing, spacing, edge_order=2)
    dw_dx, dw_dy, _dw_dz = np.gradient(w, spacing, spacing, spacing, edge_order=2)
    omega = np.stack((dw_dy - dv_dz, du_dz - dw_dx, dv_dx - du_dy), axis=-1)
    return omega, np.linalg.norm(omega, axis=-1)


def enstrophy_aspect_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict[str, float | bool]:
    mag = np.asarray(magnitude, dtype=float)
    coords = np.asarray(axis, dtype=float)
    n = coords.size
    if mag.shape != (n, n, n) or n < 3 or not np.isfinite(mag).all():
        raise ValueError("magnitude must be finite shape-(n,n,n) matching axis")
    xx, yy, zz = np.meshgrid(coords, coords, coords, indexing="ij")
    radius = np.hypot(xx, yy)
    weight = np.square(mag) * tensor_trapezoid_weights(n)
    total = float(np.sum(weight))
    if total <= np.finfo(float).tiny:
        raise ValueError("zero full-grid enstrophy")
    axial_rms = float(np.sqrt(np.sum(weight * np.square(np.abs(zz))) / total))
    radial_rms = float(np.sqrt(np.sum(weight * np.square(radius)) / total))
    if radial_rms <= np.finfo(float).tiny:
        raise ValueError("zero radial enstrophy moment")
    aspect = float(axial_rms / radial_rms)
    return {
        "full_axial_rms": axial_rms,
        "full_radial_rms": radial_rms,
        "full_aspect_ratio": aspect,
        "axially_elongated_proxy": bool(axial_rms > radial_rms),
        "full_enstrophy_trapezoid": float(total * (coords[1] - coords[0]) ** 3),
    }


def measure_field(field: VelocityField) -> dict[str, Any]:
    axis, velocity = sample_velocity_grid(field)
    spacing = float(axis[1] - axis[0])
    _omega, magnitude = cartesian_vorticity(velocity, spacing)
    metrics = enstrophy_aspect_metrics(magnitude, axis)
    return {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "spacing": spacing,
        "metrics": metrics,
        "proxy_interpretation": (
            "aspect>1 means only that this candidate's enstrophy-weighted axial extent exceeds its radial extent under this autonomous cube diagnostic; it is not an OpenAI numerical target or a correspondence verdict"
        ),
    }


def _protocol() -> dict[str, Any]:
    return {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "cartesian_curl_order": 2,
        "integration_rule": "tensor_product_trapezoid",
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
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
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected stable axial-aspect receipt schema/task")
    identity = receipt.get("stable_candidate_identity")
    materialization = receipt.get("materialization_join")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, materialization, measurement, truth)):
        raise ValueError("malformed stable axial-aspect receipt")
    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    if identity != expected_identity:
        raise ValueError("stable candidate/callable identity drift")
    if receipt.get("protocol") != _protocol():
        raise ValueError("frozen axial-aspect protocol drift")
    if receipt.get("public_observable") != PUBLIC_OBSERVABLE:
        raise ValueError("public qualitative observable drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("internal metric provenance drift")
    if receipt.get("external_method_screen") != EXTERNAL_METHOD_SCREEN:
        raise ValueError("external method screen drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    for key in ("stable_identity_receipt_sha256", "legacy_whole_candidate_identity_sha256", "materialization_evidence_sha256"):
        if not _is_sha256(materialization.get(key)):
            raise ValueError(f"invalid materialization join {key}")
    if materialization.get("included_in_measurement_binding") is not False:
        raise ValueError("materialization evidence leaked into semantic measurement binding")
    if materialization.get("stable_identity_main_merge") != STABLE_IDENTITY_MAIN_MERGE:
        raise ValueError("stable identity main-merge ancestry drift")

    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing aspect metrics")
    axial = float(metrics.get("full_axial_rms"))
    radial = float(metrics.get("full_radial_rms"))
    aspect = float(metrics.get("full_aspect_ratio"))
    if not all(np.isfinite(v) and v > 0.0 for v in (axial, radial, aspect)):
        raise ValueError("nonpositive/nonfinite aspect metric")
    if not np.isclose(aspect, axial / radial, rtol=0.0, atol=2.0e-15):
        raise ValueError("aspect ratio arithmetic drift")
    if metrics.get("axially_elongated_proxy") is not bool(axial > radial):
        raise ValueError("axial-elongation proxy drift")

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
    """Freshly measure axial aspect and bind it to the admitted stable semantic identity."""
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    stable_receipt = stable_identity.execute(constrained_root=constrained_root, bundle_dir=bundle_dir)
    stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    candidate_semantic_id = stable["candidate_semantic_identity_sha256"]
    velocity_semantic_id = stable["velocity_semantic_identity_sha256"]
    if candidate_semantic_id != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("fresh stable candidate semantic identity does not match admitted CR-A9-105 identity")
    if velocity_semantic_id != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("fresh stable callable identity does not match admitted CR-A9-105 identity")

    identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    legacy_id = str(candidate.identity_sha256)
    stable_legacy_id = stable_receipt["materialization_evidence"]["legacy_whole_candidate_identity_sha256"]
    if legacy_id != stable_legacy_id:
        raise ValueError("loaded callable materialization disagrees with stable-identity-verified bundle")

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
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_main_merge": STABLE_IDENTITY_MAIN_MERGE,
            "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
            "legacy_whole_candidate_identity_sha256": legacy_id,
            "materialization_evidence_sha256": stable_receipt["materialization_evidence_sha256"],
            "included_in_measurement_binding": False,
        },
        "protocol": _protocol(),
        "public_observable": dict(PUBLIC_OBSERVABLE),
        "metric_provenance": dict(METRIC_PROVENANCE),
        "external_method_screen": dict(EXTERNAL_METHOD_SCREEN),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "direct_contribution": (
            "freshly closes the ST052 axial-vorticity observable measurement against the admitted stable semantic candidate/callable identity without reusing the historical CR-A9-104 receipt"
        ),
        "remaining_limits": [
            "the official public source supplies no numerical axial-aspect target",
            "aspect>1 is a candidate-only proxy and does not establish visual/source correspondence",
            "no renderer/camera registration or pixel matching is used",
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
    args.report.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "candidate_semantic_identity_sha256": receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"],
                "velocity_semantic_identity_sha256": receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"],
                "legacy_whole_candidate_identity_sha256": receipt["materialization_join"]["legacy_whole_candidate_identity_sha256"],
                "full_aspect_ratio": receipt["measurement"]["metrics"]["full_aspect_ratio"],
                "axially_elongated_proxy": receipt["measurement"]["metrics"]["axially_elongated_proxy"],
                "measurement_sha256": receipt["measurement_sha256"],
                "visualization_ready": receipt["truth_boundary"]["visualization_ready"],
                "pde_validated": receipt["truth_boundary"]["pde_validated"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
