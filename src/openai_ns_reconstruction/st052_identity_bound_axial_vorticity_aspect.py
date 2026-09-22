"""Identity-bound axial/vorticity aspect diagnostic for frozen ST052-M.

This module fills one narrowly documented visualization-evidence gap: the merged
Agent-9 cylindrical morphology receipt measures swirl/inward/core trends but not
an explicit axial-vorticity aspect observable.  The diagnostic below evaluates
only the already-frozen, save/reloaded ST052-M ``velocity(x,y,z,t)`` field.
It does not modify the candidate, fit OpenAI pixels, infer a source numerical
target, or perform Navier--Stokes acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from . import st052_identity_bound_morphology_execution as identity_bridge

SCHEMA = "st052-identity-bound-axial-vorticity-aspect/v1"
TASK_ID = "CR-A9-104"
CANDIDATE_ID = identity_bridge.CANDIDATE_ID
TIME = 0.50
GRID_RESOLUTION = 25
BOX = (-2.0, 2.0)

# This is a qualitative source label only.  PR #906 remains open, so this
# diagnostic records the preregistered source-contract provenance without
# treating that PR as admitted numerical evidence.
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

# Match the already-preregistered Agent-7 coordinate definition, but do not
# import the open Agent-7 branch or its control/Jacobian machinery.  We
# independently recode only the target-free vorticity/enstrophy moment metric.
METRIC_PROVENANCE = {
    "source_pr": 1069,
    "source_head": "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7",
    "source_path": "experiments/root_st052/agent7_st052m_observable_space_five_channel_jacobian.py",
    "source_blob_sha": "862f5b25554b1a6e7c038418092e0b20d9435bef",
    "metric_origin_path": "experiments/root_st052/agent7_st052m_threshold_free_morphology.py",
    "metric_origin_blob_sha": "c072712ced553eb6aef09a5789f53d5cea1b7aa3",
    "classification": "internal_protocol_reimplementation",
    "migration_scope": "full enstrophy-weighted axial RMS, radial RMS, and their aspect ratio only",
    "difference": (
        "Agent-7 computes morphology-control sensitivities; this Agent-9 module "
        "evaluates the frozen authenticated ST052-M velocity directly and binds "
        "the measurement to the candidate/velocity identity."
    ),
}

# No new third-party method is adopted in this increment.  Preserve the same
# external-screen conclusion already used by the merged morphology engine.
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


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def tensor_trapezoid_weights(resolution: int) -> np.ndarray:
    """Product trapezoid weights; the common h^3 factor cancels in moment ratios."""
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
    """Sample only the public Cartesian velocity callable on a frozen cube."""
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
    values = np.asarray(
        field.velocity(xx.ravel(), yy.ravel(), zz.ravel(), float(time)), dtype=float
    )
    if values.shape != (n**3, 3) or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape-(n^3,3) data")
    return axis, values.reshape(n, n, n, 3)


def cartesian_vorticity(field: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Independent centered-grid Cartesian curl with second-order boundary stencils."""
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
    magnitude = np.linalg.norm(omega, axis=-1)
    return omega, magnitude


def enstrophy_aspect_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict[str, float | bool]:
    """Target-free full-grid enstrophy moments compatible with Agent-7's coordinate."""
    mag = np.asarray(magnitude, dtype=float)
    coords = np.asarray(axis, dtype=float)
    n = coords.size
    if mag.shape != (n, n, n) or n < 3 or not np.isfinite(mag).all():
        raise ValueError("magnitude must be finite shape-(n,n,n) matching axis")

    xx, yy, zz = np.meshgrid(coords, coords, coords, indexing="ij")
    radius = np.hypot(xx, yy)
    abs_z = np.abs(zz)
    weight = np.square(mag) * tensor_trapezoid_weights(n)
    total = float(np.sum(weight))
    if total <= np.finfo(float).tiny:
        raise ValueError("zero full-grid enstrophy")

    axial_rms = float(np.sqrt(np.sum(weight * np.square(abs_z)) / total))
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
            "aspect>1 means the candidate's enstrophy-weighted axial extent exceeds its radial extent "
            "under this autonomous cube diagnostic. It is not an OpenAI numerical target or a source-correspondence verdict."
        ),
    }


def _measurement_binding_sha256(identity: dict[str, Any], measurement_sha256: str) -> str:
    return _canonical_sha256(
        {
            "schema": SCHEMA,
            "candidate_identity": identity,
            "protocol": {
                "time": TIME,
                "grid_resolution": GRID_RESOLUTION,
                "box": list(BOX),
                "cartesian_curl_order": 2,
                "integration_rule": "tensor_product_trapezoid",
                "source_numeric_targets_used": False,
                "renderer_or_camera_used": False,
                "pixel_loss_used": False,
            },
            "measurement_sha256": measurement_sha256,
            "metric_provenance": METRIC_PROVENANCE,
        }
    )


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected axial-aspect receipt schema/task")
    identity = receipt.get("candidate_identity")
    protocol = receipt.get("protocol")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, protocol, measurement, truth)):
        raise ValueError("malformed axial-aspect receipt")
    if identity.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("candidate identity drift")
    for key in ("candidate_sha256", "velocity_identity_sha256"):
        value = identity.get(key)
        if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError(f"invalid {key}")
    expected_protocol = {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "cartesian_curl_order": 2,
        "integration_rule": "tensor_product_trapezoid",
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
    }
    if protocol != expected_protocol:
        raise ValueError("frozen axial-aspect protocol drift")
    if receipt.get("public_observable") != PUBLIC_OBSERVABLE:
        raise ValueError("public qualitative observable drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("internal metric provenance drift")
    if receipt.get("external_method_screen") != EXTERNAL_METHOD_SCREEN:
        raise ValueError("external method screen drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")

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
    binding = _measurement_binding_sha256(identity, measurement_sha)
    if receipt.get("candidate_measurement_binding_sha256") != binding:
        raise ValueError("candidate/measurement binding mismatch")
    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("receipt checksum mismatch")


def execute(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    """Evaluate the frozen ST052 aspect metric on the authenticated save/load identity."""
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    identity_bridge._install_constrained_package_path(constrained_root)
    acceptance = importlib.import_module(
        "openai_ns_reconstruction.constrained_st052_runtime_dependency_acceptance"
    )
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    acceptance_receipt = acceptance.audit(constrained_root)

    contract_path = constrained_root / identity_bridge.DEPENDENCY_CONTRACT_REL
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract_sha = identity_bridge._sha256_file(contract_path)
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    velocity_payload = identity_bridge.build_velocity_identity_payload(
        candidate.manifest, contract, dependency_contract_sha256=contract_sha
    )
    candidate_sha = str(candidate.identity_sha256)
    velocity_sha = identity_bridge._canonical_sha256(velocity_payload)
    identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_sha256": candidate_sha,
        "velocity_identity_sha256": velocity_sha,
    }

    measurement = measure_field(candidate)
    measurement_sha = _canonical_sha256(measurement)
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_identity": identity,
        "identity_ancestry": {
            "merged_a9_morphology_execution_commit": "bcd0982d13bca63ed9408b8224a649377edc5a59",
            "constrained_runtime_acceptance_merge": identity_bridge.CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE,
            "source_head": identity_bridge.SOURCE_HEAD,
            "source_tree": identity_bridge.SOURCE_TREE,
            "source_runtime_identity_sha256": identity_bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
            "runtime_dependency_contract_sha256": contract_sha,
            "runtime_dependency_acceptance_receipt_sha256": identity_bridge._canonical_sha256(acceptance_receipt),
        },
        "protocol": {
            "time": TIME,
            "grid_resolution": GRID_RESOLUTION,
            "box": list(BOX),
            "cartesian_curl_order": 2,
            "integration_rule": "tensor_product_trapezoid",
            "source_numeric_targets_used": False,
            "renderer_or_camera_used": False,
            "pixel_loss_used": False,
        },
        "public_observable": dict(PUBLIC_OBSERVABLE),
        "metric_provenance": dict(METRIC_PROVENANCE),
        "external_method_screen": dict(EXTERNAL_METHOD_SCREEN),
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "direct_contribution": (
            "adds the missing candidate-bound axial/vorticity aspect observable for the same frozen ST052-M [u,v,w] identity"
        ),
        "remaining_limits": [
            "the official public source supplies no numerical axial-aspect target",
            "this autonomous aspect proxy does not establish visual/source correspondence",
            "no streamline trajectory proof or camera/frame registration is inferred",
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
    print(
        json.dumps(
            {
                "candidate_sha256": receipt["candidate_identity"]["candidate_sha256"],
                "velocity_identity_sha256": receipt["candidate_identity"]["velocity_identity_sha256"],
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
