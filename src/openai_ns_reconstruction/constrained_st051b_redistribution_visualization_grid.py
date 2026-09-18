"""Governed cross-language visualization grid for the frozen ST051-B swirl child.

The exported field is the exact ST051-B reconstruction from PR #460 followed by
Agent 7 PR #469's frozen radial swirl redistribution.  The grid is a delivery
artifact for MATLAB/Python visualization; it is not Navier--Stokes validation,
OpenAI-field identification, or visual-correspondence evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.io import netcdf_file

TASK_ID = "CR-A9-053"
SCHEMA = "st051b_frozen_redistribution_velocity_netcdf_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
PARENT_PR = 460
PARENT_HEAD_SHA = "4b784f1b8457af2ead49295631d834d4e882000b"
PARENT_ID = "ST051-B"
PARENT_PARENT_ID = "ST050R-P"
PARENT_PARENT_RAW_SHA256 = "f5b5a2869991387e503a522b599af86ea814ba02bdae2f8b56fd16674f7dab19"
PARENT_REPORTED_RAW_SHA256 = "0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d"
TRANSFER_PR = 469
TRANSFER_HEAD_SHA = "95f94188270f02511c8bdf660822f77530c9983a"
TRANSFER_SOURCE_PR = 458
SOURCE_ALPHA = 2.520520814687742
GAIN = 0.025
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
REFERENCE_TIMES = np.asarray([0.25, 0.50, 0.75], dtype=np.float64)
EVALUATION_BOUNDS = (-2.0, 2.0)
DEFAULT_RESOLUTION = 33
ENERGY_QUADRATURE_ORDER = 72
EXPECTED_REFERENCE_SCALE = 1.0014791925672812
SCIPY_SOURCE_COMMIT = "b12c772edbc1fe0d3db9481cdfcc2e311569cb25"
_FALSE_TRUTH_FLAGS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "source_correspondence_verified",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)
VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


def _canonical_hash(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _git_head(repository_root: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _validate_resolution(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("resolution must be an odd integer")
    value = int(value)
    if value < 5 or value > 129 or value % 2 == 0:
        raise ValueError("resolution must be odd and lie in [5,129]")
    return value


def compact_bump(r, lo: float, hi: float) -> np.ndarray:
    """The exact C-infinity radial bump used by the frozen Agent-7 transfer."""
    r = np.asarray(r, dtype=float)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    x = (r - mid) / half
    out = np.zeros_like(r)
    mask = np.abs(x) < 1.0
    xm = x[mask]
    out[mask] = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    return out


def redistribution_profile(r) -> np.ndarray:
    return compact_bump(r, *INNER_WINDOW) - SOURCE_ALPHA * compact_bump(r, *OUTER_WINDOW)


def _transform_velocity_values(points: np.ndarray, values: np.ndarray, *, gain: float, scale: float) -> np.ndarray:
    """Apply the frozen cylindrical swirl redistribution to Cartesian values."""
    points = np.asarray(points, dtype=float)
    values = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or values.shape != points.shape:
        raise ValueError("points and velocity values must have shape (n,3)")
    if not np.isfinite(points).all() or not np.isfinite(values).all():
        raise ValueError("points and velocity values must be finite")
    out = values.copy()
    x, y = points[:, 0], points[:, 1]
    radius = np.hypot(x, y)
    mask = radius > 1e-14
    if np.any(mask):
        rx, ry = x[mask] / radius[mask], y[mask] / radius[mask]
        ux, uy = out[mask, 0].copy(), out[mask, 1].copy()
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        ut *= 1.0 + float(gain) * redistribution_profile(radius[mask])
        out[mask, 0] = rx * ur - ry * ut
        out[mask, 1] = ry * ur + rx * ut
    return float(scale) * out


def _load_parent(candidate_root: str | Path):
    """Load the exact frozen ST051-B reconstruction, rejecting source drift."""
    root = Path(candidate_root).resolve()
    if _git_head(root) != PARENT_HEAD_SHA:
        raise ValueError("ST051 source checkout is not the pinned audited head")
    experiment = root / "experiments" / "root_st051"
    recipe_path = experiment / "recipes" / f"{PARENT_ID}.json"
    if not recipe_path.is_file():
        raise ValueError("missing frozen ST051-B recipe")
    record = json.loads(recipe_path.read_text())
    if (
        record.get("parent_id") != PARENT_PARENT_ID
        or record.get("parent_original_raw_sha256") != PARENT_PARENT_RAW_SHA256
        or record.get("edge_modes") is not True
        or record.get("pde_validated") is not False
        or record.get("source_correspondence_verified") is not False
    ):
        raise ValueError("ST051-B recipe provenance/truth state drift")
    modifiers_sha = record.get("modifiers_sha256")
    if not isinstance(modifiers_sha, str) or len(modifiers_sha) != 64:
        raise ValueError("missing ST051-B modifier identity")

    previous_path = list(sys.path)
    stale_modules = (
        "replay_st051", "aligned_continuation", "audit_st051", "replay_recovery",
        "replay_st050", "replay_st048", "replay_st047", "continuation",
        "boundary_shear", "spacetime", "validate", "diagnostics",
    )
    try:
        sys.path.insert(0, str(experiment))
        for name in stale_modules:
            sys.modules.pop(name, None)
        replay = importlib.import_module("replay_st051")
        family, raw = replay.reconstruct(PARENT_ID)
    finally:
        sys.path[:] = previous_path
    return family, raw, modifiers_sha


def _parent_velocity(family, raw, points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    times = np.full(len(points), float(time), dtype=float)
    values, _pressure = family.fields(raw, points, times)
    values = np.asarray(values, dtype=np.float64)
    if values.shape != points.shape or not np.isfinite(values).all():
        raise ValueError("ST051-B returned invalid velocity values")
    return values


def _reference_energy_and_child_scale(family, raw, *, order: int = ENERGY_QUADRATURE_ORDER) -> tuple[float, float, float]:
    """Recompute the one common reference-energy scale used by PR #469."""
    nodes_r, weights_r = leggauss(int(order))
    nodes_z, weights_z = leggauss(int(order))
    radius = nodes_r + 1.0
    axial = 2.0 * nodes_z
    rr, zz = np.meshgrid(radius, axial, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    parent = _parent_velocity(family, raw, points, 0.25).reshape(rr.shape + (3,))
    weights = weights_r[:, None] * (2.0 * weights_z[None, :]) * np.pi * rr
    parent_energy = float(np.sum(weights * np.sum(parent * parent, axis=-1)))
    child = parent.copy()
    child[..., 1] *= 1.0 + GAIN * redistribution_profile(rr)
    child_energy = float(np.sum(weights * np.sum(child * child, axis=-1)))
    if not (parent_energy > 0.0 and child_energy > 0.0):
        raise ValueError("nonpositive reference energy")
    scale = float(np.sqrt(parent_energy / child_energy))
    return parent_energy, child_energy, scale


def load_st051b_redistribution_velocity(candidate_root: str | Path) -> tuple[VelocityCallable, dict]:
    """Return the exact frozen visualization child as a callable velocity."""
    family, raw, modifiers_sha = _load_parent(candidate_root)
    parent_energy, raw_child_energy, scale = _reference_energy_and_child_scale(family, raw)
    if abs(scale - EXPECTED_REFERENCE_SCALE) > 5e-12:
        raise ValueError("frozen redistribution normalization drift")

    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be finite with shape (n,3)")
        if not np.isfinite(time) or not (0.25 <= float(time) <= 0.75):
            raise ValueError("time must lie in the registered [0.25,0.75] window")
        parent = _parent_velocity(family, raw, points, float(time))
        values = _transform_velocity_values(points, parent, gain=GAIN, scale=scale)
        if not np.isfinite(values).all():
            raise ValueError("transformed velocity is nonfinite")
        return values

    metadata = {
        "parent_pr": PARENT_PR,
        "parent_head_sha": PARENT_HEAD_SHA,
        "parent_candidate_id": PARENT_ID,
        "parent_reported_raw_sha256": PARENT_REPORTED_RAW_SHA256,
        "parent_recipe_modifiers_sha256": modifiers_sha,
        "transfer_pr": TRANSFER_PR,
        "transfer_head_sha": TRANSFER_HEAD_SHA,
        "transfer_source_pr": TRANSFER_SOURCE_PR,
        "source_alpha": SOURCE_ALPHA,
        "gain": GAIN,
        "inner_window": list(INNER_WINDOW),
        "outer_window": list(OUTER_WINDOW),
        "parent_reference_energy": parent_energy,
        "raw_child_reference_energy": raw_child_energy,
        "common_reference_energy_scale": scale,
        "reconstruction_note": (
            "Mathematical ST051-B is reconstructed from the exact pinned readable recipe; "
            "the frozen PR #469 swirl transform is independently reimplemented without refitting."
        ),
        "pde_validated": False,
        "source_correspondence_verified": False,
    }
    return velocity, metadata


def _hash_array(hasher, name: str, value: np.ndarray) -> None:
    array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
    hasher.update(name.encode("ascii"))
    hasher.update(b"\0")
    hasher.update(np.asarray(array.shape, dtype="<i8").tobytes())
    hasher.update(array.tobytes(order="C"))


def _new_grid_hasher(time: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray):
    hasher = hashlib.sha256()
    for value in (SCHEMA, PARENT_HEAD_SHA, TRANSFER_HEAD_SHA, PARENT_ID, "time,x,y,z;u,v,w"):
        hasher.update(value.encode("ascii"))
        hasher.update(b"\0")
    for name, array in (("time", time), ("x", x), ("y", y), ("z", z)):
        _hash_array(hasher, name, array)
    return hasher


def _text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("ascii")
    return str(value)


def export_velocity_netcdf(
    velocity: VelocityCallable,
    path: str | os.PathLike[str],
    *,
    resolution: int = DEFAULT_RESOLUTION,
    overwrite: bool = False,
) -> dict:
    """Export one candidate-specific NetCDF3 64-bit-offset visualization grid."""
    resolution = _validate_resolution(resolution)
    target = Path(path)
    if target.suffix.lower() != ".nc":
        raise ValueError("output path must end in .nc")
    if target.exists() and not overwrite:
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)

    axis = np.linspace(*EVALUATION_BOUNDS, resolution, dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    times = REFERENCE_TIMES.copy()
    hasher = _new_grid_hasher(times, axis, axis, axis)
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
    speed_rms: list[float] = []
    speed_max: list[float] = []

    try:
        with netcdf_file(str(temporary), mode="w", version=2) as nc:
            nc.createDimension("time", len(times))
            nc.createDimension("x", resolution)
            nc.createDimension("y", resolution)
            nc.createDimension("z", resolution)
            for name, data in (("time", times), ("x", axis), ("y", axis), ("z", axis)):
                var = nc.createVariable(name, "d", (name,))
                var[:] = data
            variables = {
                name: nc.createVariable(name, "d", ("time", "x", "y", "z"))
                for name in ("u", "v", "w", "speed")
            }
            any_nonzero = False
            for index, current_time in enumerate(times):
                values = np.asarray(velocity(points, float(current_time)), dtype=np.float64)
                if values.shape != points.shape or not np.isfinite(values).all():
                    raise ValueError("velocity returned invalid values")
                values = values.reshape(resolution, resolution, resolution, 3)
                speed = np.linalg.norm(values, axis=-1)
                for component, name in enumerate(("u", "v", "w")):
                    variables[name][index] = values[..., component]
                    _hash_array(hasher, f"{name}[{index}]", values[..., component])
                variables["speed"][index] = speed
                speed_rms.append(float(np.sqrt(np.mean(speed * speed))))
                speed_max.append(float(np.max(speed)))
                any_nonzero = any_nonzero or bool(np.any(values != 0.0))
            if not any_nonzero:
                raise ValueError("exact-zero sampled velocity cannot be exported")

            digest = hasher.hexdigest()
            nc.schema = SCHEMA
            nc.task_id = TASK_ID
            nc.parent_pr = PARENT_PR
            nc.parent_head_sha = PARENT_HEAD_SHA
            nc.transfer_pr = TRANSFER_PR
            nc.transfer_head_sha = TRANSFER_HEAD_SHA
            nc.candidate = "ST051-B+frozen-swirl-redistribution"
            nc.grid_sha256 = digest
            nc.dimension_order = "time,x,y,z"
            nc.component_order = "u,v,w"
            nc.registered_box = "[-2,2]^3"
            nc.registered_time_window = "[0.25,0.75]"
            nc.delivery_semantics = "candidate_specific_visualization_grid_not_pde_evidence"
            nc.velocity_export_ready = 1
            for flag in _FALSE_TRUTH_FLAGS:
                setattr(nc, flag, 0)
        os.replace(temporary, target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return {
        "path": str(target),
        "resolution": resolution,
        "times": [float(value) for value in times],
        "grid_sha256": digest,
        "speed_rms_by_time": speed_rms,
        "speed_max_by_time": speed_max,
    }


def load_visualization_netcdf(path: str | os.PathLike[str]) -> dict:
    """Fail-closed reload and integrity validation for an exported grid."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    with netcdf_file(str(source), mode="r", mmap=False) as nc:
        expected = {
            "schema": SCHEMA,
            "task_id": TASK_ID,
            "parent_head_sha": PARENT_HEAD_SHA,
            "transfer_head_sha": TRANSFER_HEAD_SHA,
            "candidate": "ST051-B+frozen-swirl-redistribution",
            "dimension_order": "time,x,y,z",
            "component_order": "u,v,w",
            "delivery_semantics": "candidate_specific_visualization_grid_not_pde_evidence",
        }
        for key, value in expected.items():
            if _text(getattr(nc, key, "")) != value:
                raise ValueError(f"NetCDF metadata drift: {key}")
        if int(getattr(nc, "parent_pr", -1)) != PARENT_PR or int(getattr(nc, "transfer_pr", -1)) != TRANSFER_PR:
            raise ValueError("source PR metadata drift")
        if int(getattr(nc, "velocity_export_ready", 0)) != 1:
            raise ValueError("velocity export readiness missing")
        for flag in _FALSE_TRUTH_FLAGS:
            if int(getattr(nc, flag, 1)) != 0:
                raise ValueError(f"unsupported truth-state promotion: {flag}")
        required = {"time", "x", "y", "z", "u", "v", "w", "speed"}
        if set(nc.variables) != required:
            raise ValueError("unexpected NetCDF variable set")
        arrays = {name: np.array(nc.variables[name].data, dtype=np.float64, copy=True) for name in required}
        stored_digest = _text(getattr(nc, "grid_sha256", ""))

    time, x, y, z = (arrays[name] for name in ("time", "x", "y", "z"))
    if not np.array_equal(time, REFERENCE_TIMES):
        raise ValueError("reference times drift")
    resolution = _validate_resolution(len(x))
    for axis in (x, y, z):
        if axis.ndim != 1 or len(axis) != resolution or not np.isfinite(axis).all():
            raise ValueError("invalid coordinate axis")
        if axis[0] != -2.0 or axis[-1] != 2.0 or not np.all(np.diff(axis) > 0.0):
            raise ValueError("coordinate axis leaves registered box")
    expected_shape = (3, resolution, resolution, resolution)
    for name in ("u", "v", "w", "speed"):
        if arrays[name].shape != expected_shape or not np.isfinite(arrays[name]).all():
            raise ValueError(f"invalid {name} array")
    speed = np.sqrt(arrays["u"] ** 2 + arrays["v"] ** 2 + arrays["w"] ** 2)
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, float(np.max(speed)))
    if float(np.max(np.abs(arrays["speed"] - speed))) > tolerance:
        raise ValueError("stored speed is inconsistent with [u,v,w]")
    if float(np.max(speed)) == 0.0:
        raise ValueError("exact-zero grid is invalid")

    hasher = _new_grid_hasher(time, x, y, z)
    for index in range(3):
        for name in ("u", "v", "w"):
            _hash_array(hasher, f"{name}[{index}]", arrays[name][index])
    if len(stored_digest) != 64 or hasher.hexdigest() != stored_digest:
        raise ValueError("NetCDF grid checksum mismatch")
    arrays["grid_sha256"] = stored_digest
    return arrays


def export_st051b_redistribution_visualization_grid(
    candidate_root: str | Path,
    output: str | os.PathLike[str],
    receipt: str | os.PathLike[str],
    *,
    resolution: int = DEFAULT_RESOLUTION,
) -> dict:
    velocity, candidate = load_st051b_redistribution_velocity(candidate_root)
    result = export_velocity_netcdf(velocity, output, resolution=resolution)
    loaded = load_visualization_netcdf(output)
    if loaded["grid_sha256"] != result["grid_sha256"]:
        raise ValueError("post-write grid identity drift")
    truth = {
        "canonical_velocity_changed": False,
        "pressure_transferred": False,
        "forcing_transferred": False,
        "heldout_pde_residual_assessed": False,
        "production_candidate_selected": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "source_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "candidate": candidate,
        "registered_contract": {
            "nu": 0.01,
            "physical_domain": "R^3",
            "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
            "support": "r < 2 and abs(z) < 2 with smooth zero extension",
            "time_interval": [0.25, 0.75],
            "reference_energy": 1.0,
            "forcing": "restricted two-parameter family only; not transferred to this child",
            "divergence_gate": 1e-5,
            "momentum_max_gate": 1e-3,
            "momentum_L2_gate": 1e-3,
        },
        "external_method": {
            "classification": "direct migration / public API only",
            "source_repo": "scipy/scipy",
            "source_commit": SCIPY_SOURCE_COMMIT,
            "api": "scipy.io.netcdf_file (NetCDF3 64-bit-offset)",
            "license": "BSD-3-Clause",
            "copied_implementation": False,
        },
        "internal_method": {
            "classification": "direct internal method reuse / independent reimplementation",
            "source_pr": TRANSFER_PR,
            "source_head_sha": TRANSFER_HEAD_SHA,
            "migrated_scope": "frozen radial windows, alpha_C, gain, azimuthal-only redistribution, common E(.25) normalization",
            "differences": "export samples the fixed child on a Cartesian grid; it does not refit the transform or reuse visual targets",
        },
        "grid": result,
        "truth_boundary": truth,
    }
    payload["receipt_sha256"] = _canonical_hash(payload)
    target = Path(receipt)
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    export = sub.add_parser("export")
    export.add_argument("--candidate-root", required=True)
    export.add_argument("--output", required=True)
    export.add_argument("--receipt", required=True)
    export.add_argument("--resolution", type=int, default=DEFAULT_RESOLUTION)
    args = parser.parse_args(argv)
    if args.command == "export":
        payload = export_st051b_redistribution_visualization_grid(
            args.candidate_root, args.output, args.receipt, resolution=args.resolution
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
