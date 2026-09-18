"""Governed cross-language visualization grid for frozen ST048-S.

This module samples only the pinned, reconstructible ST048-S challenger and
writes a self-describing NetCDF grid for MATLAB/Python visualization.  It is a
delivery artifact, not PDE validation and not evidence of hidden OpenAI data.
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
from scipy.io import netcdf_file

TASK_ID = "CR-A9-045"
SCHEMA = "st048s_velocity_netcdf_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
SOURCE_PR = 390
SOURCE_HEAD_SHA = "97695a86f85ce68fb4ae70c41fc81c904d655183"
CANDIDATE_ID = "ST048-S"
PARENT_ID = "ST047-E"
PARENT_RAW_SHA256 = "dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0"
CANDIDATE_RAW_SHA256 = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
REFERENCE_TIMES = np.asarray([0.25, 0.50, 0.75], dtype=np.float64)
EVALUATION_BOUNDS = (-2.0, 2.0)
DEFAULT_RESOLUTION = 33
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


def _validate_resolution(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("resolution must be an odd integer")
    value = int(value)
    if value < 5 or value > 129 or value % 2 == 0:
        raise ValueError("resolution must be odd and lie in [5,129]")
    return value


def _hash_array(hasher, name: str, value: np.ndarray) -> None:
    array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
    hasher.update(name.encode("ascii"))
    hasher.update(b"\0")
    hasher.update(np.asarray(array.shape, dtype="<i8").tobytes())
    hasher.update(array.tobytes(order="C"))


def _new_grid_hasher(time: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray):
    hasher = hashlib.sha256()
    for value in (SCHEMA, SOURCE_HEAD_SHA, CANDIDATE_ID, CANDIDATE_RAW_SHA256, "time,x,y,z;u,v,w"):
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


def _git_head(repository_root: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def load_st048s_velocity(candidate_root: str | Path) -> tuple[VelocityCallable, dict]:
    """Reconstruct ST048-S only from the exact pinned PR #390 tree."""
    root = Path(candidate_root).resolve()
    if _git_head(root) != SOURCE_HEAD_SHA:
        raise ValueError("ST048 source checkout is not the pinned audited head")
    experiment = root / "experiments" / "root_st048"
    records = json.loads((experiment / "recipes.json").read_text())
    record = records.get(CANDIDATE_ID)
    if not isinstance(record, dict):
        raise ValueError("missing frozen ST048-S recipe")
    if (
        record.get("parent_id") != PARENT_ID
        or record.get("parent_original_raw_sha256") != PARENT_RAW_SHA256
        or record.get("original_local_raw_sha256") != CANDIDATE_RAW_SHA256
        or record.get("pde_validated") is not False
        or record.get("source_correspondence_verified") is not False
    ):
        raise ValueError("ST048-S recipe provenance/truth state drift")

    previous_path = list(sys.path)
    stale_modules = (
        "replay_st048", "boundary_shear", "replay_st047", "continuation",
        "structure_and_particles", "replay_st046", "acceleration_fit",
        "mechanism_audit", "localized_model", "spacetime", "validate",
    )
    try:
        sys.path.insert(0, str(experiment))
        for name in stale_modules:
            sys.modules.pop(name, None)
        replay = importlib.import_module("replay_st048")
        family, raw = replay.reconstruct(CANDIDATE_ID)
    finally:
        sys.path[:] = previous_path

    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be finite with shape (n,3)")
        times = np.full(len(points), float(time), dtype=float)
        values, _pressure = family.fields(raw, points, times)
        values = np.asarray(values, dtype=np.float64)
        if values.shape != points.shape or not np.isfinite(values).all():
            raise ValueError("ST048-S returned invalid velocity values")
        return values

    return velocity, {
        "source_pr": SOURCE_PR,
        "source_head_sha": SOURCE_HEAD_SHA,
        "candidate_id": CANDIDATE_ID,
        "parent_id": PARENT_ID,
        "parent_original_raw_sha256": PARENT_RAW_SHA256,
        "original_raw_candidate_sha256": CANDIDATE_RAW_SHA256,
        "recipe_modifiers_sha256": record.get("modifiers_sha256"),
        "reconstruction_note": (
            "Mathematical ST048-S is reconstructed from the exact pinned readable recipe. "
            "Regenerated JSON metadata is not claimed byte-identical to the archived raw child."
        ),
        "pde_validated": False,
        "source_correspondence_verified": False,
    }


def _export_velocity_netcdf(
    velocity: VelocityCallable,
    path: str | os.PathLike[str],
    *,
    resolution: int,
    overwrite: bool = False,
) -> dict:
    resolution = _validate_resolution(resolution)
    target = Path(path)
    if target.suffix.lower() != ".nc":
        raise ValueError("output path must end in .nc")
    if target.exists() and not overwrite:
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)

    axis = np.linspace(EVALUATION_BOUNDS[0], EVALUATION_BOUNDS[1], resolution, dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    time = REFERENCE_TIMES.copy()
    hasher = _new_grid_hasher(time, axis, axis, axis)
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
    speed_rms: list[float] = []
    speed_max: list[float] = []

    try:
        with netcdf_file(str(temporary), mode="w", version=2) as nc:
            nc.createDimension("time", len(time))
            nc.createDimension("x", resolution)
            nc.createDimension("y", resolution)
            nc.createDimension("z", resolution)
            for name, data in (("time", time), ("x", axis), ("y", axis), ("z", axis)):
                variable = nc.createVariable(name, "d", (name,))
                variable[:] = data
            u_var = nc.createVariable("u", "d", ("time", "x", "y", "z"))
            v_var = nc.createVariable("v", "d", ("time", "x", "y", "z"))
            w_var = nc.createVariable("w", "d", ("time", "x", "y", "z"))
            speed_var = nc.createVariable("speed", "d", ("time", "x", "y", "z"))

            any_nonzero = False
            for index, current_time in enumerate(time):
                values = np.asarray(velocity(points, float(current_time)), dtype=np.float64)
                if values.shape != points.shape or not np.isfinite(values).all():
                    raise ValueError("velocity returned invalid values")
                values = values.reshape(resolution, resolution, resolution, 3)
                speed = np.linalg.norm(values, axis=-1)
                u_var[index] = values[..., 0]
                v_var[index] = values[..., 1]
                w_var[index] = values[..., 2]
                speed_var[index] = speed
                for component, label in enumerate(("u", "v", "w")):
                    _hash_array(hasher, f"{label}[{index}]", values[..., component])
                speed_rms.append(float(np.sqrt(np.mean(speed * speed))))
                speed_max.append(float(np.max(speed)))
                any_nonzero = any_nonzero or bool(np.any(values != 0.0))
            if not any_nonzero:
                raise ValueError("exact-zero sampled velocity cannot be exported")

            digest = hasher.hexdigest()
            nc.schema = SCHEMA
            nc.task_id = TASK_ID
            nc.source_pr = SOURCE_PR
            nc.source_head_sha = SOURCE_HEAD_SHA
            nc.candidate = CANDIDATE_ID
            nc.candidate_raw_sha256 = CANDIDATE_RAW_SHA256
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
        "times": [float(v) for v in time],
        "grid_sha256": digest,
        "speed_rms_by_time": speed_rms,
        "speed_max_by_time": speed_max,
    }


def load_st048s_netcdf(path: str | os.PathLike[str]) -> dict:
    """Fail-closed validation and in-memory load of one ST048-S grid."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    with netcdf_file(str(source), mode="r", mmap=False) as nc:
        expected_attrs = {
            "schema": SCHEMA,
            "task_id": TASK_ID,
            "source_head_sha": SOURCE_HEAD_SHA,
            "candidate": CANDIDATE_ID,
            "candidate_raw_sha256": CANDIDATE_RAW_SHA256,
            "dimension_order": "time,x,y,z",
            "component_order": "u,v,w",
            "delivery_semantics": "candidate_specific_visualization_grid_not_pde_evidence",
        }
        for key, expected in expected_attrs.items():
            if _text(getattr(nc, key, "")) != expected:
                raise ValueError(f"NetCDF metadata drift: {key}")
        if int(getattr(nc, "source_pr", -1)) != SOURCE_PR or int(getattr(nc, "velocity_export_ready", 0)) != 1:
            raise ValueError("source/export metadata drift")
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
        if not np.all(np.diff(axis) > 0.0) or axis[0] != -2.0 or axis[-1] != 2.0:
            raise ValueError("coordinate axis leaves registered box")
    expected_shape = (3, resolution, resolution, resolution)
    for name in ("u", "v", "w", "speed"):
        if arrays[name].shape != expected_shape or not np.isfinite(arrays[name]).all():
            raise ValueError(f"invalid {name} array")
    recomputed_speed = np.sqrt(arrays["u"] ** 2 + arrays["v"] ** 2 + arrays["w"] ** 2)
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, float(np.max(recomputed_speed)))
    if float(np.max(np.abs(arrays["speed"] - recomputed_speed))) > tolerance:
        raise ValueError("stored speed is inconsistent with [u,v,w]")
    if float(np.max(recomputed_speed)) == 0.0:
        raise ValueError("exact-zero grid is invalid")

    hasher = _new_grid_hasher(time, x, y, z)
    for index in range(3):
        for label in ("u", "v", "w"):
            _hash_array(hasher, f"{label}[{index}]", arrays[label][index])
    if len(stored_digest) != 64 or hasher.hexdigest() != stored_digest:
        raise ValueError("NetCDF grid checksum mismatch")
    arrays["grid_sha256"] = stored_digest
    return arrays


def export_st048s_visualization_grid(
    candidate_root: str | Path,
    output: str | os.PathLike[str],
    receipt: str | os.PathLike[str],
    *,
    resolution: int = DEFAULT_RESOLUTION,
) -> dict:
    """Reconstruct pinned ST048-S, export the grid, verify it, and write a receipt."""
    velocity, metadata = load_st048s_velocity(candidate_root)
    result = _export_velocity_netcdf(velocity, output, resolution=resolution)
    loaded = load_st048s_netcdf(output)
    if loaded["grid_sha256"] != result["grid_sha256"]:
        raise ValueError("post-write grid identity drift")
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "candidate": metadata,
        "registered_contract": {
            "nu": 0.01,
            "physical_domain": "R^3",
            "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
            "support": "r < 2 and abs(z) < 2 with smooth zero extension",
            "time_interval": [0.25, 0.75],
            "reference_times": [0.25, 0.50, 0.75],
            "forcing": "unchanged preregistered restricted two-parameter family; not altered by export",
            "reference_energy": "E(0.25)=1 +/- 0.001; not altered by export",
            "validation": "independent PDE validation remains separate from this sampled visualization grid",
        },
        "grid": result,
        "matlab_python_usage": {
            "format": "NetCDF-3 64-bit-offset",
            "dimensions": "time,x,y,z",
            "variables": ["u", "v", "w", "speed"],
            "matlab": "ncinfo / ncread",
            "python": "scipy.io.netcdf_file or compatible NetCDF reader",
        },
        "external_method": {
            "classification": "direct_migration_public_api_only",
            "source_repo": "scipy/scipy",
            "source_commit": "b12c772edbc1fe0d3db9481cdfcc2e311569cb25",
            "license": "BSD-3-Clause",
            "scope": "scipy.io.netcdf_file for cross-language NetCDF serialization only",
            "copied_source_code": False,
            "dependency_delta": "none",
        },
        "truth_boundary": {
            "grid_derivatives_are_pde_evidence": False,
            "openai_numerical_velocity_used": False,
            "hidden_time_alignment_used": False,
            "camera_registration_used": False,
            "image_fit_used": False,
            "visual_acceptance_threshold_defined": False,
            "free_residual_force_used": False,
            "collapsed_velocity_accepted": False,
            **{flag: False for flag in _FALSE_TRUTH_FLAGS},
        },
    }
    payload["report_sha256"] = _canonical_hash(payload)
    receipt_path = Path(receipt)
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--resolution", type=int, default=DEFAULT_RESOLUTION)
    args = parser.parse_args()
    report = export_st048s_visualization_grid(
        args.candidate_root, args.output, args.receipt, resolution=args.resolution
    )
    print(json.dumps({
        "task_id": report["task_id"],
        "grid_sha256": report["grid"]["grid_sha256"],
        "report_sha256": report["report_sha256"],
        "speed_rms_by_time": report["grid"]["speed_rms_by_time"],
        "speed_max_by_time": report["grid"]["speed_max_by_time"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
