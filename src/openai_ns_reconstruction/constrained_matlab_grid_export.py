"""Truth-bounded MATLAB interchange for frozen Cartesian velocity grids.

This module serializes an already-evaluated candidate grid.  It does not
select a candidate, alter velocity values, fit a camera, or establish visual
or PDE validity.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
import tempfile
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.io import loadmat, savemat


_SCHEMA_VERSION = 1
_LAYOUT = "time,x,y,z,component"
_COMPONENT_ORDER = "u,v,w"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_TRUTH_KEYS = (
    "velocity_changed",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


@dataclass(frozen=True)
class MatlabVelocityGrid:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    t: np.ndarray
    velocity: np.ndarray
    candidate_sha256: str
    candidate_provenance: str
    grid_provenance: str
    per_time_max_speed: np.ndarray
    per_time_rms_speed: np.ndarray
    truth_flags: Mapping[str, bool]


@dataclass(frozen=True)
class MatlabExportReceipt:
    path: str
    candidate_sha256: str
    velocity_shape: tuple[int, int, int, int, int]
    global_max_speed: float
    per_time_max_speed: tuple[float, ...]
    roundtrip_verified: bool
    truth_flags: Mapping[str, bool]


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _require_sha256(value: str) -> str:
    value = _require_text("candidate_sha256", value)
    if _SHA256_RE.fullmatch(value) is None:
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
    return value


def _axis(name: str, values: np.ndarray, *, minimum_size: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 1 or arr.size < minimum_size:
        raise ValueError(f"{name} must be one-dimensional with at least {minimum_size} entries")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if arr.size > 1 and not np.all(np.diff(arr) > 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return np.array(arr, dtype=np.float64, copy=True)


def _readonly(array: np.ndarray) -> np.ndarray:
    out = np.array(array, copy=True)
    out.setflags(write=False)
    return out


def _mat_text(raw: dict[str, object], key: str) -> str:
    if key not in raw:
        raise ValueError(f"MAT file is missing required variable {key!r}")
    arr = np.asarray(raw[key])
    if arr.dtype.kind not in "US":
        raise ValueError(f"MAT variable {key!r} must be text")
    flat = arr.astype(str).ravel()
    if flat.size == 0:
        raise ValueError(f"MAT variable {key!r} must not be empty")
    if flat.size == 1:
        return str(flat[0])
    return "".join(str(item) for item in flat)


def _mat_scalar_int(raw: dict[str, object], key: str) -> int:
    if key not in raw:
        raise ValueError(f"MAT file is missing required variable {key!r}")
    arr = np.asarray(raw[key])
    if arr.size != 1 or not np.issubdtype(arr.dtype, np.number):
        raise ValueError(f"MAT variable {key!r} must be a numeric scalar")
    value = float(arr.reshape(-1)[0])
    if not np.isfinite(value) or value != int(value):
        raise ValueError(f"MAT variable {key!r} must be an integer scalar")
    return int(value)


def _truth_payload() -> dict[str, np.ndarray]:
    return {
        f"truth_{key}": np.array([[0]], dtype=np.uint8)
        for key in _TRUTH_KEYS
    }


def _read_truth(raw: dict[str, object]) -> Mapping[str, bool]:
    truth: dict[str, bool] = {}
    for key in _TRUTH_KEYS:
        value = _mat_scalar_int(raw, f"truth_{key}")
        if value not in (0, 1):
            raise ValueError(f"truth_{key} must be 0 or 1")
        if value != 0:
            raise ValueError(
                f"truth_{key}=1 is forbidden for this serialization-only artifact"
            )
        truth[key] = False
    return MappingProxyType(truth)


def load_velocity_grid_mat(path: str | os.PathLike[str]) -> MatlabVelocityGrid:
    """Load and validate a trusted MAT file produced by this module."""

    file_path = Path(path)
    if file_path.suffix.lower() != ".mat":
        raise ValueError("MATLAB grid path must end in .mat")
    if not file_path.is_file():
        raise FileNotFoundError(file_path)

    raw = loadmat(
        file_path,
        appendmat=False,
        squeeze_me=False,
        chars_as_strings=True,
        verify_compressed_data_integrity=True,
    )
    if _mat_scalar_int(raw, "schema_version") != _SCHEMA_VERSION:
        raise ValueError("unsupported MATLAB grid schema_version")
    if _mat_text(raw, "layout") != _LAYOUT:
        raise ValueError("MATLAB grid layout metadata is not the governed layout")
    if _mat_text(raw, "component_order") != _COMPONENT_ORDER:
        raise ValueError("MATLAB grid component order is not u,v,w")

    x = _axis("x", np.asarray(raw.get("x", [])).reshape(-1), minimum_size=2)
    y = _axis("y", np.asarray(raw.get("y", [])).reshape(-1), minimum_size=2)
    z = _axis("z", np.asarray(raw.get("z", [])).reshape(-1), minimum_size=2)
    t = _axis("t", np.asarray(raw.get("t", [])).reshape(-1), minimum_size=1)

    expected = (t.size, x.size, y.size, z.size)
    components: list[np.ndarray] = []
    for key in ("u_txyz", "v_txyz", "w_txyz"):
        if key not in raw:
            raise ValueError(f"MAT file is missing required variable {key!r}")
        component = np.asarray(raw[key], dtype=np.float64)
        if component.shape != expected:
            raise ValueError(
                f"{key} has shape {component.shape}, expected {expected}"
            )
        if not np.all(np.isfinite(component)):
            raise ValueError(f"{key} contains non-finite values")
        components.append(component)
    velocity = np.stack(components, axis=-1)

    speed = np.linalg.norm(velocity, axis=-1)
    per_time_max = np.max(speed, axis=(1, 2, 3))
    per_time_rms = np.sqrt(np.mean(speed * speed, axis=(1, 2, 3)))
    if np.any(per_time_max == 0.0):
        raise ValueError("each exported time slice must be nonzero")

    candidate_sha256 = _require_sha256(_mat_text(raw, "candidate_sha256"))
    candidate_provenance = _require_text(
        "candidate_provenance", _mat_text(raw, "candidate_provenance")
    )
    grid_provenance = _require_text(
        "grid_provenance", _mat_text(raw, "grid_provenance")
    )
    truth_flags = _read_truth(raw)

    return MatlabVelocityGrid(
        x=_readonly(x),
        y=_readonly(y),
        z=_readonly(z),
        t=_readonly(t),
        velocity=_readonly(velocity),
        candidate_sha256=candidate_sha256,
        candidate_provenance=candidate_provenance,
        grid_provenance=grid_provenance,
        per_time_max_speed=_readonly(per_time_max),
        per_time_rms_speed=_readonly(per_time_rms),
        truth_flags=truth_flags,
    )


def export_velocity_grid_mat(
    path: str | os.PathLike[str],
    *,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    velocity: np.ndarray,
    candidate_sha256: str,
    candidate_provenance: str,
    grid_provenance: str,
    overwrite: bool = False,
    do_compression: bool = True,
) -> MatlabExportReceipt:
    """Write a frozen ``(time,x,y,z,component)`` grid to MATLAB Level-5 MAT.

    Serialization preserves the supplied numeric field exactly.  Scientific
    acceptance remains external to this function.
    """

    file_path = Path(path)
    if file_path.suffix.lower() != ".mat":
        raise ValueError("MATLAB grid path must end in .mat")
    if not file_path.parent.is_dir():
        raise FileNotFoundError(file_path.parent)
    if file_path.exists() and not overwrite:
        raise FileExistsError(file_path)

    x_arr = _axis("x", x, minimum_size=2)
    y_arr = _axis("y", y, minimum_size=2)
    z_arr = _axis("z", z, minimum_size=2)
    t_arr = _axis("t", t, minimum_size=1)

    values = np.asarray(velocity, dtype=np.float64)
    expected_shape = (t_arr.size, x_arr.size, y_arr.size, z_arr.size, 3)
    if values.shape != expected_shape:
        raise ValueError(
            f"velocity has shape {values.shape}, expected {expected_shape}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity must contain only finite values")

    speed = np.linalg.norm(values, axis=-1)
    per_time_max = np.max(speed, axis=(1, 2, 3))
    per_time_rms = np.sqrt(np.mean(speed * speed, axis=(1, 2, 3)))
    if np.any(per_time_max == 0.0):
        raise ValueError("each exported time slice must be nonzero")

    sha = _require_sha256(candidate_sha256)
    candidate_provenance = _require_text(
        "candidate_provenance", candidate_provenance
    )
    grid_provenance = _require_text("grid_provenance", grid_provenance)

    payload: dict[str, object] = {
        "schema_version": np.array([[ _SCHEMA_VERSION ]], dtype=np.int32),
        "x": x_arr.reshape(-1, 1),
        "y": y_arr.reshape(-1, 1),
        "z": z_arr.reshape(-1, 1),
        "t": t_arr.reshape(-1, 1),
        "u_txyz": np.array(values[..., 0], copy=True),
        "v_txyz": np.array(values[..., 1], copy=True),
        "w_txyz": np.array(values[..., 2], copy=True),
        "layout": _LAYOUT,
        "component_order": _COMPONENT_ORDER,
        "candidate_sha256": sha,
        "candidate_provenance": candidate_provenance,
        "grid_provenance": grid_provenance,
        "per_time_max_speed": per_time_max.reshape(-1, 1),
        "per_time_rms_speed": per_time_rms.reshape(-1, 1),
        "global_max_speed": np.array([[float(np.max(per_time_max))]]),
        "scientific_claim_scope": "serialization_only",
        "mat_format": "MATLAB Level-5 via scipy.io.savemat",
    }
    payload.update(_truth_payload())

    temp_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{file_path.stem}.", suffix=".mat", dir=file_path.parent, delete=False
        ) as handle:
            temp_name = handle.name
        savemat(
            temp_name,
            payload,
            appendmat=False,
            format="5",
            long_field_names=True,
            do_compression=bool(do_compression),
            oned_as="column",
        )

        replay = load_velocity_grid_mat(temp_name)
        if not (
            np.array_equal(replay.x, x_arr)
            and np.array_equal(replay.y, y_arr)
            and np.array_equal(replay.z, z_arr)
            and np.array_equal(replay.t, t_arr)
            and np.array_equal(replay.velocity, values)
            and replay.candidate_sha256 == sha
            and replay.candidate_provenance == candidate_provenance
            and replay.grid_provenance == grid_provenance
        ):
            raise RuntimeError("MATLAB grid round-trip verification failed")

        if file_path.exists() and not overwrite:
            raise FileExistsError(file_path)
        os.replace(temp_name, file_path)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass

    truth_flags = MappingProxyType({key: False for key in _TRUTH_KEYS})
    return MatlabExportReceipt(
        path=str(file_path),
        candidate_sha256=sha,
        velocity_shape=expected_shape,
        global_max_speed=float(np.max(per_time_max)),
        per_time_max_speed=tuple(float(v) for v in per_time_max),
        roundtrip_verified=True,
        truth_flags=truth_flags,
    )
