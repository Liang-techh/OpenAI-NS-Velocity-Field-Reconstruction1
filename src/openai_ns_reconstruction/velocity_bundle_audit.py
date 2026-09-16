"""Fail-closed reproducibility audit for an exported velocity-field bundle.

The audit verifies that a saved visualization grid is exactly tied to a frozen
candidate and can be regenerated through the public ``VelocityField`` API.  It
is an artifact-integrity check only: passing does not imply PDE acceptance,
OpenAI-field identification, or paper-exact reconstruction.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


REQUIRED_GRID_KEYS = frozenset({"x", "y", "z", "times", "u", "v", "w"})
PUBLIC_METADATA_KEYS = (
    "family",
    "components",
    "coordinates",
    "units",
    "time_interval",
    "compact_support",
    "visual_correspondence",
    "pde_acceptance",
    "grid_layout",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_velocity_bundle(field, artifact_dir, *, atol: float = 5e-13) -> dict:
    """Verify a saved ``grid.npz`` against its metadata and public evaluator.

    ``field`` is intentionally duck-typed.  It must expose ``sha256``,
    ``metadata()`` and ``grid(x, y, z, times)``.  This keeps the audit separate
    from optimizer/training internals and makes it usable for future candidate
    families that obey the same public delivery contract.
    """
    artifact_dir = Path(artifact_dir)
    metadata_path = artifact_dir / "metadata.json"
    grid_path = artifact_dir / "grid.npz"
    _require(metadata_path.is_file(), "missing metadata.json")
    _require(grid_path.is_file(), "missing grid.npz")
    _require(np.isfinite(atol) and atol >= 0.0, "atol must be finite and nonnegative")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    public_metadata = field.metadata()

    for key in PUBLIC_METADATA_KEYS:
        _require(key in metadata, f"metadata missing {key}")
        _require(key in public_metadata, f"public evaluator metadata missing {key}")
        _require(metadata[key] == public_metadata[key], f"metadata mismatch for {key}")

    _require("candidate_sha256" in metadata, "metadata missing candidate_sha256")
    _require(
        metadata["candidate_sha256"] == field.sha256,
        f"candidate SHA256 mismatch: metadata={metadata['candidate_sha256']} field={field.sha256}",
    )
    actual_grid_sha256 = _sha256(grid_path)
    _require(
        metadata.get("grid_npz_sha256") == actual_grid_sha256,
        "grid.npz SHA256 mismatch",
    )

    with np.load(grid_path, allow_pickle=False) as grid:
        _require(REQUIRED_GRID_KEYS.issubset(grid.files), "grid.npz missing required arrays")
        x = np.asarray(grid["x"], dtype=float)
        y = np.asarray(grid["y"], dtype=float)
        z = np.asarray(grid["z"], dtype=float)
        times = np.asarray(grid["times"], dtype=float)
        stored = np.stack((grid["u"], grid["v"], grid["w"]), axis=-1)

    for name, axis in (("x", x), ("y", y), ("z", z), ("times", times)):
        _require(axis.ndim == 1 and axis.size > 0, f"{name} must be a nonempty 1D array")
        _require(np.all(np.isfinite(axis)), f"{name} contains nonfinite values")
        if axis.size > 1:
            _require(np.all(np.diff(axis) > 0), f"{name} must be strictly increasing")

    expected_shape = (times.size, x.size, y.size, z.size, 3)
    _require(stored.shape == expected_shape, "stored component array shape mismatch")
    _require(metadata.get("grid_shape") == list(expected_shape), "metadata grid_shape mismatch")
    _require(np.all(np.isfinite(stored)), "stored velocity grid contains nonfinite values")

    t0, t1 = map(float, metadata["time_interval"])
    _require(np.all((times >= t0) & (times <= t1)), "saved reference times leave time_interval")

    regenerated = np.asarray(field.grid(x, y, z, times), dtype=float)
    _require(regenerated.shape == stored.shape, "public evaluator regenerated wrong grid shape")
    _require(np.all(np.isfinite(regenerated)), "public evaluator regenerated nonfinite values")
    max_abs_error = float(np.max(np.abs(regenerated - stored)))
    _require(max_abs_error <= atol, "saved grid does not reproduce through public evaluator")

    return {
        "bundle_schema_version": 1,
        "status": "passed",
        "family": metadata["family"],
        "candidate_sha256": field.sha256,
        "grid_npz_sha256": actual_grid_sha256,
        "grid_shape": list(expected_shape),
        "reference_times": [float(t) for t in times],
        "max_abs_regeneration_error": max_abs_error,
        "truth_boundary": {
            "velocity_export_ready": True,
            "visualization_ready": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_correspondence_verified": False,
            "note": "artifact reproducibility only; stronger claims require independent evidence",
        },
    }
