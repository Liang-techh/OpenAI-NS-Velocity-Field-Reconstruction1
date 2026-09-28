"""Compare the fitted and solenoidal-generator time derivatives in NS momentum.

The latest ``localized_drift_1800_fit`` candidate is reconstructed as an actual
field, including its pressure corrections and two local 180-control patches.
On the independent 18,720-point grid this module computes the Cartesian
five-point spatial and physical-time jets directly.  It then replaces only
the fitted ``u_t`` by the frozen-time solenoidal-scale generator.  Pressure,
the pressure gradient, convection, and viscous Laplacian remain those of the
same actual field in both residuals.

This is a bounded reference-time residual diagnostic.  It does not certify a
forcing, PDE solution, time trajectory, or scale recursion.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from enriched_shape_replay import build_field  # noqa: E402
from full_wave_tangent import LocalPotentialField, _unpack_full  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from scale_transport_generator import (  # noqa: E402
    MEAN_ANGLES,
    VelocityAdapter,
    axisymmetric_swirl_mean,
    solenoidal_scale_generator,
)


DRIFT_PATH = ROOT / "localized_drift_1800_fit.json"
DRIFT_SOURCE_PATH = ROOT / "localized_drift_1800_fit.py"
GRID_PATH = ROOT / "full_wave_dense_tangent.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
GENERATOR_PATH = ROOT / "scale_transport_generator.py"
REPLAY_PATH = ROOT / "enriched_shape_replay.py"
MEAN_PATH = ROOT / "broad_meridional_constrained.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
OUTPUT_PATH = ROOT / "scale_generator_momentum_defect.json"
CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
CHUNK_SIZE = 768


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _metric(values, weights, points):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    square = float(np.sum(weights * np.sum(values * values, axis=1)))
    index = int(np.argmax(norms))
    return {
        "point_count": int(len(values)),
        "volume_L2": float(np.sqrt(max(square, 0.0))),
        "volume_RMS": float(np.sqrt(max(square, 0.0) / max(float(np.sum(weights)), 1.0e-300))),
        "max_norm": float(norms[index]),
        "max_point": points[index].tolist(),
        "max_index": index,
        "max_vector": values[index].tolist(),
        "physical_volume": float(np.sum(weights)),
    }


def _array_sha256(*arrays) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _load_grid():
    grid_report = json.loads(GRID_PATH.read_text(encoding="utf-8"))
    cache = grid_report.get("new_frozen_cache")
    if not isinstance(cache, dict):
        raise ValueError("full_wave_dense_tangent.json has no new_frozen_cache")
    points = np.asarray(cache["points"], dtype=float)
    weights = np.asarray(cache["weights"], dtype=float)
    if points.shape != (18720, 3) or weights.shape != (18720,):
        raise ValueError(f"Expected independent grid (18720,3)/(18720,), got {points.shape}/{weights.shape}")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("Independent grid contains invalid points or weights")
    return grid_report, points, weights


def _build_latest_field(drift_report):
    if drift_report.get("status") != "completed" or not drift_report.get("selected", {}).get("feasible", False):
        raise ValueError("Latest localized drift report is not completed and feasible")
    selected = drift_report["selected"]
    if len(selected.get("patch1_coefficients", [])) != 180 or len(selected.get("patch2_coefficients", [])) != 180:
        raise ValueError("Expected two selected 180-control local patches")
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    base, snapshot = build_field(candidate)
    install_in_field(base.mean if hasattr(base, "mean") else base)
    carrier = np.asarray(snapshot["inputs"]["wave"]["carrier"], dtype=float)
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    field = base
    patch_metadata = []
    for number, word in ((1, "first"), (2, "second")):
        control = np.asarray(selected[f"patch{number}_coefficients"], dtype=float)
        derivatives, pressures = _unpack_full(control, 9)
        field = LocalPotentialField(
            field,
            drift_report["inputs"][f"{word}_patch_center"],
            drift_report["inputs"][f"{word}_patch_widths"],
            2,
            {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier},
            tuple(np.zeros(27, dtype=complex) for _ in range(3)),
            derivatives,
            pressures,
            tau0,
        )
        patch_metadata.append({
            "patch": number,
            "center": list(drift_report["inputs"][f"{word}_patch_center"]),
            "widths": list(drift_report["inputs"][f"{word}_patch_widths"]),
            "control_count": int(len(control)),
        })
    field.nu = float(snapshot["inputs"]["mean"]["nu"])
    return field, snapshot, patch_metadata


def _find_inner_h(field, default=0.005):
    """Find the registered similarity exponent without assuming wrapper layout."""

    pending = [field]
    seen = set()
    while pending:
        node = pending.pop()
        if node is None or id(node) in seen:
            continue
        seen.add(id(node))
        value = getattr(node, "h", None)
        if value is not None:
            try:
                value = float(value)
                if np.isfinite(value):
                    return value
            except (TypeError, ValueError):
                pass
        for name in ("base", "mean", "inner", "reference", "current"):
            if hasattr(node, name):
                pending.append(getattr(node, name))
    return float(default)


def _jets_with_fitted_ut(field, points, tau, hspace, htime):
    """Exact stencil convention used by affine_momentum.jets, split by output."""

    points = np.asarray(points, dtype=float)
    velocity, pressure = field.fields(points, tau)
    velocity = np.asarray(velocity, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    gradient = np.empty((len(points), 3, 3), dtype=float)
    pressure_gradient = np.empty((len(points), 3), dtype=float)
    laplacian = np.zeros_like(velocity)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = hspace
        um2, pm2 = field.fields(points - 2.0 * offset, tau)
        um, pm = field.fields(points - offset, tau)
        up, pp = field.fields(points + offset, tau)
        up2, pp2 = field.fields(points + 2.0 * offset, tau)
        gradient[:, :, axis] = (um2 - 8.0 * um + 8.0 * up - up2) / (12.0 * hspace)
        pressure_gradient[:, axis] = (pm2 - 8.0 * pm + 8.0 * pp - pp2) / (12.0 * hspace)
        laplacian += (-up2 + 16.0 * up - 30.0 * velocity + 16.0 * um - um2) / (12.0 * hspace * hspace)
    # affine_momentum.jets uses this minus sign because physical time is -tau.
    ut_fitted = -(
        field.fields(points, tau - 2.0 * htime)[0]
        - 8.0 * field.fields(points, tau - htime)[0]
        + 8.0 * field.fields(points, tau + htime)[0]
        - field.fields(points, tau + 2.0 * htime)[0]
    ) / (12.0 * htime)
    return velocity, gradient, pressure_gradient, laplacian, np.asarray(ut_fitted, dtype=float)


def run(output_path=OUTPUT_PATH, chunk_size=CHUNK_SIZE):
    started = time.perf_counter()
    output_path = Path(output_path)
    drift_raw = DRIFT_PATH.read_bytes()
    drift_report = json.loads(drift_raw)
    grid_report, points, weights = _load_grid()
    frozen_report = json.loads(FROZEN_PATH.read_text(encoding="utf-8"))
    tau0 = float(frozen_report["inputs"]["mean"]["tau"])
    nu = float(frozen_report["inputs"]["mean"]["nu"])
    hspace = float(frozen_report["timesteps"]["hspace"])
    htime = float(frozen_report["timesteps"]["htime"])
    if tau0 != float(drift_report.get("inputs", {}).get("reference_tau", tau0)) and "reference_tau" in drift_report.get("inputs", {}):
        raise ValueError("Drift report reference tau disagrees with frozen geometry")
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Reference-time actual NS momentum diagnostic on the independent 18,720-point grid. "
            "The latest localized drift field is evaluated directly. The fitted physical-time "
            "derivative is replaced by the solenoidal scale generator while pressure, pressure "
            "gradient, convection, and viscosity remain unchanged. No forcing, PDE, trajectory, "
            "or recursive-scale acceptance is claimed."
        ),
        "sources": {
            "localized_drift_report": {"path": DRIFT_PATH.name, "sha256": _sha256(DRIFT_PATH)},
            "localized_drift_source": {"path": DRIFT_SOURCE_PATH.name, "sha256": _sha256(DRIFT_SOURCE_PATH)},
            "grid_report": {"path": GRID_PATH.name, "sha256": _sha256(GRID_PATH)},
            "grid_array_sha256": _array_sha256(points, weights),
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "generator": {"path": GENERATOR_PATH.name, "sha256": _sha256(GENERATOR_PATH)},
            "cache": {"path": CACHE_PATH.name},
            "field_replay": {"path": REPLAY_PATH.name, "sha256": _sha256(REPLAY_PATH)},
            "mean_report": {"path": MEAN_PATH.name, "sha256": _sha256(MEAN_PATH)},
            "candidate_report": {"path": CANDIDATE_PATH.name, "sha256": _sha256(CANDIDATE_PATH)},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "tau0": tau0,
            "k0": float(frozen_report["inputs"]["mean"]["k"]),
            "physical_time": "t=-tau",
            "viscosity": nu,
            "forcing": [0.0, 0.0, 0.0],
            "hspace": hspace,
            "htime": htime,
            "spatial_stencil": "five-point fourth-order Cartesian derivative",
            "time_stencil": "-(f(tau-2h)-8f(tau-h)+8f(tau+h)-f(tau+2h))/(12h)",
            "chunk_size": int(chunk_size),
            "mean_angles": MEAN_ANGLES,
            "cached_residual_used": False,
        },
        "field_provenance": {
            "latest_drift_selected_reason": drift_report["selected"].get("reason"),
            "base_candidate": CANDIDATE_PATH.name,
            "patches": [
                {
                    "patch": int(number),
                    "center": drift_report["inputs"][f"{word}_patch_center"],
                    "widths": drift_report["inputs"][f"{word}_patch_widths"],
                    "control_count": 180,
                }
                for number, word in ((1, "first"), (2, "second"))
            ],
        },
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points), "grid_sha256": report["sources"]["grid_array_sha256"]}), flush=True)

    field, snapshot, patch_metadata = _build_latest_field(drift_report)
    h = _find_inner_h(field)
    report["field_provenance"]["loaded_patch_metadata"] = patch_metadata
    report["inputs"]["similarity_exponent_h"] = h
    report["status"] = "fields_ready"
    _save(output_path, report)
    print(json.dumps({"stage": "fields_ready", "patches": patch_metadata}), flush=True)

    if int(chunk_size) <= 0:
        raise ValueError("chunk_size must be positive")
    reference = VelocityAdapter(field, tau0)
    base_residual_chunks = []
    generator_residual_chunks = []
    fitted_ut_chunks = []
    generator_ut_chunks = []
    delta_ut_chunks = []
    velocity_chunks = []
    gradient_chunks = []
    pressure_gradient_chunks = []
    laplacian_chunks = []
    for start in range(0, len(points), int(chunk_size)):
        stop = min(start + int(chunk_size), len(points))
        block_points = points[start:stop]
        block_weights = weights[start:stop]
        velocity, gradient, pressure_gradient, laplacian, ut_fitted = _jets_with_fitted_ut(
            field, block_points, tau0, hspace, htime
        )
        mean_swirl = axisymmetric_swirl_mean(reference, block_points, MEAN_ANGLES)
        ut_generator = solenoidal_scale_generator(
            block_points, velocity, gradient, tau0, h, mean_swirl
        )
        # The latest field uses the same registered h as the frozen mean. It is
        # stored in the localized report when available; fallback is the value
        # used by the validated generator artifact.
        convection = np.einsum("nij,nj->ni", gradient, velocity)
        shared_spatial = pressure_gradient - nu * laplacian + convection
        base_residual_chunks.append(ut_fitted + shared_spatial)
        generator_residual_chunks.append(ut_generator + shared_spatial)
        fitted_ut_chunks.append(ut_fitted)
        generator_ut_chunks.append(ut_generator)
        delta_ut_chunks.append(ut_generator - ut_fitted)
        velocity_chunks.append(velocity)
        gradient_chunks.append(gradient)
        pressure_gradient_chunks.append(pressure_gradient)
        laplacian_chunks.append(laplacian)
        print(json.dumps({"stage": "chunk", "start": start, "stop": stop, "chunks_done": (stop + int(chunk_size) - 1) // int(chunk_size)}), flush=True)

    base_residual = np.concatenate(base_residual_chunks, axis=0)
    generator_residual = np.concatenate(generator_residual_chunks, axis=0)
    fitted_ut = np.concatenate(fitted_ut_chunks, axis=0)
    generator_ut = np.concatenate(generator_ut_chunks, axis=0)
    delta_ut = np.concatenate(delta_ut_chunks, axis=0)
    velocity = np.concatenate(velocity_chunks, axis=0)
    gradient = np.concatenate(gradient_chunks, axis=0)
    pressure_gradient = np.concatenate(pressure_gradient_chunks, axis=0)
    laplacian = np.concatenate(laplacian_chunks, axis=0)
    np.savez_compressed(
        CACHE_PATH,
        points=points,
        weights=weights,
        velocity=velocity,
        gradient=gradient,
        pressure_gradient=pressure_gradient,
        laplacian=laplacian,
        fitted_ut=fitted_ut,
        generator_ut=generator_ut,
        delta_ut=delta_ut,
        original_residual=base_residual,
        generator_residual=generator_residual,
    )
    report["cache"] = {
        "path": CACHE_PATH.name,
        "sha256": _sha256(CACHE_PATH),
        "keys": [
            "points", "weights", "velocity", "gradient", "pressure_gradient",
            "laplacian", "fitted_ut", "generator_ut", "delta_ut",
            "original_residual", "generator_residual",
        ],
        "array_shapes": {
            "points": list(points.shape),
            "weights": list(weights.shape),
            "velocity": list(velocity.shape),
            "gradient": list(gradient.shape),
            "pressure_gradient": list(pressure_gradient.shape),
            "laplacian": list(laplacian.shape),
            "fitted_ut": list(fitted_ut.shape),
            "generator_ut": list(generator_ut.shape),
            "delta_ut": list(delta_ut.shape),
            "original_residual": list(base_residual.shape),
            "generator_residual": list(generator_residual.shape),
        },
    }
    report["fitted_ut"] = _metric(fitted_ut, weights, points)
    report["generator_ut"] = _metric(generator_ut, weights, points)
    report["delta_ut_generator_minus_fitted"] = _metric(delta_ut, weights, points)
    report["residual_original_fitted_ut"] = _metric(base_residual, weights, points)
    report["residual_replaced_generator_ut"] = _metric(generator_residual, weights, points)
    report["comparison"] = {
        "residual_volume_L2_change": float(
            report["residual_replaced_generator_ut"]["volume_L2"]
            - report["residual_original_fitted_ut"]["volume_L2"]
        ),
        "residual_volume_L2_relative_change": float(
            (report["residual_replaced_generator_ut"]["volume_L2"] - report["residual_original_fitted_ut"]["volume_L2"])
            / max(report["residual_original_fitted_ut"]["volume_L2"], 1.0e-300)
        ),
        "residual_max_change": float(
            report["residual_replaced_generator_ut"]["max_norm"]
            - report["residual_original_fitted_ut"]["max_norm"]
        ),
        "interpretation": (
            "Only u_t was replaced. A positive L2 change means the scale-compatible generator "
            "worsens this reference-time residual for the frozen pressure and spatial field; "
            "a negative change means improvement."
        ),
    }
    report["pressure_and_forcing"] = {
        "pressure_source": "field.fields(..., tau0) from latest reconstructed candidate, differentiated spatially by the same five-point stencil",
        "pressure_replaced": False,
        "forcing": [0.0, 0.0, 0.0],
        "viscous_term": "-nu*laplacian(u), unchanged between residuals",
        "convection": "(u dot grad)u, unchanged between residuals",
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({
        "stage": "completed",
        "original_l2": report["residual_original_fitted_ut"]["volume_L2"],
        "generator_l2": report["residual_replaced_generator_ut"]["volume_L2"],
        "delta_ut_l2": report["delta_ut_generator_minus_fitted"]["volume_L2"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    args = parser.parse_args()
    run(args.output, args.chunk_size)
