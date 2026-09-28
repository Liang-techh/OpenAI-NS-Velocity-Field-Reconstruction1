"""Cached multiscale momentum replay for the frozen corrected reference.

The completed refined replay stores velocity, Cartesian Jacobian, Laplacian,
pressure gradient, mean swirl, generator, and residual on 44,400 quadrature
points.  This diagnostic reuses those jets and adds only the axial second
derivative needed to anisotropically map the Laplacian.  It evaluates the
full momentum residual at scales ``s = 1, 1/2, 1/4`` on the correspondingly
mapped physical-volume quadrature.

The pressure rule is the explicit kinematic choice ``p_s=s**-1 p_0``.  This
is a scale-map diagnostic, not a claim that pressure symmetry closes the
Navier--Stokes equation.  No PDE, trajectory, or recursive-scale acceptance
is inferred from these residuals.
"""

from __future__ import annotations

import argparse
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

import scale_reference_candidate as candidate_module  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from scale_transport_generator import solenoidal_scale_generator  # noqa: E402


ACTUAL_REPLAY_PATH = ROOT / "scale_reference_trust_nonlinear_actual.json"
OUTPUT_PATH = ROOT / "scale_reference_multiscale.json"
MEAN_ANGLES = 12
SCALES = (1.0, 0.5, 0.25)
HSPACE_DEFAULT = 1.0e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _save(path: Path, report: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _metric_state() -> dict:
    return {
        "weighted_square": 0.0,
        "weight_sum": 0.0,
        "processed": 0,
        "max_norm": None,
        "max_point": None,
        "max_index": -1,
        "max_vector": None,
    }


def _metric_update(state: dict, values: np.ndarray, weights: np.ndarray,
                   points: np.ndarray, offset: int) -> None:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    state["weighted_square"] += float(np.sum(weights * np.sum(values * values, axis=1)))
    state["weight_sum"] += float(np.sum(weights))
    index = int(np.argmax(norms))
    if state["max_norm"] is None or float(norms[index]) > float(state["max_norm"]):
        state["max_norm"] = float(norms[index])
        state["max_point"] = points[index].tolist()
        state["max_index"] = int(offset + index)
        state["max_vector"] = values[index].tolist()
    state["processed"] += int(len(values))


def _metric_finish(state: dict, count: int, volume: float) -> dict:
    square = max(float(state["weighted_square"]), 0.0)
    return {
        "point_count": int(count),
        "volume_L2": float(np.sqrt(square)),
        "volume_RMS": float(np.sqrt(square / max(float(volume), 1.0e-300))),
        "max_norm": float(state["max_norm"]),
        "max_point": state["max_point"],
        "max_index": int(state["max_index"]),
        "max_vector": state["max_vector"],
        "physical_volume": float(volume),
    }


def _load_replay_report(path: Path) -> tuple[dict, Path, list[Path], float]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "completed":
        raise ValueError(f"Replay report is not completed: {report.get('status')}")
    actual_case = report.get("cases", {}).get("corrected_reference", {})
    if actual_case.get("status") != "completed":
        raise ValueError("Corrected reference case is not completed")
    cache = report.get("chunk_caches", {}).get("corrected_reference", {})
    cache_dir = Path(cache.get("directory", ""))
    if not cache_dir.exists():
        raise FileNotFoundError(f"Missing replay cache directory: {cache_dir}")
    chunk_paths = sorted(cache_dir.glob("*.npz"), key=lambda value: int(value.stem))
    expected = int(report["grid"]["point_count"])
    if len(chunk_paths) != int(actual_case.get("chunks_completed", 0)):
        raise ValueError("Replay cache chunk count does not match completed case")
    if len(chunk_paths) == 0:
        raise ValueError("Replay cache contains no chunks")
    hspace = float(report.get("inputs", {}).get("hspace", HSPACE_DEFAULT))
    return report, cache_dir, chunk_paths, hspace


def _axial_second_derivative(reference, points: np.ndarray, center_velocity: np.ndarray,
                             hspace: float) -> np.ndarray:
    offset = np.zeros(3, dtype=float)
    offset[2] = float(hspace)
    stencil_points = np.concatenate((
        points - 2.0 * offset,
        points - offset,
        points + offset,
        points + 2.0 * offset,
    ), axis=0)
    values = np.asarray(reference.velocity(stencil_points), dtype=float)
    n = len(points)
    minus2, minus1, plus1, plus2 = np.split(values, 4)
    return (-plus2 + 16.0 * plus1 - 30.0 * center_velocity + 16.0 * minus1 - minus2) / (12.0 * hspace * hspace)


def _ring_means(points: np.ndarray, velocity: np.ndarray, jacobian: np.ndarray,
                laplacian: np.ndarray, hzz: np.ndarray) -> tuple[np.ndarray, ...]:
    count = len(points)
    if count % MEAN_ANGLES != 0:
        raise ValueError(f"Chunk count {count} is not divisible by {MEAN_ANGLES}")
    radius = np.linalg.norm(points[:, :2], axis=1)
    if np.any(radius <= 0.0):
        raise ValueError("Axis points are not supported by cylindrical ring derivatives")
    e_r = np.column_stack((points[:, 0] / radius, points[:, 1] / radius, np.zeros(count)))
    e_theta = np.column_stack((-points[:, 1] / radius, points[:, 0] / radius, np.zeros(count)))
    ring_radius = radius.reshape(-1, MEAN_ANGLES)
    ring_z = points[:, 2].reshape(-1, MEAN_ANGLES)
    if float(np.max(np.ptp(ring_radius, axis=1))) > 2.0e-12 or float(np.max(np.ptp(ring_z, axis=1))) > 2.0e-12:
        raise ValueError("Cache ordering does not contain contiguous angular rings")
    swirl = np.einsum("ni,ni->n", e_theta, velocity)
    swirl_r = np.einsum("ni,nij,nj->n", e_theta, jacobian, e_r)
    swirl_z = np.einsum("ni,nij,nj->n", e_theta, jacobian, np.broadcast_to([0.0, 0.0, 1.0], (count, 3)))
    theta_laplacian = np.einsum("ni,ni->n", e_theta, laplacian)
    theta_hzz = np.einsum("ni,ni->n", e_theta, hzz)
    def broadcast_ring(values):
        return np.repeat(np.mean(np.asarray(values).reshape(-1, MEAN_ANGLES), axis=1), MEAN_ANGLES)
    return (
        broadcast_ring(swirl),
        broadcast_ring(swirl_r),
        broadcast_ring(swirl_z),
        broadcast_ring(theta_laplacian),
        broadcast_ring(theta_hzz),
    )


def _mapped_residual(points: np.ndarray, velocity: np.ndarray, jacobian: np.ndarray,
                     laplacian: np.ndarray, pressure_gradient: np.ndarray,
                     mean_swirl: np.ndarray, hzz: np.ndarray, scale: float,
                     tau0: float, h: float, nu: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    scale = float(scale)
    sr = scale ** 0.5
    sz = scale ** (0.5 - h)
    det_s = scale ** (1.5 - h)
    piola_amplitude = scale ** (0.5 - h)
    Bdiag = piola_amplitude * np.array([sr, sr, sz], dtype=float) / det_s
    inv_s = np.array([1.0 / sr, 1.0 / sr, 1.0 / sz], dtype=float)
    swirl_increment = scale ** (-0.5 - h) - scale ** (-0.5)
    radius = np.linalg.norm(points[:, :2], axis=1)
    if np.any(radius <= 0.0):
        raise ValueError("Mapped replay does not support axis points")
    e_r = np.column_stack((points[:, 0] / radius, points[:, 1] / radius, np.zeros(len(points))))
    e_theta = np.column_stack((-points[:, 1] / radius, points[:, 0] / radius, np.zeros(len(points))))
    target_points = points * np.array([sr, sr, sz], dtype=float)[None, :]
    target_radius = sr * radius
    mean, mean_r, mean_z, theta_laplacian, theta_hzz = _ring_means(
        points, velocity, jacobian, laplacian, hzz
    )
    # The cache is ordered as output-row/component-column Cartesian jets.
    mapped_velocity = velocity * Bdiag[None, :] + swirl_increment * mean[:, None] * e_theta
    mapped_jacobian = Bdiag[None, :, None] * jacobian * inv_s[None, None, :]
    grad_mean_target = np.column_stack((
        mean_r * e_r[:, 0] * inv_s[0],
        mean_r * e_r[:, 1] * inv_s[1],
        mean_z * inv_s[2],
    ))
    mapped_jacobian += swirl_increment * (
        np.einsum("ni,nj->nij", e_theta, grad_mean_target)
        - (mean / target_radius)[:, None, None] * np.einsum("ni,nj->nij", e_r, e_theta)
    )
    lap_xy = laplacian - hzz
    mapped_laplacian = Bdiag[None, :] * (lap_xy * (inv_s[None, :] ** 2) + hzz * (inv_s[2] ** 2))
    mean_lap_xy = theta_laplacian - theta_hzz
    mapped_laplacian += swirl_increment * e_theta * (
        mean_lap_xy * (inv_s[0] ** 2) + theta_hzz * (inv_s[2] ** 2)
    )[:, None]
    mapped_pressure_gradient = scale ** -1.0 * pressure_gradient * inv_s[None, :]
    generator = solenoidal_scale_generator(
        target_points,
        mapped_velocity,
        mapped_jacobian,
        tau0 * scale,
        h,
        scale ** (-0.5 - h) * mean,
    )
    convection = np.einsum("nij,nj->ni", mapped_jacobian, mapped_velocity)
    residual = generator + convection + mapped_pressure_gradient - nu * mapped_laplacian
    return target_points, residual, mapped_jacobian


def run(output_path: Path | str = OUTPUT_PATH,
        replay_path: Path | str = ACTUAL_REPLAY_PATH) -> dict:
    started = time.perf_counter()
    output_path = Path(output_path)
    replay_path = Path(replay_path)
    replay, cache_dir, chunk_paths, hspace = _load_replay_report(replay_path)
    step_name = replay["sources"]["step"]["path"]
    step_path = ROOT / step_name
    if not step_path.exists():
        raise FileNotFoundError(f"Missing replay step source: {step_path}")
    candidate = candidate_module.load_reference(step=step_path)
    try:
        install_in_field(candidate.base_field)
        grouped_backend = "installed"
    except Exception as exc:  # pragma: no cover - fallback is intentional
        grouped_backend = f"unavailable: {type(exc).__name__}: {exc}"
    tau0 = float(candidate.tau0)
    h = float(candidate.h)
    nu = float(candidate.nu)
    grid_count = int(replay["grid"]["point_count"])
    states = {str(scale): _metric_state() for scale in SCALES}
    divergence_states = {str(scale): _metric_state() for scale in SCALES}
    agreement = _metric_state()
    agreement_max_abs = 0.0
    processed = 0
    expected_begin = 0
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Cached multiscale momentum diagnostic for the frozen corrected reference. "
            "The saved Cartesian jets are mapped by the anisotropic Piola/swirl rule "
            "at scales 1, 1/2, and 1/4; only axial Hzz is recomputed with a fourth-order "
            "four-point stencil. The pressure rule is the explicit kinematic choice "
            "p_s=s^-1 p_0. This does not establish Navier--Stokes symmetry, PDE closure, "
            "finite-time evolution, or scale recursion."
        ),
        "sources": {
            "actual_replay": {"path": replay_path.name, "sha256": _sha256(replay_path)},
            "candidate_loader": {"path": replay["sources"]["candidate_loader"]["path"], "sha256": _sha256(ROOT / replay["sources"]["candidate_loader"]["path"])},
            "step": {"path": step_name, "sha256": _sha256(step_path)},
            "generator": {"path": replay["sources"]["generator"]["path"], "sha256": _sha256(ROOT / replay["sources"]["generator"]["path"])},
            "grid": dict(replay["sources"]["grid"]),
            "chunk_cache": {"directory": str(cache_dir), "fingerprint": replay["chunk_caches"]["corrected_reference"]["fingerprint"], "chunk_count": len(chunk_paths)},
        },
        "inputs": {
            "point_count": grid_count,
            "tau0": tau0,
            "similarity_exponent_h": h,
            "viscosity": nu,
            "physical_time_convention": "t=-tau",
            "scales": list(SCALES),
            "mean_angles": MEAN_ANGLES,
            "axial_stencil_step": hspace,
            "axial_stencil": "four-point fourth-order Hzz using +/-h and +/-2h; center velocity reused from cache",
            "pressure_rule": "p_s=s^-1 p_0; grad_p_s=s^-1 S^-1 grad_p_0",
            "pressure_interpretation": "explicit kinematic pressure scaling; not asserted as Navier--Stokes pressure symmetry",
            "quadrature_rule": "saved points/weights mapped by det(S); no domain narrowing",
            "grouped_backend": grouped_backend,
        },
        "reference_cache_metric": replay["cases"]["corrected_reference"]["metric"],
        "scales": {
            str(scale): {
                "status": "running",
                "spatial_matrix": [scale ** 0.5, scale ** 0.5, scale ** (0.5 - h)],
                "determinant": scale ** (1.5 - h),
                "piola_amplitude": scale ** (0.5 - h),
                "swirl_increment_factor": scale ** (-0.5 - h) - scale ** (-0.5),
                "physical_volume": float(replay["grid"]["physical_volume"] * scale ** (1.5 - h)),
            }
            for scale in SCALES
        },
        "cache_agreement": {"status": "running"},
        "chunk_progress": {"processed": 0, "total": grid_count},
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "chunks": len(chunk_paths), "point_count": grid_count}, sort_keys=True), flush=True)

    for chunk_path in chunk_paths:
        begin = int(chunk_path.stem)
        if begin != expected_begin:
            raise ValueError(f"Unexpected chunk offset {begin}; expected {expected_begin}")
        with np.load(chunk_path, allow_pickle=False) as cached:
            required = {"points", "weights", "velocity", "jacobian", "pressure_gradient", "laplacian", "mean_swirl", "generator", "residual"}
            missing = sorted(required.difference(cached.files))
            if missing:
                raise ValueError(f"Chunk {chunk_path.name} missing {missing}")
            points = np.asarray(cached["points"], dtype=float)
            weights = np.asarray(cached["weights"], dtype=float)
            velocity = np.asarray(cached["velocity"], dtype=float)
            jacobian = np.asarray(cached["jacobian"], dtype=float)
            pressure_gradient = np.asarray(cached["pressure_gradient"], dtype=float)
            laplacian = np.asarray(cached["laplacian"], dtype=float)
            cached_mean_swirl = np.asarray(cached["mean_swirl"], dtype=float)
            cached_residual = np.asarray(cached["residual"], dtype=float)
        if len(points) == 0 or len(points) != len(weights) or velocity.shape != (len(points), 3) or jacobian.shape != (len(points), 3, 3):
            raise ValueError(f"Invalid chunk shapes in {chunk_path.name}")
        hzz = _axial_second_derivative(candidate, points, velocity, hspace)
        for scale in SCALES:
            key = str(scale)
            target_points, residual, mapped_jacobian = _mapped_residual(
                points, velocity, jacobian, laplacian, pressure_gradient,
                cached_mean_swirl, hzz, scale, tau0, h, nu,
            )
            target_weights = weights * scale ** (1.5 - h)
            _metric_update(states[key], residual, target_weights, target_points, begin)
            divergence = np.trace(mapped_jacobian, axis1=1, axis2=2)
            _metric_update(
                divergence_states[key],
                np.column_stack((divergence, np.zeros((len(divergence), 2)))),
                target_weights,
                target_points,
                begin,
            )
            if scale == 1.0:
                difference = residual - cached_residual
                _metric_update(agreement, difference, weights, points, begin)
                agreement_max_abs = max(agreement_max_abs, float(np.max(np.abs(difference))))
        processed += len(points)
        expected_begin += len(points)
        report["chunk_progress"]["processed"] = int(processed)
        if processed == grid_count or processed % (768 * 8) == 0:
            _save(output_path, report)
            print(json.dumps({"stage": "chunk", "processed": processed, "total": grid_count}, sort_keys=True), flush=True)
    if processed != grid_count:
        raise ValueError(f"Only processed {processed} of {grid_count} points")

    for scale in SCALES:
        key = str(scale)
        report["scales"][key].update({
            "status": "completed",
            "residual": _metric_finish(states[key], grid_count, states[key]["weight_sum"]),
            "divergence": {
                "max_abs": float(divergence_states[key]["max_norm"]),
                "weighted_RMS": float(np.sqrt(divergence_states[key]["weighted_square"] / max(divergence_states[key]["weight_sum"], 1.0e-300))),
                "max_point": divergence_states[key]["max_point"],
                "max_index": int(divergence_states[key]["max_index"]),
            },
        })
    cached_metric = replay["cases"]["corrected_reference"]["metric"]
    agreement_metric = _metric_finish(agreement, grid_count, agreement["weight_sum"])
    report["cache_agreement"] = {
        "status": "completed",
        "scale_one_metric": report["scales"]["1.0"]["residual"],
        "cached_metric": cached_metric,
        "difference_metric": agreement_metric,
        "difference_max_abs_component": agreement_max_abs,
        "relative_volume_L2": float(agreement_metric["volume_L2"] / max(float(cached_metric["volume_L2"]), 1.0e-300)),
        "criterion": "scale 1 uses the same mapped points and weights; difference is from 12-angle mean swirl and mapped jet reconstruction",
    }
    report["status"] = "completed"
    report["status_detail"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({"stage": "completed", "elapsed_seconds": report["elapsed_seconds"], "scales": {key: report["scales"][key]["residual"]["volume_L2"] for key in report["scales"]}, "s1_difference_l2": agreement_metric["volume_L2"], "s1_difference_max": agreement_max_abs}, sort_keys=True), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--replay", type=Path, default=ACTUAL_REPLAY_PATH)
    args = parser.parse_args()
    try:
        run(args.output, args.replay)
    except Exception as exc:
        report = {}
        if args.output.exists():
            try:
                report = json.loads(args.output.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                report = {}
        report.update(status="failed", accepted=False, pde_validated=False, scale_recursion_established=False, error=f"{type(exc).__name__}: {exc}")
        _save(args.output, report)
        raise
