"""Independent reference-time replay for the frozen scale-generator step.

This diagnostic deliberately uses the 44,400-point refined grid, rather than
the 18,720 points used by the reference-velocity fit.  Both the uncorrected
reference (with the earlier pressure projection) and the corrected reference
are differentiated with a five-point Cartesian spatial stencil.  Their time
derivative is the analytic solenoidal scale generator, so this file does not
use a fitted temporal derivative.

The loader supplies the initial curl correction. Pressure composition is
checked against independently decoded prior and selected coefficients before
the grid replay. Each completed chunk is saved atomically with its numerical
jets and source fingerprint; failure reporting preserves completed cases.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Callable

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scale_reference_candidate as candidate_module  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from scale_transport_generator import (  # noqa: E402
    MEAN_ANGLES,
    axisymmetric_swirl_mean,
    solenoidal_scale_generator,
)
from supported_fourier_basis import basis_data  # noqa: E402


CANDIDATE_PATH = ROOT / "scale_reference_candidate.py"
STEP_PATH = ROOT / "scale_reference_velocity_step.json"
PRESSURE_PATH = ROOT / "scale_generator_pressure_projection.json"
GRID_PATH = ROOT / "refined_wave_momentum_cache.npz"
OUTPUT_PATH = ROOT / "scale_reference_actual_replay.json"
RECOVERED_BASELINE_PATH = ROOT / "scale_reference_baseline_recovered.json"
HSPACE_DEFAULT = 1.0e-6
CHUNK_SIZE = 768


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays: Any) -> str:
    digest = hashlib.sha256()
    for value in arrays:
        array = np.ascontiguousarray(np.asarray(value))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _save(path: Path, report: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _metric_update(state: dict[str, Any], values: np.ndarray,
                   weights: np.ndarray, points: np.ndarray) -> None:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    state["weighted_square"] += float(np.sum(weights * np.sum(values * values, axis=1)))
    state["weight_sum"] += float(np.sum(weights))
    index = int(np.argmax(norms))
    if state["max_norm"] is None or float(norms[index]) > state["max_norm"]:
        state["max_norm"] = float(norms[index])
        state["max_point"] = points[index].tolist()
        state["max_index"] = int(state["processed"] + index)
        state["max_vector"] = values[index].tolist()
    state["processed"] += int(len(values))


def _metric_finish(state: dict[str, Any], count: int, volume: float) -> dict[str, Any]:
    square = max(float(state["weighted_square"]), 0.0)
    return {
        "point_count": int(count),
        "volume_L2": float(np.sqrt(square)),
        "volume_RMS": float(np.sqrt(square / max(volume, 1.0e-300))),
        "max_norm": float(state["max_norm"]),
        "max_point": state["max_point"],
        "max_index": int(state["max_index"]),
        "max_vector": state["max_vector"],
        "physical_volume": float(volume),
    }


def _metric_state() -> dict[str, Any]:
    return {
        "weighted_square": 0.0,
        "weight_sum": 0.0,
        "processed": 0,
        "max_norm": None,
        "max_point": None,
        "max_index": -1,
        "max_vector": None,
    }


class _FrozenView:
    """Add a scalar pressure correction to an existing frozen field."""

    def __init__(self, inner: Any, tau0: float,
                 pressure_correction: Callable[[np.ndarray], np.ndarray] | None,
                 label: str):
        self.inner = inner
        self.tau0 = float(tau0)
        self.pressure_correction = pressure_correction
        self.label = label
        self.h = float(getattr(inner, "h", 0.005))
        self.nu = float(getattr(inner, "nu", 0.01))

    def fields(self, points: Any, tau: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        if tau is not None and not np.isclose(float(tau), self.tau0, rtol=0.0, atol=2.0e-15):
            raise ValueError(f"{self.label} is frozen at tau0")
        array = np.asarray(points, dtype=float)
        velocity, pressure = self.inner.fields(array, self.tau0)
        velocity = np.asarray(velocity, dtype=float)
        pressure = np.asarray(pressure, dtype=float)
        if self.pressure_correction is not None:
            pressure = pressure + self.pressure_correction(array)
        return velocity, pressure

    def velocity(self, points: Any) -> np.ndarray:
        return self.fields(points, self.tau0)[0]


def _decode_pressure_block(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    block = np.asarray(block, dtype=float)
    if block.shape != (45,):
        raise ValueError(f"pressure block must have shape (45,), got {block.shape}")
    q = 9
    cursor = 0
    modes: list[np.ndarray] = []
    modes.append(block[cursor:cursor + q].astype(complex))
    cursor += q
    for _mode in (1, 2):
        modes.append(block[cursor:cursor + 2 * q:2] + 1j * block[cursor + 1:cursor + 2 * q:2])
        cursor += 2 * q
    if cursor != len(block):
        raise ValueError("pressure layout decode failed")
    return tuple(modes)  # type: ignore[return-value]


def _load_pressure_projection(reference: Any) -> tuple[Callable[[np.ndarray], np.ndarray], dict[str, Any]]:
    report = json.loads(PRESSURE_PATH.read_text(encoding="utf-8"))
    if report.get("status") != "completed":
        raise ValueError("Pressure projection report is not completed")
    fit = report.get("fit", {})
    coefficients = np.asarray(fit.get("coefficients"), dtype=float)
    if coefficients.shape != (90,):
        raise ValueError(f"Pressure projection coefficients have shape {coefficients.shape}, expected (90,)")
    inputs = report["inputs"]
    centers = (
        np.asarray(inputs["first_patch_center"], dtype=float),
        np.asarray(inputs["second_patch_center"], dtype=float),
    )
    widths = (
        np.asarray(inputs["first_patch_widths"], dtype=float),
        np.asarray(inputs["second_patch_widths"], dtype=float),
    )
    carriers = tuple(np.asarray(value, dtype=float) for value in (
        inputs["carrier"], inputs["carrier"],
    ))
    modes = tuple(
        tuple(_decode_pressure_block(coefficients[45 * index:45 * (index + 1)]))
        for index in range(2)
    )

    def correction(points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        result = np.zeros(len(points), dtype=float)
        for patch in range(2):
            for mode in (0, 1, 2):
                carrier = np.zeros(2, dtype=float) if mode == 0 else mode * carriers[patch]
                scalar = basis_data(points, centers[patch], widths[patch], mode, 2, carrier)[1]
                result += np.einsum("nq,q->n", scalar, modes[patch][mode]).real
        return result

    metadata = {
        "path": PRESSURE_PATH.name,
        "sha256": _sha256(PRESSURE_PATH),
        "coefficient_count": 90,
        "centers": [value.tolist() for value in centers],
        "widths": [value.tolist() for value in widths],
        "carrier": np.asarray(inputs["carrier"], dtype=float).tolist(),
    }
    return correction, metadata


def _pressure_from_step_increment(candidate: Any, points: np.ndarray) -> np.ndarray:
    """Evaluate the joint-step pressure increment, excluding prior projection."""

    step = getattr(candidate, "step", {})
    selected = step.get("selected", {}) if isinstance(step, dict) else {}
    blocks = []
    for number in (1, 2):
        value = selected.get(f"pressure_coefficients_patch{number}")
        if value is None:
            value = step.get(f"pressure_coefficients_patch{number}") if isinstance(step, dict) else None
        if value is None:
            raise ValueError("Loader step has no selected pressure increment blocks")
        blocks.append(_decode_pressure_block(np.asarray(value, dtype=float)))
    result = np.zeros(len(points), dtype=float)
    for patch in range(2):
        for mode in (0, 1, 2):
            carrier = np.zeros(2, dtype=float) if mode == 0 else mode * np.asarray(candidate.carriers[patch][1], dtype=float)
            scalar = basis_data(points, candidate.center[patch], candidate.widths[patch], mode, 2, carrier)[1]
            result += np.einsum("nq,q->n", scalar, blocks[patch][mode]).real
    return result


def _pressure_composition(candidate: Any, pressure_projection: Callable[[np.ndarray], np.ndarray],
                          tau0: float) -> tuple[Any, Any, dict[str, Any]]:
    """Build baseline/corrected views and verify loader pressure semantics."""

    # Use a small off-support/interior mixture so the pressure differences are
    # not inferred from an identically zero compact basis value.
    sample = np.array([
        [0.001000, 0.0, 0.0],
        [0.001000, 0.0, 0.000010],
        [0.001180, 0.000020, 0.000005],
    ], dtype=float)
    raw_p = np.asarray(candidate.base_field.fields(sample, tau0)[1], dtype=float)
    corrected_p = np.asarray(candidate.fields(sample, tau0)[1], dtype=float)
    prior = np.asarray(pressure_projection(sample), dtype=float)
    selected = _pressure_from_step_increment(candidate, sample)
    base_reference = getattr(candidate, "base_reference", None)
    base_reference_delta = None
    if base_reference is not None:
        base_reference_delta = np.asarray(base_reference.fields(sample, tau0)[1], dtype=float) - raw_p

    tol = 2.0e-8 * max(1.0, float(np.max(np.abs(prior))), float(np.max(np.abs(selected))))
    corrected_delta = corrected_p - raw_p
    loader_has_prior = bool(np.max(np.abs(corrected_delta - (prior + selected))) <= tol)
    loader_has_selected_only = bool(np.max(np.abs(corrected_delta - selected)) <= tol)
    if loader_has_prior:
        # Build baseline from the immutable raw field and prior projection so
        # a loader's convenience ``base_reference`` view cannot double-count.
        baseline = _FrozenView(candidate.base_field, tau0, pressure_projection, "baseline reference")
        corrected = candidate
        composition = "loader corrected field includes prior projection"
    elif loader_has_selected_only:
        baseline = _FrozenView(candidate.base_field, tau0, pressure_projection, "baseline reference")
        corrected = _FrozenView(candidate, tau0, pressure_projection, "corrected reference")
        composition = "loader selected increment plus external prior projection"
    else:
        raise ValueError(
            "Ambiguous pressure composition: loader does not identify prior projection; "
            f"prior_delta={prior.tolist()}, selected_increment={selected.tolist()}, "
            f"base_reference_delta={None if base_reference_delta is None else base_reference_delta.tolist()}"
        )
    return baseline, corrected, {
        "composition": composition,
        "prior_projection_sample": prior.tolist(),
        "selected_pressure_sample": selected.tolist(),
        "base_reference_pressure_delta_sample": None if base_reference_delta is None else base_reference_delta.tolist(),
        "corrected_pressure_delta_sample": corrected_delta.tolist(),
        "comparison_tolerance": tol,
        "loader_has_prior_projection": loader_has_prior,
        "loader_base_reference_pressure_delta_is_prior_plus_increment": bool(
            base_reference_delta is not None
            and np.max(np.abs(base_reference_delta - (prior + selected))) <= tol
        ),
        "loader_corrected_delta_is_selected_increment_only": loader_has_selected_only,
        "step_selected_increment_sample": selected.tolist(),
    }


def _load_grid() -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    with np.load(GRID_PATH, allow_pickle=False) as loaded:
        if "points" not in loaded.files or "weights" not in loaded.files:
            raise ValueError("Refined cache must contain points and weights")
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError(f"Expected refined grid (44400,3)/(44400,), got {points.shape}/{weights.shape}")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("Refined grid contains invalid points or weights")
    return points, weights, {
        "point_count": int(len(points)),
        "physical_volume": float(np.sum(weights)),
        "array_sha256": _array_sha256(points, weights),
        "path": GRID_PATH.name,
        "sha256": _sha256(GRID_PATH),
    }


def _spatial_jets(field: Any, points: np.ndarray, tau0: float, hspace: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    velocity, pressure = field.fields(points, tau0)
    velocity = np.asarray(velocity, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    gradient = np.empty((len(points), 3, 3), dtype=float)
    pressure_gradient = np.empty((len(points), 3), dtype=float)
    laplacian = np.zeros_like(velocity)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = hspace
        stencil_points = np.concatenate((
            points - 2.0 * offset,
            points - offset,
            points + offset,
            points + 2.0 * offset,
        ), axis=0)
        stencil_velocity, stencil_pressure = field.fields(stencil_points, tau0)
        n = len(points)
        um2, um, up, up2 = np.split(np.asarray(stencil_velocity, dtype=float), 4)
        pm2, pm, pp, pp2 = np.split(np.asarray(stencil_pressure, dtype=float), 4)
        gradient[:, :, axis] = (um2 - 8.0 * um + 8.0 * up - up2) / (12.0 * hspace)
        pressure_gradient[:, axis] = (pm2 - 8.0 * pm + 8.0 * pp - pp2) / (12.0 * hspace)
        laplacian += (-up2 + 16.0 * up - 30.0 * velocity + 16.0 * um - um2) / (12.0 * hspace * hspace)
    return velocity, gradient, pressure_gradient, laplacian


def _run_case(field: Any, label: str, points: np.ndarray, weights: np.ndarray,
              tau0: float, h: float, nu: float, hspace: float, chunk_size: int,
              report: dict[str, Any], output_path: Path) -> dict[str, Any]:
    state = _metric_state()
    div_state = _metric_state()
    count = len(points)
    volume = float(np.sum(weights))
    fingerprint = hashlib.sha256(json.dumps(dict(sources=report['sources'],
        tau0=tau0,h=h,nu=nu,hspace=hspace,label=label,chunk_size=chunk_size),sort_keys=True).encode()).hexdigest()
    cache_dir = ROOT / '.scale_replay_cache' / fingerprint
    cache_dir.mkdir(parents=True,exist_ok=True)
    report.setdefault('chunk_caches',{})[label] = dict(directory=str(cache_dir),fingerprint=fingerprint)
    report["cases"][label] = {"status": "running", "chunks_completed": 0, "chunk_count": int((count + chunk_size - 1) // chunk_size)}
    _save(output_path, report)
    for begin in range(0, count, chunk_size):
        end = min(begin + chunk_size, count)
        p = points[begin:end]
        w = weights[begin:end]
        cache_path = cache_dir / f'{begin:06d}.npz'
        if cache_path.exists():
            with np.load(cache_path,allow_pickle=False) as cached:
                if not np.array_equal(cached['points'],p) or not np.array_equal(cached['weights'],w):
                    raise ValueError('Chunk cache grid mismatch')
                residual,jacobian = cached['residual'],cached['jacobian']
        else:
            velocity, jacobian, pressure_gradient, laplacian = _spatial_jets(field, p, tau0, hspace)
            mean_swirl = axisymmetric_swirl_mean(field, p, MEAN_ANGLES)
            generator = solenoidal_scale_generator(p, velocity, jacobian, tau0, h, mean_swirl)
            convection = np.einsum("nij,nj->ni", jacobian, velocity)
            residual = generator + convection + pressure_gradient - nu * laplacian
            temporary = cache_path.with_suffix('.tmp')
            with temporary.open('wb') as stream:
                np.savez_compressed(stream,points=p,weights=w,velocity=velocity,jacobian=jacobian,
                    pressure_gradient=pressure_gradient,laplacian=laplacian,mean_swirl=mean_swirl,
                    generator=generator,residual=residual)
            temporary.replace(cache_path)
        _metric_update(state, residual, w, p)
        divergence = np.trace(jacobian, axis1=1, axis2=2)
        div_values = np.column_stack((divergence, np.zeros((len(divergence), 2))))
        _metric_update(div_state, div_values, w, p)
        report["cases"][label]["chunks_completed"] = int(end // chunk_size + (1 if end == count else 0))
        report["cases"][label]["processed"] = int(end)
        if end == count or end % (chunk_size * 8) == 0:
            _save(output_path, report)
            print(json.dumps({"stage": "replay_chunk", "case": label, "processed": int(end), "total": int(count)}, sort_keys=True), flush=True)
    metric = _metric_finish(state, count, volume)
    div_metric = _metric_finish(div_state, count, volume)
    report["cases"][label] = {
        "status": "completed",
        "metric": metric,
        "divergence": {
            "max_abs": float(div_metric["max_norm"]),
            "weighted_RMS": float(div_metric["volume_RMS"]),
            "max_point": div_metric["max_point"],
            "max_index": div_metric["max_index"],
        },
        "chunks_completed": int((count + chunk_size - 1) // chunk_size),
        "processed": int(count),
    }
    _save(output_path, report)
    return report["cases"][label]


def run(output_path: Path | str = OUTPUT_PATH, chunk_size: int = CHUNK_SIZE,
        hspace: float = HSPACE_DEFAULT, step_path: Path | str = STEP_PATH,
        case_mode: str = "both") -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    step_path = Path(step_path)
    if case_mode not in ("both", "baseline", "corrected"):
        raise ValueError("case_mode must be both, baseline, or corrected")
    prior_report = None
    if case_mode == 'corrected':
        recovered = json.loads(RECOVERED_BASELINE_PATH.read_text(encoding='utf-8'))
        meta = recovered['inputs']
        if (meta['grid_sha256'] != _sha256(GRID_PATH)
            or meta['pressure_projection_sha256'] != _sha256(PRESSURE_PATH)
            or meta['raw_reference_sha256'] != _sha256(ROOT/meta['raw_reference_source'])
            or meta['hspace'] != hspace):
            raise ValueError('Recovered baseline does not match current grid, pressure, reference or stencil')
        if recovered['cases']['uncorrected_plus_pressure_projection']['status'] != 'completed':
            raise ValueError('Completed baseline is required before corrected replay')
        prior_report = recovered
    if case_mode == "corrected" and prior_report is None and output_path.exists():
        try:
            prior_report = json.loads(output_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            prior_report = None
    if case_mode == "corrected" and (
        prior_report is None
        or "uncorrected_plus_pressure_projection" not in prior_report.get("cases", {})
    ) and RECOVERED_BASELINE_PATH.exists():
        try:
            recovered = json.loads(RECOVERED_BASELINE_PATH.read_text(encoding="utf-8"))
            if "uncorrected_plus_pressure_projection" in recovered.get("cases", {}):
                if prior_report is None:
                    prior_report = {}
                prior_report["cases"] = dict(prior_report.get("cases", {}))
                prior_report["cases"]["uncorrected_plus_pressure_projection"] = recovered["cases"]["uncorrected_plus_pressure_projection"]
                prior_report["baseline_recovery"] = {
                    "path": RECOVERED_BASELINE_PATH.name,
                    "sha256": _sha256(RECOVERED_BASELINE_PATH),
                    "status": recovered.get("status"),
                }
        except (OSError, ValueError):
            pass
    report: dict[str, Any] = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Independent frozen-reference scale-generator momentum replay on the refined "
            "44,400-point grid. Five-point Cartesian spatial jets and analytic G(u) are used; "
            "no fitted temporal derivative, trajectory, PDE acceptance, or recursive-scale claim."
        ),
        "sources": {
            "candidate_loader": {"path": CANDIDATE_PATH.name, "sha256": _sha256(CANDIDATE_PATH)},
            "step": {"path": step_path.name, "sha256": _sha256(step_path)},
            "pressure_projection": {"path": PRESSURE_PATH.name, "sha256": _sha256(PRESSURE_PATH)},
            "grid": {"path": GRID_PATH.name, "sha256": _sha256(GRID_PATH)},
            "baseline_recovery": ({"path": RECOVERED_BASELINE_PATH.name, "sha256": _sha256(RECOVERED_BASELINE_PATH)}
                                   if RECOVERED_BASELINE_PATH.exists() else None),
            "generator": {"path": (ROOT / "scale_transport_generator.py").name,
                           "sha256": _sha256(ROOT / "scale_transport_generator.py")},
        },
        "inputs": {"chunk_size": int(chunk_size), "hspace": float(hspace), "spatial_stencil": "five-point fourth-order", "temporal_stencil": None, "case_mode": case_mode, "step_path": step_path.name},
        "cases": {},
        "status_detail": "loading",
    }
    if prior_report is not None:
        old_cases = prior_report.get("cases", {})
        if "uncorrected_plus_pressure_projection" in old_cases:
            report["cases"]["uncorrected_plus_pressure_projection"] = old_cases["uncorrected_plus_pressure_projection"]
            report["baseline_reused"] = {
                "source_report": output_path.name,
                "source_status": prior_report.get("status"),
                "source_elapsed_seconds": prior_report.get("elapsed_seconds"),
                "reason": "same raw reference and prior pressure projection; baseline was valid before corrected pressure-composition fix",
            }
        if "corrected_reference" in old_cases:
            report["superseded_inflight"] = {
                "case": old_cases["corrected_reference"],
                "reason": "prior run double-counted the pressure projection during composition inference",
            }
        if prior_report.get("baseline_recovery"):
            report["baseline_recovery"] = prior_report["baseline_recovery"]
    _save(output_path, report)
    if chunk_size <= 0 or not np.isfinite(hspace) or hspace <= 0.0:
        raise ValueError("chunk_size and hspace must be positive")
    candidate = candidate_module.load_reference(step=step_path)
    try:
        install_in_field(candidate.base_field)
    except Exception as exc:  # grouped backend is an optimization, not a prerequisite
        report["grouped_backend"] = {"status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}
    else:
        report["grouped_backend"] = {"status": "installed"}
    tau0 = float(candidate.tau0)
    h = float(candidate.h)
    nu = float(candidate.nu)
    points, weights, grid_meta = _load_grid()
    pressure_projection, pressure_meta = _load_pressure_projection(candidate)
    baseline, corrected, composition = _pressure_composition(candidate, pressure_projection, tau0)
    report["grid"] = grid_meta
    report["inputs"].update({"tau0": tau0, "similarity_exponent_h": h, "viscosity": nu})
    report["pressure_composition"] = composition
    report["status_detail"] = "fields_ready"
    report["pressure_projection"] = pressure_meta
    _save(output_path, report)
    print(json.dumps({"stage": "fields_ready", "point_count": int(len(points)), "chunk_size": int(chunk_size)}, sort_keys=True), flush=True)

    if case_mode in ("both", "baseline"):
        _run_case(baseline, "uncorrected_plus_pressure_projection", points, weights, tau0, h, nu, hspace, chunk_size, report, output_path)
    if case_mode in ("both", "corrected"):
        _run_case(corrected, "corrected_reference", points, weights, tau0, h, nu, hspace, chunk_size, report, output_path)
    if case_mode == 'baseline':
        report.update(status='completed',elapsed_seconds=time.perf_counter()-started)
        _save(output_path,report)
        return report
    base_metric = report["cases"]["uncorrected_plus_pressure_projection"]["metric"]
    corrected_metric = report["cases"]["corrected_reference"]["metric"]
    base_div = report["cases"]["uncorrected_plus_pressure_projection"]["divergence"]
    corrected_div = report["cases"]["corrected_reference"]["divergence"]
    base_l2 = float(base_metric["volume_L2"])
    base_max = float(base_metric["max_norm"])
    report["comparison"] = {
        "volume_L2_change": float(corrected_metric["volume_L2"] - base_l2),
        "volume_L2_relative_change": float((corrected_metric["volume_L2"] - base_l2) / max(base_l2, 1.0e-300)),
        "volume_L2_percent_change": float(100.0 * (corrected_metric["volume_L2"] - base_l2) / max(base_l2, 1.0e-300)),
        "max_norm_change": float(corrected_metric["max_norm"] - base_max),
        "max_norm_relative_change": float((corrected_metric["max_norm"] - base_max) / max(base_max, 1.0e-300)),
        "max_norm_percent_change": float(100.0 * (corrected_metric["max_norm"] - base_max) / max(base_max, 1.0e-300)),
        "divergence_max_abs_change": float(corrected_div["max_abs"] - base_div["max_abs"]),
        "divergence_weighted_RMS_change": float(corrected_div["weighted_RMS"] - base_div["weighted_RMS"]),
    }
    report["status"] = "completed"
    report["status_detail"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({"stage": "completed", "baseline_l2": base_l2, "corrected_l2": corrected_metric["volume_L2"], "baseline_max": base_max, "corrected_max": corrected_metric["max_norm"], "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    parser.add_argument("--hspace", type=float, default=HSPACE_DEFAULT)
    parser.add_argument("--step", type=Path, default=STEP_PATH)
    parser.add_argument("--case", choices=("both", "baseline", "corrected"), default="both")
    args = parser.parse_args()
    try:
        run(args.output, args.chunk_size, args.hspace, args.step, args.case)
    except Exception as exc:
        try:
            report = json.loads(args.output.read_text(encoding='utf-8'))
        except (OSError,ValueError):
            report = {}
        report.update(status='failed',accepted=False,pde_validated=False,
            scale_recursion_established=False,error=f'{type(exc).__name__}: {exc}')
        _save(args.output, report)
        raise


if __name__ == "__main__":
    main()
