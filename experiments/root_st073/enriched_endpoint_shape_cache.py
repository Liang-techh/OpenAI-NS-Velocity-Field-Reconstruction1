"""Degree-enriched finite-endpoint shape oracle.

This module reuses the previously frozen endpoint mean and mode-1 offset from
``wave_endpoint_shape_cache.npz``.  It only rebuilds the analytic tangent
columns for degree 3 mode 0, degree 2 mode 1, and degree 3 mode 2, giving the
264-column layout ``64 + 72 + 128``.  No mean endpoint replay is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from supported_fourier_analytic_jets import basis_jets
from wave_endpoint_shape_cache import (
    DEFAULT_CACHE as OLD_CACHE,
    FROZEN_GEOMETRY,
    ROOT,
    SHAPE_ROWS,
    _observables_and_jacobian,
    _sha256,
)


OLD_REPORT = ROOT / "wave_endpoint_shape_cache.json"
DEFAULT_CACHE = ROOT / "enriched_endpoint_shape_cache.npz"
DEFAULT_REPORT = ROOT / "enriched_endpoint_shape_cache.json"
DEFAULT_CANDIDATE = ROOT / "shape_direction_margin_tangent.json"
OLD_CANDIDATE = ROOT / "wave_moment_cone_tangent.json"
DEGREE_BY_MODE = {0: 3, 1: 2, 2: 3}
MODE_LIST = (0, 1, 2)
NU = 0.01


def _layout_for_mode(mode, degree):
    q = (degree + 1) ** 2
    velocity_q = 3 * q
    rows = []
    if mode == 0:
        for index in range(velocity_q):
            rows.append((mode, "velocity", "real", index))
        for index in range(q):
            rows.append((mode, "pressure", "real", index))
    else:
        for index in range(velocity_q):
            rows.append((mode, "velocity", "real", index))
            rows.append((mode, "velocity", "imag", index))
        for index in range(q):
            rows.append((mode, "pressure", "real", index))
            rows.append((mode, "pressure", "imag", index))
    return rows


def _all_layout():
    result = []
    for mode in MODE_LIST:
        result.extend(_layout_for_mode(mode, DEGREE_BY_MODE[mode]))
    return result


def _map_polynomial_index(index, old_degree, new_degree):
    old_q = (old_degree + 1) ** 2
    new_q = (new_degree + 1) ** 2
    component = int(index) // old_q
    polynomial = int(index) % old_q
    a, b = divmod(polynomial, old_degree + 1)
    if a > new_degree or b > new_degree:
        raise ValueError("Old polynomial index does not fit enriched degree")
    return component * new_q + a * (new_degree + 1) + b


def _old_to_new_mapping(old_layout, new_layout):
    new_index = {entry: index for index, entry in enumerate(new_layout)}
    mapping = []
    for mode, kind, component, index in old_layout:
        old_degree = 2
        new_degree = DEGREE_BY_MODE[int(mode)]
        if kind == "pressure":
            mapped_index = _map_polynomial_index(index, 0, 0)
            # Pressure indices are polynomial-only, unlike velocity block
            # indices. The helper above is used with one-component blocks.
            old_q = (old_degree + 1) ** 2
            new_q = (new_degree + 1) ** 2
            a, b = divmod(int(index), old_degree + 1)
            mapped_index = a * (new_degree + 1) + b
        else:
            mapped_index = _map_polynomial_index(index, old_degree, new_degree)
        key = (int(mode), kind, component, int(mapped_index))
        if key not in new_index:
            raise ValueError(f"Missing enriched layout entry for {key}")
        mapping.append(new_index[key])
    return np.asarray(mapping, dtype=int)


def _analytic_response(points, geometry, dtau):
    carriers = {
        0: np.zeros(2, dtype=float),
        1: np.asarray(geometry["carrier"], dtype=float),
        2: 2.0 * np.asarray(geometry["carrier"], dtype=float),
    }
    blocks_v = []
    blocks_j = []
    precision = []
    layout = []
    for mode in MODE_LIST:
        degree = DEGREE_BY_MODE[mode]
        q = (degree + 1) ** 2
        values, gradients, _, _, _ = basis_jets(
            points, geometry["center"], geometry["widths"], mode,
            degree, carriers[mode], NU
        )
        if mode == 0:
            for index in range(3 * q):
                blocks_v.append(dtau * values[:, :, index].real)
                blocks_j.append(dtau * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
        else:
            for index in range(3 * q):
                blocks_v.append(dtau * values[:, :, index].real)
                blocks_j.append(dtau * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets")
                blocks_v.append(-dtau * values[:, :, index].imag)
                blocks_j.append(-dtau * gradients[:, :, :, index].imag)
                layout.append((mode, "velocity", "imag", index))
                precision.append("analytic_basis_jets")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "imag", index))
                precision.append("zero_velocity_pressure")
    velocity = np.stack(blocks_v, axis=-1)
    gradient = np.stack(blocks_j, axis=-1)
    if velocity.shape != (len(points), 3, 264):
        raise ValueError(f"Unexpected enriched velocity response shape {velocity.shape}")
    if gradient.shape != (len(points), 3, 3, 264):
        raise ValueError(f"Unexpected enriched gradient response shape {gradient.shape}")
    return velocity, gradient, layout, precision


class EnrichedEndpointShape:
    """Evaluate endpoint shape observables for 264 enriched tangent controls."""

    def __init__(self, cache_path=DEFAULT_CACHE):
        cache_path = Path(cache_path)
        with np.load(cache_path, allow_pickle=False) as loaded:
            data = {key: loaded[key] for key in loaded.files}
        self.cache_path = cache_path
        self.points = np.asarray(data["points"], dtype=float)
        self.weights = np.asarray(data["weights"], dtype=float)
        self.velocity_offset = np.asarray(data["velocity_offset"], dtype=float)
        self.gradient_offset = np.asarray(data["gradient_offset"], dtype=float)
        self.velocity_control = np.asarray(data["velocity_control"], dtype=float)
        self.gradient_control = np.asarray(data["gradient_control"], dtype=float)
        self.reference = np.asarray(data["reference"], dtype=float)
        self.coefficients_original = np.asarray(data["coefficients_original"], dtype=float)
        self.k0 = float(np.asarray(data["k0"]).ravel()[0])
        self.k_endpoint = float(np.asarray(data["k_endpoint"]).ravel()[0])
        self.tau0 = float(np.asarray(data["tau0"]).ravel()[0])
        self.tau_endpoint = float(np.asarray(data["tau_endpoint"]).ravel()[0])
        self.dtau = float(np.asarray(data["dtau"]).ravel()[0])
        self.point_count = len(self.points)
        self.control_count = self.velocity_control.shape[-1]
        if self.control_count != 264:
            raise ValueError(f"Expected 264 controls, got {self.control_count}")

    def evaluate(self, control):
        control = np.asarray(control, dtype=float)
        if control.shape != (264,):
            raise ValueError(f"control must have shape (264,), got {control.shape}")
        velocity = self.velocity_offset + np.einsum("ncq,q->nc", self.velocity_control, control)
        gradient = self.gradient_offset + np.einsum("ncdq,q->ncd", self.gradient_control, control)
        values, jacobian, _ = _observables_and_jacobian(
            velocity, gradient, self.weights, self.points,
            self.velocity_control, self.gradient_control
        )
        return values, jacobian


def build_cache(cache_path=DEFAULT_CACHE, report_path=DEFAULT_REPORT,
                candidate_path=DEFAULT_CANDIDATE):
    started = time.perf_counter()
    old_report_raw = OLD_REPORT.read_bytes()
    old_report = json.loads(old_report_raw)
    if old_report.get("status") != "completed" or old_report.get("control_count") != 180:
        raise ValueError("Old endpoint cache report is not the expected completed 180-column cache")
    with np.load(OLD_CACHE, allow_pickle=False) as loaded:
        old = {key: loaded[key] for key in loaded.files}
    candidate_path = Path(candidate_path)
    candidate_raw = candidate_path.read_bytes()
    candidate = json.loads(candidate_raw)
    if candidate.get("status") != "completed":
        raise ValueError("Candidate must be completed")
    geometry_raw = FROZEN_GEOMETRY.read_bytes()
    geometry = json.loads(geometry_raw)["inputs"]["wave"]
    shape_rows_raw = SHAPE_ROWS.read_bytes()
    shape_rows = json.loads(shape_rows_raw)
    points = np.asarray(old["points"], dtype=float)
    weights = np.asarray(old["weights"], dtype=float)
    grid_hash = hashlib.sha256(points.tobytes() + weights.tobytes()).hexdigest()
    dtau = float(np.asarray(old["dtau"]).ravel()[0])
    velocity_control, gradient_control, new_layout, precision = _analytic_response(
        points, geometry, dtau
    )
    old_layout = [
        (int(entry["mode"]), entry["kind"], entry["component"], int(entry["index"]))
        for entry in old_report["control_layout"]
    ]
    mapping = _old_to_new_mapping(old_layout, new_layout)
    if len(mapping) != 180 or len(set(mapping.tolist())) != 180:
        raise ValueError("Old-to-new mapping is not one-to-one")

    cache_path = Path(cache_path)
    np.savez_compressed(
        cache_path,
        points=points,
        weights=weights,
        velocity_offset=np.asarray(old["velocity_offset"], dtype=float),
        gradient_offset=np.asarray(old["gradient_offset"], dtype=float),
        velocity_control=velocity_control,
        gradient_control=gradient_control,
        reference=np.asarray(old["reference"], dtype=float),
        coefficients_original=np.asarray(old["coefficients_original"], dtype=float),
        k0=np.asarray(old["k0"]), k_endpoint=np.asarray(old["k_endpoint"]),
        tau0=np.asarray(old["tau0"]), tau_endpoint=np.asarray(old["tau_endpoint"]),
        dtau=np.asarray(old["dtau"]),
        old_to_new_mapping=mapping,
        source_hash_old_cache=np.array([_sha256(OLD_CACHE)]),
        source_hash_old_report=np.array([_sha256(OLD_REPORT)]),
        source_hash_candidate=np.array([_sha256(candidate_path)]),
        source_hash_geometry=np.array([_sha256(FROZEN_GEOMETRY)]),
        source_hash_shape_rows=np.array([_sha256(SHAPE_ROWS)]),
        source_hash_grid=np.array([grid_hash]),
        source_hash_mean=old.get("source_hash_mean", np.array([old_report["mean_source_sha256"]]))
    )

    # Construct the class directly from the just-saved arrays for a single
    # immutable API/replay path.
    oracle = EnrichedEndpointShape(cache_path)
    old_oracle_module = __import__("wave_endpoint_shape_cache", fromlist=["EndpointShape"])
    old_oracle = old_oracle_module.EndpointShape(OLD_CACHE)
    candidate_checks = {}
    for label, path in (("original_wave_moment_cone", OLD_CANDIDATE),
                        ("shape_direction_margin", candidate_path)):
        raw = path.read_bytes()
        report = json.loads(raw)
        old_control = np.asarray(report["selected"]["tangent_coefficients"], dtype=float)
        enriched_control = np.zeros(264, dtype=float)
        enriched_control[mapping] = old_control
        old_values, old_jacobian = old_oracle.evaluate(old_control)
        new_values, new_jacobian = oracle.evaluate(enriched_control)
        candidate_checks[label] = dict(
            control_source=path.name,
            old_values=old_values.tolist(), new_values=new_values.tolist(),
            value_max_abs_error=float(np.max(np.abs(new_values - old_values))),
            mapped_jacobian_max_abs_error=float(
                np.max(np.abs(new_jacobian[:, mapping] - old_jacobian))
            ),
        )

    shape_control_old = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    shape_control = np.zeros(264, dtype=float)
    shape_control[mapping] = shape_control_old
    values, jacobian = oracle.evaluate(shape_control)
    scales = np.maximum(np.abs(shape_control), 1.0)
    direction = np.linspace(1.0, 2.0, 264) * scales
    direction /= np.linalg.norm(direction)
    fd_step = 1.0e3
    plus, _ = oracle.evaluate(shape_control + fd_step * direction)
    minus, _ = oracle.evaluate(shape_control - fd_step * direction)
    fd = (plus - minus) / (2.0 * fd_step)
    analytic = jacobian @ direction
    gradient_error = float(np.max(np.abs(fd - analytic)))
    gradient_relative_error = gradient_error / max(float(np.max(np.abs(fd))), 1.0e-30)

    report = dict(
        status="completed", accepted=False, pde_validated=False,
        scale_recursion_established=False, wave_integrated=False,
        source=candidate_path.name, source_sha256=_sha256(candidate_path),
        old_cache_source=OLD_CACHE.name, old_cache_sha256=_sha256(OLD_CACHE),
        old_report_source=OLD_REPORT.name, old_report_sha256=_sha256(OLD_REPORT),
        geometry_source=FROZEN_GEOMETRY.name, geometry_sha256=_sha256(FROZEN_GEOMETRY),
        shape_rows_source=SHAPE_ROWS.name, shape_rows_sha256=_sha256(SHAPE_ROWS),
        mean_source=old_report["mean_source"], mean_source_sha256=old_report["mean_source_sha256"],
        grid_hash=grid_hash, cache_path=cache_path.name,
        cache_sha256=_sha256(cache_path),
        reference_k=oracle.k0, endpoint_k=oracle.k_endpoint,
        tau_reference=oracle.tau0, tau_endpoint=oracle.tau_endpoint,
        dtau=oracle.dtau, point_count=int(oracle.point_count),
        control_count=264, mode_control_counts={"mode0_degree3": 64,
                                                  "mode1_degree2": 72,
                                                  "mode2_degree3": 128},
        pressure_column_count=64,
        precision_convention={
            "offset_reuse": "velocity_offset and gradient_offset copied bytewise from wave_endpoint_shape_cache.npz",
            "tangent_velocity_gradient": "double precision supported_fourier_analytic_jets multiplied by dtau",
            "pressure_columns": "exact zero velocity and gradient",
        },
        old_to_new_mapping=mapping.tolist(),
        old_control_layout=[dict(mode=a, kind=b, component=c, index=d) for a, b, c, d in old_layout],
        control_layout=[dict(mode=a, kind=b, component=c, index=d, precision=p)
                        for (a, b, c, d), p in zip(new_layout, precision)],
        reference_observables=np.asarray(old["reference"], dtype=float).tolist(),
        coefficients_original=np.asarray(old["coefficients_original"], dtype=float).tolist(),
        mapped_candidate_checks=candidate_checks,
        directional_gradient_check={
            "control_source": candidate_path.name, "step": fd_step,
            "analytic": analytic.tolist(), "finite_difference": fd.tolist(),
            "max_abs_error": gradient_error, "relative_error": gradient_relative_error,
        },
        scope=(
            "Degree-enriched fixed-cylinder endpoint oracle. The mean and initial "
            "wave offset are reused from the frozen 180-column cache; only analytic "
            "tangent columns are rebuilt for degree 3 mode 0, degree 2 mode 1, and "
            "degree 3 mode 2. No mean replay, PDE, trajectory, or scale-recursion claim."
        ),
        elapsed_seconds=float(time.perf_counter() - started),
    )
    report_path = Path(report_path)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "point_count": report["point_count"],
        "control_count": report["control_count"],
        "mapped_errors": {key: value["value_max_abs_error"]
                          for key, value in candidate_checks.items()},
        "gradient_relative_error": gradient_relative_error,
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    args = parser.parse_args()
    build_cache(args.cache, args.output, args.candidate)
