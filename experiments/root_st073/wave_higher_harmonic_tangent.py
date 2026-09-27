"""Bounded mode-2 polynomial enrichment for the locked wave tangent.

The cone-passing candidate in ``wave_moment_cone_tangent.json`` fixes the
instantaneous mode-1 wave and all mode-0/mode-1 tangent controls.  This
diagnostic replaces only the mode-2 derivative and pressure-gradient block
with degree 2, 3, and 4 compact Fourier bases.  The cached quadratic source
from ``wave_momentum_projection.npz`` is used on its training grid.  A dense
grid check uses the stored original-wave residual plus a directly recomputed
finite-difference delta for the locked initial wave and analytic tangent
columns for the fixed/new controls.

This is a local spatial tangent fit.  The added mode-2 derivative and
pressure terms vanish from the instantaneous velocity at the reference time;
no finite-time, cone, moment, PDE, or recursion acceptance is claimed.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field  # noqa: E402
from broad_wave_mean_fit import norms  # noqa: E402
from full_wave_tangent import (  # noqa: E402
    LocalPotentialField,
    _basis_columns,
)
from grouped_joined_field import install_in_field  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402
from wave_momentum_projection import residual_and_jacobian  # noqa: E402


CACHE_PATH = ROOT / "wave_momentum_projection.npz"
CACHE_REPORT_PATH = ROOT / "wave_momentum_projection.json"
CANDIDATE_PATH = ROOT / "wave_moment_cone_tangent.json"
WAVE_PATH = ROOT / "wave_dynamics_moment_codesign.json"
DENSE_PATH = ROOT / "full_wave_dense_tangent.json"
OUTPUT_PATH = ROOT / "wave_higher_harmonic_tangent.json"

MODE2_START = 108  # mode 0 has 36 real columns; mode 1 has 72
RCOND = 1.0e-10
DEGREES = (2, 3, 4)


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _packed_complex(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _decode_packed(value):
    value = np.asarray(value, dtype=float)
    return value[..., 0] + 1j * value[..., 1]


def _metric(residual, weights):
    result = norms(np.asarray(residual), np.asarray(weights))
    result["point_count"] = int(len(residual))
    result["physical_volume"] = float(np.sum(weights))
    return result


def _mode_block_columns(points, center, widths, mode, degree, carrier):
    """Return [velocity derivative, pressure gradient] real columns.

    Nonzero Fourier modes use the same interleaved ``[Re, -Im]`` convention as
    ``full_wave_tangent._basis_columns``.  The returned matrix has
    ``2 * (3*q + q)`` columns for a degree with ``q`` scalar polynomials.
    """

    velocity, _, pressure_gradient = basis_data(
        points, center, widths, mode, degree, carrier
    )
    blocks = []
    qv = velocity.shape[-1]
    qp = pressure_gradient.shape[-1]
    for tensor, count in ((velocity, qv), (pressure_gradient, qp)):
        for j in range(count):
            blocks.append(tensor[:, :, j].real.reshape(-1))
            blocks.append(-tensor[:, :, j].imag.reshape(-1))
    return np.stack(blocks, axis=1)


def _decode_mode_block(block, degree):
    """Decode an interleaved nonzero-mode tangent block."""

    q = (degree + 1) ** 2
    velocity_count = 3 * q
    block = np.asarray(block, dtype=float)
    expected = 2 * (velocity_count + q)
    if block.shape != (expected,):
        raise ValueError(f"expected mode block shape {(expected,)}, got {block.shape}")
    velocity = block[: 2 * velocity_count: 2] + 1j * block[1: 2 * velocity_count: 2]
    pressure_start = 2 * velocity_count
    pressure = (
        block[pressure_start:pressure_start + 2 * q:2]
        + 1j * block[pressure_start + 1:pressure_start + 2 * q:2]
    )
    return velocity, pressure


class Mode2TangentCorrection:
    """Add only a fitted mode-2 tangent correction to an existing field.

    ``base_field`` should already contain the locked wave and fixed mode-0 /
    mode-1 controls.  The velocity derivative is multiplied by
    ``tau0 - tau`` because physical time is ``t = -tau``.  Pressure is a
    constant local pressure correction, so its spatial gradient enters the
    momentum residual without that time factor.
    """

    def __init__(self, base_field, center, widths, carrier_mode1, degree,
                 derivative_coefficients, pressure_coefficients, tau0):
        self.base_field = base_field
        self.center = tuple(float(v) for v in center)
        self.widths = tuple(float(v) for v in widths)
        self.carrier_mode2 = 2.0 * np.asarray(carrier_mode1, dtype=float)
        self.degree = int(degree)
        self.derivative_coefficients = np.asarray(derivative_coefficients, dtype=complex)
        self.pressure_coefficients = np.asarray(pressure_coefficients, dtype=complex)
        self.tau0 = float(tau0)

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        tau = float(np.asarray(tau).ravel()[0])
        velocity, pressure = self.base_field.fields(points, tau)
        basis_velocity, basis_pressure, _ = basis_data(
            points, self.center, self.widths, 2, self.degree, self.carrier_mode2
        )
        physical_delta = self.tau0 - tau
        velocity = velocity + physical_delta * np.einsum(
            "niq,q->ni", basis_velocity, self.derivative_coefficients
        ).real
        pressure = pressure + np.einsum(
            "nq,q->n", basis_pressure, self.pressure_coefficients
        ).real
        return velocity, pressure


def _weighted_fit(design, base_residual, weights, rcond=RCOND):
    """Fit design * coefficient to ``-base_residual`` after volume weighting."""

    weights = np.asarray(weights, dtype=float)
    row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_design = np.asarray(design, dtype=float) * row_weight[:, None]
    weighted_rhs = -np.asarray(base_residual, dtype=float).reshape(-1) * row_weight
    scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized = weighted_design / scales[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    cutoff = float(rcond) * max(float(singular[0]), 1.0e-300)
    rank = int(np.sum(singular > cutoff))
    scaled, _, lstsq_rank, _ = np.linalg.lstsq(
        normalized, weighted_rhs, rcond=float(rcond)
    )
    coefficients = scaled / scales
    return coefficients, {
        "rcond": float(rcond),
        "column_count": int(design.shape[1]),
        "rank_from_singular_values": rank,
        "rank_from_lstsq": int(lstsq_rank),
        "singular_values": singular.tolist(),
        "cutoff": cutoff,
        "column_scales": scales.tolist(),
    }


def _candidate_mode_snapshot(candidate, wave_report):
    selected = candidate["selected"]
    tangent = np.asarray(selected["tangent_coefficients"], dtype=float)
    if tangent.shape != (180,):
        raise ValueError(f"expected 180 candidate tangent coefficients, got {tangent.shape}")
    wave = _decode_packed(wave_report["selected"]["coefficients_original"])
    if wave.shape != (27,):
        raise ValueError(f"expected 27 complex initial-wave coefficients, got {wave.shape}")
    return tangent, wave


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _build_dense_base(mean, points, center, widths, carrier, wave, tau, hspace, htime):
    """Direct five-point FD residual for the locked initial wave on new grid."""

    q = 9
    # The LocalPotentialField evaluates only degree-2 initial modes here.  All
    # tangent derivatives and pressure coefficients are zero, so this is the
    # instantaneous mean plus the locked wave state.
    initial = (np.zeros(3 * q, complex), wave, np.zeros(3 * q, complex))
    zeros_d = tuple(np.zeros(3 * q, complex) for _ in (0, 1, 2))
    zeros_p = tuple(np.zeros(q, complex) for _ in (0, 1, 2))
    field = LocalPotentialField(
        mean, center, widths, 2,
        {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier},
        initial, zeros_d, zeros_p, tau,
    )
    return momentum(jets(field, points, tau, hspace, htime))


def _independent_design(points, center, widths, carrier, degree, x_fixed, x_mode2):
    """Compose fixed mode 0/1 and a selected mode-2 correction on new grid."""

    carriers = {0: np.zeros(2), 1: np.asarray(carrier), 2: 2.0 * np.asarray(carrier)}
    fixed_all, _ = _basis_columns(points, center, widths, carriers, 2)
    fixed = fixed_all[:, :MODE2_START]
    mode2 = _mode_block_columns(points, center, widths, 2, degree, 2.0 * np.asarray(carrier))
    return fixed @ x_fixed + mode2 @ x_mode2, fixed, mode2


def run():
    started = time.perf_counter()
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    wave_report = json.loads(WAVE_PATH.read_text(encoding="utf-8"))
    dense = json.loads(DENSE_PATH.read_text(encoding="utf-8"))
    geometry_name = str(wave_report.get("source", "wave_stress_growth_codesign.json"))
    geometry_report = json.loads((ROOT / geometry_name).read_text(encoding="utf-8"))
    candidate_x, wave = _candidate_mode_snapshot(candidate, wave_report)
    dense_cache = dense["new_frozen_cache"]
    dense_points = np.asarray(dense_cache["points"], dtype=float)
    dense_weights = np.asarray(dense_cache["weights"], dtype=float)
    dense_original_residual = np.asarray(dense_cache["residual"], dtype=float)
    center = tuple(float(v) for v in geometry_report["center"])
    widths = tuple(float(v) for v in geometry_report["widths"])
    carrier = np.asarray(geometry_report["carrier"], dtype=float)
    tau = float(dense["initial_tau"])
    hspace = float(dense["timesteps"]["hspace"])
    htime = float(dense["timesteps"]["htime"])
    viscosity = float(dense["viscosity"])
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "wave_integrated": False,
        "scope": (
            "Spatial mode-2 degree enrichment only. The locked instantaneous "
            "wave and mode-0/mode-1 tangent controls are fixed; only the mode-2 "
            "derivative and pressure-gradient block is replaced. Added controls "
            "vanish from instantaneous velocity at the reference time. No finite "
            "time, cone, moment, PDE, or recursion acceptance is claimed."
        ),
        "sources": {
            "candidate": {"path": CANDIDATE_PATH.name, "sha256": _sha256(CANDIDATE_PATH)},
            "wave_initial": {"path": WAVE_PATH.name, "sha256": _sha256(WAVE_PATH)},
            "wave_geometry": {"path": geometry_name, "sha256": _sha256(ROOT / geometry_name)},
            "projection_cache": {"path": CACHE_PATH.name, "sha256": _sha256(CACHE_PATH)},
            "projection_report": {"path": CACHE_REPORT_PATH.name, "sha256": _sha256(CACHE_REPORT_PATH)},
            "dense_report": {"path": DENSE_PATH.name, "sha256": _sha256(DENSE_PATH)},
        },
        "inputs": {
            "center": list(center),
            "widths": list(widths),
            "carrier_mode1": carrier.tolist(),
            "carrier_mode2": (2.0 * carrier).tolist(),
            "tau": tau,
            "physical_time": "t=-tau",
            "viscosity": viscosity,
            "hspace": hspace,
            "htime": htime,
            "degrees": list(DEGREES),
            "rcond": RCOND,
            "candidate_tangent_layout": "mode0 [36 real], then mode1 [velocity/pressure interleaved Re,-Im], then mode2 [velocity/pressure interleaved Re,-Im]",
            "mode2_start": MODE2_START,
            "candidate_tangent_coefficients": candidate_x.tolist(),
            "locked_wave_coefficients": _packed_complex(wave),
            "locked_wave_source": WAVE_PATH.name,
            "geometry_source": geometry_name,
        },
        "dense_grid": {
            "point_count": int(len(dense_points)),
            "weight_sum": float(np.sum(dense_weights)),
            "geometry": dense.get("new_grid", {}),
            "original_frozen_residual_definition": "full_wave_dense_tangent.new_frozen_cache residual for its original wave, before tangent corrections",
            "original_frozen_residual_metric": _metric(dense_original_residual, dense_weights),
        },
        "provenance_note": (
            "The projection cache was built from wave_momentum_projection's earlier "
            "wave source, while this locked candidate uses wave_dynamics_moment_codesign. "
            "R0/A/B/L and tangent_design are geometry/mean/basis operators; evaluating "
            "them at the locked wave reproduces the candidate training residual exactly."
        ),
        "reusable_wrapper": {
            "class": "Mode2TangentCorrection",
            "velocity_rule": "(tau0 - tau) * Re(V_mode2_degree_d @ derivative_coefficients)",
            "pressure_rule": "Re(P_mode2_degree_d @ pressure_coefficients), constant in physical time",
            "instantaneous_wave_unchanged": True,
        },
    }
    _save(report)
    print(json.dumps({"stage": "initialized", "training_points": int(len(np.load(CACHE_PATH, allow_pickle=False)["points"])), "dense_points": int(len(dense_points))}), flush=True)

    with np.load(CACHE_PATH, allow_pickle=False) as cached:
        cache = {key: cached[key] for key in cached.files}
    points = np.asarray(cache["points"], dtype=float)
    weights = np.asarray(cache["weights"], dtype=float)
    tangent_design = np.asarray(cache["tangent_design"], dtype=float)
    wave_dynamic = wave
    wave_x = np.concatenate((wave_dynamic.real, wave_dynamic.imag))
    dynamic_base, _ = residual_and_jacobian(cache, wave_x)
    fixed_training = tangent_design[:, :MODE2_START]
    mode2_training_seed_design = tangent_design[:, MODE2_START:]
    fixed_residual = dynamic_base + (fixed_training @ candidate_x[:MODE2_START]).reshape(-1, 3)
    seed_mode2 = candidate_x[MODE2_START:]
    seed_residual = fixed_residual + (mode2_training_seed_design @ seed_mode2).reshape(-1, 3)
    candidate_metric = _metric(seed_residual, weights)
    report["training_source"] = {
        "cache_arrays": {key: list(value.shape) for key, value in cache.items()},
        "point_count": int(len(points)),
        "base_locked_wave_metric": _metric(dynamic_base, weights),
        "fixed_mode0_mode1_metric": _metric(fixed_residual, weights),
        "locked_candidate_seed_metric": candidate_metric,
        "candidate_report_metric": candidate["selected"]["training_momentum"],
        "candidate_metric_l2_abs_difference": float(candidate_metric["volume_L2"] - candidate["selected"]["training_momentum"]["momentum_volume_L2"]),
        "candidate_reconstruction": "R0 + L(c_locked) + (B c_locked)(A c_locked) + tangent_design @ candidate_tangent_coefficients",
    }
    if abs(report["training_source"]["candidate_metric_l2_abs_difference"]) > 1.0e-5:
        raise ValueError("cached reconstruction does not reproduce locked candidate training metric")
    report["status"] = "training_source_reconstructed"
    _save(report)
    print(json.dumps({"stage": "training_source_reconstructed", "seed_l2": candidate_metric["volume_L2"]}), flush=True)

    fit_results = {}
    for degree in DEGREES:
        design = _mode_block_columns(points, center, widths, 2, degree, 2.0 * carrier)
        if degree == 2:
            if design.shape != mode2_training_seed_design.shape:
                raise ValueError("degree-2 replacement design shape differs from cached mode-2 block")
            relative = np.linalg.norm(design - mode2_training_seed_design) / max(np.linalg.norm(mode2_training_seed_design), 1.0e-300)
        else:
            relative = None
        fitted, svd = _weighted_fit(design, fixed_residual, weights, RCOND)
        fitted_residual = fixed_residual + (design @ fitted).reshape(-1, 3)
        fit_results[str(degree)] = {
            "degree": degree,
            "q": (degree + 1) ** 2,
            "design_shape": list(design.shape),
            "degree2_design_relative_error": relative,
            "seed_coefficients": seed_mode2.tolist() if degree == 2 else None,
            "fitted_coefficients": fitted.tolist(),
            "fitted_coefficients_complex": {
                "velocity": _packed_complex(_decode_mode_block(fitted, degree)[0]),
                "pressure": _packed_complex(_decode_mode_block(fitted, degree)[1]),
            },
            "svd": svd,
            "seed_metric": candidate_metric if degree == 2 else None,
            "fitted_metric": _metric(fitted_residual, weights),
            "fitted_residual_definition": "fixed locked-wave/mode0/mode1 source plus analytic mode-2 degree-d tangent columns",
        }
        report["fits"] = fit_results
        report["status"] = "training_fits_in_progress"
        _save(report)
        print(json.dumps({"stage": "training_fit_done", "degree": degree,
                          "rank": svd["rank_from_singular_values"],
                          "seed_l2": candidate_metric["volume_L2"] if degree == 2 else None,
                          "fit_l2": fit_results[str(degree)]["fitted_metric"]["volume_L2"]}), flush=True)

    report["status"] = "training_fits_complete"
    _save(report)
    print(json.dumps({"stage": "training_fits_complete"}), flush=True)

    # The dense cache residual is for the earlier wave source.  Recompute only
    # the locked-wave base residual by direct 5-point FD jets, then add the
    # analytic tangent deltas.  This avoids a second full field replay for each
    # polynomial degree while retaining the independent spatial grid.
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    print(json.dumps({"stage": "dense_base_fd_started", "point_count": int(len(dense_points))}), flush=True)
    dense_locked_base = _build_dense_base(
        mean, dense_points, center, widths, carrier, wave_dynamic, tau, hspace, htime,
    )
    dense_base_delta = dense_locked_base - dense_original_residual
    report["dense_grid"]["locked_wave_base_fd_metric"] = _metric(dense_locked_base, dense_weights)
    report["dense_grid"]["stored_original_to_locked_base_delta_metric"] = _metric(dense_base_delta, dense_weights)
    report["dense_grid"]["independent_method"] = (
        "stored full_wave_dense_tangent.new_frozen_cache original-wave residual + "
        "direct affine_momentum.jets five-point FD delta to locked wave + analytic "
        "basis_data tangent columns"
    )
    report["status"] = "dense_base_fd_complete"
    _save(report)
    print(json.dumps({"stage": "dense_base_fd_complete", "locked_base_l2": report["dense_grid"]["locked_wave_base_fd_metric"]["volume_L2"]}), flush=True)

    dense_carriers = {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier}
    dense_fixed_all, dense_layout = _basis_columns(
        dense_points, center, widths, dense_carriers, 2
    )
    dense_fixed = dense_fixed_all[:, :MODE2_START]
    dense_fixed_delta = (dense_fixed @ candidate_x[:MODE2_START]).reshape(-1, 3)
    dense_results = {}
    for degree in DEGREES:
        fitted = np.asarray(fit_results[str(degree)]["fitted_coefficients"], dtype=float)
        mode2_design = _mode_block_columns(
            dense_points, center, widths, 2, degree, 2.0 * carrier
        )
        dense_predicted = (
            dense_locked_base
            + dense_fixed_delta
            + (mode2_design @ fitted).reshape(-1, 3)
        )
        dense_results[str(degree)] = {
            "degree": degree,
            "mode2_design_shape": list(mode2_design.shape),
            "metric": _metric(dense_predicted, dense_weights),
            "residual_definition": "direct FD locked-wave base plus fixed mode0/mode1 analytic response plus fitted mode2 analytic response",
        }
    report["dense_grid"]["fixed_tangent_layout"] = dense_layout[:MODE2_START]
    report["dense_grid"]["fits"] = dense_results
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "training": {degree: fit_results[str(degree)]["fitted_metric"]["volume_L2"] for degree in DEGREES},
        "dense": {degree: dense_results[str(degree)]["metric"]["volume_L2"] for degree in DEGREES},
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
