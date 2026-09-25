"""Sequential scale-transfer correction for the width-six ST073 modes.

This is a bounded numerical screen.  It starts from the published fixed-scale
``WideJointModes`` amplitudes, fits the right knot of each adjacent k band in
sequence, and evaluates the resulting C2 quintic scale interpolation on
off-grid holdouts.  Every residual uses the independent Cartesian finite
difference jet, including the time derivative of the interpolated amplitudes.
The mode constructors retain the existing divergence-free and interface
vanishing structure.  No PDE or continuum acceptance claim is made.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from affine_momentum import jets, momentum
from wide_modes import CachedWidth, WideJointModes


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "scale_transfer_recurrence.json"
BASE_REPORT = ROOT / "wide_collocation" / "report.json"
NU = 0.01
RATIO = 6.0
KNOTS = np.array([1.0, 3.0, 5.0], dtype=float)
NMODES = 24


def quintic_interpolate(k: float, values: np.ndarray, knots: np.ndarray = KNOTS) -> np.ndarray:
    """C2 smoothstep using only the two endpoint states of each band.

    The zero first and second derivatives at knots keep a completed earlier
    band unchanged when the next right knot is fitted.
    """
    values = np.asarray(values, dtype=float)
    if values.shape != (len(knots), NMODES):
        raise ValueError("expected one length-24 vector per knot")
    if k < knots[0] - 1e-12 or k > knots[-1] + 1e-12:
        raise ValueError(f"k={k} outside interpolation interval")
    interval = int(np.clip(np.searchsorted(knots, k, side="right") - 1, 0, len(knots) - 2))
    h = float(knots[interval + 1] - knots[interval])
    s = float(np.clip((k - knots[interval]) / h, 0.0, 1.0))
    s3, s4, s5 = s**3, s**4, s**5
    blend = 10.0 * s3 - 15.0 * s4 + 6.0 * s5
    i, j = interval, interval + 1
    return (1.0 - blend) * values[i] + blend * values[j]


class ScaleField:
    """Existing width-six modes with C2 amplitudes interpolated in k."""

    def __init__(self, base: CachedWidth, knot_amplitudes: np.ndarray):
        self.base = base
        self.knots = np.asarray(knot_amplitudes, dtype=float).copy()
        if self.knots.shape != (3, NMODES):
            raise ValueError("expected three by 24 knot amplitudes")
        self.nu = base.nu

    @staticmethod
    def k_from_tau(tau: float) -> float:
        tau = float(tau)
        if not np.isfinite(tau) or tau <= 0:
            raise ValueError("positive finite tau required")
        return float(-np.log2(2.0 * tau))

    def amplitudes(self, tau: float) -> np.ndarray:
        return quintic_interpolate(self.k_from_tau(tau), self.knots)

    def fields(self, points, tau):
        # The independent jet evaluates one physical time band at a time.
        ts = np.asarray(tau, dtype=float)
        if ts.ndim != 0 and not np.all(ts == ts.flat[0]):
            raise ValueError("ScaleField.fields expects one tau per jet call")
        scalar_tau = float(ts.flat[0])
        amplitudes = self.amplitudes(scalar_tau)
        return WideJointModes(self.base, amplitudes).fields(points, scalar_tau)


def similarity_points(base: CachedWidth, k: float, etas, ys, angles) -> np.ndarray:
    """Cartesian points in the existing inner-to-outer width-six annulus."""
    tau = 0.5 * 2.0 ** (-float(k))
    eta, y, angle = np.broadcast_arrays(np.asarray(etas), np.asarray(ys), np.asarray(angles))
    q = tau / (1.0 - eta * eta)
    ri = np.sqrt(2.0 * base.nu * q * (3.0 / 64.0))
    r = ri * (1.0 + (RATIO - 1.0) * y)
    z = np.sqrt(base.nu) * q ** base.inner.D * eta
    return np.column_stack((r * np.cos(angle), r * np.sin(angle), z))


def tau_for_k(k: float) -> float:
    return float(0.5 * 2.0 ** (-float(k)))


def training_samples(base: CachedWidth, k_values: tuple[float, ...]) -> list[dict]:
    samples = []
    etas = np.array([-0.25, 0.25])
    ys = np.array([0.25, 0.50, 0.75])
    angles = np.array([0.0, 0.5 * np.pi])
    for k in k_values:
        eta, y, angle = np.meshgrid(etas, ys, angles, indexing="ij")
        samples.append({"k": float(k), "tau": tau_for_k(k), "points": similarity_points(base, k, eta.ravel(), y.ravel(), angle.ravel())})
    return samples


def concatenate_jets(field, samples: list[dict], htime_factor: float = 0.00025):
    rows = []
    for sample in samples:
        tau = sample["tau"]
        hspace = 0.0005 * np.sqrt(field.nu * tau)
        htime = htime_factor * tau
        rows.append(jets(field, sample["points"], tau, hspace, htime))
    return tuple(np.concatenate([row[i] for row in rows], axis=0) for i in range(3))


def mode_difference_jets(base: CachedWidth, samples: list[dict], knot_index: int):
    """Jets for unit correction modes at one knot, minus the zero field."""
    zero = ScaleField(base, np.zeros((3, NMODES)))
    zero_jets = concatenate_jets(zero, samples)
    mode_rows = []
    for mode_index in range(NMODES):
        knots = np.zeros((3, NMODES))
        knots[knot_index, mode_index] = 1.0
        trial = concatenate_jets(ScaleField(base, knots), samples)
        mode_rows.append(tuple(trial[i] - zero_jets[i] for i in range(3)))
    mode_u = np.stack([row[0] for row in mode_rows], axis=0)
    mode_g = np.stack([row[1] for row in mode_rows], axis=0)
    mode_l = np.stack([row[2] for row in mode_rows], axis=0)
    return zero_jets, mode_u, mode_g, mode_l


def residual_and_jacobian(increment, fixed_jets, mode_u, mode_g, mode_l, start, scales, regularization=1e-3):
    """Full Cartesian nonlinear momentum and its exact algebraic Jacobian."""
    amplitudes = start + np.asarray(increment, dtype=float)
    fixed_u, fixed_g, fixed_l = fixed_jets
    u = fixed_u + np.einsum("kni,k->ni", mode_u, amplitudes)
    g = fixed_g + np.einsum("knij,k->nij", mode_g, amplitudes)
    linear = fixed_l + np.einsum("kni,k->ni", mode_l, amplitudes)
    residual = momentum((u, g, linear))
    jac_res = mode_l + np.einsum("knij,nj->kni", mode_g, u) + np.einsum("nij,knj->kni", g, mode_u)
    # The caller passes one scale per point, so the block split is represented
    # by a flat vector and the same normalization is applied pointwise.
    residual_flat = (residual / scales[:, None]).ravel()
    jac_flat = (jac_res / scales[None, :, None]).transpose(1, 2, 0).reshape(len(residual) * 3, NMODES)
    return np.r_[residual_flat, regularization * np.asarray(increment)], np.vstack((jac_flat, regularization * np.eye(NMODES)))


def fit_band(base: CachedWidth, knot_amplitudes: np.ndarray, right_index: int, samples: list[dict]):
    """Fit one right knot from the preceding knot state."""
    fixed_knots = np.asarray(knot_amplitudes, dtype=float).copy()
    start = fixed_knots[right_index].copy()
    fixed_knots[right_index] = 0.0
    fixed_jets = concatenate_jets(ScaleField(base, fixed_knots), samples)
    _, mode_u, mode_g, mode_l = mode_difference_jets(base, samples, right_index)

    # Pointwise scales keep the high-k band from dominating only because its
    # residual has a larger physical magnitude.  They are computed from the
    # uncorrected right-knot state, before fitting this band.
    reference_u = fixed_jets[0] + np.einsum("kni,k->ni", mode_u, start)
    reference_g = fixed_jets[1] + np.einsum("knij,k->nij", mode_g, start)
    reference_l = fixed_jets[2] + np.einsum("kni,k->ni", mode_l, start)
    reference_residual = momentum((reference_u, reference_g, reference_l))
    point_scales = np.maximum(np.linalg.norm(reference_residual, axis=1), 1e-6)

    def objective(delta):
        return residual_and_jacobian(delta, fixed_jets, mode_u, mode_g, mode_l, start, point_scales)[0]

    def jacobian(delta):
        return residual_and_jacobian(delta, fixed_jets, mode_u, mode_g, mode_l, start, point_scales)[1]

    # Keep each sequential correction bounded and preserve the existing mode
    # convention's [-4,4] amplitude envelope wherever possible.
    increment_cap = 0.35
    lower = np.maximum(-increment_cap, -4.0 - start)
    upper = np.minimum(increment_cap, 4.0 - start)
    result = least_squares(
        objective,
        np.zeros(NMODES),
        jac=jacobian,
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
        max_nfev=8,
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-8,
        verbose=0,
    )
    update = np.asarray(result.x, dtype=float)
    knot_amplitudes[right_index] = start + update
    return result, update, point_scales, objective(np.zeros(NMODES)), objective(update)


def holdout_points(base: CachedWidth, k: float):
    eta_values = np.array([-0.23, 0.17])
    g, _ = leggauss(10)
    ys = (g + 1.0) / 2.0
    angles = np.array([0.0, 0.25 * np.pi, 0.5 * np.pi, 0.75 * np.pi])
    eta, y, angle = np.meshgrid(eta_values, ys, angles, indexing="ij")
    return similarity_points(base, k, eta.ravel(), y.ravel(), angle.ravel())


def evaluate_peak(base: CachedWidth, field, k: float):
    points = holdout_points(base, k)
    tau = tau_for_k(k)
    jet = jets(field, points, tau, 0.0005 * np.sqrt(field.nu * tau), 0.00025 * tau)
    residual = momentum(jet)
    divergence = np.trace(jet[1], axis1=1, axis2=2)
    return {
        "k": float(k),
        "eta_holdout": [-0.23, 0.17],
        "points": int(len(points)),
        "peak": float(np.max(np.linalg.norm(residual, axis=1))),
        "volume_l2_sample_rms": float(np.sqrt(np.mean(np.sum(residual * residual, axis=1)))),
        "divergence_max": float(np.max(np.abs(divergence))),
    }


def evaluate_volume(base: CachedWidth, field, k: float, ny: int = 8, ne: int = 6):
    gy, wy = leggauss(ny)
    ge, we = leggauss(ne)
    eta = np.repeat(0.4 * ge, ny)
    y = np.tile((gy + 1.0) / 2.0, ne)
    tau = tau_for_k(k)
    q = tau / (1.0 - eta * eta)
    ri = np.sqrt(2.0 * base.nu * q * (3.0 / 64.0))
    r = ri * (1.0 + (RATIO - 1.0) * y)
    z = np.sqrt(base.nu) * q ** base.inner.D * eta
    points = np.column_stack((r, np.zeros_like(r), z))
    jet = jets(field, points, tau, 0.0005 * np.sqrt(field.nu * tau), 0.00025 * tau)
    residual = momentum(jet)
    divergence = np.trace(jet[1], axis1=1, axis2=2)
    ze = np.sqrt(base.nu) * q ** base.inner.D * (1.0 + 2.0 * base.inner.D * eta * eta / (1.0 - eta * eta))
    weights = np.repeat(0.4 * we, ny) * np.tile(wy / 2.0, ne) * 2.0 * np.pi * r * (RATIO - 1.0) * ri * ze
    return {
        "k": float(k),
        "volume": float(np.sum(weights)),
        "momentum_peak": float(np.max(np.linalg.norm(residual, axis=1))),
        "momentum_volume_L2": float(np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)))),
        "divergence_max": float(np.max(np.abs(divergence))),
        "quadrature": {"radial_order": ny, "eta_order": ne, "eta_half_width": 0.4},
    }


def product_rule_delta(base: CachedWidth, field: ScaleField, k: float):
    """Compare varying-scale FD time jets with the frozen-amplitude jet."""
    tau = tau_for_k(k)
    points = holdout_points(base, k)[:16]
    hspace = 0.0005 * np.sqrt(field.nu * tau)
    htime = 0.00025 * tau
    varying = jets(field, points, tau, hspace, htime)
    frozen = WideJointModes(base, field.amplitudes(tau))
    fixed = jets(frozen, points, tau, hspace, htime)
    delta = varying[2] - fixed[2]
    return {
        "k": float(k),
        "points": int(len(points)),
        "max_linear_jet_difference": float(np.max(np.linalg.norm(delta, axis=1))),
        "rms_linear_jet_difference": float(np.sqrt(np.mean(np.sum(delta * delta, axis=1)))),
        "interpretation": "Independent Cartesian FD time jets include the product-rule contribution from d a(k(tau))/dt; frozen comparison holds amplitudes at the central k.",
    }


def interface_check(base: CachedWidth, field: ScaleField, k: float):
    tau = tau_for_k(k)
    eta = np.array([-0.23, 0.17])
    q = tau / (1.0 - eta * eta)
    ri = np.sqrt(2.0 * base.nu * q * (3.0 / 64.0))
    z = np.sqrt(base.nu) * q ** base.inner.D * eta
    rows = []
    for side, factor in (("inner", 1.0), ("outer", RATIO)):
        points = np.column_stack((factor * ri, np.zeros_like(ri), z))
        u0, p0 = base.fields(points, tau)
        u1, p1 = field.fields(points, tau)
        rows.append({"side": side, "velocity_difference_max": float(np.max(np.abs(u1 - u0))), "pressure_difference_max": float(np.max(np.abs(p1 - p0)))})
    return rows


def run(output: Path = OUT) -> dict:
    started = time.time()
    if not BASE_REPORT.is_file():
        raise FileNotFoundError(f"required fixed-scale source missing: {BASE_REPORT}")
    source = json.loads(BASE_REPORT.read_text())
    initial_amplitudes = np.asarray(source["amplitudes"], dtype=float)
    if initial_amplitudes.shape != (NMODES,):
        raise ValueError("wide_collocation report does not contain 24 amplitudes")
    base = CachedWidth(6.0)
    initial_knots = np.repeat(initial_amplitudes[None, :], len(KNOTS), axis=0)
    knots = initial_knots.copy()
    band_specs = ((0, 1.5, 2.5), (1, 3.5, 4.5))
    stages = []
    for stage, (right_index, *k_values) in enumerate(band_specs, start=1):
        samples = training_samples(base, tuple(k_values))
        result, update, scales, before, after = fit_band(base, knots, right_index, samples)
        stage_row = {
            "stage": stage,
            "band": [float(KNOTS[right_index - 1]), float(KNOTS[right_index])],
            "training_k": [float(k) for k in k_values],
            "training_points": int(sum(len(sample["points"]) for sample in samples)),
            "right_knot": float(KNOTS[right_index]),
            "success": bool(result.success),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "initial_normalized_loss": float(before @ before),
            "final_normalized_loss": float(after @ after),
            "increment_max_abs": float(np.max(np.abs(update))),
            "increment_l2": float(np.linalg.norm(update)),
            "increment": update.tolist(),
            "coefficient_norm_before": float(np.linalg.norm(knots[right_index] - update)),
            "coefficient_norm_after": float(np.linalg.norm(knots[right_index])),
            "full_cartesian_components": 3,
            "fd_jet": "affine_momentum.jets fourth-order Cartesian space and centered tau finite difference; interpolation varies inside every stencil",
        }
        stages.append(stage_row)

    baseline_field = ScaleField(base, initial_knots)
    corrected_field = ScaleField(base, knots)
    holdout_ks = (2.25, 4.25)
    holdouts = []
    volumes = []
    product_rule = []
    interfaces = []
    for k in holdout_ks:
        base_peak = evaluate_peak(base, baseline_field, k)
        corrected_peak = evaluate_peak(base, corrected_field, k)
        holdouts.append({"k": float(k), "baseline": base_peak, "sequential_scale": corrected_peak, "peak_improved": corrected_peak["peak"] < base_peak["peak"], "sample_rms_improved": corrected_peak["volume_l2_sample_rms"] < base_peak["volume_l2_sample_rms"]})
        base_volume = evaluate_volume(base, baseline_field, k)
        corrected_volume = evaluate_volume(base, corrected_field, k)
        volumes.append({"k": float(k), "baseline": base_volume, "sequential_scale": corrected_volume, "volume_l2_improved": corrected_volume["momentum_volume_L2"] < base_volume["momentum_volume_L2"], "volume_peak_improved": corrected_volume["momentum_peak"] < base_volume["momentum_peak"]})
        product_rule.append(product_rule_delta(base, corrected_field, k))
        interfaces.extend(interface_check(base, corrected_field, k))

    base_peaks = [row["baseline"]["peak"] for row in holdouts]
    corrected_peaks = [row["sequential_scale"]["peak"] for row in holdouts]
    base_l2 = [row["baseline"]["momentum_volume_L2"] for row in volumes]
    corrected_l2 = [row["sequential_scale"]["momentum_volume_L2"] for row in volumes]
    report = {
        "schema": "st073_scale_transfer_recurrence_v1",
        "experiment": "sequential k-scale correction with C2 quintic amplitude interpolation",
        "status": "bounded research screen; no PDE acceptance claim",
        "source": {"fixed_scale_report": str(BASE_REPORT.relative_to(ROOT)), "width_ratio": 6.0, "mode_family": "CachedWidth(6.) + WideJointModes"},
        "scale_protocol": {"knots": KNOTS.tolist(), "bands": [[1.0, 3.0], [3.0, 5.0]], "k_definition": "k = -log2(2 tau)", "interpolation": "piecewise quintic smoothstep with zero first and second derivatives at knots; C2 across k=3 and completed earlier bands stay fixed", "sequential_rule": "stage n fits only the right knot a_(n+1) from the current left-knot state; previous knots are frozen", "active_coefficients_per_stage": NMODES, "increment_cap": 0.35, "global_amplitude_envelope": [-4.0, 4.0]},
        "initial_amplitudes": initial_amplitudes.tolist(),
        "final_knot_amplitudes": knots.tolist(),
        "stages": stages,
        "holdouts": holdouts,
        "volume_audit": volumes,
        "product_rule_fd_jets": product_rule,
        "interface_checks": interfaces,
        "aggregate": {"baseline_peak_max": float(max(base_peaks)), "sequential_peak_max": float(max(corrected_peaks)), "baseline_volume_L2_max": float(max(base_l2)), "sequential_volume_L2_max": float(max(corrected_l2)), "peak_max_improved": bool(max(corrected_peaks) < max(base_peaks)), "volume_L2_max_improved": bool(max(corrected_l2) < max(base_l2)), "peak_relative_change": float(max(corrected_peaks) / max(base_peaks) - 1.0), "volume_L2_relative_change": float(max(corrected_l2) / max(base_l2) - 1.0)},
        "runtime_seconds": time.time() - started,
        "scope": "Full Cartesian momentum is evaluated only on finite training/holdout point sets in the registered width-six annulus and eta slab. Divergence/interface checks are finite-difference diagnostics. No whole-space closure, continuum supremum, energy, or theorem transfer is asserted.",
        "pde_validated": False,
        "global_field_ready": False,
    }
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(output), "aggregate": report["aggregate"], "runtime_seconds": report["runtime_seconds"]}, indent=2), flush=True)
    return report


if __name__ == "__main__":
    run()
