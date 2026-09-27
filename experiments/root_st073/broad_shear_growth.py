"""Bounded scalar screen for a broad annular mean shear.

This experiment reuses the factor-8 degree-2 Fourier energy matrices from
``fourier_shear_codesign.json`` and adds one smooth axisymmetric shear
direction.  It searches only one scalar amplitude, with the sampled
lambda-squared cone at the 22 solve-control locations as a feasibility
screen.  The result is a candidate diagnostic; it is not a moment, PDE, or
Navier--Stokes acceptance test.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_annular_shear import BroadAnnularShear  # noqa: E402
from fourier_patch_evolution import basis_jets  # noqa: E402
from joined_field import coordinates  # noqa: E402
from midplane_resolved_feasibility import ZeroBackground  # noqa: E402
from outer_feedback_evolution import Trajectory, build_current  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


INPUT_PATH = ROOT / "fourier_patch_implicit.json"
MATRIX_PATH = ROOT / "fourier_shear_codesign.json"
MEAN_PATH = ROOT / "outer_feedback_evolution.json"
OUTPUT_PATH = ROOT / "broad_shear_growth.json"
TARGET_GROWTH = 1.0e4
WIDTH_FACTOR = 8.0
QUADRATURE_ORDER = 16
REFINED_PANEL_ORDER = 12
REFINED_Z_ORDER = 16
CONe_FLOOR_FACTOR = 0.01


def decode_complex(value):
    value = np.asarray(value, dtype=float)
    return value[..., 0] + 1j * value[..., 1]


def pack_complex(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def hermitian(matrix):
    return 0.5 * (matrix + matrix.conj().T)


def spatial_gradient(field, points, tau, h):
    """Fourth-order Cartesian gradient using only spatial field calls."""

    points = np.asarray(points, dtype=float)
    gradient = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        direction = np.zeros(3)
        direction[axis] = h
        minus2 = field.fields(points - 2.0 * direction, tau)[0]
        minus1 = field.fields(points - direction, tau)[0]
        plus1 = field.fields(points + direction, tau)[0]
        plus2 = field.fields(points + 2.0 * direction, tau)[0]
        gradient[:, :, axis] = (minus2 - 8.0 * minus1 + 8.0 * plus1 - plus2) / (12.0 * h)
    return gradient


def node_points(inner, tau):
    locations = list(dict.fromkeys(
        [(eta, y) for eta in (-0.2, 0.0, 0.2) for y in (0.5, 0.75)]
        + [(eta + de, 0.75 + dy) for eta in (-0.2, 0.2)
           for de in (-0.005, 0.0, 0.005)
           for dy in (-0.005, 0.0, 0.005)]
    ))
    points = []
    labels = []
    for eta, y in locations:
        X = inner.p.X_max * (1.0 + 15.0 * y) ** 2
        points.append(inner.from_similarity([X], [eta], tau)[0])
        labels.append({"eta": eta, "y": y})
    return np.asarray(points), labels


def cone_lambda2(points, velocity, gradient):
    """Sample the lambda-squared cone diagnostic at theta equal to zero."""

    values = []
    for point, value, jacobian in zip(points, velocity, gradient):
        radius = float(np.hypot(point[0], point[1]))
        safe = max(radius, np.finfo(float).tiny)
        ca, sa = point[0] / safe, point[1] / safe
        u_theta = -sa * value[0] + ca * value[1]
        j_theta_r = -sa * jacobian[0, 0] + ca * jacobian[1, 0]
        swirl = u_theta / safe
        shear = np.array([j_theta_r - swirl, jacobian[2, 0]], dtype=float)
        shear_norm = float(np.linalg.norm(shear))
        if shear_norm == 0.0:
            values.append(float("nan"))
            continue
        normal = shear / shear_norm
        values.append(float(-2.0 * swirl * normal[0]
                            * (2.0 * swirl * normal[0] + shear_norm)))
    return np.asarray(values, dtype=float)


def candidate_node_lambda(node_cache, amplitude):
    velocity = node_cache["baseline_velocity"] + amplitude * node_cache["unit_velocity"]
    gradient = node_cache["baseline_gradient"] + amplitude * node_cache["unit_gradient"]
    return cone_lambda2(node_cache["points"], velocity, gradient)


def evaluate_amplitude(amplitude, mode_data, node_cache, cone_floor):
    growths = {}
    vectors = {}
    totals = {}
    for mode, item in mode_data.items():
        total = hermitian(item["H0"] + amplitude * item["H_broad"])
        values, vectors_mode = np.linalg.eigh(total)
        growths[str(mode)] = float(values[-1])
        vectors[str(mode)] = vectors_mode[:, -1]
        totals[str(mode)] = total
    best_mode = max(growths, key=growths.get)
    node_lambda = candidate_node_lambda(node_cache, amplitude)
    finite_cone = bool(np.all(np.isfinite(node_lambda)))
    cone_ok = finite_cone and bool(np.all(node_lambda >= cone_floor))
    growth_ok = bool(growths[best_mode] >= TARGET_GROWTH)
    return {
        "amplitude": float(amplitude),
        "growths": growths,
        "best_mode": int(best_mode),
        "best_growth_lambda": float(growths[best_mode]),
        "node_lambda2": node_lambda,
        "min_node_lambda2": float(np.nanmin(node_lambda)),
        "cone_ok": cone_ok,
        "growth_ok": growth_ok,
        "feasible": bool(cone_ok and growth_ok),
        "vectors": vectors,
        "totals": totals,
    }


def scan_sign(sign, mode_data, node_cache, cone_floor):
    """Find the first feasible point on a signed geometric amplitude ray."""

    magnitudes = np.geomspace(1.0e-4, 1.0e7, 32)
    records = []
    previous = evaluate_amplitude(0.0, mode_data, node_cache, cone_floor)
    lower = 0.0
    upper = None
    for magnitude in magnitudes:
        amplitude = float(sign * magnitude)
        current = evaluate_amplitude(amplitude, mode_data, node_cache, cone_floor)
        records.append({
            "amplitude": current["amplitude"],
            "best_mode": current["best_mode"],
            "best_growth_lambda": current["best_growth_lambda"],
            "min_node_lambda2": current["min_node_lambda2"],
            "growth_ok": current["growth_ok"],
            "cone_ok": current["cone_ok"],
            "feasible": current["feasible"],
        })
        if current["feasible"]:
            upper = amplitude
            break
        lower = amplitude
        previous = current
    if upper is None:
        return {
            "sign": int(sign),
            "bracket_found": False,
            "records": records,
            "best_attempt": max(
                [previous] + [evaluate_amplitude(r["amplitude"], mode_data, node_cache, cone_floor)
                               for r in records],
                key=lambda item: (item["cone_ok"], item["min_node_lambda2"],
                                  item["best_growth_lambda"]),
            ),
        }

    # The search has a feasible endpoint and a non-feasible endpoint.  Refine
    # the first observed crossing; retain a monotonicity flag because the
    # cone and growth conditions are not mathematically guaranteed monotone.
    for _ in range(32):
        midpoint = 0.5 * (lower + upper)
        trial = evaluate_amplitude(midpoint, mode_data, node_cache, cone_floor)
        if trial["feasible"]:
            upper = midpoint
        else:
            lower = midpoint
    selected = evaluate_amplitude(upper, mode_data, node_cache, cone_floor)
    return {
        "sign": int(sign),
        "bracket_found": True,
        "lower_nonfeasible": float(lower),
        "upper_feasible": float(upper),
        "records": records,
        "selected": selected,
    }


def y_of_r(r, z, tau, inner):
    co = coordinates(r / np.sqrt(inner.nu), z / np.sqrt(inner.nu), tau, inner.inner.h)
    q = float(np.asarray(co["q"]).reshape(-1)[0])
    ri = np.sqrt(2.0 * inner.nu * q * inner.join_X)
    return float((r / ri - 1.0) / (inner.ratio - 1.0))


def refined_quadrature(center, widths, tau, inner, z_order=REFINED_Z_ORDER,
                       panel_order=REFINED_PANEL_ORDER):
    """Tensor quadrature with panels at the four physical cutoff roots."""

    z_nodes, z_weights = leggauss(z_order)
    points = []
    weights = []
    rlo = float(center[0] - widths[0])
    rhi = float(center[0] + widths[0])
    cut_values = (0.01, 0.10, 0.88, 0.99)
    for z_node, z_weight in zip(z_nodes, z_weights):
        z = float(center[1] + z_node * widths[1])
        cuts = [rlo, rhi]
        for target in cut_values:
            f_lo = y_of_r(rlo, z, tau, inner) - target
            f_hi = y_of_r(rhi, z, tau, inner) - target
            if f_lo == 0.0:
                cuts.append(rlo)
            elif f_hi == 0.0:
                cuts.append(rhi)
            elif f_lo * f_hi < 0.0:
                cuts.append(brentq(lambda r: y_of_r(r, z, tau, inner) - target,
                                   rlo, rhi, xtol=1.0e-14, rtol=1.0e-12))
        cuts = sorted(set(float(value) for value in cuts))
        r_nodes, r_weights = leggauss(panel_order)
        for left, right in zip(cuts[:-1], cuts[1:]):
            midpoint = 0.5 * (left + right)
            half = 0.5 * (right - left)
            for r_node, r_weight in zip(r_nodes, r_weights):
                r = midpoint + half * r_node
                points.append((r, 0.0, z))
                weights.append(z_weight * half * r_weight * r)
    points = np.asarray(points, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights /= np.sum(weights)
    return points, weights, {
        "z_order": int(z_order),
        "panel_order": int(panel_order),
        "radial_bounds": [rlo, rhi],
        "point_count": int(len(points)),
        "panels_per_z_max": int(len(cut_values) + 1),
    }


def refined_mode_matrix(points, weights, center, widths, degree, carrier, mode,
                        baseline_gradient, broad_gradient, nu, h, amplitude):
    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, mode, degree, carrier, nu, h
    )
    q = velocity.shape[-1]
    velocity_flat = velocity.reshape(-1, q)
    gradient_flat = gradient.reshape(-1, q)
    gradient_weight = np.repeat(weights, 9)
    stiffness = gradient_flat.conj().T @ (gradient_weight[:, None] * gradient_flat)
    strain0 = 0.5 * (baseline_gradient + np.swapaxes(baseline_gradient, 1, 2))
    strainb = 0.5 * (broad_gradient + np.swapaxes(broad_gradient, 1, 2))
    matrix0 = np.einsum("nia,nij,njb,n->ab", velocity.conj(), strain0, velocity, weights)
    matrixb = np.einsum("nia,nij,njb,n->ab", velocity.conj(), strainb, velocity, weights)
    weighted_velocity = velocity_flat * np.sqrt(np.repeat(weights, 3))[:, None]
    _, singular, right = np.linalg.svd(weighted_velocity, full_matrices=False)
    cutoff = max(float(singular[0]), np.finfo(float).tiny) * 1.0e-10
    keep = singular > cutoff
    whitening = right[keep].conj().T / singular[keep][None, :]
    h0 = hermitian(-nu * whitening.conj().T @ stiffness @ whitening
                   - whitening.conj().T @ matrix0 @ whitening)
    hb = hermitian(-whitening.conj().T @ matrixb @ whitening)
    total = hermitian(h0 + amplitude * hb)
    return float(np.linalg.eigvalsh(total)[-1]), int(np.sum(keep)), float(singular[0])


def run(input_path=INPUT_PATH, matrix_path=MATRIX_PATH, mean_path=MEAN_PATH,
        output_path=OUTPUT_PATH):
    started = time.perf_counter()
    report = json.loads(Path(input_path).read_text(encoding="utf-8"))
    matrix_report = json.loads(Path(matrix_path).read_text(encoding="utf-8"))
    if not report.get("state_data"):
        print("state_data absent; no broad-shear report written", flush=True)
        return None
    state_data = report["state_data"]
    tau = float(-state_data["physical_times"][0])
    center = tuple(float(value) for value in matrix_report["center"])
    widths = tuple(float(value) for value in matrix_report["widths"])
    degree = int(matrix_report["degree"])
    modes = [int(value) for value in matrix_report["modes"] if int(value) >= 1]
    carriers = {int(mode): np.asarray(value, dtype=float)
                for mode, value in matrix_report["carriers"].items()}
    _, base, current = build_current()
    mean_report = json.loads(Path(mean_path).read_text(encoding="utf-8"))
    mean = Trajectory(current, mean_report["nodes"])
    h = 0.0005 * np.sqrt(mean.nu * tau)
    points = np.asarray(matrix_report["matrix_cache"]["points"], dtype=float)
    weights = np.asarray(matrix_report["matrix_cache"]["weights"], dtype=float)
    if center[0] - widths[0] <= 0.0:
        raise ValueError("factor-8 radial support reaches the cylindrical axis")
    print(json.dumps({"stage": "mean_ready", "tau": tau,
                      "point_count": len(points)}), flush=True)

    zero = ZeroBackground(base)
    broad_unit = BroadAnnularShear(zero, 1.0)
    broad_velocity = broad_unit.fields(points, tau)[0]
    broad_gradient = spatial_gradient(broad_unit, points, tau, h)
    print(json.dumps({"stage": "broad_factor8_ready", "point_count": len(points)}), flush=True)

    node_pts, node_labels = node_points(current.inner, tau)
    node_cache = {
        "points": node_pts,
        "labels": node_labels,
        "baseline_velocity": mean.fields(node_pts, tau)[0],
        "baseline_gradient": spatial_gradient(mean, node_pts, tau, h),
        "unit_velocity": broad_unit.fields(node_pts, tau)[0],
        "unit_gradient": spatial_gradient(broad_unit, node_pts, tau, h),
    }
    baseline_node_lambda = cone_lambda2(
        node_pts, node_cache["baseline_velocity"], node_cache["baseline_gradient"]
    )
    cone_floor = CONe_FLOOR_FACTOR * baseline_node_lambda
    if not np.all(np.isfinite(baseline_node_lambda)):
        raise ValueError("baseline node cone diagnostic contains non-finite values")
    print(json.dumps({"stage": "node_cache_ready", "nodes": len(node_pts),
                      "min_baseline_lambda2": float(np.min(baseline_node_lambda))}), flush=True)

    mode_data = {}
    cache_modes = matrix_report["matrix_cache"]["modes"]
    for mode in modes:
        item = cache_modes[str(mode)]
        q = decode_complex(item["whitening"])
        h0 = decode_complex(item["H0"])
        velocity, _, _ = basis_data(
            points, center=center, widths=widths, mode=mode, degree=degree,
            carrier=carriers[mode]
        )
        strainb = 0.5 * (broad_gradient + np.swapaxes(broad_gradient, 1, 2))
        broad_strain_matrix = np.einsum(
            "nia,nij,njb,n->ab", velocity.conj(), strainb, velocity, weights
        )
        h_broad = hermitian(-q.conj().T @ broad_strain_matrix @ q)
        mode_data[str(mode)] = {
            "H0": h0,
            "H_broad": h_broad,
            "whitening": q,
            "dimension": int(q.shape[0]),
        }
        print(json.dumps({"stage": "mode_matrix_ready", "mode": mode,
                          "rank": int(q.shape[1])}), flush=True)

    scans = [scan_sign(sign, mode_data, node_cache, cone_floor) for sign in (1, -1)]
    feasible_scans = [item for item in scans if item.get("bracket_found")]
    selected_scan = None
    if feasible_scans:
        selected_scan = min(
            feasible_scans,
            key=lambda item: abs(item["selected"]["amplitude"]),
        )
    best_attempt = None
    if selected_scan is None:
        attempts = []
        for scan in scans:
            attempts.append(scan.get("best_attempt"))
        attempts = [item for item in attempts if item is not None]
        if attempts:
            best_attempt = max(
                attempts,
                key=lambda item: (bool(item["cone_ok"]), item["min_node_lambda2"],
                                  item["best_growth_lambda"]),
            )
    print(json.dumps({"stage": "scalar_search_done",
                      "selected": selected_scan is not None,
                      "selected_amplitude": None if selected_scan is None else selected_scan["selected"]["amplitude"]}), flush=True)

    selected = None
    refined = None
    higher_order = None
    if selected_scan is not None:
        selected_eval = selected_scan["selected"]
        amplitude = float(selected_eval["amplitude"])
        mode = int(selected_eval["best_mode"])
        vector = selected_eval["vectors"][str(mode)]
        potential_coefficients = mode_data[str(mode)]["whitening"] @ vector
        added_velocity = amplitude * broad_velocity
        baseline_velocity = mean.fields(points, tau)[0]
        selected = {
            "delta_state": [0.0] * 9,
            "mode": mode,
            "amplitude": amplitude,
            "growth_lambda": float(selected_eval["best_growth_lambda"]),
            "energy_rate_2lambda": 2.0 * float(selected_eval["best_growth_lambda"]),
            "growths_by_mode": selected_eval["growths"],
            "node_lambda2": selected_eval["node_lambda2"].tolist(),
            "min_node_lambda2": float(selected_eval["min_node_lambda2"]),
            "mean_increment": {
                "type": "broad_annular_shear",
                "parameters": BroadAnnularShear(zero, amplitude).parameters(),
            },
            "potential_coefficients": pack_complex(potential_coefficients),
            "field_coefficient": pack_complex(potential_coefficients),
            "potential_coefficient_dimension": int(len(potential_coefficients)),
            "added_swirl_rms": float(np.sqrt(np.sum(
                weights * np.sum(added_velocity ** 2, axis=1)
            ))),
            "added_swirl_peak": float(np.max(np.linalg.norm(added_velocity, axis=1))),
            "baseline_velocity_rms": float(np.sqrt(np.sum(
                weights * np.sum(baseline_velocity ** 2, axis=1)
            ))),
            "baseline_velocity_peak": float(np.max(np.linalg.norm(baseline_velocity, axis=1))),
        }

        print(json.dumps({"stage": "refined_quadrature_start", "amplitude": amplitude}), flush=True)
        refined_points, refined_weights, refined_meta = refined_quadrature(
            center, widths, tau, broad_unit
        )
        baseline_gradient_refined = spatial_gradient(mean, refined_points, tau, h)
        broad_gradient_refined = spatial_gradient(broad_unit, refined_points, tau, h)
        refined_growths = {}
        refined_ranks = {}
        for mode_item in modes:
            refined_growths[str(mode_item)], refined_ranks[str(mode_item)], _ = refined_mode_matrix(
                refined_points, refined_weights, center, widths, degree,
                carriers[mode_item], mode_item, baseline_gradient_refined,
                broad_gradient_refined, mean.nu, h, amplitude
            )
        refined_best_mode = max(refined_growths, key=refined_growths.get)
        refined_growth = float(refined_growths[refined_best_mode])
        fixed_growth = float(selected["growth_lambda"])
        relative_difference = abs(refined_growth - fixed_growth) / max(abs(fixed_growth), 1.0)
        refined = {
            **refined_meta,
            "growths_by_mode": refined_growths,
            "best_mode": int(refined_best_mode),
            "best_growth_lambda": refined_growth,
            "fixed_selected_mode_growth_lambda": float(refined_growths[str(mode)]),
            "factor8_growth_lambda": fixed_growth,
            "relative_difference": float(relative_difference),
            "growth_discrepancy_status": "unresolved" if relative_difference > 0.02 else "consistent",
            "ranks": refined_ranks,
        }
        print(json.dumps({"stage": "refined_quadrature_done", "point_count": len(refined_points),
                          "best_growth_lambda": refined_growth,
                          "relative_difference": relative_difference}), flush=True)

        # A fixed-amplitude, higher-panel-order replay checks whether the
        # cutoff split itself has converged.  This is deliberately one more
        # spatial screen; it does not launch a time integration.
        print(json.dumps({"stage": "higher_panel_replay_start",
                          "panel_order": 16, "amplitude": amplitude}), flush=True)
        high_points, high_weights, high_meta = refined_quadrature(
            center, widths, tau, broad_unit, panel_order=16
        )
        baseline_gradient_high = spatial_gradient(mean, high_points, tau, h)
        broad_gradient_high = spatial_gradient(broad_unit, high_points, tau, h)
        high_growths = {}
        high_ranks = {}
        for mode_item in modes:
            high_growths[str(mode_item)], high_ranks[str(mode_item)], _ = refined_mode_matrix(
                high_points, high_weights, center, widths, degree,
                carriers[mode_item], mode_item, baseline_gradient_high,
                broad_gradient_high, mean.nu, h, amplitude
            )
        high_best_mode = max(high_growths, key=high_growths.get)
        high_growth = float(high_growths[high_best_mode])
        selected_mode_high_growth = float(high_growths[str(mode)])
        high_vs_split = abs(high_growth - refined_growth) / max(abs(refined_growth), 1.0)
        high_vs_coarse = abs(high_growth - fixed_growth) / max(abs(fixed_growth), 1.0)
        higher_order = {
            **high_meta,
            "growths_by_mode": high_growths,
            "best_mode": int(high_best_mode),
            "best_growth_lambda": high_growth,
            "fixed_selected_mode_growth_lambda": selected_mode_high_growth,
            "factor8_growth_lambda": fixed_growth,
            "relative_difference_from_panel12": float(high_vs_split),
            "relative_difference_from_factor8": float(high_vs_coarse),
            "growth_convergence_status": "consistent_with_split" if high_vs_split <= 0.02 else "unresolved",
            "ranks": high_ranks,
        }
        selected["resolved_replay_growth_lambda"] = selected_mode_high_growth
        selected["resolved_replay_passes_growth_gate"] = bool(selected_mode_high_growth >= TARGET_GROWTH)
        selected["resolved_replay_passes_22_node_gate"] = bool(
            np.all(selected_eval["node_lambda2"] >= cone_floor)
        )
        print(json.dumps({"stage": "higher_panel_replay_done",
                          "point_count": len(high_points),
                          "best_growth_lambda": high_growth,
                          "relative_difference_from_panel12": high_vs_split,
                          "relative_difference_from_factor8": high_vs_coarse}), flush=True)

    baseline_matrix = {}
    for mode, item in mode_data.items():
        baseline_matrix[mode] = {
            "lambda_max": float(np.linalg.eigvalsh(item["H0"])[-1]),
            "lambda_min": float(np.linalg.eigvalsh(item["H0"])[0]),
        }
    result = {
        "accepted": False,
        "pde_validated": False,
        "scope": "Bounded one-scalar broad-annular shear growth screen only; no moment feasibility, solve_control, time integration, nonlinear transfer, or Navier--Stokes acceptance.",
        "source_report": "fourier_patch_implicit.json",
        "matrix_source": "fourier_shear_codesign.json",
        "mean_report": "outer_feedback_evolution.json",
        "initial_k": float(-np.log2(2.0 * tau)),
        "initial_physical_time": -tau,
        "center": list(center),
        "widths": list(widths),
        "width_factor": WIDTH_FACTOR,
        "degree": degree,
        "modes": modes,
        "carriers": {str(mode): carriers[mode].tolist() for mode in modes},
        "target_growth_lambda": TARGET_GROWTH,
        "cone_floor_factor": CONe_FLOOR_FACTOR,
        "broad_parameters_unit": BroadAnnularShear(zero, 1.0).parameters(),
        "baseline_diagnostic": baseline_matrix,
        "baseline_node_lambda2": baseline_node_lambda.tolist(),
        "node_labels": node_labels,
        "search": {
            "method": "signed geometric bracket followed by bisection",
            "signs": [1, -1],
            "scans": [
                {key: value for key, value in scan.items()
                 if key not in ("selected", "best_attempt")}
                for scan in scans
            ],
        },
        "best_attempt": None if best_attempt is None else {
            key: value for key, value in best_attempt.items()
            if key not in ("vectors", "totals")
        },
        "selected_candidate": selected,
        "refined_replay": refined,
        "higher_order_fixed_amplitude": higher_order,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    Path(output_path).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "saved", "output": str(output_path),
                      "selected": selected is not None}), flush=True)
    return result


if __name__ == "__main__":
    run()
