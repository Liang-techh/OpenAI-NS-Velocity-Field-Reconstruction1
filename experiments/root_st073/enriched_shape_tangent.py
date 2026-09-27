"""Constrained 236-column shape tangent with a degree-3 mode-2 block.

The first 108 controls are the frozen degree-2 mode-0 and mode-1 controls.
The degree-2 mode-2 block (72 real columns) is replaced by the degree-3
mode-2 block (128 real columns), giving 236 controls.  Four moment equalities,
81 cone inequalities, and the three signed fractional shape-rate inequalities
are refit in this enlarged basis.  The instantaneous wave is unchanged.

This is an assembled local tangent diagnostic.  It does not provide a
finite-time replay, PDE, cone continuation, or scale-recursion acceptance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import LinearConstraint, linprog, minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from constrained_tangent_projection import ConstrainedTangent, quadratic_moment_target  # noqa: E402
from full_wave_tangent import _basis_columns  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402
from vortex_state_observables import fixed_cylinder  # noqa: E402
from wave_higher_harmonic_tangent import Mode2TangentCorrection, _mode_block_columns  # noqa: E402
from wave_momentum_projection import _metric, residual_and_jacobian  # noqa: E402
from wave_shape_tangent_rows import (  # noqa: E402
    _build_field,
    _curl,
    _metric_rate,
    _state_with_fd_jets,
)


CANDIDATE_PATH = ROOT / "wave_moment_cone_tangent.json"
SHAPE_PATH = ROOT / "wave_shape_tangent_rows.json"
MOMENT_PATH = ROOT / "wave_dynamics_mean_compatibility.json"
CONE_PATH = ROOT / "wave_mean_cone_projection.json"
CACHE_PATH = ROOT / "wave_momentum_projection.npz"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
OLD_PATH = ROOT / "shape_direction_margin_tangent.json"
WRAPPER_PATH = ROOT / "wave_higher_harmonic_tangent.py"
OUTPUT_PATH = ROOT / "enriched_shape_tangent.json"

FIXED_COLUMNS = 108
MODE2_DEGREE = 3
RELATIVE_SHAPE_RATE_FLOOR = 10.0
NU = 0.01


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _packed_complex(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _save(report: dict, path: Path) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _mode2_layout(degree: int):
    q = (int(degree) + 1) ** 2
    layout = []
    for kind, count in (("velocity", 3 * q), ("pressure", q)):
        for index in range(count):
            layout.append({"mode": 2, "kind": kind, "component": "real", "index": index})
            layout.append({"mode": 2, "kind": kind, "component": "imag", "index": index})
    return layout


def _mode2_shape_rows(points, weights, geometry, reference, dt_dk, degree):
    carriers = 2.0 * np.asarray(geometry["carrier"], dtype=float)
    values, gradients, _, _, _ = basis_jets(
        points,
        geometry["center"],
        geometry["widths"],
        2,
        int(degree),
        carriers,
        NU,
    )
    curls = _curl(gradients)
    rows = []
    q = (int(degree) + 1) ** 2
    for index in range(values.shape[-1]):
        rows.append(
            _metric_rate(
                values[:, :, index].real,
                curls[:, :, index].real,
                points,
                weights,
                reference,
                dt_dk,
            )
        )
        rows.append(
            _metric_rate(
                -values[:, :, index].imag,
                -curls[:, :, index].imag,
                points,
                weights,
                reference,
                dt_dk,
            )
        )
    # Pressure does not alter the shape observables.
    rows.extend([np.zeros(3, dtype=float) for _ in range(2 * q)])
    return np.asarray(rows, dtype=float).T


def _weighted_constraint_setup(design, moment_rows, target, residual, weights):
    sqrt_weights = np.repeat(np.sqrt(np.asarray(weights, dtype=float)), 3)
    weighted_design = np.asarray(design, dtype=float) * sqrt_weights[:, None]
    projector = ConstrainedTangent(weighted_design, moment_rows)
    particular = (projector.particular @ target) / projector.scales
    reduced = projector.design @ projector.nullspace
    u, singular, vh = np.linalg.svd(reduced, full_matrices=False)
    keep = singular > 1.0e-10 * singular[0]
    q = u[:, keep]
    mapping = ((projector.nullspace @ vh[keep].T) / singular[keep]) / projector.scales[:, None]
    shifted = np.asarray(residual, dtype=float).reshape(-1) * sqrt_weights
    shifted += weighted_design @ particular
    q0 = -q.T @ shifted
    scale = max(float(np.linalg.norm(shifted)), 1.0)
    return {
        "sqrt_weights": sqrt_weights,
        "weighted_design": weighted_design,
        "projector": projector,
        "particular": particular,
        "nullspace": projector.nullspace,
        "Q": q,
        "mapping": mapping,
        "q0": q0,
        "scale": scale,
        "singular_values": singular.tolist(),
        "keep_rank": int(np.sum(keep)),
    }


def run(relative_shape_rate=RELATIVE_SHAPE_RATE_FLOOR, output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    raw_paths = {
        "candidate": CANDIDATE_PATH,
        "shape": SHAPE_PATH,
        "moments": MOMENT_PATH,
        "cones": CONE_PATH,
        "projection_cache": CACHE_PATH,
        "frozen_grid": FROZEN_PATH,
        "wrapper": WRAPPER_PATH,
    }
    raw = {name: path.read_bytes() for name, path in raw_paths.items()}
    data = {
        "candidate": json.loads(raw["candidate"]),
        "shape": json.loads(raw["shape"]),
        "moments": json.loads(raw["moments"]),
        "cones": json.loads(raw["cones"]),
        "frozen": json.loads(raw["frozen_grid"]),
    }
    if any(data[name].get("status") != "completed" for name in ("candidate", "shape", "moments", "cones")):
        raise ValueError("All frozen source reports must be completed")
    candidate = data["candidate"]
    shape = data["shape"]
    moments = data["moments"]["reusable_moment_linearization"]
    cones = data["cones"].get("reusable_cone_linearization", data["cones"])
    if shape["source_sha256"] != hashlib.sha256(raw["candidate"]).hexdigest():
        raise ValueError("Shape rows do not describe the exact frozen candidate")

    with np.load(CACHE_PATH, allow_pickle=False) as loaded:
        cache = {key: loaded[key] for key in loaded.files}
    c = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    wave_x = np.r_[c[:, 0], c[:, 1]]
    residual, _ = residual_and_jacobian(cache, wave_x)
    points = np.asarray(cache["points"], dtype=float)
    weights = np.asarray(cache["weights"], dtype=float)
    geometry = data["frozen"]["inputs"]["wave"]
    if int(geometry["degree"]) != 2 or int(geometry["mode"]) != 1:
        raise ValueError("Frozen wave geometry is not the expected degree-2 mode-1 patch")
    carrier = np.asarray(geometry["carrier"], dtype=float)

    E = np.asarray(moments["moment_rows"], dtype=float)
    target, _ = quadratic_moment_target(
        moments["baseline_moments"], moments["wave_moment_forms"], wave_x
    )
    C = np.asarray(cones["cone_control_rows"], dtype=float)
    cone_base = np.asarray(cones["cone_baseline"], dtype=float) + np.einsum(
        "i,kij,j->k", wave_x, np.asarray(cones["cone_wave_forms"], dtype=float), wave_x
    )
    cone_lower = np.asarray(cones["cone_lower"], dtype=float)

    old_report = json.loads(OLD_PATH.read_text(encoding="utf-8"))
    old_selected = old_report.get("selected", {})
    S_old = np.asarray(shape["shape_control_rows"], dtype=float)
    shape_base = np.asarray(shape["shape_baseline"], dtype=float)
    signs = np.asarray(shape["desired_signs"], dtype=float)
    reference = np.asarray(shape["reference_observables"], dtype=float)
    shape_scale = np.maximum(np.abs(reference), 1.0e-30)
    if S_old.shape != (3, 180):
        raise ValueError(f"unexpected frozen shape rows {S_old.shape}")

    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "wave_integrated": False,
        "relative_shape_rate_floor": float(relative_shape_rate),
        "scope": (
            "Fixed instantaneous wave; first 108 degree-2 mode-0/mode-1 controls "
            "retained and degree-2 mode-2 controls replaced by degree-3 mode-2 "
            "controls. Four moments, 81 cones, and three signed shape-rate rows "
            "are assembled. No finite-time, PDE, or recursion acceptance."
        ),
        "sources": {
            name: {"path": path.name, "sha256": _sha256(path)}
            for name, path in raw_paths.items()
        },
        "inputs": {
            "training_point_count": int(len(points)),
            "training_weight_sum": float(np.sum(weights)),
            "wave_mode": int(geometry["mode"]),
            "wave_degree": int(geometry["degree"]),
            "center": list(geometry["center"]),
            "widths": list(geometry["widths"]),
            "carrier_mode1": carrier.tolist(),
            "carrier_mode2": (2.0 * carrier).tolist(),
            "old_control_count": 180,
            "fixed_control_count": FIXED_COLUMNS,
            "new_mode2_degree": MODE2_DEGREE,
            "new_mode2_q": (MODE2_DEGREE + 1) ** 2,
            "new_control_count": FIXED_COLUMNS + 2 * (3 + 1) * (MODE2_DEGREE + 1) ** 2,
            "layout": "first 108 frozen degree-2 mode-0/mode-1 columns, then degree-3 mode-2 velocity [Re,-Im] and pressure [Re,-Im] columns",
            "wave_coefficients_original": candidate["selected"]["coefficients_original"],
            "baseline_reference_observables": reference.tolist(),
            "desired_signs": signs.tolist(),
        },
        "old_degree2_comparison": {
            "path": OLD_PATH.name,
            "status": old_report.get("status"),
            "training_momentum": old_selected.get("training_momentum"),
            "assembled_moment_max_abs": old_selected.get("assembled_moment_max_abs"),
            "assembled_cone_min_margin": old_selected.get("assembled_cone_min_margin"),
            "assembled_shape_rates": old_selected.get("assembled_shape_rates"),
            "signed_relative_shape_rates": old_selected.get("signed_relative_shape_rates"),
        },
    }
    _save(report, output_path)
    print(json.dumps({"stage": "sources_loaded", "training_points": len(points)}), flush=True)

    # Recompute the reference shape state once, as requested, and derive only
    # the new mode-2 shape rows.  The frozen baseline/reference/signs remain
    # authoritative because the instantaneous velocity is unchanged.
    snapshot = data["frozen"]
    field, base, dynamic_report, field_geometry, _, mean_report = _build_field(candidate, snapshot)
    k0 = float(dynamic_report["k"])
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    shape_points, shape_weights, shape_domain = fixed_cylinder(base.inner, k0, angles=12)
    print(json.dumps({"stage": "shape_reference_started", "point_count": len(shape_points)}), flush=True)
    velocity_ref, omega_ref, density_ref, angular_speed_ref, reference_detail = _state_with_fd_jets(
        field, shape_points, shape_weights, tau0
    )
    reference_state = dict(reference_detail)
    reference_state.update({
        "density": density_ref,
        "omega": omega_ref,
        "angular_speed": angular_speed_ref,
    })
    reference_error = {
        key: float(reference_state[key] - shape["reference_observables_detail"][key])
        for key in (
            "enstrophy_radial_rms",
            "enstrophy_aspect_ratio",
            "enstrophy_weighted_angular_speed",
        )
    }
    dt_dk = float(shape["dt_dk"])
    S_mode2_degree2_recomputed = _mode2_shape_rows(
        shape_points, shape_weights, field_geometry, reference_state, dt_dk, 2
    )
    S_mode2_degree3 = _mode2_shape_rows(
        shape_points, shape_weights, field_geometry, reference_state, dt_dk, MODE2_DEGREE
    )
    mode2_old = S_old[:, FIXED_COLUMNS:]
    degree2_shape_relative_error = float(
        np.linalg.norm(S_mode2_degree2_recomputed - mode2_old)
        / max(np.linalg.norm(mode2_old), 1.0e-300)
    )
    S_enriched = np.concatenate((S_old[:, :FIXED_COLUMNS], S_mode2_degree3), axis=1)
    report["shape_rows"] = {
        "point_count": int(len(shape_points)),
        "domain": shape_domain,
        "reference_recomputed": {
            key: float(reference_state[key])
            for key in (
                "enstrophy_radial_rms",
                "enstrophy_aspect_ratio",
                "enstrophy_weighted_angular_speed",
            )
        },
        "reference_error_vs_frozen": reference_error,
        "degree2_mode2_row_relative_error_vs_frozen": degree2_shape_relative_error,
        "old_rows_shape": list(S_old.shape),
        "new_mode2_rows_shape": list(S_mode2_degree3.shape),
        "new_rows_shape": list(S_enriched.shape),
        "dt_dk": dt_dk,
        "baseline_source": "wave_shape_tangent_rows.json unchanged because instantaneous velocity is identical",
        "rows_mode2_degree3": S_mode2_degree3.tolist(),
        "control_layout": [
            *shape["control_layout"][:FIXED_COLUMNS],
            *_mode2_layout(MODE2_DEGREE),
        ],
    }
    if len(report["shape_rows"]["control_layout"]) != FIXED_COLUMNS + S_mode2_degree3.shape[1]:
        raise ValueError("enriched shape control layout length mismatch")
    report["status"] = "shape_rows_assembled"
    _save(report, output_path)
    print(json.dumps({
        "stage": "shape_rows_assembled",
        "new_rows": list(S_enriched.shape),
        "degree2_row_relative_error": degree2_shape_relative_error,
        "reference_error": reference_error,
    }), flush=True)

    # Build the 236-column momentum operator.  The first 108 columns are
    # exactly the cached degree-2 operator; only the mode-2 block is enriched.
    D_fixed = np.asarray(cache["tangent_design"], dtype=float)[:, :FIXED_COLUMNS]
    D_mode2_degree3 = _mode_block_columns(
        points,
        field_geometry["center"],
        field_geometry["widths"],
        2,
        MODE2_DEGREE,
        2.0 * carrier,
    )
    D = np.concatenate((D_fixed, D_mode2_degree3), axis=1)
    ncontrol = D.shape[1]
    if ncontrol != 236:
        raise ValueError(f"expected 236 enriched controls, got {ncontrol}")
    zeros_pad = np.zeros((E.shape[0], ncontrol - E.shape[1]))
    E_enriched = np.concatenate((E, zeros_pad), axis=1)
    C_enriched = np.concatenate((C, np.zeros((C.shape[0], ncontrol - C.shape[1]))), axis=1)
    target_setup = _weighted_constraint_setup(
        D, E_enriched, target, residual, weights
    )
    S = S_enriched
    C_all = np.vstack((C_enriched, signs[:, None] * S / shape_scale[:, None]))
    base_all = np.r_[cone_base, signs * shape_base / shape_scale]
    lower_all = np.r_[cone_lower, np.full(3, float(relative_shape_rate))]
    mapping = target_setup["mapping"]
    scale = target_setup["scale"]
    matrix = C_all @ mapping * scale
    row_norms = np.linalg.norm(matrix, axis=1)
    responsive = row_norms > max(row_norms.max() * 1.0e-12, 1.0e-14)
    safety = np.r_[np.full(len(C_enriched), 1.0e-4), np.zeros(3)] * responsive
    rhs = lower_all + safety - base_all - C_all @ target_setup["particular"]
    report["operator"] = {
        "control_count": ncontrol,
        "moment_rows_shape": list(E_enriched.shape),
        "cone_rows_shape": list(C_enriched.shape),
        "shape_rows_shape": list(S.shape),
        "moment_padding": "first 108 cached columns retained; 128 degree-3 mode-2 columns are zero in moment rows",
        "cone_padding": "first 108 cached columns retained; 128 degree-3 mode-2 columns are zero in cone rows",
        "moment_rank": target_setup["projector"].metadata["moment_rank"],
        "objective_rank": target_setup["keep_rank"],
        "null_design_singular_values": target_setup["singular_values"],
        "responsive_constraint_count": int(np.sum(responsive)),
        "fixed_constraint_count": int(np.sum(~responsive)),
    }
    report["status"] = "qp_assembled"
    _save(report, output_path)
    print(json.dumps({
        "stage": "qp_assembled",
        "controls": ncontrol,
        "objective_rank": target_setup["keep_rank"],
        "constraints": int(len(rhs)),
        "responsive": int(np.sum(responsive)),
    }), flush=True)

    if np.any(rhs[~responsive] > 1.0e-7):
        report.update(
            status="fixed_row_failure",
            fixed_row_rhs=rhs[~responsive].tolist(),
            assembled_feasible=False,
        )
        _save(report, output_path)
        return report

    matrix_n = matrix[responsive] / row_norms[responsive, None]
    rhs_n = rhs[responsive] / row_norms[responsive]
    z0 = target_setup["q0"] / scale
    phase = linprog(
        np.r_[np.zeros(len(z0)), 1.0],
        A_ub=np.column_stack((-matrix_n, -np.ones(len(rhs_n)))),
        b_ub=-rhs_n,
        bounds=[(None, None)] * len(z0) + [(0, None)],
        method="highs",
    )
    report["phase_one"] = {
        "success": bool(phase.success),
        "status": int(phase.status),
        "message": str(phase.message),
        "slack": float(phase.x[-1]) if phase.x is not None else None,
        "failure_classification": (
            "linear_program_infeasible" if phase.status == 2 else "solver_failure_or_unresolved"
        ) if not phase.success else "feasible",
    }
    if not phase.success or phase.x[-1] > 1.0e-8:
        report.update(
            status="restricted_basis_infeasible_or_unresolved",
            assembled_feasible=False,
        )
        _save(report, output_path)
        print(json.dumps(report["phase_one"]), flush=True)
        return report

    fit = minimize(
        lambda z: 0.5 * np.sum((z - z0) ** 2),
        phase.x[:-1],
        jac=lambda z: z - z0,
        method="SLSQP",
        constraints=[LinearConstraint(matrix_n, rhs_n, np.inf)],
        options=dict(maxiter=500, ftol=1.0e-12),
    )
    z = phase.x[:-1]
    if np.min(matrix_n @ fit.x - rhs_n) >= -1.0e-9 and np.sum((fit.x - z0) ** 2) <= np.sum((z - z0) ** 2):
        z = fit.x
    control = target_setup["particular"] + mapping @ (scale * z)
    corrected = residual + (D @ control).reshape(-1, 3)
    cone_margins = cone_base + C_enriched @ control - cone_lower
    shape_rates = shape_base + S @ control
    signed_shape = signs * shape_rates / shape_scale
    report["optimizer"] = {
        "success": bool(fit.success),
        "status": int(fit.status),
        "message": str(fit.message),
        "selected_phase_one_or_slsqp": "slsqp" if np.allclose(z, fit.x) else "phase_one",
    }
    report["selected"] = {
        "coefficients_original": candidate["selected"]["coefficients_original"],
        "coefficients_whitened": candidate["selected"]["coefficients_whitened"],
        "tangent_coefficients": control.tolist(),
        "tangent_coefficient_count": int(len(control)),
        "tangent_layout": report["shape_rows"]["control_layout"],
        "mode2_degree3_block": control[FIXED_COLUMNS:].tolist(),
        "mode2_degree3_block_complex": None,
        "growth_lambda": candidate["selected"]["growth_lambda"],
        "training_momentum": _metric(corrected, weights),
        "assembled_moment_max_abs": float(np.max(np.abs(E_enriched @ control - target))),
        "assembled_cone_min_margin": float(np.min(cone_margins)),
        "assembled_cone_location_pass_count": int(np.sum(np.all(cone_margins.reshape(-1, 3) >= -1.0e-7, axis=1))),
        "assembled_cone_inequality_pass_count": int(np.sum(cone_margins >= -1.0e-7)),
        "assembled_shape_rates": shape_rates.tolist(),
        "signed_relative_shape_rates": signed_shape.tolist(),
        "relative_shape_floor_min_margin": float(np.min(signed_shape - relative_shape_rate)),
    }
    # Decode the selected degree-3 mode-2 block without passing it to any
    # legacy 180-column replay.
    block = control[FIXED_COLUMNS:]
    q3 = (MODE2_DEGREE + 1) ** 2
    velocity_count = 3 * q3
    velocity_block = block[:2 * velocity_count:2] + 1j * block[1:2 * velocity_count:2]
    pressure_start = 2 * velocity_count
    pressure_block = (
        block[pressure_start:pressure_start + 2 * q3:2]
        + 1j * block[pressure_start + 1:pressure_start + 2 * q3:2]
    )
    report["selected"]["mode2_degree3_block_complex"] = {
        "velocity": _packed_complex(velocity_block),
        "pressure": _packed_complex(pressure_block),
    }
    report["selected"]["assembled_feasible"] = bool(
        np.max(np.abs(E_enriched @ control - target)) < 1.0e-5
        and np.min(cone_margins) >= -1.0e-7
        and np.min(signed_shape - relative_shape_rate) >= -1.0e-7
    )
    report["selected"]["legacy_180_replay_compatible"] = False
    report["selected"]["wrapper"] = {
        "class": "Mode2TangentCorrection",
        "module": WRAPPER_PATH.name,
        "degree": MODE2_DEGREE,
        "base_field_requirement": "fixed wave plus first-108 mode-0/mode-1 controls",
    }
    report["assembled_feasible"] = report["selected"]["assembled_feasible"]
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(json.dumps({
        "stage": "completed",
        "assembled_feasible": report["assembled_feasible"],
        "training_l2": report["selected"]["training_momentum"]["momentum_volume_L2"],
        "moment_error": report["selected"]["assembled_moment_max_abs"],
        "cone_min_margin": report["selected"]["assembled_cone_min_margin"],
        "shape_min_margin": report["selected"]["relative_shape_floor_min_margin"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--relative-shape-rate", type=float, default=RELATIVE_SHAPE_RATE_FLOOR)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.relative_shape_rate, args.output)
