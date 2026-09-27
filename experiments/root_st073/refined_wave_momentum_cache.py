"""Higher-quadrature 264-column cache for the current enriched wave candidate.

The spatial grid is a fresh order-20, 12-angle split grid.  The frozen base
residual is evaluated once with the current initial wave and no tangent
controls using the existing five-point finite-difference jet.  The 264-column
tangent matrix contains degree-3 real mode 0, degree-2 complex mode 1, and
degree-3 complex mode 2 responses.  No optimization or constraint refit is
performed here.
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


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field  # noqa: E402
from full_wave_tangent import LocalPotentialField, _grid  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from wave_higher_harmonic_tangent import _mode_block_columns  # noqa: E402
from wave_momentum_projection import _metric  # noqa: E402


FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
MEAN_PATH = ROOT / "broad_meridional_constrained.json"
OUTPUT_PATH = ROOT / "refined_wave_momentum_cache.json"
NPZ_PATH = ROOT / "refined_wave_momentum_cache.npz"

ORDER = 20
ANGLES = 12
SHIFT = 0.23
MODE0_DEGREE = 3
MODE1_DEGREE = 2
MODE2_DEGREE = 3


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict) -> None:
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _packed_complex(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1)


def _unpack_packed(value):
    value = np.asarray(value, dtype=float)
    if value.shape != (27, 2):
        raise ValueError(f"Expected packed current wave shape (27, 2), got {value.shape}")
    return value[:, 0] + 1j * value[:, 1]


def _max_residual_detail(residual, points, weights):
    residual = np.asarray(residual, dtype=float)
    norms = np.linalg.norm(residual, axis=1)
    index = int(np.argmax(norms))
    return {
        "index": index,
        "point": np.asarray(points[index], dtype=float).tolist(),
        "weight": float(weights[index]),
        "vector": residual[index].tolist(),
        "norm": float(norms[index]),
    }


def _mode0_real_columns(points, geometry):
    """Degree-3 mode-0 real velocity/pressure columns in canonical order."""

    full = _mode_block_columns(
        points,
        geometry["center"],
        geometry["widths"],
        0,
        MODE0_DEGREE,
        np.zeros(2),
    )
    return full[:, ::2]


def _mixed_tangent_design(points, geometry):
    mode0 = _mode0_real_columns(points, geometry)
    mode1 = _mode_block_columns(
        points,
        geometry["center"],
        geometry["widths"],
        1,
        MODE1_DEGREE,
        np.asarray(geometry["carrier"], dtype=float),
    )
    mode2 = _mode_block_columns(
        points,
        geometry["center"],
        geometry["widths"],
        2,
        MODE2_DEGREE,
        2.0 * np.asarray(geometry["carrier"], dtype=float),
    )
    design = np.column_stack((mode0, mode1, mode2))
    if design.shape[1] != 264:
        raise ValueError(f"Expected 264 tangent columns, got {design.shape}")
    return design, {
        "mode0_degree3_real": list(mode0.shape),
        "mode1_degree2_complex": list(mode1.shape),
        "mode2_degree3_complex": list(mode2.shape),
        "combined": list(design.shape),
        "layout": (
            "mode0 degree-3 real [64], mode1 degree-2 complex [72], "
            "mode2 degree-3 complex [128]"
        ),
    }


def run(output_path=OUTPUT_PATH, npz_path=NPZ_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    npz_path = Path(npz_path)
    frozen_raw = FROZEN_PATH.read_bytes()
    candidate_raw = CANDIDATE_PATH.read_bytes()
    mean_raw = MEAN_PATH.read_bytes()
    frozen = json.loads(frozen_raw)
    candidate = json.loads(candidate_raw)
    if frozen.get("status") not in ("completed", None):
        raise ValueError("Frozen geometry report is not complete")
    if candidate.get("status") != "completed":
        raise ValueError("Current endpoint candidate is not complete")
    selected = candidate["selected"]
    control = np.asarray(selected["tangent_coefficients"], dtype=float)
    if control.shape != (264,):
        raise ValueError(f"Expected current 264-control candidate, got {control.shape}")
    wave_coefficients = _unpack_packed(selected["coefficients_original"])
    geometry = frozen["inputs"]["wave"]
    mean_inputs = frozen["inputs"]["mean"]
    timesteps = frozen["timesteps"]
    tau = float(mean_inputs["tau"])
    hspace = float(timesteps["hspace"])
    htime = float(timesteps["htime"])
    nu = float(mean_inputs["nu"])
    breaks = list(frozen["inputs"]["geometry"]["radial_breaks"])
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "wave_integrated": False,
        "scope": (
            "Fresh order-20, 12-angle higher-quadrature cache for the current "
            "264-control enriched candidate. The base residual is a single "
            "finite-difference full-field replay with no tangent controls; the "
            "264-column tangent response is analytic basis_data. No optimization, "
            "finite-time, PDE, or recursion acceptance is claimed."
        ),
        "sources": {
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "current_candidate": {"path": CANDIDATE_PATH.name, "sha256": _sha256(CANDIDATE_PATH)},
            "mean": {"path": MEAN_PATH.name, "sha256": _sha256(MEAN_PATH)},
        },
        "inputs": {
            "order": ORDER,
            "angles": ANGLES,
            "angle_shift": SHIFT,
            "tau": tau,
            "physical_time": "t=-tau",
            "nu": nu,
            "hspace": hspace,
            "htime": htime,
            "radial_breaks": breaks,
            "geometry_source": frozen["inputs"]["wave"],
            "mean_source": mean_inputs,
            "control_count": 264,
            "degrees": {"mode0": MODE0_DEGREE, "mode1": MODE1_DEGREE, "mode2": MODE2_DEGREE},
            "coefficients_original": _packed_complex(wave_coefficients).tolist(),
            "coefficient_layout": "packed [[Re(c_j), Im(c_j)]], j=0..26; current endpoint candidate",
            "tangent_control_layout": "mode0 degree-3 real [64], mode1 degree-2 complex [72], mode2 degree-3 complex [128]",
            "tangent_coefficients": control.tolist(),
        },
        "old_training_reference": {
            "source": CANDIDATE_PATH.name,
            "point_count": selected["training_momentum"]["point_count"],
            "candidate_training_momentum": selected["training_momentum"],
            "candidate_training_definition": "previous order-9 projection grid with fixed candidate tangent response",
        },
    }
    _save(report)
    print(json.dumps({"stage": "sources_loaded", "order": ORDER, "angles": ANGLES}), flush=True)

    mean, mean_report = load_saved_field()
    replacements = install_in_field(mean)
    points, weights, grid_metadata = _grid(
        mean,
        geometry["center"],
        geometry["widths"],
        tau,
        breaks,
        ORDER,
        ANGLES,
        SHIFT,
    )
    report["grid"] = {
        **grid_metadata,
        "weights_sum": float(np.sum(weights)),
        "point_count": int(len(points)),
        "grouped_backend_replacements": int(replacements),
    }
    report["status"] = "grid_built"
    _save(report)
    print(json.dumps({"stage": "grid_built", "point_count": len(points), "weights_sum": float(np.sum(weights))}), flush=True)

    # Current wave at the reference time, with zero tangent controls.  The
    # LocalPotentialField uses the frozen geometry but the current endpoint
    # wave coefficient snapshot from enriched_mean_endpoint_tangent.json.
    q2 = (int(geometry["degree"]) + 1) ** 2
    initial = (np.zeros(3 * q2, complex), wave_coefficients, np.zeros(3 * q2, complex))
    zero_derivatives = tuple(np.zeros(3 * q2, complex) for _ in (0, 1, 2))
    zero_pressures = tuple(np.zeros(q2, complex) for _ in (0, 1, 2))
    field = LocalPotentialField(
        mean,
        geometry["center"],
        geometry["widths"],
        int(geometry["degree"]),
        {0: np.zeros(2), 1: np.asarray(geometry["carrier"], dtype=float),
         2: 2.0 * np.asarray(geometry["carrier"], dtype=float)},
        initial,
        zero_derivatives,
        zero_pressures,
        tau,
    )
    print(json.dumps({"stage": "base_fd_started", "point_count": len(points)}), flush=True)
    base_residual = momentum(jets(field, points, tau, hspace, htime))
    report["base_replay"] = {
        "method": "affine_momentum.jets five-point spatial/temporal FD on LocalPotentialField",
        "tangent_controls": "zero",
        "metric": _metric(base_residual, weights),
    }
    report["status"] = "base_fd_complete"
    _save(report)
    print(json.dumps({"stage": "base_fd_complete", "l2": report["base_replay"]["metric"]["momentum_volume_L2"]}), flush=True)

    print(json.dumps({"stage": "tangent_design_started", "point_count": len(points)}), flush=True)
    tangent_design, layout = _mixed_tangent_design(points, geometry)
    candidate_residual = base_residual + (tangent_design @ control).reshape(-1, 3)
    report["tangent_design"] = {
        "shape": list(tangent_design.shape),
        "layout": layout,
        "method": "supported_fourier_basis.basis_data analytic velocity and pressure-gradient responses",
        "candidate_metric": _metric(candidate_residual, weights),
        "candidate_response_metric": _metric((tangent_design @ control).reshape(-1, 3), weights),
        "candidate_source_training_metric": selected["training_momentum"],
        "base_max_residual": _max_residual_detail(base_residual, points, weights),
        "candidate_max_residual": _max_residual_detail(candidate_residual, points, weights),
    }
    arrays = {
        "points": np.asarray(points, dtype=float),
        "weights": np.asarray(weights, dtype=float),
        "residual": np.asarray(base_residual, dtype=float),
        "tangent_design": np.asarray(tangent_design, dtype=float),
        "coefficients_original": _packed_complex(wave_coefficients),
        "tangent_coefficients": np.asarray(control, dtype=float),
        "center": np.asarray(geometry["center"], dtype=float),
        "widths": np.asarray(geometry["widths"], dtype=float),
        "carrier": np.asarray(geometry["carrier"], dtype=float),
        "tau": np.asarray([tau], dtype=float),
        "hspace": np.asarray([hspace], dtype=float),
        "htime": np.asarray([htime], dtype=float),
        "viscosity": np.asarray([nu], dtype=float),
        "orders": np.asarray([ORDER, ANGLES], dtype=int),
    }
    np.savez_compressed(npz_path, **arrays)
    report["npz"] = {
        "path": npz_path.name,
        "sha256": _sha256(npz_path),
        "size_bytes": int(npz_path.stat().st_size),
        "keys": {key: list(value.shape) for key, value in arrays.items()},
        "required_keys": ["points", "weights", "residual", "tangent_design", "coefficients_original"],
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "point_count": len(points),
        "tangent_shape": list(tangent_design.shape),
        "base_l2": report["base_replay"]["metric"]["momentum_volume_L2"],
        "candidate_l2": report["tangent_design"]["candidate_metric"]["momentum_volume_L2"],
        "npz": npz_path.name,
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--npz", type=Path, default=NPZ_PATH)
    args = parser.parse_args()
    run(args.output, args.npz)
