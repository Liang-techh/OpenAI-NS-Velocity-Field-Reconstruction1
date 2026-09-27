"""Actual compatibility replay for the frozen degree-3 mean candidate.

Only the first 64 entries of the 264-control candidate are applied to the
axisymmetric mean: 48 degree-3 mode-0 velocity derivatives and 16 pressure
coefficients.  The saved nonaxisymmetric wave is added through its complete
angular-mean force.  Moment and cone matrices are loaded from the frozen
degree-3 row report and independently checked against the actual replay.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_constrained import load_saved_field as load_mean  # noqa: E402
from broad_shear_dynamic_control import load_saved_field as load_dynamic  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from joint_wave_mean_fit import _wave_metadata  # noqa: E402
from midplane_integrated_moment_balance import evaluate as replay_moments  # noqa: E402
from wave_dynamics_mean_compatibility import (  # noqa: E402
    MeanMode0Field,
    _cone_replay,
    _wave_moment_column,
)


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
REPLAY_PATH = ROOT / "enriched_mean_endpoint_replay.json"
ROWS_PATH = ROOT / "enriched_mean_constraint_rows.json"
OUTPUT_PATH = ROOT / "enriched_mean_compatibility.json"
MODE0_DEGREE = 3
MODE0_Q = (MODE0_DEGREE + 1) ** 2
MODE0_CONTROL_COUNT = 4 * MODE0_Q
MOMENT_ORDER = 96
MOMENT_ZFACTOR = 0.002
CONE_ORDER = 64


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict, path: Path) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _quadratic(forms, vector):
    return np.einsum("i,rij,j->r", vector, forms, vector)


def _resolve_wave(candidate: dict, candidate_path: Path):
    if all(key in candidate for key in ("center", "widths", "carrier", "mode", "degree")):
        wave = _wave_metadata(candidate_path)
    else:
        source_wave_path = ROOT / candidate["source"]
        wave = _wave_metadata(source_wave_path)
        wave["source"] = (
            f"{candidate_path.name}:selected + {source_wave_path.name}:geometry"
        )
    packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    if packed.shape != (27, 2):
        raise ValueError(f"expected wave coefficient shape (27, 2), got {packed.shape}")
    wave["potential_coefficients"] = packed.tolist()
    wave["coefficients"] = packed[:, 0] + 1j * packed[:, 1]
    return wave


def run(candidate_path=CANDIDATE_PATH, replay_path=REPLAY_PATH,
        output_path=OUTPUT_PATH):
    started = time.perf_counter()
    candidate_path = Path(candidate_path)
    replay_path = Path(replay_path)
    output_path = Path(output_path)
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    rows_report = json.loads(ROWS_PATH.read_text(encoding="utf-8"))
    saved = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    if candidate.get("status") != "completed":
        raise ValueError("candidate must be completed")
    if replay.get("status") != "completed":
        raise ValueError("independent replay must be completed")
    if rows_report.get("status") != "completed":
        raise ValueError("degree-3 row report must be completed")

    tangent = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    if tangent.shape != (264,):
        raise ValueError(f"expected 264 candidate controls, got {tangent.shape}")
    mode0_derivative = tangent[: 3 * MODE0_Q].astype(complex)
    mode0_pressure = tangent[3 * MODE0_Q : MODE0_CONTROL_COUNT].astype(complex)
    wave = _resolve_wave(candidate, candidate_path)

    mean, mean_report = load_mean(MEAN_PATH)
    grouped_replacements = int(install_in_field(mean))
    dynamic, dynamic_report = load_dynamic(DYNAMIC_PATH)
    install_in_field(dynamic)
    k = float(saved["k"])
    tau = float(saved["tau"])
    nu = float(saved["nu"])
    if not np.isclose(k, float(dynamic_report["k"]), rtol=1.0e-12, atol=0.0):
        raise ValueError("candidate and dynamic reference k differ")
    hspace = 5.0e-4 * np.sqrt(nu * tau)
    htime = 1.0e-4 * tau
    breaks = [float(value) for value in saved["quadrature"]["radial_split_breaks"]]
    field = MeanMode0Field(
        mean,
        wave["center"],
        wave["widths"],
        MODE0_DEGREE,
        mode0_derivative,
        mode0_pressure,
        tau,
    )

    report = {
        "status": "initialized",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "k": k,
        "tau": tau,
        "physical_time": "t=-tau",
        "nu": nu,
        "grouped_backend_replacements": grouped_replacements,
        "finite_difference_steps": {"hspace": hspace, "htime": htime},
        "sources": {
            "mean": {"path": MEAN_PATH.name, "sha256": _sha256(MEAN_PATH)},
            "dynamic": {"path": DYNAMIC_PATH.name, "sha256": _sha256(DYNAMIC_PATH)},
            "candidate": {"path": candidate_path.name, "sha256": _sha256(candidate_path)},
            "replay": {
                "path": replay_path.name,
                "sha256": _sha256(replay_path),
                "status": replay.get("status"),
            },
            "constraint_rows": {"path": ROWS_PATH.name, "sha256": _sha256(ROWS_PATH)},
        },
        "candidate_selected_metrics": candidate.get("selected", {}),
        "wave": {
            key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in wave.items()
            if key not in {"coefficients"}
        },
        "tangent_source": {
            "layout": "enriched_mean_endpoint_tangent: mode-0 degree-3 first 64 controls; remaining 200 controls omitted from angular mean",
            "selected_count": int(len(tangent)),
            "mode0_degree": MODE0_DEGREE,
            "mode0_control_count": MODE0_CONTROL_COUNT,
            "mode0_derivative": mode0_derivative.real.tolist(),
            "mode0_pressure": mode0_pressure.real.tolist(),
            "mode0_derivative_l2": float(np.linalg.norm(mode0_derivative)),
            "mode0_pressure_l2": float(np.linalg.norm(mode0_pressure)),
            "ignored_tail_l2": float(np.linalg.norm(tangent[MODE0_CONTROL_COUNT:])),
        },
        "quadrature": {
            "moment_order": MOMENT_ORDER,
            "moment_zfactor": MOMENT_ZFACTOR,
            "cone_order": CONE_ORDER,
            "radial_split_breaks": breaks,
            "cone_geometry": "saved 27 locations and H transforms; grouped radial panels",
        },
        "scope": (
            "Instantaneous actual order-96 moment and order-64 sampled cone replay "
            "for the degree-3 mode-0 mean controls plus the complete saved wave force. "
            "No optimization, oscillatory full-PDE validation, finite-time trajectory, "
            "continuum cone, or scale-recursion claim."
        ),
    }
    _save(report, output_path)
    print(
        json.dumps(
            {
                "stage": "initialized",
                "grouped_backend_replacements": grouped_replacements,
                "mode0_derivative_l2": report["tangent_source"]["mode0_derivative_l2"],
            }
        ),
        flush=True,
    )

    moment_rows = np.asarray(rows_report["moment_rows"], dtype=float)
    wave_moment_forms = np.asarray(rows_report["wave_moment_forms"], dtype=float)
    if moment_rows.shape != (4, MODE0_CONTROL_COUNT):
        raise ValueError(f"unexpected moment row shape {moment_rows.shape}")
    if wave_moment_forms.shape != (4, 54, 54):
        raise ValueError(f"unexpected wave moment form shape {wave_moment_forms.shape}")
    report["reusable_moment_linearization"] = {
        "moment_rows": moment_rows.tolist(),
        "baseline_moments": rows_report["baseline_moments"],
        "wave_moment_forms": wave_moment_forms.tolist(),
        "row_order": ["eta=-0.2 theta", "eta=-0.2 axial", "eta=+0.2 theta", "eta=+0.2 axial"],
        "source": ROWS_PATH.name,
    }
    report["status"] = "forms_loaded"
    report["elapsed_forms_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(
        json.dumps(
            {
                "stage": "forms_loaded",
                "moment_rows": list(moment_rows.shape),
                "wave_forms": list(wave_moment_forms.shape),
                "elapsed_seconds": report["elapsed_forms_seconds"],
            }
        ),
        flush=True,
    )

    mean_moments = np.asarray(
        replay_moments(
            field,
            dynamic.inner,
            k,
            MOMENT_ORDER,
            MOMENT_ZFACTOR,
            radial_breaks=breaks,
        ),
        dtype=float,
    )
    wave_moment, _ = _wave_moment_column(
        dynamic.inner, k, breaks, wave, nu, hspace, MOMENT_ORDER
    )
    joint_moments = mean_moments + wave_moment
    baseline_moments = np.asarray(rows_report["baseline_moments"], dtype=float)
    linearized_mean = baseline_moments + moment_rows @ tangent[:MODE0_CONTROL_COUNT]
    wave_x = np.r_[wave["coefficients"].real, wave["coefficients"].imag]
    quadratic_wave_moment = _quadratic(wave_moment_forms, wave_x)
    report["moment_replay"] = {
        "order": MOMENT_ORDER,
        "mean_mode0_corrected": mean_moments.tolist(),
        "wave_force_column": wave_moment.tolist(),
        "joint_angular_mean": joint_moments.tolist(),
        "mean_mode0_max_abs": float(np.max(np.abs(mean_moments))),
        "joint_max_abs": float(np.max(np.abs(joint_moments))),
        "baseline_moments": baseline_moments.tolist(),
        "linearized_mean_moments": linearized_mean.tolist(),
        "linearized_mean_reconstruction_max_abs": float(
            np.max(np.abs(linearized_mean - mean_moments))
        ),
        "quadratic_wave_moments": quadratic_wave_moment.tolist(),
        "direct_wave_moments": wave_moment.tolist(),
        "wave_reconstruction_max_abs": float(
            np.max(np.abs(quadratic_wave_moment - wave_moment))
        ),
        "saved_mean_reference": (
            json.loads((ROOT / "broad_meridional_constrained_replay.json").read_text())
            .get("moment_replay")
            if (ROOT / "broad_meridional_constrained_replay.json").exists()
            else None
        ),
        "source": "midplane_integrated_moment_balance.evaluate plus wave_mean_flux.single_mode_force radial panels",
    }
    report["status"] = "moment_complete"
    report["elapsed_moment_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(
        json.dumps(
            {
                "stage": "moment_complete",
                "joint_max_abs": report["moment_replay"]["joint_max_abs"],
                "linearized_mean_reconstruction_max_abs": report["moment_replay"]["linearized_mean_reconstruction_max_abs"],
                "wave_reconstruction_max_abs": report["moment_replay"]["wave_reconstruction_max_abs"],
                "elapsed_seconds": report["elapsed_moment_seconds"],
            }
        ),
        flush=True,
    )

    cone = _cone_replay(
        saved,
        dynamic.inner,
        k,
        breaks,
        field,
        wave,
        nu,
        hspace,
        htime,
    )
    cone_rows = np.asarray(rows_report["cone_control_rows"], dtype=float)
    cone_baseline = np.asarray(rows_report["cone_baseline"], dtype=float)
    cone_lower = np.asarray(rows_report["cone_lower"], dtype=float)
    cone_wave_forms = np.asarray(rows_report["cone_wave_forms"], dtype=float)
    predicted_margin = (
        cone_baseline
        + cone_rows @ tangent[:MODE0_CONTROL_COUNT]
        + _quadratic(cone_wave_forms, wave_x)
        - cone_lower
    )
    actual_margin = np.asarray(
        [value for row in cone["rows"] for value in row["total_margin"]],
        dtype=float,
    )
    report["cone_replay"] = cone
    report["cone_reconstruction"] = {
        "predicted_margin_min": float(np.min(predicted_margin)),
        "actual_margin_min": float(np.min(actual_margin)),
        "predicted_vs_actual_max_abs": float(
            np.max(np.abs(predicted_margin - actual_margin))
        ),
        "predicted_pass_count": int(np.sum(predicted_margin >= -1.0e-7)),
        "actual_pass_count": cone["inequality_pass_count"],
        "row_count": cone["inequality_count"],
        "source": ROWS_PATH.name,
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(
        json.dumps(
            {
                "stage": "completed",
                "joint_moment_max": report["moment_replay"]["joint_max_abs"],
                "cone_min_margin": cone["minimum_margin"],
                "cone_pass_count": cone["pass_count"],
                "cone_row_count": cone["row_count"],
                "cone_inequality_pass_count": cone["inequality_pass_count"],
                "cone_point_count": cone["point_count"],
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", default=str(CANDIDATE_PATH))
    parser.add_argument("--replay", default=str(REPLAY_PATH))
    parser.add_argument("--output", default=str(OUTPUT_PATH))
    args = parser.parse_args()
    run(args.candidate, args.replay, args.output)
