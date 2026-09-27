"""Reusable cone rows for the instantaneous mean/wave tangent QP.

The exported inequality convention is

    cone_baseline + cone_control_rows @ tangent_x
        + cone_wave_forms(x_wave) >= cone_lower,

where the 81 rows are three H-transformed inequalities at each of the 27
saved locations.  ``tangent_x`` uses the full 180-column tangent layout; only
the first 36 mode-0 derivative/pressure columns are nonzero at the reference
time.  ``x_wave`` uses the all-real/all-imag 54-column complex mode-1 layout.

This is a matrix construction and independent snapshot cross-check.  It does
not optimize, integrate, or establish a continuum cone or PDE solution.
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

from broad_meridional_constrained import load_saved_field  # noqa: E402
from broad_shear_dynamic_control import load_saved_field as load_dynamic_field  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from joint_wave_mean_fit import (  # noqa: E402
    _grouped_cone_panels,
    _wave_metadata,
)
from supported_fourier_basis import basis_data  # noqa: E402
from fourier_patch_evolution import basis_jets  # noqa: E402


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
SAVED_REPLAY_PATH = ROOT / "broad_meridional_constrained_replay.json"
CURRENT_CANDIDATE_PATH = ROOT / "wave_dynamics_moment_codesign.json"
CURRENT_COMPAT_PATH = ROOT / "wave_dynamics_moment_compatibility.json"
PREVIOUS_CANDIDATE_PATH = ROOT / "wave_dynamics_codesign.json"
PREVIOUS_COMPAT_PATH = ROOT / "wave_dynamics_mean_compatibility.json"
OUTPUT_PATH = ROOT / "wave_mean_cone_projection.json"

CONE_ORDER = 64
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _candidate_wave(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if "selected_candidate" in data:
        wave = _wave_metadata(path)
    else:
        source_path = ROOT / data["source"]
        wave = _wave_metadata(source_path)
        wave["source"] = f"{Path(path).name}:selected + {source_path.name}:geometry"
        wave["potential_coefficients"] = data["selected"]["coefficients_original"]
    packed = np.asarray(wave["potential_coefficients"], dtype=float)
    if packed.shape != (27, 2):
        raise ValueError(f"expected 27 packed wave coefficients, got {packed.shape}")
    wave["coefficients"] = packed[:, 0] + 1j * packed[:, 1]
    return data, wave


def _cone_geometry(saved, inner, k, breaks):
    points, panels = _grouped_cone_panels(saved, inner, k, breaks, CONE_ORDER)
    diagnostics = saved["cone_problem"]["diagnostics"]
    H = np.asarray([row["H"] for row in diagnostics], dtype=float)
    return points, panels, H


def _mode0_rows(points, panels, H, center, widths):
    """Build target response and H-transformed 180-column tangent rows."""

    target_response = np.zeros((len(panels), 2, 36), dtype=float)
    for index, (sl, radii, weights, radius) in enumerate(panels):
        velocity, _, pressure_gradient = basis_data(
            points[sl], center, widths, 0, 2, (0.0, 0.0)
        )
        response = np.concatenate((velocity.real, pressure_gradient.real), axis=2)
        target_response[index, 0] = -np.einsum(
            "n,n,nq->q", weights, radii**2 / radius**2, response[:, 1, :]
        )
        target_response[index, 1] = -np.einsum(
            "n,n,nq->q", weights, radii / radius, response[:, 2, :]
        )
    transformed = np.einsum("rab,rbi->rai", H, target_response)
    rows = np.zeros((len(panels) * 3, 180), dtype=float)
    rows[:, :36] = transformed.reshape(len(panels) * 3, 36)
    return target_response, rows


def _wave_forms(points, panels, H, wave, nu, hspace):
    """Build 81 symmetric quadratic forms for the complete wave force."""

    forms = np.zeros((len(panels) * 3, 54, 54), dtype=float)
    for index, (sl, radii, weights, radius) in enumerate(panels):
        velocity, gradient, _, _, _ = basis_jets(
            points[sl], wave["center"], wave["widths"], wave["mode"],
            wave["degree"], wave["carrier"], nu, hspace,
        )
        velocity_real = np.concatenate((velocity, 1j * velocity), axis=-1)
        gradient_real = np.concatenate((gradient, 1j * gradient), axis=-1)
        bilinear = 0.5 * np.einsum(
            "ncdi,ndj->ncij", gradient_real, velocity_real.conj()
        ).real
        target_forms = np.stack((
            -np.tensordot(weights * radii**2 / radius**2,
                          bilinear[:, 1], axes=(0, 0)),
            -np.tensordot(weights * radii / radius,
                          bilinear[:, 2], axes=(0, 0)),
        ), axis=0)
        target_forms = 0.5 * (
            target_forms + np.swapaxes(target_forms, 1, 2)
        )
        forms[index * 3:(index + 1) * 3] = np.einsum(
            "rab,rbij->raij", H[index:index + 1], target_forms[None, ...]
        )[0]
    forms = 0.5 * (forms + np.swapaxes(forms, 1, 2))
    return forms


def _wave_vector(wave):
    return np.r_[wave["coefficients"].real, wave["coefficients"].imag]


def _quadratic(forms, x):
    return np.einsum("i,rij,j->r", x, forms, x)


def _targets(rows, key):
    return np.asarray([row[key] for row in rows], dtype=float)


def _snapshot_crosscheck(name, report, candidate, baseline_target,
                         target_response, wave_forms, H, lower):
    tangent = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    wave_x = _wave_vector(_candidate_wave_from_candidate(candidate)[1])
    predicted_mean_target = baseline_target + np.einsum(
        "rbi,i->rb", target_response, tangent[:36]
    )
    actual_mean_target = _targets(report["cone_replay"]["rows"], "mean_target")
    predicted_wave_target = None
    actual_wave_target = _targets(report["cone_replay"]["rows"], "wave_target")
    predicted_margin = (
        np.einsum("rab,rb->ra", H, predicted_mean_target).reshape(-1)
        + _quadratic(wave_forms, wave_x) - lower
    )
    actual_margin = np.asarray(
        [value for row in report["cone_replay"]["rows"]
         for value in row["total_margin"]], dtype=float
    )
    predicted_wave_target = actual_wave_target
    return {
        "tangent_l2": float(np.linalg.norm(tangent)),
        "wave_coefficient_l2": float(np.linalg.norm(wave_x)),
        "mean_target_reconstruction_max_abs": float(
            np.max(np.abs(predicted_mean_target - actual_mean_target))
        ),
        "wave_target_reconstruction_reference": predicted_wave_target.tolist(),
        "total_margin_reconstruction_max_abs": float(
            np.max(np.abs(predicted_margin - actual_margin))
        ),
        "actual_mean_target": actual_mean_target.tolist(),
        "predicted_mean_target": predicted_mean_target.tolist(),
        "actual_total_margin": actual_margin.tolist(),
        "predicted_total_margin": predicted_margin.tolist(),
    }


def _candidate_wave_from_candidate(candidate):
    """Resolve a candidate wrapper to its source geometry and selected c."""

    path = ROOT / candidate["source"]
    wave = _wave_metadata(path)
    wave["potential_coefficients"] = candidate["selected"]["coefficients_original"]
    packed = np.asarray(wave["potential_coefficients"], dtype=float)
    wave["coefficients"] = packed[:, 0] + 1j * packed[:, 1]
    return candidate, wave


def run():
    started = time.perf_counter()
    saved = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    current_candidate, current_wave = _candidate_wave(CURRENT_CANDIDATE_PATH)
    previous_candidate, previous_wave = _candidate_wave(PREVIOUS_CANDIDATE_PATH)
    current_compat = json.loads(CURRENT_COMPAT_PATH.read_text(encoding="utf-8"))
    previous_compat = json.loads(PREVIOUS_COMPAT_PATH.read_text(encoding="utf-8"))
    mean, mean_report = load_saved_field(MEAN_PATH)
    grouped_replacements = int(install_in_field(mean))
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k = float(saved["k"])
    tau = float(saved["tau"])
    nu = float(saved["nu"])
    hspace = 5.0e-4 * np.sqrt(nu * tau)
    breaks = [float(v) for v in saved["quadrature"]["radial_split_breaks"]]
    points, panels, H = _cone_geometry(saved, dynamic.inner, k, breaks)
    target_response, control_rows = _mode0_rows(
        points, panels, H, current_wave["center"], current_wave["widths"]
    )
    current_mean_target = _targets(current_compat["cone_replay"]["rows"], "mean_target")
    current_tangent = np.asarray(
        current_candidate["selected"]["tangent_coefficients"], dtype=float
    )
    baseline_target = current_mean_target - np.einsum(
        "rbi,i->rb", target_response, current_tangent[:36]
    )
    cone_baseline = np.einsum("rab,rb->ra", H, baseline_target).reshape(-1)
    cone_lower = np.tile(MARGIN, len(panels))
    wave_forms = _wave_forms(
        points, panels, H, current_wave, nu, hspace
    )
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "k": k,
        "tau": tau,
        "nu": nu,
        "physical_time": "t=-tau",
        "grouped_backend_replacements": grouped_replacements,
        "sources": {
            "mean": {"path": MEAN_PATH.name, "sha256": _sha256(MEAN_PATH)},
            "dynamic": {"path": DYNAMIC_PATH.name, "sha256": _sha256(DYNAMIC_PATH)},
            "current_candidate": {"path": CURRENT_CANDIDATE_PATH.name,
                                   "sha256": _sha256(CURRENT_CANDIDATE_PATH)},
            "current_compatibility": {"path": CURRENT_COMPAT_PATH.name,
                                       "sha256": _sha256(CURRENT_COMPAT_PATH)},
            "previous_candidate": {"path": PREVIOUS_CANDIDATE_PATH.name,
                                    "sha256": _sha256(PREVIOUS_CANDIDATE_PATH)},
            "previous_compatibility": {"path": PREVIOUS_COMPAT_PATH.name,
                                        "sha256": _sha256(PREVIOUS_COMPAT_PATH)},
            "saved_cone_replay": {"path": SAVED_REPLAY_PATH.name,
                                   "sha256": _sha256(SAVED_REPLAY_PATH)},
        },
        "geometry": {
            "location_count": len(panels),
            "inequality_count": int(len(cone_lower)),
            "cone_order": CONE_ORDER,
            "grouped_point_count": int(len(points)),
            "radial_split_breaks": breaks,
            "physical_volume_panels": "saved H locations with grouped radial quadrature",
        },
        "wave": {
            "mode": current_wave["mode"],
            "degree": current_wave["degree"],
            "center": np.asarray(current_wave["center"], dtype=float).tolist(),
            "widths": np.asarray(current_wave["widths"], dtype=float).tolist(),
            "carrier": np.asarray(current_wave["carrier"], dtype=float).tolist(),
            "current_coefficients_original": current_wave["potential_coefficients"],
            "previous_coefficients_original": previous_wave["potential_coefficients"],
            "coefficient_layout": "x_wave=[Re(c_0..c_26), Im(c_0..c_26)]",
        },
        "row_layout": {
            "cone_control_rows": "81x180; rows grouped by location then H row; only tangent columns 0:36 are nonzero",
            "tangent_layout": "full_wave_tangent 180 real controls: mode0 derivative 27, mode0 pressure 9, then interleaved complex mode1/2 controls",
            "wave_forms": "81 symmetric 54x54 forms; x.T @ cone_wave_forms[row] @ x",
            "inequality": "cone_baseline + cone_control_rows @ tangent_x + quadratic_wave >= cone_lower",
            "cone_lower_per_location": MARGIN.tolist(),
        },
        "cone_baseline": cone_baseline.tolist(),
        "baseline_target": baseline_target.tolist(),
        "cone_control_rows": control_rows.tolist(),
        "cone_wave_forms": wave_forms.tolist(),
        "cone_H": H.tolist(),
        "cone_lower": cone_lower.tolist(),
        "cross_checks": {},
        "scope": (
            "Reusable instantaneous mean/wave cone projection for the next "
            "tangent QP. Rows use saved 27 H transforms and cone64 grouped "
            "radial quadrature. No optimization, finite-time evolution, "
            "continuum cone, or PDE acceptance is claimed."
        ),
    }
    report["cross_checks"]["current_moment_candidate"] = _snapshot_crosscheck(
        "current", current_compat, current_candidate, baseline_target,
        target_response, wave_forms, H, cone_lower,
    )
    report["cross_checks"]["previous_unconstrained_candidate"] = _snapshot_crosscheck(
        "previous", previous_compat, previous_candidate, baseline_target,
        target_response, wave_forms, H, cone_lower,
    )
    # Compatibility snapshots store two-component wave targets, whereas this
    # projection stores their H-transformed 81-row values.  Project the direct
    # snapshots explicitly for a transparent reconstruction.
    for label, report_key, compat, candidate in (
        ("current", "current_moment_candidate", current_compat, current_candidate),
        ("previous", "previous_unconstrained_candidate", previous_compat, previous_candidate),
    ):
        wave_x = _wave_vector(_candidate_wave_from_candidate(candidate)[1])
        direct_target = _targets(compat["cone_replay"]["rows"], "wave_target")
        projected_direct = np.einsum("rab,rb->ra", H, direct_target).reshape(-1)
        report["cross_checks"][report_key]["wave_projected_direct_vs_quadratic_max_abs"] = float(
            np.max(np.abs(projected_direct - _quadratic(wave_forms, wave_x)))
        )
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "rows": list(control_rows.shape),
        "wave_forms": list(wave_forms.shape),
        "current_mean_reconstruction": report["cross_checks"]["current_moment_candidate"]["mean_target_reconstruction_max_abs"],
        "previous_mean_reconstruction": report["cross_checks"]["previous_unconstrained_candidate"]["mean_target_reconstruction_max_abs"],
        "current_wave_reconstruction": report["cross_checks"]["current_moment_candidate"]["wave_projected_direct_vs_quadratic_max_abs"],
        "previous_wave_reconstruction": report["cross_checks"]["previous_unconstrained_candidate"]["wave_projected_direct_vs_quadratic_max_abs"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
