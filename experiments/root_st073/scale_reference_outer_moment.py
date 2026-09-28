"""Compact outer swirl neutralizer for the measured global axial moment.

The correction is

    delta u = C B((r-r_c)/w_r) B(z/w_z) (-y, x, 0),

with the C-infinity bump ``B(x)=exp(-1/(1-x^2))`` on ``|x|<1``.  It is an
axisymmetric exact curl, hence divergence-free, and its amplitude is chosen
to cancel the measured global ``M_z`` of the current reference.  The annulus
is disjoint from the inner wave/local correction supports.  It lies inside the
registered outer collar envelope, so this is a moment-neutrality diagnostic,
not a replacement for a globally constrained assembly.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from global_collar_tangent import BOXES  # noqa: E402


OUTPUT_PATH = ROOT / "scale_reference_outer_moment.json"
MOMENT_REPORT_PATH = ROOT / "scale_reference_angular_momentum.json"
REFERENCE_REPORT_PATH = ROOT / "scale_reference_trust_nonlinear_fit.json"
DEGREE = 2
RADIUS_CENTER = 0.0075
RADIUS_HALF_WIDTH = 0.001
AXIAL_CENTER = 0.0
AXIAL_HALF_WIDTH = 0.0008
NORMALIZATION_ORDERS = (4, 8, 16, 32, 64)
CALIBRATION_ORDER = 32
VALIDATION_ORDER = 64


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _bump(value: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return B, dB/dx, and d2B/dx2, zero outside the open unit interval."""

    value = np.asarray(value, dtype=float)
    bump = np.zeros_like(value)
    first = np.zeros_like(value)
    second = np.zeros_like(value)
    active = np.abs(value) < 1.0
    if np.any(active):
        x = value[active]
        one_minus = 1.0 - x * x
        exponent = np.exp(-1.0 / one_minus)
        q1 = -2.0 * x / one_minus**2
        q2 = -2.0 / one_minus**2 - 8.0 * x * x / one_minus**3
        bump[active] = exponent
        first[active] = exponent * q1
        second[active] = exponent * (q2 + q1 * q1)
    return bump, first, second


class OuterMomentNeutralizer:
    """Analytic compact axisymmetric swirl correction."""

    def __init__(self, amplitude: float, radius_center: float = RADIUS_CENTER,
                 radius_half_width: float = RADIUS_HALF_WIDTH,
                 axial_center: float = AXIAL_CENTER,
                 axial_half_width: float = AXIAL_HALF_WIDTH):
        self.amplitude = float(amplitude)
        self.radius_center = float(radius_center)
        self.radius_half_width = float(radius_half_width)
        self.axial_center = float(axial_center)
        self.axial_half_width = float(axial_half_width)
        if self.radius_center <= self.radius_half_width:
            raise ValueError("The outer annulus must stay away from the axis")
        if self.radius_half_width <= 0.0 or self.axial_half_width <= 0.0:
            raise ValueError("Bump half-widths must be positive")

    def _factors(self, points: np.ndarray) -> tuple[np.ndarray, ...]:
        x, y, z = np.asarray(points, dtype=float).T
        radius = np.hypot(x, y)
        radial = (radius - self.radius_center) / self.radius_half_width
        axial = (z - self.axial_center) / self.axial_half_width
        br, dbr, ddbr = _bump(radial)
        bz, dbz, ddbz = _bump(axial)
        f = self.amplitude * br * bz
        fr = self.amplitude * dbr * bz / self.radius_half_width
        frr = self.amplitude * ddbr * bz / self.radius_half_width**2
        fz = self.amplitude * br * dbz / self.axial_half_width
        fzz = self.amplitude * br * ddbz / self.axial_half_width**2
        return radius, f, fr, frr, fz, fzz

    def velocity(self, points: Any) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (N,3)")
        x, y, _ = points.T
        _, f, _, _, _, _ = self._factors(points)
        return np.column_stack((-f * y, f * x, np.zeros(len(points))))

    def jacobian(self, points: Any) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (N,3)")
        x, y, _ = points.T
        radius, f, fr, _, fz, _ = self._factors(points)
        inv_radius = 1.0 / np.maximum(radius, 1.0e-300)
        result = np.zeros((len(points), 3, 3), dtype=float)
        result[:, 0, 0] = -fr * x * y * inv_radius
        result[:, 0, 1] = -fr * y * y * inv_radius - f
        result[:, 0, 2] = -fz * y
        result[:, 1, 0] = fr * x * x * inv_radius + f
        result[:, 1, 1] = fr * x * y * inv_radius
        result[:, 1, 2] = fz * x
        return result

    def laplacian(self, points: Any) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (N,3)")
        x, y, _ = points.T
        radius, _, fr, frr, _, fzz = self._factors(points)
        scalar_factor = frr + fzz + 3.0 * fr / np.maximum(radius, 1.0e-300)
        return np.column_stack((-y * scalar_factor, x * scalar_factor, np.zeros(len(points))))


def _normalization(order: int) -> dict[str, float]:
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    radial_x = nodes
    radial_r = RADIUS_CENTER + RADIUS_HALF_WIDTH * radial_x
    br, _, _ = _bump(radial_x)
    axial_y = nodes
    bz, _, _ = _bump(axial_y)
    radial_one = RADIUS_HALF_WIDTH * float(np.sum(weights * radial_r**3 * br))
    axial_one = AXIAL_HALF_WIDTH * float(np.sum(weights * bz))
    radial_two = RADIUS_HALF_WIDTH * float(np.sum(weights * radial_r**3 * br**2))
    axial_two = AXIAL_HALF_WIDTH * float(np.sum(weights * bz**2))
    denominator = 2.0 * np.pi * radial_one * axial_one
    kinetic_factor = np.pi * radial_two * axial_two
    return {
        "order": int(order),
        "radial_integral_B": radial_one,
        "axial_integral_B": axial_one,
        "radial_integral_B2": radial_two,
        "axial_integral_B2": axial_two,
        "moment_denominator_for_C": denominator,
        "kinetic_factor_for_C2": kinetic_factor,
    }


def _max_speed(amplitude: float) -> float:
    # The axial bump is maximal at z=0. A dense radial scan is enough to make
    # this reported maximum a deterministic bound for the diagnostic.
    radial_x = np.linspace(-1.0, 1.0, 200001)
    radial_r = RADIUS_CENTER + RADIUS_HALF_WIDTH * radial_x
    br, _, _ = _bump(radial_x)
    return float(abs(amplitude) * np.max(radial_r * br) * np.exp(-1.0))


def _support_check() -> dict[str, Any]:
    outer_r = [RADIUS_CENTER - RADIUS_HALF_WIDTH, RADIUS_CENTER + RADIUS_HALF_WIDTH]
    outer_z = [AXIAL_CENTER - AXIAL_HALF_WIDTH, AXIAL_CENTER + AXIAL_HALF_WIDTH]
    geometry = json.loads((ROOT / "full_wave_frozen_cache.json").read_text(encoding="utf-8"))["inputs"]["wave"]
    wave_r = [float(geometry["center"][0] - geometry["widths"][0]), float(geometry["center"][0] + geometry["widths"][0])]
    wave_z = [float(geometry["center"][1] - geometry["widths"][1]), float(geometry["center"][1] + geometry["widths"][1])]
    disjoint_wave = outer_r[0] > wave_r[1] or outer_r[1] < wave_r[0] or outer_z[1] < wave_z[0] or outer_z[0] > wave_z[1]
    collar_overlap = any(
        outer_r[0] < box[1] and outer_r[1] > box[0] and outer_z[0] < box[3] and outer_z[1] > box[2]
        for box in BOXES
    )
    return {
        "outer_radial_support": outer_r,
        "outer_axial_support": outer_z,
        "inner_wave_radial_support": wave_r,
        "inner_wave_axial_support": wave_z,
        "disjoint_from_inner_wave_support": bool(disjoint_wave),
        "overlaps_registered_outer_collar_envelope": bool(collar_overlap),
    }


def _analytic_checks(neutralizer: OuterMomentNeutralizer) -> dict[str, Any]:
    points = np.asarray([
        [0.0075, 0.0, 0.0],
        [0.0072, 0.0003, 0.0001],
        [0.0068, -0.0002, -0.0003],
    ], dtype=float)
    jacobian = neutralizer.jacobian(points)
    h_space = 1.0e-8
    finite_difference_jacobian = np.zeros_like(jacobian)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = h_space
        finite_difference_jacobian[:, :, axis] = (
            neutralizer.velocity(points + offset) - neutralizer.velocity(points - offset)
        ) / (2.0 * h_space)
    laplacian = neutralizer.laplacian(points)
    h_laplacian = 2.0e-7
    finite_difference_laplacian = np.zeros_like(laplacian)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = h_laplacian
        finite_difference_laplacian += (
            neutralizer.velocity(points + offset)
            - 2.0 * neutralizer.velocity(points)
            + neutralizer.velocity(points - offset)
        ) / h_laplacian**2
    jacobian_scale = max(float(np.max(np.abs(finite_difference_jacobian))), 1.0)
    laplacian_scale = max(float(np.max(np.abs(finite_difference_laplacian))), 1.0)
    return {
        "points": points.tolist(),
        "jacobian_step": h_space,
        "laplacian_step": h_laplacian,
        "jacobian_max_abs_error": float(np.max(np.abs(jacobian - finite_difference_jacobian))),
        "jacobian_relative_error": float(np.max(np.abs(jacobian - finite_difference_jacobian)) / jacobian_scale),
        "laplacian_max_abs_error": float(np.max(np.abs(laplacian - finite_difference_laplacian))),
        "laplacian_relative_error": float(np.max(np.abs(laplacian - finite_difference_laplacian)) / laplacian_scale),
        "analytic_divergence_max_abs": float(np.max(np.abs(np.trace(jacobian, axis1=1, axis2=2)))),
    }


def run(output_path: Path = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    moment_report = json.loads(MOMENT_REPORT_PATH.read_text(encoding="utf-8"))
    if moment_report.get("status") != "completed":
        raise ValueError("Angular-momentum source report is not completed")
    reference_mz = float(moment_report["angular_momentum"]["M_z_reference"])
    normalization_rows = [_normalization(order) for order in NORMALIZATION_ORDERS]
    selected = next(row for row in normalization_rows if row["order"] == CALIBRATION_ORDER)
    amplitude = -reference_mz / selected["moment_denominator_for_C"]
    neutralizer = OuterMomentNeutralizer(amplitude)
    cancellation_rows = []
    for row in normalization_rows:
        corrected = reference_mz + amplitude * row["moment_denominator_for_C"]
        cancellation_rows.append({
            "order": row["order"],
            "correction_M_z": float(amplitude * row["moment_denominator_for_C"]),
            "corrected_M_z": float(corrected),
            "relative_corrected_to_reference": float(abs(corrected) / max(abs(reference_mz), 1.0e-300)),
        })
    validation = next(row for row in normalization_rows if row["order"] == VALIDATION_ORDER)
    correction_self_energy = amplitude * amplitude * selected["kinetic_factor_for_C2"]
    source_quadrature = moment_report.get("quadrature", {})
    source_mz_uncertainty = source_quadrature.get("coarse_to_fine_relative_difference")
    support = _support_check()
    report: dict[str, Any] = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Compact outer axisymmetric swirl calibration for necessary global axial "
            "angular-momentum neutrality. The correction is divergence-free and smooth "
            "by construction; no Navier--Stokes or recursive-scale acceptance is claimed."
        ),
        "sources": {
            "angular_momentum_report": {"path": MOMENT_REPORT_PATH.name, "sha256": _sha256(MOMENT_REPORT_PATH)},
            "reference_report": {"path": REFERENCE_REPORT_PATH.name, "sha256": _sha256(REFERENCE_REPORT_PATH)},
            "frozen_geometry": {"path": "full_wave_frozen_cache.json", "sha256": _sha256(ROOT / "full_wave_frozen_cache.json")},
            "collar_geometry": {"path": "global_collar_tangent.py", "sha256": _sha256(ROOT / "global_collar_tangent.py")},
        },
        "parameters": {
            "radius_center": RADIUS_CENTER,
            "radius_half_width": RADIUS_HALF_WIDTH,
            "axial_center": AXIAL_CENTER,
            "axial_half_width": AXIAL_HALF_WIDTH,
            "bump": "exp(-1/(1-x^2)) for |x|<1, zero otherwise",
            "amplitude_C": amplitude,
            "reference_M_z": reference_mz,
            "calibration_normalization_order": CALIBRATION_ORDER,
            "independent_validation_normalization_order": VALIDATION_ORDER,
        },
        "support": support,
        "normalization": {
            "orders": list(NORMALIZATION_ORDERS),
            "rows": normalization_rows,
            "cancellation_rows": cancellation_rows,
            "calibration_order": CALIBRATION_ORDER,
            "independent_validation_order": VALIDATION_ORDER,
            "calibration_identity_note": (
                "Cancellation at the calibration order is an identity by construction; "
                "only the independent validation order tests quadrature transfer."
            ),
            "independent_validation_corrected_M_z": float(
                reference_mz + amplitude * validation["moment_denominator_for_C"]
            ),
            "independent_validation_relative_residual": float(
                abs(reference_mz + amplitude * validation["moment_denominator_for_C"])
                / max(abs(reference_mz), 1.0e-300)
            ),
            "source_M_z_quadrature_uncertainty_relative": source_mz_uncertainty,
            "source_uncertainty_note": (
                "The reference M_z inherits the coarse-to-fine discrepancy reported by "
                "scale_reference_angular_momentum.json; outer normalization quadrature "
                "does not reduce that source-field uncertainty."
            ),
            "moment_formula": "delta M_z = C * 2*pi * integral(r^3 B_r dr) * integral(B_z dz)",
            "correction_self_energy_formula": "delta E_self = pi*C^2*integral(r^3 B_r^2 dr)*integral(B_z^2 dz)",
        },
        "correction": {
            "reference_M_z": reference_mz,
            "calibration_order": CALIBRATION_ORDER,
            "correction_M_z_calibrated": float(amplitude * selected["moment_denominator_for_C"]),
            "corrected_M_z_calibrated": float(reference_mz + amplitude * selected["moment_denominator_for_C"]),
            "independent_validation_order": VALIDATION_ORDER,
            "corrected_M_z_validation": float(reference_mz + amplitude * validation["moment_denominator_for_C"]),
            "correction_self_energy": float(correction_self_energy),
            "total_energy_cross_term_uncomputed": True,
            "total_energy_note": "Reported energy is integral |delta u|^2/2 only; the cross term with the existing field is not computed.",
            "max_speed": _max_speed(amplitude),
            "velocity_formula": "C*B_r*B_z*(-y,x,0)",
            "jacobian_exposed": True,
            "laplacian_exposed": True,
            "divergence": "identically zero for the axisymmetric swirl form",
        },
        "neutrality_scope": {
            "moment_condition_only": True,
            "global_neutrality_certificate": False,
            "note": "The independent order-64 normalization validates the bump integral, while source M_z uncertainty remains inherited from its order-4/order-8 field quadrature.",
        },
        "analytic_checks": _analytic_checks(neutralizer),
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "reference_M_z": reference_mz,
        "amplitude_C": amplitude,
        "corrected_M_z_calibrated": report["correction"]["corrected_M_z_calibrated"],
        "corrected_M_z_validation": report["correction"]["corrected_M_z_validation"],
        "correction_self_energy": correction_self_energy,
        "max_speed": report["correction"]["max_speed"],
    }), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
