"""Bounded full-momentum fit for a compact meridional time-slope.

The saved broad-shear control is evaluated at one instantaneous reference
scale ``k0``.  The velocity directions in this file are multiplied by
``k-k0``; consequently their instantaneous velocity is zero while their
physical-time derivative is nonzero because ``t=-tau`` and
``dk/dt=1/(tau*log(2))``.  The fit includes the complete Cartesian

    R = partial_t u + (u dot grad) u + grad p - nu Delta u

residual, with pressure gradients evaluated by the same fourth-order
Cartesian finite-difference jet as the velocity terms.

This is an instantaneous, finite-grid, unconstrained exploratory diagnostic.
It does not establish the moment/cone constraints, a continuum estimate,
scale recursion, forcing, or Navier--Stokes acceptance.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_shear_dynamic_control import (  # noqa: E402
    BroadShearSlope,
    load_saved_field as load_dynamic_field,
)
from broad_shear_pressure import BroadCentrifugalPressure  # noqa: E402
from fourier_shear_feasibility import ControlledMean  # noqa: E402
from midplane_outer_slope_patch_repair import P_WINDOWS  # noqa: E402
from midplane_resolved_feasibility import ZeroBackground  # noqa: E402
from outer_pressure_modes import OuterPressure  # noqa: E402
from outer_swirl_slope import OuterSwirlSlope  # noqa: E402
from separated_moment_modes import (  # noqa: E402
    RADIAL_WINDOWS_THREE,
    SeparatedMomentModes,
)


DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
OUTPUT_PATH = ROOT / "broad_meridional_momentum.json"

# The pressure windows retain the existing compact pressure directions and add
# one window exactly over the broad-shear compact centrifugal cutoff.
# Extend the separated streamfunction directions into the broad onset and
# outer collar.  The three middle windows retain the existing exact basis
# geometry; the two extra compact windows expose the residual where the
# broad-shear cutoff actually starts and ends.
MERIDIONAL_WINDOWS = ((0.01, 0.12),) + tuple(RADIAL_WINDOWS_THREE) + ((0.72, 0.99),)
# Scalar pressure modes cover the same inner/onset pieces as well as the
# existing outer directions and the dedicated centrifugal cutoff window.
PRESSURE_WINDOWS = ((0.01, 0.12), (0.12, 0.38)) + tuple(P_WINDOWS[:2]) + ((0.72, 0.93),) + (P_WINDOWS[2],)
AXIAL_POWERS = (0, 1, 2)
REFERENCE_ORDER = 8
# Split at broad-shear, streamfunction and pressure cutoffs.  Five Gauss
# nodes per interval give 55 distinct radial training nodes; the holdout uses
# a different order on the same physical domain.
RADIAL_QUADRATURE_BREAKS = (
    0.01, 0.10, 0.12, 0.38, 0.40, 0.58, 0.60, 0.72,
    0.88, 0.92, 0.93, 0.99,
)


class MeridionalSlope:
    """One exact streamfunction mode with zero instantaneous velocity."""

    def __init__(self, unit, k0):
        self.unit = unit
        self.k0 = float(k0)
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(unit, name))

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        k = -np.log2(2.0 * float(np.asarray(tau).ravel()[0]))
        velocity, pressure = self.unit.fields(points, tau)
        return (k - self.k0) * velocity, (k - self.k0) * pressure


class CompactPressureDirection:
    """Pressure-only direction equal to the compact centrifugal primitive."""

    def __init__(self, primitive):
        self.primitive = primitive
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(primitive, name))

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        return np.zeros_like(points), self.primitive.increment(points, tau)


class CorrectedField:
    """Saved compact-pressure baseline plus serializable linear directions."""

    def __init__(self, baseline, modes, coefficients):
        self.baseline = baseline
        self.modes = tuple(modes)
        self.coefficients = np.asarray(coefficients, dtype=float)
        if self.coefficients.shape != (len(self.modes),):
            raise ValueError("Coefficient count does not match correction modes")
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(baseline, name))

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        velocity, pressure = self.baseline.fields(points, tau)
        velocity = np.asarray(velocity, dtype=float).copy()
        pressure = np.asarray(pressure, dtype=float).copy()
        for coefficient, mode in zip(self.coefficients, self.modes):
            if coefficient == 0.0:
                continue
            delta_velocity, delta_pressure = mode.fields(points, tau)
            velocity += coefficient * delta_velocity
            pressure += coefficient * delta_pressure
        return velocity, pressure


def _dynamic_reference(dynamic, dynamic_report):
    """Construct the no-broad reference used by the centrifugal primitive.

    The dynamic report has 19 controls: nine pressure controls, nine swirl
    state controls, and one broad-shear slope.  ``ControlledMean`` consumes
    only the first 18, so passing all 19 would create a ten-component state
    and fail during the pressure primitive's radial quadrature.
    """

    controls = np.asarray(dynamic_report["control"], dtype=float)
    if controls.shape != (19,):
        raise ValueError("Expected nineteen broad-shear dynamic controls")
    return ControlledMean(
        dynamic.base.base,
        np.asarray(dynamic_report["state"], dtype=float),
        controls[:18],
        float(dynamic_report["k"]),
    )


def build_compact_baseline(dynamic, dynamic_report, *, order=REFERENCE_ORDER):
    """Add the compact centrifugal primitive to the saved dynamic field."""

    increment = dynamic_report.get("mean_increment")
    if not increment or increment.get("type") != "broad_annular_shear":
        raise ValueError("Saved field has no broad-annular-shear increment")
    parameters = dict(increment["parameters"])
    amplitude = float(parameters["amplitude"])
    reference = _dynamic_reference(dynamic, dynamic_report)
    return BroadCentrifugalPressure(
        dynamic,
        reference,
        lambda tau: amplitude,
        parameters,
        dynamic_report["radial_breaks"],
        datum="compact",
        order=order,
        compact_window=(0.72, 0.93),
    )


def build_modes(dynamic, k0):
    """Build meridional, existing swirl, broad-shear, and pressure directions."""

    zero = ZeroBackground(dynamic.base)
    modes = []
    names = []

    # A single-knot unit mode is enough because the outer factor supplies the
    # desired k-slope.  The second knot keeps SeparatedMomentModes smooth when
    # the finite time jet samples a small neighborhood around k0.
    for radial_index, window in enumerate(MERIDIONAL_WINDOWS):
        for axial_index, power in enumerate(AXIAL_POWERS):
            amplitudes = np.zeros((2, 2, len(MERIDIONAL_WINDOWS), len(AXIAL_POWERS)))
            amplitudes[0, 1, radial_index, axial_index] = 1.0
            unit = SeparatedMomentModes(
                zero,
                amplitudes,
                windows=MERIDIONAL_WINDOWS,
                knots=(float(k0), float(k0) + 8.0),
                axial_powers=AXIAL_POWERS,
            )
            modes.append(MeridionalSlope(unit, k0))
            names.append(f"meridional_poloidal_window_{window[0]:.2f}_{window[1]:.2f}_eta^{power}")

    # These are the existing compact swirl time directions used by the
    # broad-shear control screen.  They are retained as candidate directions,
    # but no old moment or cone row is imposed in this residual fit.
    for index in range(9):
        modes.append(OuterSwirlSlope(zero, np.eye(9)[index], k0, P_WINDOWS))
        names.append(f"existing_outer_swirl_slope_{index}")

    params = dict(dynamic.base.base.parameters()) if hasattr(dynamic.base.base, "parameters") else {}
    # The saved broad-shear parameters are attached by the dynamic report and
    # are filled by the caller after construction when this direction is used.
    modes.append(BroadShearSlope(zero, k0, **params))
    names.append("existing_broad_shear_slope")

    # Scalar pressure modes have zero velocity and therefore enter the fit
    # through their complete Cartesian gradient, including radial cutoffs.
    for index in range(len(PRESSURE_WINDOWS) * 3):
        modes.append(
            OuterPressure(
                zero,
                np.eye(len(PRESSURE_WINDOWS) * 3)[index],
                windows=PRESSURE_WINDOWS,
            )
        )
        window = PRESSURE_WINDOWS[index // 3]
        power = index % 3
        names.append(f"pressure_window_{window[0]:.2f}_{window[1]:.2f}_eta^{power}")

    return modes, names


def build_modes_with_parameters(dynamic, dynamic_report, k0, compact_pressure=None):
    """Build directions while taking broad-shear shape parameters from JSON."""

    modes, names = build_modes(dynamic, k0)
    # ``build_modes`` creates the broad direction last among velocity modes.
    broad_index = len(MERIDIONAL_WINDOWS) * len(AXIAL_POWERS) + 9
    params = dict(dynamic_report["mean_increment"]["parameters"])
    params.pop("amplitude", None)
    zero = ZeroBackground(dynamic.base)
    modes[broad_index] = BroadShearSlope(zero, k0, **params)
    if compact_pressure is not None:
        modes.append(CompactPressureDirection(compact_pressure))
        names.append("compact_centrifugal_pressure_primitive")
    return modes, names


def quadrature_nodes(inner, tau, eta_order, radial_order, angle_order,
                     *, eta_shift=0.0, radial_shift=0.0,
                     eta_bounds=(-0.5, 0.5), radial_bounds=(0.01, 0.99),
                     angle_shift=0.0, radial_breaks=None):
    """Return Cartesian annulus nodes and physical-volume weights.

    The integration coordinates are ``eta`` and ``y`` with
    ``r = ri*(1+(ratio-1)*y)`` and ``z = sqrt(nu)*q**D*eta``, where
    ``q=tau/(1-eta**2)`` and ``ri=sqrt(2*nu*q*Xmax)``.  Thus the cylindrical
    volume Jacobian is ``2*pi*r*dr/dy*dz/deta``.  The angle quadrature sums
    to ``2*pi`` exactly for the evenly spaced nodes.
    """

    ge, we = leggauss(int(eta_order))
    elo, ehi = eta_bounds
    ylo, yhi = radial_bounds
    eta = (ehi - elo) * (ge + 1.0) / 2.0 + elo + float(eta_shift)
    if radial_breaks is None:
        gy, wy = leggauss(int(radial_order))
        y = (yhi - ylo) * (gy + 1.0) / 2.0 + ylo + float(radial_shift)
        radial_weights = (yhi - ylo) * wy / 2.0
    else:
        breaks = np.asarray(radial_breaks, dtype=float)
        if (breaks[0] != ylo or breaks[-1] != yhi
                or np.any(np.diff(breaks) <= 0.0)):
            raise ValueError("Radial split endpoints must match radial_bounds")
        gy, gw = leggauss(int(radial_order))
        y_parts = []
        w_parts = []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            y_parts.append((hi - lo) * (gy + 1.0) / 2.0 + lo)
            w_parts.append((hi - lo) * gw / 2.0)
        y = np.concatenate(y_parts) + float(radial_shift)
        radial_weights = np.concatenate(w_parts)
    if np.any(np.abs(eta) >= 1.0) or np.any(y <= 0.0) or np.any(y >= 1.0):
        raise ValueError("Shifted quadrature leaves the registered annulus")
    eta = np.asarray(eta, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = (float(angle_shift) + 2.0 * np.pi * np.arange(int(angle_order))
             / float(angle_order)) % (2.0 * np.pi)
    E, Y, T = np.meshgrid(eta, y, theta, indexing="ij")
    eta_weights = we * (ehi - elo) / 2.0
    WE, WY, WT = np.meshgrid(eta_weights, radial_weights,
                              np.full(len(theta), 2.0 * np.pi / len(theta)),
                              indexing="ij")
    Xmax = float(inner.p.X_max)
    ratio = float(getattr(inner, "outer_ratio", 16.0))
    q = float(tau) / (1.0 - E * E)
    ri = np.sqrt(2.0 * inner.nu * q * Xmax)
    radial_factor = 1.0 + (ratio - 1.0) * Y
    radius = ri * radial_factor
    X = Xmax * radial_factor * radial_factor
    points = inner.from_similarity(X.ravel(), E.ravel(), tau, T.ravel())
    D = 0.5 - float(inner.h)
    z_eta = np.sqrt(inner.nu) * q**D * (
        1.0 + (1.0 - 2.0 * inner.h) * E * E / (1.0 - E * E)
    )
    dr_dy = (ratio - 1.0) * ri
    weights = (radius * dr_dy * z_eta * WE * WY * WT).ravel()
    return points, weights


def evaluate_field(field, points, tau, weights):
    """Evaluate complete finite-difference momentum and divergence metrics."""

    hspace = 5.0e-4 * np.sqrt(field.nu * tau)
    htime = 1.0e-4 * tau
    jet = jets(field, points, tau, hspace, htime)
    residual = momentum(jet)
    divergence = np.trace(jet[1], axis1=1, axis2=2)
    norms = np.linalg.norm(residual, axis=1)
    return {
        "point_count": int(len(points)),
        "momentum_max": float(np.max(norms)),
        "momentum_volume_L2": float(np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)))),
        "momentum_volume_RMS": float(np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)) / np.sum(weights))),
        "divergence_max": float(np.max(np.abs(divergence))),
        "divergence_volume_L2": float(np.sqrt(np.sum(weights * divergence * divergence))),
        "residual_components": residual.tolist(),
    }


def fit_linear_directions(baseline, modes, points, weights, tau,
                          *, raw_reference=None):
    """Weighted complete-vector least squares at fixed instantaneous time."""

    hspace = 5.0e-4 * np.sqrt(baseline.nu * tau)
    htime = 1.0e-4 * tau
    baseline_jet = jets(baseline, points, tau, hspace, htime)
    baseline_residual = momentum(baseline_jet)
    raw_residual = None
    if raw_reference is not None:
        # The compact pressure direction is pressure-only.  At k0 the
        # baseline and raw field have identical velocity, spatial velocity
        # gradient and time derivative, so their complete residual difference
        # is exactly the primitive's Cartesian pressure gradient.  Reusing
        # this difference avoids a second expensive primitive quadrature.
        raw_residual = momentum(jets(raw_reference, points, tau, hspace, htime))
    columns = []
    for index, mode in enumerate(modes):
        if isinstance(mode, CompactPressureDirection) and raw_residual is not None:
            mode_residual = baseline_residual - raw_residual
        else:
            mode_residual = momentum(jets(mode, points, tau, hspace, htime))
        columns.append(mode_residual)
        print(json.dumps({"stage": "direction_jets", "index": index,
                          "point_count": len(points)}), flush=True)
    columns = np.stack(columns, axis=-1)  # point, component, direction
    sqrt_weights = np.sqrt(np.asarray(weights, dtype=float))
    weighted_columns = columns * sqrt_weights[:, None, None]
    weighted_target = -baseline_residual * sqrt_weights[:, None]
    flattened_columns = weighted_columns.reshape(-1, len(modes))
    flattened_target = weighted_target.reshape(-1)
    column_scales = np.linalg.norm(flattened_columns, axis=0)
    safe_scales = np.maximum(column_scales, 1.0e-30)
    scaled = flattened_columns / safe_scales[None, :]
    # Small ridge in dimensionless columns prevents an unidentifiable
    # direction from exploding while leaving the full residual objective
    # dominant.  No force or residual shortcut is introduced.
    ridge = 1.0e-8
    augmented = np.vstack((scaled, np.sqrt(ridge) * np.eye(len(modes))))
    target = np.concatenate((flattened_target, np.zeros(len(modes))))
    scaled_coefficients, _, rank, singular = np.linalg.lstsq(
        augmented, target, rcond=1.0e-12
    )
    coefficients = scaled_coefficients / safe_scales
    predicted = baseline_residual + np.einsum("ncp,p->nc", columns, coefficients)
    return coefficients, baseline_residual, columns, predicted, {
        "rank": int(rank),
        "singular_values": singular.tolist(),
        "column_scales": column_scales.tolist(),
        "ridge": ridge,
        "weighted_training_l2_before": float(np.sqrt(np.sum(weights[:, None] * baseline_residual**2))),
        "weighted_training_l2_after_linear": float(np.sqrt(np.sum(weights[:, None] * predicted**2))),
        "linear_model_max_before": float(np.max(np.linalg.norm(baseline_residual, axis=1))),
        "linear_model_max_after": float(np.max(np.linalg.norm(predicted, axis=1))),
    }


def load_saved_field(path=OUTPUT_PATH):
    """Reconstruct the compact meridional candidate from its JSON coefficients."""

    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if "coefficients" not in report:
        raise ValueError("Meridional momentum experiment has no fitted coefficients")
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    baseline = build_compact_baseline(dynamic, dynamic_report,
                                      order=int(report.get("compact_pressure_order", REFERENCE_ORDER)))
    modes, _ = build_modes_with_parameters(dynamic, dynamic_report,
                                           float(report["k"]),
                                           compact_pressure=baseline)
    return CorrectedField(baseline, modes, report["coefficients"]), report


def run():
    started = time.perf_counter()
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k0 = float(dynamic_report["k"])
    tau0 = 0.5 * 2.0**(-k0)
    baseline = build_compact_baseline(dynamic, dynamic_report)
    modes, mode_names = build_modes_with_parameters(dynamic, dynamic_report, k0,
                                                    compact_pressure=baseline)

    # Moderate tensor quadrature keeps the full-vector fit and replay bounded
    # while retaining the physical Jacobian and angular integration.
    train_points, train_weights = quadrature_nodes(
        dynamic.inner, tau0, 3, 5, 1,
        radial_bounds=(0.01, 0.99), radial_breaks=RADIAL_QUADRATURE_BREAKS,
    )
    holdout_points, holdout_weights = quadrature_nodes(
        dynamic.inner, tau0, 4, 4, 1,
        radial_bounds=(0.01, 0.99), radial_breaks=RADIAL_QUADRATURE_BREAKS,
        angle_shift=np.pi / 7.0,
    )

    print(json.dumps({"stage": "fit_started", "k": k0, "tau": tau0,
                      "train_points": len(train_points),
                      "holdout_points": len(holdout_points),
                      "mode_count": len(modes)}), flush=True)
    coefficients, baseline_residual, columns, predicted, fit_info = fit_linear_directions(
        baseline, modes, train_points, train_weights, tau0,
        raw_reference=dynamic,
    )
    candidate = CorrectedField(baseline, modes, coefficients)
    # Direct jet replay is retained separately from the cached linear model.
    raw_train = evaluate_field(dynamic, train_points, tau0, train_weights)
    compact_train = evaluate_field(baseline, train_points, tau0, train_weights)
    corrected_train = evaluate_field(candidate, train_points, tau0, train_weights)
    raw_holdout = evaluate_field(dynamic, holdout_points, tau0, holdout_weights)
    compact_holdout = evaluate_field(baseline, holdout_points, tau0, holdout_weights)
    corrected_holdout = evaluate_field(candidate, holdout_points, tau0, holdout_weights)

    report = {
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "status": "instantaneous_unconstrained_full_momentum_fit",
        "source_dynamic_report": DYNAMIC_PATH.name,
        "field_loader": "load_saved_field",
        "k": k0,
        "tau": tau0,
        "physical_time": "t=-tau",
        "dk_dt": float(1.0 / (tau0 * np.log(2.0))),
        "instantaneous_slope_velocity_zero": True,
        "compact_pressure_order": REFERENCE_ORDER,
        "compact_pressure": {
            "datum": "compact",
            "cutoff_window": [0.72, 0.93],
            "reference_control_slice": "control[:18] (pressure 9 + state 9); broad slope control retained separately",
            "source_wrapper": "broad_shear_pressure.BroadCentrifugalPressure",
            "radial_cutoff_gradient_included": True,
            "time_scope": "Only tau0 is evaluated; the callable amplitude is held fixed because finite-time evolution is outside scope.",
        },
        "residual_definition": "R=partial_t u+(u dot grad)u+grad p-nu Delta u; fourth-order Cartesian FD jets; t=-tau",
        "basis": {
            "mode_count": len(modes),
            "mode_names": mode_names,
            "meridional_windows": [list(v) for v in MERIDIONAL_WINDOWS],
            "axial_powers": list(AXIAL_POWERS),
            "pressure_windows": [list(v) for v in PRESSURE_WINDOWS],
            "existing_swirl_direction_count": 9,
            "broad_shear_direction_count": 1,
            "pressure_direction_count": len(PRESSURE_WINDOWS) * 3,
            "compact_pressure_direction_count": 1,
            "construction": "SeparatedMomentModes streamfunction radial/axial bases multiplied by (k-k0); pressure modes are compact scalar bumps",
        },
        "coefficients": coefficients.tolist(),
        "fit": fit_info,
        "training": {
            "eta_order": 3,
            "radial_order_per_split_interval": 5,
            "radial_split_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "angle_order": 1,
            "eta_bounds": [-0.5, 0.5],
            "radial_fraction_bounds": [0.01, 0.99],
            "point_count": len(train_points),
            "physical_volume": float(np.sum(train_weights)),
            "volume_weight": "r*(dr/dy)*(dz/deta) dy deta dtheta; mapped Gauss weights in y/eta and angular weights summing to 2*pi",
        },
        "holdout": {
            "independent": True,
            "eta_order": 4,
            "radial_order_per_split_interval": 4,
            "radial_split_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "eta_shift": 0.0,
            "radial_fraction_shift": 0.0,
            "angle_shift": float(np.pi / 7.0),
            "point_count": len(holdout_points),
            "physical_volume": float(np.sum(holdout_weights)),
            "volume_weight": "same domain and physical Jacobian, different Gauss orders in eta/y and a rotated angular node",
        },
        "baseline_raw_loader": {
            "training": raw_train,
            "holdout": raw_holdout,
        },
        "baseline_compact_pressure": {
            "training": compact_train,
            "holdout": compact_holdout,
        },
        "corrected": {
            "training": corrected_train,
            "holdout": corrected_holdout,
        },
        "compatibility": {
            "moment_constraints_used": False,
            "cone_constraints_used": False,
            "existing_dynamic_moment_and_cone_constraints_preserved": False,
            "status": "unconstrained_exploratory",
            "lost_compatibility": "The complete-vector volume fit did not impose the saved four moment equations or outer/wave cone rows; those constraints must be re-audited before any promotion.",
        },
        "scope": "One instantaneous reference tau with finite training and shifted holdout quadrature. Full vector residual, pressure gradients, compact radial cutoff and divergence are evaluated. No forcing=R shortcut, PDE acceptance, finite-time evolution, continuum bound, or scale recursion claim.",
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "stage": "fit_complete",
        "status": report["status"],
        "baseline_compact_max": compact_train["momentum_max"],
        "corrected_train_max": corrected_train["momentum_max"],
        "corrected_holdout_max": corrected_holdout["momentum_max"],
        "corrected_train_l2": corrected_train["momentum_volume_L2"],
        "corrected_holdout_l2": corrected_holdout["momentum_volume_L2"],
        "elapsed_seconds": report["elapsed_seconds"],
        "output": str(OUTPUT_PATH),
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
