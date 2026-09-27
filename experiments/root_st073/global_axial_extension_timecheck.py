"""Bounded time and streamfunction checks for the global localization prototype.

This companion diagnostic leaves ``global_axial_extension.py`` frozen.  It
checks the analytic axisymmetric streamfunction reconstruction against the
actual mean field at the reference scale and at a small positive ``delta k``.
It also verifies the complete frozen wave support and plateau equality at the
endpoint ``k0 + 1e-6``.  These are construction checks only; no PDE,
finite-time, or scale-recursion acceptance is inferred.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_constrained import load_saved_field as load_mean  # noqa: E402
from enriched_shape_replay import Mode0TangentCorrection  # noqa: E402
from full_wave_tangent import LocalPotentialField, _unpack_full  # noqa: E402
from global_axial_extension import (  # noqa: E402
    CANDIDATE_PATH,
    OUTPUT_PATH as GLOBAL_OUTPUT_PATH,
    SNAPSHOT_PATH,
    GlobalAxialExtension,
    _streamfunction,
    build_candidate,
)
from grouped_joined_field import install_in_field  # noqa: E402
from joined_field import coordinates  # noqa: E402
from wave_higher_harmonic_tangent import (  # noqa: E402
    Mode2TangentCorrection,
    _decode_mode_block,
)


OUTPUT_PATH = ROOT / "global_axial_extension_timecheck.json"
GLOBAL_SOURCE_PATH = ROOT / "global_axial_extension.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.isfinite(points).all():
        raise ValueError("points must be finite")
    return points


def _walk_nodes(root):
    """Yield wrapper nodes once through the known field containers."""

    seen = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node is None or id(node) in seen:
            continue
        seen.add(id(node))
        yield node
        if isinstance(node, dict):
            stack.extend(node.values())
            continue
        if isinstance(node, (list, tuple, set)):
            stack.extend(node)
            continue
        for name in (
            "base", "current", "modified", "reference", "unit", "compact",
            "baseline", "values", "pressures", "modes", "field", "joined",
        ):
            if hasattr(node, name):
                value = getattr(node, name)
                if isinstance(value, (list, tuple, set, dict)):
                    stack.extend(value if not isinstance(value, dict) else value.values())
                else:
                    stack.append(value)


def _axisymmetric_points(mean, tau):
    """Use one inner, one bridge, and one outer-collar point at eta = .23."""

    eta = 0.23
    q = float(tau) / (1.0 - eta**2)
    ri = math.sqrt(2.0 * float(mean.nu) * q * float(mean.join_X))
    z = float(math.sqrt(mean.nu) * q ** (0.5 - mean.inner.h) * eta)
    radii = {
        "inner": 0.50 * ri,
        "bridge": ri * (1.0 + 0.50 * (float(mean.ratio) - 1.0)),
        "outer_collar": ri * (1.0 + 0.90 * (float(mean.ratio) - 1.0)),
    }
    labels = list(radii)
    points = np.asarray([[radii[label], 0.0, z] for label in labels], dtype=float)
    return labels, points


def _five_point_derivative(field, points, tau, axis, h):
    points = _points(points)
    offsets = (-2.0, -1.0, 1.0, 2.0)
    values = []
    for offset in offsets:
        shifted = points.copy()
        shifted[:, axis] += offset * h
        values.append(_streamfunction(field, shifted, tau))
    return (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h)


def _streamfunction_velocity_rows(mean, labels, points, tau, steps):
    actual, _ = mean.fields(points, tau)
    radius = np.hypot(points[:, 0], points[:, 1])
    rows = []
    for h in steps:
        dpsi_dr = _five_point_derivative(mean, points, tau, axis=0, h=float(h))
        dpsi_dz = _five_point_derivative(mean, points, tau, axis=2, h=float(h))
        predicted = np.column_stack((-dpsi_dz / radius, dpsi_dr / radius))
        observed = np.column_stack((actual[:, 0], actual[:, 2]))
        error = predicted - observed
        rows.append(
            {
                "space_step": float(h),
                "actual_meridional": observed.tolist(),
                "streamfunction_meridional": predicted.tolist(),
                "absolute_error": error.tolist(),
                "max_absolute_error": float(np.max(np.abs(error))),
                "relative_l2_error": float(
                    np.linalg.norm(error) / max(np.linalg.norm(observed), 1.0e-30)
                ),
            }
        )
    return rows


def _support_check(localized, geometry, tau):
    center = np.asarray(geometry["center"], dtype=float)
    widths = np.asarray(geometry["widths"], dtype=float)
    radial = np.asarray([center[0] - widths[0], center[0] + widths[0]], dtype=float)
    axial = np.asarray([center[1] - widths[1], center[1] + widths[1]], dtype=float)
    points = np.asarray(
        [[r, 0.0, z] for r in radial for z in axial], dtype=float
    )
    coord = coordinates(
        np.hypot(points[:, 0], points[:, 1]) / math.sqrt(localized.nu),
        points[:, 2] / math.sqrt(localized.nu),
        float(tau),
        localized.inner.h,
    )
    eta = np.asarray(coord["eta"], dtype=float)
    r1, r2 = localized.radial_radii(float(tau))
    eta_inside = bool(np.max(np.abs(eta)) <= localized.eta_flat + 1.0e-14)
    radial_inside = bool(np.max(np.abs(radial)) <= r1 + 1.0e-14)
    return {
        "tau": float(tau),
        "wave_support_points": points.tolist(),
        "eta_range": [float(np.min(eta)), float(np.max(eta))],
        "radial_range": [float(np.min(radial)), float(np.max(radial))],
        "eta_flat": float(localized.eta_flat),
        "eta_outer": float(localized.eta_outer),
        "radial_plateau_end": float(r1),
        "radial_cutoff_end": float(r2),
        "inside_eta_plateau": eta_inside,
        "inside_radial_plateau": radial_inside,
        "full_support_inside_plateau": bool(eta_inside and radial_inside),
    }


def _build_original(candidate, snapshot, tau0):
    """Build the same frozen candidate with its original, unlocalized mean."""

    original_mean, _ = load_mean()
    install_in_field(original_mean)
    geometry = snapshot["inputs"]["wave"]
    control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    wave = packed[:, 0] + 1j * packed[:, 1]
    legacy = np.r_[np.zeros(36), control[64:136], np.zeros(72)]
    derivatives, pressures = _unpack_full(legacy, 9)
    field = LocalPotentialField(
        original_mean,
        geometry["center"],
        geometry["widths"],
        2,
        {0: np.zeros(2), 1: np.asarray(geometry["carrier"]), 2: 2.0 * np.asarray(geometry["carrier"])},
        (np.zeros(27, complex), wave, np.zeros(27, complex)),
        derivatives,
        pressures,
        tau0,
    )
    field = Mode0TangentCorrection(field, geometry, control[:64], tau0)
    mode2_velocity, mode2_pressure = _decode_mode_block(control[-128:], 3)
    field = Mode2TangentCorrection(
        field,
        geometry["center"],
        geometry["widths"],
        np.asarray(geometry["carrier"]),
        3,
        mode2_velocity,
        mode2_pressure,
        tau0,
    )
    field.nu = original_mean.nu
    return field


def _broad_shear_audit(mean, points, tau, k0, k1):
    entries = []
    for node in _walk_nodes(mean):
        if type(node).__name__ != "BroadShearSlope":
            continue
        base = node.base
        base_psi = _streamfunction(base, points, tau)
        adapter_psi = _streamfunction(node, points, tau)
        expected_psi = (float(k1) - float(k0)) * base_psi
        base_velocity, _ = base.fields(points, tau)
        base_meridional = np.column_stack((base_velocity[:, 0], base_velocity[:, 2]))
        entries.append(
            {
                "base_type": type(base).__name__,
                "base_lower_type": type(getattr(base, "base", None)).__name__,
                "base_streamfunction_max_abs": float(np.max(np.abs(base_psi))),
                "adapter_streamfunction_max_abs": float(np.max(np.abs(adapter_psi))),
                "factor_corrected_expected_max_abs": float(np.max(np.abs(expected_psi))),
                "adapter_minus_factor_corrected_max_abs": float(
                    np.max(np.abs(adapter_psi - expected_psi))
                ),
                "base_meridional_velocity_max_abs": float(np.max(np.abs(base_meridional))),
                "current_base_is_meridionally_zero": bool(
                    np.max(np.abs(base_meridional)) <= 1.0e-14
                ),
                "factor_at_test_delta_k": float(k1 - k0),
                "source_fix_required_for_current_candidate": bool(
                    np.max(np.abs(base_psi)) > 1.0e-14
                    or np.max(np.abs(base_meridional)) > 1.0e-14
                ),
            }
        )
    return entries


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    full, localized, mean, candidate, snapshot, mean_report = build_candidate()
    k0 = float(snapshot["inputs"]["mean"]["k"])
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    delta_k = 1.0e-3
    k1 = k0 + delta_k
    tau1 = 0.5 * 2.0 ** (-k1)
    endpoint_delta_k = 1.0e-6
    k_endpoint = k0 + endpoint_delta_k
    tau_endpoint = 0.5 * 2.0 ** (-k_endpoint)
    steps_by_state = {}
    point_rows = {}
    all_labels = None
    all_points = []
    for label, k, tau in (("k0", k0, tau0), ("k0_plus_1e-3", k1, tau1)):
        labels, points = _axisymmetric_points(mean, tau)
        all_labels = labels
        all_points.append(points)
        h0 = 5.0e-4 * math.sqrt(float(mean.nu) * tau)
        steps = (h0, 0.5 * h0)
        steps_by_state[label] = [float(value) for value in steps]
        point_rows[label] = {
            "k": float(k),
            "tau": float(tau),
            "delta_k_from_k0": float(k - k0),
            "labels": labels,
            "points": points.tolist(),
            "rows": _streamfunction_velocity_rows(mean, labels, points, tau, steps),
            "psi_values": _streamfunction(mean, points, tau).tolist(),
            "velocity_max_abs": float(np.max(np.abs(mean.fields(points, tau)[0]))),
        }

    audit_points = np.vstack(all_points)
    broad_audit = _broad_shear_audit(mean, audit_points, tau1, k0, k1)
    geometry = snapshot["inputs"]["wave"]
    endpoint_support = _support_check(localized, geometry, tau_endpoint)
    endpoint_points = np.asarray(endpoint_support["wave_support_points"], dtype=float)
    endpoint_reference = np.vstack(
        (
            endpoint_points,
            np.asarray([geometry["center"][0], 0.0, geometry["center"][1]], dtype=float),
        )
    )
    # Keep the original wave constructor's reference tau at k0; evaluate both
    # localized and unlocalized fields at the endpoint tau_endpoint below.
    original = _build_original(candidate, snapshot, tau0)
    old_u, old_p = original.fields(endpoint_reference, tau_endpoint)
    new_u, new_p = full.fields(endpoint_reference, tau_endpoint)
    endpoint_plateau = {
        "k": float(k_endpoint),
        "tau": float(tau_endpoint),
        "delta_k_from_k0": float(endpoint_delta_k),
        "point_count": int(len(endpoint_reference)),
        "velocity_max_abs_difference": float(np.max(np.abs(old_u - new_u))),
        "pressure_max_abs_difference": float(np.max(np.abs(old_p - new_p))),
        "equality_pass": bool(
            np.max(np.abs(old_u - new_u)) <= 1.0e-10
            and np.max(np.abs(old_p - new_p)) <= 1.0e-10
        ),
        "support_containment_pass": bool(endpoint_support["full_support_inside_plateau"]),
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "finite_time_evolution_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Bounded reconstruction check for the frozen global axial localization "
            "at k0 and k0+1e-3. Five-point streamfunction derivatives are compared "
            "with actual mean meridional velocities; endpoint support and plateau "
            "equality are checked at k0+1e-6. No PDE or finite-time acceptance is claimed."
        ),
        "sources": {
            "candidate": {"path": CANDIDATE_PATH.name, "sha256": _sha256(CANDIDATE_PATH)},
            "frozen_cache": {"path": SNAPSHOT_PATH.name, "sha256": _sha256(SNAPSHOT_PATH)},
            "global_source": {
                "path": GLOBAL_SOURCE_PATH.name,
                "sha256": _sha256(GLOBAL_SOURCE_PATH),
            },
            "timecheck_source": {
                "path": Path(__file__).name,
                "sha256": _sha256(Path(__file__).resolve()),
            },
            "global_report": {"path": GLOBAL_OUTPUT_PATH.name, "sha256": _sha256(GLOBAL_OUTPUT_PATH)},
            "mean_source": {
                "path": mean_report.get("source_dynamic_report", "broad_shear_dynamic_control.json"),
                "coefficients": mean_report.get("coefficients"),
            },
        },
        "states": {
            "reference": {"k": float(k0), "tau": float(tau0)},
            "slope_test": {"k": float(k1), "tau": float(tau1), "delta_k": float(delta_k)},
            "physical_time": "t=-tau",
            "dk_dt_at_slope_test": float(1.0 / (tau1 * math.log(2.0))),
            "steps": steps_by_state,
        },
        "streamfunction_reconstruction": {
            "formula": "u_r=-partial_z psi/r, u_z=partial_r psi/r",
            "derivative": "five-point centered Cartesian radial/axial stencil",
            "points": point_rows,
            "point_types": all_labels,
        },
        "broad_shear_adapter_audit": {
            "nodes": broad_audit,
            "summary": (
                "The adapter now includes the (k-k0) factor used by BroadShearSlope.fields. "
                "The current candidate also has a BroadAnnularShear over "
                "ZeroBackground lower base, so the corrected factor is numerically "
                "zero in its meridional streamfunction."
            ),
            "source_modified": True,
        },
        "endpoint_support": endpoint_support,
        "endpoint_plateau_check": endpoint_plateau,
        "limitations": {
            "registered_slab": "The localization remains a tau-dependent similarity/radial cutoff on the registered finite slab.",
            "time": "The k0+1e-3 state is an independent instantaneous reconstruction check, not an integrated trajectory.",
            "pde": "No Navier--Stokes closure, force balance, or scale-recursion acceptance is asserted.",
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(output_path),
                "streamfunction_max_error": max(
                    row["max_absolute_error"]
                    for state in point_rows.values()
                    for row in state["rows"]
                ),
                "endpoint_support_inside_plateau": endpoint_support["full_support_inside_plateau"],
                "endpoint_plateau_velocity_error": endpoint_plateau["velocity_max_abs_difference"],
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
