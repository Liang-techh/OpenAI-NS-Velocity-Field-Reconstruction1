"""Joint linear control screen for a broad-shear k slope.

The fixed broad annular shear and the nine state values are taken from the
saved growth candidate.  At the candidate scale this script solves for nine
compact pressure coefficients, nine outer swirl k-slopes, and one broad
annular-shear k-slope.  Every slope direction has zero instantaneous
velocity, so the moment and cone rows are linear in these nineteen controls.

This is an instantaneous finite-grid diagnostic.  It does not certify a
continuum cone, an evolved wave, a pressure primitive, or Navier--Stokes
acceptance.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_annular_shear import BroadAnnularShear  # noqa: E402
from grouped_outer_cache import outer_cache  # noqa: E402
from midplane_integrated_moment_balance import evaluate  # noqa: E402
from midplane_outer_residual_source import outer_cones  # noqa: E402
from midplane_resolved_feasibility import ZeroBackground, _integrated_moment_coefficients  # noqa: E402
from midplane_outer_slope_patch_repair import P_WINDOWS  # noqa: E402
from outer_feedback_evolution import (  # noqa: E402
    SwirlValue,
    Trajectory,
    build_current,
    choose_control,
)
from outer_pressure_modes import OuterPressure  # noqa: E402
from outer_swirl_slope import OuterSwirlSlope  # noqa: E402
from radial_continuation import ROOT as PROJECT_ROOT  # noqa: E402
from joined_field import coordinates  # noqa: E402


GROWTH_PATH = ROOT / "broad_shear_growth.json"
FEASIBILITY_PATH = ROOT / "fourier_shear_feasibility.json"
MEAN_PATH = ROOT / "outer_feedback_evolution.json"
OUTPUT_PATH = ROOT / "broad_shear_dynamic_control.json"
ORDER = 48
REPLAY_MOMENT_ORDER = 96
REPLAY_CONE_ORDER = 64
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)


def old_outer_locations():
    locations = ([(eta, y) for eta in (-0.2, 0.0, 0.2)
                  for y in (0.5, 0.75)]
                 + [(eta + de, 0.75 + dy) for eta in (-0.2, 0.2)
                    for de in (-0.005, 0.0, 0.005)
                    for dy in (-0.005, 0.0, 0.005)])
    return list(dict.fromkeys(locations))


def wave_locations(wave_report):
    center = wave_report["center_similarity"]
    eta = float(center["eta"])
    y = float(center["y"])
    return [(eta, y), (eta - 0.03, y), (eta + 0.03, y),
            (eta, y - 0.07), (eta, y + 0.07)]


class BroadShearSlope:
    """Unit broad shear multiplied by k minus k0.

    Since physical time is t equal to minus tau and k equals minus log2 of
    2 tau, this gives dk/dt equal to 1 divided by tau log 2 in ``jets``.
    """

    def __init__(self, base, k0, *, power=2.0, onset=(0.01, 0.10),
                 offset=(0.88, 0.99), reference_y=0.325):
        self.k0 = float(k0)
        self.base = BroadAnnularShear(
            base, 1.0, power=power, onset=onset, offset=offset,
            reference_y=reference_y,
        )
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(self.base, name))

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        k = -np.log2(2.0 * float(np.asarray(tau).ravel()[0]))
        u, p = self.base.fields(points, tau)
        factor = k - self.k0
        return factor * u, factor * p


class DynamicControlledMean:
    """Callable field loader for a feasible nineteen-control result."""

    def __init__(self, base, amplitude, state, control, k0):
        self.base = base
        self.amplitude = float(amplitude)
        self.state = np.asarray(state, dtype=float)
        self.control = np.asarray(control, dtype=float)
        self.k0 = float(k0)
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(base, name))

    def fields(self, points, tau):
        tau = float(np.asarray(tau).ravel()[0])
        k = -np.log2(2.0 * tau)
        delta_k = k - self.k0
        broad = BroadAnnularShear(
            self.base,
            self.amplitude + delta_k * self.control[18],
            power=2.0,
            onset=(0.01, 0.10),
            offset=(0.88, 0.99),
            reference_y=0.325,
        )
        state = self.state + delta_k * self.control[9:18]
        field = OuterPressure(
            SwirlValue(broad, state), self.control[:9], windows=P_WINDOWS
        )
        return field.fields(points, tau)


def load_saved_field(path=OUTPUT_PATH):
    """Reload a feasible dynamic mean diagnostic as a callable field."""

    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if "control" not in report:
        raise ValueError("The dynamic control solve did not produce a field.")
    _, _, base = build_current()
    return DynamicControlledMean(
        base,
        report["fixed_broad_amplitude"],
        report["state"],
        report["control"],
        report["k"],
    ), report


def as_float_list(values):
    return np.asarray(values).tolist()


def cone_linear_rows(cache, labels, kind, center_offset=0):
    """Build the three normalized cone inequalities at each cached center."""

    baseline = cache["baseline"]
    unit_linear = np.asarray(cache["modes"][2])
    residual = momentum(baseline)
    rows = []
    rhs = []
    diagnostics = []
    for index, (panel, label) in enumerate(zip(cache["panels"], labels)):
        center_index = index + int(center_offset)
        sl, radii, weights, radius = panel
        u = baseline[0][center_index]
        gradient = baseline[1][center_index]
        F = float(u[1] / radius)
        shear = np.array([gradient[1, 0] - F, gradient[2, 0]], dtype=float)
        shear_norm = float(np.linalg.norm(shear))
        if shear_norm <= 1.0e-30:
            raise RuntimeError(f"Degenerate cone shear at {kind} node {index}")
        normal = shear / shear_norm
        tangent = np.array([-normal[1], normal[0]])
        lam = float(-2.0 * F * normal[0]
                    * (2.0 * F * normal[0] + shear_norm))
        if lam <= 0.0:
            raise RuntimeError(f"Nonpositive baseline lambda2 at {kind} node {index}: {lam}")
        target = np.array([
            -np.dot(weights * radii**2, residual[sl, 1]) / radius**2,
            -np.dot(weights * radii, residual[sl, 2]) / radius,
        ])
        unit_target = np.stack([
            -np.einsum("n,pn->p", weights * radii**2 / radius**2,
                        unit_linear[:, sl, 1]),
            -np.einsum("n,pn->p", weights * radii / radius,
                        unit_linear[:, sl, 2]),
        ])
        scale = max(abs(float(target @ normal)), 1.0)
        L = np.sqrt(lam)
        G = 0.95 * abs(2.0 * F * normal[0])
        if G <= 1.0e-30:
            raise RuntimeError(f"Degenerate cone G at {kind} node {index}")
        H = np.stack([
            -normal / scale,
            (-G * normal - L * tangent) / (G * scale),
            (-G * normal + L * tangent) / (G * scale),
        ])
        rows.extend(H @ unit_target)
        rhs.extend(H @ target - MARGIN)
        diagnostics.append({
            "index": int(index),
            "kind": kind,
            "label": label,
            "physical_center": np.asarray(cache["centers"][center_index], dtype=float).tolist(),
            "F": F,
            "lambda_squared": lam,
            "target_dot_N": float(target @ normal),
            "target_dot_K": float(target @ tangent),
            "source_T": target.tolist(),
            "H": H.tolist(),
            "normal": normal.tolist(),
            "baseline_margin_rows": (H @ target - MARGIN).tolist(),
            "unit_target": unit_target.tolist(),
        })
    return np.asarray(rows, dtype=float), np.asarray(rhs, dtype=float), diagnostics


def make_report(k, state, amplitude, descriptor, radial_breaks):
    return {
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "status": "initializing",
        "source": "broad_shear_growth.json",
        "mean_source": "outer_feedback_evolution.json",
        "k": float(k),
        "tau": float(0.5 * 2.0 ** (-k)),
        "dk_dt": float(1.0 / ((0.5 * 2.0 ** (-k)) * np.log(2.0))),
        "state": np.asarray(state, dtype=float).tolist(),
        "fixed_broad_amplitude": float(amplitude),
        "mean_increment": descriptor,
        "radial_breaks": list(radial_breaks),
        "unknowns": [*(f"pressure_{i}" for i in range(9)),
                     *(f"outer_swirl_k_slope_{i}" for i in range(9)),
                     "broad_shear_k_slope"],
        "scope": "Instantaneous linear nineteen-control broad-shear slope screen with finite moment and cone samples; no continuum, evolved wave, pressure primitive, or Navier--Stokes acceptance.",
    }


def run():
    started = time.perf_counter()
    growth = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    candidate = growth.get("selected_candidate")
    if candidate is None:
        raise ValueError("broad_shear_growth.json has no selected candidate")
    growth_descriptor = candidate.get("mean_increment")
    if not growth_descriptor or growth_descriptor.get("type") != "broad_annular_shear":
        raise ValueError("Expected broad_annular_shear selected candidate")
    params = dict(growth_descriptor["parameters"])
    amplitude = float(params.pop("amplitude"))
    _, base, current = build_current()
    mean_data = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    original = Trajectory(current, mean_data["nodes"])
    k = float(growth["initial_k"])
    delta_state = np.asarray(candidate.get("delta_state", np.zeros(9)), dtype=float)
    if delta_state.shape != (9,):
        raise ValueError("Growth candidate delta_state must have length 9")
    state = np.asarray(original.state(k), dtype=float) + delta_state
    descriptor = {
        "type": "broad_annular_shear",
        "parameters": {"amplitude": amplitude, **params},
    }
    feasibility = json.loads(FEASIBILITY_PATH.read_text(encoding="utf-8"))
    radial_breaks = list(feasibility.get("radial_breaks", ()))
    if not radial_breaks:
        raise ValueError("fourier_shear_feasibility.json has no radial_breaks")
    report = make_report(k, state, amplitude, descriptor, radial_breaks)
    output_path = OUTPUT_PATH

    def save():
        report["elapsed_seconds"] = float(time.perf_counter() - started)
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    save()
    print(json.dumps({"stage": "inputs_ready", "k": k,
                      "amplitude": amplitude, "state_dimension": len(state)}), flush=True)

    fixed_mean = BroadAnnularShear(current, amplitude, **params)
    field = SwirlValue(fixed_mean, state)
    zero = ZeroBackground(base)
    units = [OuterPressure(zero, np.eye(9)[index], windows=P_WINDOWS)
             for index in range(9)]
    units += [OuterSwirlSlope(zero, np.eye(9)[index], k, P_WINDOWS)
              for index in range(9)]
    units.append(BroadShearSlope(zero, k, **params))
    names = report["unknowns"]
    unit_check_tau = 0.5 * 2.0 ** (-k)
    unit_check_points = np.asarray([
        fixed_mean.inner.from_similarity(
            [fixed_mean.inner.p.X_max * (1.0 + 15.0 * y) ** 2], [eta], unit_check_tau
        )[0]
        for eta, y in old_outer_locations()[:6]
    ])
    unit_values = np.stack([unit.fields(unit_check_points, unit_check_tau)[0]
                            for unit in units])
    report["unit_zero_velocity_max"] = float(np.max(np.abs(unit_values)))
    report["broad_slope_definition"] = {
        "physical_time": "t=-tau",
        "dk_dt": report["dk_dt"],
        "factor": "k-k0",
        "instantaneous_velocity_max": float(np.max(np.abs(unit_values[-1]))),
    }
    save()
    print(json.dumps({"stage": "units_ready", "unknown_count": len(units),
                      "unit_zero_velocity_max": report["unit_zero_velocity_max"]}), flush=True)

    # The moment problem is assembled at the same order as the existing
    # outer-control solve.  Its quadratic array is retained as evidence and
    # checked to be numerically zero because all nineteen units have zero
    # instantaneous velocity.
    from adaptive_bridge_moment_fit import moment_slices  # noqa: E402

    moment_data = moment_slices(
        fixed_mean.inner, base, field, orders=(k,), n=ORDER,
        unit_fields=units, unit_fields_are_deltas=True,
        radial_breaks=radial_breaks,
    )[0]
    m, E, Q = _integrated_moment_coefficients(moment_data)
    report["moment_problem"] = {
        "order": ORDER,
        "equation_count": int(len(m)),
        "offset": m.tolist(),
        "linear": E.tolist(),
        "quadratic": Q.tolist(),
        "max_quadratic_abs": float(np.max(np.abs(Q))),
        "unit_moment_panels": len(moment_data["panels"]),
        "physical_panels": [
            {
                "eta": float(eta),
                "R": float(fixed_mean.inner.from_similarity(
                    [fixed_mean.inner.p.X_max * 16.0**2], [eta],
                    0.5 * 2.0 ** (-k)
                )[0, 0]),
                "z": float(fixed_mean.inner.from_similarity(
                    [fixed_mean.inner.p.X_max * 16.0**2], [eta],
                    0.5 * 2.0 ** (-k)
                )[0, 2]),
            }
            for eta in (-0.2, 0.2)
        ],
    }
    print(json.dumps({"stage": "moment_problem_ready",
                      "equation_count": len(m),
                      "max_quadratic_abs": float(np.max(np.abs(Q)))}), flush=True)

    old_locations = old_outer_locations()
    wave_report = json.loads((ROOT / "broad_shear_wave_cone.json").read_text(encoding="utf-8"))
    wave_locs = wave_locations(wave_report)
    locations = old_locations + wave_locs
    labels = ([{"eta": float(eta), "y": float(y), "set": "old_outer"}
               for eta, y in old_locations]
              + [{"eta": float(eta), "y": float(y), "set": "wave_center"}
                 for eta, y in wave_locs])
    cache = outer_cache(
        field, units, k=k, order=ORDER, locations=locations,
        radial_breaks=radial_breaks,
    )
    report["cache_metadata"] = {
        "point_count": int(cache["point_count"]),
        "quadrature_point_count": int(cache["quadrature_point_count"]),
        "group_count": int(cache["group_count"]),
        "location_count": len(locations),
        "old_outer_count": len(old_locations),
        "wave_center_count": len(wave_locs),
    }
    print(json.dumps({"stage": "grouped_cache_ready", **report["cache_metadata"]}), flush=True)

    old_cache = {
        **cache,
        "panels": cache["panels"][:len(old_locations)],
    }
    wave_start = len(old_locations)
    wave_cache = {
        **cache,
        "panels": cache["panels"][wave_start:],
    }
    A_old, b_old, old_diag = cone_linear_rows(old_cache, labels[:len(old_locations)], "old_outer")
    A_wave, b_wave, wave_diag = cone_linear_rows(
        wave_cache, labels[len(old_locations):], "wave_center",
        center_offset=wave_start,
    )
    A = np.vstack([A_old, A_wave])
    b = np.concatenate([b_old, b_wave])
    reference = np.r_[np.asarray(original.pressure(k), dtype=float),
                      np.asarray(original.state(k, 1), dtype=float), 0.0]
    equation_scale = np.maximum(np.max(np.abs(E), axis=1), 1.0e-12)
    En = E / equation_scale[:, None]
    mn = m / equation_scale
    report["cone_problem"] = {
        "old_outer_node_count": len(old_diag),
        "wave_center_node_count": len(wave_diag),
        "old_outer_constraint_count": int(len(b_old)),
        "wave_center_constraint_count": int(len(b_wave)),
        "diagnostics": old_diag + wave_diag,
    }
    report["linear_problem"] = {
        "unknown_names": names,
        "A": A.tolist(),
        "b": b.tolist(),
        "E": E.tolist(),
        "m": m.tolist(),
        "equation_scale": equation_scale.tolist(),
        "En": En.tolist(),
        "mn": mn.tolist(),
        "reference": reference.tolist(),
        "shape": {"inequalities": list(A.shape), "equalities": list(E.shape)},
    }
    report["status"] = "linear_problem_assembled"
    save()
    print(json.dumps({"stage": "linear_problem_saved",
                      "inequalities": list(A.shape), "equalities": list(E.shape),
                      "output": str(output_path)}), flush=True)

    lp = linprog(
        np.zeros(len(units)), A_ub=-A, b_ub=b,
        A_eq=En, b_eq=-mn,
        bounds=[(None, None)] * len(units), method="highs",
    )
    report["linear_program"] = {
        "success": bool(lp.success),
        "status": int(lp.status),
        "message": str(lp.message),
        "nit": int(getattr(lp, "nit", -1)),
        "equality_residual_max": None if lp.x is None else float(np.max(np.abs(En @ lp.x + mn))),
        "minimum_constraint": None if lp.x is None else float(np.min(A @ lp.x + b)),
    }
    save()
    print(json.dumps({"stage": "linear_program_done", **report["linear_program"]}), flush=True)
    if not lp.success:
        report["status"] = "linear_program_infeasible" if lp.status == 2 else "linear_program_failed"
        report["note"] = (
            "LP status 2 is a finite sampled infeasibility result, not a continuum impossibility proof."
            if lp.status == 2 else
            "The LP solver did not return an optimum; this is a numerical or bounded-solver failure, not a continuum impossibility proof."
        )
        save()
        return report

    try:
        control, selection = choose_control(En, mn, A, b, reference, lp.x)
    except RuntimeError as error:
        report["status"] = "control_projection_failed"
        report["projection_error"] = str(error)
        report["note"] = "The feasible LP point was retained; projection failure is numerical, not a continuum impossibility proof."
        save()
        return report
    control = np.asarray(control, dtype=float)
    report["control"] = control.tolist()
    report["control_selection"] = selection
    report["control_residuals"] = {
        "moment_scaled_max": float(np.max(np.abs(En @ control + mn))),
        "minimum_constraint": float(np.min(A @ control + b)),
        "moment_raw_max": float(np.max(np.abs(E @ control + m))),
    }
    report["status"] = "control_selected"
    save()
    print(json.dumps({"stage": "control_selected", **report["control_residuals"]}), flush=True)

    dynamic = DynamicControlledMean(current, amplitude, state, control, k)
    report["field_loader"] = "load_saved_field"
    moments = evaluate(dynamic, fixed_mean.inner, k, REPLAY_MOMENT_ORDER,
                       0.002, radial_breaks=radial_breaks)
    report["moment_replay"] = {
        "order": REPLAY_MOMENT_ORDER,
        "values": moments.tolist(),
        "max_abs": float(np.max(np.abs(moments))),
    }
    save()
    print(json.dumps({"stage": "moment_replay", **report["moment_replay"]}), flush=True)

    replay_outer = old_locations[:6]
    replay_wave = wave_locs
    outer_rows = outer_cones(
        dynamic, k, order=REPLAY_CONE_ORDER,
        radial_breaks=radial_breaks, locations=replay_outer,
    )
    wave_rows = outer_cones(
        dynamic, k, order=REPLAY_CONE_ORDER,
        radial_breaks=radial_breaks, locations=replay_wave,
    )
    report["cone_replay"] = {
        "order": REPLAY_CONE_ORDER,
        "outer6": outer_rows,
        "wave5": wave_rows,
        "outer6_pass_count": int(sum(row["cone_pass"] for row in outer_rows)),
        "wave5_pass_count": int(sum(row["cone_pass"] for row in wave_rows)),
    }
    save()
    print(json.dumps({"stage": "cone_replay",
                      "outer6_pass_count": report["cone_replay"]["outer6_pass_count"],
                      "wave5_pass_count": report["cone_replay"]["wave5_pass_count"]}), flush=True)

    tau = 0.5 * 2.0 ** (-k)
    radial = np.linspace(0.01, 0.99, 61)
    points = np.concatenate([
        fixed_mean.inner.from_similarity(
            fixed_mean.inner.p.X_max * (1.0 + 15.0 * radial) ** 2,
            np.full(len(radial), eta), tau
        ) for eta in (-0.3, -0.1, 0.0, 0.1, 0.3)
    ])
    residual = momentum(jets(
        dynamic, points, tau,
        0.0005 * np.sqrt(dynamic.nu * tau), 0.0001 * tau,
    ))
    report["sampled_momentum"] = {
        "point_count": int(len(points)),
        "peak": float(np.max(np.linalg.norm(residual, axis=1))),
        "sample_rms": float(np.sqrt(np.mean(np.sum(residual ** 2, axis=1)))),
    }
    report["status"] = "replayed"
    report["replay_pass_summary"] = {
        "moment_max_below_1e-3": bool(report["moment_replay"]["max_abs"] < 1.0e-3),
        "outer6_all_pass": bool(report["cone_replay"]["outer6_pass_count"] == 6),
        "wave5_all_pass": bool(report["cone_replay"]["wave5_pass_count"] == 5),
    }
    save()
    print(json.dumps({"stage": "replayed", "status": report["status"],
                      "replay_pass_summary": report["replay_pass_summary"],
                      "sampled_momentum": report["sampled_momentum"],
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    run()
