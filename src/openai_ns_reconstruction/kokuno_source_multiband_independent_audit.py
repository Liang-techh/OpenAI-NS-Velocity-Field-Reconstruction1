"""Independent Agent-4 re-audit of Agent-2 #461 source-scheduled multi-band curls.

This is the exact scientific contract of Agent-4 #451: fresh manufactured
source-compatible data, an independent Cartesian FD4 curl oracle, centered
FD6 coefficient derivatives, the same held-out seed, resolution ladder,
mutation tests, and frozen local guards.  The only intended dependency change
is Agent-2 #461's deterministic aggregation repair.  This remains a structural
preflight, not the formal full-domain Navier-Stokes gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import ceil
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_source_multiband_real_pair_family import KokunoSourceMultiBandRealPairFamily

SCHEMA = "kokuno-agent4-source-multiband-independent-audit-v2"
TASK_ID = "KOKUNO-A4-SOURCE-MULTIBAND-INDEPENDENT-REAUDIT-028"
BASE_PR = 461
BASE_HEAD = "aafb67af44c40794701ed0d66dbf0a08c1e45fb3"
PREVIOUS_AGENT4_PR = 451
PREVIOUS_AGENT4_HEAD = "3e8170bde17cfeee5189e453f304388e65564931"
SEED = 9173171
FD4_STEPS = (0.032, 0.016, 0.008)
COEFF_FD6_STEP = 2.0e-4
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5
FROZEN_GUARDS = {
    "source_schedule_relative_error_max": 5.0e-14,
    "finest_cartesian_curl_relative_rms_max": 2.0e-5,
    "curl_refinement_ratio_min": 4.0,
    "finest_normalized_divergence_rms_max": 2.0e-5,
    "finest_normalized_divergence_point_max": 8.0e-5,
    "divergence_refinement_ratio_min": 3.0,
    "shared_first_band_schedule_mutation_relative_rms_min": 0.10,
    "label_schedule_swap_mutation_relative_rms_min": 0.10,
}
CASES = (
    {"name": "low", "h": 0.0025, "ells": (6, 8, 10)},
    {"name": "mid", "h": 0.0060, "ells": (40, 42, 44)},
    {"name": "k_transition", "h": 0.0090, "ells": (221, 223, 225)},
)


def _canonical(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _schedule(ell: int, h: float) -> tuple[float, float, int]:
    q = float(2.0 ** (-int(ell)))
    eps = float(q**h)
    return q, eps, int(ceil(eps ** -0.5))


def _fd6(fn: Callable[[np.ndarray], np.ndarray], x: np.ndarray, h: float) -> np.ndarray:
    return (
        -fn(x - 3 * h)
        + 9 * fn(x - 2 * h)
        - 45 * fn(x - h)
        + 45 * fn(x + h)
        - 9 * fn(x + 2 * h)
        + fn(x + 3 * h)
    ) / (60 * h)


def _fd4_axis(
    fn: Callable[[np.ndarray], np.ndarray], p: np.ndarray, axis: int, h: float
) -> np.ndarray:
    m2 = p.copy()
    m1 = p.copy()
    p1 = p.copy()
    p2 = p.copy()
    m2[axis] -= 2 * h
    m1[axis] -= h
    p1[axis] += h
    p2[axis] += 2 * h
    return (fn(m2) - 8 * fn(m1) + 8 * fn(p1) - fn(p2)) / (12 * h)


def _partition(r: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = 0.72 + 0.11 * r - 0.07 * z + 0.025 * r * z
    b = 0.43 - 0.06 * r + 0.09 * z + 0.02 * r * r
    ar, az = 0.11 + 0.025 * z, -0.07 + 0.025 * r
    br, bz = -0.06 + 0.04 * r, np.full_like(r, 0.09)
    ca, sa, cb, sb = np.cos(a), np.sin(a), np.cos(b), np.sin(b)
    eta = np.stack((ca, sa * cb, sa * sb), axis=-1)
    dr = np.stack(
        (-sa * ar, ca * ar * cb - sa * sb * br, ca * ar * sb + sa * cb * br),
        axis=-1,
    )
    dz = np.stack(
        (-sa * az, ca * az * cb - sa * sb * bz, ca * az * sb + sa * cb * bz),
        axis=-1,
    )
    return eta, dr, dz


def _phase_n_raw(
    case: dict[str, Any], j: int, r: np.ndarray, th: np.ndarray, z: np.ndarray
):
    ell = case["ells"][j]
    k = _schedule(ell, case["h"])[2]
    p = (1, 2, 1)[j] / k
    a, c = 0.83 + 0.09 * j, 0.47 - 0.06 * j
    d, e, f = 0.035 * (j + 1), 0.012 * (j + 1), -0.010 * (j + 1)
    phase = p * th + a * r + c * z + d * r * z + e * r * r + f * z * z
    n = np.stack(
        (a + d * z + 2 * e * r, np.full_like(r, p) / r, c + d * r + 2 * f * z),
        axis=-1,
    )
    raw = np.stack(
        (
            (0.61 + 0.025 * r - 0.017 * z) + 1j * (0.08 + 0.013 * r + 0.007 * z),
            (-0.31 + 0.021 * z + 0.011 * r) + 1j * (0.12 - 0.009 * r + 0.006 * z),
            (0.24 + 0.018 * r * z - 0.008 * r) + 1j * (-0.07 + 0.010 * z - 0.005 * r),
        ),
        axis=-1,
    )
    return phase, n, raw


def _transverse_C(
    case: dict[str, Any],
    j: int,
    r: np.ndarray,
    th: np.ndarray,
    z: np.ndarray,
    k_override=None,
):
    phase, n, raw = _phase_n_raw(case, j, r, th, z)
    norm2 = np.sum(n * n, axis=-1)
    t = raw - n * (np.sum(n * raw, axis=-1) / norm2)[..., None]
    k = _schedule(case["ells"][j], case["h"])[2] if k_override is None else int(k_override)
    C = 1j * np.cross(n, t) / (k * norm2[..., None])
    return phase, n, t, C


def _inputs(case: dict[str, Any], xyz: np.ndarray) -> dict[str, Any]:
    pts = np.asarray(xyz, float)
    if pts.ndim == 1:
        pts = pts[None, :]
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r, th = np.hypot(x, y), np.arctan2(y, x)
    if np.any(r <= 0):
        raise ValueError("audit point crossed cylindrical axis")
    eta, dre, dze = _partition(r, z)
    phases, ns, ts, drs, dzs = [], [], [], [], []
    for j, _ell in enumerate(case["ells"]):
        ph, n, t, _C = _transverse_C(case, j, r, th, z)
        dr = _fd6(lambda rr: _transverse_C(case, j, rr, th, z)[3], r, COEFF_FD6_STEP)
        dz = _fd6(lambda zz: _transverse_C(case, j, r, th, zz)[3], z, COEFF_FD6_STEP)
        phases.append(ph)
        ns.append(n)
        ts.append(t)
        drs.append(dr)
        dzs.append(dz)
    return {
        "R": r,
        "theta": th,
        "phase": np.stack(phases, -1),
        "n_phi": np.stack(ns, -2),
        "t_plus": np.stack(ts, -2),
        "D_r_C_plus": np.stack(drs, -2),
        "D_z_C_plus": np.stack(dzs, -2),
        "eta": eta,
        "D_r_eta": dre,
        "D_z_eta": dze,
        "beta_labels": tuple((int(ell), (j, 0, 0)) for j, ell in enumerate(case["ells"])),
    }


def _cyl_to_cart(v: np.ndarray, th: float) -> np.ndarray:
    c, s = np.cos(th), np.sin(th)
    return np.asarray((v[0] * c - v[1] * s, v[0] * s + v[1] * c, v[2]), complex)


def _potential(case: dict[str, Any], j: int, p: np.ndarray, k_override=None) -> np.ndarray:
    x, y, z = np.asarray(p, float)
    r, th = float(np.hypot(x, y)), float(np.arctan2(y, x))
    if r <= 0:
        raise ValueError("Cartesian stencil crossed axis")
    ph, _n, _t, C = _transverse_C(
        case,
        j,
        np.asarray([r]),
        np.asarray([th]),
        np.asarray([z]),
        k_override,
    )
    eta = _partition(np.asarray([r]), np.asarray([z]))[0][0, j]
    k = _schedule(case["ells"][j], case["h"])[2] if k_override is None else int(k_override)
    return _cyl_to_cart(eta * C[0] * np.exp(1j * k * ph[0]), th)


def _curl_fd4(fn: Callable[[np.ndarray], np.ndarray], p: np.ndarray, h: float) -> np.ndarray:
    J = np.empty((3, 3), complex)
    for a in range(3):
        J[:, a] = _fd4_axis(fn, p, a, h)
    return np.asarray(
        (J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]), complex
    )


def _reference(
    case: dict[str, Any], p: np.ndarray, h: float, mode: str = "correct"
) -> np.ndarray:
    out = []
    A = 0.5 + case["h"]
    first_q, _e, first_k = _schedule(case["ells"][0], case["h"])
    sched = tuple(reversed(case["ells"])) if mode == "reversed" else case["ells"]
    for j, _ell in enumerate(case["ells"]):
        q, _eps, k = _schedule(sched[j], case["h"])
        if mode == "shared_first":
            q, k = first_q, first_k
        curl = _curl_fd4(lambda x, jj=j, kk=k: _potential(case, jj, x, kk), p, h)
        out.append(q ** (-A) * 2.0 * curl.real)
    return np.stack(out, 0)


def _public(case: dict[str, Any], p: np.ndarray) -> dict[str, Any]:
    return KokunoSourceMultiBandRealPairFamily(h=case["h"]).physical_family(**_inputs(case, p))


def _samples() -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    probes = []
    for case in CASES:
        r = np.concatenate((rng.uniform(0.42, 1.30, 8), rng.uniform(0.14, 0.24, 4)))
        th = rng.uniform(-2.7, 2.7, r.size)
        z = rng.uniform(-0.55, 0.55, r.size)
        for i in range(r.size):
            probes.append(
                {
                    "case": case,
                    "region": "axis_near" if i >= 8 else "off_grid",
                    "xyz": np.asarray(
                        (r[i] * np.cos(th[i]), r[i] * np.sin(th[i]), z[i]), float
                    ),
                }
            )
    return probes


def _rms(v) -> float:
    a = np.asarray(v, float)
    return float(np.sqrt(np.mean(a * a)))


def _curl_metrics(probes, h: float) -> dict[str, Any]:
    err, ref, point, bycase = [], [], [], {}
    for probe in probes:
        case, p = probe["case"], probe["xyz"]
        u = _public(case, p)["velocity_physical_cartesian_by_beta"][0]
        v = _reference(case, p, h)
        e, n = float(np.linalg.norm(u - v)), float(np.linalg.norm(v))
        err.append(e)
        ref.append(n)
        point.append(e / max(n, 1e-300))
        bycase.setdefault(case["name"], []).append((e, n))
    return {
        "step": h,
        "curl_relative_rms": _rms(err) / max(_rms(ref), 1e-300),
        "curl_point_relative_max": max(point),
        "parameter_cases": {
            k: _rms([x[0] for x in v]) / max(_rms([x[1] for x in v]), 1e-300)
            for k, v in bycase.items()
        },
    }


def _total(case, p):
    return np.asarray(_public(case, p)["velocity_physical_cartesian_total"][0], float)


def _div_metrics(probes, h: float) -> dict[str, Any]:
    div, grad, pts = [], [], []
    regions = {"off_grid": [], "axis_near": []}
    for probe in probes:
        case, p = probe["case"], probe["xyz"]
        J = np.empty((3, 3), float)
        for a in range(3):
            J[:, a] = _fd4_axis(lambda x: _total(case, x), p, a, h)
        d, g = float(abs(np.trace(J))), float(np.linalg.norm(J))
        n = d / max(g, 1e-300)
        div.append(d)
        grad.append(g)
        pts.append(n)
        regions[probe["region"]].append(n)
    return {
        "step": h,
        "normalized_divergence_rms": _rms(div) / max(_rms(grad), 1e-300),
        "normalized_divergence_point_max": max(pts),
        "regions": {k: {"normalized_point_max": max(v)} for k, v in regions.items()},
    }


def _schedule_error(probes) -> float:
    worst = 0.0
    for probe in probes[::12]:
        case = probe["case"]
        out = _public(case, probe["xyz"])
        eq = np.asarray([_schedule(e, case["h"])[0] for e in case["ells"]])
        ee = np.asarray([_schedule(e, case["h"])[1] for e in case["ells"]])
        for actual, expected in (
            (out["Q_source_by_label"], eq),
            (out["epsilon_source_by_label"], ee),
        ):
            worst = max(
                worst,
                float(
                    np.max(
                        np.abs(np.asarray(actual) - expected)
                        / np.maximum(np.abs(expected), 1e-300)
                    )
                ),
            )
    return worst


def _consistent_permutation(case, p) -> tuple[bool, float | None, bool | None]:
    contract = KokunoSourceMultiBandRealPairFamily(h=case["h"])
    base = _inputs(case, p)
    base_out = contract.physical_family(**base)
    baseline = base_out["velocity_physical_cartesian_total"][0]
    baseline_by_band = base_out["velocity_physical_cartesian_by_band"][0]
    perm = np.asarray((2, 1, 0))
    data = dict(base)
    for key in ("phase", "eta", "D_r_eta", "D_z_eta"):
        data[key] = np.take(data[key], perm, axis=-1)
    for key in ("n_phi", "t_plus", "D_r_C_plus", "D_z_C_plus"):
        data[key] = np.take(data[key], perm, axis=-2)
    data["beta_labels"] = tuple(base["beta_labels"][int(i)] for i in perm)
    try:
        moved_out = contract.physical_family(**data)
    except RuntimeError:
        return False, None, None
    moved = moved_out["velocity_physical_cartesian_total"][0]
    moved_by_band = moved_out["velocity_physical_cartesian_by_band"][0]
    relative = float(np.linalg.norm(moved - baseline)) / max(
        float(np.linalg.norm(baseline)), 1e-300
    )
    return True, relative, bool(np.array_equal(moved_by_band, baseline_by_band) and np.array_equal(moved, baseline))


def _mutations(probes) -> dict[str, Any]:
    h = FD4_STEPS[-1]
    sd, sr, ld, lr, perm = [], [], [], [], []
    failures = 0
    bitwise = []
    for probe in probes:
        case, p = probe["case"], probe["xyz"]
        correct = np.sum(_reference(case, p, h), axis=0)
        shared = np.sum(_reference(case, p, h, "shared_first"), axis=0)
        swapped = np.sum(_reference(case, p, h, "reversed"), axis=0)
        sd.append(float(np.linalg.norm(shared - correct)))
        sr.append(float(np.linalg.norm(correct)))
        ld.append(float(np.linalg.norm(swapped - correct)))
        lr.append(float(np.linalg.norm(correct)))
        ok, rel, exact = _consistent_permutation(case, p)
        failures += int(not ok)
        if rel is not None:
            perm.append(rel)
        if exact is not None:
            bitwise.append(exact)
    return {
        "shared_first_band_schedule_relative_rms": _rms(sd) / max(_rms(sr), 1e-300),
        "label_schedule_swap_relative_rms": _rms(ld) / max(_rms(lr), 1e-300),
        "consistent_public_label_permutation_failures": failures,
        "consistent_public_label_permutation_relative_max_when_evaluable": max(perm) if perm else None,
        "consistent_public_label_permutation_bitwise_all": bool(bitwise) and all(bitwise),
        "consistent_public_label_permutation_passed": failures == 0,
    }


def _ratio(a, b):
    return float(a / max(b, 1e-300))


def build_report() -> dict[str, Any]:
    probes = _samples()
    curl = [_curl_metrics(probes, h) for h in FD4_STEPS]
    div = [_div_metrics(probes, h) for h in FD4_STEPS]
    cr = [_ratio(curl[i]["curl_relative_rms"], curl[i + 1]["curl_relative_rms"]) for i in range(2)]
    dr = [
        _ratio(div[i]["normalized_divergence_rms"], div[i + 1]["normalized_divergence_rms"])
        for i in range(2)
    ]
    se = _schedule_error(probes)
    mut = _mutations(probes)
    checks = {
        "source_schedule": se <= FROZEN_GUARDS["source_schedule_relative_error_max"],
        "finest_cartesian_curl": curl[-1]["curl_relative_rms"] <= FROZEN_GUARDS["finest_cartesian_curl_relative_rms_max"],
        "curl_refinement": min(cr) >= FROZEN_GUARDS["curl_refinement_ratio_min"],
        "finest_normalized_divergence_rms": div[-1]["normalized_divergence_rms"] <= FROZEN_GUARDS["finest_normalized_divergence_rms_max"],
        "finest_normalized_divergence_point_max": div[-1]["normalized_divergence_point_max"] <= FROZEN_GUARDS["finest_normalized_divergence_point_max"],
        "divergence_refinement": min(dr) >= FROZEN_GUARDS["divergence_refinement_ratio_min"],
        "shared_schedule_mutation": mut["shared_first_band_schedule_relative_rms"] >= FROZEN_GUARDS["shared_first_band_schedule_mutation_relative_rms_min"],
        "label_swap_mutation": mut["label_schedule_swap_relative_rms"] >= FROZEN_GUARDS["label_schedule_swap_mutation_relative_rms_min"],
        "consistent_public_label_permutation": mut["consistent_public_label_permutation_passed"],
    }
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent2_parent_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "public_artifact": "KokunoSourceMultiBandRealPairFamily",
            "previous_agent4_pr": PREVIOUS_AGENT4_PR,
            "previous_agent4_head": PREVIOUS_AGENT4_HEAD,
        },
        "seed": SEED,
        "parameter_cases": [
            {
                **c,
                "schedules": [
                    {"ell": e, "Q": _schedule(e, c["h"])[0], "epsilon": _schedule(e, c["h"])[1], "k": _schedule(e, c["h"])[2]}
                    for e in c["ells"]
                ],
            }
            for c in CASES
        ],
        "sampling": {
            "total_points": len(probes),
            "off_grid_points": sum(p["region"] == "off_grid" for p in probes),
            "axis_near_points": sum(p["region"] == "axis_near" for p in probes),
            "minimum_nominal_radius": min(float(np.hypot(p["xyz"][0], p["xyz"][1])) for p in probes),
            "largest_cartesian_stencil_offset": 2 * FD4_STEPS[0],
        },
        "independent_operator": {
            "reference": "fresh Cartesian FD4 curl of independently reconstructed localized vector potential",
            "coefficient_derivatives_supplied_to_public_path": "fresh centered FD6 on independently reconstructed C_plus",
            "fd4_steps": list(FD4_STEPS),
            "coefficient_derivative_step": COEFF_FD6_STEP,
            "training_tensor_or_loss_read": False,
            "agent2_complete_curl_helper_used_by_reference": False,
            "agent2_residual_operator_reused": False,
            "free_forcing_used": False,
        },
        "source_schedule_relative_error_max": se,
        "curl_resolution_ladder": curl,
        "curl_refinement_ratios": cr,
        "divergence_resolution_ladder": div,
        "divergence_refinement_ratios": dr,
        "mutation": mut,
        "frozen_local_guards": dict(FROZEN_GUARDS),
        "checks": checks,
        "local_structural_preflight_passed": bool(all(checks.values())),
        "formal_project_gates": {
            "normalized_momentum_max_L2": FORMAL_MOMENTUM_GATE,
            "divergence_max_L2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
        },
        "truth_boundary": {
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "source_actual_partition_labels_instantiated": False,
            "public_xyz_t_velocity_correction_materialized": False,
            "global_leading_velocity_pressure_available": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "after_correction_global_velocity_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "limitations": [
            "Fresh phase/pulse/partition data are manufactured, not recovered source data.",
            "#461 still lacks public source-bound velocity(x,y,z,t) and actual positive-order/background binding.",
            "No leading pressure or materialized correction composite exists, so formal NS momentum is not run.",
            "A PASS here only closes the #451 deterministic aggregation/permutation seam under the unchanged frozen audit.",
        ],
    }
    payload["sha256"] = hashlib.sha256(_canonical(payload).encode()).hexdigest()
    return payload


def write_report(path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing audit report: {target}")
    target.write_text(json.dumps(build_report(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/kokuno_agent4/source_multiband_independent_audit.json")
    args = parser.parse_args(argv)
    path = write_report(args.output)
    report = json.loads(path.read_text())
    print(
        json.dumps(
            {
                "local_structural_preflight_passed": report["local_structural_preflight_passed"],
                "checks": report["checks"],
                "finest_curl_relative_rms": report["curl_resolution_ladder"][-1]["curl_relative_rms"],
                "finest_normalized_divergence_rms": report["divergence_resolution_ladder"][-1]["normalized_divergence_rms"],
                "finest_normalized_divergence_point_max": report["divergence_resolution_ladder"][-1]["normalized_divergence_point_max"],
                "mutation": report["mutation"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
