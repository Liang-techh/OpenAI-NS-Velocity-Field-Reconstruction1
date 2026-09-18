"""Independent Agent-4 audit of the source-scheduled Kokuno multi-band family.

Agent 2 PR #450 is the first current handoff that evaluates at least two
source-scheduled dyadic bands with a distinct Q_ell/epsilon_ell for every beta
label. Its focused regression compares against the same complete-curl helper
used by the implementation. This module instead treats the public
``physical_family`` output as a black box and compares it with a separately
written Cartesian vector-potential curl oracle on fresh held-out data.

The manufactured phase, transverse pulse and slow partition are deliberately
fresh and analytic. They are source-compatible test data, not recovered
positive-order/background or auxiliary-torus source data. Consequently this
is a local structural preflight only; it does not run the formal full-domain
Navier--Stokes gate.
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

SCHEMA = "kokuno-agent4-source-multiband-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-MULTIBAND-INDEPENDENT-AUDIT-027"
BASE_PR = 450
BASE_HEAD = "e2e583a0388a177b5b02a51b4d1f55a8a7da1c68"
SEED = 9173171
FD4_STEPS = (0.032, 0.016, 0.008)
COEFFICIENT_DERIVATIVE_STEP = 2.0e-4
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# Frozen before execution. These are local implementation-consistency guards,
# not replacements for the formal PDE thresholds above.
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
    # ell*h crosses 2 here, so ceil(epsilon^-1/2) changes inside the family.
    {"name": "k_transition", "h": 0.0090, "ells": (221, 223, 225)},
)


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _schedule(ell: int, h: float) -> tuple[float, float, int]:
    q = float(2.0 ** (-int(ell)))
    epsilon = float(q**h)
    return q, epsilon, int(ceil(epsilon ** -0.5))


def _fd6_first(fn: Callable[[np.ndarray], np.ndarray], x: np.ndarray, h: float) -> np.ndarray:
    return (
        -fn(x - 3.0 * h)
        + 9.0 * fn(x - 2.0 * h)
        - 45.0 * fn(x - h)
        + 45.0 * fn(x + h)
        - 9.0 * fn(x + 2.0 * h)
        + fn(x + 3.0 * h)
    ) / (60.0 * h)


def _fd4_axis(
    fn: Callable[[np.ndarray], np.ndarray], point: np.ndarray, axis: int, step: float
) -> np.ndarray:
    m2 = np.asarray(point, dtype=float).copy()
    m1 = m2.copy()
    p1 = m2.copy()
    p2 = m2.copy()
    m2[axis] -= 2.0 * step
    m1[axis] -= step
    p1[axis] += step
    p2[axis] += 2.0 * step
    return (fn(m2) - 8.0 * fn(m1) + 8.0 * fn(p1) - fn(p2)) / (12.0 * step)


def _partition(r: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alpha = 0.72 + 0.11 * r - 0.07 * z + 0.025 * r * z
    beta = 0.43 - 0.06 * r + 0.09 * z + 0.02 * r * r
    ar = 0.11 + 0.025 * z
    az = -0.07 + 0.025 * r
    br = -0.06 + 0.04 * r
    bz = np.full_like(r, 0.09)

    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    eta = np.stack((ca, sa * cb, sa * sb), axis=-1)
    dr = np.stack(
        (
            -sa * ar,
            ca * ar * cb - sa * sb * br,
            ca * ar * sb + sa * cb * br,
        ),
        axis=-1,
    )
    dz = np.stack(
        (
            -sa * az,
            ca * az * cb - sa * sb * bz,
            ca * az * sb + sa * cb * bz,
        ),
        axis=-1,
    )
    return eta, dr, dz


def _phase_n_raw(
    ell: int,
    h: float,
    label_index: int,
    r: np.ndarray,
    theta: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _, _, k = _schedule(ell, h)
    angular_integer = (1, 2, 1)[label_index]
    p = angular_integer / k
    a = 0.83 + 0.09 * label_index
    c = 0.47 - 0.06 * label_index
    d = 0.035 * (label_index + 1)
    e = 0.012 * (label_index + 1)
    f = -0.010 * (label_index + 1)
    phase = p * theta + a * r + c * z + d * r * z + e * r * r + f * z * z
    n = np.stack(
        (
            a + d * z + 2.0 * e * r,
            np.full_like(r, p) / r,
            c + d * r + 2.0 * f * z,
        ),
        axis=-1,
    )
    raw = np.stack(
        (
            (0.61 + 0.025 * r - 0.017 * z)
            + 1j * (0.08 + 0.013 * r + 0.007 * z),
            (-0.31 + 0.021 * z + 0.011 * r)
            + 1j * (0.12 - 0.009 * r + 0.006 * z),
            (0.24 + 0.018 * r * z - 0.008 * r)
            + 1j * (-0.07 + 0.010 * z - 0.005 * r),
        ),
        axis=-1,
    )
    return phase, n, raw


def _transverse_and_coefficient(
    ell: int,
    h: float,
    label_index: int,
    r: np.ndarray,
    theta: np.ndarray,
    z: np.ndarray,
    *,
    k_override: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phase, n, raw = _phase_n_raw(ell, h, label_index, r, theta, z)
    norm_sq = np.sum(n * n, axis=-1)
    projection = np.sum(n * raw, axis=-1) / norm_sq
    t = raw - n * projection[..., None]
    k = _schedule(ell, h)[2] if k_override is None else int(k_override)
    coefficient = 1j * np.cross(n, t) / (k * norm_sq[..., None])
    return phase, n, t, coefficient


def _coefficient_derivatives(
    ell: int,
    h: float,
    label_index: int,
    r: np.ndarray,
    theta: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    step = COEFFICIENT_DERIVATIVE_STEP
    dr = _fd6_first(
        lambda rr: _transverse_and_coefficient(ell, h, label_index, rr, theta, z)[3],
        r,
        step,
    )
    dz = _fd6_first(
        lambda zz: _transverse_and_coefficient(ell, h, label_index, r, theta, zz)[3],
        z,
        step,
    )
    return dr, dz


def _inputs(case: dict[str, Any], xyz: np.ndarray) -> dict[str, Any]:
    points = np.asarray(xyz, dtype=float)
    if points.ndim == 1:
        points = points[None, :]
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    r = np.hypot(x, y)
    if np.any(r <= 0.0):
        raise ValueError("audit points must remain away from the cylindrical axis")
    theta = np.arctan2(y, x)
    eta, dr_eta, dz_eta = _partition(r, z)

    phase_cols = []
    n_cols = []
    t_cols = []
    dr_cols = []
    dz_cols = []
    for j, ell in enumerate(case["ells"]):
        phase, n, t, _ = _transverse_and_coefficient(ell, case["h"], j, r, theta, z)
        dr_c, dz_c = _coefficient_derivatives(ell, case["h"], j, r, theta, z)
        phase_cols.append(phase)
        n_cols.append(n)
        t_cols.append(t)
        dr_cols.append(dr_c)
        dz_cols.append(dz_c)

    return {
        "R": r,
        "theta": theta,
        "phase": np.stack(phase_cols, axis=-1),
        "n_phi": np.stack(n_cols, axis=-2),
        "t_plus": np.stack(t_cols, axis=-2),
        "D_r_C_plus": np.stack(dr_cols, axis=-2),
        "D_z_C_plus": np.stack(dz_cols, axis=-2),
        "eta": eta,
        "D_r_eta": dr_eta,
        "D_z_eta": dz_eta,
        "beta_labels": tuple((int(ell), (j, 0, 0)) for j, ell in enumerate(case["ells"])),
    }


def _cyl_to_cart(vector: np.ndarray, theta: float) -> np.ndarray:
    v = np.asarray(vector, dtype=np.complex128)
    c, s = np.cos(theta), np.sin(theta)
    return np.asarray((v[0] * c - v[1] * s, v[0] * s + v[1] * c, v[2]))


def _direct_plus_potential(
    case: dict[str, Any], label_index: int, point: np.ndarray, *, k_override: int | None = None
) -> np.ndarray:
    x, y, z = np.asarray(point, dtype=float)
    r = float(np.hypot(x, y))
    if r <= 0.0:
        raise ValueError("Cartesian curl stencil crossed the axis")
    theta = float(np.arctan2(y, x))
    ell = int(case["ells"][label_index])
    phase, _, _, coefficient = _transverse_and_coefficient(
        ell,
        case["h"],
        label_index,
        np.asarray([r]),
        np.asarray([theta]),
        np.asarray([z]),
        k_override=k_override,
    )
    eta = _partition(np.asarray([r]), np.asarray([z]))[0][0, label_index]
    k = _schedule(ell, case["h"])[2] if k_override is None else int(k_override)
    potential_cyl = eta * coefficient[0] * np.exp(1j * k * phase[0])
    return _cyl_to_cart(potential_cyl, theta)


def _cartesian_curl_fd4(
    fn: Callable[[np.ndarray], np.ndarray], point: np.ndarray, step: float
) -> np.ndarray:
    jac = np.empty((3, 3), dtype=np.complex128)
    for axis in range(3):
        jac[:, axis] = _fd4_axis(fn, point, axis, step)
    return np.asarray(
        (
            jac[2, 1] - jac[1, 2],
            jac[0, 2] - jac[2, 0],
            jac[1, 0] - jac[0, 1],
        ),
        dtype=np.complex128,
    )


def _reference_by_beta(
    case: dict[str, Any], point: np.ndarray, step: float, *, shared_first_schedule: bool = False
) -> np.ndarray:
    output = []
    first_q, _, first_k = _schedule(case["ells"][0], case["h"])
    A = 0.5 + case["h"]
    for j, ell in enumerate(case["ells"]):
        q, _, k = _schedule(ell, case["h"])
        use_q = first_q if shared_first_schedule else q
        use_k = first_k if shared_first_schedule else k
        curl = _cartesian_curl_fd4(
            lambda p, jj=j, kk=use_k: _direct_plus_potential(case, jj, p, k_override=kk),
            point,
            step,
        )
        output.append((use_q ** (-A)) * (2.0 * curl.real))
    return np.stack(output, axis=0)


def _public(
    case: dict[str, Any], xyz: np.ndarray, *, labels: tuple[Any, ...] | None = None
) -> dict[str, Any]:
    contract = KokunoSourceMultiBandRealPairFamily(h=case["h"])
    data = _inputs(case, xyz)
    if labels is not None:
        data["beta_labels"] = labels
    return contract.physical_family(**data)


def _sample_points() -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    probes: list[dict[str, Any]] = []
    for case in CASES:
        radii = np.concatenate((rng.uniform(0.42, 1.30, 8), rng.uniform(0.14, 0.24, 4)))
        theta = rng.uniform(-2.7, 2.7, radii.size)
        z = rng.uniform(-0.55, 0.55, radii.size)
        for idx in range(radii.size):
            probes.append(
                {
                    "case": case,
                    "region": "axis_near" if idx >= 8 else "off_grid",
                    "xyz": np.asarray(
                        (
                            radii[idx] * np.cos(theta[idx]),
                            radii[idx] * np.sin(theta[idx]),
                            z[idx],
                        ),
                        dtype=float,
                    ),
                }
            )
    return probes


def _rms(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _evaluate_curl_step(probes: list[dict[str, Any]], step: float) -> dict[str, Any]:
    errors: list[float] = []
    refs: list[float] = []
    point_rel: list[float] = []
    case_rel: dict[str, list[tuple[float, float]]] = {}
    for probe in probes:
        case = probe["case"]
        point = probe["xyz"]
        public = _public(case, point)["velocity_physical_cartesian_by_beta"][0]
        reference = _reference_by_beta(case, point, step)
        diff = float(np.linalg.norm(public - reference))
        ref = float(np.linalg.norm(reference))
        errors.append(diff)
        refs.append(ref)
        point_rel.append(diff / max(ref, 1.0e-300))
        case_rel.setdefault(case["name"], []).append((diff, ref))
    result: dict[str, Any] = {
        "step": step,
        "curl_error_rms": _rms(errors),
        "curl_reference_rms": _rms(refs),
        "curl_relative_rms": _rms(errors) / max(_rms(refs), 1.0e-300),
        "curl_point_relative_max": float(max(point_rel)),
        "parameter_cases": {},
    }
    for name, pairs in case_rel.items():
        e = [p[0] for p in pairs]
        r = [p[1] for p in pairs]
        result["parameter_cases"][name] = _rms(e) / max(_rms(r), 1.0e-300)
    return result


def _public_total(case: dict[str, Any], point: np.ndarray) -> np.ndarray:
    return np.asarray(_public(case, point)["velocity_physical_cartesian_total"][0], dtype=float)


def _divergence_metrics(probes: list[dict[str, Any]], step: float) -> dict[str, Any]:
    divergences: list[float] = []
    grad_norms: list[float] = []
    normalized_points: list[float] = []
    regions: dict[str, list[float]] = {"off_grid": [], "axis_near": []}
    for probe in probes:
        case = probe["case"]
        point = probe["xyz"]
        jac = np.empty((3, 3), dtype=float)
        for axis in range(3):
            jac[:, axis] = _fd4_axis(lambda p: _public_total(case, p), point, axis, step)
        divergence = float(abs(np.trace(jac)))
        grad_norm = float(np.linalg.norm(jac))
        normalized = divergence / max(grad_norm, 1.0e-300)
        divergences.append(divergence)
        grad_norms.append(grad_norm)
        normalized_points.append(normalized)
        regions[probe["region"]].append(normalized)
    return {
        "step": step,
        "divergence_abs_rms": _rms(divergences),
        "gradient_frobenius_rms": _rms(grad_norms),
        "normalized_divergence_rms": _rms(divergences) / max(_rms(grad_norms), 1.0e-300),
        "normalized_divergence_point_max": float(max(normalized_points)),
        "regions": {
            name: {"normalized_point_max": float(max(values))}
            for name, values in regions.items()
        },
    }


def _schedule_error(probes: list[dict[str, Any]]) -> float:
    worst = 0.0
    for probe in probes[::12]:
        case = probe["case"]
        out = _public(case, probe["xyz"])
        expected_q = np.asarray([_schedule(ell, case["h"])[0] for ell in case["ells"]])
        expected_eps = np.asarray([_schedule(ell, case["h"])[1] for ell in case["ells"]])
        for actual, expected in (
            (np.asarray(out["Q_source_by_label"]), expected_q),
            (np.asarray(out["epsilon_source_by_label"]), expected_eps),
        ):
            worst = max(
                worst,
                float(
                    np.max(
                        np.abs(actual - expected) / np.maximum(np.abs(expected), 1.0e-300)
                    )
                ),
            )
    return worst


def _mutation_metrics(probes: list[dict[str, Any]]) -> dict[str, float]:
    step = FD4_STEPS[-1]
    shared_diffs: list[float] = []
    shared_refs: list[float] = []
    swap_diffs: list[float] = []
    swap_refs: list[float] = []
    for probe in probes:
        case = probe["case"]
        point = probe["xyz"]
        correct_ref = np.sum(_reference_by_beta(case, point, step), axis=0)
        shared = np.sum(
            _reference_by_beta(case, point, step, shared_first_schedule=True), axis=0
        )
        shared_diffs.append(float(np.linalg.norm(shared - correct_ref)))
        shared_refs.append(float(np.linalg.norm(correct_ref)))

        base = _public(case, point)["velocity_physical_cartesian_total"][0]
        labels = tuple((int(ell), (j, 0, 0)) for j, ell in enumerate(case["ells"]))
        swapped = (labels[-1], labels[1], labels[0])
        mutated = _public(case, point, labels=swapped)[
            "velocity_physical_cartesian_total"
        ][0]
        swap_diffs.append(float(np.linalg.norm(mutated - base)))
        swap_refs.append(float(np.linalg.norm(base)))
    return {
        "shared_first_band_schedule_relative_rms": _rms(shared_diffs)
        / max(_rms(shared_refs), 1.0e-300),
        "label_schedule_swap_relative_rms": _rms(swap_diffs)
        / max(_rms(swap_refs), 1.0e-300),
    }


def _ratio(coarse: float, fine: float) -> float:
    return float(coarse / max(fine, 1.0e-300))


def build_report() -> dict[str, Any]:
    probes = _sample_points()
    curl_ladder = [_evaluate_curl_step(probes, step) for step in FD4_STEPS]
    div_ladder = [_divergence_metrics(probes, step) for step in FD4_STEPS]
    curl_ratios = [
        _ratio(curl_ladder[i]["curl_relative_rms"], curl_ladder[i + 1]["curl_relative_rms"])
        for i in range(2)
    ]
    div_ratios = [
        _ratio(
            div_ladder[i]["normalized_divergence_rms"],
            div_ladder[i + 1]["normalized_divergence_rms"],
        )
        for i in range(2)
    ]
    schedule_error = _schedule_error(probes)
    mutation = _mutation_metrics(probes)
    finest_curl = curl_ladder[-1]
    finest_div = div_ladder[-1]
    checks = {
        "source_schedule": schedule_error
        <= FROZEN_GUARDS["source_schedule_relative_error_max"],
        "finest_cartesian_curl": finest_curl["curl_relative_rms"]
        <= FROZEN_GUARDS["finest_cartesian_curl_relative_rms_max"],
        "curl_refinement": min(curl_ratios)
        >= FROZEN_GUARDS["curl_refinement_ratio_min"],
        "finest_normalized_divergence_rms": finest_div["normalized_divergence_rms"]
        <= FROZEN_GUARDS["finest_normalized_divergence_rms_max"],
        "finest_normalized_divergence_point_max": finest_div[
            "normalized_divergence_point_max"
        ]
        <= FROZEN_GUARDS["finest_normalized_divergence_point_max"],
        "divergence_refinement": min(div_ratios)
        >= FROZEN_GUARDS["divergence_refinement_ratio_min"],
        "shared_schedule_mutation": mutation[
            "shared_first_band_schedule_relative_rms"
        ]
        >= FROZEN_GUARDS["shared_first_band_schedule_mutation_relative_rms_min"],
        "label_swap_mutation": mutation["label_schedule_swap_relative_rms"]
        >= FROZEN_GUARDS["label_schedule_swap_mutation_relative_rms_min"],
    }
    passed = bool(all(checks.values()))
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent2_parent_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "public_artifact": "KokunoSourceMultiBandRealPairFamily",
        },
        "seed": SEED,
        "parameter_cases": [
            {
                **case,
                "schedules": [
                    {
                        "ell": ell,
                        "Q": _schedule(ell, case["h"])[0],
                        "epsilon": _schedule(ell, case["h"])[1],
                        "k": _schedule(ell, case["h"])[2],
                    }
                    for ell in case["ells"]
                ],
            }
            for case in CASES
        ],
        "sampling": {
            "total_points": len(probes),
            "off_grid_points": sum(p["region"] == "off_grid" for p in probes),
            "axis_near_points": sum(p["region"] == "axis_near" for p in probes),
            "minimum_nominal_radius": min(
                float(np.hypot(p["xyz"][0], p["xyz"][1])) for p in probes
            ),
            "largest_cartesian_stencil_offset": 2.0 * FD4_STEPS[0],
        },
        "independent_operator": {
            "reference": (
                "fresh Cartesian FD4 curl of independently reconstructed localized vector "
                "potential; no Agent-2 complete-curl helper"
            ),
            "coefficient_derivatives_supplied_to_public_path": (
                "fresh centered FD6 on independently reconstructed C_plus"
            ),
            "fd4_steps": list(FD4_STEPS),
            "coefficient_derivative_step": COEFFICIENT_DERIVATIVE_STEP,
            "training_tensor_or_loss_read": False,
            "agent2_complete_curl_helper_used_by_reference": False,
            "agent2_residual_operator_reused": False,
            "free_forcing_used": False,
        },
        "source_schedule_relative_error_max": schedule_error,
        "curl_resolution_ladder": curl_ladder,
        "curl_refinement_ratios": curl_ratios,
        "divergence_resolution_ladder": div_ladder,
        "divergence_refinement_ratios": div_ratios,
        "mutation": mutation,
        "frozen_local_guards": dict(FROZEN_GUARDS),
        "checks": checks,
        "local_structural_preflight_passed": passed,
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
            "The phase, transverse pulse and three-band squared partition are fresh manufactured source-compatible data, not recovered source data.",
            "Agent 2 #450 still has no public source-bound velocity(x,y,z,t) and no actual positive-order/background binding.",
            "No leading pressure or materialized correction composite is available, so the formal NS momentum gate is not run.",
        ],
    }
    payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return payload


def write_report(path: str | Path) -> Path:
    target = Path(path)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing audit report: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    report = build_report()
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/source_multiband_independent_audit.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    report = json.loads(path.read_text(encoding="utf-8"))
    print(
        json.dumps(
            {
                "local_structural_preflight_passed": report[
                    "local_structural_preflight_passed"
                ],
                "finest_curl_relative_rms": report["curl_resolution_ladder"][-1][
                    "curl_relative_rms"
                ],
                "finest_normalized_divergence_rms": report[
                    "divergence_resolution_ladder"
                ][-1]["normalized_divergence_rms"],
                "finest_normalized_divergence_point_max": report[
                    "divergence_resolution_ladder"
                ][-1]["normalized_divergence_point_max"],
                "mutation": report["mutation"],
            },
            sort_keys=True,
        )
    )
    return 0 if report["local_structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
