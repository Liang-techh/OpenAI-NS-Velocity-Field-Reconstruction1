"""Independent numerical audit of the source-compatible PA.10 axis domain.

This is the repaired Agent-4 audit for Agent-1 #644.  The first A4 attempt
(#646, head ab567b5..., workflow 35447156458) failed before producing a
scientific receipt because an 8192-point *uniform* random cloud happened to
contain no point in the narrow set K_delta={|Z_*|<=delta_*}.  That was a
sampling-design failure, not evidence that the source PA.8 implication failed.

The repaired protocol keeps every source parameter and scientific threshold
fixed.  It retains the three dense independent real-axis ladders and adds a
stratified off-grid check: the finest independent grid is used only to locate a
point deep inside each rare region, then a fresh seeded local cloud is drawn
around that point.  The public upstream certificate helpers are never called.

The result is still only a leading-side prerequisite audit.  It is not a
Navier--Stokes momentum residual and cannot set ``pde_validated``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Literal

import numpy as np

from .kokuno_pa10_source_axis_domain import KokunoPA10SourceCompatibleAxisDomain

SCHEMA = "kokuno-agent4-pa10-source-axis-independent-audit-v2"
AGENT1_PR = 644
AGENT1_HEAD = "1bd16378651f6441ec71c6cf1ae4c3abccb70cc5"
PRIOR_A4_PR = 646
PRIOR_A4_HEAD = "ab567b5325c54f676d35dd2a70d5215105dc0266"
PRIOR_FAILED_WORKFLOW = 35447156458

REAL_LEVELS = (4097, 8193, 16385)
TARGETED_OFFGRID_SEED = 9173341
TARGETED_OFFGRID_COUNT = 256
TARGETED_DRAW_BATCH = 8192
TARGETED_MAX_SHRINK_LEVELS = 16
OFFGRID_CELL_DISTANCE_FLOOR = 1.0e-10
COMPLEX_TUBE_MUTATION_FACTOR = 100.0
G_BOUND_MUTATION_FACTOR = 0.99


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _segment_distance(z: complex, lo: float, hi: float) -> float:
    x = float(np.real(z))
    y = float(np.imag(z))
    if x < lo:
        return math.hypot(lo - x, y)
    if x > hi:
        return math.hypot(x - hi, y)
    return abs(y)


def _H(eta: np.ndarray | complex, *, D: float, j0: float):
    return j0 + (D + 4.0) * eta - j0 * eta * eta - 4.0 * eta * eta * eta


def _Z(
    eta: np.ndarray,
    *,
    A: float,
    D: float,
    j0: float,
    pressure_square: float,
) -> np.ndarray:
    e = np.asarray(eta, dtype=float)
    d = 1.0 - e * e
    U = 4.0 * e + j0
    H = _H(e, D=D, j0=j0)
    B = 1.0 - 2.0 * e * U
    den = 1.0 + e * e
    Pi = -pressure_square / den**2
    Pi_eta = 4.0 * pressure_square * e / den**3
    return -A * B * U - 4.0 * H - d * Pi_eta + 4.0 * A * e * Pi


def _H_shift_roots(*, D: float, j0: float, target: complex) -> np.ndarray:
    coeff = np.asarray([-4.0, -j0, D + 4.0, j0 - target], dtype=complex)
    return np.roots(coeff)


def _independent_g_norm(x: float) -> float:
    return (1.0 + x) / (1.0 - x) ** 3


def _point_satisfies_pa8_separation(
    H: float, Z: float, *, sigma: float, delta: float
) -> bool:
    if abs(Z) <= delta:
        return abs(H) > 10.0 * sigma
    return True


def _targeted_offgrid_region(
    *,
    region: Literal["K_delta", "H_small"],
    lo: float,
    hi: float,
    datum: KokunoPA10SourceCompatibleAxisDomain,
    dense_eta: np.ndarray,
    dense_H: np.ndarray,
    dense_Z: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build a non-vacuous fresh off-grid cloud for one narrow PA.8 region.

    The dense grid is used only to locate a deep interior anchor by minimizing
    |Z| (for K_delta) or |H| (for H_small).  No upstream finite-cover or axis
    certificate is consulted.  Around the anchor, the sampling half-width is
    reduced geometrically until a fresh seeded cloud contains at least the
    preregistered target number of region points.  Exact dense-grid nodes are
    excluded explicitly.
    """

    if dense_eta.ndim != 1 or dense_eta.size != REAL_LEVELS[-1]:
        raise ValueError("targeted off-grid locator requires the frozen finest grid")
    spacing = float((hi - lo) / (dense_eta.size - 1))
    if region == "K_delta":
        anchor_index = int(np.argmin(np.abs(dense_Z)))
    elif region == "H_small":
        anchor_index = int(np.argmin(np.abs(dense_H)))
    else:  # pragma: no cover - Literal plus defensive fail-closed branch
        raise ValueError("unknown targeted region")
    anchor = float(dense_eta[anchor_index])

    rng = np.random.default_rng(seed)
    selected: np.ndarray | None = None
    selected_level: int | None = None
    selected_halfwidth: float | None = None
    selected_acceptance_count: int | None = None
    selected_min_cell_distance: float | None = None

    for level in range(TARGETED_MAX_SHRINK_LEVELS):
        halfwidth = spacing * (0.5**level)
        candidates = anchor + rng.uniform(
            -halfwidth, halfwidth, size=TARGETED_DRAW_BATCH
        )
        candidates = candidates[(candidates > lo) & (candidates < hi)]
        if candidates.size == 0:
            continue
        cell_coordinate = (candidates - lo) / spacing
        offgrid_cell_distance = np.abs(cell_coordinate - np.rint(cell_coordinate))
        offgrid = offgrid_cell_distance > OFFGRID_CELL_DISTANCE_FLOOR

        H = np.asarray(_H(candidates, D=datum.D, j0=datum.j0), dtype=float)
        Z = _Z(
            candidates,
            A=datum.A,
            D=datum.D,
            j0=datum.j0,
            pressure_square=datum.pressure_square,
        )
        if region == "K_delta":
            inside = np.abs(Z) <= datum.delta_star
        else:
            inside = np.abs(H) <= 10.0 * datum.sigma_star
        accepted = np.unique(candidates[inside & offgrid])
        if accepted.size >= TARGETED_OFFGRID_COUNT:
            selected = accepted[:TARGETED_OFFGRID_COUNT]
            selected_level = level
            selected_halfwidth = halfwidth
            selected_acceptance_count = int(accepted.size)
            selected_cells = (selected - lo) / spacing
            selected_min_cell_distance = float(
                np.min(np.abs(selected_cells - np.rint(selected_cells)))
            )
            break

    if selected is None:
        raise RuntimeError(
            f"could not materialize {TARGETED_OFFGRID_COUNT} off-grid points in {region}"
        )

    return selected, {
        "region": region,
        "seed": seed,
        "target_count": TARGETED_OFFGRID_COUNT,
        "selected_count": int(selected.size),
        "anchor_eta": anchor,
        "anchor_dense_index": anchor_index,
        "dense_spacing": spacing,
        "selected_shrink_level": int(selected_level),
        "selected_halfwidth": float(selected_halfwidth),
        "accepted_before_truncation": int(selected_acceptance_count),
        "minimum_distance_from_dense_grid_in_cells": float(selected_min_cell_distance),
        "all_points_strictly_off_dense_grid": bool(
            selected_min_cell_distance > OFFGRID_CELL_DISTANCE_FLOOR
        ),
    }


def run_audit() -> dict[str, Any]:
    datum = KokunoPA10SourceCompatibleAxisDomain()
    lo = -1.0 - datum.enlarged_real_margin
    hi = 1.0 + datum.enlarged_real_margin
    sigma = float(datum.sigma_star)
    delta = float(datum.delta_star)

    real_levels: list[dict[str, Any]] = []
    finest_eta: np.ndarray | None = None
    finest_H: np.ndarray | None = None
    finest_Z: np.ndarray | None = None
    for count in REAL_LEVELS:
        eta = np.linspace(lo, hi, count, dtype=float)
        H = np.asarray(_H(eta, D=datum.D, j0=datum.j0), dtype=float)
        Z = _Z(
            eta,
            A=datum.A,
            D=datum.D,
            j0=datum.j0,
            pressure_square=datum.pressure_square,
        )
        K = np.abs(Z) <= delta
        Hsmall = np.abs(H) <= 10.0 * sigma
        chi = H * H / (H * H + sigma * sigma)
        real_levels.append(
            {
                "count": count,
                "K_count": int(np.count_nonzero(K)),
                "Hsmall_count": int(np.count_nonzero(Hsmall)),
                "min_abs_H_on_K": float(np.min(np.abs(H[K]))) if np.any(K) else None,
                "min_chi_on_K": float(np.min(chi[K])) if np.any(K) else None,
                "min_Z_on_Hsmall": float(np.min(Z[Hsmall])) if np.any(Hsmall) else None,
                "max_abs_Z_on_K": float(np.max(np.abs(Z[K]))) if np.any(K) else None,
            }
        )
        if count == REAL_LEVELS[-1]:
            finest_eta, finest_H, finest_Z = eta, H, Z

    if finest_eta is None or finest_H is None or finest_Z is None:
        raise RuntimeError("finest real-axis grid was not materialized")
    finest = real_levels[-1]
    real_nonvacuous = int(finest["K_count"]) > 0 and int(finest["Hsmall_count"]) > 0
    real_pa8_pass = bool(
        real_nonvacuous
        and float(finest["min_abs_H_on_K"]) > 10.0 * sigma
        and float(finest["min_chi_on_K"]) > 0.99
        and float(finest["min_Z_on_Hsmall"]) > 2.0 * delta
    )

    K_eta, K_meta = _targeted_offgrid_region(
        region="K_delta",
        lo=lo,
        hi=hi,
        datum=datum,
        dense_eta=finest_eta,
        dense_H=finest_H,
        dense_Z=finest_Z,
        seed=TARGETED_OFFGRID_SEED,
    )
    H_eta, H_meta = _targeted_offgrid_region(
        region="H_small",
        lo=lo,
        hi=hi,
        datum=datum,
        dense_eta=finest_eta,
        dense_H=finest_H,
        dense_Z=finest_Z,
        seed=TARGETED_OFFGRID_SEED + 1,
    )

    K_H = np.asarray(_H(K_eta, D=datum.D, j0=datum.j0), dtype=float)
    K_Z = _Z(
        K_eta,
        A=datum.A,
        D=datum.D,
        j0=datum.j0,
        pressure_square=datum.pressure_square,
    )
    K_chi = K_H * K_H / (K_H * K_H + sigma * sigma)
    H_H = np.asarray(_H(H_eta, D=datum.D, j0=datum.j0), dtype=float)
    H_Z = _Z(
        H_eta,
        A=datum.A,
        D=datum.D,
        j0=datum.j0,
        pressure_square=datum.pressure_square,
    )
    targeted_pass = bool(
        K_eta.size == TARGETED_OFFGRID_COUNT
        and H_eta.size == TARGETED_OFFGRID_COUNT
        and np.all(np.abs(K_Z) <= delta)
        and np.min(np.abs(K_H)) > 10.0 * sigma
        and np.min(K_chi) > 0.99
        and np.all(np.abs(H_H) <= 10.0 * sigma)
        and np.min(H_Z) > 2.0 * delta
        and K_meta["all_points_strictly_off_dense_grid"]
        and H_meta["all_points_strictly_off_dense_grid"]
    )

    h_roots = np.concatenate(
        [
            _H_shift_roots(D=datum.D, j0=datum.j0, target=1j * sigma),
            _H_shift_roots(D=datum.D, j0=datum.j0, target=-1j * sigma),
        ]
    )
    h_root_distances = np.asarray(
        [_segment_distance(complex(z), lo, hi) for z in h_roots], dtype=float
    )
    min_h_root_distance = float(np.min(h_root_distances))
    tube = float(datum.complex_tube_radius)

    L_root = 1.0 / math.sqrt(2.0 * datum.h)
    L_root_distance = min(
        _segment_distance(complex(L_root, 0.0), lo, hi),
        _segment_distance(complex(-L_root, 0.0), lo, hi),
    )
    pi_root_distance = min(
        _segment_distance(1j, lo, hi),
        _segment_distance(-1j, lo, hi),
    )
    complex_pass = bool(
        min_h_root_distance > tube
        and L_root_distance > tube
        and pi_root_distance > tube
    )

    x = float(datum.coefficient_rho / datum.cauchy_radius)
    exact_g_series = _independent_g_norm(x)
    public_g = float(datum.pressure_operator_inputs()["g_norm_upper"])
    partials = {
        str(N): float(sum((beta + 1) ** 2 * x**beta for beta in range(N + 1)))
        for N in (8, 16, 32, 64)
    }
    g_pass = bool(public_g + 2.0e-15 >= exact_g_series)

    mutation = {
        "complex_tube_x100_detected": bool(
            min_h_root_distance <= COMPLEX_TUBE_MUTATION_FACTOR * tube
        ),
        "g_bound_0p99_detected": bool(
            G_BOUND_MUTATION_FACTOR * public_g < exact_g_series
        ),
        "synthetic_PA8_bad_point_detected": bool(
            not _point_satisfies_pa8_separation(
                0.0, 0.0, sigma=sigma, delta=delta
            )
        ),
    }

    failed: list[str] = []
    if not real_pa8_pass:
        failed.append("independent_real_axis_PA8")
    if not targeted_pass:
        failed.append("targeted_offgrid_PA8")
    if not complex_pass:
        failed.append("independent_complex_zero_exclusion")
    if not g_pass:
        failed.append("independent_g_coefficient_norm")
    if not all(mutation.values()):
        failed.append("mutation_detection")

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "agent1_pr": AGENT1_PR,
            "agent1_head": AGENT1_HEAD,
            "source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "source_path": "navier-stokes/navier_stokes_workbench.tex",
        },
        "protocol_repair": {
            "prior_agent4_pr": PRIOR_A4_PR,
            "prior_agent4_head": PRIOR_A4_HEAD,
            "prior_failed_workflow": PRIOR_FAILED_WORKFLOW,
            "prior_failure_stage": "focused warnings-as-errors regression before receipt generation",
            "prior_failure": "uniform 8192-point off-grid cloud contained zero K_delta samples",
            "repair": "stratified fresh off-grid sampling localized by the independent finest grid",
            "source_parameters_changed": False,
            "scientific_thresholds_changed": False,
            "final_project_gates_changed": False,
        },
        "frozen_protocol": {
            "real_axis_levels": list(REAL_LEVELS),
            "targeted_offgrid_seed": TARGETED_OFFGRID_SEED,
            "targeted_offgrid_count_per_region": TARGETED_OFFGRID_COUNT,
            "targeted_draw_batch": TARGETED_DRAW_BATCH,
            "targeted_max_shrink_levels": TARGETED_MAX_SHRINK_LEVELS,
            "complex_zero_method": "numpy polynomial roots of H_*(z)=+/- i sigma_*",
            "complex_tube_mutation_factor": COMPLEX_TUBE_MUTATION_FACTOR,
            "g_bound_mutation_factor": G_BOUND_MUTATION_FACTOR,
            "final_project_gates_unchanged": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
            },
        },
        "public_parameters": {
            "sigma_star": sigma,
            "delta_star": delta,
            "enlarged_real_margin": float(datum.enlarged_real_margin),
            "complex_tube_radius": tube,
            "cauchy_radius": float(datum.cauchy_radius),
            "coefficient_rho": float(datum.coefficient_rho),
            "h": float(datum.h),
            "A": float(datum.A),
            "D": float(datum.D),
            "j0": float(datum.j0),
            "pressure_square": float(datum.pressure_square),
        },
        "real_axis": {
            "levels": real_levels,
            "nonvacuous_K_and_Hsmall": real_nonvacuous,
            "PA8_independent_sampling_passed": real_pa8_pass,
        },
        "targeted_offgrid": {
            "K_delta": {
                **K_meta,
                "max_abs_Z": float(np.max(np.abs(K_Z))),
                "min_abs_H": float(np.min(np.abs(K_H))),
                "min_chi": float(np.min(K_chi)),
            },
            "H_small": {
                **H_meta,
                "max_abs_H": float(np.max(np.abs(H_H))),
                "min_Z": float(np.min(H_Z)),
            },
            "targeted_offgrid_PA8_passed": targeted_pass,
        },
        "complex_zero_exclusion": {
            "H_plusminus_i_sigma_root_distances": [float(v) for v in h_root_distances],
            "min_H_plusminus_i_sigma_root_distance": min_h_root_distance,
            "min_H_root_distance_over_tube": min_h_root_distance / tube,
            "L_root_distance": float(L_root_distance),
            "Pi_denominator_root_distance": float(pi_root_distance),
            "independent_zero_exclusion_passed": complex_pass,
        },
        "g_coefficient_norm": {
            "rho_over_cauchy": x,
            "independent_closed_form": exact_g_series,
            "public_upper": public_g,
            "public_to_independent_ratio": public_g / exact_g_series,
            "partial_sums": partials,
            "independent_norm_domination_passed": g_pass,
        },
        "mutation": mutation,
        "failed_guards": failed,
        "source_axis_domain_independent_preflight_passed": not failed,
        "truth_boundary": {
            "source_axis_domain_independently_audited": not failed,
            "source_rho_machine_bound_independently_audited": g_pass and complex_pass,
            "source_g_coefficient_norm_machine_bound_independently_audited": g_pass,
            "source_Phi_radius_one_ball_norm_machine_bound": False,
            "source_Phi_radius_one_ball_lipschitz_machine_bound": False,
            "source_R1_R2_machine_bound": False,
            "source_operator_M_K_machine_bound": False,
            "global_pressure_matched": False,
            "global_leading_profile_reconstructed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    identity = dict(report)
    report["receipt_sha256"] = hashlib.sha256(
        _canonical_json(identity).encode("utf-8")
    ).hexdigest()
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = run_audit()
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
