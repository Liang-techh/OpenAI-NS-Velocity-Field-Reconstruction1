from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.kokuno_signed_physical_covariance_rank_screen import (
    KokunoSignedPhysicalCovarianceRankScreen,
)
from openai_ns_reconstruction.kokuno_source_signed_complete_curl_family import (
    KokunoSourceSignedCompleteCurlFamily,
)


def supplied_source_compatible_signed_family() -> tuple[dict, np.ndarray]:
    """Predeclared #482 supplied-data screen; no post-result phase tuning."""
    radial = np.linspace(0.78, 1.10, 5, dtype=float)
    quadrature_phase = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False, dtype=float)
    R, theta = np.meshgrid(radial, quadrature_phase, indexing="ij")

    alpha = 0.55 + 0.04 * (R - 0.90)
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.04 * tangent
    D_z_eta = np.zeros_like(eta)

    phase = np.empty(R.shape + (2, 2), dtype=float)
    phase[..., 0, 0] = quadrature_phase[None, :] + 0.17 + 0.20 * R
    phase[..., 0, 1] = quadrature_phase[None, :] - 0.23 + 0.20 * R
    phase[..., 1, 0] = quadrature_phase[None, :] + 0.31 - 0.10 * R
    phase[..., 1, 1] = quadrature_phase[None, :] + 0.07 - 0.10 * R

    n_phi = np.zeros(R.shape + (2, 2, 3), dtype=float)
    n_phi[...] = np.array([1.0, 0.0, 1.0])
    t_plus = np.empty(R.shape + (2, 2, 3), dtype=np.complex128)
    prototypes = np.array(
        [
            [[1.0 + 0.2j, 0.0 + 0.3j, -1.0 - 0.2j],
             [0.0 + 0.0j, 1.0 + 0.1j, 0.0 + 0.0j]],
            [[0.83 + 0.166j, 0.0 + 0.249j, -0.83 - 0.166j],
             [0.0 + 0.0j, 0.83 + 0.083j, 0.0 + 0.0j]],
        ],
        dtype=np.complex128,
    )
    t_plus[...] = prototypes
    zeros = np.zeros_like(t_plus)

    A_c = np.array([2.0, 2.2])
    u_star = np.array([3.0, 2.8])
    h_plus = np.array([1.1, 0.95])
    h_minus = np.array([0.9, 1.05])
    p = np.array([0.40, 0.36])
    q = np.array([0.08, -0.06])

    out = KokunoSourceSignedCompleteCurlFamily(h=0.005).physical_family(
        R=R,
        theta=theta,
        phase=phase,
        n_phi=n_phi,
        t_plus_prototype=t_plus,
        D_r_C_plus_prototype=zeros,
        D_z_C_plus_prototype=zeros,
        eta=eta,
        D_r_eta=D_r_eta,
        D_z_eta=D_z_eta,
        beta_labels=((5, "a"), (6, "b")),
        A_c=A_c,
        u_star=u_star,
        h_plus=h_plus,
        h_minus=h_minus,
        T_N=-A_c * p,
        T_K=u_star * q,
        D_r_A_c=np.array([0.03, -0.02]),
        D_r_u_star=np.array([0.04, 0.01]),
        D_r_h_plus=np.array([0.05, -0.03]),
        D_r_h_minus=np.array([-0.02, 0.04]),
        D_r_T_N=np.array([0.12, -0.09]),
        D_r_T_K=np.array([0.08, 0.06]),
        D_z_A_c=np.array([-0.02, 0.01]),
        D_z_u_star=np.array([0.03, -0.04]),
        D_z_h_plus=np.array([-0.04, 0.02]),
        D_z_h_minus=np.array([0.03, -0.01]),
        D_z_T_N=np.array([-0.10, 0.07]),
        D_z_T_K=np.array([0.06, -0.05]),
        direction_gap_eta=0.1,
    )
    return out, radial


def _raw_sign_velocity_rank(physical: dict) -> dict:
    columns = np.asarray(
        physical["velocity_physical_cartesian_by_beta_sign"], dtype=float
    )
    sign_columns = np.sum(columns, axis=-3)
    matrix = np.moveaxis(sign_columns, -2, 0).reshape(2, -1)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    ratio = float(singular_values[1] / singular_values[0]) if singular_values[0] > 0.0 else 0.0
    return {
        "raw_physical_sign_column_norms": np.linalg.norm(matrix, axis=1).tolist(),
        "raw_physical_sign_singular_values": singular_values.tolist(),
        "raw_physical_sign_smallest_singular_ratio": ratio,
        "raw_physical_sign_rank_two_at_1e8": bool(
            singular_values[0] > 0.0 and singular_values[1] > 1.0e-8 * singular_values[0]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    physical, radial = supplied_source_compatible_signed_family()
    raw_rank = _raw_sign_velocity_rank(physical)
    result = KokunoSignedPhysicalCovarianceRankScreen().evaluate(
        physical, averaging_axes=(1,)
    )
    payload = {
        "schema": "kokuno-a3-signed-physical-covariance-rank-v1",
        "agent2_parent": "PR #482 supplied-data signed complete-curl physical family",
        "radial_nodes": radial.tolist(),
        "phase_samples": 16,
        "beta_labels": [list(v) for v in result["beta_labels"]],
        "sigma_labels": list(result["sigma_labels"]),
        **raw_rank,
        "rank_two_cells": result["rank_two_cells"],
        "total_cells": result["total_cells"],
        "rank_two_fraction": result["rank_two_fraction"],
        "minimum_smallest_singular_value": result["minimum_smallest_singular_value"],
        "minimum_singular_value_ratio": result["minimum_singular_value_ratio"],
        "minimum_signed_response_novelty": result["minimum_signed_response_novelty"],
        "singular_values": np.asarray(result["singular_values"], dtype=float).tolist(),
        "covariance_response_jacobian": np.asarray(
            result["covariance_response_jacobian"], dtype=float
        ).tolist(),
        "maximum_total_reconstruction_abs_error": result[
            "maximum_total_reconstruction_abs_error"
        ],
        "maximum_sign_tangent_sum_abs_error": result[
            "maximum_sign_tangent_sum_abs_error"
        ],
        "maximum_quadratic_homogeneity_relative_error": result[
            "maximum_quadratic_homogeneity_relative_error"
        ],
        "cross_sign_and_cross_beta_terms_retained": True,
        "local_supplied_signed_physical_covariance_rank_two": result[
            "local_supplied_signed_physical_covariance_rank_two"
        ],
        "reference_covariance_rank_two": True,
        "physical_complete_curl_covariance_rank_two_assessed": True,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_signed_auxiliary_rectangles_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "source_actual_partition_labels_instantiated": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "real_candidate_defect_consumed": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_rerun_allowed": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "provenance_note": (
            "The sign-resolved complete-curl implementation is Agent-2 #482 and retains the "
            "corrected-reader sigma axis and Q/epsilon schedule. This deterministic receipt uses "
            "predeclared repository-supplied source-compatible phase/prototype/background inputs; "
            "it is not recovered Kokuno positive-order/background or actual signed rectangles."
        ),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
