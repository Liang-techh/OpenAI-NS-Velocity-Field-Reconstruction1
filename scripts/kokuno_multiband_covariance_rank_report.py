from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.kokuno_multiband_covariance_rank_screen import (
    KokunoMultiBandCovarianceRankScreen,
)
from openai_ns_reconstruction.kokuno_source_multiband_real_pair_family import (
    KokunoSourceMultiBandRealPairFamily,
)


def supplied_source_compatible_family() -> dict:
    radial = np.linspace(0.82, 1.18, 5, dtype=float)
    angle = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False, dtype=float)
    R, theta = np.meshgrid(radial, angle, indexing="ij")

    alpha = 0.23 + 0.15 * R
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.15 * tangent
    D_z_eta = np.zeros_like(eta)

    phase = np.stack(
        (0.91 * R + 0.70 * theta + 0.13, 1.07 * R - 0.50 * theta - 0.27),
        axis=-1,
    )
    n_phi = np.zeros(R.shape + (2, 3), dtype=float)
    n_phi[..., 0] = 1.0

    t_plus = np.empty(R.shape + (2, 3), dtype=np.complex128)
    t_plus[..., 0, :] = np.array([0.0, 0.24 + 0.07j, -0.11 + 0.04j])
    t_plus[..., 1, :] = np.array([0.0, -0.18 + 0.05j, 0.16 - 0.03j])
    D_r_C_plus = np.zeros_like(t_plus)
    D_z_C_plus = np.zeros_like(t_plus)

    family = KokunoSourceMultiBandRealPairFamily(h=0.004)
    out = family.physical_family(
        R=R,
        theta=theta,
        phase=phase,
        n_phi=n_phi,
        t_plus=t_plus,
        D_r_C_plus=D_r_C_plus,
        D_z_C_plus=D_z_C_plus,
        eta=eta,
        D_r_eta=D_r_eta,
        D_z_eta=D_z_eta,
        beta_labels=((5, (0, 0, 0)), (6, (1, 0, 0))),
    )
    out["_report_radial_nodes"] = radial
    return out


def _velocity_band_rank(family: dict) -> dict:
    by_band = np.asarray(family["velocity_physical_cartesian_by_band"], dtype=float)
    if by_band.shape[-2] != 2 or by_band.shape[-1] != 3:
        raise RuntimeError("receipt expects exactly two Cartesian band columns")
    matrix = np.moveaxis(by_band, -2, 0).reshape(2, -1)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    ratio = float(singular_values[1] / singular_values[0]) if singular_values[0] > 0.0 else 0.0
    return {
        "band_column_norms": np.asarray(family["band_column_norms"], dtype=float).tolist(),
        "velocity_band_singular_values": singular_values.tolist(),
        "velocity_band_smallest_singular_ratio": ratio,
        "velocity_band_rank_two_at_1e-8": bool(
            singular_values[0] > 0.0 and singular_values[1] > 1.0e-8 * singular_values[0]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    family = supplied_source_compatible_family()
    radial = np.asarray(family.pop("_report_radial_nodes"), dtype=float)
    velocity_rank = _velocity_band_rank(family)
    result = KokunoMultiBandCovarianceRankScreen().evaluate(family, averaging_axes=(1,))
    payload = {
        "schema": "kokuno-a3-supplied-multiband-covariance-rank-v1",
        "radial_nodes": radial.tolist(),
        "active_ell_bands": list(result["active_ell_bands"]),
        **velocity_rank,
        "rank_two_cells": result["rank_two_cells"],
        "total_cells": result["total_cells"],
        "rank_two_fraction": result["rank_two_fraction"],
        "minimum_smallest_singular_value": result["minimum_smallest_singular_value"],
        "minimum_singular_value_ratio": result["minimum_singular_value_ratio"],
        "minimum_band_response_novelty": result["minimum_band_response_novelty"],
        "singular_values": np.asarray(result["singular_values"], dtype=float).tolist(),
        "local_supplied_family_covariance_rank_two": result[
            "local_supplied_family_covariance_rank_two"
        ],
        "actual_source_mode_family_bound": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "finite_correction_cycle_rerun_allowed": False,
        "full_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "provenance_note": (
            "Agent-2 #461 source-scheduled supplied-mode interface; phase/pulse data in this receipt "
            "are repository-manufactured source-compatible inputs, not recovered Kokuno modes"
        ),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
