from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_multiband_covariance_rank_screen import (
    KokunoMultiBandCovarianceRankScreen,
)


def _physical_family(*, duplicate: bool = False):
    phi = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    radial_scale = np.array([1.0, 0.8, 1.2], dtype=float)[:, None]
    w0 = np.zeros((3, 8, 3), dtype=float)
    w1 = np.zeros_like(w0)
    w0[..., 0] = radial_scale * np.cos(phi)
    w0[..., 1] = 0.7 * radial_scale * np.sin(phi)
    w0[..., 2] = 0.2 * radial_scale * np.cos(phi)
    if duplicate:
        w1 = 0.35 * w0
    else:
        w1[..., 0] = 0.6 * radial_scale * np.sin(phi)
        w1[..., 1] = -0.9 * radial_scale * np.cos(phi)
        w1[..., 2] = 0.4 * radial_scale * np.sin(phi)
    by_beta = np.stack((w0, w1), axis=-2)
    return {
        "beta_labels": ((5, (0, 0, 0)), (6, (1, 0, 0))),
        "active_ell_bands": (5, 6),
        "velocity_physical_cylindrical_by_beta": by_beta,
        "velocity_physical_cylindrical_total": np.sum(by_beta, axis=-2),
        "deterministic_permutation_invariant_aggregation": True,
    }


def test_distinct_band_tangents_give_local_rank_two_without_source_promotion():
    result = KokunoMultiBandCovarianceRankScreen().evaluate(
        _physical_family(), averaging_axes=(1,)
    )
    assert result["local_supplied_family_covariance_rank_two"] is True
    assert result["rank_two_cells"] == 3
    assert result["total_cells"] == 3
    assert result["rank_two_fraction"] == 1.0
    assert result["minimum_smallest_singular_value"] > 0.0
    assert result["minimum_singular_value_ratio"] > 1.0e-2
    assert result["minimum_band_response_novelty"] > 0.5
    assert result["repository_coefficient_unit_mapping_available"] is True
    assert result["source_coefficient_unit_mapping_available"] is False
    assert result["actual_source_mode_family_bound"] is False
    assert result["genuinely_independent_second_covariance_column_ready"] is False
    assert result["finite_correction_cycle_rerun_allowed"] is False
    assert result["full_ns_residual_assessed"] is False
    assert result["residual_reduction_claimed"] is False
    assert result["pde_validated"] is False


def test_duplicate_band_velocity_fails_covariance_rank_two():
    result = KokunoMultiBandCovarianceRankScreen().evaluate(
        _physical_family(duplicate=True), averaging_axes=(1,)
    )
    assert result["local_supplied_family_covariance_rank_two"] is False
    assert result["rank_two_cells"] == 0
    assert result["rank_two_fraction"] == 0.0
    assert result["minimum_smallest_singular_value"] < 1.0e-12
    assert result["genuinely_independent_second_covariance_column_ready"] is False


def test_deterministic_agent2_aggregation_is_required():
    data = _physical_family()
    data["deterministic_permutation_invariant_aggregation"] = False
    with pytest.raises(ValueError, match="deterministic aggregation"):
        KokunoMultiBandCovarianceRankScreen().evaluate(data, averaging_axes=(1,))


def test_averaging_axes_are_explicit_and_fail_closed():
    screen = KokunoMultiBandCovarianceRankScreen()
    with pytest.raises(ValueError, match="nonempty"):
        screen.evaluate(_physical_family(), averaging_axes=())
    with pytest.raises(ValueError, match="outside"):
        screen.evaluate(_physical_family(), averaging_axes=(2,))
    with pytest.raises(ValueError, match="duplicates"):
        screen.evaluate(_physical_family(), averaging_axes=(1, 1))
