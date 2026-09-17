import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_interior_compact_swirl_capacity import (
    COMPACT_MODE,
    TRUTH_BOUNDARY,
    audit_interior_compact_swirl_capacity,
)


def test_interior_compact_swirl_is_independent_localized_and_truth_bounded():
    report = audit_interior_compact_swirl_capacity()

    assert report["normalized_response_rank"] == 3
    assert np.isfinite(report["normalized_response_condition_number"])
    assert report["normalized_response_condition_number"] < 1.0e3
    assert report["finite_difference_refinement_relative_change"] < 1.0e-9
    assert report["compact_novelty_outside_F10_F20_span"] > 1.0e-3

    overall = report["overall_response"][COMPACT_MODE]
    assert overall["swirl_fraction"] > 1.0 - 1.0e-12
    assert overall["total_rms"] > 0.0

    regional = report["regional_response"]
    assert regional["core"][COMPACT_MODE]["total_rms"] > 0.0
    assert regional["interior"][COMPACT_MODE]["total_rms"] > 0.0
    assert regional["radial_flank"][COMPACT_MODE]["total_rms"] < 1.0e-14
    assert regional["axial_flank"][COMPACT_MODE]["total_rms"] < 1.0e-14
    assert report["compact_collar_response_energy_fraction"] < 1.0e-24

    visual = report["visualization_fingerprint_leverage"]
    assert visual["max_relative_physical_plateau_collar_energy_change"] < 1.0e-12
    assert visual["max_fractional_weighted_radial_rms_leverage"] > 1.0e-6

    assert report["parameter_growth"] == {
        "fixed_spatial_basis_shapes_added_if_materialized": 1,
        "scalar_coefficients_added_if_materialized": 1,
        "shape_parameters_fitted": 0,
        "temporal_degrees_added": 0,
    }
    assert report["basis_definition"]["residual_fitted_shape_parameters"] == 0
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["perturbed_candidate_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
