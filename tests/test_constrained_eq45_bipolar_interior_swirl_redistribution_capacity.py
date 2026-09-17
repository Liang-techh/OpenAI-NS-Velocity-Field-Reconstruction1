import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_interior_swirl_redistribution_capacity import (
    ORTHOGONALIZER,
    REDISTRIBUTION_MODE,
    TRUTH_BOUNDARY,
    audit_interior_swirl_redistribution_capacity,
)


def test_interior_swirl_redistribution_is_independent_localized_and_truth_bounded():
    report = audit_interior_swirl_redistribution_capacity()

    assert report["normalized_response_rank"] == 4
    assert np.isfinite(report["normalized_response_condition_number"])
    assert report["normalized_response_condition_number"] < 1.0e4
    assert report["finite_difference_refinement_relative_change"] < 1.0e-9
    assert report["redistribution_novelty_outside_F10_F20_compact_span"] > 1.0e-3

    overall = report["overall_response"][REDISTRIBUTION_MODE]
    assert overall["swirl_fraction"] > 1.0 - 1.0e-12
    assert overall["total_rms"] > 0.0

    regional = report["regional_response"]
    assert regional["core"][REDISTRIBUTION_MODE]["total_rms"] > 0.0
    assert regional["interior"][REDISTRIBUTION_MODE]["total_rms"] > 0.0
    assert regional["radial_flank"][REDISTRIBUTION_MODE]["total_rms"] < 1.0e-14
    assert regional["axial_flank"][REDISTRIBUTION_MODE]["total_rms"] < 1.0e-14
    assert report["redistribution_collar_response_energy_fraction"] < 1.0e-24

    basis = report["basis_definition"]
    assert basis["orthogonalizer"] == ORTHOGONALIZER == 6.5
    assert basis["analytic_parent_inner_product"] == 0.0
    assert basis["residual_fitted_shape_parameters"] == 0
    assert basis["public_image_fitted_shape_parameters"] == 0

    signs = report["sign_localization"]
    assert np.isclose(signs["analytic_sign_change_r_over_Rp"], np.sqrt(2.0 / 13.0))
    assert signs["inner_signed_swirl_mean"] > 0.0
    assert signs["shoulder_signed_swirl_mean"] < 0.0
    assert signs["inner_total_rms"] > 0.0
    assert signs["shoulder_total_rms"] > 0.0

    visual = report["visualization_fingerprint_leverage"]
    assert visual["coefficient_is_candidate_bound"] is False
    assert visual["max_relative_physical_plateau_collar_energy_change"] < 1.0e-12
    assert visual["plus_minus_radial_centroid_span"] > 1.0e-6
    assert visual["plus_minus_inner_energy_fraction_span"] > 1.0e-6

    assert report["parameter_growth"] == {
        "fixed_spatial_basis_shapes_added_if_materialized": 1,
        "scalar_coefficients_added_if_materialized": 1,
        "shape_parameters_fitted": 0,
        "temporal_degrees_added": 0,
    }
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["redistribution_coefficient_selected"] is False
    assert report["truth_boundary"]["perturbed_candidate_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
