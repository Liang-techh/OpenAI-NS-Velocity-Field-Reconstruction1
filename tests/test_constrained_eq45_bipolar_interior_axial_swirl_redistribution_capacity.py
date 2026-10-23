import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_interior_axial_swirl_redistribution_capacity import (
    AXIAL_REDISTRIBUTION_MODE,
    TRUTH_BOUNDARY,
    audit_interior_axial_swirl_redistribution_capacity,
)


def test_axial_interior_swirl_redistribution_is_independent_localized_and_truth_bounded():
    report = audit_interior_axial_swirl_redistribution_capacity()

    assert report["normalized_response_rank"] == 4
    assert np.isfinite(report["normalized_response_condition_number"])
    assert report["normalized_response_condition_number"] < 1.0e3
    assert report["finite_difference_refinement_relative_change"] < 1.0e-9
    assert report["axial_redistribution_novelty_outside_F10_F20_compact_span"] > 1.0e-3

    overall = report["overall_response"][AXIAL_REDISTRIBUTION_MODE]
    assert overall["swirl_fraction"] > 1.0 - 1.0e-12
    assert overall["total_rms"] > 0.0

    signs = report["sign_localization"]
    assert signs["center_signed_swirl_mean"] > 0.0
    assert signs["axial_shoulder_signed_swirl_mean"] < 0.0
    assert signs["center_total_rms"] > 0.0
    assert signs["axial_shoulder_total_rms"] > 0.0

    regional = report["regional_response"]
    assert regional["center"][AXIAL_REDISTRIBUTION_MODE]["total_rms"] > 0.0
    assert regional["axial_shoulder"][AXIAL_REDISTRIBUTION_MODE]["total_rms"] > 0.0
    assert regional["radial_flank"][AXIAL_REDISTRIBUTION_MODE]["total_rms"] < 1.0e-14
    assert regional["axial_flank"][AXIAL_REDISTRIBUTION_MODE]["total_rms"] < 1.0e-14
    assert report["axial_redistribution_collar_response_energy_fraction"] < 1.0e-24

    visual = report["visualization_fingerprint_leverage"]
    assert visual["max_relative_physical_plateau_collar_energy_change"] < 1.0e-12
    assert visual["plus_minus_abs_z_centroid_span"] > 1.0e-6
    assert visual["plus_minus_central_energy_fraction_span"] > 1.0e-6

    assert report["parameter_growth"] == {
        "fixed_spatial_basis_shapes_added_if_materialized": 1,
        "scalar_coefficients_added_if_materialized": 1,
        "shape_parameters_fitted": 0,
        "temporal_degrees_added": 0,
    }
    assert report["basis_definition"]["residual_fitted_shape_parameters"] == 0
    assert report["basis_definition"]["public_image_fitted_shape_parameters"] == 0
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["perturbed_candidate_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
