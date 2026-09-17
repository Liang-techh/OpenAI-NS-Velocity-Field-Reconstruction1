import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_interior_swirl_2d_jacobian import (
    AXIAL_REDISTRIBUTION_MODE,
    RADIAL_REDISTRIBUTION_MODE,
    TRUTH_BOUNDARY,
    audit_interior_swirl_2d_jacobian,
)


def test_joint_interior_swirl_controls_are_independent_and_truth_bounded():
    report = audit_interior_swirl_2d_jacobian()

    velocity = report["public_velocity_response"]
    pair = velocity["radial_axial_pair"]
    assert velocity["normalized_rank"] == 5
    assert np.isfinite(velocity["normalized_condition_number"])
    assert velocity["normalized_condition_number"] < 1.0e3
    assert velocity["finite_difference_refinement_relative_change"] < 1.0e-9
    assert velocity["radial_novelty_outside_F10_F20_compact_span"] > 1.0e-3
    assert velocity["axial_novelty_outside_F10_F20_compact_radial_span"] > 1.0e-3
    assert pair["normalized_rank"] == 2
    assert np.isfinite(pair["normalized_condition_number"])
    assert pair["normalized_condition_number"] < 1.0e2
    assert abs(pair["cosine"]) < 0.99

    morphology = report["morphology_jacobian"]
    jacobian = np.asarray(morphology["fine_jacobian"])
    assert morphology["control_order"] == [RADIAL_REDISTRIBUTION_MODE, AXIAL_REDISTRIBUTION_MODE]
    assert jacobian.shape == (4, 2)
    assert morphology["normalized_rank"] == 2
    assert np.isfinite(morphology["normalized_condition_number"])
    assert morphology["normalized_condition_number"] < 1.0e3
    assert morphology["finite_difference_refinement_relative_change"] < 0.1
    assert morphology["centroid_rank"] == 2
    assert abs(jacobian[0, 0]) > 1.0e-5
    assert abs(jacobian[1, 1]) > 1.0e-7

    joint = report["joint_visualization_fingerprint"]
    assert joint["coefficient_is_candidate_bound"] is False
    assert joint["radial_centroid_span"] > 1.0e-4
    assert joint["abs_z_centroid_span"] > 1.0e-6
    assert joint["max_relative_physical_plateau_collar_energy_change"] < 1.0e-12

    contract = report["basis_contract"]
    assert contract["analytic_parent_radial_inner_product"] == 0.0
    assert contract["analytic_parent_axial_inner_product"] == 0.0
    assert contract["analytic_radial_axial_inner_product"] == 0.0
    assert contract["residual_fitted_shape_parameters"] == 0
    assert contract["public_image_fitted_shape_parameters"] == 0

    assert report["parameter_growth_this_increment"] == {
        "new_spatial_basis_shapes_added": 0,
        "new_scalar_coefficients_selected": 0,
        "shape_parameters_fitted": 0,
        "temporal_degrees_added": 0,
    }
    assert report["routing"]["rank_obstruction_to_two_mode_interior_swirl_control_observed"] is False
    assert report["routing"]["mixed_rz_swirl_mode_justified_by_capacity_rank_alone"] is False
    assert report["routing"]["stop_interior_swirl_basis_growth_on_capacity_ground"] is True

    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["new_basis_mode_added"] is False
    assert report["truth_boundary"]["perturbed_candidate_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
