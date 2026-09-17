import numpy as np
from openai_ns_reconstruction.constrained_eq45_bipolar_f20_capacity import TRUTH_BOUNDARY, audit_bipolar_f20_capacity


def test_bipolar_f20_adds_stable_outer_radial_swirl_control_only():
    report = audit_bipolar_f20_capacity()
    assert report["embedding_max_abs_velocity_error"] < 1e-12
    assert report["response_rank"] == 2
    assert np.isfinite(report["response_condition_number"])
    assert report["response_condition_number"] < 100.0
    assert report["finite_difference_refinement_relative_change"] < 1e-9
    assert report["F20_novelty_outside_F10_span"] > 0.10
    loc = report["localization_at_t_050"]
    assert loc["F10"]["swirl_fraction"] > 1.0 - 1e-10
    assert loc["F20"]["swirl_fraction"] > 1.0 - 1e-10
    assert report["F20_vs_F10_outer_selectivity_ratio"] > 1.5
    assert report["parameter_growth"] == {"spatial_modes_added": 1, "scalar_coefficients_added_if_selected": 1}
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
