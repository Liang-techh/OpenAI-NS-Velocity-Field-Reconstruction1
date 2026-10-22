import math

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_axial_envelope_flattening_capacity import (
    _flattened_envelope_jets,
    audit_axial_envelope_flattening_capacity,
)
from openai_ns_reconstruction.constrained_eq45_profile_basis import Eq45CompactProfileBasis


def test_alpha_zero_replays_governed_profile_jets():
    basis = Eq45CompactProfileBasis.seed()
    x = np.array([0.0, 0.4, 1.2, 2.8, 3.9])
    eta = np.array([0.0, 0.2, -0.45, 0.72, -0.95])
    expected = basis.evaluate(x, eta).stacked()
    actual = _flattened_envelope_jets(basis, x, eta, 0.0).stacked()
    assert np.max(np.abs(actual - expected)) < 5e-14


def test_axial_envelope_flattening_audit_preserves_support_and_reports_capacity():
    report = audit_axial_envelope_flattening_capacity(
        alpha_steps=(0.02, 0.01), trial_alpha=1.0, morphology_grid_size=33
    )

    growth = report["representation_increment"]
    assert growth["new_spatial_basis_shapes_added"] == 0
    assert growth["diagnostic_geometry_parameters_added"] == 1
    assert growth["new_production_fit_parameters_added"] == 0
    assert growth["production_parameter_selected"] is False

    structure = report["structure_checks"]
    assert structure["alpha_zero_parent_replay_max_abs_velocity_error"] < 1e-10
    assert structure["outside_physical_support_max_abs_velocity"] < 1e-12
    assert structure["cartesian_fd_divergence_max_abs"] < 2e-5
    assert structure["representative_core_signs_all_pass"] is True

    velocity = report["velocity_space_diagnostics"]
    assert velocity["rank_condition"]["normalized_rank"] == 5
    assert velocity["rank_condition"]["normalized_condition_number"] < 10.0
    assert velocity["finite_difference_refinement_relative_change"] < 1e-3
    assert velocity["envelope_novelty_outside_existing_four_span"] > 0.40
    assert velocity["envelope_response_rms_per_public_velocity_component"] > 0.0

    delta = report["trial_delta"]
    assert delta["axial_q90_over_Zp_delta"] > 0.07
    assert delta["axial_q99_over_Zp_delta"] > 0.07
    assert delta["axial_rms_over_Zp_delta"] > delta["radial_rms_over_Rp_delta"] > 0.0
    assert len(report["time_slice_morphology"]) == 3
    assert all(row["delta"]["axial_q90_over_Zp_delta"] > 0.05 for row in report["time_slice_morphology"])

    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False


def test_axial_envelope_flattening_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="alpha"):
        audit_axial_envelope_flattening_capacity(trial_alpha=1.2)
    with pytest.raises(ValueError, match="alpha_steps"):
        audit_axial_envelope_flattening_capacity(alpha_steps=(0.01, 0.02))
    with pytest.raises(ValueError, match="morphology_grid_size"):
        audit_axial_envelope_flattening_capacity(morphology_grid_size=24)
