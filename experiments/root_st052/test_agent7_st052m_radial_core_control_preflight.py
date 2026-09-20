from __future__ import annotations

import math

import numpy as np

import agent7_st052m_radial_core_control_preflight as audit


def test_frozen_family_preserves_live_shape_and_support() -> None:
    audit._assert_source_lock()
    assert audit.K_VALUES == (0.50, 0.75, 1.00, 1.50, 2.00)
    assert audit.LIVE_K == 1.0
    q = np.linspace(0.0, audit.Q_OUTER, 101)
    live_value, live_deriv = audit.radial_profile_k(q, 1.0)
    source_value, source_deriv = audit.reservoir._radial_profile_and_derivative(q)
    assert np.array_equal(live_value, source_value)
    assert np.array_equal(live_deriv, source_deriv)
    for k in audit.K_VALUES:
        assert audit.support_check(k)["outside_or_boundary_max_abs_Cz_over_Z"] == 0.0


def test_root_formula_matches_live_response_and_is_monotone() -> None:
    roots = []
    for k in audit.K_VALUES:
        row = audit.root_geometry(k)
        assert row["analytic_numeric_q_mismatch"] <= audit.ROOT_TOL
        assert 0.0 < row["analytic_q_root"] < audit.Q_OUTER
        roots.append(row["r_star"])
    assert all(a > b for a, b in zip(roots, roots[1:]))
    expected_live = math.sqrt((13.0 + 5.0 * math.sqrt(17.0)) / 30.0)
    assert abs(audit.root_geometry(1.0)["r_star"] - expected_live) <= 5.0e-14


def test_shape_control_is_independent_and_well_conditioned() -> None:
    for k in audit.K_VALUES:
        row = audit.response_geometry(k)
        if k == 1.0:
            assert row["identical_to_live_shape"] is True
            assert row["two_column_rank_vs_live_k1"] == 1
            assert row["two_column_condition_vs_live_k1"] is None
        else:
            assert row["identical_to_live_shape"] is False
            assert row["two_column_rank_vs_live_k1"] == 2
            assert row["two_column_condition_vs_live_k1"] <= audit.CONDITION_MAX
            assert abs(row["response_cosine_vs_live_k1"]) < 1.0


def test_flux_identity_and_live_sensitivity() -> None:
    for k in audit.K_VALUES:
        row = audit.flux_refinement(k)
        assert row["passes"] is True
        assert row["finest_abs_signed_flux"] <= audit.FLUX_TOL
        assert row["medium_to_fine_change"] <= audit.REFINEMENT_TOL
    sensitivity = audit.sensitivity_at_live()
    assert sensitivity["analytic_dr_star_dk"] < 0.0
    assert sensitivity["relative_mismatch"] <= audit.SENSITIVITY_REL_TOL


def test_truth_boundary_and_report(tmp_path) -> None:
    report = audit.run(tmp_path / "report.json")
    assert report["task_id"] == audit.TASK_ID
    assert report["prereg_issue"] == 800
    assert report["source_parent"] == {"pr": 794, "head": audit.SOURCE_PARENT_HEAD}
    assert report["family"]["candidate_selection"] is None
    assert report["family"]["public_openai_core_radius_target"] is None
    assert report["capacity_conclusion"]["same_dimension_core_thickness_control_demonstrated"] is True
    assert report["capacity_conclusion"]["strict_monotone_core_radius_control_on_frozen_grid"] is True
    assert report["capacity_conclusion"]["second_poloidal_basis_justified_by_this_audit"] is False
    assert report["capacity_conclusion"]["candidate_k_selected"] is False
    assert report["capacity_conclusion"]["total_child_visual_improvement_established"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["parameter_fit_performed"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["public_image_numeric_target_used"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["openai_field_identified"] is False
