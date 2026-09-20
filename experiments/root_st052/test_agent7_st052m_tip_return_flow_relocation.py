from __future__ import annotations

import math

import agent7_st052m_tip_return_flow_relocation as audit


def test_preregistered_minimal_integer_rule_selects_p9():
    p, q_star, z_star = audit.select_inner_exponent()
    assert p == 9
    assert q_star > 1.0
    assert q_star < 3.24
    assert z_star >= 1.60
    _, previous_z = audit.physical_sign_change(8)
    assert previous_z < 1.60


def test_p9_sign_polynomial_and_physical_root_are_frozen():
    coeffs = audit.sign_polynomial_coefficients(9, 4)
    assert [str(v) for v in coeffs] == ["-27", "1764/25", "-81/25"]
    roots = audit.sign_roots(9, 4)
    assert math.isclose(roots[0], 0.04675485317423917, rel_tol=0.0, abs_tol=2e-15)
    assert math.isclose(roots[1], 2.566578480159094, rel_tol=0.0, abs_tol=2e-15)
    q_star, z_star = audit.physical_sign_change(9, 4)
    assert q_star == roots[1]
    assert math.isclose(z_star, 1.602054456052944, rel_tol=0.0, abs_tol=2e-15)


def test_same_dimension_reshape_moves_return_flow_outward_without_removing_it():
    report = audit.evaluate()
    selected = report["selected_envelope"]
    comparison = report["baseline_comparison"]
    scope = report["structural_scope"]

    assert selected["inner_exponent_p"] == 9
    assert selected["outer_exponent_m"] == 4
    assert selected["physical_abs_z_sign_change"] >= 1.60
    assert comparison["reshaped_return_collar_thickness"] < comparison["baseline_return_collar_thickness"]
    assert comparison["return_collar_reduction_fraction"] > 0.40
    assert scope["compact_lobe_integral_Cr_exact"] == 0.0
    assert scope["return_flow_removed"] is False
    assert scope["return_flow_relocated_outward"] is True
    assert scope["full_total_child_tip_inwardness_established"] is False


def test_preregistered_spot_signs_and_support_isolation():
    report = audit.evaluate()
    checks = report["implementation_spot_checks"]
    assert all(v < 0.0 for v in checks["inward_radial_components"].values())
    assert all(v > 0.0 for v in checks["outward_return_radial_components"].values())
    assert all(v == 0.0 for v in checks["shoulder_correction_norms"].values())
    assert checks["outside_support_correction_max_abs"] == 0.0


def test_truth_boundary_keeps_this_as_representation_preflight_only():
    report = audit.evaluate()
    assert report["prereg_issue"] == 731
    assert report["source_parent"] == {
        "pr": 723,
        "head": "1b16147b93e7b5279c833e8629e10aff1ef4f5a6",
    }
    assert report["basis_dimension_changed"] is False
    assert report["new_spatial_basis_direction_added"] is False
    assert report["new_temporal_basis_added"] is False
    assert report["nonlinear_child_evaluated"] is False
    assert report["coefficient_alpha_fitted"] is False
    assert report["trajectory_result_used"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["decision"]["same_dimension_envelope_reshape_justified_for_one_nonlinear_replay"] is True
    assert report["decision"]["second_poloidal_channel_justified_now"] is False
    assert report["decision"]["velocity_changed_this_increment"] is False
    assert report["decision"]["closer_visualization_delivery_this_increment"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False
