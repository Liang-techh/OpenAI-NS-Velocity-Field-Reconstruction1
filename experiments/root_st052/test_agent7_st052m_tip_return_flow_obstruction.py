from __future__ import annotations

import math

import agent7_st052m_tip_return_flow_obstruction as audit


def test_exact_sign_polynomial_and_physical_root():
    c2, c1, c0 = audit.sign_polynomial_coefficients()
    assert str(c2) == "-17"
    assert str(c1) == "954/25"
    assert str(c0) == "-81/25"

    roots = audit.sign_roots()
    assert roots[0] < float(audit.AXIAL_Q_A)
    assert float(audit.AXIAL_Q_A) < roots[1] < float(audit.AXIAL_Q_B)
    q_star, z_star = audit.physical_sign_change()
    assert q_star == roots[1]
    assert math.isclose(q_star, 2.1563200050974434, rel_tol=0.0, abs_tol=2e-15)
    assert math.isclose(z_star, 1.4684413522839252, rel_tol=0.0, abs_tol=2e-15)


def test_frozen_implementation_has_inward_core_and_outer_return_flow():
    report = audit.evaluate()
    analytic = report["analytic_result"]
    assert analytic["single_compact_tip_channel_strictly_inward_everywhere_possible"] is False
    assert analytic["return_flow_required_by_compact_lobe_identity"] is True
    assert analytic["compact_lobe_integral_Cr_exact"] == 0.0

    checks = report["implementation_spot_checks"]
    assert all(v < 0.0 for v in checks["inward_radial_components"].values())
    assert all(v > 0.0 for v in checks["outward_return_radial_components"].values())
    assert checks["shoulder_plus_abs_z_correction_norm"] == 0.0
    assert checks["shoulder_minus_abs_z_correction_norm"] == 0.0


def test_fresh_holdout_seed_is_local_inward_evidence_not_full_support_claim():
    report = audit.evaluate()
    holdout = report["fresh_holdout_scope"]
    assert holdout["tip_seed_abs_z"] == 1.25
    assert holdout["tip_seed_inside_analytic_inward_subband"] is True
    assert holdout["full_tip_support_inwardness_inferred_from_holdout"] is False


def test_truth_boundary_and_no_basis_growth():
    report = audit.evaluate()
    assert report["new_spatial_basis_added"] is False
    assert report["new_temporal_basis_added"] is False
    assert report["source_701_candidate_retuned"] is False
    assert report["trajectory_result_used"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["decision"]["generic_swirl_or_time_basis_growth_justified"] is False
    assert report["decision"]["velocity_changed_this_increment"] is False
    assert report["decision"]["closer_visualization_delivery_this_increment"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False
