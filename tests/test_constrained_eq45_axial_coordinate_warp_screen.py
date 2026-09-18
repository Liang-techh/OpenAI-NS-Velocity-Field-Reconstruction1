import numpy as np

from openai_ns_reconstruction.constrained_eq45_axial_coordinate_warp_screen import (
    PHYSICAL_Z_CUT,
    _warp_z_and_jacobian,
    audit_axial_coordinate_warp_screen,
)


def test_axial_warp_is_identity_at_zero_and_support_boundary():
    z = np.linspace(-PHYSICAL_Z_CUT, PHYSICAL_Z_CUT, 4001)
    mapped0, jac0 = _warp_z_and_jacobian(z, 0.0)
    np.testing.assert_allclose(mapped0, z, atol=0.0, rtol=0.0)
    np.testing.assert_allclose(jac0, 1.0, atol=0.0, rtol=0.0)

    mapped, jac = _warp_z_and_jacobian(z, 0.5)
    assert mapped[0] == -PHYSICAL_Z_CUT
    assert mapped[-1] == PHYSICAL_Z_CUT
    assert np.all(np.diff(mapped) > 0.0)
    assert float(np.min(jac)) >= 0.5 - 1e-14
    assert np.all(np.isfinite(jac))


def test_axial_warp_report_keeps_truth_boundary_and_energy_preflight():
    report = audit_axial_coordinate_warp_screen(
        betas=(0.0, 0.2),
        quadrature_orders=(16, 24),
        morphology_times=(0.25, 0.5, 0.75),
        morphology_grid_size=17,
    )
    assert report["task_id"] == "CR003-AXIAL-COORDINATE-WARP-SCREEN-056"
    assert report["screen_protocol"]["new_basis_shapes_added"] == 0
    assert report["screen_protocol"]["production_beta_selected"] is False
    assert report["truth_boundary"]["canonical_velocity_changed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["held_out_pde_residual"]["evaluated"] is False

    baseline, trial = report["rows"]
    assert baseline["beta"] == 0.0
    assert baseline["baseline_parent_replay_without_common_scale"] is True
    assert abs(baseline["common_velocity_scale"] - 1.0) < 1e-15
    assert trial["beta"] == 0.2
    assert trial["reference_energy_gate_pass"] is True
    assert trial["validation_energy_range_all_pass"] is True
    assert trial["structure_checks"]["outside_physical_support_max_abs_velocity"] == 0.0
    assert trial["structure_checks"]["map_monotone_on_support"] is True
    assert trial["structure_checks"]["representative_core_signs_all_pass"] is True
    assert trial["piola_divergence_identity_used"] is True
