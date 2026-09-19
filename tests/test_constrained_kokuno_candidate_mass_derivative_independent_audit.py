from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_candidate_mass_derivative_independent_audit import (
    FD6_STEPS,
    GUARDS,
    H_VALUES,
    PULSE_RESOLUTIONS,
    SEED,
    _geometry,
    _independent_mass_derivatives,
    _materialized_mass,
    run_audit,
)


def test_independent_mass_derivative_matches_black_box_value_difference() -> None:
    h = H_VALUES[1]
    geometry = _geometry(h)
    resolution = PULSE_RESOLUTIONS[-1]
    s_r, s_z = 0.031, -0.027
    D_r_h, D_z_h = _independent_mass_derivatives(geometry, resolution, s_r, s_z)
    step = 1.0e-5
    plus_r = np.asarray(
        _materialized_mass(geometry, h, resolution, s_r + step, s_z)["h_sigma_by_beta_sign"]
    )
    minus_r = np.asarray(
        _materialized_mass(geometry, h, resolution, s_r - step, s_z)["h_sigma_by_beta_sign"]
    )
    plus_z = np.asarray(
        _materialized_mass(geometry, h, resolution, s_r, s_z + step)["h_sigma_by_beta_sign"]
    )
    minus_z = np.asarray(
        _materialized_mass(geometry, h, resolution, s_r, s_z - step)["h_sigma_by_beta_sign"]
    )
    fd_r = (plus_r - minus_r) / (2.0 * step)
    fd_z = (plus_z - minus_z) / (2.0 * step)
    np.testing.assert_allclose(D_r_h, fd_r, rtol=2e-9, atol=2e-12)
    np.testing.assert_allclose(D_z_h, fd_z, rtol=2e-9, atol=2e-12)


def test_frozen_protocol_and_truth_boundary() -> None:
    report = run_audit()
    assert report["seed"] == SEED
    assert tuple(report["metrics"]["fd6_steps"]) == FD6_STEPS
    assert tuple(report["metrics"]["pulse_resolutions"]) == PULSE_RESOLUTIONS
    assert tuple(report["metrics"]["h_values"]) == H_VALUES
    assert report["guards"] == GUARDS
    assert report["independence"]["parent_directional_derivative_helper_used_as_oracle"] is False
    assert report["independence"]["training_tensor_or_loss_used"] is False
    assert report["independence"]["pressure_or_forcing_fit_used"] is False
    assert report["local_structural_preflight_passed"] is True
    assert all(report["checks"].values())
    truth = report["truth_boundary"]
    assert truth["candidate_mass_derivative_bridge_independently_audited"] is True
    assert truth["public_xyz_t_velocity_osc_provider_available"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
