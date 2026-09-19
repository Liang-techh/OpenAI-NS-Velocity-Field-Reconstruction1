import hashlib
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_source_axis_domain import (
    KokunoPA10SourceCompatibleAxisDomain,
)


def test_default_choice_certifies_source_PA8_on_enlarged_interval():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    cert = datum.real_interval_certificate()
    assert cert["source_PA8_implication_machine_certified"] is True
    assert cert["failed_intervals"] == 0
    assert cert["H_large_intervals"] > 0
    assert cert["Z_positive_exclusion_intervals"] > 0
    assert cert["required_H_abs_lower_on_K"] >= 10.0 * datum.sigma_star
    assert cert["PA8_chi_lower_from_dichotomy"] > 0.99
    assert cert["minimum_certified_H_abs_lower_on_H_branch"] > 10.0 * datum.sigma_star
    assert cert["minimum_certified_Z_lower_on_H_small_branch"] > 2.0 * datum.delta_star
    lo, hi = cert["enlarged_real_interval"]
    assert lo < -1.0
    assert hi > 1.0


def test_new_sigma_is_not_laundered_old_selected_sigma():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    assert datum.sigma_star == 1.0e-3
    assert datum.old_selected_sigma_star == 0.5
    assert datum.sigma_star != datum.old_selected_sigma_star
    truth = datum.truth_boundary
    assert truth["selected_source_compatible_sigma_PA8_machine_certified"] is True
    assert truth["old_selected_sigma_star_equals_source_compatible_choice"] is False
    assert truth["old_selected_local_pressure_audit_transfers_to_new_sigma"] is False
    assert truth["source_hidden_sigma_star_recovered"] is False


def test_complex_tube_cauchy_radius_and_g_norm_are_machine_bounded():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    cert = datum.complex_domain_certificate()
    assert cert["source_complex_neighborhood_machine_certified"] is True
    assert cert["H_plusminus_i_sigma_abs_lower"] > 0.0
    assert cert["H_square_plus_sigma_square_abs_lower"] > 0.0
    assert cert["L_abs_lower"] > 0.0
    assert cert["Pi_0_one_plus_z_squared_abs_lower"] > 0.0
    assert 0.0 < cert["coefficient_rho"] < cert["cauchy_radius"] < cert["complex_tube_radius"]
    assert cert["source_normalized_g_complex_sup_upper"] == 1.0
    assert 1.0 < cert["source_normalized_g_coefficient_norm_upper"] < 2.0
    assert cert["source_C_numeric_value_materialized"] is False


def test_pressure_operator_handoff_exposes_only_rho_and_g_not_missing_Phi():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    inputs = datum.pressure_operator_inputs()
    assert inputs["rho"] == datum.coefficient_rho
    assert math.isfinite(inputs["g_norm_upper"])
    assert inputs["g_norm_upper"] > 1.0
    assert inputs["Phi_radius_one_ball_norm"] is None
    assert inputs["Phi_radius_one_ball_lipschitz"] is None
    assert inputs["ready_for_pressure_operator_after_Phi_bound"] is True
    assert inputs["ready_for_source_pressure_ball_now"] is False


def test_axis_state_is_vectorized_nontrivial_and_finite():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    eta = np.asarray([-1.0, -0.25, 0.0, 0.25, 1.0])
    state = datum.axis_state(eta)
    for key in ("H_star", "Z_star", "Pi_0", "chi", "zeta_star"):
        assert state[key].shape == eta.shape
        assert np.all(np.isfinite(state[key]))
    assert np.max(np.abs(state["Z_star"])) > 1.0
    assert np.max(state["chi"]) > 0.99
    assert np.ptp(state["zeta_star"]) > 0.0
    with pytest.raises(ValueError):
        datum.axis_state(np.asarray([1.1]))


def test_truth_boundary_keeps_remaining_leading_and_PDE_gates_closed():
    truth = KokunoPA10SourceCompatibleAxisDomain().truth_boundary
    assert truth["selected_source_compatible_coefficient_rho_machine_bound"] is True
    assert truth["source_normalized_g_coefficient_norm_upper_executable"] is True
    assert truth["source_Phi_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_report_is_deterministic_hash_bound_and_saveable(tmp_path):
    datum = KokunoPA10SourceCompatibleAxisDomain()
    first = datum.report()
    second = datum.report()
    assert first == second
    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == expected
    path = tmp_path / "source_axis_domain.json"
    saved = datum.save_report(path)
    assert json.loads(path.read_text()) == saved


def test_invalid_domain_inputs_fail_closed():
    with pytest.raises(ValueError):
        KokunoPA10SourceCompatibleAxisDomain(sigma_star=0.5)
    with pytest.raises(ValueError):
        KokunoPA10SourceCompatibleAxisDomain(complex_tube_radius=1e-5, cauchy_radius=1e-5)
    with pytest.raises(ValueError):
        KokunoPA10SourceCompatibleAxisDomain(cauchy_radius=1e-5, coefficient_rho=1e-5)
    with pytest.raises(ValueError):
        KokunoPA10SourceCompatibleAxisDomain(partition_intervals=100)
