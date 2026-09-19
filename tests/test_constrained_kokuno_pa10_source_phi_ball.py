import hashlib
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_source_axis_domain import (
    KokunoPA10SourceCompatibleAxisDomain,
)
from openai_ns_reconstruction.kokuno_pa10_source_phi_ball import (
    KokunoPA10SourcePhiBallBounds,
)


def test_chi_bound_uses_certified_tube_and_is_finite():
    bound = KokunoPA10SourcePhiBallBounds()
    cert = bound.chi_certificate()
    assert 0.0 < cert["relative_H_variation_e"] < 1.0
    assert cert["single_factor_ratio_upper"] > 1.0
    assert 1.0 < cert["chi_complex_sup_upper"] < 10.0
    assert 1.0 < cert["cauchy_weight_sum_upper"] < 2.0
    assert cert["chi_multiplier_norm_upper"] > cert["chi_complex_sup_upper"]
    assert cert["chi_multiplier_norm_upper"] < 20.0


def test_phi_center_and_radius_one_ball_are_machine_bounded():
    bound = KokunoPA10SourcePhiBallBounds()
    cert = bound.phi_ball_certificate()
    assert math.isfinite(cert["inverse_one_plus_T_absolute_series_upper"])
    assert cert["inverse_series_terms_before_tail_majorant"] > 0
    assert cert["Phi0_coefficient_norm_upper"] >= 1.0
    assert cert["Phi_radius_one_ball_norm_upper"] > cert["Phi0_coefficient_norm_upper"]
    assert cert["Phi_radius_one_ball_norm_upper"] <= math.nextafter(
        cert["Phi0_coefficient_norm_upper"] + 1.0, math.inf
    )
    assert cert["Phi_radius_one_ball_lipschitz_upper"] == 1.0


def test_pressure_handoff_is_no_longer_missing_Phi_inputs():
    bound = KokunoPA10SourcePhiBallBounds()
    inputs = bound.completed_pressure_operator_inputs()
    assert inputs["rho"] == 1.0e-6
    assert 1.0 < inputs["g_norm_upper"] < 2.0
    assert math.isfinite(inputs["Phi_radius_one_ball_norm"])
    assert inputs["Phi_radius_one_ball_norm"] > 1.0
    assert inputs["Phi_radius_one_ball_lipschitz"] == 1.0
    assert inputs["ready_for_source_compatible_pressure_operator"] is True
    assert inputs["mixed_Y_Phi_Y_still_required_for_R1"] is True
    assert inputs["mixed_Y_u_Y_still_required_for_R2"] is True


def test_existing_pressure_algebra_accepts_source_compatible_inputs():
    pressure = KokunoPA10SourcePhiBallBounds().pressure_ball_envelope()
    for name in ("integrand_g2_Phi2", "p", "p_eta", "Y_p_Y"):
        assert math.isfinite(pressure[name]["norm"])
        assert math.isfinite(pressure[name]["lipschitz"])
        assert pressure[name]["norm"] > 0.0
        assert pressure[name]["lipschitz"] > 0.0
    assert pressure["p_eta"]["norm"] > pressure["p"]["norm"]
    assert pressure["p_eta"]["lipschitz"] > pressure["p"]["lipschitz"]


def test_truth_boundary_promotes_only_phi_pressure_inputs_not_R1_R2_MK():
    truth = KokunoPA10SourcePhiBallBounds().truth_boundary
    assert truth["source_compatible_chi_multiplier_norm_machine_bound"] is True
    assert truth["source_compatible_Phi_radius_one_ball_norm_machine_bound"] is True
    assert truth["source_compatible_Phi_radius_one_ball_lipschitz_machine_bound"] is True
    assert truth["source_compatible_pressure_operator_inputs_complete"] is True
    assert truth["source_axis_domain_independent_agent4_admission_required"] is True
    assert truth["source_Phi_ball_independent_agent4_admission_required"] is True
    assert truth["source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_R2_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_report_is_deterministic_hash_bound_and_saveable(tmp_path):
    bound = KokunoPA10SourcePhiBallBounds()
    first = bound.report()
    second = bound.report()
    assert first == second
    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == expected
    path = tmp_path / "source_phi_ball.json"
    saved = bound.save_report(path)
    assert json.loads(path.read_text()) == saved


def test_rejects_wrong_domain_type():
    with pytest.raises(TypeError):
        KokunoPA10SourcePhiBallBounds(domain=object())


def test_source_compatible_domain_identity_is_preserved():
    domain = KokunoPA10SourceCompatibleAxisDomain()
    bound = KokunoPA10SourcePhiBallBounds(domain=domain)
    report = bound.report()
    assert report["source_compatible_axis_identity"]["axis_receipt_sha256"] == domain.sha256
    assert report["source_compatible_axis_identity"]["sigma_star"] == 1.0e-3
    assert report["integration_boundary"]["old_sigma_0p5_local_pressure_interface_reused"] is False
