from __future__ import annotations

import copy
import hashlib
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_chi_multiplier import (
    KokunoPA10SelectedChiMultiplier,
)


def test_selected_complex_tube_is_certified_fail_closed() -> None:
    certifier = KokunoPA10SelectedChiMultiplier()
    certificate = certifier.certificate()

    assert certificate["selected_axis_complex_domain_certified"] is True
    assert certificate["analytic_tube_radius"] == pytest.approx(5.0e-3)
    assert certificate["coefficient_rho"] == pytest.approx(5.0e-4)
    assert certificate["H_plusminus_i_sigma_factor_abs_lower"] > 0.0
    assert certificate["chi_denominator_abs_lower"] > 0.0
    assert certificate["L_abs_lower"] > 0.0
    assert certificate["one_plus_eta_squared_abs_lower"] > 0.0


def test_tube_bounds_dominate_representative_complex_samples() -> None:
    certifier = KokunoPA10SelectedChiMultiplier()
    certificate = certifier.certificate()
    R = certifier.analytic_tube_radius

    # This is a regression/sanity replay, not the proof: the proof itself is
    # the analytic triangle/Cauchy bound implemented in certificate().
    for real in (-1.0, -0.75, -0.25, 0.0, 0.25, 0.75, 1.0):
        for imag in (-R, -0.5 * R, 0.0, 0.5 * R, R):
            z = complex(real, imag)
            H = certifier.H_star_complex(z)
            chi = certifier.chi_complex(z)
            assert abs(H) <= certificate["H_tube_abs_upper"]
            assert abs(chi) <= certificate["chi_sup_abs_upper"]


def test_multiplier_bound_uses_source_weight_cauchy_sum() -> None:
    certifier = KokunoPA10SelectedChiMultiplier()
    certificate = certifier.certificate()
    x = certifier.coefficient_rho / certifier.analytic_tube_radius
    expected_series_factor = (1.0 + x) / ((1.0 - x) ** 3)
    assert certificate["coefficient_multiplier_series_factor"] == pytest.approx(
        expected_series_factor
    )
    assert certificate["selected_multiplier_norm_chi_upper"] == pytest.approx(
        certificate["chi_sup_abs_upper"] * expected_series_factor
    )
    assert math.isfinite(certifier.multiplier_norm_chi_upper)
    assert certifier.multiplier_norm_chi_upper > 1.0


def test_selected_multiplier_bound_feeds_569_operator_primitives() -> None:
    certifier = KokunoPA10SelectedChiMultiplier()
    receipt = certifier.operator_input_receipt()
    assert receipt["coefficient_rho"] == pytest.approx(certifier.coefficient_rho)
    assert receipt["selected_multiplier_norm_chi_upper"] == pytest.approx(
        certifier.multiplier_norm_chi_upper
    )
    assert receipt["mixed_derivative_template_constant"] > 0.0
    assert math.isfinite(receipt["inverse_one_plus_T_absolute_bound"])
    assert receipt["inverse_one_plus_T_absolute_bound"] > 1.0


def test_invalid_selected_radii_fail_closed() -> None:
    with pytest.raises(ValueError, match="analytic_tube_radius"):
        KokunoPA10SelectedChiMultiplier(analytic_tube_radius=0.0)
    with pytest.raises(ValueError, match="analytic_tube_radius"):
        KokunoPA10SelectedChiMultiplier(analytic_tube_radius=0.25)
    with pytest.raises(ValueError, match="coefficient_rho_fraction"):
        KokunoPA10SelectedChiMultiplier(coefficient_rho_fraction=0.0)
    with pytest.raises(ValueError, match="coefficient_rho_fraction"):
        KokunoPA10SelectedChiMultiplier(coefficient_rho_fraction=1.0)


def test_truth_boundary_does_not_promote_selected_bound_to_source_certificate() -> None:
    truth = KokunoPA10SelectedChiMultiplier().report()["truth_boundary"]
    assert truth["selected_axis_chi_complex_tube_machine_bound"] is True
    assert truth["selected_coefficient_rho_executable"] is True
    assert truth["selected_chi_multiplier_norm_machine_bound"] is True
    assert truth["source_sigma_star_admissibility_verified"] is False
    assert truth["source_enlarged_real_interval_I_recovered"] is False
    assert truth["source_rho_machine_bound"] is False
    assert truth["source_chi_multiplier_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path) -> None:
    certifier = KokunoPA10SelectedChiMultiplier()
    payload = certifier.to_payload()
    replay = KokunoPA10SelectedChiMultiplier.from_payload(payload)
    assert replay.sha256 == certifier.sha256
    assert replay.to_payload() == payload

    path = tmp_path / "selected_chi_multiplier.json"
    certifier.save_json(path)
    loaded = KokunoPA10SelectedChiMultiplier.load_json(path)
    assert loaded.sha256 == certifier.sha256
    assert json.loads(path.read_text())["sha256"] == certifier.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_chi_multiplier_norm_machine_bound"] = True
    body = copy.deepcopy(tampered)
    body.pop("sha256")
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    tampered["sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoPA10SelectedChiMultiplier.from_payload(tampered)
