from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_phi0_profile_contract import (
    DEFAULT_SERIES_TERMS,
    SCHEMA,
    KokunoPA10Phi0ProfileContract,
)


@pytest.fixture(scope="module")
def profiles() -> KokunoPA10Phi0ProfileContract:
    return KokunoPA10Phi0ProfileContract()


def test_phi0_matches_source_series_shape_and_is_vectorized(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    Y = np.asarray([[0.0], [0.7], [4.1]])
    eta = np.asarray([[-1.0, -0.3, 0.0, 0.4, 1.0]])
    values = profiles.values(Y, eta)
    axis = profiles.axis_profiles.values(Y, eta)

    assert values["Phi_0"].shape == (3, 5)
    np.testing.assert_array_equal(values["chi"], axis["chi"])
    np.testing.assert_array_equal(values["u_0"], axis["u_0"])
    np.testing.assert_allclose(
        values["f0_argument"], values["Y"] * values["chi"], rtol=0.0, atol=0.0
    )
    assert np.all(values["f0_tail_upper"] >= 0.0)
    assert float(np.max(values["f0_tail_upper"])) < 1.0e-140


def test_phi0_axis_is_regular_and_source_initial_coefficients_are_exact(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    eta = np.linspace(-1.0, 1.0, 31)
    values = profiles.values(0.0, eta)
    deriv = profiles.derivatives(0.0, eta)

    np.testing.assert_array_equal(values["Phi_0"], np.ones_like(eta))
    np.testing.assert_array_equal(values["f0_argument"], np.zeros_like(eta))
    np.testing.assert_array_equal(values["f0_prime"], np.full_like(eta, -0.25))
    np.testing.assert_allclose(values["f0_second"], np.full_like(eta, 1.0 / 24.0))
    np.testing.assert_allclose(deriv["Phi_0_Y"], -0.25 * values["chi"])
    np.testing.assert_allclose(deriv["Phi_0_eta"], np.zeros_like(eta), atol=0.0)
    np.testing.assert_allclose(
        deriv["Phi_0_YY"], values["chi"] ** 2 / 24.0, rtol=2.0e-15, atol=0.0
    )


def test_source_f0_ode_recurrence_holds_on_full_argument_range(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    z = np.linspace(0.0, 4.1, 257)
    f0, fp, fpp = profiles.f0_triplet(z)
    residual = 2.0 * (z * fpp + 2.0 * fp) + f0
    assert float(np.max(np.abs(residual))) < 2.0e-14
    assert np.all(np.isfinite(f0))
    assert np.any(np.abs(f0 - 1.0) > 0.1)


def test_analytic_Y_eta_derivatives_match_independent_centered_differences(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    Y = 1.37
    eta = 0.23
    h = 2.0e-6
    analytic = profiles.derivatives(Y, eta)

    y_plus = profiles.values(Y + h, eta)["Phi_0"]
    y_minus = profiles.values(Y - h, eta)["Phi_0"]
    eta_plus = profiles.values(Y, eta + h)["Phi_0"]
    eta_minus = profiles.values(Y, eta - h)["Phi_0"]
    dy = (y_plus - y_minus) / (2.0 * h)
    deta = (eta_plus - eta_minus) / (2.0 * h)

    np.testing.assert_allclose(analytic["Phi_0_Y"], dy, rtol=2.0e-8, atol=2.0e-10)
    np.testing.assert_allclose(analytic["Phi_0_eta"], deta, rtol=3.0e-7, atol=2.0e-9)

    yp_ep = profiles.values(Y + h, eta + h)["Phi_0"]
    yp_em = profiles.values(Y + h, eta - h)["Phi_0"]
    ym_ep = profiles.values(Y - h, eta + h)["Phi_0"]
    ym_em = profiles.values(Y - h, eta - h)["Phi_0"]
    mixed = (yp_ep - yp_em - ym_ep + ym_em) / (4.0 * h * h)
    np.testing.assert_allclose(
        analytic["Phi_0_Y_eta"], mixed, rtol=2.0e-4, atol=2.0e-6
    )


def test_configuration_roundtrip_is_deterministic(
    profiles: KokunoPA10Phi0ProfileContract, tmp_path
) -> None:
    path = tmp_path / "phi0_profile_config.json"
    saved = profiles.save_configuration(path)
    assert saved["schema"] == SCHEMA
    assert saved["series_terms"] == DEFAULT_SERIES_TERMS
    assert json.loads(path.read_text()) == saved

    restored = KokunoPA10Phi0ProfileContract.load_configuration(path)
    assert restored.configuration() == profiles.configuration()
    before = profiles.values(np.asarray([0.0, 1.1, 4.1]), np.asarray([0.0, 0.2, -0.7]))
    after = restored.values(np.asarray([0.0, 1.1, 4.1]), np.asarray([0.0, 0.2, -0.7]))
    for name in before:
        np.testing.assert_array_equal(before[name], after[name])


def test_fail_closed_on_invalid_series_or_source_argument(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    with pytest.raises(ValueError):
        KokunoPA10Phi0ProfileContract(series_terms=15)
    with pytest.raises(TypeError):
        KokunoPA10Phi0ProfileContract(series_terms=True)
    with pytest.raises(ValueError):
        profiles.f0_triplet(-1.0e-8)


def test_report_keeps_coordinate_and_science_boundary_fail_closed(
    profiles: KokunoPA10Phi0ProfileContract,
) -> None:
    report = profiles.report()
    assert report == profiles.report()
    assert report["receipt_sha256"] == profiles.sha256
    assert report["source"]["corrected_release_date"] == "2026-09-09"
    assert report["coordinate_contract"]["native_similarity_coordinates"] == ["Y", "eta"]
    assert report["coordinate_contract"]["repository_X_identification_asserted"] is False
    assert report["coordinate_contract"]["physical_F_from_Phi_mapping_resolved"] is False

    checks = report["machine_checks"]
    assert checks["Phi0_axis_value_exact_one_on_probe"] is True
    assert checks["Phi0_nontrivial_off_axis_on_probe"] is True
    assert checks["f0_source_ode_max_abs_error"] < 2.0e-14
    assert checks["all_profile_values_finite"] is True
    assert report["numerical_realization"]["maximum_probe_tail_upper"] < 1.0e-140

    truth = report["truth_boundary"]
    assert truth["source_native_phi0_contraction_center_executable"] is True
    assert truth["source_native_phi0_analytic_derivatives_executable"] is True
    assert truth["source_phi0_is_contraction_center_not_final_corrected_Phi"] is True
    assert truth["repository_X_identified_with_source_Y"] is False
    assert truth["physical_F_from_Phi_mapping_resolved"] is False
    assert truth["fixed_point_correction_materialized"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["unified_cartesian_velocity_export_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
