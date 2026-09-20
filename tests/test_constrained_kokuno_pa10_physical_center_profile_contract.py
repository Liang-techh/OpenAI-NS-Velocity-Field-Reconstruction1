from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_physical_center_profile_contract import (
    KokunoPA10PhysicalCenterProfileContract,
)


def _relative_max(a: np.ndarray, b: np.ndarray) -> float:
    scale = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(np.abs(a - b) / scale))


def test_source_mapping_is_executable_vectorized_and_nontrivial() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    assert profile.Lambda == pytest.approx(10.0)
    assert profile.source_X_interval == pytest.approx((0.0, 0.41))

    X = np.array([[0.0], [0.07], [0.23], [0.41]])
    eta = np.array([[-0.75, -0.2, 0.35, 0.8]])
    values = profile.values(X, eta)
    assert values["F_0"].shape == (4, 4)
    assert np.array_equal(values["Y"], profile.Lambda * values["X"])
    assert np.all(np.isfinite(values["F_0"]))
    assert np.all(np.isfinite(values["U_0"]))
    assert np.all(np.isfinite(values["v_0"]))
    assert np.any(np.abs(values["F_0"]) > 1e-12)
    assert np.any(np.abs(values["U_0"]) > 1e-12)

    center = profile.center_profiles.values(values["Y"], values["eta"])
    axis = profile.axis_profiles.values(values["Y"], values["eta"])
    assert np.allclose(values["F_0"], values["g"] * center["Phi_0"], rtol=0, atol=0)
    assert np.allclose(
        values["U_0"], axis["U_star"] + center["u_0"] / profile.Lambda,
        rtol=0,
        atol=0,
    )
    assert np.all(values["E_0"][0] == 0.0)
    assert np.all(values["V_0"][0] == 0.0)


def test_phi_star_primitive_differentiates_to_public_zeta() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    eta = np.array([-0.82, -0.21, 0.17, 0.71])
    step = 2.0e-6
    fd = (
        profile.zeta_primitive(eta + step) - profile.zeta_primitive(eta - step)
    ) / (2.0 * step)
    public = profile.axis_profiles.values(np.zeros_like(eta), eta)["zeta_star"]
    assert _relative_max(fd, public) < 2.0e-8
    assert profile.zeta_primitive(np.array([0.0]))[0] == pytest.approx(0.0, abs=2e-15)
    assert np.all(profile.phi_star(eta) > 0.0)


def test_analytic_first_derivatives_match_centered_differences() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    X = np.array([0.05, 0.13, 0.29])
    eta = np.array([-0.6, 0.2, 0.72])
    deriv = profile.derivatives(X, eta)

    hx = 2.0e-6
    plus_x = profile.values(X + hx, eta)
    minus_x = profile.values(X - hx, eta)
    for key, derivative_key in (("F_0", "F_0_X"), ("U_0", "U_0_X"), ("v_0", "v_0_X")):
        fd = (plus_x[key] - minus_x[key]) / (2.0 * hx)
        assert _relative_max(fd, deriv[derivative_key]) < 3.0e-7

    he = 1.0e-6
    plus_e = profile.values(X, eta + he)
    minus_e = profile.values(X, eta - he)
    for key, derivative_key in (("F_0", "F_0_eta"), ("U_0", "U_0_eta")):
        fd = (plus_e[key] - minus_e[key]) / (2.0 * he)
        assert _relative_max(fd, deriv[derivative_key]) < 2.0e-6


def test_center_incompressibility_identity_and_axis_regular_V0() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    X = np.linspace(0.0, profile.source_X_interval[1], 23)
    eta = np.linspace(-0.95, 0.95, 21)
    XX, EE = np.meshgrid(X, eta, indexing="ij")
    values = profile.values(XX, EE)
    derivatives = profile.derivatives(XX, EE)
    axis_values = profile.axis_profiles.values(values["Y"], values["eta"])
    lhs = values["v_0"] + values["X"] * derivatives["v_0_X"]
    rhs = (
        2.0 * float(profile.axis_profiles.domain.A) * values["eta"] * values["U_0"]
        - axis_values["d"] * derivatives["U_0_eta"]
        + 2.0 * values["eta"] * values["X"] * derivatives["U_0_X"]
    ) / axis_values["L"]
    scale = np.maximum(1.0, np.maximum(np.abs(lhs), np.abs(rhs)))
    assert float(np.max(np.abs(lhs - rhs) / scale)) < 5.0e-15

    axis = profile.values(np.zeros_like(eta), eta)
    assert np.all(axis["V_0"] == 0.0)
    assert np.all(np.isfinite(axis["v_0"]))


def test_configuration_roundtrip_and_truth_boundary(tmp_path) -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    path = tmp_path / "physical_center.json"
    saved = profile.save_configuration(path)
    loaded = KokunoPA10PhysicalCenterProfileContract.load_configuration(path)
    assert loaded.configuration() == saved

    X = np.array([0.0, 0.1, 0.31])
    eta = np.array([-0.5, 0.0, 0.5])
    original = profile.values(X, eta)
    replay = loaded.values(X, eta)
    for key in ("F_0", "U_0", "v_0", "V_0", "phi_star"):
        assert np.array_equal(original[key], replay[key])

    truth = profile.truth_boundary
    assert truth["source_Y_equals_Lambda_X_mapping_resolved"] is True
    assert truth["physical_F_from_Phi_mapping_resolved"] is True
    assert truth["source_complex_C_normalization_certified"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["cartesian_spacetime_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_report_exposes_current_C_normalization_barrier_without_promoting_it() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    report = profile.report()
    checks = report["machine_checks"]
    assert checks["Y_equals_Lambda_X_exact_on_probe"] is True
    assert checks["axis_V0_exact_zero_on_probe"] is True
    assert checks["profiles_nontrivial"] is True
    assert checks["all_profile_values_finite"] is True
    assert np.isfinite(checks["incompressibility_scalar_identity_max_abs_defect"])
    assert checks["real_phi_star_probe_max"] > profile.C
    assert checks["configured_C_exceeds_real_probe_phi_star"] is False
    assert report["normalization_limitation"]["complex_domain_condition_certified_here"] is False
    assert report["normalization_limitation"]["real_probe_is_not_complex_domain_proof"] is True
    json.dumps(report, allow_nan=False)


def test_invalid_domain_inputs_fail_closed() -> None:
    profile = KokunoPA10PhysicalCenterProfileContract()
    with pytest.raises(ValueError):
        profile.values(-1e-12, 0.0)
    with pytest.raises(ValueError):
        profile.values(profile.source_X_interval[1] + 1e-9, 0.0)
    with pytest.raises(ValueError):
        profile.values(0.1, profile.eta_interval[1] + 1e-9)
    with pytest.raises(ValueError):
        profile.zeta_primitive(np.nan)
