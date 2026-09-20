from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_leading_axis_profile_contract import (
    SCHEMA,
    KokunoPA10LeadingAxisProfileContract,
)


@pytest.fixture(scope="module")
def profiles() -> KokunoPA10LeadingAxisProfileContract:
    return KokunoPA10LeadingAxisProfileContract()


def test_profile_contract_reuses_source_axis_state_and_is_vectorized(
    profiles: KokunoPA10LeadingAxisProfileContract,
) -> None:
    Y = np.asarray([[0.0], [1.0], [4.1]])
    eta = np.asarray([[-1.0, -0.25, 0.0, 0.5, 1.0]])
    values = profiles.values(Y, eta)
    state = profiles.domain.axis_state(np.broadcast_arrays(Y, eta)[1])

    assert values["U_star"].shape == (3, 5)
    for name in (
        "d", "L", "U_star", "H_star", "Z_star", "Pi_0", "Pi_0_eta",
        "chi", "zeta_star",
    ):
        np.testing.assert_array_equal(values[name], state[name])

    D = profiles.domain.D
    expected_W = (
        1.0
        - 4.0 * values["d"]
        - 2.0 * D * values["eta"] * values["U_star"]
    )
    np.testing.assert_allclose(values["W_star"], expected_W, rtol=0.0, atol=0.0)


def test_u0_center_is_axis_regular_nontrivial_and_radially_linear(
    profiles: KokunoPA10LeadingAxisProfileContract,
) -> None:
    eta = np.linspace(-1.0, 1.0, 17)
    axis = profiles.values(0.0, eta)
    off_axis = profiles.values(1.0, eta)
    farther = profiles.values(2.0, eta)

    np.testing.assert_array_equal(axis["u_0"], np.zeros_like(eta))
    assert np.any(np.abs(off_axis["u_0"]) > 0.0)
    np.testing.assert_allclose(farther["u_0"], 2.0 * off_axis["u_0"])
    np.testing.assert_allclose(off_axis["u_0"], off_axis["u_0_Y"])


def test_analytic_eta_derivatives_match_centered_differences(
    profiles: KokunoPA10LeadingAxisProfileContract,
) -> None:
    Y = 1.1
    eta = 0.2
    eps = 1.0e-6
    deriv = profiles.derivatives(Y, eta)
    plus = profiles.values(Y, eta + eps)
    minus = profiles.values(Y, eta - eps)

    pairs = {
        "U_star": "U_star_eta",
        "H_star": "H_star_eta",
        "W_star": "W_star_eta",
        "Z_star": "Z_star_eta",
        "chi": "chi_eta",
        "zeta_star": "zeta_star_eta",
        "u_0": "u_0_eta",
    }
    for value_name, derivative_name in pairs.items():
        numeric = (plus[value_name] - minus[value_name]) / (2.0 * eps)
        np.testing.assert_allclose(
            deriv[derivative_name], numeric, rtol=3.0e-5, atol=1.0e-7
        )

    np.testing.assert_allclose(deriv["u_0_Y"], plus["u_0_Y"], rtol=5.0e-5)


def test_configuration_roundtrip_preserves_values_and_source_choices(
    profiles: KokunoPA10LeadingAxisProfileContract, tmp_path
) -> None:
    path = tmp_path / "leading_axis_profiles.json"
    saved = profiles.save_configuration(path)
    assert saved["schema"] == SCHEMA
    assert json.loads(path.read_text()) == saved

    restored = KokunoPA10LeadingAxisProfileContract.load_configuration(path)
    assert restored.configuration() == profiles.configuration()
    first = profiles.values(np.asarray([0.0, 1.0]), np.asarray([0.0, 0.5]))
    second = restored.values(np.asarray([0.0, 1.0]), np.asarray([0.0, 0.5]))
    for name in first:
        np.testing.assert_array_equal(first[name], second[name])


def test_report_keeps_mapping_and_science_gates_fail_closed(
    profiles: KokunoPA10LeadingAxisProfileContract,
) -> None:
    report = profiles.report()
    assert report == profiles.report()
    assert report["receipt_sha256"] == profiles.sha256
    assert report["coordinate_contract"]["native_similarity_coordinates"] == ["Y", "eta"]
    assert report["coordinate_contract"]["repository_X_identification_asserted"] is False
    checks = report["machine_checks"]
    assert checks["axis_regular_u0_exact_zero_on_probe"] is True
    assert checks["u0_nontrivial_off_axis_on_probe"] is True
    assert checks["H_zeta_over_L_equals_minus_chi_abs_error"] < 1.0e-12
    assert checks["all_profile_values_finite"] is True
    assert checks["all_analytic_derivatives_finite"] is True

    truth = report["truth_boundary"]
    assert truth["source_native_leading_axis_profile_api_executable"] is True
    assert truth["source_native_profile_analytic_derivatives_executable"] is True
    assert truth["source_native_profile_configuration_serializable"] is True
    assert truth["repository_X_identified_with_source_Y"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["unified_cartesian_velocity_export_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
