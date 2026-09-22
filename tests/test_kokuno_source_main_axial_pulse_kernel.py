import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_main_axial_pulse_kernel import (
    MAIN_XI_MAX,
    PARENT_EXACT_HEAD,
    SOURCE_COMMIT,
    KokunoSourceMainAxialPulseKernel,
    main_kernel_R0,
    main_kernel_R0_prime,
    phi,
    principal_Kb,
    principal_amplitude,
    smooth_sigma,
    smooth_sigma_prime,
)


def test_public_cutoff_phi_and_principal_amplitude_contract():
    z = np.array([-1.0, 0.0, 0.5, 1.0, 2.0])
    np.testing.assert_array_equal(smooth_sigma(z), np.array([0.0, 0.0, 0.5, 1.0, 1.0]))
    np.testing.assert_array_equal(smooth_sigma_prime(np.array([-1.0, 0.0, 1.0, 2.0])), 0.0)

    x = np.array([0.02, 0.5, 10.0])
    np.testing.assert_allclose(phi(x), x - 0.01, rtol=0.0, atol=2e-15)
    assert main_kernel_R0(0.0) == pytest.approx(0.0, abs=0.0)
    assert main_kernel_R0(MAIN_XI_MAX) == pytest.approx(0.0, abs=0.0)
    assert float(main_kernel_R0(5.0)) > 0.0

    kb = principal_Kb()
    amp = principal_amplitude()
    assert 0.2 < kb < 0.25
    assert 0.9 < amp < 1.2
    assert 4.0 * kb * amp * amp == pytest.approx(1.0 - math.exp(-26.0), rel=2e-15, abs=2e-15)


def test_r0_analytic_derivative_replays_centered_difference():
    for xi in (0.01, 2.0, 10.35):
        h = 2e-6
        fd = (float(main_kernel_R0(xi + h)) - float(main_kernel_R0(xi - h))) / (2.0 * h)
        analytic = float(main_kernel_R0_prime(xi))
        assert analytic == pytest.approx(fd, rel=2e-6, abs=2e-7)


def test_vectorized_profile_has_exact_entry_seam_and_nontrivial_axial_pulse():
    kernel = KokunoSourceMainAxialPulseKernel()
    lam = 1.0
    xi = np.array([0.0, 0.5, 5.0, 10.5, 11.0])
    X = np.exp(xi / lam)
    eta = np.linspace(-0.2, 0.2, xi.size)
    E0 = 1.0 / (1.0 + eta * eta)
    values = kernel.profile_values(X, eta, X_p=1.0, E_entry=E0, lambda_value=lam)

    assert values["E_source_main_pulse_principal"][0] == pytest.approx(E0[0], rel=0.0, abs=0.0)
    assert values["U_source_main_pulse_principal"][0] == pytest.approx(0.0, abs=0.0)
    assert values["U_source_main_pulse_principal"][-1] == pytest.approx(0.0, abs=0.0)
    assert np.max(np.abs(values["U_source_main_pulse_principal"][1:-1])) > 0.0
    np.testing.assert_allclose(
        values["F_source_main_pulse_principal"],
        values["E_source_main_pulse_principal"] / np.sqrt(2.0 * X),
        rtol=2e-15,
        atol=0.0,
    )


def test_profile_radial_derivatives_are_analytic_not_training_loss():
    kernel = KokunoSourceMainAxialPulseKernel()
    lam = 1.0
    eta = 0.17
    E0 = 1.0 / (1.0 + eta * eta)
    for xi in (2.0, 10.35):
        X = math.exp(xi)
        h = 2e-6 * X
        deriv = kernel.radial_derivatives(X, eta, X_p=1.0, E_entry=E0, lambda_value=lam)
        plus = kernel.profile_values(X + h, eta, X_p=1.0, E_entry=E0, lambda_value=lam)
        minus = kernel.profile_values(X - h, eta, X_p=1.0, E_entry=E0, lambda_value=lam)
        for value_key, derivative_key in (
            ("E_source_main_pulse_principal", "E_X_source_main_pulse_principal"),
            ("F_source_main_pulse_principal", "F_X_source_main_pulse_principal"),
            ("U_source_main_pulse_principal", "U_X_source_main_pulse_principal"),
        ):
            fd = (float(plus[value_key]) - float(minus[value_key])) / (2.0 * h)
            assert float(deriv[derivative_key]) == pytest.approx(fd, rel=3e-6, abs=1e-13)


def test_eta_derivative_consumes_only_entry_jet_for_autonomous_principal_amp():
    kernel = KokunoSourceMainAxialPulseKernel()
    eta = np.array([-0.2, 0.0, 0.25])
    E0 = 1.0 / (1.0 + eta * eta)
    E0_eta = -2.0 * eta / (1.0 + eta * eta) ** 2
    X = np.exp(np.array([1.0, 3.0, 8.0]))
    out = kernel.eta_derivatives(
        X,
        eta,
        X_p=1.0,
        E_entry=E0,
        E_entry_eta=E0_eta,
        lambda_value=1.0,
    )
    ratio = E0_eta / E0
    np.testing.assert_allclose(
        out["E_eta_source_main_pulse_principal"],
        ratio * out["E_source_main_pulse_principal"],
        rtol=3e-15,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        out["U_eta_source_main_pulse_principal"],
        ratio * out["U_source_main_pulse_principal"],
        rtol=3e-15,
        atol=1e-15,
    )


def test_domain_serialization_and_truth_boundary_fail_closed():
    kernel = KokunoSourceMainAxialPulseKernel()
    with pytest.raises(ValueError):
        kernel.profile_values(math.exp(11.01), 0.0, X_p=1.0, E_entry=1.0, lambda_value=1.0)
    with pytest.raises(ValueError):
        kernel.profile_values(0.99, 0.0, X_p=1.0, E_entry=1.0, lambda_value=1.0)

    config = kernel.to_configuration()
    restored = KokunoSourceMainAxialPulseKernel.from_configuration(config)
    assert restored.semantic_sha256 == kernel.semantic_sha256
    assert config["parent_exact_head"] == PARENT_EXACT_HEAD
    assert config["source"]["commit"] == SOURCE_COMMIT

    mutated = copy.deepcopy(config)
    mutated["truth_boundary"]["source_exact_amplitude_root_materialized"] = True
    with pytest.raises(ValueError):
        KokunoSourceMainAxialPulseKernel.from_configuration(mutated)

    truth = kernel.truth_boundary
    assert truth["source_main_pulse_kernel_formula_materialized"] is True
    assert truth["repository_autonomous_principal_amplitude_materialized"] is True
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["source_pulse_end_MJ_corrections_materialized"] is False
    assert truth["current_cartesian_pulse_velocity_composed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
