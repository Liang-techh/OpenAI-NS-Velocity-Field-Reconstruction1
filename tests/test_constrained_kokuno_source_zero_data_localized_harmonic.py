from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_sensitivity import (
    SourceModeForcingDirectionalJet,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizationJet,
    materialize_source_zero_data_localized_harmonic,
    source_zero_data_localized_harmonic_contract,
)


EPSILON = 0.25
P = 0.0
P_Z = 0.0
X_0 = 1.0
K = 2.0
M = 1
PULSE_V = 0.30
THETA = 0.41


def _zero_background() -> SourceBackgroundPhaseJet:
    return SourceBackgroundPhaseJet(
        f=0.0,
        g=0.0,
        f_r=0.0,
        g_r=0.0,
        f_z=0.0,
        g_z=0.0,
        f_rr=0.0,
        g_rr=0.0,
        f_rz=0.0,
        g_rz=0.0,
        f_zz=0.0,
        g_zz=0.0,
    )


def _forcing(_: float) -> SourceModeForcingDirectionalJet:
    value = np.asarray((0.0 + 0.0j, 0.7 + 0.2j, -0.35 + 0.1j))
    zero = np.zeros(3, dtype=np.complex128)
    return SourceModeForcingDirectionalJet(value=value, dr_value=zero, dz_value=zero)


def _eta_data(radius, z_normalized):
    radius = np.asarray(radius, dtype=float)
    z_normalized = np.asarray(z_normalized, dtype=float)
    eta = np.exp(-0.7 * (radius - 1.05) ** 2 - 0.4 * z_normalized**2)
    dr_eta = -1.4 * (radius - 1.05) * eta
    # D_z = epsilon partial_Z on the auxiliary-independent specialization.
    dz_eta = EPSILON * (-0.8 * z_normalized) * eta
    return eta, dr_eta, dz_eta


def _run(radius=1.22, z_normalized=0.23, *, pulse_cutoff=0.8, eta_override=None):
    if eta_override is None:
        eta, dr_eta, dz_eta = _eta_data(radius, z_normalized)
    else:
        eta, dr_eta, dz_eta = eta_override
    return materialize_source_zero_data_localized_harmonic(
        PULSE_V,
        radius,
        THETA,
        z_normalized,
        _zero_background(),
        _forcing,
        SourceLocalizationJet(
            eta=eta,
            dr_eta=dr_eta,
            dz_eta=dz_eta,
            pulse_cutoff=pulse_cutoff,
        ),
        epsilon=EPSILON,
        p=P,
        p_z=P_Z,
        x_0=X_0,
        k=K,
        m=M,
    )


def test_zero_data_chain_materializes_nontrivial_localized_complete_curl():
    out = _run()
    assert out.amplitude_sensitivity.panels == (32, 64, 128)
    assert np.max(np.abs(out.amplitude_sensitivity.amplitude_jet.t_m)) > 1.0e-8
    assert np.max(np.abs(out.localized_coefficient)) > 1.0e-8
    assert np.max(np.abs(out.complex_velocity_cylindrical)) > 1.0e-8
    np.testing.assert_allclose(
        out.real_pair_velocity_cylindrical,
        2.0 * np.real(out.complex_velocity_cylindrical),
        rtol=0.0,
        atol=0.0,
    )


def test_product_rule_keeps_support_gradient_curl_terms():
    radius = 1.22
    z_normalized = 0.23
    eta, dr_eta, dz_eta = _eta_data(radius, z_normalized)
    localized = _run(radius, z_normalized, pulse_cutoff=0.8)
    base = _run(
        radius,
        z_normalized,
        pulse_cutoff=1.0,
        eta_override=(1.0, 0.0, 0.0),
    )
    naive = 0.8 * eta * base.complex_velocity_cylindrical
    assert np.linalg.norm(localized.complex_velocity_cylindrical - naive) > 1.0e-6
    np.testing.assert_allclose(
        localized.dr_localization_weight,
        0.8 * dr_eta,
        rtol=0.0,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(
        localized.dz_localization_weight,
        0.8 * dz_eta,
        rtol=0.0,
        atol=1.0e-15,
    )


def test_independent_fd_curl_of_localized_vector_potential_matches_velocity():
    radius = 1.22
    z_normalized = 0.23
    h = 2.0e-5
    center = _run(radius, z_normalized)
    rp = _run(radius + h, z_normalized)
    rm = _run(radius - h, z_normalized)
    zp = _run(radius, z_normalized + h)
    zm = _run(radius, z_normalized - h)

    dr_a = (rp.complex_vector_potential - rm.complex_vector_potential) / (2.0 * h)
    dz_a = EPSILON * (
        zp.complex_vector_potential - zm.complex_vector_potential
    ) / (2.0 * h)
    a = center.complex_vector_potential
    # p=0, so there is no theta dependence in this manufactured phase/support.
    fd_curl = np.asarray(
        (
            -dz_a[1],
            dz_a[0] - dr_a[2],
            dr_a[1] + a[1] / radius,
        )
    )
    np.testing.assert_allclose(
        center.complex_velocity_cylindrical,
        fd_curl,
        rtol=3.0e-6,
        atol=3.0e-8,
    )


def test_chart_cartesian_rotation_is_only_a_basis_rotation():
    out = _run()
    u_r, u_theta, u_z = out.real_pair_velocity_cylindrical
    c = np.cos(THETA)
    s = np.sin(THETA)
    expected = np.asarray((u_r * c - u_theta * s, u_r * s + u_theta * c, u_z))
    np.testing.assert_allclose(
        out.real_pair_velocity_chart_cartesian,
        expected,
        rtol=0.0,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        np.linalg.norm(out.real_pair_velocity_chart_cartesian),
        np.linalg.norm(out.real_pair_velocity_cylindrical),
        rtol=0.0,
        atol=2.0e-14,
    )


def test_zero_cutoff_and_zero_support_jet_fail_to_exact_zero_not_posthoc_mask():
    cutoff_zero = _run(pulse_cutoff=0.0)
    support_zero = _run(
        pulse_cutoff=1.0,
        eta_override=(0.0, 0.0, 0.0),
    )
    for out in (cutoff_zero, support_zero):
        assert np.array_equal(out.localized_coefficient, np.zeros_like(out.localized_coefficient))
        assert np.array_equal(
            out.complex_vector_potential, np.zeros_like(out.complex_vector_potential)
        )
        assert np.array_equal(
            out.complex_velocity_cylindrical, np.zeros_like(out.complex_velocity_cylindrical)
        )
        assert np.array_equal(
            out.real_pair_velocity_chart_cartesian,
            np.zeros_like(out.real_pair_velocity_chart_cartesian),
        )


def test_batch_evaluation_preserves_trailing_three_vector_shape():
    radius = np.asarray((1.16, 1.24, 1.31))
    z = np.asarray((-0.18, 0.04, 0.21))
    eta, dr_eta, dz_eta = _eta_data(radius, z)
    out = materialize_source_zero_data_localized_harmonic(
        PULSE_V,
        radius,
        THETA,
        z,
        _zero_background(),
        _forcing,
        SourceLocalizationJet(eta, dr_eta, dz_eta, 0.75),
        epsilon=EPSILON,
        p=P,
        p_z=P_Z,
        x_0=X_0,
        k=K,
        m=M,
    )
    assert out.complex_velocity_cylindrical.shape == (3, 3)
    assert out.real_pair_velocity_chart_cartesian.shape == (3, 3)
    assert np.all(np.isfinite(out.real_pair_velocity_chart_cartesian))


@pytest.mark.parametrize(
    "localization,m,match",
    [
        (SourceLocalizationJet(1.0001, 0.0, 0.0, 1.0), 1, "squared partition"),
        (SourceLocalizationJet(0.5, 0.0, 0.0, -0.1), 1, "pulse_cutoff"),
        (SourceLocalizationJet(0.5, 0.0, 0.0, 1.1), 1, "pulse_cutoff"),
        (SourceLocalizationJet(0.5, 0.0, 0.0, 1.0), 0, "positive integer"),
        (SourceLocalizationJet(0.5, 0.0, 0.0, 1.0), -1, "positive"),
    ],
)
def test_invalid_localization_or_pair_representative_fails_closed(localization, m, match):
    with pytest.raises(ValueError, match=match):
        materialize_source_zero_data_localized_harmonic(
            PULSE_V,
            1.2,
            THETA,
            0.1,
            _zero_background(),
            _forcing,
            localization,
            epsilon=EPSILON,
            p=P,
            p_z=P_Z,
            x_0=X_0,
            k=K,
            m=m,
        )


def test_public_surface_exposes_no_integrator_or_residual_tuning_knob():
    params = inspect.signature(materialize_source_zero_data_localized_harmonic).parameters
    forbidden = {
        "panels",
        "steps",
        "rtol",
        "atol",
        "max_step",
        "residual",
        "defect",
        "pressure",
        "gain",
        "threshold",
    }
    assert forbidden.isdisjoint(params)


def test_contract_keeps_upstream_inputs_and_scientific_promotions_fail_closed():
    contract = source_zero_data_localized_harmonic_contract()
    assert contract["source_amplitude_directional_jet_consumed"] is True
    assert contract["potential_level_product_rule_localization_executable"] is True
    assert contract["cutoff_gradient_complete_curl_terms_retained"] is True
    assert contract["real_conjugate_pair_velocity_materialized_in_source_chart"] is True
    assert contract["caller_supplies_mode_forcing_directional_jet"] is True
    assert contract["caller_supplies_corrected_background"] is True
    assert contract["caller_supplies_source_partition_jet"] is True
    assert contract["caller_supplies_pulse_cutoff_value"] is True
    assert contract["concrete_source_partition_bumps_reconstructed"] is False
    assert contract["source_forcing_provider_materialized"] is False
    assert contract["actual_corrected_background_provider_materialized"] is False
    assert contract["source_physical_Q_scaling_applied"] is False
    assert contract["project_domain_coordinate_map_applied"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False
