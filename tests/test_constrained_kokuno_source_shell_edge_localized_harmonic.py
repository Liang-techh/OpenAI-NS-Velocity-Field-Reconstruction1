from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
)
from openai_ns_reconstruction.kokuno_source_shell_edge_localized_harmonic import (
    SourceShellCoordinateJet,
    apply_source_shell_edge_to_localization,
    materialize_source_zero_data_shell_localized_harmonic,
    source_shell_edge_jet,
    source_shell_edge_localized_harmonic_contract,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_sensitivity import (
    SourceModeForcingDirectionalJet,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizationJet,
)


def _background() -> SourceBackgroundPhaseJet:
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


def _forcing(_v: float) -> SourceModeForcingDirectionalJet:
    value = np.asarray((0.0 + 0.0j, 0.7 + 0.2j, -0.35 + 0.1j))
    zero = np.zeros(3, dtype=np.complex128)
    return SourceModeForcingDirectionalJet(value, zero, zero)


def _formula(x: float, x_a: float, x_b: float, a_a: float, a_b: float) -> float:
    if not x_a < x < x_b:
        return 0.0
    d_a = np.log(x / x_a)
    d_b = np.log(x_b / x)
    return float(np.exp(-0.5 * (a_a / d_a**2 + a_b / d_b**2)))


def test_shell_edge_matches_source_formula_and_independent_fd_derivative() -> None:
    x = 2.0
    x_a, x_b = 1.0, 4.0
    a_a, a_b = 1.3, 0.7
    jet = source_shell_edge_jet(
        SourceShellCoordinateJet(x=x, dr_x=0.23, dz_x=-0.17),
        x_a=x_a,
        x_b=x_b,
        a_a=a_a,
        a_b=a_b,
    )
    assert float(jet.weight) == pytest.approx(_formula(x, x_a, x_b, a_a, a_b), abs=2e-15)

    h = 2.0e-6
    fd = (
        _formula(x + h, x_a, x_b, a_a, a_b)
        - _formula(x - h, x_a, x_b, a_a, a_b)
    ) / (2.0 * h)
    assert float(jet.dx_weight) == pytest.approx(fd, rel=2e-9, abs=2e-11)
    assert float(jet.dr_weight) == pytest.approx(fd * 0.23, rel=2e-9, abs=2e-11)
    assert float(jet.dz_weight) == pytest.approx(fd * -0.17, rel=2e-9, abs=2e-11)


def test_shell_edge_batch_and_smooth_zero_extension() -> None:
    x = np.asarray((0.7, 1.0, 1.000000001, 2.0, 3.999999999, 4.0, 5.0))
    jet = source_shell_edge_jet(
        SourceShellCoordinateJet(x=x, dr_x=np.ones_like(x), dz_x=-2.0),
        x_a=1.0,
        x_b=4.0,
        a_a=1.0,
        a_b=1.0,
    )
    assert jet.weight.shape == x.shape
    assert np.array_equal(jet.inside_shell, (x > 1.0) & (x < 4.0))
    assert jet.weight[0] == 0.0
    assert jet.weight[1] == 0.0
    assert jet.weight[-2] == 0.0
    assert jet.weight[-1] == 0.0
    assert jet.dx_weight[0] == 0.0
    assert jet.dx_weight[1] == 0.0
    assert jet.dx_weight[-2] == 0.0
    assert jet.dx_weight[-1] == 0.0
    assert np.all(np.isfinite(jet.weight))
    assert np.all(np.isfinite(jet.dx_weight))
    assert np.all(np.isfinite(jet.dr_weight))
    assert np.all(np.isfinite(jet.dz_weight))


def test_shell_edge_is_multiplied_at_partition_jet_level() -> None:
    base = SourceLocalizationJet(
        eta=np.asarray((0.8, 0.6)),
        dr_eta=np.asarray((0.05, -0.03)),
        dz_eta=np.asarray((-0.04, 0.02)),
        pulse_cutoff=0.73,
    )
    shell = source_shell_edge_jet(
        SourceShellCoordinateJet(
            x=np.asarray((1.7, 2.6)),
            dr_x=np.asarray((0.11, -0.08)),
            dz_x=np.asarray((-0.07, 0.12)),
        ),
        x_a=1.0,
        x_b=4.0,
        a_a=0.9,
        a_b=1.1,
    )
    combined = apply_source_shell_edge_to_localization(base, shell)
    np.testing.assert_allclose(combined.eta, np.asarray(base.eta) * shell.weight, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        combined.dr_eta,
        np.asarray(base.dr_eta) * shell.weight + np.asarray(base.eta) * shell.dr_weight,
        rtol=2e-15,
        atol=2e-15,
    )
    np.testing.assert_allclose(
        combined.dz_eta,
        np.asarray(base.dz_eta) * shell.weight + np.asarray(base.eta) * shell.dz_weight,
        rtol=2e-15,
        atol=2e-15,
    )
    assert combined.pulse_cutoff == base.pulse_cutoff


def test_wrapper_consumes_shell_before_complete_curl_and_zeroes_outside() -> None:
    common = dict(
        pulse_v=0.30,
        radius=1.22,
        theta=0.41,
        z_normalized=0.23,
        background=_background(),
        forcing_jet_m=_forcing,
        localization_without_shell_edge=SourceLocalizationJet(
            eta=0.8,
            dr_eta=0.04,
            dz_eta=-0.03,
            pulse_cutoff=0.75,
        ),
        x_a=1.0,
        x_b=4.0,
        a_a=0.8,
        a_b=0.9,
        epsilon=0.25,
        p=0.0,
        p_z=0.0,
        x_0=1.0,
        k=2.0,
        m=1,
    )
    interior = materialize_source_zero_data_shell_localized_harmonic(
        shell_coordinate=SourceShellCoordinateJet(x=2.0, dr_x=0.14, dz_x=-0.09),
        **common,
    )
    outside = materialize_source_zero_data_shell_localized_harmonic(
        shell_coordinate=SourceShellCoordinateJet(x=4.2, dr_x=0.14, dz_x=-0.09),
        **common,
    )
    assert float(interior.shell_edge.weight) > 0.0
    assert np.max(np.abs(interior.harmonic.complex_velocity_cylindrical)) > 1.0e-12
    assert float(outside.shell_edge.weight) == 0.0
    assert np.max(np.abs(outside.harmonic.complex_velocity_cylindrical)) == 0.0
    assert np.max(np.abs(outside.harmonic.real_pair_velocity_chart_cartesian)) == 0.0


def test_invalid_shell_inputs_fail_closed() -> None:
    coordinate = SourceShellCoordinateJet(x=2.0, dr_x=0.0, dz_x=0.0)
    for kwargs in (
        dict(x_a=0.0, x_b=4.0, a_a=1.0, a_b=1.0),
        dict(x_a=1.0, x_b=1.0, a_a=1.0, a_b=1.0),
        dict(x_a=1.0, x_b=4.0, a_a=0.0, a_b=1.0),
        dict(x_a=1.0, x_b=4.0, a_a=1.0, a_b=np.inf),
    ):
        with pytest.raises(ValueError):
            source_shell_edge_jet(coordinate, **kwargs)
    with pytest.raises(ValueError):
        source_shell_edge_jet(
            SourceShellCoordinateJet(x=np.nan, dr_x=0.0, dz_x=0.0),
            x_a=1.0,
            x_b=4.0,
            a_a=1.0,
            a_b=1.0,
        )
    with pytest.raises(ValueError):
        source_shell_edge_jet(
            SourceShellCoordinateJet(x=-2.0, dr_x=0.0, dz_x=0.0),
            x_a=1.0,
            x_b=4.0,
            a_a=1.0,
            a_b=1.0,
        )


def test_contract_preserves_source_realization_boundary() -> None:
    contract = source_shell_edge_localized_harmonic_contract()
    assert contract["task"] == "K2-OSC-099"
    assert contract["source_sqrt_zeta_shell_edge_materialized"] is True
    assert contract["source_shell_analytic_x_derivative_materialized"] is True
    assert contract["source_shell_normalized_directional_pullback_materialized"] is True
    assert contract["source_shell_smooth_zero_extension_materialized"] is True
    assert contract["source_shell_edge_consumed_before_complete_curl"] is True
    assert contract["caller_supplies_shell_coordinate_directional_jet"] is True
    assert contract["caller_supplies_remaining_partition_jet"] is True
    assert contract["caller_supplies_pulse_cutoff_value"] is True
    assert contract["concrete_source_chi_partition_bumps_reconstructed"] is False
    assert contract["source_pulse_cutoff_provider_materialized"] is False
    assert contract["source_forcing_provider_materialized"] is False
    assert contract["actual_corrected_background_provider_materialized"] is False
    assert contract["source_physical_Q_scaling_applied"] is False
    assert contract["project_domain_coordinate_map_applied"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_public_surfaces_have_no_residual_or_accuracy_tuning_knobs() -> None:
    forbidden = {
        "residual",
        "target",
        "gain",
        "forcing",
        "pressure",
        "rtol",
        "atol",
        "steps",
        "panels",
        "max_step",
    }
    for function in (
        source_shell_edge_jet,
        apply_source_shell_edge_to_localization,
        materialize_source_zero_data_shell_localized_harmonic,
    ):
        parameters = set(inspect.signature(function).parameters)
        assert forbidden.isdisjoint(parameters)
