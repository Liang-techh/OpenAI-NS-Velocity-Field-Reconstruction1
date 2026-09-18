from __future__ import annotations

import json
from math import pi

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_phase import (
    KokunoOscillatoryPhaseContract,
)
from openai_ns_reconstruction.kokuno_source_support_localized_curl import (
    KokunoSourceSupportLocalizedCurl,
)
from openai_ns_reconstruction.kokuno_source_wave_shell import (
    KokunoSourceWaveShellLocalizedCurl,
)


def _phase_contract(*, epsilon=0.2, m=2, j=1):
    k = int(np.ceil(epsilon ** -0.5))
    return KokunoOscillatoryPhaseContract(
        p=j / k,
        p_z=0.62,
        x0=1.08,
        epsilon=epsilon,
        m=m,
    )


def _transverse_mode(phase_contract, R):
    n = np.asarray(
        [1.03, phase_contract.p / R, 0.57],
        dtype=float,
    )
    raw = np.asarray([0.71 + 0.04j, -0.24 + 0.03j, 0.19 - 0.02j])
    t = raw - n * (np.dot(n, raw) / np.dot(n, n))
    return n, t


def _mode_inputs(phase_contract, R):
    n, t = _transverse_mode(phase_contract, R)
    D_r_C = np.asarray([0.021 + 0.013j, -0.018 + 0.007j, 0.011 - 0.009j])
    D_z_C = np.asarray([-0.014 + 0.008j, 0.009 - 0.004j, 0.017 + 0.012j])
    return n, t, D_r_C, D_z_C


def test_source_grid_makes_fast_harmonic_integer_and_2pi_periodic():
    phase_contract = _phase_contract(epsilon=0.2, m=2, j=-2)
    shell = KokunoSourceWaveShellLocalizedCurl(Lambda=16.0)
    assert phase_contract.k == 3
    assert phase_contract.j == -2
    assert shell.azimuthal_mode_number(phase_contract) == -4
    assert phase_contract.k_m * phase_contract.p == pytest.approx(-4.0, abs=2e-12)

    R = 0.72
    n, t, D_r_C, D_z_C = _mode_inputs(phase_contract, R)
    kwargs = dict(
        phase_contract=phase_contract,
        R=R,
        s_Q=0.5,
        n_phi=n,
        t_m=t,
        D_r_C_m=D_r_C,
        D_z_C_m=D_z_C,
        eta=0.61,
        D_r_eta=0.07,
        D_z_eta=-0.015,
    )
    base_phase = 0.37
    a = shell.localized_mode(phase=base_phase, **kwargs)
    b = shell.localized_mode(phase=base_phase + 2.0 * pi * phase_contract.p, **kwargs)
    np.testing.assert_allclose(a["vector_potential"], b["vector_potential"], rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(a["velocity"], b["velocity"], rtol=2e-12, atol=2e-12)

    # Agent-4's earlier manufactured p=1.35 is useful as a generic local
    # structural probe, but it is not on the source angular grid for eps=.2:
    # k*p=3*1.35=4.05 is not an integer.  The source phase contract therefore
    # fails closed before this source-bound wrapper can consume it.
    with pytest.raises(ValueError, match="p must satisfy p=j/k"):
        KokunoOscillatoryPhaseContract(p=1.35, p_z=0.62, x0=1.08, epsilon=0.2, m=1)


def test_source_inner_annulus_is_an_explicit_axis_guard_and_batches():
    shell = KokunoSourceWaveShellLocalizedCurl(Lambda=16.0)
    assert shell.X_a == pytest.approx(0.25)
    assert shell.axis_lower_bound == pytest.approx(0.5)

    edge = shell.source_shell_coordinates(0.5, 0.5)
    assert edge["X"] == pytest.approx(0.25)
    batch = shell.source_shell_coordinates(
        np.asarray([0.5, 0.8, 1.1]),
        np.asarray([0.5, 0.8, 2.0]),
    )
    np.testing.assert_allclose(batch["X"], np.asarray([0.25, 0.4, 0.3025]))

    with pytest.raises(ValueError, match="below source active annulus"):
        shell.source_shell_coordinates(0.499, 0.5)
    with pytest.raises(ValueError, match="below source active annulus"):
        shell.source_shell_coordinates(0.99, 2.0)
    with pytest.raises(ValueError, match="1/2<=s_Q<=2"):
        shell.source_shell_coordinates(0.9, 0.49)
    with pytest.raises(ValueError, match="R>0"):
        shell.source_shell_coordinates(0.0, 0.5)


def test_wrapper_delegates_exactly_after_source_domain_checks():
    phase_contract = _phase_contract(epsilon=0.25, m=-1, j=2)
    shell = KokunoSourceWaveShellLocalizedCurl(Lambda=25.0)
    R = 0.73
    s_Q = 0.6
    phase = -0.29
    n, t, D_r_C, D_z_C = _mode_inputs(phase_contract, R)
    eta = 0.58
    D_r_eta = 0.06
    D_z_eta = -0.012

    wrapped = shell.localized_mode(
        phase_contract,
        R,
        s_Q,
        phase,
        n,
        t,
        D_r_C,
        D_z_C,
        eta,
        D_r_eta,
        D_z_eta,
    )
    raw = KokunoSourceSupportLocalizedCurl(
        epsilon=phase_contract.epsilon,
        m=phase_contract.m,
    ).localized_mode(
        R,
        phase,
        n,
        t,
        D_r_C,
        D_z_C,
        eta,
        D_r_eta,
        D_z_eta,
    )
    for key in (
        "C_local",
        "D_r_C_local",
        "D_z_C_local",
        "leading_local",
        "remainder_local",
        "vector_potential",
        "velocity",
    ):
        np.testing.assert_allclose(wrapped[key], raw[key], rtol=0.0, atol=0.0)
    assert wrapped["azimuthal_mode_number"] == phase_contract.m * phase_contract.j
    assert wrapped["X"] == pytest.approx(R * R / (2.0 * s_Q))

    bad_n = n.copy()
    bad_n[1] *= 1.01
    with pytest.raises(ValueError, match="theta component"):
        shell.localized_mode(
            phase_contract,
            R,
            s_Q,
            phase,
            bad_n,
            t,
            D_r_C,
            D_z_C,
            eta,
            D_r_eta,
            D_z_eta,
        )


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    shell = KokunoSourceWaveShellLocalizedCurl(Lambda=64.0)
    path = shell.save_json(tmp_path / "source_wave_shell.json")
    loaded = KokunoSourceWaveShellLocalizedCurl.load_json(path)
    assert loaded.sha256 == shell.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["source_angular_label_grid_enforced"] is True
    assert truth["source_inner_annulus_guard_executable"] is True
    assert truth["source_axis_separation_guard_executable"] is True
    assert truth["source_actual_Lambda_recovered"] is False
    assert truth["concrete_source_partition_bumps_reconstructed"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceWaveShellLocalizedCurl.from_payload(payload)
