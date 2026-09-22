from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_fixedq_preflight import (
    MAP_RELATIVE_GATE,
    Q_RELATIVE_SPREAD_GATE,
    CurrentI4RF30FixedQPreflightError,
    _canonical_band,
    materialize_current_i4_rf30_fixedq_preflight,
    truth_boundary,
)


class _FakeLeading:
    A = 0.6
    D = 0.4
    X_I4_start = 0.5
    X_I4_end = 1.5
    semantic_sha256 = "1" * 64

    @staticmethod
    def similarity_coordinates(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if np.max(np.abs(z)) != 0.0:
            raise ValueError("fake mechanics backend only implements z=0")
        q = 1.0 - t
        eta = np.zeros_like(q)
        X = (x * x + y * y) / (2.0 * q)
        return {"q": q, "eta": eta, "X": X}

    @staticmethod
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((x, -y, x + z), axis=-1)


class _FakeComposite:
    leading_backend = _FakeLeading()


def _fake_differential(points, time):
    p = np.asarray(points, dtype=float)
    if p.shape[-1] != 3:
        raise ValueError
    x = p[..., 0]
    y = p[..., 1]
    z = p[..., 2]
    velocity = np.stack((x, -y, x + z), axis=-1)
    jac = np.zeros(p.shape[:-1] + (3, 3), dtype=float)
    jac[..., 0, 0] = 1.0
    jac[..., 1, 1] = -1.0
    jac[..., 2, 0] = 1.0
    jac[..., 2, 2] = 1.0
    return SimpleNamespace(
        velocity=velocity,
        velocity_jacobian=jac,
        interior_mask=np.ones(p.shape[:-1], dtype=bool),
    )


def _backend() -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=_FakeComposite(),
        differential_function=_fake_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256="3" * 64,
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def test_canonical_band_is_frozen_by_q_only() -> None:
    ell, Q, s_q = _canonical_band(0.5)
    assert ell == 1
    assert Q == 0.5
    assert s_q == 1.0

    ell, Q, s_q = _canonical_band(0.3)
    assert ell == 2
    assert Q == 0.25
    assert 1.0 <= s_q < 2.0

    with pytest.raises(CurrentI4RF30FixedQPreflightError):
        _canonical_band(0.0)
    with pytest.raises(CurrentI4RF30FixedQPreflightError):
        _canonical_band(1.1)


def test_fixed_q_map_and_observable_covariance_are_candidate_bound_but_not_rf30_haar() -> None:
    radius = np.asarray((0.8, 0.9, 1.0), dtype=float)
    receipt = materialize_current_i4_rf30_fixedq_preflight(
        _backend(), radius, 0.0, 0.5
    )

    assert receipt.Q == 0.5
    assert receipt.s_Q == 1.0
    assert receipt.all_strict_I4 is True
    assert receipt.fixed_q_map_closure_relative_max <= MAP_RELATIVE_GATE
    assert receipt.q_relative_spread <= Q_RELATIVE_SPREAD_GATE
    assert receipt.finest_oscillatory_interior_fraction == 1.0
    assert receipt.finest_covariance_trace_rms > 0.0
    assert receipt.normalized_physical_azimuthal_mean_used is True
    assert receipt.normalized_source_auxiliary_t2_haar_mean_used is False
    assert receipt.source_auxiliary_t2_provider_available is False
    assert receipt.rf30_repository_candidate_state_authorized is False

    Q = receipt.Q
    A = receipt.A
    expected_wrr = 0.5 * (Q ** (2.0 * A)) * radius * radius
    expected_wtt = expected_wrr
    expected_wzz = expected_wrr
    np.testing.assert_allclose(receipt.mean_W_rr_physical_theta, expected_wrr, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(receipt.mean_W_thetatheta_physical_theta, expected_wtt, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(receipt.mean_W_zz_physical_theta, expected_wzz, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(receipt.mean_W_zr_physical_theta, 0.0, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(receipt.mean_W_ztheta_physical_theta, 0.0, rtol=0.0, atol=2.0e-14)

    R = radius / np.sqrt(Q)
    expected_dR_wrr = (Q ** (2.0 * A + 1.0)) * R
    np.testing.assert_allclose(
        receipt.dR_mean_W_rr_physical_theta,
        expected_dR_wrr,
        rtol=0.0,
        atol=3.0e-14,
    )
    np.testing.assert_allclose(receipt.dZ_mean_W_zr_physical_theta, 0.0, rtol=0.0, atol=2.0e-14)
    assert max(receipt.successive_covariance_relative_differences) <= 2.0e-14
    assert len(receipt.preflight_sha256) == 64


def test_preflight_rejects_non_exact_backend_type_and_bad_radial_grid() -> None:
    with pytest.raises(TypeError):
        materialize_current_i4_rf30_fixedq_preflight(object(), (0.8, 0.9, 1.0), 0.0, 0.5)

    with pytest.raises(CurrentI4RF30FixedQPreflightError):
        materialize_current_i4_rf30_fixedq_preflight(
            _backend(), (0.8, 0.8, 1.0), 0.0, 0.5
        )


def test_truth_boundary_never_promotes_physical_theta_mean_to_source_haar_or_residual() -> None:
    truth = truth_boundary()
    assert truth["fixed_q_coordinate_map_materialized"] is True
    assert truth["actual_a2_oscillatory_jacobian_consumed"] is True
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["physical_azimuthal_covariance_materialized"] is True
    assert truth["normalized_source_auxiliary_t2_haar_mean_used"] is False
    assert truth["physical_theta_mean_promoted_to_source_auxiliary_haar"] is False
    assert truth["rf30_repository_candidate_state_authorized"] is False
    assert truth["rf30_defect_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["forbidden_scientific_controls_exposed"] is False
