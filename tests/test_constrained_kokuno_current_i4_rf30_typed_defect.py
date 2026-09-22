from __future__ import annotations

from dataclasses import replace
import math
from types import SimpleNamespace

import numpy as np

from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_fixedq_preflight import (
    materialize_current_i4_rf30_fixedq_preflight,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_typed_defect import (
    materialize_current_i4_rf30_typed_defect,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_rf30_auxiliary_t2_haar_bridge import (
    AuxiliaryT2ProviderMetadata,
    AuxiliaryT2WaveSamples,
)


class _FakeSchedule:
    lambda_outer = 0.05
    h = 0.005
    log_c_patch = math.log(0.8)

    @staticmethod
    def reserved_log_intervals():
        return {
            "I1": (-6.0, -5.0),
            "I2": (-4.0, -3.0),
            "I3": (-2.0, -1.0),
            "I4": (math.log(0.5), math.log(1.5)),
        }


class _FakeLeading:
    A = 0.505
    D = 0.495
    X_I4_start = 0.5
    X_I4_end = 1.5
    semantic_sha256 = "1" * 64
    outer_schedule = _FakeSchedule()

    @staticmethod
    def similarity_coordinates(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if np.max(np.abs(z)) != 0.0:
            raise ValueError("fixture implements z=0 only")
        q = 1.0 - t
        return {
            "q": q,
            "eta": np.zeros_like(q),
            "X": (x * x + y * y) / (2.0 * q),
        }

    @staticmethod
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((x, -y, x + z), axis=-1)


class _FakeComposite:
    leading_backend = _FakeLeading()


def _fake_differential(points, time):
    p = np.asarray(points, dtype=float)
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


class _TrigonometricAuxT2Provider:
    def __init__(self, metadata: AuxiliaryT2ProviderMetadata):
        self.metadata = metadata

    def sample_auxiliary_wave(self, *, R, Z, T, y1, y2):
        R = np.asarray(R, dtype=float)[:, None, None]
        y1 = np.asarray(y1, dtype=float)[None, :, :]
        y2 = np.asarray(y2, dtype=float)[None, :, :]
        c1 = np.cos(2.0 * np.pi * y1)
        s2 = np.sin(2.0 * np.pi * y2)

        A = 1.0 + R
        B = 1.0 + 2.0 * R + float(Z)
        C = 2.0 - R
        wr = A * c1
        wt = C * s2
        wz = B * c1
        zero = np.zeros_like(wr)
        return AuxiliaryT2WaveSamples(
            w_r=wr,
            w_theta=wt,
            w_z=wz,
            dR_w_r=c1 + zero,
            dR_w_theta=-s2 + zero,
            dR_w_z=2.0 * c1 + zero,
            dZ_w_r=zero,
            dZ_w_theta=zero,
            dZ_w_z=c1 + zero,
        )


def _metadata(*, repository_candidate: bool = False) -> AuxiliaryT2ProviderMetadata:
    preflight = materialize_current_i4_rf30_fixedq_preflight(
        _backend(), np.asarray((0.8, 0.9, 1.0)), 0.0, 0.5
    )
    return AuxiliaryT2ProviderMetadata(
        candidate_semantic_sha256=preflight.candidate_semantic_sha256,
        oscillatory_runtime_sha256=preflight.oscillatory_runtime_sha256,
        leading_semantic_sha256=preflight.leading_semantic_sha256,
        provider_semantic_sha256="5" * 64,
        source_blob_sha1="6" * 40,
        provider_kind="repository_candidate" if repository_candidate else "mechanics_fixture",
        source_auxiliary_t2_field_materialized=True,
        normalized_haar_measure_total_mass_one=True,
        recomputed_from_actual_candidate=repository_candidate,
        surrogate_or_preaveraged_covariance_used=False,
        heldout_data_used=False,
        residual_as_forcing_used=False,
    )


def _trapz(values: np.ndarray, R: np.ndarray) -> float:
    return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(R)))


def test_auxiliary_haar_state_materializes_rf30_defect_without_manual_PJ() -> None:
    radius = np.asarray((0.8, 0.9, 1.0), dtype=float)
    result = materialize_current_i4_rf30_typed_defect(
        _backend(),
        _TrigonometricAuxT2Provider(_metadata()),
        radius,
        0.0,
        0.5,
    )

    state = result.state
    R = np.asarray(state.radius_R, dtype=float)
    Wrr = np.asarray(state.mean_W_rr, dtype=float)
    Wtt = np.asarray(state.mean_W_thetatheta, dtype=float)
    Wzt = np.asarray(state.mean_W_ztheta, dtype=float)
    Wzz = np.asarray(state.mean_W_zz, dtype=float)
    dR_Wrr = np.asarray(state.dR_mean_W_rr, dtype=float)
    dZ_Wzr = np.asarray(state.dZ_mean_W_zr, dtype=float)
    g_r = -(dR_Wrr + Wrr / R) - dZ_Wzr + Wtt / R
    expected_P = _trapz(g_r, R)
    expected_Jtheta = _trapz((R**2) * Wzt, R)
    expected_Jz = _trapz(R * Wzz, R) - 0.5 * _trapz((R**2) * g_r, R)

    np.testing.assert_allclose(result.defect.P, expected_P, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(result.defect.J_theta, expected_Jtheta, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(result.defect.J_z, expected_Jz, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(result.defect.radial_source_g_r, g_r, rtol=0.0, atol=2e-13)

    lam = _FakeSchedule.lambda_outer
    h = _FakeSchedule.h
    a = (2.0 ** (0.5 + lam)) * math.exp(_FakeSchedule.log_c_patch)
    expected_V = a * np.power(R, -1.0 - 2.0 * lam)
    np.testing.assert_allclose(state.base_V, expected_V, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(state.base_G, 0.0, rtol=0.0, atol=0.0)

    assert result.source_background_from_exact_a1_i4_schedule is True
    assert result.normalized_source_auxiliary_t2_haar_mean_used is True
    assert result.typed_rf30_state_materialized is True
    assert result.rf30_defect_materialized is True
    assert result.repository_candidate_defect_evidence is False
    assert result.candidate.evidence_kind == "mechanics-only"
    assert result.rf31_five_row_system_materialized is False
    assert result.correction_applied is False
    assert result.heldout_ns_residual_assessed is False
    assert len(result.source_chart_sha256) == 64
    assert len(result.receipt_sha256) == 64


def test_unpinned_repository_provider_still_cannot_promote_rf30_defect() -> None:
    provider = _TrigonometricAuxT2Provider(_metadata(repository_candidate=True))
    result = materialize_current_i4_rf30_typed_defect(
        _backend(), provider, (0.8, 0.9, 1.0), 0.0, 0.5
    )
    assert result.provider_kind == "repository_candidate"
    assert result.state.actual_candidate_recomputed is True
    assert result.state.normalized_haar_mean_used is True
    assert result.repository_candidate_defect_evidence is False
    assert result.candidate.evidence_kind == "mechanics-only"


def test_provider_forbidden_provenance_is_rejected_upstream() -> None:
    provider = _TrigonometricAuxT2Provider(
        replace(_metadata(), heldout_data_used=True)
    )
    try:
        materialize_current_i4_rf30_typed_defect(
            _backend(), provider, (0.8, 0.9, 1.0), 0.0, 0.5
        )
    except Exception as exc:
        assert "held-out" in str(exc)
    else:
        raise AssertionError("held-out-constructed auxiliary provider must fail closed")


def test_truth_boundary_keeps_rf30_defect_separate_from_correction_and_ns_claims() -> None:
    truth = truth_boundary()
    assert truth["normalized_auxiliary_t2_haar_bridge_consumed"] is True
    assert truth["physical_theta_mean_promoted_to_source_auxiliary_haar"] is False
    assert truth["source_i4_background_derived_from_exact_a1_schedule"] is True
    assert truth["caller_supplied_base_background_allowed"] is False
    assert truth["caller_supplied_rf30_defect_allowed"] is False
    assert truth["typed_rf30_state_adapter_executable"] is True
    assert truth["rf30_defect_operator_reused_from_1126"] is True
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["repository_provider_blob_pinned_in_parent"] is False
    assert truth["current_i4_repository_candidate_rf30_defect_materialized"] is False
    assert truth["rf31_five_row_system_materialized"] is False
    assert truth["correction_applied"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
