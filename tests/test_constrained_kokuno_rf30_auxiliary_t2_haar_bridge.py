from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_fixedq_preflight import (
    materialize_current_i4_rf30_fixedq_preflight,
)
from openai_ns_reconstruction.kokuno_rf30_auxiliary_t2_haar_bridge import (
    TORUS_ORDERS,
    AuxiliaryT2HaarBridgeError,
    AuxiliaryT2ProviderMetadata,
    AuxiliaryT2WaveSamples,
    materialize_current_i4_rf30_auxiliary_haar_bridge,
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


def _metadata() -> AuxiliaryT2ProviderMetadata:
    preflight = materialize_current_i4_rf30_fixedq_preflight(
        _backend(), np.asarray((0.8, 0.9, 1.0)), 0.0, 0.5
    )
    return AuxiliaryT2ProviderMetadata(
        candidate_semantic_sha256=preflight.candidate_semantic_sha256,
        oscillatory_runtime_sha256=preflight.oscillatory_runtime_sha256,
        leading_semantic_sha256=preflight.leading_semantic_sha256,
        provider_semantic_sha256="5" * 64,
        source_blob_sha1="6" * 40,
        provider_kind="mechanics_fixture",
        source_auxiliary_t2_field_materialized=True,
        normalized_haar_measure_total_mass_one=True,
        recomputed_from_actual_candidate=False,
        surrogate_or_preaveraged_covariance_used=False,
        heldout_data_used=False,
        residual_as_forcing_used=False,
    )


def test_auxiliary_t2_normalized_haar_covariance_and_derivatives_are_constructed_in_bridge() -> None:
    radius = np.asarray((0.8, 0.9, 1.0), dtype=float)
    receipt = materialize_current_i4_rf30_auxiliary_haar_bridge(
        _backend(),
        _TrigonometricAuxT2Provider(_metadata()),
        radius,
        0.0,
        0.5,
    )

    R = np.asarray(receipt.radius_R)
    A = 1.0 + R
    B = 1.0 + 2.0 * R + receipt.Z
    C = 2.0 - R

    np.testing.assert_allclose(receipt.mean_W_rr_auxiliary_haar, 0.5 * A * A, atol=2e-14)
    np.testing.assert_allclose(receipt.mean_W_zr_auxiliary_haar, 0.5 * A * B, atol=2e-14)
    np.testing.assert_allclose(
        receipt.mean_W_thetatheta_auxiliary_haar, 0.5 * C * C, atol=2e-14
    )
    np.testing.assert_allclose(receipt.mean_W_ztheta_auxiliary_haar, 0.0, atol=2e-14)
    np.testing.assert_allclose(receipt.mean_W_zz_auxiliary_haar, 0.5 * B * B, atol=2e-14)
    np.testing.assert_allclose(receipt.dR_mean_W_rr_auxiliary_haar, A, atol=2e-14)
    np.testing.assert_allclose(receipt.dZ_mean_W_zr_auxiliary_haar, 0.5 * A, atol=2e-14)
    np.testing.assert_allclose(receipt.mean_w_r_auxiliary_haar, 0.0, atol=2e-14)
    np.testing.assert_allclose(receipt.mean_w_theta_auxiliary_haar, 0.0, atol=2e-14)
    np.testing.assert_allclose(receipt.mean_w_z_auxiliary_haar, 0.0, atol=2e-14)

    assert receipt.torus_orders == TORUS_ORDERS
    assert max(receipt.successive_covariance_relative_differences) <= 2.0e-14
    assert receipt.normalized_source_auxiliary_t2_haar_mean_used is True
    assert receipt.covariance_formed_in_bridge_from_raw_wave_samples is True
    assert receipt.covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives is True
    assert receipt.provider_checksum_pinned_for_repository_candidate is False
    assert receipt.rf30_covariance_bridge_materialized is True
    assert receipt.rf30_repository_candidate_state_authorized is False
    assert len(receipt.receipt_sha256) == 64


@pytest.mark.parametrize(
    "field,value",
    [
        ("surrogate_or_preaveraged_covariance_used", True),
        ("heldout_data_used", True),
        ("residual_as_forcing_used", True),
        ("normalized_haar_measure_total_mass_one", False),
        ("source_auxiliary_t2_field_materialized", False),
    ],
)
def test_bridge_fails_closed_on_forbidden_provider_provenance(field, value) -> None:
    metadata = replace(_metadata(), **{field: value})
    with pytest.raises(AuxiliaryT2HaarBridgeError):
        materialize_current_i4_rf30_auxiliary_haar_bridge(
            _backend(),
            _TrigonometricAuxT2Provider(metadata),
            (0.8, 0.9, 1.0),
            0.0,
            0.5,
        )


def test_bridge_fails_closed_on_candidate_identity_drift() -> None:
    metadata = replace(_metadata(), candidate_semantic_sha256="a" * 64)
    with pytest.raises(AuxiliaryT2HaarBridgeError):
        materialize_current_i4_rf30_auxiliary_haar_bridge(
            _backend(),
            _TrigonometricAuxT2Provider(metadata),
            (0.8, 0.9, 1.0),
            0.0,
            0.5,
        )


def test_unpinned_repository_provider_cannot_become_scientific_evidence() -> None:
    metadata = replace(
        _metadata(),
        provider_kind="repository_candidate",
        recomputed_from_actual_candidate=True,
    )
    receipt = materialize_current_i4_rf30_auxiliary_haar_bridge(
        _backend(),
        _TrigonometricAuxT2Provider(metadata),
        (0.8, 0.9, 1.0),
        0.0,
        0.5,
    )
    assert receipt.normalized_source_auxiliary_t2_haar_mean_used is True
    assert receipt.provider_checksum_pinned_for_repository_candidate is False
    assert receipt.rf30_repository_candidate_state_authorized is False


def test_truth_boundary_keeps_operator_mechanics_separate_from_rf30_and_ns_claims() -> None:
    truth = truth_boundary()
    assert truth["normalized_source_auxiliary_t2_haar_operator_executable"] is True
    assert truth["haar_measure_total_mass_one"] is True
    assert truth["covariance_formed_in_bridge_from_raw_wave_samples"] is True
    assert truth["covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives"] is True
    assert truth["preaveraged_covariance_input_exposed"] is False
    assert truth["physical_theta_mean_promoted_to_source_auxiliary_haar"] is False
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["repository_provider_blob_pinned"] is False
    assert truth["rf30_repository_candidate_state_authorized"] is False
    assert truth["rf30_defect_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["forbidden_scientific_controls_exposed"] is False
