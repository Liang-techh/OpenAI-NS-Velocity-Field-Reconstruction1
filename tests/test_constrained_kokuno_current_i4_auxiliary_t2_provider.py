from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i4_auxiliary_t2_provider import (
    CurrentI4AuxiliaryT2ProviderError,
    KokunoCurrentI4AuxiliaryT2Provider,
    _runtime_payload_sha,
    bind_current_i4_auxiliary_t2_provider,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_fixedq_preflight import (
    materialize_current_i4_rf30_fixedq_preflight,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import default_field
from openai_ns_reconstruction.kokuno_rf30_auxiliary_t2_haar_bridge import (
    AuxiliaryT2ProviderMetadata,
    AuxiliaryT2WaveSamples,
    materialize_current_i4_rf30_auxiliary_haar_bridge,
)


class _FakeLeading:
    A = 0.6
    D = 0.4
    X_I4_start = 0.10
    X_I4_end = 0.40
    semantic_sha256 = "1" * 64

    @staticmethod
    def similarity_coordinates(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        q = 1.0 - t
        eta = z / np.power(q, _FakeLeading.D)
        X = (x * x + y * y) / (2.0 * q)
        return {"q": q, "eta": eta, "X": X}

    @staticmethod
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((0.0 * x, 0.0 * y, 0.0 * z), axis=-1)


class _FakeComposite:
    leading_backend = _FakeLeading()


def _fake_differential(points, time):
    p = np.asarray(points, dtype=float)
    velocity = np.zeros_like(p)
    jac = np.zeros(p.shape[:-1] + (3, 3), dtype=float)
    return SimpleNamespace(
        velocity=velocity,
        velocity_jacobian=jac,
        interior_mask=np.ones(p.shape[:-1], dtype=bool),
    )


def _backend() -> ExactCurrentI4NonlinearBackend:
    runtime_sha = _runtime_payload_sha(default_field())
    return ExactCurrentI4NonlinearBackend(
        composite_field=_FakeComposite(),
        differential_function=_fake_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256=runtime_sha,
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _inputs():
    radius = np.asarray((0.40, 0.45, 0.50), dtype=float)
    return radius, 0.0, 0.5


def _provider() -> KokunoCurrentI4AuxiliaryT2Provider:
    radius, z, t = _inputs()
    return bind_current_i4_auxiliary_t2_provider(_backend(), radius, z, t)


def _cylindrical_total(field, r, z, t, theta):
    rr = np.asarray(r, dtype=float)[:, None, None]
    tt = np.asarray(theta, dtype=float)[None, :, :]
    x = rr * np.cos(tt)
    y = rr * np.sin(tt)
    out = field.evaluate(x, y, np.full_like(x, z), np.full_like(x, t))
    v = np.asarray(out["velocity_cartesian_total"], dtype=float)
    c = np.cos(tt)
    s = np.sin(tt)
    return np.stack(
        (v[..., 0] * c + v[..., 1] * s,
         -v[..., 0] * s + v[..., 1] * c,
         v[..., 2]),
        axis=-1,
    )


def test_provider_metadata_is_directly_typed_for_a3_and_keeps_truth_boundary() -> None:
    provider = _provider()
    metadata = provider.metadata
    assert isinstance(metadata, AuxiliaryT2ProviderMetadata)
    assert metadata.provider_kind == "repository_candidate"
    assert metadata.source_auxiliary_t2_field_materialized is True
    assert metadata.normalized_haar_measure_total_mass_one is True
    assert metadata.recomputed_from_actual_candidate is True
    assert metadata.surrogate_or_preaveraged_covariance_used is False
    assert metadata.heldout_data_used is False
    assert metadata.residual_as_forcing_used is False
    assert len(metadata.provider_semantic_sha256) == 64
    assert len(metadata.source_blob_sha1) == 40

    truth = truth_boundary()
    assert truth["same_identity_current_i4_candidate_bound"] is True
    assert truth["actual_sign_resolved_complete_curl_runtime_reused"] is True
    assert truth["raw_auxiliary_t2_wave_materialized"] is True
    assert truth["preaveraged_covariance_supplied"] is False
    assert truth["physical_theta_mean_relabelled_as_auxiliary_haar"] is False
    assert truth["auxiliary_sign_phase_lift_repository_autonomous"] is True
    assert truth["source_exact_auxiliary_t2_mode_family_recovered"] is False
    assert truth["source_exact_auxiliary_phase_assignment_recovered"] is False
    assert truth["source_exact_slow_derivatives_claimed"] is False
    assert truth["agent3_rf30_defect_materialized_here"] is False
    assert truth["agent3_mean_correction_materialized_here"] is False
    assert truth["forcing_or_pressure_added"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False


def test_auxiliary_t2_diagonal_exactly_replays_actual_sign_summed_candidate() -> None:
    provider = _provider()
    chart = provider.chart
    y = np.arange(8, dtype=float) / 8.0
    diagonal = np.broadcast_to(y[None, :], (8, 8))
    R = np.asarray(chart.radius_R, dtype=float)
    lifted = provider._raw_wave(R, chart.Z, chart.T, diagonal, diagonal)

    r = np.sqrt(chart.Q) * R
    z = (chart.Q ** chart.D) * chart.Z
    t = 1.0 - chart.Q * chart.T
    physical = _cylindrical_total(
        default_field(), r, z, t, 2.0 * np.pi * diagonal
    )
    expected = (chart.Q ** chart.A) * physical
    np.testing.assert_allclose(lifted, expected, rtol=0.0, atol=3.0e-12)


def test_raw_samples_and_fd4_slow_derivatives_are_finite_and_shape_exact() -> None:
    provider = _provider()
    chart = provider.chart
    y = np.arange(6, dtype=float) / 6.0
    y1, y2 = np.meshgrid(y, y, indexing="ij")
    sample = provider.sample_auxiliary_wave(
        R=np.asarray(chart.radius_R), Z=chart.Z, T=chart.T, y1=y1, y2=y2
    )
    assert isinstance(sample, AuxiliaryT2WaveSamples)
    expected = (len(chart.radius_R), 6, 6)
    for name in (
        "w_r", "w_theta", "w_z",
        "dR_w_r", "dR_w_theta", "dR_w_z",
        "dZ_w_r", "dZ_w_theta", "dZ_w_z",
    ):
        value = np.asarray(getattr(sample, name), dtype=float)
        assert value.shape == expected
        assert np.all(np.isfinite(value))
    wave = np.stack((sample.w_r, sample.w_theta, sample.w_z), axis=-1)
    assert float(np.sqrt(np.mean(wave * wave))) > 1.0e-12


def test_a3_haar_bridge_accepts_provider_but_stays_unpinned_and_unauthorized() -> None:
    backend = _backend()
    radius, z, t = _inputs()
    provider = bind_current_i4_auxiliary_t2_provider(backend, radius, z, t)
    receipt = materialize_current_i4_rf30_auxiliary_haar_bridge(
        backend, provider, radius, z, t
    )
    assert receipt.provider_kind == "repository_candidate"
    assert receipt.provider_semantic_sha256 == provider.semantic_sha256
    assert receipt.provider_source_blob_sha1 == provider.metadata.source_blob_sha1
    assert receipt.normalized_source_auxiliary_t2_haar_mean_used is True
    assert receipt.covariance_formed_in_bridge_from_raw_wave_samples is True
    assert receipt.covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives is True
    assert receipt.provider_checksum_pinned_for_repository_candidate is False
    assert receipt.rf30_repository_candidate_state_authorized is False
    assert np.all(np.isfinite(receipt.mean_W_rr_auxiliary_haar))
    assert np.all(np.isfinite(receipt.mean_W_zz_auxiliary_haar))


def test_provider_fails_closed_on_chart_or_torus_drift() -> None:
    provider = _provider()
    chart = provider.chart
    y = np.arange(4, dtype=float) / 4.0
    y1, y2 = np.meshgrid(y, y, indexing="ij")
    with pytest.raises(CurrentI4AuxiliaryT2ProviderError, match="R drifted"):
        provider.sample_auxiliary_wave(
            R=np.asarray(chart.radius_R) * 1.01,
            Z=chart.Z,
            T=chart.T,
            y1=y1,
            y2=y2,
        )
    bad = y1.copy()
    bad[0, 0] = 1.0
    with pytest.raises(CurrentI4AuxiliaryT2ProviderError, match="\[0,1\)"):
        provider.sample_auxiliary_wave(
            R=np.asarray(chart.radius_R), Z=chart.Z, T=chart.T, y1=bad, y2=y2
        )


def test_provider_requires_strict_i4_preflight() -> None:
    backend = _backend()
    preflight = materialize_current_i4_rf30_fixedq_preflight(
        backend, np.asarray((0.70, 0.75, 0.80)), 0.0, 0.5
    )
    assert preflight.all_strict_I4 is False
    with pytest.raises(CurrentI4AuxiliaryT2ProviderError, match="strict current-I4"):
        KokunoCurrentI4AuxiliaryT2Provider(backend, preflight)
