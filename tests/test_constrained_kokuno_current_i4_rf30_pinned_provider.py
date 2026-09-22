from __future__ import annotations

from dataclasses import replace
import math
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i4_auxiliary_t2_provider import (
    _runtime_payload_sha,
    bind_current_i4_auxiliary_t2_provider,
)
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_pinned_provider import (
    PINNED_AGENT2_HEAD,
    PINNED_AGENT2_PR,
    PINNED_AGENT2_PROVIDER_BLOB,
    CurrentI4RF30PinnedProviderError,
    materialize_current_i4_rf30_pinned_provider_evidence,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import default_field


class _Schedule:
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


class _Leading:
    A = 0.505
    D = 0.495
    X_I4_start = 0.5
    X_I4_end = 1.5
    semantic_sha256 = "1" * 64
    outer_schedule = _Schedule()

    @staticmethod
    def similarity_coordinates(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        q = 1.0 - t
        eta = z / np.power(q, _Leading.D)
        X = (x * x + y * y) / (2.0 * q)
        return {"q": q, "eta": eta, "X": X}

    @staticmethod
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((0.0 * x, 0.0 * y, 0.0 * z), axis=-1)


class _Composite:
    leading_backend = _Leading()


def _fake_differential(points, time):
    p = np.asarray(points, dtype=float)
    return SimpleNamespace(
        velocity=np.zeros_like(p),
        velocity_jacobian=np.zeros(p.shape[:-1] + (3, 3), dtype=float),
        interior_mask=np.ones(p.shape[:-1], dtype=bool),
    )


def _backend() -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=_Composite(),
        differential_function=_fake_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256=_runtime_payload_sha(default_field()),
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _inputs():
    return np.asarray((0.80, 0.90, 1.00), dtype=float), 0.0, 0.5


class _MetadataDriftProvider:
    def __init__(self, provider):
        self._provider = provider
        self.metadata = replace(provider.metadata, source_blob_sha1="0" * 40)

    def sample_auxiliary_wave(self, **kwargs):
        return self._provider.sample_auxiliary_wave(**kwargs)


def test_exact_a2_provider_promotes_parent_rf30_defect_without_manual_defect_input() -> None:
    backend = _backend()
    radius, z, t = _inputs()
    provider = bind_current_i4_auxiliary_t2_provider(backend, radius, z, t)
    result = materialize_current_i4_rf30_pinned_provider_evidence(
        backend, provider, radius, z, t
    )

    assert PINNED_AGENT2_PR == 1198
    assert PINNED_AGENT2_HEAD == "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
    assert PINNED_AGENT2_PROVIDER_BLOB == "96168ac6583ad5e71aa044bd558c6feeadf051f4"
    assert provider.metadata.source_blob_sha1 == PINNED_AGENT2_PROVIDER_BLOB
    assert result.provider_source_blob_sha1 == PINNED_AGENT2_PROVIDER_BLOB
    assert result.exact_provider_blob_pinned is True
    assert result.same_candidate_identity_preserved is True
    assert result.actual_candidate_recomputed is True
    assert result.normalized_source_auxiliary_t2_haar_mean_used is True
    assert result.source_fixed_q_chart_used is True
    assert result.surrogate_defect_used is False
    assert result.residual_as_forcing_shortcut_used is False
    assert result.heldout_samples_used_to_construct_state is False
    assert result.repository_candidate_rf30_defect_evidence is True
    assert math.isfinite(result.defect_tuple_l2)
    assert result.defect_tuple_l2 >= 0.0
    assert len(result.receipt_sha256) == 64

    assert result.rf31_five_row_system_materialized is False
    assert result.rf34_rf39_correction_materialized is False
    assert result.rf44_rf49_remainder_recomputed is False
    assert result.cartesian_correction_velocity_materialized is False
    assert result.finite_correction_cycle_run is False
    assert result.heldout_ns_residual_assessed is False
    assert result.pde_validated is False


def test_provider_blob_drift_fails_before_rf30_promotion() -> None:
    backend = _backend()
    radius, z, t = _inputs()
    provider = bind_current_i4_auxiliary_t2_provider(backend, radius, z, t)
    with pytest.raises(CurrentI4RF30PinnedProviderError, match="exact pinned implementation"):
        materialize_current_i4_rf30_pinned_provider_evidence(
            backend, _MetadataDriftProvider(provider), radius, z, t
        )


def test_truth_boundary_stops_at_repository_candidate_rf30_evidence() -> None:
    truth = truth_boundary()
    assert truth["exact_provider_blob_pin_enforced"] is True
    assert truth["parent_haar_and_rf30_operators_reused_without_formula_duplication"] is True
    assert truth["repository_candidate_rf30_promotion_adapter_executable"] is True
    assert truth["source_exact_auxiliary_mode_recovery_claimed"] is False
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["rf31_five_row_system_materialized"] is False
    assert truth["rf34_rf39_correction_materialized"] is False
    assert truth["rf44_rf49_remainder_recomputed"] is False
    assert truth["cartesian_correction_velocity_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["pde_validated"] is False
