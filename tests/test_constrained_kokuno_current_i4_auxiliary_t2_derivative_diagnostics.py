from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i4_auxiliary_t2_derivative_diagnostics import (
    AUDIT_RELATIVE_STEPS,
    AUDIT_TORUS_ORDER,
    AuxiliaryT2DerivativeDiagnosticError,
    materialize_current_i4_auxiliary_t2_derivative_diagnostics,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_current_i4_auxiliary_t2_provider import (
    FD4_RELATIVE_STEP,
    PARENT_AGENT3_HEAD,
    KokunoCurrentI4AuxiliaryT2Provider,
    _runtime_payload_sha,
    bind_current_i4_auxiliary_t2_provider,
)
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import default_field


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
    return SimpleNamespace(
        velocity=np.zeros_like(p),
        velocity_jacobian=np.zeros(p.shape[:-1] + (3, 3), dtype=float),
        interior_mask=np.ones(p.shape[:-1], dtype=bool),
    )


def _backend() -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=_FakeComposite(),
        differential_function=_fake_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256=_runtime_payload_sha(default_field()),
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _provider() -> KokunoCurrentI4AuxiliaryT2Provider:
    return bind_current_i4_auxiliary_t2_provider(
        _backend(), np.asarray((0.40, 0.45, 0.50), dtype=float), 0.0, 0.5
    )


def test_three_resolution_diagnostic_is_identity_bound_and_truthful() -> None:
    receipt = materialize_current_i4_auxiliary_t2_derivative_diagnostics(_provider())
    assert receipt.parent_agent2_pr == 1198
    assert receipt.parent_agent2_head == "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
    assert receipt.parent_provider_source_blob_sha1 == "96168ac6583ad5e71aa044bd558c6feeadf051f4"
    assert receipt.provider_source_blob_sha1 == receipt.parent_provider_source_blob_sha1
    assert receipt.relative_steps == AUDIT_RELATIVE_STEPS
    assert receipt.relative_steps[-1] == FD4_RELATIVE_STEP
    assert receipt.torus_order == AUDIT_TORUS_ORDER
    assert receipt.raw_wave_rms > 1.0e-12
    assert len(receipt.dR_rms_by_step) == 3
    assert len(receipt.dZ_rms_by_step) == 3
    assert len(receipt.dR_successive_relative_rms_difference) == 2
    assert len(receipt.dZ_successive_relative_rms_difference) == 2
    assert np.all(np.isfinite(receipt.dR_rms_by_step))
    assert np.all(np.isfinite(receipt.dZ_rms_by_step))
    assert np.all(np.isfinite(receipt.dR_successive_relative_rms_difference))
    assert np.all(np.isfinite(receipt.dZ_successive_relative_rms_difference))

    truth = truth_boundary()
    assert truth["three_resolution_slow_derivative_diagnostic_materialized"] is True
    assert truth["audit_derivatives_recomputed_from_raw_wave"] is True
    assert truth["production_private_fd4_helper_used_by_audit_operator"] is False
    assert truth["offgrid_auxiliary_phase_sampling_used"] is True
    assert truth["production_fd4_finest_step_replayed"] is True
    assert truth["numerical_stability_observations_are_pde_acceptance_gates"] is False
    assert truth["source_exact_slow_derivatives_claimed"] is False
    assert truth["source_exact_auxiliary_phase_assignment_recovered"] is False
    assert truth["agent3_provider_pin_modified_here"] is False
    assert truth["agent3_rf30_defect_materialized_here"] is False
    assert truth["agent3_mean_correction_materialized_here"] is False
    assert truth["forcing_or_pressure_added"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False


def test_independent_finest_fd4_replays_production_derivative_operator() -> None:
    receipt = materialize_current_i4_auxiliary_t2_derivative_diagnostics(_provider())
    assert receipt.production_derivative_operator_consistent is True
    assert receipt.production_dR_replay_max_abs <= 5.0e-11
    assert receipt.production_dZ_replay_max_abs <= 5.0e-11
    assert receipt.production_dR_replay_relative_rms <= 5.0e-11
    assert receipt.production_dZ_replay_relative_rms <= 5.0e-11


def test_receipt_is_deterministic_for_same_exact_provider() -> None:
    provider = _provider()
    first = materialize_current_i4_auxiliary_t2_derivative_diagnostics(provider)
    second = materialize_current_i4_auxiliary_t2_derivative_diagnostics(provider)
    assert first.receipt_sha256 == second.receipt_sha256
    assert first.to_dict() == second.to_dict()


def test_corrupted_public_derivative_is_observed_as_operator_mismatch(monkeypatch) -> None:
    provider = _provider()
    original = provider.sample_auxiliary_wave

    def corrupted(**kwargs):
        sample = original(**kwargs)
        bad = 1.1 * np.asarray(sample.dR_w_r, dtype=float) + 1.0e-3
        return replace(sample, dR_w_r=bad)

    monkeypatch.setattr(provider, "sample_auxiliary_wave", corrupted)
    receipt = materialize_current_i4_auxiliary_t2_derivative_diagnostics(provider)
    assert receipt.production_derivative_operator_consistent is False
    assert receipt.production_dR_replay_max_abs > 1.0e-6


def test_diagnostic_rejects_nonexact_provider_type() -> None:
    with pytest.raises(TypeError, match="KokunoCurrentI4AuxiliaryT2Provider"):
        materialize_current_i4_auxiliary_t2_derivative_diagnostics(object())


def test_public_surface_exposes_no_scientific_or_fd_tuning_knobs() -> None:
    truth = truth_boundary()
    assert truth["public_parameters"] == ("provider",)
    assert truth["forbidden_scientific_controls_exposed"] is False
    # Parent remains tied to A3's exact current-I4 consumer identity; this audit
    # does not rewrite or authorize that lane.
    assert PARENT_AGENT3_HEAD == "21b0ade966fda3b7344273c7bff2efe7e7b9a93e"
