from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_i4_rf34_rf39_correction as module
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30PreMeanState,
)


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
    outer_schedule = _Schedule()


class _Composite:
    leading_backend = _Leading()


def _unused_differential(points, time):
    raise AssertionError("A2 differential path must not be reimplemented/called in adapter test")


def _backend() -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=_Composite(),
        differential_function=_unused_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256="3" * 64,
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _preflight():
    R = np.linspace(1.001, math.sqrt(3.0) - 0.001, 4001)
    X = 0.5 * R * R
    return SimpleNamespace(
        candidate_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256="3" * 64,
        leading_semantic_sha256="1" * 64,
        preflight_sha256="a" * 64,
        Q=0.5,
        s_Q=1.0,
        A=0.505,
        D=0.495,
        eta=0.0,
        Z=0.0,
        T=1.0,
        radius_R=tuple(float(v) for v in R),
        X=tuple(float(v) for v in X),
        all_strict_I4=True,
    )


def _typed(preflight, *, repository: bool = False):
    chart_id, chart_sha = module._chart_identity(preflight)
    R = np.asarray(preflight.radius_R, dtype=float)
    lam = _Schedule.lambda_outer
    a = (2.0 ** (0.5 + lam)) * math.exp(_Schedule.log_c_patch)
    V = a * R ** (-1.0 - 2.0 * lam)
    bump = np.exp(-((R - 1.34) / 0.11) ** 2)
    zeros = np.zeros_like(R)
    state = RF30PreMeanState(
        source_candidate_id=f"current-I4-A2-{module.AGENT2_COMPOSITE_PR}-A1-{module.AGENT1_LEADING_PR}",
        source_candidate_sha256=preflight.candidate_semantic_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        radius_R=tuple(float(v) for v in R),
        base_V=tuple(float(v) for v in V),
        base_G=tuple(float(v) for v in zeros),
        mean_W_rr=tuple(float(v) for v in zeros),
        mean_W_zr=tuple(float(v) for v in zeros),
        mean_W_thetatheta=tuple(float(v) for v in (0.04 * R * bump)),
        mean_W_ztheta=tuple(float(v) for v in (-0.025 * bump)),
        mean_W_zz=tuple(float(v) for v in (0.03 * bump)),
        dR_mean_W_rr=tuple(float(v) for v in zeros),
        dZ_mean_W_zr=tuple(float(v) for v in zeros),
        actual_candidate_recomputed=repository,
        oscillatory_covariance_recomputed=True,
        normalized_haar_mean_used=True,
        source_fixed_q_chart_used=True,
        surrogate_defect_used=False,
        residual_as_forcing_shortcut_used=False,
        heldout_samples_used_to_construct_state=False,
    )
    identity = CandidateIdentity(
        candidate_id=state.source_candidate_id,
        candidate_sha256=state.source_candidate_sha256,
        evidence_kind="repository-candidate" if repository else "mechanics-only",
    )
    return SimpleNamespace(
        candidate=identity,
        state=state,
        fixed_q_preflight_sha256=preflight.preflight_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        receipt_sha256="c" * 64,
        repository_candidate_defect_evidence=repository,
    )


def _install_fixture(monkeypatch, *, repository: bool = False, events=None):
    preflight = _preflight()

    def fake_preflight(backend, radius, z, t):
        return preflight

    def fake_typed(backend, provider, radius, z, t):
        if events is not None:
            events.append("defect")
        return _typed(preflight, repository=repository)

    monkeypatch.setattr(
        module, "materialize_current_i4_rf30_fixedq_preflight", fake_preflight
    )
    monkeypatch.setattr(
        module, "materialize_current_i4_rf30_typed_defect", fake_typed
    )
    return preflight


def test_freezes_source_i4_basis_before_defect_and_reuses_rf34_rf39(monkeypatch) -> None:
    events = []
    preflight = _install_fixture(monkeypatch, events=events)
    original_basis = module._canonical_chart_basis

    def recorded_basis(context, radius_R):
        events.append("basis")
        return original_basis(context, radius_R)

    monkeypatch.setattr(module, "_canonical_chart_basis", recorded_basis)
    result = module.materialize_current_i4_rf34_rf39_correction(
        _backend(), object(), np.asarray((1.1, 1.2, 1.3)), 0.0, 0.5
    )

    assert events[:2] == ["basis", "defect"]
    assert result.basis_frozen_before_defect_evaluation is True
    assert result.patch_context_frozen_before_defect_evaluation is True
    assert result.current_i4_rf31_system_materialized is True
    assert result.current_i4_rf34_rf39_correction_materialized is True
    assert result.correction.source_rf34_rf39_correction_materialized is True
    assert result.correction.compact_support_preserved is True
    assert result.correction.two_zero_moments_preserved is True
    assert result.correction.five_row_closure_max_relative <= 2.0e-7
    assert result.correction.coefficient_direct_solve_max_relative <= 2.0e-7
    assert result.correction.correction_chart_max > 0.0
    assert result.repository_provider_blob_pinned is False
    assert result.repository_candidate_correction_evidence is False
    assert result.correction_applied_to_candidate is False
    assert result.rf44_rf49_nonlinear_remainder_recomputed is False
    assert result.cartesian_correction_velocity_materialized is False
    assert result.heldout_ns_residual_assessed is False
    assert result.pde_validated is False

    assert math.isclose(result.patch_context.mean_patch_x_min, 1.0)
    assert math.isclose(result.patch_context.mean_patch_x_max, math.sqrt(3.0))
    assert result.patch_context.source_exact_bump_claimed is False
    assert result.patch_context.heldout_samples_used_to_choose_context is False
    assert result.source_chart_sha256 == module._chart_identity(preflight)[1]
    assert len(result.receipt_sha256) == 64


def test_unpinned_parent_cannot_promote_repository_candidate_correction(monkeypatch) -> None:
    _install_fixture(monkeypatch, repository=True)
    with pytest.raises(
        module.CurrentI4RF34RF39CorrectionError,
        match="provider blob is unpinned",
    ):
        module.materialize_current_i4_rf34_rf39_correction(
            _backend(), object(), (1.1, 1.2, 1.3), 0.0, 0.5
        )


def test_typed_chart_drift_fails_closed(monkeypatch) -> None:
    preflight = _preflight()

    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf30_fixedq_preflight",
        lambda *args, **kwargs: preflight,
    )
    typed = _typed(preflight)
    typed.source_chart_sha256 = "d" * 64
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf30_typed_defect",
        lambda *args, **kwargs: typed,
    )

    with pytest.raises(
        module.CurrentI4RF34RF39CorrectionError,
        match="source chart hash drifted",
    ):
        module.materialize_current_i4_rf34_rf39_correction(
            _backend(), object(), (1.1, 1.2, 1.3), 0.0, 0.5
        )


def test_truth_boundary_keeps_compact_correction_separate_from_application_and_ns() -> None:
    truth = module.truth_boundary()
    assert truth["current_i4_rf30_typed_defect_consumed"] is True
    assert truth["rf31_basis_derived_from_source_i4_schedule_before_defect"] is True
    assert truth["rf34_rf39_source_specific_compact_realization_reused"] is True
    assert truth["caller_supplied_rf31_basis_allowed"] is False
    assert truth["caller_supplied_defect_allowed"] is False
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["repository_provider_blob_pinned_in_parent"] is False
    assert truth["current_i4_repository_candidate_rf34_rf39_correction_materialized"] is False
    assert truth["correction_applied_to_candidate"] is False
    assert truth["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert truth["cartesian_correction_velocity_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
