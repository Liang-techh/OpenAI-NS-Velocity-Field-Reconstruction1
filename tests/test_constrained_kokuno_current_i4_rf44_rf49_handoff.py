from __future__ import annotations

from dataclasses import asdict
from types import SimpleNamespace

import pytest

import openai_ns_reconstruction.kokuno_current_i4_rf44_rf49_handoff as module
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_rf30_rf31_typed_mean_correction import CandidateIdentity
from openai_ns_reconstruction.kokuno_rf34_rf39_compact_mean_correction import (
    RF34RF39CorrectionReceipt,
)


def _backend() -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=SimpleNamespace(leading_backend=object()),
        differential_function=lambda *args, **kwargs: None,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256="3" * 64,
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _identity() -> CandidateIdentity:
    return CandidateIdentity(
        candidate_id="current-I4-A2-1080-A1-1079",
        candidate_sha256="7" * 64,
        evidence_kind="mechanics-only",
    )


def _correction(identity: CandidateIdentity) -> RF34RF39CorrectionReceipt:
    return RF34RF39CorrectionReceipt(
        candidate=identity,
        parent_system_sha256="a" * 64,
        source_chart_id="chart-current-i4",
        source_chart_sha256="b" * 64,
        patch_context_sha256="c" * 64,
        autonomous_bump_family_sha256="d" * 64,
        bump_log_centers=(0.0, 1.0, 2.0),
        bump_log_halfwidth=0.1,
        profile_moments=(("m0", 1.0),),
        profile_targets_Pq_Jthetaq_Jzq=(1.0, 2.0, 3.0),
        coefficients_u0_u1_u2_s0_s1=(1.0, 2.0, 3.0, 4.0, 5.0),
        direct_linear_solve_coefficients=(1.0, 2.0, 3.0, 4.0, 5.0),
        delta_v_chart=(0.1, 0.2, 0.3),
        gamma_d_chart=(0.01, 0.02, 0.03),
        basis_match_max_abs=0.0,
        reserved_base_G_max_abs=0.0,
        reserved_base_V_max_relative_error=0.0,
        five_row_closure_max_relative=0.0,
        coefficient_direct_solve_max_relative=0.0,
        correction_chart_l2=0.1,
        correction_chart_max=0.2,
        repository_candidate_correction_evidence=False,
        compact_support_preserved=True,
        two_zero_moments_preserved=True,
        source_rf34_rf39_correction_materialized=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        complete_ns_defect=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
    )


class _Provider:
    metadata = SimpleNamespace(provider_semantic_sha256="9" * 64)

    def sample_auxiliary_wave(self, **kwargs):
        raise AssertionError("raw-wave path is monkeypatched in the handoff unit test")

    def rf44_pre_update_state(self, **kwargs):
        raise AssertionError("delegated RF44 method is consumed by the reused #1142 operator")

    def rf44_post_update_state(self, **kwargs):
        raise AssertionError("delegated RF44 method is consumed by the reused #1142 operator")

    def rf44_correction_increment(self, **kwargs):
        raise AssertionError("delegated RF44 method is consumed by the reused #1142 operator")


def _install_handoff_fixture(monkeypatch, *, parent_sha_mismatch: bool = False):
    identity = _identity()
    correction = _correction(identity)
    current = SimpleNamespace(
        candidate=identity,
        correction=correction,
        source_chart_id=correction.source_chart_id,
        source_chart_sha256=correction.source_chart_sha256,
        receipt_sha256="e" * 64,
        repository_candidate_correction_evidence=False,
    )
    preflight = object()
    context = object()
    basis = object()
    typed = SimpleNamespace(
        candidate=identity,
        state=object(),
        source_chart_id=correction.source_chart_id,
        source_chart_sha256=correction.source_chart_sha256,
        receipt_sha256="f" * 64,
    )

    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_correction",
        lambda *args, **kwargs: current,
    )
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf30_fixedq_preflight",
        lambda *args, **kwargs: preflight,
    )
    monkeypatch.setattr(
        module,
        "_freeze_patch_context_and_basis",
        lambda *args, **kwargs: (identity, context, basis),
    )
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf30_typed_defect",
        lambda *args, **kwargs: typed,
    )
    monkeypatch.setattr(
        module,
        "_bind_typed_state",
        lambda *args, **kwargs: (identity, basis, context),
    )

    expected_sha = module._rf44_sha(asdict(correction))
    observed_sha = "0" * 64 if parent_sha_mismatch else expected_sha
    recomputed = SimpleNamespace(
        parent_correction_sha256=observed_sha,
        candidate=identity,
        source_chart_id=correction.source_chart_id,
        source_chart_sha256=correction.source_chart_sha256,
        repository_candidate_recompute_evidence=False,
        pre_defect=SimpleNamespace(state_sha256="1" * 64),
        post_defect=SimpleNamespace(state_sha256="2" * 64),
        defect_gain_ratio_2=0.75,
    )
    monkeypatch.setattr(
        module,
        "materialize_rf44_rf49_postupdate_recompute",
        lambda *args, **kwargs: recomputed,
    )
    return current, recomputed


def test_routes_exact_current_i4_correction_into_existing_rf44_rf49(monkeypatch) -> None:
    current, recomputed = _install_handoff_fixture(monkeypatch)
    provider = _Provider()
    receipt = module.materialize_current_i4_rf44_rf49_handoff(
        _backend(), provider, (1.1, 1.2, 1.3), 0.0, 0.5
    )

    assert receipt.candidate == current.candidate
    assert receipt.current_i4_rf34_rf39_receipt_sha256 == current.receipt_sha256
    assert receipt.source_chart_id == current.source_chart_id
    assert receipt.source_chart_sha256 == current.source_chart_sha256
    assert receipt.provider_semantic_sha256 == "9" * 64
    assert receipt.rf44_rf49 is recomputed
    assert receipt.rf44_rf49.defect_gain_ratio_2 == 0.75
    assert receipt.exact_rf34_rf39_replay_verified is True
    assert receipt.rf44_state_provider_same_object_as_auxiliary_wave_provider is True
    assert receipt.rf44_rf49_operator_recompute_executed is True
    assert receipt.repository_provider_blob_pinned is False
    assert receipt.current_i4_repository_candidate_rf44_rf49_remainder_recomputed is False
    assert receipt.chart_correction_applied_in_recompute is True
    assert receipt.cartesian_correction_velocity_materialized is False
    assert receipt.finite_correction_cycle_run is False
    assert receipt.heldout_ns_residual_assessed is False
    assert receipt.residual_reduction_claimed is False
    assert receipt.pde_validated is False
    assert len(receipt.receipt_sha256) == 64


def test_parent_correction_sha_mismatch_fails_closed(monkeypatch) -> None:
    _install_handoff_fixture(monkeypatch, parent_sha_mismatch=True)
    with pytest.raises(
        module.CurrentI4RF44RF49HandoffError,
        match="exact #1190 RF34--RF39 correction",
    ):
        module.materialize_current_i4_rf44_rf49_handoff(
            _backend(), _Provider(), (1.1, 1.2, 1.3), 0.0, 0.5
        )


def test_missing_rf44_provider_method_fails_before_parent_execution(monkeypatch) -> None:
    called = []
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_correction",
        lambda *args, **kwargs: called.append("parent"),
    )
    with pytest.raises(
        module.CurrentI4RF44RF49HandoffError,
        match="rf44_pre_update_state",
    ):
        module.materialize_current_i4_rf44_rf49_handoff(
            _backend(), SimpleNamespace(metadata=SimpleNamespace()), (1.1,), 0.0, 0.5
        )
    assert called == []


def test_truth_boundary_keeps_recompute_distinct_from_scientific_admission() -> None:
    truth = module.truth_boundary()
    assert truth["current_i4_rf34_rf39_correction_consumed"] is True
    assert truth["existing_rf44_rf49_operator_reused_without_formula_duplication"] is True
    assert truth["same_provider_object_required_for_auxiliary_t2_and_rf44_states"] is True
    assert truth["caller_supplied_defect_allowed"] is False
    assert truth["caller_supplied_correction_allowed"] is False
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["repository_provider_blob_pinned_in_parent"] is False
    assert truth["rf44_rf49_operator_recompute_executable"] is True
    assert truth["current_i4_repository_candidate_rf44_rf49_remainder_recomputed"] is False
    assert truth["cartesian_correction_velocity_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["residual_defined_free_forcing_allowed"] is False
    assert truth["pde_validated"] is False
