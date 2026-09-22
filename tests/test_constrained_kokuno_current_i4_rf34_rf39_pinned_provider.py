from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import pytest

import openai_ns_reconstruction.kokuno_current_i4_rf34_rf39_pinned_provider as module


CANDIDATE_ID = "current-I4-A2-1080-A1-1079"
CANDIDATE_SHA = "2" * 64
CHART_ID = "current-I4-fixed-Q-aaaaaaaaaaaaaaaa"
CHART_SHA = "a" * 64
TYPED_SHA = "b" * 64
PINNED_SHA = "c" * 64
PARENT_CORRECTION_SHA = "d" * 64
SYSTEM_SHA = "e" * 64


def _pinned():
    return SimpleNamespace(
        candidate_id=CANDIDATE_ID,
        candidate_sha256=CANDIDATE_SHA,
        source_chart_id=CHART_ID,
        source_chart_sha256=CHART_SHA,
        provider_semantic_sha256="3" * 64,
        provider_source_blob_sha1=module.PINNED_AGENT2_PROVIDER_BLOB,
        parent_typed_rf30_receipt_sha256=TYPED_SHA,
        receipt_sha256=PINNED_SHA,
        exact_provider_blob_pinned=True,
        repository_candidate_rf30_defect_evidence=True,
    )


def _compact_correction():
    return SimpleNamespace(
        parent_system_sha256=SYSTEM_SHA,
        source_chart_sha256=CHART_SHA,
        patch_context_sha256="f" * 64,
        autonomous_bump_family_sha256="1" * 64,
        coefficients_u0_u1_u2_s0_s1=(1.0, -0.4, 0.2, 0.08, -0.03),
        delta_v_chart=(0.0, 0.01, -0.02, 0.0),
        gamma_d_chart=(0.0, -0.004, 0.003, 0.0),
        correction_chart_l2=0.023,
        correction_chart_max=0.02,
        five_row_closure_max_relative=5.0e-10,
        coefficient_direct_solve_max_relative=6.0e-10,
        source_rf34_rf39_correction_materialized=True,
        compact_support_preserved=True,
        two_zero_moments_preserved=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
    )


def _parent():
    return SimpleNamespace(
        candidate=SimpleNamespace(
            candidate_id=CANDIDATE_ID,
            candidate_sha256=CANDIDATE_SHA,
            evidence_kind="mechanics-only",
        ),
        typed_rf30_receipt_sha256=TYPED_SHA,
        source_chart_id=CHART_ID,
        source_chart_sha256=CHART_SHA,
        correction=_compact_correction(),
        basis_frozen_before_defect_evaluation=True,
        patch_context_frozen_before_defect_evaluation=True,
        current_i4_rf31_system_materialized=True,
        current_i4_rf34_rf39_correction_materialized=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256=PARENT_CORRECTION_SHA,
    )


def _install(monkeypatch, *, pinned=None, parent=None):
    pinned = _pinned() if pinned is None else pinned
    parent = _parent() if parent is None else parent
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf30_pinned_provider_evidence",
        lambda *args, **kwargs: pinned,
    )
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_correction",
        lambda *args, **kwargs: parent,
    )
    return pinned, parent


def test_exact_pinned_rf30_receipt_promotes_existing_compact_correction(monkeypatch) -> None:
    pinned, parent = _install(monkeypatch)
    result = module.materialize_current_i4_rf34_rf39_pinned_provider_correction(
        object(), object(), (0.8, 0.9, 1.0), 0.0, 0.5
    )

    assert result.candidate_id == CANDIDATE_ID
    assert result.candidate_sha256 == CANDIDATE_SHA
    assert result.source_chart_sha256 == CHART_SHA
    assert result.provider_source_blob_sha1 == module.PINNED_AGENT2_PROVIDER_BLOB
    assert result.pinned_rf30_receipt_sha256 == pinned.receipt_sha256
    assert result.shared_typed_rf30_receipt_sha256 == TYPED_SHA
    assert result.parent_rf34_rf39_receipt_sha256 == parent.receipt_sha256
    assert result.rf31_system_sha256 == SYSTEM_SHA
    assert len(result.correction_payload_sha256) == 64
    assert len(result.receipt_sha256) == 64

    assert result.exact_provider_blob_pinned is True
    assert result.same_typed_rf30_receipt_preserved is True
    assert result.same_candidate_identity_preserved is True
    assert result.same_source_chart_identity_preserved is True
    assert result.basis_frozen_before_defect_evaluation is True
    assert result.patch_context_frozen_before_defect_evaluation is True
    assert result.compact_support_preserved is True
    assert result.two_zero_moments_preserved is True
    assert result.repository_candidate_rf31_system_evidence is True
    assert result.repository_candidate_rf34_rf39_correction_evidence is True

    assert result.coefficients_u0_u1_u2_s0_s1 == (1.0, -0.4, 0.2, 0.08, -0.03)
    assert result.correction_chart_l2 == 0.023
    assert result.correction_chart_max == 0.02
    assert result.five_row_closure_max_relative == 5.0e-10
    assert result.coefficient_direct_solve_max_relative == 6.0e-10

    assert result.correction_applied_to_candidate is False
    assert result.rf44_rf49_nonlinear_remainder_recomputed is False
    assert result.cartesian_correction_velocity_materialized is False
    assert result.finite_correction_cycle_run is False
    assert result.heldout_ns_residual_assessed is False
    assert result.pde_validated is False


def test_typed_rf30_receipt_mismatch_fails_before_promotion(monkeypatch) -> None:
    parent = _parent()
    parent.typed_rf30_receipt_sha256 = "0" * 64
    _install(monkeypatch, parent=parent)

    with pytest.raises(
        module.CurrentI4PinnedRF34RF39Error,
        match="exact typed RF30 receipt",
    ):
        module.materialize_current_i4_rf34_rf39_pinned_provider_correction(
            object(), object(), (0.8, 0.9, 1.0), 0.0, 0.5
        )


def test_downstream_application_or_rf44_claim_fails_closed(monkeypatch) -> None:
    parent = _parent()
    parent.correction_applied_to_candidate = True
    _install(monkeypatch, parent=parent)
    with pytest.raises(
        module.CurrentI4PinnedRF34RF39Error,
        match="unexpectedly applies correction",
    ):
        module.materialize_current_i4_rf34_rf39_pinned_provider_correction(
            object(), object(), (0.8, 0.9, 1.0), 0.0, 0.5
        )

    parent = _parent()
    parent.correction.rf44_rf49_nonlinear_remainder_recomputed = True
    _install(monkeypatch, parent=parent)
    with pytest.raises(
        module.CurrentI4PinnedRF34RF39Error,
        match="RF44-RF49 recomputation",
    ):
        module.materialize_current_i4_rf34_rf39_pinned_provider_correction(
            object(), object(), (0.8, 0.9, 1.0), 0.0, 0.5
        )


def test_truth_boundary_promotes_only_correction_evidence_not_ns() -> None:
    truth = module.truth_boundary()
    assert truth["exact_agent2_provider_blob_pin_reused"] is True
    assert truth["same_typed_rf30_receipt_required_for_promotion"] is True
    assert truth["rf30_rf31_rf34_rf39_formulas_reused_without_duplication"] is True
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["repository_candidate_rf31_system_evidence_materializable"] is True
    assert truth["repository_candidate_rf34_rf39_correction_evidence_materializable"] is True
    assert truth["correction_applied_to_candidate"] is False
    assert truth["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert truth["cartesian_correction_velocity_materialized"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["pde_validated"] is False
