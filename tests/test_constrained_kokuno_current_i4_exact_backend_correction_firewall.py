from __future__ import annotations

from types import SimpleNamespace

import pytest

import openai_ns_reconstruction.kokuno_current_i4_exact_backend_correction_firewall as module
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)


def _fake_differential(points, time):
    return None


class _FakeComposite:
    pass


def _backend(*, composite_semantic="2" * 64) -> ExactCurrentI4NonlinearBackend:
    return ExactCurrentI4NonlinearBackend(
        composite_field=_FakeComposite(),
        differential_function=_fake_differential,
        composite_semantic_sha256=composite_semantic,
        oscillatory_runtime_sha256="3" * 64,
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _parent_receipt():
    return SimpleNamespace(
        candidate_id="current-I4-A2-1080-A1-1079",
        candidate_sha256="5" * 64,
        source_chart_id="current-I4-fixed-Q-aaaaaaaaaaaaaaaa",
        source_chart_sha256="6" * 64,
        correction_payload_sha256="7" * 64,
        coefficients_u0_u1_u2_s0_s1=(0.8, -0.3, 0.11, 0.05, -0.02),
        correction_chart_l2=0.019,
        correction_chart_max=0.014,
        five_row_closure_max_relative=7.0e-10,
        coefficient_direct_solve_max_relative=8.0e-10,
        exact_provider_blob_pinned=True,
        same_typed_rf30_receipt_preserved=True,
        same_candidate_identity_preserved=True,
        same_source_chart_identity_preserved=True,
        repository_candidate_rf31_system_evidence=True,
        repository_candidate_rf34_rf39_correction_evidence=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        finite_correction_cycle_run=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256="8" * 64,
    )


def test_directly_constructed_fixture_backend_fails_before_parent_correction(monkeypatch) -> None:
    called = {"parent": False}

    def _parent(*args, **kwargs):
        called["parent"] = True
        return _parent_receipt()

    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_pinned_provider_correction",
        _parent,
    )

    with pytest.raises(
        module.CurrentI4ExactBackendFirewallError,
        match="exact A2 #1080/#960 backend rebind failed",
    ):
        module.materialize_current_i4_rf34_rf39_exact_backend_correction(
            _backend(), object(), (0.8, 0.9, 1.0), 0.0, 0.5
        )
    assert called["parent"] is False


def test_exact_rebind_is_used_for_parent_and_promotes_only_scoped_correction(monkeypatch) -> None:
    original = _backend()
    rebound = _backend()
    calls = {"bound": 0, "parent_backend": None}

    def _bind(cls, composite_field, differential_function):
        calls["bound"] += 1
        assert composite_field is original.composite_field
        assert differential_function is original.differential_function
        return rebound

    def _parent(backend, provider, radius, z, t):
        calls["parent_backend"] = backend
        return _parent_receipt()

    monkeypatch.setattr(
        module.ExactCurrentI4NonlinearBackend,
        "bind",
        classmethod(_bind),
    )
    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_pinned_provider_correction",
        _parent,
    )

    result = module.materialize_current_i4_rf34_rf39_exact_backend_correction(
        original, object(), (0.8, 0.9, 1.0), 0.0, 0.5
    )

    assert calls["bound"] == 1
    assert calls["parent_backend"] is rebound
    assert result.exact_backend_rebind_executed is True
    assert result.concrete_a2_composite_identity_authenticated is True
    assert result.concrete_a2_differential_identity_authenticated is True
    assert result.caller_dataclass_identity_replayed_exactly is True
    assert result.pinned_auxiliary_t2_provider_consumed is True
    assert result.repository_candidate_rf30_defect_evidence is True
    assert result.repository_candidate_rf31_system_evidence is True
    assert result.repository_candidate_rf34_rf39_correction_evidence is True
    assert result.coefficients_u0_u1_u2_s0_s1 == (0.8, -0.3, 0.11, 0.05, -0.02)
    assert result.correction_chart_l2 == 0.019
    assert result.correction_applied_to_candidate is False
    assert result.rf44_rf49_nonlinear_remainder_recomputed is False
    assert result.cartesian_correction_velocity_materialized is False
    assert result.finite_correction_cycle_run is False
    assert result.heldout_ns_residual_assessed is False
    assert result.pde_validated is False
    assert len(result.exact_backend_receipt_sha256) == 64
    assert len(result.receipt_sha256) == 64


def test_rebound_identity_mismatch_fails_before_parent_correction(monkeypatch) -> None:
    original = _backend()
    rebound = _backend(composite_semantic="9" * 64)
    called = {"parent": False}

    monkeypatch.setattr(
        module.ExactCurrentI4NonlinearBackend,
        "bind",
        classmethod(lambda cls, composite_field, differential_function: rebound),
    )

    def _parent(*args, **kwargs):
        called["parent"] = True
        return _parent_receipt()

    monkeypatch.setattr(
        module,
        "materialize_current_i4_rf34_rf39_pinned_provider_correction",
        _parent,
    )

    with pytest.raises(
        module.CurrentI4ExactBackendFirewallError,
        match="identity disagrees",
    ):
        module.materialize_current_i4_rf34_rf39_exact_backend_correction(
            original, object(), (0.8, 0.9, 1.0), 0.0, 0.5
        )
    assert called["parent"] is False


def test_truth_boundary_rejects_fixture_promotion_and_keeps_ns_false() -> None:
    truth = module.truth_boundary()
    assert truth["parent_pr"] == 1210
    assert truth["parent_head"] == "e75940c33127b0725ddc703eaae97590158f2a56"
    assert truth["parent_blob"] == "b99b3bee2b2dad6003c11acae7946c9be53347fb"
    assert truth["exact_backend_rebind_required"] is True
    assert truth["exact_backend_bind_authenticates_concrete_module_class_function_and_blobs"] is True
    assert truth["manual_ExactCurrentI4NonlinearBackend_dataclass_construction_is_sufficient"] is False
    assert truth["zero_or_synthetic_backend_fixture_is_scientific_correction_evidence"] is False
    assert truth["parent_1210_promotion_without_exact_rebind_is_sufficient"] is False
    assert truth["parent_rf30_rf31_rf34_rf39_formulas_changed"] is False
    assert truth["agent2_curl_or_jacobian_reimplemented"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["correction_applied_to_candidate"] is False
    assert truth["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["pde_validated"] is False
