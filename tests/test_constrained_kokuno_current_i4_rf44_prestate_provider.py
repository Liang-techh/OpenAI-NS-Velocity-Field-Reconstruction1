from __future__ import annotations

from dataclasses import replace
import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_i4_rf44_prestate_provider as module
from openai_ns_reconstruction.kokuno_current_i4_auxiliary_t2_provider import (
    _runtime_payload_sha,
)
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from openai_ns_reconstruction.kokuno_current_i4_rf30_typed_defect import (
    CurrentI4RF30TypedDefectReceipt,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import default_field
from openai_ns_reconstruction.kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30DefectTuple,
    RF30PreMeanState,
)
from openai_ns_reconstruction.kokuno_rf34_rf39_compact_mean_correction import (
    RF34RF39CorrectionReceipt,
)
from openai_ns_reconstruction.kokuno_rf44_rf49_postupdate_recompute import RF44MeanState


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
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        # A simple nonzero radial physical leading field.  On theta=0,
        # u_r=u_x=2r, which makes the Q^A normalization independently checkable.
        return np.stack((2.0 * x, 2.0 * y, 0.0 * z), axis=-1)


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
    return ExactCurrentI4NonlinearBackend(
        composite_field=_FakeComposite(),
        differential_function=_fake_differential,
        composite_semantic_sha256="2" * 64,
        oscillatory_runtime_sha256=_runtime_payload_sha(default_field()),
        differential_semantic_sha256="4" * 64,
        composite_source_blob=AGENT2_COMPOSITE_SOURCE_BLOB,
        differential_source_blob=AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    )


def _inputs():
    return np.asarray((0.40, 0.45, 0.50), dtype=float), 0.0, 0.5


def _provider() -> module.KokunoCurrentI4RF44PreStateProvider:
    radius, z, t = _inputs()
    return module.bind_current_i4_rf44_prestate_provider(_backend(), radius, z, t)


def _identity(provider) -> CandidateIdentity:
    return CandidateIdentity(
        candidate_id="current-I4-A2-1080-A1-1079",
        candidate_sha256=provider._preflight.candidate_semantic_sha256,
        evidence_kind="repository-candidate",
    )


def _correction(identity: CandidateIdentity, chart_id: str, chart_sha: str) -> RF34RF39CorrectionReceipt:
    return RF34RF39CorrectionReceipt(
        candidate=identity,
        parent_system_sha256="a" * 64,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
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
        repository_candidate_correction_evidence=True,
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


def _typed_receipt(provider) -> CurrentI4RF30TypedDefectReceipt:
    identity = _identity(provider)
    R = np.asarray(provider.chart.radius_R, dtype=float)
    chart_id = "current-I4-fixed-Q-test"
    chart_sha = "b" * 64
    Wrr = 0.4 + 0.1 * R
    Wzr = -0.03 + 0.02 * R
    Wtt = 0.2 + 0.05 * R
    Wzt = 0.07 - 0.01 * R
    Wzz = 0.3 + 0.04 * R
    dR = np.full_like(R, 0.1)
    dZ = np.full_like(R, -0.025)
    V = 1.0 + 0.2 * R
    G = np.zeros_like(R)
    state = RF30PreMeanState(
        source_candidate_id=identity.candidate_id,
        source_candidate_sha256=identity.candidate_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        radius_R=tuple(float(v) for v in R),
        base_V=tuple(float(v) for v in V),
        base_G=tuple(float(v) for v in G),
        mean_W_rr=tuple(float(v) for v in Wrr),
        mean_W_zr=tuple(float(v) for v in Wzr),
        mean_W_thetatheta=tuple(float(v) for v in Wtt),
        mean_W_ztheta=tuple(float(v) for v in Wzt),
        mean_W_zz=tuple(float(v) for v in Wzz),
        dR_mean_W_rr=tuple(float(v) for v in dR),
        dZ_mean_W_zr=tuple(float(v) for v in dZ),
        actual_candidate_recomputed=True,
        oscillatory_covariance_recomputed=True,
        normalized_haar_mean_used=True,
        source_fixed_q_chart_used=True,
        surrogate_defect_used=False,
        residual_as_forcing_shortcut_used=False,
        heldout_samples_used_to_construct_state=False,
    )
    defect = RF30DefectTuple(
        P=0.0,
        J_theta=0.0,
        J_z=0.0,
        target_rows=(0.0, 0.0, 0.0, 0.0, 0.0),
        radial_source_g_r=tuple(0.0 for _ in R),
        state_sha256="e" * 64,
    )
    return CurrentI4RF30TypedDefectReceipt(
        candidate=identity,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        fixed_q_preflight_sha256=provider._preflight.preflight_sha256,
        auxiliary_haar_receipt_sha256="a" * 64,
        provider_semantic_sha256=provider.metadata.provider_semantic_sha256,
        provider_source_blob_sha1=provider.metadata.source_blob_sha1,
        provider_kind="repository_candidate",
        state=state,
        defect=defect,
        lambda_value=0.05,
        h_value=0.005,
        log_c_patch=0.0,
        base_V_min_abs=float(np.min(np.abs(V))),
        base_V_max_abs=float(np.max(np.abs(V))),
        base_G_max_abs=0.0,
        source_background_from_exact_a1_i4_schedule=True,
        normalized_source_auxiliary_t2_haar_mean_used=True,
        typed_rf30_state_materialized=True,
        rf30_defect_materialized=True,
        repository_candidate_defect_evidence=True,
        rf31_five_row_system_materialized=False,
        correction_applied=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256="f" * 64,
    )


def _handoff(provider, typed):
    correction = _correction(typed.candidate, typed.source_chart_id, typed.source_chart_sha256)
    current = SimpleNamespace(
        candidate=typed.candidate,
        source_chart_id=typed.source_chart_id,
        source_chart_sha256=typed.source_chart_sha256,
    )
    return current, correction


def test_raw_provider_identity_and_samples_are_delegated_unchanged() -> None:
    provider = _provider()
    raw = provider._raw
    assert provider.metadata == raw.metadata
    assert provider.semantic_sha256 == raw.semantic_sha256
    assert provider.metadata.source_blob_sha1 == module.PARENT_AGENT2_RAW_PROVIDER_BLOB

    y = np.arange(4, dtype=float) / 4.0
    y1, y2 = np.meshgrid(y, y, indexing="ij")
    kwargs = dict(
        R=np.asarray(provider.chart.radius_R),
        Z=provider.chart.Z,
        T=provider.chart.T,
        y1=y1,
        y2=y2,
    )
    delegated = provider.sample_auxiliary_wave(**kwargs)
    direct = raw.sample_auxiliary_wave(**kwargs)
    for name in (
        "w_r", "w_theta", "w_z",
        "dR_w_r", "dR_w_theta", "dR_w_z",
        "dZ_w_r", "dZ_w_theta", "dZ_w_z",
    ):
        np.testing.assert_array_equal(getattr(delegated, name), getattr(direct, name))


def test_prestate_uses_actual_leading_b_and_rf44_reduces_to_rf30(monkeypatch) -> None:
    provider = _provider()
    typed = _typed_receipt(provider)
    monkeypatch.setattr(module, "materialize_current_i4_rf30_typed_defect", lambda *a, **k: typed)
    current, correction = _handoff(provider, typed)

    pre = provider.rf44_pre_update_state(current=current, correction=correction)
    assert isinstance(pre, RF44MeanState)
    p = provider._preflight
    radius_phys = np.asarray(p.radius_physical, dtype=float)
    expected_b = (p.Q ** p.A) * (2.0 * radius_phys)
    np.testing.assert_allclose(pre.base_b, expected_b, rtol=0.0, atol=2.0e-15)

    n = len(pre.radius_R)
    zeros = np.zeros(n)
    for name in ("beta", "gamma", "v", "tstar_beta", "epsilon_laplacian_beta"):
        np.testing.assert_array_equal(np.asarray(getattr(pre, name), dtype=float), zeros)
    assert pre.dR_rr_bundle == typed.state.dR_mean_W_rr
    assert pre.dZ_zr_bundle == typed.state.dZ_mean_W_zr
    assert pre.mean_W_rr == typed.state.mean_W_rr
    assert pre.mean_W_zr == typed.state.mean_W_zr
    assert pre.mean_W_thetatheta == typed.state.mean_W_thetatheta
    assert pre.mean_W_ztheta == typed.state.mean_W_ztheta
    assert pre.mean_W_zz == typed.state.mean_W_zz
    assert pre.actual_candidate_recomputed is True
    assert pre.normalized_haar_mean_used is True
    assert pre.source_fixed_q_chart_used is True
    assert pre.operator_terms_recomputed_from_state is True
    assert pre.correction_applied is False

    # Implementation-distinct algebra check: specialize RF44 manually and
    # compare it with the RF30 pre-mean expression, rather than calling either
    # production defect helper.
    R = np.asarray(pre.radius_R, dtype=float)
    b = np.asarray(pre.base_b, dtype=float)
    V = np.asarray(pre.base_V, dtype=float)
    beta = np.asarray(pre.beta, dtype=float)
    gamma = np.asarray(pre.gamma, dtype=float)
    v = np.asarray(pre.v, dtype=float)
    Wrr = np.asarray(pre.mean_W_rr, dtype=float)
    Wzr = np.asarray(pre.mean_W_zr, dtype=float)
    Wtt = np.asarray(pre.mean_W_thetatheta, dtype=float)
    rf44 = (
        -np.asarray(pre.tstar_beta, dtype=float)
        -np.asarray(pre.dR_rr_bundle, dtype=float)
        -(2.0 * b * beta + beta * beta + Wrr) / R
        -np.asarray(pre.dZ_zr_bundle, dtype=float)
        +(2.0 * V * v + v * v + Wtt) / R
        +np.asarray(pre.epsilon_laplacian_beta, dtype=float)
    )
    rf30 = (
        -np.asarray(typed.state.dR_mean_W_rr, dtype=float)
        -Wrr / R
        -np.asarray(typed.state.dZ_mean_W_zr, dtype=float)
        +Wtt / R
    )
    np.testing.assert_allclose(rf44, rf30, rtol=0.0, atol=0.0)


def test_prestate_does_not_consume_mean_correction_coefficients(monkeypatch) -> None:
    provider = _provider()
    typed = _typed_receipt(provider)
    monkeypatch.setattr(module, "materialize_current_i4_rf30_typed_defect", lambda *a, **k: typed)
    current, correction = _handoff(provider, typed)
    first = provider.rf44_pre_update_state(current=current, correction=correction)
    changed = replace(
        correction,
        delta_v_chart=(9.0, 8.0, 7.0),
        gamma_d_chart=(-3.0, -2.0, -1.0),
        coefficients_u0_u1_u2_s0_s1=(10.0, 11.0, 12.0, 13.0, 14.0),
    )
    second = provider.rf44_pre_update_state(current=current, correction=changed)
    assert first == second


def test_candidate_or_chart_drift_fails_closed(monkeypatch) -> None:
    provider = _provider()
    typed = _typed_receipt(provider)
    monkeypatch.setattr(module, "materialize_current_i4_rf30_typed_defect", lambda *a, **k: typed)
    current, correction = _handoff(provider, typed)

    bad_current = SimpleNamespace(
        candidate=typed.candidate,
        source_chart_id=typed.source_chart_id + "-drift",
        source_chart_sha256=typed.source_chart_sha256,
    )
    with pytest.raises(module.CurrentI4RF44PreStateProviderError, match="source-chart identity drifted"):
        provider.rf44_pre_update_state(current=bad_current, correction=correction)

    bad_identity = CandidateIdentity(
        candidate_id=typed.candidate.candidate_id + "-drift",
        candidate_sha256=typed.candidate.candidate_sha256,
        evidence_kind="repository-candidate",
    )
    bad_correction = replace(correction, candidate=bad_identity)
    with pytest.raises(module.CurrentI4RF44PreStateProviderError, match="candidate identity drifted"):
        provider.rf44_pre_update_state(current=current, correction=bad_correction)


def test_surrogate_or_heldout_rf30_state_fails_before_promotion(monkeypatch) -> None:
    provider = _provider()
    typed = _typed_receipt(provider)
    bad_state = replace(typed.state, surrogate_defect_used=True)
    bad_typed = replace(typed, state=bad_state)
    monkeypatch.setattr(module, "materialize_current_i4_rf30_typed_defect", lambda *a, **k: bad_typed)
    current, correction = _handoff(provider, typed)
    with pytest.raises(
        module.CurrentI4RF44PreStateProviderError,
        match="surrogate/held-out/residual-as-forcing",
    ):
        provider.rf44_pre_update_state(current=current, correction=correction)


def test_prestate_identity_is_separate_and_downstream_methods_remain_absent() -> None:
    provider = _provider()
    assert len(provider.rf44_prestate_semantic_sha256) == 64
    assert len(provider.rf44_prestate_source_blob_sha1) == 40
    assert provider.rf44_prestate_source_blob_sha1 != provider.metadata.source_blob_sha1
    assert not hasattr(provider, "rf44_post_update_state")
    assert not hasattr(provider, "rf44_correction_increment")

    assert tuple(inspect.signature(module.bind_current_i4_rf44_prestate_provider).parameters) == (
        "backend", "radius", "z", "t"
    )
    truth = module.truth_boundary()
    assert truth["exact_a2_1198_raw_provider_delegated_unchanged"] is True
    assert truth["actual_complete_curl_wave_reused"] is True
    assert truth["normalized_source_auxiliary_t2_haar_state_reused"] is True
    assert truth["current_i4_leading_radial_background_recomputed"] is True
    assert truth["physical_to_normalized_velocity_scaling_Q_to_A_applied"] is True
    assert truth["rf30_preupdate_tuple_beta_gamma_v_zero"] is True
    assert truth["rf44_preupdate_state_materialized"] is True
    assert truth["rf44_preupdate_reduces_exactly_to_rf30_operator"] is True
    assert truth["rf44_prestate_has_separate_semantic_identity"] is True
    assert truth["raw_provider_blob_pin_relabelled_as_rf44_prestate_pin"] is False
    assert truth["auxiliary_sign_phase_lift_repository_autonomous"] is True
    assert truth["source_exact_auxiliary_phase_assignment_recovered"] is False
    assert truth["source_exact_slow_derivatives_claimed"] is False
    assert truth["source_exact_rf44_prestate_claimed"] is False
    assert truth["rf44_postupdate_state_materialized"] is False
    assert truth["rf44_correction_increment_materialized"] is False
    assert truth["agent3_mean_correction_solved_here"] is False
    assert truth["cartesian_correction_velocity_materialized"] is False
    assert truth["rf44_rf49_repository_remainder_recomputed"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
    assert truth["residual_defined_free_forcing_allowed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
