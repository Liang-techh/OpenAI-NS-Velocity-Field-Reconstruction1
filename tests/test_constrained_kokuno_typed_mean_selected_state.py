from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_typed_mean_damped_cycle_bridge as bridge
import openai_ns_reconstruction.kokuno_typed_mean_selected_state as selected
from openai_ns_reconstruction.kokuno_finite_correction_cycle_contract import CorrectionFieldProvider
from openai_ns_reconstruction.kokuno_finite_correction_cycle_ledger import CycleStateProviders
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)


def _before() -> CycleStateProviders:
    identity = CycleIdentity("typed-selected-state-test", 0, "scale-1")
    return CycleStateProviders(
        VectorFieldProvider(identity, "analytic:u=t*y e_x", lambda x, y, z, t: (t * y, 0.0, 0.0)),
        VectorFieldProvider(identity, "analytic:u_t=y e_x", lambda x, y, z, t: (y, 0.0, 0.0)),
        ScalarFieldProvider(identity, "analytic:p=0", lambda x, y, z, t: 0.0),
        RestrictedForcingProvider(
            identity,
            "analytic:fixed-zero-forcing",
            "regression zero forcing fixed independently of residual",
            lambda x, y, z, t: (0.0, 0.0, 0.0),
        ),
    )


def _points():
    held_in = [
        (-0.31, 0.27, 0.14, 0.37),
        (0.22, -0.19, -0.33, 0.63),
        (0.17, 0.41, -0.12, 0.46),
    ]
    held_out = [
        (0.11, 0.36, -0.21, 0.41),
        (-0.28, -0.17, 0.29, 0.59),
        (0.33, -0.44, 0.18, 0.52),
    ]
    update = [
        (0.09, 0.24, -0.16, 0.35),
        (-0.21, -0.32, 0.27, 0.55),
        (0.25, 0.38, 0.11, 0.67),
    ]
    return held_in, held_out, update


def _patch_upstream(monkeypatch: pytest.MonkeyPatch, before: CycleStateProviders) -> None:
    monkeypatch.setattr(
        bridge,
        "bind_agent2_complete_curl_backend",
        lambda backend, identity: SimpleNamespace(source_agent2_complete_curl_certified=False),
    )

    def fake_materialize(*args, **kwargs):
        to_identity = kwargs.get("to_identity")
        if to_identity is None:
            to_identity = args[9]
        correction = CorrectionFieldProvider(
            before.identity,
            to_identity,
            "analytic:typed-delta-a->external-complete-curl:oversized-shear",
            lambda x, y, z, t: (-2.4 * t * y, 0.0, 0.0),
            lambda x, y, z, t: (-2.4 * y, 0.0, 0.0),
        )
        return SimpleNamespace(
            identity=before.identity,
            to_identity=to_identity,
            correction_provider=correction,
            adapter=SimpleNamespace(source_agent2_complete_curl_certified=False),
            to_receipt=lambda: {"regression": True},
        )

    monkeypatch.setattr(bridge, "materialize_typed_mean_velocity_correction_handoff", fake_materialize)


def _run(monkeypatch: pytest.MonkeyPatch):
    before = _before()
    _patch_upstream(monkeypatch, before)
    held_in, held_out, update = _points()
    result = selected.materialize_selected_typed_mean_damped_step(
        before,
        geometry=object(),
        reference_scale=object(),
        reference_operator=object(),
        base_amplitudes=object(),
        agent2_backend=object(),
        held_in_points=held_in,
        held_out_points=held_out,
        update_check_points=update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    return before, result


def test_selected_state_is_exact_admitted_damped_transition(monkeypatch: pytest.MonkeyPatch) -> None:
    before, result = _run(monkeypatch)
    search = result.bridge_report.damped_search
    assert search.selected_alpha == 0.5
    assert search.heldout_used_for_alpha_selection is False
    assert search.heldout_trial_evaluation_count == 1
    assert search.selected_step_admitted is True
    assert result.before_identity == before.identity
    assert result.after_identity == search.selected_transition.to_identity
    assert result.correction.from_identity == before.identity
    assert result.correction.to_identity == result.after_identity
    assert result.correction.source_ref == search.selected_transition.correction_source_ref
    assert result.after_state.pressure.source_ref == before.pressure.source_ref
    assert result.after_state.restricted_forcing.source_ref == before.restricted_forcing.source_ref
    assert (
        result.after_state.restricted_forcing.restriction_receipt
        == before.restricted_forcing.restriction_receipt
    )

    point = (0.19, 0.31, -0.27, 0.47)
    before_u = np.asarray(before.velocity.evaluator(*point), dtype=float)
    delta_u = np.asarray(result.correction.velocity_evaluator(*point), dtype=float)
    after_u = np.asarray(result.after_state.velocity.evaluator(*point), dtype=float)
    before_ut = np.asarray(before.velocity_dt.evaluator(*point), dtype=float)
    delta_ut = np.asarray(result.correction.velocity_dt_evaluator(*point), dtype=float)
    after_ut = np.asarray(result.after_state.velocity_dt.evaluator(*point), dtype=float)
    assert np.array_equal(after_u, before_u + delta_u)
    assert np.array_equal(after_ut, before_ut + delta_ut)
    assert np.isclose(after_u[0] / before_u[0], -0.2, atol=2e-15)
    assert np.isclose(after_ut[0] / before_ut[0], -0.2, atol=2e-15)


def test_selected_state_replay_rejects_provenance_drift(monkeypatch: pytest.MonkeyPatch) -> None:
    before = _before()
    _patch_upstream(monkeypatch, before)
    held_in, held_out, update = _points()
    original = selected.damping._build_trial
    calls = 0

    def mutated_build(*args, **kwargs):
        nonlocal calls
        calls += 1
        state, correction = original(*args, **kwargs)
        if calls == 5:  # four frozen-grid trials, then the selected-state replay
            correction = CorrectionFieldProvider(
                correction.from_identity,
                correction.to_identity,
                correction.source_ref + ":mutated",
                correction.velocity_evaluator,
                correction.velocity_dt_evaluator,
            )
        return state, correction

    monkeypatch.setattr(selected.damping, "_build_trial", mutated_build)
    with pytest.raises(ValueError, match="provenance drifted"):
        selected.materialize_selected_typed_mean_damped_step(
            before,
            geometry=object(),
            reference_scale=object(),
            reference_operator=object(),
            base_amplitudes=object(),
            agent2_backend=object(),
            held_in_points=held_in,
            held_out_points=held_out,
            update_check_points=update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )
    assert calls == 5


def test_public_selected_state_signature_has_no_selection_or_surrogate_inputs() -> None:
    parameters = set(
        inspect.signature(selected.materialize_selected_typed_mean_damped_step).parameters
    )
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "rhs",
        "delta_y",
        "delta_a",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "alpha",
        "damping_grid",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = selected.truth_boundary()
    assert boundary["selected_damped_after_state_materialized"] is True
    assert boundary["selected_correction_provider_materialized"] is True
    assert boundary["selected_alpha_caller_supplied"] is False
    assert boundary["heldout_used_for_alpha_selection"] is False
    assert boundary["heldout_evaluated_again_during_state_materialization"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["pde_validated"] is False
