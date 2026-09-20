from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_typed_mean_damped_cycle_bridge as bridge
from openai_ns_reconstruction.kokuno_finite_correction_cycle_contract import CorrectionFieldProvider
from openai_ns_reconstruction.kokuno_finite_correction_cycle_ledger import CycleStateProviders
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)


def _before() -> CycleStateProviders:
    identity = CycleIdentity("typed-damped-cycle-bridge-test", 0, "scale-1")
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


def _patch_upstream(monkeypatch: pytest.MonkeyPatch, before: CycleStateProviders, *, wrong_identity: bool = False) -> None:
    monkeypatch.setattr(
        bridge,
        "bind_agent2_complete_curl_backend",
        lambda backend, identity: SimpleNamespace(source_agent2_complete_curl_certified=False),
    )

    def fake_materialize(*args, **kwargs):
        identity = before.identity
        if wrong_identity:
            identity = CycleIdentity(identity.cycle_id, identity.cycle_index, "wrong")
        to_identity = kwargs.get("to_identity")
        if to_identity is None:
            # positional slot immediately before keyword-only viscosity/spatial_step
            to_identity = args[9]
        correction = CorrectionFieldProvider(
            before.identity,
            to_identity,
            "analytic:typed-delta-a->external-complete-curl:oversized-shear",
            lambda x, y, z, t: (-2.4 * t * y, 0.0, 0.0),
            lambda x, y, z, t: (-2.4 * y, 0.0, 0.0),
        )
        return SimpleNamespace(
            identity=identity,
            to_identity=to_identity,
            correction_provider=correction,
            adapter=SimpleNamespace(source_agent2_complete_curl_certified=False),
            to_receipt=lambda: {"regression": True},
        )

    monkeypatch.setattr(bridge, "materialize_typed_mean_velocity_correction_handoff", fake_materialize)


def test_bridge_reuses_heldin_only_damping_and_one_shot_heldout(monkeypatch: pytest.MonkeyPatch) -> None:
    before = _before()
    _patch_upstream(monkeypatch, before)
    held_in, held_out, update = _points()
    report = bridge.materialize_and_admit_typed_mean_damped_step(
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
    search = report.damped_search
    assert search.selected_alpha == 0.5
    assert search.heldout_used_for_alpha_selection is False
    assert search.heldout_trial_evaluation_count == 1
    assert search.selected_step_admitted is True
    assert np.isclose(search.selected_transition.held_in_rms_ratio, 0.2, atol=2e-9)
    assert np.isclose(search.selected_transition.held_out_rms_ratio, 0.2, atol=2e-9)
    assert report.handoff.adapter.source_agent2_complete_curl_certified is False


def test_bridge_rejects_handoff_identity_mutation(monkeypatch: pytest.MonkeyPatch) -> None:
    before = _before()
    _patch_upstream(monkeypatch, before, wrong_identity=True)
    held_in, held_out, update = _points()
    with pytest.raises(ValueError, match="changed the before-state identity"):
        bridge.materialize_and_admit_typed_mean_damped_step(
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


def test_public_bridge_signature_has_no_surrogate_or_retuning_inputs() -> None:
    parameters = set(inspect.signature(bridge.materialize_and_admit_typed_mean_damped_step).parameters)
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
        "damping_grid",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = bridge.truth_boundary()
    assert boundary["typed_mean_radial_chain_recomputed_internally"] is True
    assert boundary["agent2_complete_curl_bound_without_reimplementation"] is True
    assert boundary["heldin_only_frozen_grid_damping_reused"] is True
    assert boundary["heldout_used_for_alpha_selection"] is False
    assert boundary["heldout_evaluated_once_after_selection"] is True
    assert boundary["caller_supplied_damping_grid_allowed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["pde_validated"] is False
