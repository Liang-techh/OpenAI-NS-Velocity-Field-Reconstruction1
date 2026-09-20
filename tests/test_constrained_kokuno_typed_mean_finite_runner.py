from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_typed_mean_finite_runner as runner
from openai_ns_reconstruction.kokuno_finite_correction_cycle_contract import (
    CorrectionFieldProvider,
)
from openai_ns_reconstruction.kokuno_finite_correction_cycle_ledger import (
    CycleStateProviders,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from openai_ns_reconstruction.kokuno_typed_mean_selected_state import (
    TypedSelectedDampedStep,
)


def _state(scale: float, index: int) -> CycleStateProviders:
    identity = CycleIdentity(
        "typed-mean-finite-runner-test",
        index,
        f"scale-{scale:g}",
    )
    return CycleStateProviders(
        VectorFieldProvider(
            identity,
            f"analytic:u={scale:g}*t*y e_x",
            lambda x, y, z, t, s=scale: (s * t * y, 0.0, 0.0),
        ),
        VectorFieldProvider(
            identity,
            f"analytic:u_t={scale:g}*y e_x",
            lambda x, y, z, t, s=scale: (s * y, 0.0, 0.0),
        ),
        ScalarFieldProvider(identity, "analytic:p=0", lambda x, y, z, t: 0.0),
        RestrictedForcingProvider(
            identity,
            "analytic:fixed-zero-forcing",
            "regression zero forcing fixed independently of residual",
            lambda x, y, z, t: (0.0, 0.0, 0.0),
        ),
    )


def _selected_step(
    before: CycleStateProviders,
    after_scale: float,
    *,
    heldout_used: bool = False,
    heldout_count: int = 1,
) -> TypedSelectedDampedStep:
    before_scale = float(before.velocity.source_ref.split("=")[1].split("*")[0])
    after = _state(after_scale, before.identity.cycle_index + 1)
    delta = after_scale - before_scale
    correction = CorrectionFieldProvider(
        before.identity,
        after.identity,
        f"analytic:delta-scale={delta:g}:alpha=0.5",
        lambda x, y, z, t, d=delta: (d * t * y, 0.0, 0.0),
        lambda x, y, z, t, d=delta: (d * y, 0.0, 0.0),
    )
    transition = SimpleNamespace(
        from_identity=before.identity,
        to_identity=after.identity,
        correction_nontrivial=True,
        finite_step_gain_guard_passed=True,
        held_in_before_rms=abs(before_scale),
        held_in_after_rms=abs(after_scale),
        held_out_before_rms=abs(before_scale),
        held_out_after_rms=abs(after_scale),
        held_in_before_max=2.0 * abs(before_scale),
        held_in_after_max=2.0 * abs(after_scale),
        held_out_before_max=3.0 * abs(before_scale),
        held_out_after_max=3.0 * abs(after_scale),
        held_in_rms_ratio=abs(after_scale / before_scale),
        held_out_rms_ratio=abs(after_scale / before_scale),
        held_in_max_ratio=abs(after_scale / before_scale),
        held_out_max_ratio=abs(after_scale / before_scale),
        correction_vector_rms=abs(delta),
        correction_divergence_max=0.0,
        to_receipt=lambda: {
            "from_cycle_index": before.identity.cycle_index,
            "to_cycle_index": after.identity.cycle_index,
        },
    )
    search = SimpleNamespace(
        selected_alpha=0.5,
        heldout_used_for_alpha_selection=heldout_used,
        heldout_trial_evaluation_count=heldout_count,
        selected_transition=transition,
        selected_step_admitted=True,
    )
    bridge = SimpleNamespace(
        damped_search=search,
        to_receipt=lambda: {
            "selected_alpha": 0.5,
            "heldout_trial_evaluation_count": heldout_count,
        },
    )
    return TypedSelectedDampedStep(
        before_identity=before.identity,
        after_state=after,
        correction=correction,
        bridge_report=bridge,
    )


def _points():
    return (
        [(-0.31, 0.27, 0.14, 0.37)],
        [(0.11, 0.36, -0.21, 0.41)],
        [(0.09, 0.24, -0.16, 0.35)],
    )


def _dummy_inputs() -> runner.TypedMeanStageInputs:
    return runner.TypedMeanStageInputs(
        geometry=object(),
        reference_scale=object(),
        reference_operator=object(),
        base_amplitudes=object(),
        agent2_backend=object(),
    )


def test_runner_feeds_exact_selected_state_into_next_typed_stage(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    initial = _state(1.0, 0)
    factory_indices: list[int] = []
    materializer_indices: list[int] = []

    def factory(state: CycleStateProviders) -> runner.TypedMeanStageInputs:
        factory_indices.append(state.identity.cycle_index)
        return _dummy_inputs()

    def fake_materialize(
        before: CycleStateProviders,
        **kwargs: object,
    ) -> TypedSelectedDampedStep:
        materializer_indices.append(before.identity.cycle_index)
        after_scale = 0.5 if before.identity.cycle_index == 0 else 0.25
        return _selected_step(before, after_scale)

    monkeypatch.setattr(
        runner,
        "materialize_selected_typed_mean_damped_step",
        fake_materialize,
    )
    held_in, held_out, update = _points()
    report = runner.materialize_typed_mean_finite_cycle(
        initial,
        factory,
        2,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=0.005,
    )

    assert factory_indices == [0, 1]
    assert materializer_indices == [0, 1]
    assert report.step_count == 2
    assert report.final_state.identity.cycle_index == 2
    assert np.allclose(
        report.final_state.velocity.evaluator(0.2, 0.4, 0.1, 0.5),
        (0.05, 0.0, 0.0),
    )
    assert report.selected_alphas == (0.5, 0.5)
    assert report.held_in_rms_contraction_factors == (0.5, 0.5)
    assert report.held_out_rms_contraction_factors == (0.5, 0.5)
    assert report.cumulative_held_in_rms_ratio == 0.25
    assert report.cumulative_held_out_rms_ratio == 0.25
    assert report.cumulative_held_in_max_ratio == 0.25
    assert report.cumulative_held_out_max_ratio == 0.25
    assert report.cumulative_correction_rms_budget == 0.75
    assert report.maximum_correction_rms == 0.5
    assert report.maximum_correction_divergence_max == 0.0
    assert report.heldout_trial_evaluation_count_total == 2
    assert report.heldout_used_for_alpha_selection is False
    assert report.every_step_nontrivial is True
    assert report.every_step_actual_defect_gain_guard_passed is True
    assert report.source_gain_bounds_passed is True
    assert [stage.stage_index for stage in report.source_stages] == [0, 1]
    assert min(stage.wave_residual_gain for stage in report.source_stages) == pytest.approx(
        0.39999
    )
    assert report.admitted is True


def test_runner_rejects_heldout_selection_contamination(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    initial = _state(1.0, 0)

    monkeypatch.setattr(
        runner,
        "materialize_selected_typed_mean_damped_step",
        lambda before, **kwargs: _selected_step(
            before,
            0.5,
            heldout_used=True,
        ),
    )
    held_in, held_out, update = _points()
    with pytest.raises(ValueError, match="held-out data contaminated"):
        runner.materialize_typed_mean_finite_cycle(
            initial,
            lambda state: _dummy_inputs(),
            1,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_runner_rejects_duplicate_heldout_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    initial = _state(1.0, 0)
    monkeypatch.setattr(
        runner,
        "materialize_selected_typed_mean_damped_step",
        lambda before, **kwargs: _selected_step(
            before,
            0.5,
            heldout_count=2,
        ),
    )
    held_in, held_out, update = _points()
    with pytest.raises(ValueError, match="exactly once"):
        runner.materialize_typed_mean_finite_cycle(
            initial,
            lambda state: _dummy_inputs(),
            1,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_runner_requires_stage_zero_and_positive_finite_length() -> None:
    held_in, held_out, update = _points()
    with pytest.raises(ValueError, match="at least one"):
        runner.materialize_typed_mean_finite_cycle(
            _state(1.0, 0),
            lambda state: _dummy_inputs(),
            0,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=0.005,
        )
    with pytest.raises(ValueError, match="cycle index 0"):
        runner.materialize_typed_mean_finite_cycle(
            _state(0.5, 1),
            lambda state: _dummy_inputs(),
            1,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_public_runner_signature_has_no_surrogate_or_retuning_inputs() -> None:
    parameters = set(
        inspect.signature(runner.materialize_typed_mean_finite_cycle).parameters
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
    boundary = runner.truth_boundary()
    assert boundary["typed_mean_finite_runner_executable"] is True
    assert boundary["selected_after_state_fed_directly_to_next_stage"] is True
    assert boundary["actual_defect_recomputed_inside_each_typed_step"] is True
    assert boundary["runner_reaudits_heldout_transitions"] is False
    assert boundary["heldout_evaluated_once_per_selected_step"] is True
    assert boundary["heldout_used_for_alpha_selection"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["pde_validated"] is False
