from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_typed_delta_a_complete_curl_backend import (
    TypedDeltaACompleteCurlBackend,
    truth_boundary,
    verification_receipt,
)


def _backend() -> tuple[TypedDeltaACompleteCurlBackend, np.ndarray]:
    radii = np.linspace(0.24, 1.26, 17)
    return (
        TypedDeltaACompleteCurlBackend(
            radii=tuple(float(v) for v in radii),
            reference_time=0.5,
            provenance="focused typed-adapter regression",
        ),
        radii,
    )


def _amplitude(radii: np.ndarray, scale: float = 1.0):
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 6
    delta_a = np.stack((scale * 0.01 * bump, -scale * 0.007 * bump), axis=-1)
    delta_a[[0, -1], :] = 0.0
    return SimpleNamespace(
        delta_a=delta_a,
        base_amplitudes=SimpleNamespace(radii=tuple(float(v) for v in radii)),
    )


def test_adapter_backend_matches_agent3_external_callable_shape() -> None:
    backend, radii = _backend()
    amplitude = _amplitude(radii)
    for evaluator in (backend.velocity_evaluator, backend.velocity_dt_evaluator):
        parameters = tuple(inspect.signature(evaluator).parameters)
        assert parameters == ("amplitude_evaluation", "x", "y", "z", "t")
        value = np.asarray(evaluator(amplitude, 0.51, -0.16, 0.23, 0.5), dtype=float)
        assert value.shape == (3,)
        assert np.all(np.isfinite(value))

    kwargs = backend.agent3_adapter_kwargs(identity="cycle-token")
    assert set(kwargs) == {
        "identity",
        "radii",
        "producer_kind",
        "provenance",
        "velocity_evaluator",
        "velocity_dt_evaluator",
        "source_agent2_complete_curl_certified",
    }
    assert kwargs["identity"] == "cycle-token"
    assert kwargs["radii"] == tuple(float(v) for v in radii)
    assert callable(kwargs["velocity_evaluator"])
    assert callable(kwargs["velocity_dt_evaluator"])
    assert kwargs["source_agent2_complete_curl_certified"] is False


def test_backend_fail_closes_on_grid_shape_and_compactness_drift() -> None:
    backend, radii = _backend()
    good = _amplitude(radii)

    wrong_grid = SimpleNamespace(
        delta_a=good.delta_a,
        base_amplitudes=SimpleNamespace(radii=tuple(float(v) for v in radii + 1.0e-7)),
    )
    with pytest.raises(ValueError, match="radial grid"):
        backend.velocity_evaluator(wrong_grid, 0.5, 0.1, 0.0, 0.5)

    bad_shape = SimpleNamespace(
        delta_a=np.zeros((len(radii), 3)),
        base_amplitudes=good.base_amplitudes,
    )
    with pytest.raises(ValueError, match="shape"):
        backend.velocity_evaluator(bad_shape, 0.5, 0.1, 0.0, 0.5)

    noncompact = _amplitude(radii)
    noncompact.delta_a[0, 0] = 1.0e-3
    with pytest.raises(ValueError, match="endpoints"):
        backend.velocity_evaluator(noncompact, 0.5, 0.1, 0.0, 0.5)


def test_zero_typed_delta_a_maps_to_exact_zero_correction() -> None:
    backend, radii = _backend()
    amplitude = _amplitude(radii, scale=0.0)
    points = ((0.44, 0.11, -0.3, 0.4), (0.73, -0.19, 0.2, 0.5), (0.96, 0.08, 0.5, 0.6))
    for point in points:
        value = backend.velocity_evaluator(amplitude, *point)
        value_dt = backend.velocity_dt_evaluator(amplitude, *point)
        potential = backend.vector_potential_evaluator(amplitude, *point)
        assert np.array_equal(value, np.zeros(3))
        assert np.array_equal(value_dt, np.zeros(3))
        assert np.array_equal(potential, np.zeros(3))


def test_truth_boundary_forbids_premature_scientific_promotion() -> None:
    boundary = truth_boundary()
    assert boundary["typed_delta_a_consumed_only_after_agent3_materialization"] is True
    assert boundary["agent3_mean_radial_chain_reimplemented"] is False
    assert boundary["vector_potential_first_complete_curl_reused"] is True
    assert boundary["forbidden_constructor_parameters_absent"] is True
    assert boundary["source_agent2_complete_curl_certified_for_full_candidate"] is False
    assert boundary["source_mean_amplitude_differential_certified_for_full_candidate"] is False
    assert boundary["independent_agent4_correction_audit_required"] is True
    assert boundary["real_full_candidate_bound"] is False
    assert boundary["heldout_ns_momentum_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["paper_exact"] is False
    assert boundary["pde_validated"] is False


def test_preregistered_adapter_receipt_is_exact_replay_only() -> None:
    receipt = verification_receipt()
    assert receipt["sample_count"] == 12
    assert receipt["failed_guards"] == []
    replay = receipt["exact_replay"]
    assert replay["velocity_max_abs_difference"] == 0.0
    assert replay["velocity_dt_max_abs_difference"] == 0.0
    assert replay["vector_potential_max_abs_difference"] == 0.0
    payload = receipt["agent3_adapter_payload"]
    assert payload["velocity_callable"] is True
    assert payload["velocity_dt_callable"] is True
    assert payload["source_agent2_complete_curl_certified"] is False
