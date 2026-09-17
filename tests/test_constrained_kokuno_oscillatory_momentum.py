import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_oscillatory_momentum import (
    DEFAULT_STEPS,
    OscillatoryMomentumIncrementContract,
    screen_frozen_amplitudes,
)


def _affine_base(points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        (
            0.20 + 0.10 * points[:, 1] + 0.02 * time,
            -0.10 + 0.05 * points[:, 0],
            0.07 - 0.03 * points[:, 0],
        )
    )


def _points() -> np.ndarray:
    return np.random.default_rng(9172792).uniform(-0.35, 0.35, size=(32, 3))


def test_analytic_increment_matches_independent_full_composite_with_refinement():
    contract = OscillatoryMomentumIncrementContract()
    correction = KokunoCompleteCurlCorrection()
    errors = []
    for step in DEFAULT_STEPS:
        sample = contract.evaluate(_affine_base, correction, _points(), 0.5, step=step)
        errors.append(sample.metrics()["direct_vs_hybrid_rms"])
        assert sample.metrics()["analytic_correction_divergence_max"] < 1.0e-12

    assert errors[1] < 0.4 * errors[0]
    assert errors[2] < 0.4 * errors[1]
    assert errors[2] < 5.0e-5


def test_zero_amplitude_replays_base_without_fake_improvement():
    contract = OscillatoryMomentumIncrementContract()
    correction = KokunoCompleteCurlCorrection(amplitude=0.0)
    sample = contract.evaluate(_affine_base, correction, _points(), 0.5, step=0.005)
    assert np.max(np.abs(sample.analytic_increment)) == 0.0
    assert np.max(np.abs(sample.direct_composite_momentum - sample.base_momentum)) < 1.0e-13
    assert sample.metrics()["composite_over_base_rms"] == pytest.approx(1.0, abs=1.0e-13)


def test_frozen_screen_keeps_predeclared_grid_and_does_not_select():
    amplitudes = (0.0, 0.125, 0.25)
    rows = screen_frozen_amplitudes(
        _affine_base,
        _points(),
        amplitudes=amplitudes,
        times=(0.5,),
        steps=(0.01,),
    )
    assert [row["amplitude"] for row in rows] == list(amplitudes)
    assert all(row["time"] == 0.5 and row["step"] == 0.01 for row in rows)
    assert all(np.isfinite(list(row.values())).all() for row in rows)


def test_screen_rejects_amplitude_outside_existing_bound():
    with pytest.raises(ValueError, match="existing"):
        screen_frozen_amplitudes(
            _affine_base,
            _points(),
            amplitudes=(2.1,),
            times=(0.5,),
            steps=(0.01,),
        )


def test_contract_rejects_bad_inputs():
    contract = OscillatoryMomentumIncrementContract()
    correction = KokunoCompleteCurlCorrection()
    with pytest.raises(ValueError, match="points"):
        contract.evaluate(_affine_base, correction, np.zeros((0, 3)), 0.5, step=0.01)
    with pytest.raises(ValueError, match="time"):
        contract.evaluate(_affine_base, correction, _points(), 0.8, step=0.01)
    with pytest.raises(ValueError, match="step"):
        contract.evaluate(_affine_base, correction, _points(), 0.5, step=0.0)
