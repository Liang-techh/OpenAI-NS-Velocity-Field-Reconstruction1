from __future__ import annotations

import pytest

from openai_ns_reconstruction import constrained_st052m_linear_ramp_fixed_render as mod


def test_frozen_time_activation() -> None:
    assert mod.expected_activation(0.25) == 0.0
    assert mod.expected_activation(0.375) == 0.25
    assert mod.expected_activation(0.50) == 0.50
    assert mod.expected_activation(0.625) == 0.75
    assert mod.expected_activation(0.75) == 1.0
    with pytest.raises(ValueError):
        mod.expected_activation(0.249)
    with pytest.raises(ValueError):
        mod.expected_activation(0.751)


def test_effect_retention_is_descriptive_only() -> None:
    assert mod.effect_retention(0.01, 0.02) == 0.5
    assert mod.effect_retention(-0.01, -0.02) == 0.5
    assert mod.effect_retention(0.0, 0.0) is None


def test_truth_boundary_is_fail_closed() -> None:
    assert mod.TASK_ID == "CR-A9-066"
    assert mod.TEMPORAL_HEAD == "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    assert mod.STATIC_RENDER_HEAD == "ee80de11775bd45157d114000378b6ee3422ff46"
    assert mod.REFERENCE_TIMES == (0.25, 0.375, 0.50, 0.625, 0.75)
    assert mod.ENDPOINT_IDENTITY_ABS_MAX == 1.0e-11
    assert not any(mod.TRUTH.values())
