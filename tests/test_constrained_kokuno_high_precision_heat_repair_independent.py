from __future__ import annotations

import math

import mpmath as mp

from openai_ns_reconstruction.kokuno_high_precision_heat_repair_independent import (
    FORMAL_DIVERGENCE_GATE,
    FORMAL_MOMENTUM_GATE,
    _independent_base_moment_strings,
    _max_relative_moment_delta,
)


def test_tanh_transform_is_exact_for_autonomous_bump() -> None:
    with mp.workdps(80):
        for y_text in ("-2.1", "-0.7", "0", "0.9", "2.3"):
            y = mp.mpf(y_text)
            s = mp.tanh(y)
            direct = mp.exp(1 - 1 / (1 - s * s))
            transformed = mp.exp(-(mp.sinh(y) ** 2))
            assert abs(direct - transformed) < mp.mpf("1e-70")


def test_independent_real_line_trapezoid_refines_without_construction_path() -> None:
    # Low-cost regression of the alternate operator itself.  The production
    # construction uses a different tanh-sinh transform and is not imported.
    lambda_text = "0.2"
    coarse = _independent_base_moment_strings(lambda_text, 96, 64)
    fine = _independent_base_moment_strings(lambda_text, 96, 128)
    delta = _max_relative_moment_delta(coarse, fine, precision_digits=96)
    assert delta < mp.mpf("1e-30")


def test_formal_project_gates_remain_fixed() -> None:
    assert math.isclose(FORMAL_MOMENTUM_GATE, 1.0e-3, rel_tol=0.0, abs_tol=0.0)
    assert math.isclose(FORMAL_DIVERGENCE_GATE, 1.0e-5, rel_tol=0.0, abs_tol=0.0)
