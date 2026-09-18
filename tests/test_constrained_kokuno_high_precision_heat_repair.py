from __future__ import annotations

from decimal import Decimal, localcontext

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_high_precision_heat_repair import (
    KokunoHighPrecisionHierarchicalHeatRepair,
)
from openai_ns_reconstruction.kokuno_hierarchical_heat_repair import (
    KokunoHierarchicalHeatRepair,
)
from openai_ns_reconstruction.kokuno_hierarchical_heat_repair_independent import (
    _map_decimal,
    independent_moment_tensors,
)


def _max_relative_tensor_error(actual, expected) -> float:
    errors = []
    for a, e in zip(np.ravel(actual), np.ravel(expected)):
        if e != 0.0:
            errors.append(abs(float(a) - float(e)) / abs(float(e)))
        else:
            errors.append(abs(float(a)))
    return max(errors)


def test_high_precision_tensor_matches_refined_independent_moments() -> None:
    repair = KokunoHighPrecisionHierarchicalHeatRepair()
    eta = 0.2
    f_eta = float(repair.outer_schedule.source_f(eta))
    with localcontext() as ctx:
        ctx.prec = 96
        linear, quadratic = repair._discrete_tensors(f_eta)

    linear_float = np.asarray([[float(value) for value in row] for row in linear])
    quadratic_float = np.asarray(
        [
            [[float(value) for value in row] for row in plane]
            for plane in quadratic
        ]
    )
    independent = independent_moment_tensors(
        lambda_outer=repair.outer_schedule.lambda_outer,
        f_eta=f_eta,
        resolution=4096,
    )

    assert _max_relative_tensor_error(linear_float, independent.linear) < 2.0e-12
    assert _max_relative_tensor_error(quadratic_float, independent.quadratic) < 2.0e-12


def test_high_precision_tensor_removes_agent4_dominant_channel_rejection() -> None:
    eta = 0.2
    repair = KokunoHighPrecisionHierarchicalHeatRepair()
    solution = repair.solve(eta)
    coefficients = solution.coefficient_decimals()
    targets = tuple(Decimal(value) for value in solution.target)
    independent = independent_moment_tensors(
        lambda_outer=repair.outer_schedule.lambda_outer,
        f_eta=float(repair.outer_schedule.source_f(eta)),
        resolution=4096,
    )
    achieved = _map_decimal(coefficients, independent, solution.precision_digits)
    independent_i_relative = abs(achieved[2] - targets[2]) / abs(targets[2])

    legacy = KokunoHierarchicalHeatRepair()
    legacy_solution = legacy.solve(eta)
    legacy_achieved = _map_decimal(
        legacy_solution.coefficient_decimals(),
        independent,
        legacy_solution.precision_digits,
    )
    legacy_target = Decimal(legacy_solution.target[2])
    legacy_i_relative = abs(legacy_achieved[2] - legacy_target) / abs(legacy_target)

    assert independent_i_relative < Decimal("1e-10")
    assert independent_i_relative < legacy_i_relative * Decimal("1e-3")
    # The production solve still closes all three rows against its own
    # high-precision tensor.  Tiny-channel independent certification needs an
    # alternate high-precision validator and is deliberately not claimed here.
    assert solution.max_relative_residual < 1.0e-40


def test_high_precision_payload_roundtrip_and_truth_boundary() -> None:
    repair = KokunoHighPrecisionHierarchicalHeatRepair()
    payload = repair.to_payload()
    truth = payload["truth_boundary"]
    assert truth["repository_float64_moment_tensor_used"] is False
    assert truth["high_precision_continuous_moment_tensor_used"] is True
    assert truth["independent_high_precision_moment_audit_completed"] is False
    assert truth["continuous_source_moment_compensation_certified"] is False
    assert truth["heat_compensation_completed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    restored = KokunoHighPrecisionHierarchicalHeatRepair.from_payload(payload)
    assert restored.to_payload() == payload

    tampered = dict(payload)
    tampered["truth_boundary"] = dict(payload["truth_boundary"])
    tampered["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="hash or content mismatch"):
        KokunoHighPrecisionHierarchicalHeatRepair.from_payload(tampered)
