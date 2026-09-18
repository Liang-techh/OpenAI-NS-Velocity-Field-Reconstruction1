from __future__ import annotations

import copy
from decimal import Decimal
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_hierarchical_heat_repair import (
    KokunoHierarchicalHeatRepair,
)


def test_default_signed_log_target_closes_on_frozen_discrete_three_bump_map() -> None:
    repair = KokunoHierarchicalHeatRepair()
    solution = repair.solve(0.2)

    assert solution.target_sign == (1, -1, 1)
    assert solution.precision_digits > 450
    assert solution.newton_steps >= 1
    assert solution.max_relative_residual < 1.0e-40
    assert solution.term_labels[:3] == (
        "linear_C_p",
        "linear_S",
        "linear_I_sub",
    )
    assert any(label.startswith("newton_") for label in solution.term_labels)

    target = [Decimal(value) for value in solution.target]
    achieved = [Decimal(value) for value in solution.achieved]
    for expected, actual in zip(target, achieved):
        assert expected != 0
        assert abs((actual - expected) / expected) < Decimal("1e-40")


def test_hierarchy_preserves_channels_that_parent_float64_target_cannot_materialize() -> None:
    repair = KokunoHierarchicalHeatRepair()
    report = repair.precision_report(0.2)
    solution = repair.solve(0.2)

    assert report["target"]["underflow_channels"] == ["C_p", "S"]
    assert report["target"]["span_decades"] > 400.0
    assert report["float64_velocity_preserves_all_repair_channels"] is False
    assert report["continuous_moment_compensation_certified"] is False

    # All three source targets remain explicitly nonzero in Decimal, including
    # the two that ordinary float64 would silently erase.
    assert all(Decimal(value) != 0 for value in solution.target)
    assert all(Decimal(value) != 0 for value in solution.coefficients)


def test_local_I2_decimal_velocity_carries_nonzero_repair_while_float_view_is_compatible() -> None:
    repair = KokunoHierarchicalHeatRepair()
    outer = repair.outer_schedule
    X_star = outer.repair_scales()["X_star"]
    r = math.sqrt(X_star)

    base = np.asarray(outer.patch_velocity(r, 0.0, 0.0, 0.5), dtype=float)
    corrected_decimal = repair.velocity_decimal(r, 0.0, 0.0, 0.5)
    corrected_float = repair.velocity(r, 0.0, 0.0, 0.5)

    assert corrected_float.shape == (3,)
    assert np.all(np.isfinite(corrected_float))
    assert corrected_decimal[0] == Decimal("0.0")
    assert corrected_decimal[2] == Decimal("0.0")
    delta_theta = corrected_decimal[1] - Decimal(repr(float(base[1])))
    assert delta_theta != 0
    assert abs(delta_theta) < abs(Decimal(repr(float(base[1]))))
    # Compatibility conversion is allowed to lose the sub-ulp repair; the
    # authoritative local repaired value is the Decimal evaluator above.
    assert corrected_float[1] == pytest.approx(base[1], rel=0.0, abs=0.0)


def test_float_velocity_wrapper_is_vectorized_on_the_reserved_patch() -> None:
    repair = KokunoHierarchicalHeatRepair()
    X_star = repair.outer_schedule.repair_scales()["X_star"]
    r = math.sqrt(X_star)

    velocity = repair.velocity(
        np.asarray([r, 0.9 * r]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.5, 0.5]),
    )
    assert velocity.shape == (2, 3)
    assert np.all(np.isfinite(velocity))
    assert np.all(velocity[:, 0] == 0.0)
    assert np.all(velocity[:, 2] == 0.0)
    assert np.all(velocity[:, 1] > 0.0)


def test_eta_endpoints_have_exact_zero_repair() -> None:
    repair = KokunoHierarchicalHeatRepair()
    for eta in (-1.0, 1.0):
        solution = repair.solve(eta)
        assert solution.target_sign == (0, 0, 0)
        assert solution.coefficients == ("0", "0", "0")
        assert solution.max_relative_residual == 0.0
        assert solution.newton_steps == 0


def test_serialization_round_trip_and_truth_promotion_fail_closed(tmp_path) -> None:
    repair = KokunoHierarchicalHeatRepair()
    path = tmp_path / "hierarchical_heat_repair.json"
    repair.save_json(path)
    loaded = KokunoHierarchicalHeatRepair.load_json(path)

    assert loaded.sha256 == repair.sha256
    assert loaded.to_payload() == repair.to_payload()

    payload = copy.deepcopy(repair.to_payload())
    truth = payload["truth_boundary"]
    assert truth["actual_heat_discrepancy_applied_to_I2_hierarchical_representation"] is True
    assert truth["float64_preserves_all_repair_channels"] is False
    assert truth["continuous_source_moment_compensation_certified"] is False
    assert truth["heat_compensation_completed"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    payload["truth_boundary"]["heat_compensation_completed"] = True
    with pytest.raises(ValueError, match="truth boundary mismatch"):
        KokunoHierarchicalHeatRepair.from_payload(payload)
