from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_autonomous_oscillatory_mean_debt as debt_module
from openai_ns_reconstruction.kokuno_autonomous_oscillatory_mean_debt import (
    EXPECTED_AUTONOMOUS_MISSING_WEIGHT,
    materialize_autonomous_oscillatory_mean_debt,
)


@pytest.fixture(scope="module")
def report() -> dict[str, object]:
    return materialize_autonomous_oscillatory_mean_debt()


def test_public_api_accepts_no_defect_stress_target_or_gain() -> None:
    params = set(inspect.signature(materialize_autonomous_oscillatory_mean_debt).parameters)
    assert not params.intersection(
        {"defect", "residual", "target", "stress", "requested_stress", "gain", "factor"}
    )


def test_replays_frozen_autonomous_factor_without_theorem_promotion(
    report: dict[str, object],
) -> None:
    factor = report["autonomous_factor"]
    assert factor["origin_pr"] == 514
    assert factor["origin_head"] == "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8"
    assert factor["prepared_N"] == 5
    assert factor["band"] == 5
    assert factor["coordinate_q"] == 1.3
    assert factor["physical_q"] == pytest.approx(0.040625, rel=0.0, abs=1.0e-18)
    assert factor["active_mask_indices"] == (4, 5)
    assert factor["autonomous_missing_weight"] == pytest.approx(
        EXPECTED_AUTONOMOUS_MISSING_WEIGHT, rel=0.0, abs=5.0e-15
    )
    assert factor["squared_partition_sum"] == pytest.approx(1.0, abs=3.0e-14)
    assert factor["profile_is_repository_autonomous"] is True
    assert factor["formal_theorem_machine_bump_identity_claimed"] is False
    assert factor["formal_missing_weight_equality_claimed"] is False
    assert factor["theorem_missing_weight_replaced"] is False
    assert factor["theorem_missing_weight_materialized"] is False


def test_component_debt_is_exact_frozen_scalar_multiple_of_actual_stress(
    report: dict[str, object],
) -> None:
    factor = report["autonomous_factor"]
    stress = report["oscillatory_requested_stress"]
    debt = report["oscillatory_component_mean_debt"]
    weight = float(factor["autonomous_missing_weight"])
    theta = np.asarray(stress["theta_e2_values"], dtype=float)
    axial = np.asarray(stress["axial_e1_values"], dtype=float)
    debt_theta = np.asarray(debt["theta_e2_values"], dtype=float)
    debt_axial = np.asarray(debt["axial_e1_values"], dtype=float)
    np.testing.assert_allclose(debt_theta, -weight * theta, rtol=0.0, atol=1.0e-13)
    np.testing.assert_allclose(debt_axial, -weight * axial, rtol=0.0, atol=1.0e-13)
    assert debt["identity_closure_max_abs"] <= 1.0e-13
    assert debt["debt_rms"] == pytest.approx(
        weight * float(stress["requested_stress_rms"]), rel=2.0e-15, abs=1.0e-15
    )
    assert math.isfinite(float(debt["debt_rms"]))
    assert float(debt["debt_rms"]) > 0.0


def test_truth_boundary_promotes_only_oscillatory_component_debt(
    report: dict[str, object],
) -> None:
    anti = report["anti_surrogate_contract"]
    truth = report["truth_boundary"]
    assert anti["caller_supplied_defect_or_stress_parameters"] == []
    assert anti["surrogate_defect_used"] is False
    assert anti["actual_oscillatory_self_defect_recomputed_internally"] is True
    assert anti["autonomous_factor_frozen_before_this_debt_evaluation"] is True
    assert anti["threshold_or_gain_tuned_after_observation"] is False
    assert truth["real_oscillatory_self_defect_component_consumed"] is True
    assert truth["oscillatory_requested_stress_component_materialized"] is True
    assert truth["autonomous_finite_head_factor_applied"] is True
    assert truth["oscillatory_component_finite_head_mean_debt_materialized"] is True
    assert truth["formal_missing_weight_materialized"] is False
    assert truth["formal_missing_weight_equality_claimed"] is False
    assert truth["full_same_cycle_composite_requested_stress_materialized"] is False
    assert truth["agent1_leading_cross_terms_included"] is False
    assert truth["matched_pressure_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["candidate_finite_head_mean_debt_materialized"] is False
    assert truth["real_full_candidate_defect_consumed"] is False
    assert truth["signed_mean_inverse_input_ready"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["finite_correction_cycle_rerun_allowed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["blowup_proved"] is False


def test_fail_closed_if_upstream_stress_launders_full_composite(monkeypatch: pytest.MonkeyPatch) -> None:
    actual = debt_module.materialize_actual_oscillatory_mean_stress

    def forged(**kwargs: object) -> dict[str, object]:
        result = actual(**kwargs)
        result["truth_boundary"]["full_same_cycle_composite_requested_stress_materialized"] = True
        return result

    monkeypatch.setattr(debt_module, "materialize_actual_oscillatory_mean_stress", forged)
    with pytest.raises(RuntimeError, match="component-only lane"):
        debt_module.materialize_autonomous_oscillatory_mean_debt()
