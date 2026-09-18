from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy_v2 import (
    KokunoHeatReplacementDiscrepancyV2,
)
from openai_ns_reconstruction.kokuno_log_rescaled_heat_repair_target import (
    KokunoLogRescaledHeatRepairTarget,
)
from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)
from openai_ns_reconstruction.kokuno_terminal_tail_schedule import (
    KokunoTerminalTailSchedule,
)


def _moderate_schedule() -> KokunoTerminalTailSchedule:
    outer = KokunoOuterReservedPatchSchedule(
        M_d=0.1,
        lambda_outer=0.65,
        h=0.005,
    )
    return KokunoTerminalTailSchedule(outer_schedule=outer)


def test_materializable_schedule_matches_corrected_v2_physical_evaluator() -> None:
    schedule = _moderate_schedule()
    target = KokunoLogRescaledHeatRepairTarget(terminal_schedule=schedule)
    heat = KokunoHeatReplacementDiscrepancyV2(
        **schedule.materialize_heat_inputs()
    )
    scales = schedule.outer_schedule.repair_scales()

    eta = np.asarray([0.0, 0.2, 0.5, 0.9])
    expected = heat.normalized_repair_target(eta, **scales)
    actual = target.materialize_target(eta)

    np.testing.assert_allclose(actual, expected, rtol=4.0e-12, atol=0.0)
    assert np.all(actual[:, 0] > 0.0)
    assert np.all(actual[:, 1] < 0.0)
    assert np.all(actual[:, 2] > 0.0)


def test_default_existence_schedule_retains_nonzero_underflowed_channels() -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    encoded = target.signed_log_target(0.2)
    lower = math.log(np.nextafter(0.0, 1.0))

    assert encoded["sign"].tolist() == [1, -1, 1]
    assert encoded["log_abs"][0] < lower
    assert encoded["log_abs"][1] < lower
    assert encoded["log_abs"][2] > lower

    with pytest.raises(OverflowError, match="underflows float64"):
        target.materialize_target(0.2)


def test_default_precision_report_exposes_float64_cancellation_barrier() -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    report = target.precision_report(0.2)

    assert report["underflow_channels"] == ["C_p", "S"]
    assert report["span_decades"] > 400.0
    assert report["estimated_required_decimal_digits"] > 400
    assert report["float64_decimal_digits"] < 20
    assert report["float64_direct_repair_safe"] is False

    log10_abs = report["log10_abs"]
    assert log10_abs[0] == pytest.approx(-664.7, abs=0.2)
    assert log10_abs[1] == pytest.approx(-464.7, abs=0.2)
    assert log10_abs[2] == pytest.approx(-243.3, abs=0.2)


def test_rescaled_target_keeps_channel_logs_when_common_mantissa_underflows() -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    rescaled = target.rescaled_target(0.2)

    assert rescaled["log_reference"] == pytest.approx(rescaled["log_abs"][2])
    assert rescaled["mantissa"][2] == pytest.approx(1.0)
    assert rescaled["mantissa"][0] == 0.0
    assert rescaled["sign"].tolist() == [1, -1, 1]
    assert np.all(np.isfinite(rescaled["log_abs"]))


def test_eta_endpoints_are_exact_zero_without_runtime_warning() -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    encoded = target.signed_log_target(np.asarray([-1.0, 1.0]))
    rescaled = target.rescaled_target(np.asarray([-1.0, 1.0]))

    assert np.all(encoded["sign"] == 0)
    assert np.all(np.isneginf(encoded["log_abs"]))
    assert np.all(rescaled["mantissa"] == 0.0)
    assert np.all(np.isneginf(rescaled["log_reference"]))


def test_serialization_round_trip_and_truth_mutation_fail_closed(tmp_path) -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    path = tmp_path / "log_target.json"
    target.save_json(path)
    loaded = KokunoLogRescaledHeatRepairTarget.load_json(path)

    assert loaded.sha256 == target.sha256
    assert loaded.to_payload() == target.to_payload()

    payload = copy.deepcopy(target.to_payload())
    assert payload["truth_boundary"]["heat_compensation_completed"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    payload["truth_boundary"]["heat_compensation_completed"] = True
    with pytest.raises(ValueError, match="truth boundary mismatch"):
        KokunoLogRescaledHeatRepairTarget.from_payload(payload)


def test_vectorized_signed_log_shape_and_source_vs_autonomous_truth() -> None:
    target = KokunoLogRescaledHeatRepairTarget()
    encoded = target.signed_log_target(np.asarray([[0.0, 0.25], [0.5, 0.75]]))
    truth = target.to_payload()["truth_boundary"]

    assert encoded["sign"].shape == (2, 2, 3)
    assert encoded["log_abs"].shape == (2, 2, 3)
    assert truth["corrected_v2_heat_target_formula_executable"] is True
    assert truth["signed_log_underflow_safe_target"] is True
    assert truth["float64_three_bump_inverse_certified"] is False
    assert truth["actual_heat_discrepancy_applied_to_patch"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["paper_exact"] is False
