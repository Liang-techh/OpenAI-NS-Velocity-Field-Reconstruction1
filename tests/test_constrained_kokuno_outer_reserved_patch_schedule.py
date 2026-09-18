import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)


def test_source_scale_relations_and_reserved_interval_order():
    schedule = KokunoOuterReservedPatchSchedule()

    assert schedule.T_d == pytest.approx(math.exp(schedule.M_d) + 10.0)
    assert schedule.log_P_star > schedule.T_d
    assert schedule.T_w == pytest.approx(60.0 * math.log(1.0 / schedule.lambda_outer))
    assert schedule.log_X_R == pytest.approx(
        math.log(110.0) + 10.0 * (math.log(schedule.C) + schedule.log_P_star)
    )
    assert schedule.log_X_w == pytest.approx(schedule.log_X_R + schedule.T_d + 2.0)
    assert schedule.log_P1 == pytest.approx(schedule.log_P_star - 0.2)
    assert schedule.log_e_w == pytest.approx(
        schedule.log_P1 - 0.5 * schedule.T_d - 0.5 - 0.5 * schedule.lambda_outer
    )
    assert schedule.log_c_patch == pytest.approx(
        schedule.log_e_w + (0.5 + schedule.lambda_outer) * schedule.log_X_w
    )

    intervals = schedule.reserved_log_intervals()
    assert intervals["I1"][1] == pytest.approx(intervals["I2"][0])
    assert intervals["I2"][1] < intervals["I3"][0]
    assert intervals["I3"][1] < intervals["I4"][0]

    i2_low, i2_high = intervals["I2"]
    assert i2_low < schedule.log_X_star < i2_high
    # The complete autonomous three-bump support from the existing repair
    # primitive remains strictly inside the source I2 interval.
    assert i2_low < schedule.log_X_star + math.log(0.80)
    assert schedule.log_X_star + math.log(1.70) < i2_high


def test_X_star_e_star_are_source_background_scales_not_independent_guesses():
    schedule = KokunoOuterReservedPatchSchedule()
    scales = schedule.repair_scales()
    X_star = scales["X_star"]
    e_star = scales["e_star"]

    assert math.log(X_star) == pytest.approx(schedule.log_X_star)
    assert math.log(e_star) == pytest.approx(schedule.log_e_star)
    assert schedule.log_e_star == pytest.approx(
        schedule.log_c_patch - (0.5 + schedule.lambda_outer) * schedule.log_X_star
    )

    eta = 0.2
    values = schedule.patch_profile_values(X_star, eta)
    source_f = 1.0 / (1.0 + eta * eta)
    assert values["x_star"] == pytest.approx(1.0)
    assert values["E"] == pytest.approx(e_star * source_f)
    assert values["U"] == pytest.approx(0.0)
    assert values["v0"] == pytest.approx(0.0)
    assert values["Pi_X"] == pytest.approx(values["F"] ** 2)
    assert values["E_X"] == pytest.approx(
        (-0.5 - schedule.lambda_outer) * values["E"] / X_star
    )
    assert values["F_X"] == pytest.approx(
        -(1.0 + schedule.lambda_outer) * values["F"] / X_star
    )
    assert np.isfinite(values["E_eta"])
    assert np.isfinite(values["F_eta"])


def test_reserved_patch_velocity_is_vectorized_finite_and_nontrivial():
    schedule = KokunoOuterReservedPatchSchedule()
    X_star = schedule.repair_scales()["X_star"]

    # At z=0,t=1/2 the native source coordinates have q=1/2 and X=r^2.
    # Choosing r=sqrt(X_*) therefore probes the exact centre of I2.
    r = math.sqrt(X_star)
    velocity = schedule.patch_velocity(
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

    profile = schedule.patch_profile_values(X_star, 0.0)
    expected_u_theta = (0.5 ** (-0.5 - schedule.h)) * profile["E"]
    assert velocity[0, 1] == pytest.approx(expected_u_theta, rel=2.0e-12)


def test_log_scale_report_remains_available_when_physical_scale_overflows():
    schedule = KokunoOuterReservedPatchSchedule(M_d=6.0)
    report = schedule.log_scale_report()

    assert math.isfinite(report["log_X_star"])
    assert math.isfinite(report["log_e_star"])
    assert report["X_star"] is None
    with pytest.raises(OverflowError):
        schedule.repair_scales()


def test_serialization_round_trip_and_truth_boundary_are_fail_closed(tmp_path):
    schedule = KokunoOuterReservedPatchSchedule()
    path = tmp_path / "outer_reserved_patch.json"
    schedule.save(path)
    loaded = KokunoOuterReservedPatchSchedule.load(path)

    assert loaded == schedule
    assert loaded.sha256 == schedule.sha256
    payload = schedule.to_payload()
    assert payload["truth_boundary"]["source_X_star_e_star_relation_executable"] is True
    assert payload["truth_boundary"]["source_hidden_numeric_choices_recovered"] is False
    assert payload["truth_boundary"]["actual_heat_discrepancy_applied_to_patch"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False

    tampered = copy.deepcopy(payload)
    tampered["derived_log_scales"]["log_X_star"] += 1.0
    with pytest.raises(ValueError):
        KokunoOuterReservedPatchSchedule.from_payload(tampered)
