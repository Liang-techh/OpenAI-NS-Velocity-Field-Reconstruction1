import copy
from decimal import Decimal
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i2_heat_repair import (
    I2_AUDIT_HEAD,
    I2_AUDIT_PR,
    KokunoPA16CurrentCartesianI2HeatRepair,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianI2HeatRepair()


def _active_i2_X(candidate, x_star=1.24):
    return math.exp(candidate.i2.outer_schedule.log_X_star + math.log(x_star))


def test_exact_frozen_i1_parent_replay_before_and_at_exit(candidate):
    eta = np.asarray([-0.43, 0.0, 0.61])
    X = np.asarray(
        [0.75 * candidate.X_I1_end, candidate.X_I1_end, candidate.X_I1_end]
    )
    child = candidate.similarity_profile_values(X, eta)
    parent = candidate.parent.similarity_profile_values(X, eta)
    mapping = {
        "F_current_i2_repair": "F_current_i1_repair",
        "U_current_i2_repair": "U_current_i1_repair",
        "E_current_i2_repair": "E_current_i1_repair",
        "M_over_X_current_i2_repair": "M_over_X_current_i1_repair",
        "M_eta_over_X_current_i2_repair": "M_eta_over_X_current_i1_repair",
        "v0_current_i2_repair": "v0_current_i1_repair",
    }
    for child_key, parent_key in mapping.items():
        np.testing.assert_array_equal(child[child_key], parent[parent_key])


def test_i2_dependency_is_same_outer_schedule_and_audited_component(candidate):
    current_schedule = candidate.parent.pre_i1.modulation.outer_schedule
    assert candidate.i2.outer_schedule.to_payload() == current_schedule.to_payload()
    assert candidate.log_X_I1_end <= candidate.log_X_I2_start
    assert candidate.log_X_I2_start < candidate.log_X_I2_end <= candidate.log_X_I3_start
    assert candidate.i2.sha256 == candidate.configuration()["bound_high_precision_i2"]["sha256"]
    assert candidate.configuration()["bound_high_precision_i2"]["independent_audit_pr"] == I2_AUDIT_PR == 339
    assert candidate.configuration()["bound_high_precision_i2"]["independent_audit_head"] == I2_AUDIT_HEAD


def test_actual_i1_exit_memory_is_carried_as_constant_physical_M(candidate):
    eta = np.asarray([-0.35, 0.0, 0.42])
    exit_total = candidate.parent.similarity_profile_values(
        np.full(eta.shape, candidate.X_I1_end), eta
    )
    exit_base = candidate.current.similarity_profile_values(
        np.full(eta.shape, candidate.X_I1_end), eta
    )
    dm_exit = np.asarray(exit_total["M_over_X_current_i1_repair"]) - np.asarray(
        exit_base["M_over_X_current_rf40_power_law"]
    )
    dm_eta_exit = np.asarray(
        exit_total["M_eta_over_X_current_i1_repair"]
    ) - np.asarray(exit_base["M_eta_over_X_current_rf40_power_law"])

    X = np.full(eta.shape, candidate.X_I2_start)
    child = candidate.similarity_profile_values(X, eta)
    np.testing.assert_allclose(
        child["delta_M_over_X_from_I1_exit"],
        (candidate.X_I1_end / X) * dm_exit,
        rtol=5e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        child["delta_M_eta_over_X_from_I1_exit"],
        (candidate.X_I1_end / X) * dm_eta_exit,
        rtol=5e-14,
        atol=0.0,
    )


def test_active_i2_uses_pure_swirl_heat_delta_without_resetting_memory(candidate):
    X = _active_i2_X(candidate)
    assert candidate.X_I2_start < X < candidate.X_I2_end
    eta = np.asarray(0.0)
    child = candidate.similarity_profile_values(np.asarray(X), eta)
    base = candidate.current.similarity_profile_values(np.asarray(X), eta)

    # The heat patch changes E/F only.  U, hence physical prefix M, is not
    # changed by the I2 local repair; the already-existing I1 remainder is kept.
    np.testing.assert_array_equal(
        child["U_current_i2_repair"], base["U_current_rf40_power_law"]
    )
    dm, dm_eta = candidate._carried_i1_memory(np.asarray(X), eta)
    np.testing.assert_allclose(
        child["M_over_X_current_i2_repair"],
        np.asarray(base["M_over_X_current_rf40_power_law"]) + dm,
        rtol=5e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        child["M_eta_over_X_current_i2_repair"],
        np.asarray(base["M_eta_over_X_current_rf40_power_law"]) + dm_eta,
        rtol=5e-14,
        atol=0.0,
    )

    receipt = candidate.i2_decimal_correction_report(X, 0.0)
    assert Decimal(receipt["delta_E_decimal"]) != 0
    assert Decimal(receipt["delta_F_decimal"]) != 0
    assert receipt["precision_digits"] >= 128
    assert receipt["solution_sha256"]


def test_i2_analytic_decimal_radial_derivative_matches_centered_difference(candidate):
    X = _active_i2_X(candidate)
    eta = 0.0
    eps = 2.0e-6
    dE, dEX, _, _ = candidate._i2_delta_decimal(X, eta)
    plus = candidate._i2_delta_decimal(X * (1.0 + eps), eta)[0]
    minus = candidate._i2_delta_decimal(X * (1.0 - eps), eta)[0]
    finite_difference = (plus - minus) / Decimal(repr(2.0 * eps * X))
    assert dE != 0
    assert dEX != 0
    relative = abs((finite_difference - dEX) / dEX)
    assert float(relative) < 2.0e-5


def test_physical_M_identity_remains_constant_across_active_i2(candidate):
    eta = np.asarray(0.0)
    X0 = _active_i2_X(candidate, 1.10)
    X1 = _active_i2_X(candidate, 1.40)
    p0 = candidate.similarity_profile_values(np.asarray(X0), eta)
    p1 = candidate.similarity_profile_values(np.asarray(X1), eta)
    M0 = X0 * float(np.asarray(p0["M_over_X_current_i2_repair"]))
    M1 = X1 * float(np.asarray(p1["M_over_X_current_i2_repair"]))
    scale = max(1.0, abs(M0), abs(M1))
    assert abs(M1 - M0) / scale < 3.0e-12
    assert float(np.asarray(p0["U_current_i2_repair"])) == 0.0
    assert float(np.asarray(p1["U_current_i2_repair"])) == 0.0


def test_vectorized_cartesian_velocity_and_axis_regularity(candidate):
    # q=1, eta=0 gives t=0,z=0 and X=r^2/2, allowing a direct active-I2
    # Cartesian probe without solving an inverse geometry by hand.
    X = _active_i2_X(candidate)
    r = math.sqrt(2.0 * X)
    points = candidate.velocity(
        np.asarray([0.0, r]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.5, 0.0]),
    )
    assert points.shape == (2, 3)
    assert np.all(np.isfinite(points))
    assert points[0, 0] == 0.0
    assert points[0, 1] == 0.0


def test_radial_derivative_api_is_finite_and_i2_delta_is_analytic(candidate):
    X = _active_i2_X(candidate)
    out = candidate.similarity_radial_derivatives(
        np.asarray([X, X * 1.01]), np.asarray([0.0, 0.0])
    )
    for key in (
        "F_current_i2_repair_X",
        "U_current_i2_repair_X",
        "E_current_i2_repair_X",
        "delta_F_high_precision_I2_X_float_view",
        "delta_E_high_precision_I2_X_float_view",
    ):
        assert np.all(np.isfinite(out[key]))


def test_save_load_binds_exact_parent_and_immutable_i2_dependency(candidate, tmp_path):
    path = tmp_path / "current_i2.json"
    payload = candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianI2HeatRepair.load_configuration(path)
    assert restored.configuration() == payload == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_high_precision_i2"]["mutable"] = True
    with pytest.raises(ValueError, match="immutable"):
        KokunoPA16CurrentCartesianI2HeatRepair.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_high_precision_i2"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="semantic SHA"):
        KokunoPA16CurrentCartesianI2HeatRepair.from_configuration(mutated)


def test_fails_closed_before_i3_and_preserves_truth_boundary(candidate):
    with pytest.raises(ValueError, match="current I2-repair domain"):
        candidate.similarity_profile_values(
            math.nextafter(candidate.X_I2_end, math.inf), 0.0
        )

    truth = candidate.truth_boundary
    assert truth["current_I2_overlay_applied"] is True
    assert truth["current_I1_exit_M_memory_carried_through_I2"] is True
    assert truth["independent_local_I2_moment_audit_exists"] is True
    for key in (
        "source_hidden_loop_parameters_recovered",
        "paper_exact_I2_bump_parameters_recovered",
        "current_I3_overlay_applied",
        "current_I4_overlay_applied",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        assert truth[key] is False
