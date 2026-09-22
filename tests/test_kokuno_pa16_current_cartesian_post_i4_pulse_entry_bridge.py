import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_post_i4_pulse_entry_bridge import (
    PARENT_EXACT_HEAD,
    SOURCE_BLOB,
    KokunoPA16CurrentCartesianPostI4PulseEntryBridge,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianPostI4PulseEntryBridge()


def _bridge_midpoint(candidate):
    return math.exp(0.5 * (candidate.log_X_I4_end + candidate.log_X_p))


def test_exact_i4_parent_replay_through_i4_exit(candidate):
    eta = np.asarray([-0.37, 0.0, 0.51])
    X = np.asarray(
        [0.97 * candidate.X_I4_end, candidate.X_I4_end, candidate.X_I4_end]
    )
    child = candidate.similarity_profile_values(X, eta)
    parent = candidate.parent.similarity_profile_values(X, eta)
    mapping = {
        "F_current_leading_to_pulse_entry": "F_current_leading_through_I4",
        "U_current_leading_to_pulse_entry": "U_current_leading_through_I4",
        "E_current_leading_to_pulse_entry": "E_current_leading_through_I4",
        "M_over_X_current_leading_to_pulse_entry": "M_over_X_current_leading_through_I4",
        "M_eta_over_X_current_leading_to_pulse_entry": "M_eta_over_X_current_leading_through_I4",
        "v0_current_leading_to_pulse_entry": "v0_current_leading_through_I4",
    }
    for child_key, parent_key in mapping.items():
        np.testing.assert_array_equal(child[child_key], parent[parent_key])


def test_corrected_schedule_has_exact_three_log_unit_bridge_to_pulse_entry(candidate):
    assert PARENT_EXACT_HEAD == "b06742ca6e189499192ede3cce40f62cdc1e35ca"
    assert math.isclose(
        candidate.log_X_p - candidate.log_X_I4_end,
        3.0,
        rel_tol=0.0,
        abs_tol=4e-12,
    )
    assert math.isclose(
        candidate.log_X_p,
        candidate.current.log_X_4,
        rel_tol=0.0,
        abs_tol=2e-12,
    )
    assert candidate.X_p == candidate.current.X_4
    cfg = candidate.configuration()["bound_source_scope"]
    assert cfg["source_blob"] == SOURCE_BLOB
    assert cfg["bridge_role"] == "post_I4_constant_lambda_to_pulse_entry"
    assert cfg["pulse_entry_identity"] == "X_p_equals_current_X_4"
    assert cfg["source_axial_pulse_materialized"] is False
    assert cfg["mutable"] is False


def test_post_i4_bridge_reuses_exact_current_rf40_power_law(candidate):
    X = np.asarray(
        [
            math.exp(candidate.log_X_I4_end + 0.25),
            _bridge_midpoint(candidate),
            math.exp(candidate.log_X_p - 0.2),
            candidate.X_p,
        ]
    )
    eta = np.asarray([-0.41, 0.0, 0.37, 0.58])
    child = candidate.similarity_profile_values(X, eta)
    base = candidate.current.similarity_profile_values(X, eta)
    np.testing.assert_array_equal(
        child["F_current_leading_to_pulse_entry"],
        base["F_current_rf40_power_law"],
    )
    np.testing.assert_array_equal(
        child["U_current_leading_to_pulse_entry"],
        base["U_current_rf40_power_law"],
    )
    np.testing.assert_array_equal(
        child["E_current_leading_to_pulse_entry"],
        base["E_current_rf40_power_law"],
    )
    assert np.all(child["U_current_leading_to_pulse_entry"] == 0.0)
    assert np.all(child["region"] == "post_I4_pre_pulse_constant_lambda_bridge")


def test_actual_i4_exit_memory_is_carried_as_one_over_X(candidate):
    eta = np.asarray([-0.33, 0.0, 0.47])
    seam_total = candidate.parent.similarity_profile_values(
        np.full(eta.shape, candidate.X_I4_end), eta
    )
    seam_base = candidate.current.similarity_profile_values(
        np.full(eta.shape, candidate.X_I4_end), eta
    )
    dm = np.asarray(seam_total["M_over_X_current_leading_through_I4"]) - np.asarray(
        seam_base["M_over_X_current_rf40_power_law"]
    )
    dm_eta = np.asarray(
        seam_total["M_eta_over_X_current_leading_through_I4"]
    ) - np.asarray(seam_base["M_eta_over_X_current_rf40_power_law"])

    X = np.full(eta.shape, _bridge_midpoint(candidate))
    child = candidate.similarity_profile_values(X, eta)
    factor = candidate.X_I4_end / X
    np.testing.assert_allclose(
        child["delta_M_over_X_carried_from_I4_exit"],
        factor * dm,
        rtol=5e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        child["delta_M_eta_over_X_carried_from_I4_exit"],
        factor * dm_eta,
        rtol=5e-14,
        atol=0.0,
    )


def test_physical_M_stays_constant_from_i4_exit_to_pulse_entry(candidate):
    eta = np.asarray(0.23)
    X0 = math.exp(candidate.log_X_I4_end + 0.4)
    X1 = _bridge_midpoint(candidate)
    X2 = candidate.X_p
    physical_M = []
    for X in (X0, X1, X2):
        p = candidate.similarity_profile_values(np.asarray(X), eta)
        assert float(np.asarray(p["U_current_leading_to_pulse_entry"])) == 0.0
        physical_M.append(
            X * float(np.asarray(p["M_over_X_current_leading_to_pulse_entry"]))
        )
    scale = max(1.0, *(abs(v) for v in physical_M))
    assert (max(physical_M) - min(physical_M)) / scale < 6e-12


def test_bridge_radial_derivatives_reuse_power_law_and_match_centered_X(candidate):
    X = _bridge_midpoint(candidate)
    eta = 0.19
    out = candidate.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    base = candidate.current.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    np.testing.assert_array_equal(
        out["F_current_leading_to_pulse_entry_X"],
        base["F_current_rf40_power_law_X"],
    )
    np.testing.assert_array_equal(
        out["U_current_leading_to_pulse_entry_X"],
        base["U_current_rf40_power_law_X"],
    )
    np.testing.assert_array_equal(
        out["E_current_leading_to_pulse_entry_X"],
        base["E_current_rf40_power_law_X"],
    )

    eps = 2e-6
    plus = candidate.similarity_profile_values(X * (1 + eps), eta)
    minus = candidate.similarity_profile_values(X * (1 - eps), eta)
    fd = (
        float(np.asarray(plus["F_current_leading_to_pulse_entry"]))
        - float(np.asarray(minus["F_current_leading_to_pulse_entry"]))
    ) / (2 * eps * X)
    analytic = float(np.asarray(out["F_current_leading_to_pulse_entry_X"]))
    assert abs(fd - analytic) / max(abs(analytic), 1e-300) < 2e-5


def test_vectorized_cartesian_velocity_reaches_pre_pulse_bridge_and_axis_is_regular(candidate):
    X = _bridge_midpoint(candidate)
    r = math.sqrt(2.0 * X)
    values = candidate.velocity(
        np.asarray([0.0, r]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.5, 0.0]),
    )
    assert values.shape == (2, 3)
    assert np.all(np.isfinite(values))
    assert values[0, 0] == 0.0
    assert values[0, 1] == 0.0
    assert float(np.linalg.norm(values[1])) > 0.0


def test_save_load_binds_parent_and_pre_pulse_source_scope(candidate, tmp_path):
    path = tmp_path / "current_post_i4_pulse_entry_bridge.json"
    payload = candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianPostI4PulseEntryBridge.load_configuration(path)
    assert restored.configuration() == payload == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["source_axial_pulse_materialized"] = True
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianPostI4PulseEntryBridge.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["pulse_entry_identity"] = "caller_selected"
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianPostI4PulseEntryBridge.from_configuration(mutated)


def test_fails_closed_after_pulse_entry_and_does_not_promote_pulse_or_pde(candidate):
    with pytest.raises(ValueError, match="post-I4-to-pulse-entry domain"):
        candidate.similarity_profile_values(math.nextafter(candidate.X_p, math.inf), 0.0)

    truth = candidate.truth_boundary
    assert truth["source_post_I4_pre_pulse_constant_lambda_bridge_identified"] is True
    assert truth["source_post_I4_pre_pulse_bridge_log_length_is_three"] is True
    assert truth["current_leading_preserved_from_I4_exit_to_pulse_entry"] is True
    assert truth["current_cartesian_velocity_executable_to_pulse_entry"] is True
    for key in (
        "source_positive_order_I3_profiles_materialized",
        "I3_positive_order_correction_materialized",
        "source_I4_mean_correction_materialized",
        "current_I4_mean_correction_materialized",
        "source_axial_pulse_materialized",
        "source_pulse_amplitude_root_materialized",
        "source_pulse_M_J_end_correction_materialized",
        "source_terminal_tail_schedule_bound_into_current_velocity",
        "source_exterior_heat_replacement_materialized",
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
