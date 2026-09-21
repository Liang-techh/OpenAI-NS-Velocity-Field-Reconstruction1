import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i4_leading_preservation import (
    PARENT_EXACT_HEAD,
    SOURCE_BLOB,
    KokunoPA16CurrentCartesianI4LeadingPreservation,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianI4LeadingPreservation()


def _i4_midpoint(candidate):
    return math.exp(0.5 * (candidate.log_X_I4_start + candidate.log_X_I4_end))


def test_exact_i3_parent_replay_through_i3_exit(candidate):
    eta = np.asarray([-0.37, 0.0, 0.51])
    X = np.asarray(
        [0.97 * candidate.X_I3_end, candidate.X_I3_end, candidate.X_I3_end]
    )
    child = candidate.similarity_profile_values(X, eta)
    parent = candidate.parent.similarity_profile_values(X, eta)
    mapping = {
        "F_current_leading_through_I4": "F_current_leading_through_I3",
        "U_current_leading_through_I4": "U_current_leading_through_I3",
        "E_current_leading_through_I4": "E_current_leading_through_I3",
        "M_over_X_current_leading_through_I4": "M_over_X_current_leading_through_I3",
        "M_eta_over_X_current_leading_through_I4": "M_eta_over_X_current_leading_through_I3",
        "v0_current_leading_through_I4": "v0_current_leading_through_I3",
    }
    for child_key, parent_key in mapping.items():
        np.testing.assert_array_equal(child[child_key], parent[parent_key])


def test_reserved_schedule_orders_i3_i4_and_source_scope_is_frozen(candidate):
    assert candidate.log_X_I3_end < candidate.log_X_I4_start
    assert candidate.log_X_I4_start < candidate.log_X_I4_end
    assert candidate.X_I4_end <= candidate.current.X_4
    cfg = candidate.configuration()["bound_source_scope"]
    assert cfg["source_blob"] == SOURCE_BLOB
    assert cfg["I4_role"] == "mean_correction_not_leading_E0_U0"
    assert cfg["leading_action"] == "preserve_RF40_power_law"
    assert cfg["mutable"] is False
    assert PARENT_EXACT_HEAD == "5fb7b583062b4a991db86e39fdcdb231022b70c9"


def test_leading_F_U_E_on_I4_are_exact_current_rf40_power_law(candidate):
    X = np.asarray(
        [
            math.exp(candidate.log_X_I4_start + 0.1),
            _i4_midpoint(candidate),
            math.exp(candidate.log_X_I4_end - 0.1),
        ]
    )
    eta = np.asarray([-0.41, 0.0, 0.58])
    child = candidate.similarity_profile_values(X, eta)
    base = candidate.current.similarity_profile_values(X, eta)
    np.testing.assert_array_equal(
        child["F_current_leading_through_I4"], base["F_current_rf40_power_law"]
    )
    np.testing.assert_array_equal(
        child["U_current_leading_through_I4"], base["U_current_rf40_power_law"]
    )
    np.testing.assert_array_equal(
        child["E_current_leading_through_I4"], base["E_current_rf40_power_law"]
    )
    assert np.all(child["U_current_leading_through_I4"] == 0.0)
    assert np.all(child["region"] == "reserved_I4_leading_preserved")


def test_actual_i3_exit_primitive_memory_is_carried_as_one_over_X(candidate):
    eta = np.asarray([-0.33, 0.0, 0.47])
    seam_total = candidate.parent.similarity_profile_values(
        np.full(eta.shape, candidate.X_I3_end), eta
    )
    seam_base = candidate.current.similarity_profile_values(
        np.full(eta.shape, candidate.X_I3_end), eta
    )
    dm = np.asarray(seam_total["M_over_X_current_leading_through_I3"]) - np.asarray(
        seam_base["M_over_X_current_rf40_power_law"]
    )
    dm_eta = np.asarray(
        seam_total["M_eta_over_X_current_leading_through_I3"]
    ) - np.asarray(seam_base["M_eta_over_X_current_rf40_power_law"])

    X = np.full(eta.shape, _i4_midpoint(candidate))
    child = candidate.similarity_profile_values(X, eta)
    factor = candidate.X_I3_end / X
    np.testing.assert_allclose(
        child["delta_M_over_X_carried_from_I3_exit"],
        factor * dm,
        rtol=5e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        child["delta_M_eta_over_X_carried_from_I3_exit"],
        factor * dm_eta,
        rtol=5e-14,
        atol=0.0,
    )


def test_physical_M_is_constant_from_post_i3_gap_through_i4(candidate):
    eta = np.asarray(0.23)
    X0 = math.exp(0.5 * (candidate.log_X_I3_end + candidate.log_X_I4_start))
    X1 = _i4_midpoint(candidate)
    p0 = candidate.similarity_profile_values(np.asarray(X0), eta)
    p1 = candidate.similarity_profile_values(np.asarray(X1), eta)
    M0 = X0 * float(np.asarray(p0["M_over_X_current_leading_through_I4"]))
    M1 = X1 * float(np.asarray(p1["M_over_X_current_leading_through_I4"]))
    scale = max(1.0, abs(M0), abs(M1))
    assert abs(M1 - M0) / scale < 4e-12
    assert float(np.asarray(p0["U_current_leading_through_I4"])) == 0.0
    assert float(np.asarray(p1["U_current_leading_through_I4"])) == 0.0


def test_i4_radial_derivative_reuses_power_law_and_matches_centered_X(candidate):
    X = _i4_midpoint(candidate)
    eta = 0.19
    out = candidate.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    base = candidate.current.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    np.testing.assert_array_equal(
        out["F_current_leading_through_I4_X"], base["F_current_rf40_power_law_X"]
    )
    np.testing.assert_array_equal(
        out["U_current_leading_through_I4_X"], base["U_current_rf40_power_law_X"]
    )
    np.testing.assert_array_equal(
        out["E_current_leading_through_I4_X"], base["E_current_rf40_power_law_X"]
    )

    eps = 2e-6
    plus = candidate.similarity_profile_values(X * (1 + eps), eta)
    minus = candidate.similarity_profile_values(X * (1 - eps), eta)
    fd = (
        float(np.asarray(plus["F_current_leading_through_I4"]))
        - float(np.asarray(minus["F_current_leading_through_I4"]))
    ) / (2 * eps * X)
    analytic = float(np.asarray(out["F_current_leading_through_I4_X"]))
    assert abs(fd - analytic) / max(abs(analytic), 1e-300) < 2e-5


def test_vectorized_cartesian_velocity_reaches_i4_and_axis_remains_regular(candidate):
    X = _i4_midpoint(candidate)
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


def test_save_load_binds_parent_and_source_scope(candidate, tmp_path):
    path = tmp_path / "current_i4_leading.json"
    payload = candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianI4LeadingPreservation.load_configuration(path)
    assert restored.configuration() == payload == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["I4_role"] = "leading_patch"
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianI4LeadingPreservation.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["mutable"] = True
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianI4LeadingPreservation.from_configuration(mutated)


def test_fails_closed_after_i4_and_does_not_promote_mean_or_pde(candidate):
    with pytest.raises(ValueError, match="leading-through-I4 domain"):
        candidate.similarity_profile_values(
            math.nextafter(candidate.X_I4_end, math.inf), 0.0
        )

    truth = candidate.truth_boundary
    assert truth["source_I4_is_reserved_mean_correction_interval"] is True
    assert truth["source_leading_E0_U0_unchanged_on_I4"] is True
    assert truth["source_I4_is_not_a_leading_E0_U0_patch"] is True
    assert truth["current_leading_preserved_through_I4"] is True
    assert truth["leading_I4_overlay_invented"] is False
    for key in (
        "source_positive_order_I3_profiles_materialized",
        "I3_positive_order_correction_materialized",
        "source_I4_mean_correction_materialized",
        "current_I4_mean_correction_materialized",
        "source_terminal_tail_schedule_bound_into_current_velocity",
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
