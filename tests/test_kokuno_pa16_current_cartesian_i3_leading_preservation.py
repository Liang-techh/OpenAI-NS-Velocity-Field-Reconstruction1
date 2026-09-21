import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i3_leading_preservation import (
    PARENT_EXACT_HEAD,
    SOURCE_BLOB,
    KokunoPA16CurrentCartesianI3LeadingPreservation,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianI3LeadingPreservation()


def _i3_midpoint(candidate):
    return math.exp(0.5 * (candidate.log_X_I3_start + candidate.log_X_I3_end))


def test_exact_i2_parent_replay_through_i2_exit(candidate):
    eta = np.asarray([-0.37, 0.0, 0.51])
    X = np.asarray(
        [0.97 * candidate.X_I2_end, candidate.X_I2_end, candidate.X_I2_end]
    )
    child = candidate.similarity_profile_values(X, eta)
    parent = candidate.parent.similarity_profile_values(X, eta)
    mapping = {
        "F_current_leading_through_I3": "F_current_i2_repair",
        "U_current_leading_through_I3": "U_current_i2_repair",
        "E_current_leading_through_I3": "E_current_i2_repair",
        "M_over_X_current_leading_through_I3": "M_over_X_current_i2_repair",
        "M_eta_over_X_current_leading_through_I3": "M_eta_over_X_current_i2_repair",
        "v0_current_leading_through_I3": "v0_current_i2_repair",
    }
    for child_key, parent_key in mapping.items():
        np.testing.assert_array_equal(child[child_key], parent[parent_key])


def test_reserved_schedule_orders_i2_i3_i4_and_source_scope_is_frozen(candidate):
    assert candidate.log_X_I2_end < candidate.log_X_I3_start
    assert candidate.log_X_I3_start < candidate.log_X_I3_end
    assert candidate.log_X_I3_end < candidate.log_X_I4_start
    assert candidate.X_I3_end <= candidate.current.X_4
    cfg = candidate.configuration()["bound_source_scope"]
    assert cfg["source_blob"] == SOURCE_BLOB
    assert cfg["I3_role"] == "positive_order_not_leading"
    assert cfg["leading_action"] == "preserve_RF40_power_law"
    assert cfg["mutable"] is False
    assert PARENT_EXACT_HEAD == "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"


def test_leading_F_U_E_on_I3_are_exact_current_rf40_power_law(candidate):
    X = np.asarray([
        math.exp(candidate.log_X_I3_start + 0.1),
        _i3_midpoint(candidate),
        math.exp(candidate.log_X_I3_end - 0.1),
    ])
    eta = np.asarray([-0.41, 0.0, 0.58])
    child = candidate.similarity_profile_values(X, eta)
    base = candidate.current.similarity_profile_values(X, eta)
    np.testing.assert_array_equal(
        child["F_current_leading_through_I3"], base["F_current_rf40_power_law"]
    )
    np.testing.assert_array_equal(
        child["U_current_leading_through_I3"], base["U_current_rf40_power_law"]
    )
    np.testing.assert_array_equal(
        child["E_current_leading_through_I3"], base["E_current_rf40_power_law"]
    )
    assert np.all(child["U_current_leading_through_I3"] == 0.0)
    assert np.all(child["region"] == "reserved_I3_leading_preserved")


def test_actual_i2_exit_primitive_memory_is_carried_as_one_over_X(candidate):
    eta = np.asarray([-0.33, 0.0, 0.47])
    seam_total = candidate.parent.similarity_profile_values(
        np.full(eta.shape, candidate.X_I2_end), eta
    )
    seam_base = candidate.current.similarity_profile_values(
        np.full(eta.shape, candidate.X_I2_end), eta
    )
    dm = np.asarray(seam_total["M_over_X_current_i2_repair"]) - np.asarray(
        seam_base["M_over_X_current_rf40_power_law"]
    )
    dm_eta = np.asarray(seam_total["M_eta_over_X_current_i2_repair"]) - np.asarray(
        seam_base["M_eta_over_X_current_rf40_power_law"]
    )

    X = np.full(eta.shape, _i3_midpoint(candidate))
    child = candidate.similarity_profile_values(X, eta)
    factor = candidate.X_I2_end / X
    np.testing.assert_allclose(
        child["delta_M_over_X_carried_from_I2_exit"], factor * dm, rtol=5e-14, atol=0.0
    )
    np.testing.assert_allclose(
        child["delta_M_eta_over_X_carried_from_I2_exit"], factor * dm_eta, rtol=5e-14, atol=0.0
    )


def test_physical_M_is_constant_from_post_i2_gap_through_i3(candidate):
    eta = np.asarray(0.23)
    X0 = math.exp(0.5 * (candidate.log_X_I2_end + candidate.log_X_I3_start))
    X1 = _i3_midpoint(candidate)
    p0 = candidate.similarity_profile_values(np.asarray(X0), eta)
    p1 = candidate.similarity_profile_values(np.asarray(X1), eta)
    M0 = X0 * float(np.asarray(p0["M_over_X_current_leading_through_I3"]))
    M1 = X1 * float(np.asarray(p1["M_over_X_current_leading_through_I3"]))
    scale = max(1.0, abs(M0), abs(M1))
    assert abs(M1 - M0) / scale < 4e-12
    assert float(np.asarray(p0["U_current_leading_through_I3"])) == 0.0
    assert float(np.asarray(p1["U_current_leading_through_I3"])) == 0.0


def test_i3_radial_derivative_reuses_public_power_law_and_matches_centered_X(candidate):
    X = _i3_midpoint(candidate)
    eta = 0.19
    out = candidate.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    base = candidate.current.similarity_radial_derivatives(np.asarray(X), np.asarray(eta))
    np.testing.assert_array_equal(
        out["F_current_leading_through_I3_X"], base["F_current_rf40_power_law_X"]
    )
    np.testing.assert_array_equal(
        out["U_current_leading_through_I3_X"], base["U_current_rf40_power_law_X"]
    )
    np.testing.assert_array_equal(
        out["E_current_leading_through_I3_X"], base["E_current_rf40_power_law_X"]
    )

    eps = 2e-6
    plus = candidate.similarity_profile_values(X * (1 + eps), eta)
    minus = candidate.similarity_profile_values(X * (1 - eps), eta)
    fd = (
        float(np.asarray(plus["F_current_leading_through_I3"]))
        - float(np.asarray(minus["F_current_leading_through_I3"]))
    ) / (2 * eps * X)
    analytic = float(np.asarray(out["F_current_leading_through_I3_X"]))
    assert abs(fd - analytic) / max(abs(analytic), 1e-300) < 2e-5


def test_vectorized_cartesian_velocity_reaches_i3_and_axis_remains_regular(candidate):
    X = _i3_midpoint(candidate)
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
    path = tmp_path / "current_i3_leading.json"
    payload = candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianI3LeadingPreservation.load_configuration(path)
    assert restored.configuration() == payload == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["I3_role"] = "leading_patch"
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianI3LeadingPreservation.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_source_scope"]["mutable"] = True
    with pytest.raises(ValueError, match="source-scope"):
        KokunoPA16CurrentCartesianI3LeadingPreservation.from_configuration(mutated)


def test_fails_closed_after_i3_and_does_not_promote_positive_order_or_pde(candidate):
    with pytest.raises(ValueError, match="leading-through-I3 domain"):
        candidate.similarity_profile_values(
            math.nextafter(candidate.X_I3_end, math.inf), 0.0
        )

    truth = candidate.truth_boundary
    assert truth["source_I3_is_positive_order_patch"] is True
    assert truth["source_I3_is_not_a_leading_E0_U0_patch"] is True
    assert truth["current_leading_preserved_through_I3"] is True
    assert truth["leading_I3_overlay_invented"] is False
    for key in (
        "source_positive_order_I3_profiles_materialized",
        "I3_positive_order_correction_materialized",
        "current_I4_mean_correction_materialized",
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
