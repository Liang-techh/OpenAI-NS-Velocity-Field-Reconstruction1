import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_autonomous_modulation import (
    KokunoPA16CurrentCartesianAutonomousModulation,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianAutonomousModulation()


def test_support_is_current_power_stage_and_stops_at_i1(candidate):
    left, right = candidate.support_log_interval
    assert candidate.outer_base.log_X_w < left < right
    assert right < candidate.log_X_I1_start < candidate.current.log_X_4
    assert candidate.truth_boundary[
        "repository_autonomous_modulation_bound_to_current_lineage"
    ]
    assert not candidate.truth_boundary["source_admissible_loop_reconstructed"]
    assert not candidate.truth_boundary["I1_repair_applied_to_current_modulation"]


def test_inherited_current_candidate_replays_exactly_before_support(candidate):
    left, _ = candidate.support_log_interval
    y = left - 0.75
    X = math.exp(y)
    eta = 0.23
    new = candidate.similarity_profile_values(X, eta)
    old = candidate.current.similarity_profile_values(X, eta)

    pairs = (
        ("F_current_autonomous_modulation", "F_current_rf40_power_law"),
        ("U_current_autonomous_modulation", "U_current_rf40_power_law"),
        ("E_current_autonomous_modulation", "E_current_rf40_power_law"),
        ("M_over_X_current_autonomous_modulation", "M_over_X_current_rf40_power_law"),
        ("M_eta_over_X_current_autonomous_modulation", "M_eta_over_X_current_rf40_power_law"),
        ("v0_current_autonomous_modulation", "v0_current_rf40_power_law"),
    )
    for new_key, old_key in pairs:
        assert np.array_equal(np.asarray(new[new_key]), np.asarray(old[old_key]))
    assert float(np.asarray(new["delta_M_over_X_autonomous_modulation"])) == 0.0


def test_active_profile_reuses_existing_source_form_modulation(candidate):
    left, right = candidate.support_log_interval
    y = left + 0.371 * (right - left)
    X = math.exp(y)
    eta = 0.2
    new = candidate.similarity_profile_values(X, eta)
    selected = candidate.modulation.profile_values_logX(y, eta)
    base = candidate.current.similarity_profile_values(X, eta)

    assert new["region"].item() == "selected_autonomous_modulation"
    for new_key, selected_key in (
        ("F_current_autonomous_modulation", "F"),
        ("U_current_autonomous_modulation", "U"),
        ("E_current_autonomous_modulation", "E"),
    ):
        assert float(np.asarray(new[new_key])) == pytest.approx(
            float(np.asarray(selected[selected_key])), rel=2.0e-14, abs=0.0
        )

    # Nontriviality is structural; it is not obtained by shrinking the field.
    delta = abs(
        float(np.asarray(new["F_current_autonomous_modulation"]))
        - float(np.asarray(base["F_current_rf40_power_law"]))
    ) + abs(
        float(np.asarray(new["U_current_autonomous_modulation"]))
        - float(np.asarray(base["U_current_rf40_power_law"]))
    )
    assert delta > 0.0


def test_modulation_prefix_memory_persists_and_decays_after_support(candidate):
    _, right = candidate.support_log_interval
    eta = 0.27
    y1 = right + 0.8
    y2 = right + 1.8
    assert y2 < candidate.log_X_I1_start
    p1 = candidate.similarity_profile_values(math.exp(y1), eta)
    p2 = candidate.similarity_profile_values(math.exp(y2), eta)
    d1 = float(np.asarray(p1["delta_M_over_X_autonomous_modulation"]))
    d2 = float(np.asarray(p2["delta_M_over_X_autonomous_modulation"]))
    de1 = float(np.asarray(p1["delta_M_eta_over_X_autonomous_modulation"]))
    de2 = float(np.asarray(p2["delta_M_eta_over_X_autonomous_modulation"]))
    assert d1 != 0.0
    assert de1 != 0.0
    factor = math.exp(-(y2 - y1))
    assert d2 == pytest.approx(d1 * factor, rel=3.0e-11, abs=1.0e-300)
    assert de2 == pytest.approx(de1 * factor, rel=3.0e-11, abs=1.0e-300)


def test_current_primitive_obeys_logX_incompressibility_ode_inside_modulation(candidate):
    left, right = candidate.support_log_interval
    y = left + 0.43 * (right - left)
    eta = 0.31
    h = 2.0e-6
    centre = candidate.similarity_profile_values(math.exp(y), eta)
    plus = candidate.similarity_profile_values(math.exp(y + h), eta)
    minus = candidate.similarity_profile_values(math.exp(y - h), eta)
    m = float(np.asarray(centre["M_over_X_current_autonomous_modulation"]))
    u = float(np.asarray(centre["U_current_autonomous_modulation"]))
    fd = (
        float(np.asarray(plus["M_over_X_current_autonomous_modulation"]))
        - float(np.asarray(minus["M_over_X_current_autonomous_modulation"]))
    ) / (2.0 * h)
    assert fd == pytest.approx(u - m, rel=2.0e-5, abs=1.0e-40)


def test_radial_derivatives_replay_with_independent_centered_X_path(candidate):
    left, right = candidate.support_log_interval
    y = left + 0.57 * (right - left)
    X = math.exp(y)
    eta = -0.22
    h = 2.0e-6
    xp = X * math.exp(h)
    xm = X * math.exp(-h)
    deriv = candidate.similarity_radial_derivatives(X, eta)
    plus = candidate.similarity_profile_values(xp, eta)
    minus = candidate.similarity_profile_values(xm, eta)
    denominator = xp - xm
    for dkey, vkey in (
        ("F_current_autonomous_modulation_X", "F_current_autonomous_modulation"),
        ("U_current_autonomous_modulation_X", "U_current_autonomous_modulation"),
        ("E_current_autonomous_modulation_X", "E_current_autonomous_modulation"),
    ):
        fd = (
            float(np.asarray(plus[vkey])) - float(np.asarray(minus[vkey]))
        ) / denominator
        exact = float(np.asarray(deriv[dkey]))
        assert fd == pytest.approx(exact, rel=8.0e-5, abs=1.0e-300)


def test_vectorized_cartesian_interface_and_axis_regularity(candidate):
    axis = candidate.velocity(
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.2]),
        np.asarray([0.5, 0.5]),
    )
    assert axis.shape == (2, 3)
    assert np.all(np.isfinite(axis))
    assert np.array_equal(axis[:, :2], np.zeros((2, 2)))

    vals = candidate.velocity(
        np.asarray([0.05, 0.07]),
        np.asarray([0.02, -0.01]),
        np.asarray([0.0, 0.1]),
        np.asarray([0.2, 0.2]),
    )
    assert vals.shape == (2, 3)
    assert np.all(np.isfinite(vals))


def test_i1_and_later_domain_fails_closed(candidate):
    eta = 0.1
    with pytest.raises(ValueError):
        candidate.similarity_profile_values(
            math.exp(candidate.log_X_I1_start + 1.0e-5), eta
        )


def test_configuration_roundtrip_and_truth_boundary(candidate, tmp_path):
    path = tmp_path / "current_autonomous_modulation.json"
    candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianAutonomousModulation.load_configuration(path)
    assert restored.configuration() == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    payload = copy.deepcopy(candidate.configuration())
    payload["schema"] = "wrong"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianAutonomousModulation.from_configuration(payload)

    truth = candidate.truth_boundary
    assert truth["repository_autonomous_modulation_bound_to_current_lineage"] is True
    assert truth["cone_modulation_completed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
