import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i1_repair import (
    KokunoPA16CurrentCartesianI1Repair,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianI1Repair()


def _assert_close(a, b, *, rtol=2.0e-11, atol=2.0e-12):
    np.testing.assert_allclose(
        np.asarray(a, dtype=float),
        np.asarray(b, dtype=float),
        rtol=rtol,
        atol=atol,
    )


def test_current_modulation_identity_and_pa17_closure(candidate):
    assert candidate.i1_family.modulation.sha256 == candidate.pre_i1.modulation.sha256
    report = candidate.i1_family.closure_report(np.asarray([0.0, 0.2, 0.6, 1.0]))
    assert report["passed"] is True
    assert report["max_abs_residual"] <= candidate.closure_tolerance
    family = candidate.i1_family.family_report()
    assert family["max_abs_coefficient"] <= family["coefficient_limit"] + 1.0e-14


def test_exact_pre_i1_replay_and_entry_handoff(candidate):
    _, support_right = candidate.pre_i1.support_log_interval
    log_x = 0.5 * (support_right + candidate.log_X_I1_start)
    X = math.exp(log_x)
    eta = np.asarray([-0.55, 0.0, 0.35])
    new = candidate.similarity_profile_values(X, eta)
    old = candidate.pre_i1.similarity_profile_values(X, eta)
    for new_key, old_key in (
        ("F_current_i1_repair", "F_current_autonomous_modulation"),
        ("U_current_i1_repair", "U_current_autonomous_modulation"),
        ("E_current_i1_repair", "E_current_autonomous_modulation"),
        ("M_over_X_current_i1_repair", "M_over_X_current_autonomous_modulation"),
        (
            "M_eta_over_X_current_i1_repair",
            "M_eta_over_X_current_autonomous_modulation",
        ),
        ("v0_current_i1_repair", "v0_current_autonomous_modulation"),
    ):
        _assert_close(new[new_key], old[old_key], rtol=0.0, atol=0.0)

    entry_new = candidate.similarity_profile_values(candidate.X_I1_start, eta)
    entry_old = candidate.pre_i1.similarity_profile_values(
        candidate.X_I1_start, eta
    )
    for new_key, old_key in (
        ("F_current_i1_repair", "F_current_autonomous_modulation"),
        ("U_current_i1_repair", "U_current_autonomous_modulation"),
        ("E_current_i1_repair", "E_current_autonomous_modulation"),
        ("M_over_X_current_i1_repair", "M_over_X_current_autonomous_modulation"),
        ("v0_current_i1_repair", "v0_current_autonomous_modulation"),
    ):
        _assert_close(entry_new[new_key], entry_old[old_key], rtol=0.0, atol=0.0)


def test_active_i1_profile_is_base_plus_repair_and_current_memory(candidate):
    X = candidate.i1_family.repair.X_1
    eta = np.asarray([-0.4, 0.25, 0.7])
    got = candidate.similarity_profile_values(X, eta)
    base = candidate.current.similarity_profile_values(X, eta)
    repair = candidate.i1_family.profile_correction_logX(math.log(X), eta)
    dm, dm_eta = candidate._carried_modulation_memory_ratio(
        np.full(eta.shape, X), eta
    )

    _assert_close(
        got["F_current_i1_repair"],
        np.asarray(base["F_current_rf40_power_law"]) + repair["delta_F"],
    )
    _assert_close(
        got["U_current_i1_repair"],
        np.asarray(base["U_current_rf40_power_law"]) + repair["delta_U"],
    )
    _assert_close(
        got["E_current_i1_repair"],
        np.asarray(base["E_current_rf40_power_law"]) + repair["delta_E"],
    )
    expected_m = (
        np.asarray(base["M_over_X_current_rf40_power_law"])
        + dm
        + np.asarray(repair["delta_M"]) / X
    )
    expected_m_eta = (
        np.asarray(base["M_eta_over_X_current_rf40_power_law"])
        + dm_eta
        + np.asarray(repair["delta_M_eta"]) / X
    )
    _assert_close(got["M_over_X_current_i1_repair"], expected_m)
    _assert_close(got["M_eta_over_X_current_i1_repair"], expected_m_eta)
    assert np.max(np.abs(np.asarray(repair["delta_U"]))) > 0.0
    assert np.all(np.asarray(got["F_current_i1_repair"]) > 0.0)
    assert np.all(np.asarray(got["E_current_i1_repair"]) >= 0.0)


def test_current_prefix_M_obeys_dM_dX_equals_U_inside_i1(candidate):
    X = candidate.i1_family.repair.X_1
    eta = 0.31
    eps = 2.0e-6
    Xm = X * (1.0 - eps)
    Xp = X * (1.0 + eps)
    pm = candidate.similarity_profile_values(Xm, eta)
    pp = candidate.similarity_profile_values(Xp, eta)
    p0 = candidate.similarity_profile_values(X, eta)
    Mm = Xm * float(np.asarray(pm["M_over_X_current_i1_repair"]))
    Mp = Xp * float(np.asarray(pp["M_over_X_current_i1_repair"]))
    dM = (Mp - Mm) / (Xp - Xm)
    U = float(np.asarray(p0["U_current_i1_repair"]))
    assert dM == pytest.approx(U, rel=2.0e-5, abs=2.0e-7)


def test_radial_derivative_replay_inside_i1(candidate):
    X = candidate.i1_family.repair.X_1
    eta = -0.28
    eps = 2.0e-6
    pm = candidate.similarity_profile_values(X * (1.0 - eps), eta)
    pp = candidate.similarity_profile_values(X * (1.0 + eps), eta)
    fd = (
        float(np.asarray(pp["F_current_i1_repair"]))
        - float(np.asarray(pm["F_current_i1_repair"]))
    ) / (2.0 * eps * X)
    analytic = candidate.similarity_radial_derivatives(X, eta)
    fx = float(np.asarray(analytic["F_current_i1_repair_X"]))
    assert fd == pytest.approx(fx, rel=3.0e-5, abs=1.0e-18)


def test_exit_report_carries_actual_current_memory(candidate):
    eta = np.asarray([0.0, 0.35, 0.8])
    report = candidate.i1_exit_report(eta)
    assert report["family_closure"]["passed"] is True
    incoming = np.asarray(
        report["incoming_autonomous_M_over_X_at_I1_start"], dtype=float
    )
    outgoing = np.asarray(
        report["outgoing_overlay_M_over_X_at_I1_exit"], dtype=float
    )
    assert np.all(np.isfinite(incoming))
    assert np.all(np.isfinite(outgoing))
    assert np.max(np.abs(incoming)) > 0.0
    assert report["independent_PDE_validation"] is False


def test_cartesian_vectorization_axis_regularity_and_domain_guard(candidate):
    q = 0.5
    X_active = candidate.i1_family.repair.X_1
    X_pre = math.exp(
        0.5
        * (
            candidate.pre_i1.support_log_interval[1]
            + candidate.log_X_I1_start
        )
    )
    radius = np.sqrt(2.0 * q * np.asarray([X_pre, X_active]))
    velocity = candidate.velocity(radius, np.zeros(2), np.zeros(2), np.full(2, 0.5))
    assert velocity.shape == (2, 3)
    assert np.all(np.isfinite(velocity))

    axis = candidate.velocity(0.0, 0.0, 0.0, 0.5)
    assert float(np.asarray(axis)[0]) == 0.0
    assert float(np.asarray(axis)[1]) == 0.0
    assert np.all(np.isfinite(axis))

    with pytest.raises(ValueError):
        candidate.similarity_profile_values(candidate.X_I1_end * 1.000001, 0.0)


def test_configuration_roundtrip_semantic_dependency_and_truth(candidate, tmp_path, monkeypatch):
    path = tmp_path / "current_i1_repair.json"
    candidate.save_configuration(path)
    restored = KokunoPA16CurrentCartesianI1Repair.load_configuration(path)
    assert restored.configuration() == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256

    original_semantic = candidate.semantic_sha256
    modulation_type = type(candidate.pre_i1.modulation)
    monkeypatch.setattr(
        modulation_type,
        "sha256",
        property(lambda self: "f" * 64),
    )
    assert candidate.semantic_sha256 != original_semantic

    truth = candidate.truth_boundary
    assert truth["I1_repair_applied_to_current_modulation"] is True
    assert truth[
        "child_semantic_identity_binds_parent_and_modulation_dependency"
    ] is True
    assert truth[
        "parent_pre_I1_semantic_dependency_binding_repaired_by_this_increment"
    ] is False
    for key in (
        "source_hidden_loop_parameters_recovered",
        "I2_overlay_applied_to_current_lineage",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        assert truth[key] is False

    payload = copy.deepcopy(candidate.configuration())
    payload["schema"] = "wrong"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianI1Repair.from_configuration(payload)
