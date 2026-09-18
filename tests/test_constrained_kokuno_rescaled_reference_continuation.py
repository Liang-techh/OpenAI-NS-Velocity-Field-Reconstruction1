import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)


def test_large_pressure_reference_transition_fits_source_core_and_freezes():
    candidate = KokunoSourceRescaledReferenceContinuation()
    assert candidate.pressure_scale > 1.0e13
    assert candidate.rescaling_lambda > 1.0e26
    assert candidate.X2 < candidate.max_source_transition_X

    eta = candidate.seed.phase_stationary_eta
    natural_X = 0.75 * candidate.X1
    natural = candidate.profile_values(natural_X, eta)
    seed = candidate.seed.profile_values(natural_X, eta)
    assert float(natural["log_F"]) == pytest.approx(float(seed["log_F"]), abs=1.0e-13)
    assert float(natural["F"]) == pytest.approx(float(seed["F"]), rel=2.0e-15)
    assert float(natural["U"]) == pytest.approx(float(seed["U"]), rel=2.0e-15)

    at_exit = candidate.profile_values(candidate.X2, eta)
    after_exit = candidate.profile_values(3.0 * candidate.X2, eta)
    assert float(after_exit["log_F"]) == pytest.approx(float(at_exit["log_F"]), abs=2.0e-14)
    assert float(after_exit["U"]) == pytest.approx(float(at_exit["U"]), rel=2.0e-14, abs=1.0e-14)
    assert float(after_exit["D_X_log_F"]) == pytest.approx(0.0, abs=0.0)
    assert float(after_exit["D_X_U"]) == pytest.approx(0.0, abs=0.0)


def test_transition_radial_slopes_follow_public_flat_step_without_raw_taylor_state():
    candidate = KokunoSourceRescaledReferenceContinuation()
    eta = 0.25
    y = 1.5 * candidate.log_transition_width
    X = candidate.X0 * math.exp(y)
    values = candidate.profile_values(X, eta)
    s = (y - candidate.log_transition_width) / candidate.log_transition_width

    from openai_ns_reconstruction.kokuno_reference_continuation import source_smooth_step

    gate = 1.0 - float(source_smooth_step(np.asarray(s)))
    expected_log_slope = gate * float(candidate._natural_log_slope(np.asarray(X), eta))
    expected_DU = gate * float(candidate._natural_DU(np.asarray(X), eta))
    assert float(values["D_X_log_F"]) == pytest.approx(expected_log_slope, rel=3.0e-14, abs=1.0e-15)
    assert float(values["D_X_U"]) == pytest.approx(expected_DU, rel=3.0e-14, abs=1.0e-15)

    eps = 2.0e-6
    plus = candidate.profile_values(X * math.exp(eps), eta)
    minus = candidate.profile_values(X * math.exp(-eps), eta)
    finite_DU = (float(plus["U"]) - float(minus["U"])) / (2.0 * eps)
    assert float(values["D_X_U"]) == pytest.approx(finite_DU, rel=2.0e-6, abs=2.0e-10)


def test_global_frozen_reference_has_prefix_velocity_and_pressure_identity():
    candidate = KokunoSourceRescaledReferenceContinuation()
    eta = candidate.seed.phase_stationary_eta
    X = np.array([candidate.X2, 1.0, 110.0])
    values = candidate.profile_values(X, np.full(X.shape, eta))
    for key in ("F", "U", "v0", "Pi", "M", "Pi_X", "log_F"):
        assert np.all(np.isfinite(values[key])), key
    assert np.allclose(values["Pi_X"], values["F"] ** 2, rtol=0.0, atol=0.0)
    assert np.max(np.abs(values["U"])) > 1.0e-3
    assert np.max(np.abs(values["v0"])) > 1.0e-3
    assert np.all(values["Pi"] < 0.0)


def test_native_velocity_is_vectorized_nontrivial_and_axis_regular():
    candidate = KokunoSourceRescaledReferenceContinuation()
    axis = candidate.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis.shape == (3,)
    assert np.all(np.isfinite(axis))
    assert axis[0] == pytest.approx(0.0, abs=0.0)
    assert axis[1] == pytest.approx(0.0, abs=0.0)
    assert abs(axis[2]) > 1.0e-3

    eta = candidate.seed.phase_stationary_eta
    q = 1.0 / (1.0 - eta * eta)
    z = q ** candidate.D * eta
    r1 = math.sqrt(2.0 * q * 1.0)
    r2 = math.sqrt(2.0 * q * 110.0)
    velocity = candidate.velocity(np.array([r1, r2]), 0.0, z, 0.0)
    assert velocity.shape == (2, 3)
    assert np.all(np.isfinite(velocity))
    assert np.all(np.linalg.norm(velocity, axis=-1) > 1.0e-3)


def test_serialization_replays_and_truth_boundary_fails_closed(tmp_path):
    candidate = KokunoSourceRescaledReferenceContinuation()
    payload = candidate.to_payload()
    replay = KokunoSourceRescaledReferenceContinuation.from_payload(payload)
    assert replay.sha256 == candidate.sha256
    assert replay.report()["all_sample_fields_finite"] is True

    path = tmp_path / "rescaled-reference.json"
    candidate.save_json(path)
    assert KokunoSourceRescaledReferenceContinuation.load_json(path).sha256 == candidate.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_fixed_point_solved"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoSourceRescaledReferenceContinuation.from_payload(tampered)


def test_report_keeps_selected_stage_separate_from_global_source_claims():
    candidate = KokunoSourceRescaledReferenceContinuation()
    report = candidate.report()
    truth = report["truth_boundary"]
    assert report["source_transition_fits_rescaled_core_domain"] is True
    assert truth["selected_source_pressure_datum_carried_into_reference_continuation"] is True
    assert truth["source_reference_continuation_executable_at_selected_pressure_scale"] is True
    assert truth["source_fixed_point_solved"] is False
    assert truth["source_appendix_B_activation_executed_at_selected_pressure_scale"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
