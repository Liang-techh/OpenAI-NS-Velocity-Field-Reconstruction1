import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_reference_continuation import source_smooth_step
from openai_ns_reconstruction.kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)


@pytest.fixture(scope="module")
def candidate():
    # Eight-point source-y quadrature is enough for the regression guards and
    # keeps this stacked constrained lane from re-solving the same smooth
    # transition at production order in every test.
    return KokunoSourceRescaledReferenceContinuation(quadrature_points=8)


def test_large_pressure_reference_transition_fits_source_core_and_freezes(candidate):
    assert candidate.pressure_scale > 1.0e13
    assert candidate.rescaling_lambda > 1.0e26
    assert candidate.X2 < candidate.max_source_transition_X

    eta = candidate.seed.phase_stationary_eta
    natural_X = 0.75 * candidate.X1
    natural = candidate._continued_scalar(natural_X, eta)
    seed = candidate.seed.profile_values(natural_X, eta)
    assert natural["log_F"] == pytest.approx(float(seed["log_F"]), abs=1.0e-13)
    assert natural["F"] == pytest.approx(float(seed["F"]), rel=2.0e-15)
    assert natural["U"] == pytest.approx(float(seed["U"]), rel=2.0e-15)

    at_exit = candidate._continued_scalar(candidate.X2, eta)
    after_exit = candidate._continued_scalar(3.0 * candidate.X2, eta)
    assert after_exit["log_F"] == pytest.approx(at_exit["log_F"], abs=2.0e-14)
    assert after_exit["U"] == pytest.approx(at_exit["U"], rel=2.0e-14, abs=1.0e-14)
    assert after_exit["D_X_log_F"] == pytest.approx(0.0, abs=0.0)
    assert after_exit["D_X_U"] == pytest.approx(0.0, abs=0.0)


def test_transition_radial_slopes_follow_public_flat_step(candidate):
    eta = 0.25
    y = 1.5 * candidate.log_transition_width
    X = candidate.X0 * math.exp(y)
    values = candidate._continued_scalar(X, eta)
    # At the selected source pressure scale X0 is ~1e-27.  Recover the source
    # logarithmic coordinate from the represented X, exactly as the public
    # formula requires, rather than comparing against the pre-roundtrip y.
    source_y = math.log(X / candidate.X0)
    s = (source_y - candidate.log_transition_width) / candidate.log_transition_width
    gate = 1.0 - float(source_smooth_step(np.asarray(s)))
    expected_log_slope = gate * float(candidate._natural_log_slope(np.asarray(X), eta))
    expected_DU = gate * float(candidate._natural_DU(np.asarray(X), eta))
    assert values["D_X_log_F"] == pytest.approx(expected_log_slope, rel=3.0e-14, abs=1.0e-15)
    assert values["D_X_U"] == pytest.approx(expected_DU, rel=3.0e-14, abs=1.0e-15)

    eps = 2.0e-6
    plus = candidate._continued_scalar(X * math.exp(eps), eta)
    minus = candidate._continued_scalar(X * math.exp(-eps), eta)
    finite_DU = (plus["U"] - minus["U"]) / (2.0 * eps)
    assert values["D_X_U"] == pytest.approx(finite_DU, rel=2.0e-6, abs=2.0e-10)


def test_frozen_reference_carries_prefix_velocity_and_pressure_identity(candidate):
    eta = candidate.seed.phase_stationary_eta
    values = candidate.profile_values(110.0, eta)
    for key in ("F", "U", "v0", "Pi", "M", "Pi_X", "log_F"):
        assert np.all(np.isfinite(values[key])), key
    assert float(values["Pi_X"]) == pytest.approx(float(values["F"]) ** 2, abs=0.0)
    assert abs(float(values["U"])) > 1.0e-3
    assert abs(float(values["v0"])) > 1.0e-3
    assert float(values["Pi"]) < 0.0


def test_native_velocity_is_nontrivial_and_axis_regular(candidate):
    axis = candidate.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis.shape == (3,)
    assert np.all(np.isfinite(axis))
    assert axis[0] == pytest.approx(0.0, abs=0.0)
    assert axis[1] == pytest.approx(0.0, abs=0.0)
    assert abs(axis[2]) > 1.0e-3

    eta = candidate.seed.phase_stationary_eta
    q = 1.0 / (1.0 - eta * eta)
    z = q ** candidate.D * eta
    r = math.sqrt(2.0 * q * 110.0)
    velocity = candidate.velocity(r, 0.0, z, 0.0)
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity) > 1.0e-3


def test_serialization_report_and_truth_boundary_fail_closed(candidate, tmp_path):
    payload = candidate.to_payload()
    replay = KokunoSourceRescaledReferenceContinuation.from_payload(payload)
    assert replay.sha256 == candidate.sha256

    report = candidate.report()
    truth = report["truth_boundary"]
    assert report["all_sample_fields_finite"] is True
    assert report["source_transition_fits_rescaled_core_domain"] is True
    assert truth["selected_source_pressure_datum_carried_into_reference_continuation"] is True
    assert truth["source_reference_continuation_executable_at_selected_pressure_scale"] is True
    assert truth["source_fixed_point_solved"] is False
    assert truth["source_appendix_B_activation_executed_at_selected_pressure_scale"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = tmp_path / "rescaled-reference.json"
    candidate.save_json(path)
    assert KokunoSourceRescaledReferenceContinuation.load_json(path).sha256 == candidate.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_fixed_point_solved"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoSourceRescaledReferenceContinuation.from_payload(tampered)
