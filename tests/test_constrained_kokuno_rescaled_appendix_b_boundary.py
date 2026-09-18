import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rescaled_appendix_b_boundary import (
    A_FINAL,
    X_I,
    KokunoSourceRescaledAppendixBBoundary,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoSourceRescaledAppendixBBoundary(max_step=0.5)


def test_selected_pressure_path_hands_off_exactly_from_rescaled_reference(candidate):
    eta = candidate.reference.seed.phase_stationary_eta
    selected = candidate.profile_values(candidate.X0, eta)
    reference = candidate.reference.profile_values(candidate.X0, eta)
    for key in ("log_F", "F", "U", "v0", "Pi", "M", "D_X_log_F", "D_X_U"):
        assert float(selected[key]) == pytest.approx(
            float(reference[key]), rel=2.0e-13, abs=2.0e-13
        )
    assert candidate.reference.pressure_scale > 1.0e13
    assert candidate.reference.rescaling_lambda > 1.0e26


def test_xi_boundary_is_log_stable_and_reaches_public_final_slopes(candidate):
    eta = candidate.reference.seed.phase_stationary_eta
    boundary = candidate.boundary_values(np.asarray([eta]))
    for value in boundary.values():
        assert np.all(np.isfinite(value))
    assert float(boundary["radial_log_slope_F_i"][0]) == pytest.approx(
        -0.5 * A_FINAL, abs=2.0e-12
    )
    assert float(boundary["radial_log_slope_U_i"][0]) == pytest.approx(0.0, abs=1.0e-12)
    assert abs(float(boundary["G_i"][0])) > 1.0e-6
    assert candidate.log_C > 1.0e10
    reconstructed = (
        candidate.log_C
        + 0.5 * math.log(2.0 * X_I)
        + float(boundary["log_F_i"][0])
    )
    assert float(boundary["ell_i"][0]) == pytest.approx(reconstructed, rel=0.0, abs=1.0e-8)


def test_selected_profile_carries_pressure_identity_and_prefix_incompressibility(candidate):
    eta = candidate.reference.seed.phase_stationary_eta
    X = 10.0
    values = candidate.profile_values(X, eta)
    assert float(values["Pi_X"]) == pytest.approx(
        float(values["F"]) ** 2, rel=2.0e-15, abs=0.0
    )
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * candidate.h * eta * eta
    expected_v0 = (
        2.0 * eta * float(values["U"])
        - 2.0 * candidate.D * eta * float(values["M"]) / X
        - d * float(values["M_eta"]) / X
    ) / L
    assert float(values["v0"]) == pytest.approx(expected_v0, rel=5.0e-13, abs=5.0e-13)
    assert np.all(
        np.isfinite(
            [
                float(values["log_F"]),
                float(values["U"]),
                float(values["v0"]),
                float(values["Pi"]),
                float(values["p1_reference"]),
                float(values["n_s_reference"]),
            ]
        )
    )


def test_native_velocity_is_nontrivial_on_selected_activation_path(candidate):
    axis = candidate.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis.shape == (3,)
    assert np.all(np.isfinite(axis))
    assert axis[0] == pytest.approx(0.0, abs=0.0)
    assert axis[1] == pytest.approx(0.0, abs=0.0)

    eta = candidate.reference.seed.phase_stationary_eta
    X = 10.0
    q = 1.0 / (1.0 - eta * eta)
    z = q ** candidate.D * eta
    r = math.sqrt(2.0 * q * X)
    velocity = candidate.velocity(r, 0.0, z, 0.0)
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity) > 1.0e-6


def test_serialization_and_truth_boundary_fail_closed(candidate, tmp_path):
    payload = candidate.to_payload()
    replay = KokunoSourceRescaledAppendixBBoundary.from_payload(payload)
    assert replay.sha256 == candidate.sha256

    truth = payload["truth_boundary"]
    assert truth["selected_pressure_appendix_B_activation_executable"] is True
    assert truth["selected_pressure_appendix_B_Xi_profile_executable"] is True
    assert truth["source_fixed_point_solved"] is False
    assert truth["actual_source_appendix_B_trajectory_reconstructed"] is False
    assert truth["actual_source_incoming_five_moments_bound"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = tmp_path / "rescaled-appendix-b.json"
    candidate.save_json(path)
    assert KokunoSourceRescaledAppendixBBoundary.load_json(path).sha256 == candidate.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_fixed_point_solved"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoSourceRescaledAppendixBBoundary.from_payload(tampered)
