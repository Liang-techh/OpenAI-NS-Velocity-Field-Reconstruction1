import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_heat_exterior_profile import (
    KokunoHeatExteriorProfile,
)


def test_source_endpoint_derivatives_and_heat_ode():
    profile = KokunoHeatExteriorProfile()
    h = profile.h

    np.testing.assert_allclose(profile.heat_factor(0.0), 1.0, rtol=0.0, atol=3e-14)
    np.testing.assert_allclose(
        profile.heat_factor(0.0, 1), -h * (1.0 + h), rtol=0.0, atol=3e-14
    )
    np.testing.assert_allclose(
        profile.heat_factor(0.0, 2),
        h * (1.0 + h) ** 2 * (2.0 + h),
        rtol=0.0,
        atol=5e-14,
    )

    Z = np.array([0.0, 0.1, 1.0, 5.0, 10.0])
    residual = profile.heat_ode_residual(Z)
    assert np.max(np.abs(residual)) < 2e-7
    assert np.all(profile.heat_factor(Z) > 0.0)


def test_profile_contract_and_partials_match_finite_difference():
    profile = KokunoHeatExteriorProfile()
    X = 3.0
    eta = 0.35
    values = profile.profile_values(X, eta)

    np.testing.assert_allclose(values["E"], math.sqrt(2.0 * X) * values["F"])
    np.testing.assert_allclose(values["Pi_X"], values["F"] ** 2)
    np.testing.assert_array_equal(values["U"], 0.0)
    np.testing.assert_array_equal(values["v0"], 0.0)

    eps = 1e-5
    F_x_fd = (
        profile.profile_values(X + eps, eta)["F"]
        - profile.profile_values(X - eps, eta)["F"]
    ) / (2.0 * eps)
    F_eta_fd = (
        profile.profile_values(X, eta + eps)["F"]
        - profile.profile_values(X, eta - eps)["F"]
    ) / (2.0 * eps)
    np.testing.assert_allclose(values["F_X"], F_x_fd, rtol=2e-8, atol=2e-11)
    np.testing.assert_allclose(values["F_eta"], F_eta_fd, rtol=2e-8, atol=2e-11)


def test_native_similarity_velocity_replays_direct_physical_heat_field():
    profile = KokunoHeatExteriorProfile()
    x = np.array([2.0, 2.0, 2.0])
    y = np.array([0.75, 0.75, 0.75])
    z = np.array([-0.5, 0.0, 0.5])
    t = 0.4

    native = profile.velocity(x, y, z, t)
    direct = profile.direct_heat_velocity(x, y, t)
    np.testing.assert_allclose(native, direct, rtol=5e-12, atol=5e-13)
    np.testing.assert_allclose(
        native, np.broadcast_to(native[1], native.shape), rtol=5e-12, atol=5e-13
    )
    assert np.linalg.norm(native[0]) > 0.1

    points = np.stack((x, y, z), axis=-1)
    np.testing.assert_allclose(
        profile.at_points(points, t), native, rtol=0.0, atol=0.0
    )


def test_domain_parameter_and_nontriviality_guards():
    profile = KokunoHeatExteriorProfile()
    with pytest.raises(ValueError, match="X > 0"):
        profile.profile_values(0.0, 0.0)
    with pytest.raises(ValueError, match="not an axis field"):
        profile.velocity(0.0, 0.0, 0.0, 0.5)
    with pytest.raises(ValueError, match="Z >= 0"):
        profile.heat_factor(-0.1)
    with pytest.raises(ValueError, match="0 < c_inf"):
        KokunoHeatExteriorProfile(c_inf=0.0)
    with pytest.raises(ValueError, match="0 < h < 1e-2"):
        KokunoHeatExteriorProfile(h=0.01)
    with pytest.raises(TypeError, match="quadrature_order"):
        KokunoHeatExteriorProfile(quadrature_order=32.0)


def test_fail_closed_json_roundtrip(tmp_path):
    profile = KokunoHeatExteriorProfile(h=0.004, c_inf=1.25, quadrature_order=64)
    path = profile.save_json(tmp_path / "heat.json")
    loaded = KokunoHeatExteriorProfile.load_json(path)
    assert loaded.sha256 == profile.sha256

    points = np.array([[2.0, 0.5, -0.25], [2.25, -0.4, 0.3]])
    np.testing.assert_allclose(
        loaded.at_points(points, 0.45), profile.at_points(points, 0.45)
    )

    payload = profile.to_payload()
    bad = copy.deepcopy(payload)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoHeatExteriorProfile.from_payload(bad)

    bad = copy.deepcopy(payload)
    bad["parameters"]["c_inf_origin"] = "recovered_hidden_parameter"
    with pytest.raises(ValueError, match="amplitude provenance"):
        KokunoHeatExteriorProfile.from_payload(bad)
