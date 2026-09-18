from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_bounded_family_coefficients import (
    KokunoBoundedFamilyCoefficientCoordinates,
)
from openai_ns_reconstruction.kokuno_source_localized_real_pair_family import (
    KokunoSourceLocalizedRealPairFamily,
)


def _direct_family():
    velocity = np.array(
        [
            [[0.10, 0.20, -0.05], [0.03, -0.04, 0.07], [-0.02, 0.06, 0.08], [0.01, -0.03, 0.02]],
            [[0.12, 0.18, -0.04], [0.02, -0.05, 0.05], [-0.01, 0.04, 0.09], [0.03, -0.02, 0.01]],
        ],
        dtype=float,
    )
    labels = [(5, (0, 0, 0)), (5, (1, 0, 0)), (6, (0, 0, 0)), (6, (1, 0, 0))]
    theta = np.array([0.2, -0.4])
    return velocity, labels, theta


def _agent2_physical_family():
    R = np.array([0.82, 1.01, 1.19])
    theta = np.array([0.17, -0.31, 0.42])
    alpha = 0.23 + 0.15 * R
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.15 * tangent
    D_z_eta = np.zeros_like(eta)
    phase = R[:, None] + np.array([0.13, -0.27])[None, :]
    n_phi = np.zeros((3, 2, 3), dtype=float)
    n_phi[..., 0] = 1.0
    t_plus = np.empty((3, 2, 3), dtype=np.complex128)
    t_plus[:, 0, :] = np.array([0.0, 0.24 + 0.07j, -0.11 + 0.04j])
    t_plus[:, 1, :] = np.array([0.0, -0.18 + 0.05j, 0.16 - 0.03j])
    D_r_C_plus = np.zeros_like(t_plus)
    D_z_C_plus = np.zeros_like(t_plus)
    Q = np.array([[0.5, 0.25], [0.5, 0.25], [0.5, 0.25]])
    labels = [(5, (0, 0, 0)), (6, (1, 0, 0))]
    family = KokunoSourceLocalizedRealPairFamily(epsilon=0.2, h=0.004)
    result = family.physical_family(
        Q,
        R,
        theta,
        phase,
        n_phi,
        t_plus,
        D_r_C_plus,
        D_z_C_plus,
        eta,
        D_r_eta,
        D_z_eta,
        labels,
    )
    return result, theta


def test_two_coordinate_family_is_exact_and_label_count_independent():
    velocity, labels, theta = _direct_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.125)
    out = contract.evaluate(
        velocity,
        labels,
        delta_common=0.03,
        delta_band=-0.02,
        theta=theta,
    )

    np.testing.assert_allclose(out["band_contrast"], [-1.0, -1.0, 1.0, 1.0])
    assert out["parameter_count"] == 2
    assert out["parameter_count"] < len(labels)
    assert out["minimum_label_multiplier"] >= 1.0 - contract.max_l1_update

    g = out["band_contrast"]
    multipliers = 1.0 + 0.03 - 0.02 * g
    manual = np.sum(velocity * multipliers[None, :, None], axis=1)
    np.testing.assert_allclose(out["velocity_physical_cylindrical_total_modulated"], manual)

    expected_common = np.sum(velocity, axis=1)
    expected_band = np.sum(velocity * g[None, :, None], axis=1)
    np.testing.assert_allclose(
        out["velocity_parameter_tangent_cylindrical"][:, 0, :], expected_common
    )
    np.testing.assert_allclose(
        out["velocity_parameter_tangent_cylindrical"][:, 1, :], expected_band
    )
    assert out["velocity_parameter_tangent_cartesian"].shape == (2, 2, 3)
    assert out["repository_coefficient_unit_mapping_available"] is True
    assert out["source_coefficient_unit_mapping_available"] is False


def test_parameter_tangents_match_independent_centered_difference():
    velocity, labels, _ = _direct_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.2)
    h = 1.0e-6
    base = contract.evaluate(velocity, labels)
    tangents = base["velocity_parameter_tangent_cylindrical"]

    plus = contract.evaluate(velocity, labels, delta_common=h)[
        "velocity_physical_cylindrical_total_modulated"
    ]
    minus = contract.evaluate(velocity, labels, delta_common=-h)[
        "velocity_physical_cylindrical_total_modulated"
    ]
    fd_common = (plus - minus) / (2.0 * h)
    np.testing.assert_allclose(fd_common, tangents[:, 0, :], rtol=2e-10, atol=2e-11)

    plus = contract.evaluate(velocity, labels, delta_band=h)[
        "velocity_physical_cylindrical_total_modulated"
    ]
    minus = contract.evaluate(velocity, labels, delta_band=-h)[
        "velocity_physical_cylindrical_total_modulated"
    ]
    fd_band = (plus - minus) / (2.0 * h)
    np.testing.assert_allclose(fd_band, tangents[:, 1, :], rtol=2e-10, atol=2e-11)


def test_consumes_agent2_physical_family_without_recomputing_curl():
    physical, theta = _agent2_physical_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.125)
    out = contract.consume_physical_family(physical, theta=theta)
    np.testing.assert_allclose(
        out["velocity_physical_cylindrical_total_base"],
        physical["velocity_physical_cylindrical_total"],
        rtol=0.0,
        atol=3e-12,
    )
    np.testing.assert_allclose(
        out["velocity_physical_cartesian_total_base"],
        physical["velocity_physical_cartesian_total"],
        rtol=0.0,
        atol=3e-12,
    )
    assert out["consumed_agent2_physical_family_contract"] is True
    assert out["band_contrast_active"] is True
    assert np.linalg.norm(out["velocity_parameter_tangent_cylindrical"][:, 1, :]) > 0.0


def test_aggregate_bound_and_single_band_fail_closed():
    velocity, labels, _ = _direct_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.1)
    with pytest.raises(ValueError, match="aggregate L1"):
        contract.evaluate(velocity, labels, delta_common=0.07, delta_band=0.04)

    one_band = [(5, (i, 0, 0)) for i in range(4)]
    out = contract.evaluate(velocity, one_band)
    assert out["band_contrast_active"] is False
    np.testing.assert_array_equal(out["band_contrast"], np.zeros(4))
    np.testing.assert_array_equal(out["velocity_parameter_tangent_cylindrical"][:, 1, :], 0.0)
    with pytest.raises(ValueError, match="unavailable"):
        contract.evaluate(velocity, one_band, delta_band=0.01)


def test_duplicate_labels_and_corrupt_declared_total_are_rejected():
    velocity, labels, _ = _direct_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates()
    duplicate = [labels[0], labels[0], labels[2], labels[3]]
    with pytest.raises(ValueError, match="unique"):
        contract.evaluate(velocity, duplicate)

    physical, _ = _agent2_physical_family()
    corrupted = dict(physical)
    corrupted["velocity_physical_cylindrical_total"] = (
        np.asarray(physical["velocity_physical_cylindrical_total"]) + 1.0e-3
    )
    with pytest.raises(ValueError, match="does not equal"):
        contract.consume_physical_family(corrupted)


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.125)
    path = contract.save_json(tmp_path / "bounded_family_coefficients.json")
    loaded = KokunoBoundedFamilyCoefficientCoordinates.load_json(path)
    assert loaded.sha256 == contract.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["repository_coefficient_unit_mapping_available"] is True
    assert truth["low_dimensional_bounded_parameterization_available"] is True
    assert truth["source_coefficient_unit_mapping_available"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["pde_validated"] is False

    payload.pop("sha256")
    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoBoundedFamilyCoefficientCoordinates.from_payload(payload)
