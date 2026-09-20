from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_velocity import (
    KokunoPA10CartesianCenterVelocity,
)
from openai_ns_reconstruction.kokuno_pa10_source_c_normalized_physical_center import (
    KokunoPA10SourceCNormalizedPhysicalCenter,
)


def _interior_source_points(field: KokunoPA10CartesianCenterVelocity):
    _, xmax = field.source_X_interval
    X = np.asarray([0.0, 0.11 * xmax, 0.27 * xmax, 0.53 * xmax])
    eta = np.asarray([0.0, 0.18, -0.31, 0.42])
    t = np.asarray([0.27, 0.41, 0.59, 0.73])
    theta = np.asarray([0.0, 0.35, 1.2, 2.1])
    return X, eta, t, theta


def test_source_coordinate_roundtrip_and_public_cartesian_formula():
    field = KokunoPA10CartesianCenterVelocity()
    X, eta, t, theta = _interior_source_points(field)
    physical = field.cartesian_from_similarity(X, eta, t, theta)
    coords = field.similarity_coordinates(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    np.testing.assert_allclose(coords["X"], X, rtol=2e-13, atol=1e-15)
    np.testing.assert_allclose(coords["eta"], eta, rtol=2e-13, atol=2e-15)
    np.testing.assert_allclose(coords["q"], physical["q"], rtol=2e-13, atol=2e-15)

    values = field.values(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    profiles = field.source_center.values(X, eta)
    q = physical["q"]
    expected_u = (
        profiles["v_0"] / (2.0 * q) * physical["x"]
        - np.power(q, -field.A - 0.5) * profiles["F_0"] * physical["y"]
    )
    expected_v = (
        profiles["v_0"] / (2.0 * q) * physical["y"]
        + np.power(q, -field.A - 0.5) * profiles["F_0"] * physical["x"]
    )
    expected_w = np.power(q, -field.A) * profiles["U_0"]
    np.testing.assert_allclose(values["u"], expected_u, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(values["v"], expected_v, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(values["w"], expected_w, rtol=2e-14, atol=0.0)

    velocity = field.velocity(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    assert velocity.shape == (4, 3)
    assert np.all(np.isfinite(velocity))
    assert np.any(np.abs(velocity) > 0.0)


def test_axis_regular_and_rotation_equivariant():
    field = KokunoPA10CartesianCenterVelocity()
    axis = field.cartesian_from_similarity(
        np.zeros(3),
        np.asarray([-0.35, 0.0, 0.35]),
        np.asarray([0.3, 0.5, 0.7]),
    )
    axis_velocity = field.velocity(axis["x"], axis["y"], axis["z"], axis["t"])
    assert np.all(axis_velocity[:, :2] == 0.0)
    assert np.all(np.isfinite(axis_velocity[:, 2]))

    _, xmax = field.source_X_interval
    X = 0.3 * xmax
    eta = 0.23
    t = 0.52
    base = field.cartesian_from_similarity(X, eta, t, 0.0)
    angle = 0.73
    rotated = field.cartesian_from_similarity(X, eta, t, angle)
    v0 = field.velocity(base["x"], base["y"], base["z"], base["t"])
    v1 = field.velocity(
        rotated["x"], rotated["y"], rotated["z"], rotated["t"]
    )
    c, s = math.cos(angle), math.sin(angle)
    expected = np.asarray(
        [c * v0[0] - s * v0[1], s * v0[0] + c * v0[1], v0[2]]
    )
    np.testing.assert_allclose(v1, expected, rtol=3e-13, atol=1e-12)


def test_field_identity_excludes_certificate_execution_settings():
    parent_a = KokunoPA10SourceCNormalizedPhysicalCenter(chunk_size=2**15)
    parent_b = KokunoPA10SourceCNormalizedPhysicalCenter(
        physical_profiles=parent_a.physical_profiles,
        chunk_size=2**14,
    )
    assert parent_a.configuration() != parent_b.configuration()
    field_a = KokunoPA10CartesianCenterVelocity(source_center=parent_a)
    field_b = KokunoPA10CartesianCenterVelocity(source_center=parent_b)
    assert field_a.field_configuration() == field_b.field_configuration()
    assert field_a.field_sha256 == field_b.field_sha256
    assert field_a.evidence_configuration() != field_b.evidence_configuration()

    X, eta, t, theta = _interior_source_points(field_a)
    physical = field_a.cartesian_from_similarity(X, eta, t, theta)
    va = field_a.velocity(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    vb = field_b.velocity(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    np.testing.assert_array_equal(va, vb)


def test_configuration_roundtrip_and_fail_closed_domain():
    field = KokunoPA10CartesianCenterVelocity()
    payload = field.field_configuration()
    replay = KokunoPA10CartesianCenterVelocity.from_configuration(
        copy.deepcopy(payload)
    )
    assert replay.field_configuration() == payload
    assert replay.field_sha256 == field.field_sha256

    _, xmax = field.source_X_interval
    physical = field.cartesian_from_similarity(0.4 * xmax, 0.2, 0.5, 0.4)
    np.testing.assert_array_equal(
        replay.velocity(physical["x"], physical["y"], physical["z"], physical["t"]),
        field.velocity(physical["x"], physical["y"], physical["z"], physical["t"]),
    )

    outside = field.cartesian_from_similarity(xmax, 0.0, 0.5, 0.0)
    with pytest.raises(ValueError, match="outside the source inner X interval"):
        field.velocity(
            1.01 * outside["x"],
            1.01 * outside["y"],
            outside["z"],
            outside["t"],
        )
    with pytest.raises(ValueError, match="registered interval"):
        field.velocity(0.0, 0.0, 0.0, 0.9)


def test_report_truth_boundary_stays_center_only():
    field = KokunoPA10CartesianCenterVelocity()
    report = field.report()
    checks = report["machine_checks"]
    assert checks["coordinate_roundtrip_max_abs_X"] < 1e-12
    assert checks["coordinate_roundtrip_max_abs_eta"] < 1e-12
    assert checks["q_equation_relative_residual_max"] < 1e-14
    assert checks["axis_transverse_velocity_exact_zero"] is True
    assert checks["velocity_nontrivial_on_probe"] is True
    assert checks["all_probe_values_finite"] is True

    truth = report["truth_boundary"]
    assert truth["inner_cartesian_spacetime_center_velocity_executable"] is True
    assert truth["field_semantic_identity_separate_from_evidence_receipt_identity"] is True
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["fixed_point_correction_materialized"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["unified_cartesian_velocity_export_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert report["scientific_gates"]["momentum_max_l2"] == 1e-3
    assert report["scientific_gates"]["divergence_max_l2"] == 1e-5
    assert report["scientific_gates"]["free_residual_defined_forcing_forbidden"] is True
