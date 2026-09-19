from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_candidate_velocity import (
    KokunoPublicCandidateOscillatoryVelocity,
    velocity_osc,
)


def test_public_provider_is_nonzero_time_varying_and_aggregates_by_sign_beta() -> None:
    field = KokunoPublicCandidateOscillatoryVelocity()
    out = field.evaluate(
        x=np.asarray((0.72, 0.83, 0.64)),
        y=np.asarray((0.11, -0.16, 0.19)),
        z=np.asarray((0.08, -0.14, 0.21)),
        t=np.asarray((0.31, 0.47, 0.69)),
    )
    total = np.asarray(out["velocity_cartesian_total"])
    by_beta = np.asarray(out["velocity_cartesian_by_beta"])
    by_sign = np.asarray(out["velocity_cartesian_by_beta_sign"])
    assert total.shape == (3, 3)
    assert by_beta.shape[:2] == (3, len(out["beta_labels"]))
    assert by_sign.shape[:3] == (3, len(out["beta_labels"]), 2)
    assert np.all(np.isfinite(total))
    assert np.linalg.norm(total) > 0.0
    np.testing.assert_allclose(np.sum(by_sign, axis=-2), by_beta, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(np.sum(by_beta, axis=-2), total, rtol=0.0, atol=2e-12)
    assert out["public_xyz_t_velocity_correction_materialized"] is True
    assert out["pde_validated"] is False
    assert out["paper_exact"] is False

    early = field.velocity_osc(0.78, 0.12, 0.09, 0.30)
    late = field.velocity_osc(0.78, 0.12, 0.09, 0.68)
    assert early.shape == (3,)
    assert late.shape == (3,)
    assert np.linalg.norm(early - late) > 1e-10
    np.testing.assert_array_equal(velocity_osc(0.78, 0.12, 0.09, 0.30), early)


def test_axis_and_outside_annulus_are_exact_zero_without_cylindrical_evaluation() -> None:
    field = KokunoPublicCandidateOscillatoryVelocity()
    inner = field.radial_inner
    outer = field.radial_outer
    x = np.asarray((0.0, 0.5 * inner, outer + 0.05, 0.76))
    y = np.zeros_like(x)
    z = np.asarray((0.0, 0.2, -0.1, 0.3))
    t = np.full_like(x, 0.5)
    out = field.evaluate(x, y, z, t)
    velocity = np.asarray(out["velocity_cartesian_total"])
    np.testing.assert_array_equal(velocity[:3], np.zeros((3, 3)))
    assert np.linalg.norm(velocity[3]) > 0.0
    np.testing.assert_array_equal(out["support_mask"], np.asarray((False, False, False, True)))


def test_analytic_candidate_mode_coefficient_radial_derivative_matches_fd4() -> None:
    field = KokunoPublicCandidateOscillatoryVelocity()
    R0 = 0.78
    theta = 0.31
    Z = 0.12
    time = 0.43
    h = 2.0e-5

    def C_at(R: float) -> np.ndarray:
        return np.asarray(
            field._inside_family(
                np.asarray((R,)), np.asarray((theta,)), np.asarray((Z,)), np.asarray((time,))
            )["candidate_mode_C_plus_prototype"]
        )[0]

    center = field._inside_family(
        np.asarray((R0,)), np.asarray((theta,)), np.asarray((Z,)), np.asarray((time,))
    )
    analytic = np.asarray(center["candidate_mode_D_r_C_plus_prototype"])[0]
    fd4 = (-C_at(R0 + 2 * h) + 8 * C_at(R0 + h) - 8 * C_at(R0 - h) + C_at(R0 - 2 * h)) / (12 * h)
    np.testing.assert_allclose(analytic, fd4, rtol=3e-7, atol=3e-8)
    np.testing.assert_array_equal(
        np.asarray(center["candidate_mode_D_z_C_plus_prototype"]),
        np.zeros_like(np.asarray(center["candidate_mode_D_z_C_plus_prototype"])),
    )


def test_serialization_freezes_numerical_realization_and_truth_boundary(tmp_path) -> None:
    field = KokunoPublicCandidateOscillatoryVelocity()
    payload = field.to_payload()
    truth = payload["truth_boundary"]
    assert truth["public_xyz_t_velocity_correction_materialized"] is True
    assert truth["candidate_autonomous_background_bound"] is True
    assert truth["candidate_autonomous_signed_mode_bound"] is True
    assert truth["candidate_autonomous_partition_bound"] is True
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert len(payload["frozen_realization"]["beta_labels"]) == len(field.beta_labels)
    assert len(payload["frozen_realization"]["h_sigma_by_beta_sign"]) == len(field.beta_labels)

    path = field.save_json(tmp_path / "kokuno_public_candidate.json")
    loaded = KokunoPublicCandidateOscillatoryVelocity.load_json(path)
    assert loaded.sha256 == field.sha256
    np.testing.assert_array_equal(
        loaded.velocity_osc(0.81, -0.13, 0.17, 0.52),
        field.velocity_osc(0.81, -0.13, 0.17, 0.52),
    )


def test_bounds_and_registered_time_fail_closed() -> None:
    with pytest.raises(ValueError, match="strictly away"):
        KokunoPublicCandidateOscillatoryVelocity(radial_center=0.4, radial_halfwidth=0.5)
    with pytest.raises(ValueError, match="strict positive amplitude gap"):
        KokunoPublicCandidateOscillatoryVelocity(normal_target_ratio=0.10, cross_target_ratio=0.10)
    field = KokunoPublicCandidateOscillatoryVelocity()
    with pytest.raises(ValueError, match="registered candidate interval"):
        field.velocity_osc(0.8, 0.1, 0.0, field.time_max + 0.01)
