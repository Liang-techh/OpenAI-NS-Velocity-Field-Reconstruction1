import json
import numpy as np

from openai_ns_reconstruction.kokuno_public_carrier_resolved_velocity import (
    KokunoCarrierResolvedCandidateOscillatoryVelocity,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import (
    KokunoPublicZPullbackCandidateOscillatoryVelocity,
    velocity_osc,
)


def _phase_frame_at_public_z(field, public_z):
    z = np.asarray(public_z, dtype=float)
    if z.ndim != 1:
        raise ValueError("test public_z must be one-dimensional")
    labels = field.beta_labels
    n_beta = len(labels)
    geometry = field._static["geometry"]
    L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
    R = np.full(z.shape, 0.73)
    theta = np.full(z.shape, 0.31)
    R_beta = np.broadcast_to(R[:, None], (z.size, n_beta))
    v = np.stack((0.25 * L_s, 0.75 * L_s), axis=-1)
    return field.bridge.phase_family.phase_frame(
        R=R_beta,
        Z=field._source_Z_by_beta(z),
        theta=theta,
        v=v,
        beta_labels=labels,
        L_s=L_s,
        **field._background(),
    )


def test_bandwise_public_z_pullback_makes_source_nphi_axial_component_the_public_phase_derivative():
    field = KokunoPublicZPullbackCandidateOscillatoryVelocity()
    z = np.asarray((0.17, -0.23), dtype=float)
    step = 2.0e-6

    center = _phase_frame_at_public_z(field, z)
    plus = _phase_frame_at_public_z(field, z + step)
    minus = _phase_frame_at_public_z(field, z - step)
    fd_public_z = (np.asarray(plus["phase"]) - np.asarray(minus["phase"])) / (2.0 * step)
    nphi_z = np.asarray(center["n_phi"])[..., 2]

    # This is the coordinate seam being repaired: with Z_beta=epsilon_beta*z,
    # ordinary public partial_z Phi equals the source-normalized D_z Phi=n_Phi,z.
    assert np.allclose(fd_public_z, nphi_z, rtol=2.0e-9, atol=2.0e-9)

    epsilon = np.asarray(field._static["geometry"]["epsilon_by_beta"], dtype=float)
    source_z = field._source_Z_by_beta(z)
    assert np.allclose(source_z, z[:, None] * epsilon[None, :], rtol=0.0, atol=0.0)

    # The old direct Z=z identification instead differentiates the source phase
    # as p_z/epsilon in public z.  Check the exact expected mismatch without
    # using any Agent-4 divergence implementation as an oracle.
    labels = field.beta_labels
    n_beta = len(labels)
    geometry = field._static["geometry"]
    L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
    R_beta = np.broadcast_to(np.full((z.size, 1), 0.73), (z.size, n_beta))
    theta = np.full(z.shape, 0.31)
    v = np.stack((0.25 * L_s, 0.75 * L_s), axis=-1)

    def direct_frame(zz):
        return field.bridge.phase_family.phase_frame(
            R=R_beta,
            Z=np.broadcast_to(np.asarray(zz)[:, None], (z.size, n_beta)),
            theta=theta,
            v=v,
            beta_labels=labels,
            L_s=L_s,
            **field._background(),
        )

    old_plus = direct_frame(z + step)
    old_minus = direct_frame(z - step)
    old_fd = (np.asarray(old_plus["phase"]) - np.asarray(old_minus["phase"])) / (2.0 * step)
    assert np.allclose(
        old_fd * epsilon[None, :, None], nphi_z, rtol=2.0e-9, atol=2.0e-9
    )
    assert float(np.max(np.abs(old_fd - nphi_z))) > 1.0e-7


def test_public_z_pullback_provider_remains_nontrivial_finite_compact_and_aggregated():
    field = KokunoPublicZPullbackCandidateOscillatoryVelocity()
    points = np.asarray(
        (
            (0.62, 0.13, 0.17, 0.37),
            (0.71, -0.28, -0.21, 0.49),
            (0.88, 0.19, 0.31, 0.61),
        ),
        dtype=float,
    )
    result = field.evaluate(points[:, 0], points[:, 1], points[:, 2], points[:, 3])
    vel = np.asarray(result["velocity_cartesian_total"])
    by_beta = np.asarray(result["velocity_cartesian_by_beta"])
    by_beta_sign = np.asarray(result["velocity_cartesian_by_beta_sign"])

    assert vel.shape == (3, 3)
    assert np.all(np.isfinite(vel))
    assert np.linalg.norm(vel) > 1.0e-8
    assert np.allclose(np.sum(by_beta, axis=-2), vel, rtol=2.0e-13, atol=2.0e-13)
    assert np.allclose(np.sum(by_beta_sign, axis=-2), by_beta, rtol=2.0e-13, atol=2.0e-13)
    assert result["public_source_D_z_pullback_applied"] is True
    assert "Z_beta=epsilon_beta*z" in result["coordinate_contract"]

    exterior = velocity_osc(
        np.asarray((0.0, 1.50, 0.72, 0.72)),
        np.asarray((0.0, 0.0, 0.0, 0.0)),
        np.asarray((0.0, 0.0, 2.05, -2.05)),
        np.asarray((0.41, 0.43, 0.47, 0.53)),
    )
    assert np.array_equal(exterior, np.zeros_like(exterior))


def test_public_z_pullback_payload_preserves_parent_background_and_truth_boundary(tmp_path):
    parent = KokunoCarrierResolvedCandidateOscillatoryVelocity()
    field = KokunoPublicZPullbackCandidateOscillatoryVelocity()
    parent_payload = parent.to_payload()
    payload = field.to_payload()
    receipt = payload["frozen_realization"]["public_z_pullback"]

    assert payload["schema"] == "kokuno-public-z-pullback-oscillatory-velocity-v1"
    assert payload["frozen_realization"]["background"] == parent_payload["frozen_realization"]["background"]
    assert receipt["coordinate_pullback"] == "Z_beta=epsilon_beta*z; partial_z=D_z(source)"
    assert receipt["source_formulas_changed"] is False
    assert receipt["autonomous_background_changed_from_parent"] is False
    assert receipt["vector_potential_complete_curl_path_changed"] is False
    assert receipt["public_coordinate_pullback_changed"] is True
    assert receipt["independent_agent4_reaudit_required"] is True

    truth = payload["truth_boundary"]
    assert truth["public_source_D_z_pullback_applied"] is True
    assert truth["public_coordinate_curl_consistency_targeted"] is True
    assert truth["independent_fd4_divergence_reaudit_required"] is True
    assert truth["independent_fd4_divergence_passed"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = field.save_json(tmp_path / "public-z-pullback.json")
    raw = json.loads(path.read_text())
    loaded = KokunoPublicZPullbackCandidateOscillatoryVelocity.load_json(path)
    assert raw["sha256"] == field.sha256
    assert loaded.sha256 == field.sha256
    assert loaded.to_payload() == field.to_payload()


def test_pullback_changes_only_public_coordinate_realization_not_frozen_parameters():
    parent = KokunoCarrierResolvedCandidateOscillatoryVelocity()
    field = KokunoPublicZPullbackCandidateOscillatoryVelocity()

    assert field.h == parent.h
    assert field.radial_center == parent.radial_center
    assert field.radial_halfwidth == parent.radial_halfwidth
    assert field.axial_center == parent.axial_center
    assert field.axial_halfwidth == parent.axial_halfwidth
    assert field.normal_target_ratio == parent.normal_target_ratio
    assert field.cross_target_ratio == parent.cross_target_ratio
    assert field.time_modulation == parent.time_modulation
    assert field.mode_imaginary_ratio == parent.mode_imaginary_ratio
    assert field._background().keys() == parent._background().keys()
    for key in field._background():
        assert np.array_equal(field._background()[key], parent._background()[key])
