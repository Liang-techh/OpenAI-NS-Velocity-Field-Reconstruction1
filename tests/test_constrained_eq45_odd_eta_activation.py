from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_odd_eta_activation import (
    ODD_ETA_MODE,
    activate_eq45_odd_eta01,
)
from openai_ns_reconstruction.constrained_eq45_profile_basis import Eq45CompactProfileBasis
from openai_ns_reconstruction.eq45_delivery import Eq45DeliveryField


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _seed() -> Eq45VelocityCandidate:
    return Eq45VelocityCandidate.load_json(SEED)


def test_zero_odd_eta_activation_is_exact_function_preserving_reparameterization():
    candidate = _seed()
    updated = activate_eq45_odd_eta01(candidate, phi01=0.0, swirl01=0.0)

    assert updated.sha256 == candidate.sha256
    assert updated.to_dict() == candidate.to_dict()

    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.3, -0.1, 0.2],
            [-0.5, 0.2, -0.3],
            [0.7, 0.0, 0.4],
        ]
    )
    np.testing.assert_allclose(
        updated.at_points(points, 0.5),
        candidate.at_points(points, 0.5),
        rtol=0.0,
        atol=0.0,
    )


def test_phi01_breaks_only_midplane_poloidal_reflection_channel_and_roundtrips(tmp_path):
    candidate = _seed()
    basis = candidate.profile_basis
    index = basis.mode_indices.index(ODD_ETA_MODE)
    assert basis.phi_coefficients[index] == 0.0

    updated = activate_eq45_odd_eta01(candidate, phi01=0.2)
    assert updated.profile_basis.mode_indices == basis.mode_indices
    assert updated.profile_basis.radial_degree == basis.radial_degree
    assert updated.profile_basis.eta_degree == basis.eta_degree
    assert updated.profile_basis.x_cut == basis.x_cut
    assert updated.profile_basis.eta_cut == basis.eta_cut
    assert updated.profile_basis.cutoff_power == basis.cutoff_power
    assert updated.profile_basis.coefficient_limit == basis.coefficient_limit
    assert updated.profile_basis.swirl_coefficients == basis.swirl_coefficients

    changed = [
        i
        for i, (before, after) in enumerate(
            zip(basis.phi_coefficients, updated.profile_basis.phi_coefficients)
        )
        if before != after
    ]
    assert changed == [index]
    assert updated.profile_basis.phi_coefficients[index] == pytest.approx(0.2)

    midplane = np.array([[0.45, 0.0, 0.0], [0.7, 0.0, 0.0]])
    baseline = candidate.at_points(midplane, 0.5)
    activated = updated.at_points(midplane, 0.5)
    np.testing.assert_allclose(baseline[:, 0], 0.0, rtol=0.0, atol=1e-14)
    assert np.max(np.abs(activated[:, 0])) > 1e-5

    saved = tmp_path / "eq45_odd_phi.json"
    updated.save_json(saved)
    reloaded = Eq45VelocityCandidate.load_json(saved)
    assert reloaded.sha256 == updated.sha256

    axes = np.array([-0.3, 0.0, 0.3])
    times = np.array([0.25, 0.5, 0.75])
    direct_grid = updated.grid(axes, axes, axes, times)
    delivery_grid = Eq45DeliveryField(candidate=reloaded).grid(axes, axes, axes, times)
    np.testing.assert_allclose(delivery_grid, direct_grid, rtol=0.0, atol=0.0)

    truth = reloaded.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_swirl01_independently_adds_z_odd_swirl_without_changing_phi():
    candidate = _seed()
    basis = candidate.profile_basis
    index = basis.mode_indices.index(ODD_ETA_MODE)
    assert basis.swirl_coefficients[index] == 0.0

    updated = activate_eq45_odd_eta01(candidate, swirl01=-0.3)
    assert updated.profile_basis.phi_coefficients == basis.phi_coefficients
    changed = [
        i
        for i, (before, after) in enumerate(
            zip(basis.swirl_coefficients, updated.profile_basis.swirl_coefficients)
        )
        if before != after
    ]
    assert changed == [index]

    mirrored = np.array([[0.5, 0.0, 0.25], [0.5, 0.0, -0.25]])
    baseline = candidate.at_points(mirrored, 0.5)
    activated = updated.at_points(mirrored, 0.5)
    np.testing.assert_allclose(baseline[0, 1], baseline[1, 1], rtol=0.0, atol=1e-14)
    assert abs(float(activated[0, 1] - activated[1, 1])) > 1e-5


def test_odd_eta_activation_fails_closed_on_invalid_or_unavailable_controls():
    candidate = _seed()
    limit = candidate.profile_basis.coefficient_limit

    with pytest.raises(ValueError, match="provide phi01 and/or swirl01"):
        activate_eq45_odd_eta01(candidate)
    with pytest.raises(ValueError, match="phi01 must be finite"):
        activate_eq45_odd_eta01(candidate, phi01=np.nan)
    with pytest.raises(ValueError, match="exceeds coefficient_limit"):
        activate_eq45_odd_eta01(candidate, swirl01=limit + 0.01)
    with pytest.raises(ValueError, match="finite scalar or None"):
        activate_eq45_odd_eta01(candidate, phi01=True)
    with pytest.raises(TypeError, match="Eq45VelocityCandidate"):
        activate_eq45_odd_eta01(object(), phi01=0.1)

    no_odd_basis = Eq45CompactProfileBasis(
        radial_degree=0,
        eta_degree=0,
        phi_coefficients=(1.0,),
        swirl_coefficients=(0.5,),
    )
    no_odd_candidate = Eq45VelocityCandidate(profile_basis=no_odd_basis)
    with pytest.raises(ValueError, match=r"does not contain mode \(0, 1\)"):
        activate_eq45_odd_eta01(no_odd_candidate, phi01=0.1)
