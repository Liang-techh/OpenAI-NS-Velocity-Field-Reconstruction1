from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import (
    Eq45VelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_poloidal_streamfunction import (
    eq45_streamfunction_poloidal_jet,
)
from openai_ns_reconstruction.constrained_eq45_profile_basis import (
    Eq45CompactProfileBasis,
)
from openai_ns_reconstruction.constrained_eq45_velocity import Eq45VelocityBackbone


ROOT = Path(__file__).resolve().parents[1]
SEED_ARTIFACT = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"


def test_seed_is_nontrivial_vectorized_axis_regular_and_time_bounded():
    candidate = Eq45VelocityCandidate.seed()
    points = np.array(
        [
            [0.0, 0.0, -0.9],
            [0.35, 0.0, -0.35],
            [0.25, -0.40, 0.0],
            [0.0, 0.0, 0.45],
            [-0.55, 0.30, 0.8],
        ]
    )
    times = np.array([0.25, 0.375, 0.50, 0.625, 0.75])
    velocity = candidate.velocity(points, times)

    assert velocity.shape == points.shape
    assert np.all(np.isfinite(velocity))
    assert np.max(np.linalg.norm(velocity, axis=-1)) > 1e-6
    assert np.array_equal(velocity[[0, 3], :2], np.zeros((2, 2)))
    assert np.array_equal(candidate.at_points(points, times), velocity)

    with pytest.raises(ValueError, match="declared interval"):
        candidate.velocity(points[:1], 0.249)
    with pytest.raises(ValueError, match="declared interval"):
        candidate.velocity(points[:1], 0.751)


def test_composition_matches_profile_adapter_and_cartesian_backbone():
    candidate = Eq45VelocityCandidate.seed()
    points = np.array(
        [
            [0.15, 0.20, -0.60],
            [0.75, -0.10, -0.15],
            [-0.55, 0.45, 0.30],
            [0.20, -0.85, 0.70],
        ]
    )
    times = np.array([0.30, 0.42, 0.58, 0.72])

    backbone = Eq45VelocityBackbone(candidate.profile_values, h=candidate.h)
    coordinates = backbone.coordinates(points, times)
    basis_jets = candidate.profile_basis.evaluate(coordinates.X, coordinates.eta)
    poloidal = eq45_streamfunction_poloidal_jet(
        coordinates.X,
        coordinates.eta,
        candidate.h,
        basis_jets.phi,
        basis_jets.phi_x,
        basis_jets.phi_eta,
        basis_jets.phi_xx,
        basis_jets.phi_xeta,
    )
    expected_profiles = np.stack(
        (poloidal.v0, basis_jets.swirl, poloidal.U), axis=-1
    )

    assert np.array_equal(
        candidate.profile_values(coordinates.X, coordinates.eta), expected_profiles
    )
    assert np.array_equal(candidate.velocity(points, times), backbone.velocity(points, times))


def test_seed_artifact_round_trip_and_xyz_interface(tmp_path):
    candidate = Eq45VelocityCandidate.seed()
    checked = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)
    assert checked.to_dict() == candidate.to_dict()
    assert checked.sha256 == candidate.sha256

    target = tmp_path / "eq45_candidate.json"
    candidate.save_json(target)
    restored = Eq45VelocityCandidate.load_json(target)
    assert restored.to_dict() == candidate.to_dict()
    assert restored.sha256 == candidate.sha256

    x = np.array([0.0, 0.2, -0.4])
    y = np.array([0.0, -0.3, 0.5])
    z = np.array([-0.4, 0.0, 0.6])
    t = np.array([0.25, 0.5, 0.75])
    via_xyz = restored.velocity_xyz(x, y, z, t)
    via_points = restored.velocity(np.stack((x, y, z), axis=-1), t)
    assert np.array_equal(via_xyz, via_points)

    truth = restored.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_fail_closed_zero_profile_and_claim_tampering():
    seed_basis = Eq45CompactProfileBasis.seed()
    zeros = (0.0,) * seed_basis.mode_count
    zero_basis = Eq45CompactProfileBasis(
        radial_degree=seed_basis.radial_degree,
        eta_degree=seed_basis.eta_degree,
        phi_coefficients=zeros,
        swirl_coefficients=zeros,
        x_cut=seed_basis.x_cut,
        eta_cut=seed_basis.eta_cut,
        cutoff_power=seed_basis.cutoff_power,
        coefficient_limit=seed_basis.coefficient_limit,
    )
    with pytest.raises(ValueError, match="nontrivial"):
        Eq45VelocityCandidate(profile_basis=zero_basis)

    payload = Eq45VelocityCandidate.seed().to_dict()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45VelocityCandidate.from_dict(payload)

    payload = Eq45VelocityCandidate.seed().to_dict()
    payload["classification"]["openai_hidden_profiles"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification"):
        Eq45VelocityCandidate.from_dict(payload)

    payload = Eq45VelocityCandidate.seed().to_dict()
    payload["schema"] = "paper_exact_velocity"
    with pytest.raises(ValueError, match="schema"):
        Eq45VelocityCandidate.from_dict(payload)
