from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_axial_correction import (
    TARGET_MODE,
    lift_eq45_axial_eta4,
)
from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate


def _coefficient_map(candidate: Eq45VelocityCandidate, channel: str) -> dict[tuple[int, int], float]:
    basis = candidate.profile_basis
    values = getattr(basis, channel)
    return dict(zip(basis.mode_indices, values, strict=True))


def test_zero_lift_preserves_existing_modes_and_public_velocity_exactly() -> None:
    base = Eq45VelocityCandidate.seed()
    lifted = lift_eq45_axial_eta4(base)

    assert lifted.profile_basis.radial_degree == base.profile_basis.radial_degree
    assert lifted.profile_basis.eta_degree == 4

    base_phi = _coefficient_map(base, "phi_coefficients")
    base_swirl = _coefficient_map(base, "swirl_coefficients")
    lifted_phi = _coefficient_map(lifted, "phi_coefficients")
    lifted_swirl = _coefficient_map(lifted, "swirl_coefficients")

    for mode, value in base_phi.items():
        assert lifted_phi[mode] == value
    for mode, value in base_swirl.items():
        assert lifted_swirl[mode] == value
    for mode in set(lifted_phi) - set(base_phi):
        assert lifted_phi[mode] == 0.0
        assert lifted_swirl[mode] == 0.0

    points = np.array(
        [
            [0.10, 0.00, 0.10],
            [0.20, -0.10, -0.15],
            [-0.18, 0.22, 0.05],
            [0.00, 0.00, 0.20],
        ],
        dtype=float,
    )
    times = np.array([0.25, 0.40, 0.60, 0.75], dtype=float)
    np.testing.assert_array_equal(lifted.at_points(points, times), base.at_points(points, times))


def test_nonzero_lift_sets_only_screened_eta4_modes_and_round_trips(tmp_path) -> None:
    base = Eq45VelocityCandidate.seed()
    lifted = lift_eq45_axial_eta4(base, phi_eta4=0.20, swirl_eta4=-0.15)

    base_modes = set(base.profile_basis.mode_indices)
    phi = _coefficient_map(lifted, "phi_coefficients")
    swirl = _coefficient_map(lifted, "swirl_coefficients")
    new_modes = set(lifted.profile_basis.mode_indices) - base_modes

    assert phi[TARGET_MODE] == pytest.approx(0.20)
    assert swirl[TARGET_MODE] == pytest.approx(-0.15)
    for mode in new_modes - {TARGET_MODE}:
        assert phi[mode] == 0.0
        assert swirl[mode] == 0.0

    path = tmp_path / "lifted.json"
    lifted.save_json(path)
    loaded = Eq45VelocityCandidate.load_json(path)
    assert loaded.to_dict() == lifted.to_dict()
    assert loaded.sha256 == lifted.sha256


def test_phi_eta4_is_a_real_velocity_direction_without_claim_promotion() -> None:
    base = Eq45VelocityCandidate.seed()
    lifted = lift_eq45_axial_eta4(base, phi_eta4=0.10)

    points = np.array(
        [
            [0.30, 0.10, 0.20],
            [0.30, 0.10, -0.20],
            [-0.25, 0.20, 0.30],
            [-0.25, 0.20, -0.30],
        ],
        dtype=float,
    )
    base_velocity = base.at_points(points, 0.5)
    lifted_velocity = lifted.at_points(points, 0.5)
    assert float(np.linalg.norm(lifted_velocity - base_velocity)) > 1e-8

    truth = lifted.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_lift_rejects_invalid_or_existing_eta4_channels() -> None:
    base = Eq45VelocityCandidate.seed()
    limit = base.profile_basis.coefficient_limit

    with pytest.raises(ValueError, match="finite"):
        lift_eq45_axial_eta4(base, phi_eta4=float("nan"))
    with pytest.raises(ValueError, match="coefficient_limit"):
        lift_eq45_axial_eta4(base, swirl_eta4=limit + 1.0)

    already_lifted = lift_eq45_axial_eta4(base)
    with pytest.raises(ValueError, match="already contains"):
        lift_eq45_axial_eta4(already_lifted, phi_eta4=0.1)
    with pytest.raises(TypeError, match="Eq45VelocityCandidate"):
        lift_eq45_axial_eta4(object())  # type: ignore[arg-type]
