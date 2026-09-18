import numpy as np
import pytest

from openai_ns_reconstruction.constrained_st051b_piola_capacity import (
    BETA,
    CRITERIA,
    apply_piola_grid,
    decision,
    sampled_kinetic_energy,
    warp_z_and_jacobian,
)


def test_piola_map_identity_endpoints_and_orientation():
    z = np.linspace(-2.0, 2.0, 4001)
    mapped0, jac0 = warp_z_and_jacobian(z, 0.0)
    assert np.array_equal(mapped0, z)
    assert np.array_equal(jac0, np.ones_like(z))

    mapped, jac = warp_z_and_jacobian(z, BETA)
    assert mapped[0] == pytest.approx(-2.0)
    assert mapped[-1] == pytest.approx(2.0)
    assert np.min(jac) > 0.0
    assert np.all(np.diff(mapped) > 0.0)
    with pytest.raises(ValueError):
        warp_z_and_jacobian(z, -0.01)


def test_zero_beta_grid_replays_and_energy_scale_is_one():
    axis = np.linspace(-2.0, 2.0, 9)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    one = np.stack((0.1 * x, -0.2 * y, 0.05 * z), axis=-1)
    fields = np.stack((one, 1.1 * one, 0.9 * one), axis=0)
    warped, scale, structure = apply_piola_grid(axis, fields, beta=0.0)
    assert scale == pytest.approx(1.0, abs=1e-14)
    assert np.max(np.abs(warped - fields)) < 1e-14
    assert structure["map_monotone"]
    assert sampled_kinetic_energy(axis, warped[0]) == pytest.approx(
        sampled_kinetic_energy(axis, fields[0]), rel=1e-14, abs=1e-14
    )


def test_decision_requires_axial_gain_and_independent_response():
    fine = []
    for time in (0.25, 0.50, 0.75):
        fine.append(
            {
                "time": time,
                "relative": {
                    "mean_axial_span": 0.011,
                    "mean_turns": 0.003,
                    "mean_radial_span": 0.02,
                },
            }
        )
    response = {
        "normalized_rank": 2,
        "normalized_condition_number": 1.2,
        "normalized_column_cosine": 0.1,
    }
    structure = {"map_monotone": True, "minimum_coordinate_jacobian": 0.9}
    assert decision(fine, response, structure)["useful_axial_capacity_transfer"]

    fine[1]["relative"]["mean_axial_span"] = CRITERIA["minimum_axial_span_gain"] - 1e-6
    assert not decision(fine, response, structure)["useful_axial_capacity_transfer"]
