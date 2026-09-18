import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)
from openai_ns_reconstruction.kokuno_source_localized_real_pair_family import (
    KokunoSourceLocalizedRealPairFamily,
)


def _inputs():
    y = np.array([5.31, 6.22])
    q = 2.0 ** (-y)
    D_r_q = q * np.array([1.0e-4, -8.0e-5])
    D_z_q = q * np.array([-7.0e-5, 9.0e-5])
    slow = np.array(
        [
            [1.10e-4, -0.70e-4, 0.40e-4],
            [0.80e-4, 1.30e-4, -1.00e-4],
        ]
    )
    D_r_slow = np.array(
        [
            [1.00e-8, -2.00e-8, 1.50e-8],
            [-1.50e-8, 0.70e-8, 1.20e-8],
        ]
    )
    D_z_slow = np.array(
        [
            [-0.80e-8, 0.40e-8, -0.60e-8],
            [1.10e-8, -1.30e-8, 0.90e-8],
        ]
    )
    return q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow


def _shifted(realization, step, direction):
    q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow = _inputs()
    if direction == "r":
        q = q + step * D_r_q
        slow = slow + step * D_r_slow
    elif direction == "z":
        q = q + step * D_z_q
        slow = slow + step * D_z_slow
    else:
        raise ValueError(direction)
    return realization.evaluate(q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow)


def test_source_compatible_partition_closes_and_derivatives_match_fd4():
    realization = KokunoSourceCompatiblePartitionRealization()
    base = realization.evaluate(*_inputs())

    assert np.max(np.abs(base["partition_error"])) < 1.0e-12
    assert np.max(np.abs(base["radial_closure"])) < 1.0e-12
    assert np.max(np.abs(base["axial_closure"])) < 1.0e-12
    assert base["source_actual_partition_recovered"] is False
    assert base["autonomous_realization"] is True

    h = 1.0e-4
    for direction, analytic_key in (("r", "D_r_eta"), ("z", "D_z_eta")):
        fm2 = _shifted(realization, -2.0 * h, direction)
        fm1 = _shifted(realization, -h, direction)
        fp1 = _shifted(realization, h, direction)
        fp2 = _shifted(realization, 2.0 * h, direction)
        assert fm2["beta_labels"] == fm1["beta_labels"] == base["beta_labels"]
        assert base["beta_labels"] == fp1["beta_labels"] == fp2["beta_labels"]
        fd4 = (fm2["eta"] - 8.0 * fm1["eta"] + 8.0 * fp1["eta"] - fp2["eta"]) / (12.0 * h)
        np.testing.assert_allclose(fd4, base[analytic_key], rtol=2.0e-7, atol=3.0e-10)


def test_source_mesh_overlap_and_agent2_family_integration():
    realization = KokunoSourceCompatiblePartitionRealization()
    data = realization.evaluate(*_inputs())

    for ell in data["ell_labels"]:
        assert data["mesh_by_ell"][ell] == pytest.approx(float(ell) ** -6, rel=0.0, abs=0.0)
    for label in data["beta_labels"]:
        assert len(label) == 2
        ell, a = label
        assert isinstance(ell, int)
        assert len(a) == 3

    for sample in range(data["eta"].shape[0]):
        active = np.flatnonzero(np.abs(data["eta"][sample]) > 1.0e-14)
        active_ells = sorted({data["beta_labels"][i][0] for i in active})
        assert active_ells[-1] - active_ells[0] <= 1
        assert active_ells[-1] - active_ells[0] <= 2  # source overlap bound

    # The new realization is directly consumable by the #380 supplied-data
    # localized real-pair family.  The mode data below are manufactured only
    # to test the software seam; they are not source pulse/background data.
    family_shape = data["eta"].shape
    n_phi = np.zeros(family_shape + (3,))
    n_phi[..., 0] = 1.0
    t_plus = np.zeros(family_shape + (3,), dtype=np.complex128)
    t_plus[..., 1] = 1.0
    zeros = np.zeros_like(t_plus)
    phase = np.zeros(family_shape)
    Q = np.full(family_shape, 0.5)
    R = np.array([1.1, 1.2])
    theta = np.array([0.2, -0.4])

    composed = KokunoSourceLocalizedRealPairFamily().physical_family(
        Q,
        R,
        theta,
        phase,
        n_phi,
        t_plus,
        zeros,
        zeros,
        data["eta"],
        data["D_r_eta"],
        data["D_z_eta"],
        data["beta_labels"],
    )
    assert composed["velocity_physical_cartesian_total"].shape == (2, 3)
    assert np.all(np.isfinite(composed["velocity_physical_cartesian_total"]))
    np.testing.assert_allclose(
        composed["velocity_cylindrical_by_beta"],
        2.0 * composed["plus_velocity_by_beta"].real,
        rtol=0.0,
        atol=3.0e-12,
    )


def test_fail_closed_partition_guard_and_truth_boundary_round_trip(tmp_path):
    realization = KokunoSourceCompatiblePartitionRealization()
    q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow = _inputs()

    with pytest.raises(ValueError, match="ell_min"):
        realization.evaluate(
            np.array([2.0 ** -4.8]),
            np.array([0.0]),
            np.array([0.0]),
            np.zeros((1, 3)),
            np.zeros((1, 3)),
            np.zeros((1, 3)),
        )
    with pytest.raises(ValueError, match=r"q.shape\+\(3,\)"):
        realization.evaluate(q, D_r_q, D_z_q, slow[:, :2], D_r_slow[:, :2], D_z_slow[:, :2])

    data = realization.evaluate(q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow)
    corrupted = data["eta"].copy()
    corrupted[..., 0] += 1.0e-4
    with pytest.raises(ValueError, match="squared partition"):
        realization.partition_contract.validate_family(
            corrupted,
            data["D_r_eta"],
            data["D_z_eta"],
            data["beta_labels"],
        )

    target = tmp_path / "partition.json"
    realization.save_json(target)
    loaded = KokunoSourceCompatiblePartitionRealization.load_json(target)
    assert loaded == realization
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["concrete_source_partition_bumps_reconstructed"] is False
    assert payload["truth_boundary"]["source_actual_partition_labels_instantiated"] is False
    assert payload["truth_boundary"]["public_xyz_t_velocity_correction_materialized"] is False

    payload["truth_boundary"]["concrete_source_partition_bumps_reconstructed"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceCompatiblePartitionRealization.from_payload(payload)
