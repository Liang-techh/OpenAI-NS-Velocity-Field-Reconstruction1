from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_partition_family import (
    KokunoSourcePartitionFamilyContract,
)


def _manufactured_family(R, Z):
    R, Z = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(Z, dtype=float))
    alpha = 0.31 + 0.27 * R - 0.19 * Z + 0.04 * R * Z
    D_r_alpha = 0.27 + 0.04 * Z
    # Manufactured source-normalized axial derivative for epsilon=.2.
    D_z_alpha = 0.2 * (-0.19 + 0.04 * R)
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = D_r_alpha[..., None] * tangent
    D_z_eta = D_z_alpha[..., None] * tangent
    return eta, D_r_eta, D_z_eta


def _fd4(fn, x, h):
    return (fn(x - 2 * h) - 8 * fn(x - h) + 8 * fn(x + h) - fn(x + 2 * h)) / (12 * h)


def test_squared_partition_and_derived_closures_hold_on_batch():
    contract = KokunoSourcePartitionFamilyContract()
    R = np.array([0.72, 0.91, 1.18, 1.37])
    Z = np.array([-0.16, 0.02, 0.13, -0.07])
    eta, D_r_eta, D_z_eta = _manufactured_family(R, Z)
    out = contract.validate_family(
        eta,
        D_r_eta,
        D_z_eta,
        beta_labels=[(12, (3, 4)), (12, (4, 4))],
    )
    np.testing.assert_allclose(out["squared_sum"], 1.0, rtol=0.0, atol=3e-16)
    np.testing.assert_allclose(out["radial_closure"], 0.0, rtol=0.0, atol=3e-16)
    np.testing.assert_allclose(out["axial_closure"], 0.0, rtol=0.0, atol=3e-16)
    assert out["beta_labels"] == ((12, (3, 4)), (12, (4, 4)))


def test_derivative_closure_matches_independent_fd_of_partition_sum():
    contract = KokunoSourcePartitionFamilyContract()
    R, Z = 0.93, -0.11
    eta, D_r_eta, D_z_eta = _manufactured_family(R, Z)
    out = contract.validate_family(
        eta,
        D_r_eta,
        D_z_eta,
        beta_labels=[(17, "a0"), (17, "a1")],
    )

    h = 2.0e-4
    radial_fd = _fd4(
        lambda RR: np.sum(_manufactured_family(RR, Z)[0] ** 2), R, h
    )
    axial_fd = 0.2 * _fd4(
        lambda ZZ: np.sum(_manufactured_family(R, ZZ)[0] ** 2), Z, h
    )
    np.testing.assert_allclose(2.0 * out["radial_closure"], radial_fd, atol=2e-12, rtol=0.0)
    np.testing.assert_allclose(2.0 * out["axial_closure"], axial_fd, atol=2e-12, rtol=0.0)


def test_bad_partition_or_derivative_fails_closed():
    contract = KokunoSourcePartitionFamilyContract()
    eta, D_r_eta, D_z_eta = _manufactured_family(0.8, 0.1)
    labels = [(9, "left"), (9, "right")]

    bad_eta = eta.copy()
    bad_eta[0] += 1.0e-3
    with pytest.raises(ValueError, match="squared partition"):
        contract.validate_family(bad_eta, D_r_eta, D_z_eta, labels)

    bad_D_r = D_r_eta.copy()
    bad_D_r[0] += 2.0e-5
    with pytest.raises(ValueError, match="D_r_eta"):
        contract.validate_family(eta, bad_D_r, D_z_eta, labels)

    bad_D_z = D_z_eta.copy()
    bad_D_z[1] -= 2.0e-5
    with pytest.raises(ValueError, match="D_z_eta"):
        contract.validate_family(eta, D_r_eta, bad_D_z, labels)


def test_partition_axis_counts_beta_once_not_plus_minus_gamma():
    contract = KokunoSourcePartitionFamilyContract()
    eta, D_r_eta, D_z_eta = _manufactured_family(1.02, -0.04)

    with pytest.raises(ValueError, match="do not duplicate gamma"):
        contract.validate_family(
            eta,
            D_r_eta,
            D_z_eta,
            beta_labels=[(14, "a0", +1), (14, "a1", -1)],
        )

    with pytest.raises(ValueError, match="unique"):
        contract.validate_family(
            eta,
            D_r_eta,
            D_z_eta,
            beta_labels=[(14, "a0"), (14, "a0")],
        )

    duplicated_eta = np.repeat(eta / np.sqrt(2.0), 2)
    duplicated_D_r = np.repeat(D_r_eta / np.sqrt(2.0), 2)
    duplicated_D_z = np.repeat(D_z_eta / np.sqrt(2.0), 2)
    # Even if numerical weights are rescaled to hide double counting, gamma
    # labels are rejected because the source squared sum is indexed by beta.
    with pytest.raises(ValueError, match="do not duplicate gamma"):
        contract.validate_family(
            duplicated_eta,
            duplicated_D_r,
            duplicated_D_z,
            beta_labels=[
                (14, "a0", +1), (14, "a0", -1),
                (14, "a1", +1), (14, "a1", -1),
            ],
        )


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    contract = KokunoSourcePartitionFamilyContract(
        partition_atol=1.0e-11,
        derivative_rtol=3.0e-9,
    )
    path = contract.save_json(tmp_path / "source_partition_family.json")
    loaded = KokunoSourcePartitionFamilyContract.load_json(path)
    assert loaded.sha256 == contract.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["source_squared_partition_family_structure_executable"] is True
    assert truth["source_beta_vs_gamma_counting_guard_executable"] is True
    assert truth["derived_partition_D_r_closure_guard_executable"] is True
    assert truth["derived_partition_D_z_closure_guard_executable"] is True
    assert truth["concrete_source_partition_bumps_reconstructed"] is False
    assert truth["source_actual_partition_labels_instantiated"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourcePartitionFamilyContract.from_payload(payload)
