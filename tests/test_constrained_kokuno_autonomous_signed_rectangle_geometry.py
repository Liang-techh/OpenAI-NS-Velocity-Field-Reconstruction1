from __future__ import annotations

import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)
from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering
from openai_ns_reconstruction.kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)


def _partition_result() -> dict:
    realization = KokunoSourceCompatiblePartitionRealization(ell_min=5)
    return realization.evaluate(
        q=np.asarray(2.0**-5.5),
        D_r_q=np.asarray(0.0),
        D_z_q=np.asarray(0.0),
        slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
    )


def test_receipt_keeps_autonomous_witness_separate_from_source_data() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry()
    receipt = geometry.receipt()
    witness = geometry.witness

    assert receipt["schema"] == "kokuno-autonomous-signed-rectangle-geometry-v1"
    assert receipt["source"]["release_date"] == "2026-09-09"
    assert witness.injectivity_guard_passed
    assert witness.separation_guard_passed
    assert witness.d_c_squared_exact > 0
    assert receipt["witness"]["center_binding"].startswith("repository_autonomous")
    assert receipt["witness"]["r0_binding"].startswith("repository_autonomous")

    truth = receipt["truth_boundary"]
    assert truth["autonomous_source_compatible_signed_rectangle_witness_instantiated"] is True
    assert truth["source_rectangle_centers_recovered"] is False
    assert truth["source_rectangle_radius_r0_recovered"] is False
    assert truth["source_actual_rectangle_labels_instantiated"] is False
    assert truth["source_actual_partition_labels_instantiated"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_partition_labels_receive_disjoint_signed_centers_and_source_scales() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry(h=0.005)
    partition = _partition_result()
    bound = geometry.instantiate(partition)

    labels = tuple(partition["beta_labels"])
    assert bound["beta_labels"] == labels
    assert len(bound["gamma_labels"]) == 2 * len(labels)
    assert bound["center_index_by_beta_sign"].shape == (len(labels), 2)
    assert bound["center_float_by_beta_sign"].shape == (len(labels), 2, 2)
    assert np.all(bound["center_index_by_beta_sign"][:, 0] != bound["center_index_by_beta_sign"][:, 1])
    assert bound["sign_labels"] == ("sigma_plus", "sigma_minus")

    for j, beta in enumerate(labels):
        ell = int(beta[0])
        band = KokunoSourceBandCovering(ell=ell, h=geometry.h)
        assert bound["Q_by_beta"][j] == pytest.approx(2.0**-ell, rel=2e-15)
        assert bound["epsilon_by_beta"][j] == pytest.approx(band.epsilon, rel=2e-15)
        assert bound["covering_level_by_beta"][j] == band.covering_level
        assert bound["L_s_by_beta"][j] == pytest.approx(
            2.0 * geometry.witness.r0 / band.c_i, rel=2e-15
        )
        assert bound["gamma_labels"][2 * j] == (ell, tuple(beta[1]), "sigma_plus")
        assert bound["gamma_labels"][2 * j + 1] == (ell, tuple(beta[1]), "sigma_minus")

    assert np.allclose(np.sum(bound["eta"] ** 2, axis=-1), 1.0, rtol=0.0, atol=2e-12)
    assert bound["rectangle_geometry_is_autonomous"] is True
    assert bound["source_actual_rectangle_labels_instantiated"] is False


def test_assignment_is_semantically_permutation_invariant() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry()
    partition = _partition_result()
    first = geometry.instantiate(partition)

    n = len(partition["beta_labels"])
    permutation = np.arange(n - 1, -1, -1)
    permuted = dict(partition)
    permuted["eta"] = np.asarray(partition["eta"])[..., permutation]
    permuted["D_r_eta"] = np.asarray(partition["D_r_eta"])[..., permutation]
    permuted["D_z_eta"] = np.asarray(partition["D_z_eta"])[..., permutation]
    permuted["beta_labels"] = tuple(partition["beta_labels"][j] for j in permutation)
    second = geometry.instantiate(permuted)

    def by_label(result: dict) -> dict:
        return {
            label: (
                tuple(int(x) for x in result["center_index_by_beta_sign"][j]),
                float(result["Q_by_beta"][j]),
                float(result["epsilon_by_beta"][j]),
                float(result["L_s_by_beta"][j]),
                int(result["covering_level_by_beta"][j]),
            )
            for j, label in enumerate(result["beta_labels"])
        }

    assert by_label(first) == by_label(second)


def test_exact_witness_certifies_every_transformed_center_constraint() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry()
    rows = geometry.witness.constraint_rows()
    assert rows
    assert all(row["separated"] for row in rows)
    assert min(row["distance_squared_numerator"] / row["distance_squared_denominator"] for row in rows) > 0.0
    assert 4.0 * geometry.witness.r0 * geometry.witness.C_v < 1.0
    assert (
        2.0
        * geometry.witness.r0
        * geometry.witness.C_v
        * (geometry.witness.C_J + 1.0)
        < geometry.witness.d_c
    )


def test_forged_or_incompatible_partition_fails_closed() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry()
    partition = _partition_result()

    duplicate = dict(partition)
    labels = list(partition["beta_labels"])
    labels[-1] = labels[0]
    duplicate["beta_labels"] = tuple(labels)
    with pytest.raises(ValueError, match="unique"):
        geometry.instantiate(duplicate)

    broken_sum = dict(partition)
    broken_sum["eta"] = np.asarray(partition["eta"]) * 0.9
    with pytest.raises(ValueError, match="sum_beta"):
        geometry.instantiate(broken_sum)

    wrong_shape = dict(partition)
    wrong_shape["D_r_eta"] = np.asarray(partition["D_r_eta"])[..., :-1]
    with pytest.raises(ValueError, match="share"):
        geometry.instantiate(wrong_shape)

    below = {
        "eta": np.asarray((1.0,)),
        "D_r_eta": np.asarray((0.0,)),
        "D_z_eta": np.asarray((0.0,)),
        "beta_labels": ((4, (0, 0, 0)),),
    }
    with pytest.raises(ValueError, match="exceed 4"):
        geometry.instantiate(below)


def test_constructor_fails_closed_on_invalid_autonomous_choices() -> None:
    with pytest.raises(ValueError, match="at least two"):
        KokunoAutonomousSignedRectangleGeometry(color_count=1)
    with pytest.raises(ValueError, match="ell0"):
        KokunoAutonomousSignedRectangleGeometry(ell0=4)
    with pytest.raises(ValueError, match="0<h<1/100"):
        KokunoAutonomousSignedRectangleGeometry(h=0.02)
    with pytest.raises(ValueError, match="safety"):
        KokunoAutonomousSignedRectangleGeometry(safety=1.0)


def test_default_witness_matches_source_covering_delta_bound_not_a_hidden_constant() -> None:
    geometry = KokunoAutonomousSignedRectangleGeometry(h=0.005, ell0=5)
    expected = math.ceil(
        KokunoSourceBandCovering.interacting_level_difference_bound(h=geometry.h, ell0=geometry.ell0)
    )
    assert geometry.witness.delta_max == expected
    assert geometry.receipt()["autonomous_conventions"]["signed_pair"].startswith("sigma_plus")
