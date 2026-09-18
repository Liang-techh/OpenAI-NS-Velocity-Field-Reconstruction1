from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_localized_real_pair_family import (
    KokunoSourceLocalizedRealPairFamily,
)
from openai_ns_reconstruction.kokuno_source_support_localized_curl import (
    KokunoSourceSupportLocalizedCurl,
)


def _manufactured_family():
    R = np.array([0.82, 1.01, 1.19])
    theta = np.array([0.17, -0.31, 0.42])
    alpha = 0.23 + 0.15 * R
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.15 * tangent
    D_z_eta = np.zeros_like(eta)
    phase = R[:, None] + np.array([0.13, -0.27])[None, :]

    n_phi = np.zeros((3, 2, 3), dtype=float)
    n_phi[..., 0] = 1.0
    t_plus = np.empty((3, 2, 3), dtype=np.complex128)
    t_plus[:, 0, :] = np.array([0.0, 0.24 + 0.07j, -0.11 + 0.04j])
    t_plus[:, 1, :] = np.array([0.0, -0.18 + 0.05j, 0.16 - 0.03j])
    D_r_C_plus = np.zeros_like(t_plus)
    D_z_C_plus = np.zeros_like(t_plus)
    Q = np.array([[0.5, 0.25], [0.5, 0.25], [0.5, 0.25]])
    labels = [(0, "cell-a"), (1, "cell-b")]
    return {
        "R": R,
        "theta": theta,
        "eta": eta,
        "D_r_eta": D_r_eta,
        "D_z_eta": D_z_eta,
        "phase": phase,
        "n_phi": n_phi,
        "t_plus": t_plus,
        "D_r_C_plus": D_r_C_plus,
        "D_z_C_plus": D_z_C_plus,
        "Q": Q,
        "labels": labels,
    }


def test_family_matches_explicit_per_beta_localized_pair_sum():
    family = KokunoSourceLocalizedRealPairFamily(epsilon=0.2, h=0.004)
    d = _manufactured_family()
    out = family.physical_family(
        d["Q"], d["R"], d["theta"], d["phase"], d["n_phi"],
        d["t_plus"], d["D_r_C_plus"], d["D_z_C_plus"],
        d["eta"], d["D_r_eta"], d["D_z_eta"], d["labels"],
    )

    plus = KokunoSourceSupportLocalizedCurl(epsilon=0.2, m=1)
    minus = KokunoSourceSupportLocalizedCurl(epsilon=0.2, m=-1)
    manual_cyl = np.zeros((3, 3), dtype=float)
    for b in range(2):
        up = plus.localized_mode(
            d["R"], d["phase"][:, b], d["n_phi"][:, b, :],
            d["t_plus"][:, b, :], d["D_r_C_plus"][:, b, :],
            d["D_z_C_plus"][:, b, :], d["eta"][:, b],
            d["D_r_eta"][:, b], d["D_z_eta"][:, b],
        )["velocity"]
        um = minus.localized_mode(
            d["R"], d["phase"][:, b], d["n_phi"][:, b, :],
            np.conjugate(d["t_plus"][:, b, :]),
            np.conjugate(d["D_r_C_plus"][:, b, :]),
            np.conjugate(d["D_z_C_plus"][:, b, :]),
            d["eta"][:, b], d["D_r_eta"][:, b], d["D_z_eta"][:, b],
        )["velocity"]
        manual_cyl += d["Q"][:, b, None] ** (-family.A) * (up + um).real

    c = np.cos(d["theta"])
    s = np.sin(d["theta"])
    manual_cart = np.stack(
        (
            manual_cyl[:, 0] * c - manual_cyl[:, 1] * s,
            manual_cyl[:, 0] * s + manual_cyl[:, 1] * c,
            manual_cyl[:, 2],
        ),
        axis=-1,
    )
    np.testing.assert_allclose(
        out["velocity_physical_cylindrical_total"], manual_cyl, rtol=2e-13, atol=2e-13
    )
    np.testing.assert_allclose(
        out["velocity_physical_cartesian_total"], manual_cart, rtol=2e-13, atol=2e-13
    )


def test_family_is_real_and_preserves_whole_partition_diagnostics():
    family = KokunoSourceLocalizedRealPairFamily(epsilon=0.2)
    d = _manufactured_family()
    out = family.normalized_family(
        d["R"], d["phase"], d["n_phi"], d["t_plus"],
        d["D_r_C_plus"], d["D_z_C_plus"], d["eta"],
        d["D_r_eta"], d["D_z_eta"], d["labels"],
    )
    np.testing.assert_allclose(
        out["velocity_cylindrical_by_beta"],
        2.0 * out["plus_velocity_by_beta"].real,
        rtol=0.0,
        atol=3e-12,
    )
    assert np.max(np.abs(out["partition_error"])) < 1e-14
    assert np.max(np.abs(out["partition_radial_closure"])) < 1e-14
    assert np.max(np.abs(out["partition_axial_closure"])) < 1e-14
    assert np.linalg.norm(out["velocity_cylindrical_total"]) > 1e-3


def test_family_fails_closed_on_partition_or_label_misuse():
    family = KokunoSourceLocalizedRealPairFamily(epsilon=0.2)
    d = _manufactured_family()
    bad_eta = d["eta"].copy()
    bad_eta[:, 0] *= 0.9
    with pytest.raises(ValueError, match="squared partition"):
        family.normalized_family(
            d["R"], d["phase"], d["n_phi"], d["t_plus"],
            d["D_r_C_plus"], d["D_z_C_plus"], bad_eta,
            d["D_r_eta"], d["D_z_eta"], d["labels"],
        )
    with pytest.raises(ValueError, match="do not duplicate"):
        family.normalized_family(
            d["R"], d["phase"], d["n_phi"], d["t_plus"],
            d["D_r_C_plus"], d["D_z_C_plus"], d["eta"],
            d["D_r_eta"], d["D_z_eta"], [(0, "cell-a", "+"), (1, "cell-b", "-")],
        )
    with pytest.raises(ValueError, match="Q_beta"):
        family.physical_family(
            np.array([[1.2, 0.25]] * 3), d["R"], d["theta"], d["phase"],
            d["n_phi"], d["t_plus"], d["D_r_C_plus"], d["D_z_C_plus"],
            d["eta"], d["D_r_eta"], d["D_z_eta"], d["labels"],
        )


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    family = KokunoSourceLocalizedRealPairFamily(epsilon=0.125, h=0.003)
    path = family.save_json(tmp_path / "localized_real_pair_family.json")
    loaded = KokunoSourceLocalizedRealPairFamily.load_json(path)
    assert loaded.sha256 == family.sha256
    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["supplied_mode_family_composition_executable"] is True
    assert truth["source_real_conjugate_harmonic_reconstruction_executable"] is True
    assert truth["concrete_source_partition_bumps_reconstructed"] is False
    assert truth["source_actual_background_path_instantiated"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["pde_validated"] is False
    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceLocalizedRealPairFamily.from_payload(payload)
