from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)
from openai_ns_reconstruction.kokuno_candidate_mass_bound_phase_curl import (
    KokunoCandidateMassBoundPhaseCurlFamily,
)
from openai_ns_reconstruction.kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)
from openai_ns_reconstruction.kokuno_source_signed_covariance_mass import (
    KokunoSourceSignedCovarianceMass,
)


def _geometry_and_mass(h: float = 0.005) -> tuple[dict, dict]:
    partition = KokunoSourceCompatiblePartitionRealization(ell_min=5).evaluate(
        q=np.asarray(2.0**-5.5),
        D_r_q=np.asarray(0.0),
        D_z_q=np.asarray(0.0),
        slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
    )
    geometry = KokunoAutonomousSignedRectangleGeometry(h=h).instantiate(partition)
    v_by_beta = []
    psi_by_beta = []
    x_by_beta_sign = []
    for L_s in geometry["L_s_by_beta"]:
        v = np.linspace(0.0, float(L_s), 129)
        v_by_beta.append(v)
        psi_by_beta.append(np.ones_like(v))
        x_by_beta_sign.append(
            np.stack((np.full_like(v, 2.0), np.full_like(v, 3.0)), axis=-1)
        )
    mass = KokunoSourceSignedCovarianceMass(h=h).materialize_from_autonomous_geometry(
        geometry,
        v_by_beta=tuple(v_by_beta),
        psi_by_beta=tuple(psi_by_beta),
        x_by_beta_sign=tuple(x_by_beta_sign),
        xi=np.linspace(-1.0, 1.0, 257),
        chi_g=np.ones(257),
        pulse_binding="repository_autonomous_source_compatible",
    )
    return geometry, mass


def _phase_and_curl_inputs(family: KokunoCandidateMassBoundPhaseCurlFamily, geometry: dict) -> dict:
    labels = tuple(geometry["beta_labels"])
    n = len(labels)
    L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
    R = np.linspace(0.72, 0.92, n)
    Z = np.linspace(0.08, -0.06, n)
    u_star = np.linspace(2.0, 2.2, n)
    phase_inputs = {
        "R": R,
        "Z": Z,
        "theta": 0.37,
        "v": np.stack((0.25 * L_s, 0.75 * L_s), axis=-1),
        "beta_labels": labels,
        "R0": np.full(n, 0.82),
        "F0": np.full(n, 1.2),
        "a": np.full(n, 3.2),
        "b_s": np.full(n, 0.2),
        "u_star": u_star,
        "L_s": L_s,
        "F": np.full(n, 1.1),
        "G": np.full(n, 0.24),
        "F_R": np.full(n, 0.08),
        "G_R": np.full(n, -0.05),
        "F_Z": np.full(n, 0.04),
        "G_Z": np.full(n, 0.03),
    }
    frame = family.phase_family.phase_frame(**phase_inputs)
    n_phi = frame["n_phi"]
    base = np.stack(
        (n_phi[..., 1], -n_phi[..., 0], np.zeros_like(n_phi[..., 0])), axis=-1
    )
    t_plus = (1.0 + 0.17j) * base.astype(np.complex128)
    A_c = -frame["c0_by_beta"] * np.sqrt(1.0 + u_star * u_star)
    eta = np.full(n, 1.0 / np.sqrt(float(n)))
    return {
        **phase_inputs,
        "t_plus_prototype": t_plus,
        "D_r_C_plus_prototype": np.zeros_like(t_plus),
        "D_z_C_plus_prototype": np.zeros_like(t_plus),
        "eta": eta,
        "D_r_eta": np.zeros(n),
        "D_z_eta": np.zeros(n),
        "A_c": A_c,
        "T_N": -0.40 * A_c,
        "T_K": 0.05 * u_star,
    }


def test_mass_bridge_matches_direct_phase_curl_with_same_h_and_derivatives() -> None:
    geometry, mass = _geometry_and_mass()
    bridge = KokunoCandidateMassBoundPhaseCurlFamily()
    inputs = _phase_and_curl_inputs(bridge, geometry)
    masses = np.asarray(mass["h_sigma_by_beta_sign"], dtype=float)
    D_r_h = 1.0e-4 * masses
    D_z_h = -2.0e-4 * masses

    bridged = bridge.physical_family_from_materialized_mass(
        mass,
        D_r_h_by_beta_sign=D_r_h,
        D_z_h_by_beta_sign=D_z_h,
        h_derivative_binding="candidate_derived_from_same_pulse_model",
        **inputs,
    )
    direct = bridge.phase_family.physical_family(
        **inputs,
        h_plus=masses[:, 0],
        h_minus=masses[:, 1],
        D_r_h_plus=D_r_h[:, 0],
        D_r_h_minus=D_r_h[:, 1],
        D_z_h_plus=D_z_h[:, 0],
        D_z_h_minus=D_z_h[:, 1],
    )

    np.testing.assert_array_equal(
        bridged["velocity_physical_cartesian_by_beta_sign"],
        direct["velocity_physical_cartesian_by_beta_sign"],
    )
    np.testing.assert_array_equal(
        bridged["velocity_physical_cartesian_by_beta"],
        direct["velocity_physical_cartesian_by_beta"],
    )
    np.testing.assert_array_equal(
        bridged["velocity_physical_cartesian_total"],
        direct["velocity_physical_cartesian_total"],
    )
    np.testing.assert_array_equal(
        bridged["candidate_velocity_osc_cartesian_total"],
        bridged["velocity_physical_cartesian_total"],
    )
    assert np.linalg.norm(bridged["candidate_velocity_osc_cartesian_total"]) > 0.0
    assert bridged["candidate_h_sigma_bound_into_complete_curl"] is True
    assert bridged["h_sigma_derivative_binding"] == "candidate_derived_from_same_pulse_model"
    assert bridged["public_xyz_t_velocity_correction_materialized"] is False
    assert bridged["paper_exact"] is False
    assert bridged["pde_validated"] is False


def test_frozen_local_binding_requires_exact_zero_mass_derivatives() -> None:
    geometry, mass = _geometry_and_mass()
    bridge = KokunoCandidateMassBoundPhaseCurlFamily()
    inputs = _phase_and_curl_inputs(bridge, geometry)
    n = len(geometry["beta_labels"])

    out = bridge.physical_family_from_materialized_mass(
        mass,
        D_r_h_by_beta_sign=np.zeros((n, 2)),
        D_z_h_by_beta_sign=np.zeros((n, 2)),
        h_derivative_binding="explicit_frozen_local_zero",
        **inputs,
    )
    assert np.all(np.isfinite(out["candidate_velocity_osc_cartesian_total"]))

    bad = np.zeros((n, 2))
    bad[0, 0] = 1e-12
    with pytest.raises(ValueError, match="exactly zero"):
        bridge.physical_family_from_materialized_mass(
            mass,
            D_r_h_by_beta_sign=bad,
            D_z_h_by_beta_sign=np.zeros((n, 2)),
            h_derivative_binding="explicit_frozen_local_zero",
            **inputs,
        )


def test_mass_schedule_labels_and_provenance_fail_closed() -> None:
    geometry, mass = _geometry_and_mass()
    bridge = KokunoCandidateMassBoundPhaseCurlFamily()
    inputs = _phase_and_curl_inputs(bridge, geometry)
    n = len(geometry["beta_labels"])
    zeros = np.zeros((n, 2))

    tampered_q = dict(mass)
    tampered_q["Q_by_beta"] = np.asarray(mass["Q_by_beta"], dtype=float).copy()
    tampered_q["Q_by_beta"][0] *= 1.01
    with pytest.raises(ValueError, match="source band schedule"):
        bridge.physical_family_from_materialized_mass(
            tampered_q,
            D_r_h_by_beta_sign=zeros,
            D_z_h_by_beta_sign=zeros,
            h_derivative_binding="explicit_frozen_local_zero",
            **inputs,
        )

    tampered_truth = dict(mass)
    tampered_truth["paper_exact"] = True
    with pytest.raises(ValueError, match="paper_exact"):
        bridge.physical_family_from_materialized_mass(
            tampered_truth,
            D_r_h_by_beta_sign=zeros,
            D_z_h_by_beta_sign=zeros,
            h_derivative_binding="explicit_frozen_local_zero",
            **inputs,
        )

    bad_inputs = dict(inputs)
    bad_labels = list(inputs["beta_labels"])
    bad_labels[0] = (int(bad_labels[0][0]), "changed-label")
    bad_inputs["beta_labels"] = tuple(bad_labels)
    with pytest.raises(ValueError, match="exactly match"):
        bridge.physical_family_from_materialized_mass(
            mass,
            D_r_h_by_beta_sign=zeros,
            D_z_h_by_beta_sign=zeros,
            h_derivative_binding="explicit_frozen_local_zero",
            **bad_inputs,
        )


def test_mass_values_cannot_be_overridden_by_caller() -> None:
    geometry, mass = _geometry_and_mass()
    bridge = KokunoCandidateMassBoundPhaseCurlFamily()
    inputs = _phase_and_curl_inputs(bridge, geometry)
    n = len(geometry["beta_labels"])
    inputs["h_plus"] = np.ones(n)
    with pytest.raises(ValueError, match="owned by this bridge"):
        bridge.physical_family_from_materialized_mass(
            mass,
            D_r_h_by_beta_sign=np.zeros((n, 2)),
            D_z_h_by_beta_sign=np.zeros((n, 2)),
            h_derivative_binding="explicit_frozen_local_zero",
            **inputs,
        )


def test_receipt_keeps_candidate_and_source_truth_separate() -> None:
    receipt = KokunoCandidateMassBoundPhaseCurlFamily().receipt()
    assert receipt["schema"] == "kokuno-candidate-mass-bound-phase-curl-v1"
    assert "explicit candidate inputs" in receipt["numerical_boundary"]["h_sigma_derivatives"]
    truth = receipt["truth_boundary"]
    assert truth["candidate_h_sigma_bound_into_complete_curl"] is True
    assert truth["candidate_velocity_osc_arrays_materialized"] is True
    assert truth["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
