from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)
from openai_ns_reconstruction.kokuno_physical_evaluation import SOURCE_B_G
from openai_ns_reconstruction.kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)
from openai_ns_reconstruction.kokuno_source_signed_covariance_mass import (
    KokunoSourceSignedCovarianceMass,
)


def _geometry_result(h: float = 0.005) -> dict:
    partition = KokunoSourceCompatiblePartitionRealization(ell_min=5).evaluate(
        q=np.asarray(2.0**-5.5),
        D_r_q=np.asarray(0.0),
        D_z_q=np.asarray(0.0),
        slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
    )
    return KokunoAutonomousSignedRectangleGeometry(h=h).instantiate(partition)


def _constant_pulses(geometry: dict) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    grids = []
    cutoffs = []
    pulses = []
    for L_s in geometry["L_s_by_beta"]:
        v = np.linspace(0.0, float(L_s), 257)
        grids.append(v)
        cutoffs.append(np.ones_like(v))
        pulses.append(np.stack((np.full_like(v, 2.0), np.full_like(v, 3.0)), axis=-1))
    return tuple(grids), tuple(cutoffs), tuple(pulses)


def test_displayed_h_sigma_formula_matches_independent_constant_integral() -> None:
    h = 0.005
    geometry = _geometry_result(h)
    v, psi, x = _constant_pulses(geometry)
    xi = np.linspace(-1.0, 1.0, 513)
    chi = np.ones_like(xi)

    evaluator = KokunoSourceSignedCovarianceMass(h=h)
    out = evaluator.materialize_from_autonomous_geometry(
        geometry,
        v_by_beta=v,
        psi_by_beta=psi,
        x_by_beta_sign=x,
        xi=xi,
        chi_g=chi,
        pulse_binding="repository_autonomous_source_compatible",
    )

    D_g = 1.0 + SOURCE_B_G**2  # half factor times integral_{-1}^{1} 1 dxi
    r0 = float(geometry["r0"])
    expected_plus = D_g * 8.0 * r0
    expected_minus = D_g * 18.0 * r0

    assert out["D_g"] == pytest.approx(D_g, rel=2e-15)
    assert np.allclose(out["h_sigma_by_beta_sign"][:, 0], expected_plus, rtol=2e-13, atol=0.0)
    assert np.allclose(out["h_sigma_by_beta_sign"][:, 1], expected_minus, rtol=2e-13, atol=0.0)

    for j, L_s in enumerate(geometry["L_s_by_beta"]):
        assert out["I0_by_beta_sign"][j, 0] == pytest.approx(4.0 * float(L_s), rel=2e-14)
        assert out["I0_by_beta_sign"][j, 1] == pytest.approx(9.0 * float(L_s), rel=2e-14)


def test_candidate_masses_feed_existing_signed_reference_columns_without_rank_promotion() -> None:
    geometry = _geometry_result()
    v, psi, x = _constant_pulses(geometry)
    evaluator = KokunoSourceSignedCovarianceMass()
    out = evaluator.materialize_from_autonomous_geometry(
        geometry,
        v_by_beta=v,
        psi_by_beta=psi,
        x_by_beta_sign=x,
        xi=np.asarray((-1.0, 0.0, 1.0)),
        chi_g=np.asarray((1.0, 1.0, 1.0)),
        pulse_binding="candidate_derived_from_executable_background",
    )
    columns = evaluator.reference_columns(out, A_c=1.7, u_star=0.8)
    H = columns["reference_column_matrix_by_beta"]
    masses = out["h_sigma_by_beta_sign"]

    assert H.shape == (len(geometry["beta_labels"]), 2, 2)
    assert np.allclose(H[:, 0, 0], -1.7 * masses[:, 0], rtol=0.0, atol=1e-15)
    assert np.allclose(H[:, 0, 1], -1.7 * masses[:, 1], rtol=0.0, atol=1e-15)
    assert np.allclose(H[:, 1, 0], -0.8 * masses[:, 0], rtol=0.0, atol=1e-15)
    assert np.allclose(H[:, 1, 1], +0.8 * masses[:, 1], rtol=0.0, atol=1e-15)
    assert columns["reference_formula_only"] is True
    assert columns["actual_complete_curl_covariance_rank_assessed"] is False
    assert columns["actual_source_mode_family_bound"] is False


def test_geometry_schedule_is_revalidated_before_mass_materialization() -> None:
    geometry = _geometry_result()
    v, psi, x = _constant_pulses(geometry)
    evaluator = KokunoSourceSignedCovarianceMass()

    tampered = dict(geometry)
    tampered["L_s_by_beta"] = np.asarray(geometry["L_s_by_beta"], dtype=float).copy()
    tampered["L_s_by_beta"][0] *= 1.01
    with pytest.raises(ValueError, match="L_s"):
        evaluator.materialize_from_autonomous_geometry(
            tampered,
            v_by_beta=v,
            psi_by_beta=psi,
            x_by_beta_sign=x,
            xi=np.asarray((-1.0, 0.0, 1.0)),
            chi_g=np.ones(3),
            pulse_binding="repository_autonomous_source_compatible",
        )


def test_pulse_and_cutoff_inputs_fail_closed() -> None:
    geometry = _geometry_result()
    v, psi, x = _constant_pulses(geometry)
    evaluator = KokunoSourceSignedCovarianceMass()

    bad_psi = list(psi)
    bad_psi[0] = bad_psi[0].copy()
    bad_psi[0][7] = 1.2
    with pytest.raises(ValueError, match="cutoff range"):
        evaluator.materialize_from_autonomous_geometry(
            geometry,
            v_by_beta=v,
            psi_by_beta=tuple(bad_psi),
            x_by_beta_sign=x,
            xi=np.asarray((-1.0, 0.0, 1.0)),
            chi_g=np.ones(3),
            pulse_binding="repository_autonomous_source_compatible",
        )

    bad_x = list(x)
    bad_x[0] = bad_x[0].copy()
    bad_x[0][3, 0] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        evaluator.materialize_from_autonomous_geometry(
            geometry,
            v_by_beta=v,
            psi_by_beta=psi,
            x_by_beta_sign=tuple(bad_x),
            xi=np.asarray((-1.0, 0.0, 1.0)),
            chi_g=np.ones(3),
            pulse_binding="repository_autonomous_source_compatible",
        )

    bad_grid = list(v)
    bad_grid[0] = bad_grid[0].copy()
    bad_grid[0][-1] *= 0.99
    with pytest.raises(ValueError, match="end at geometry L_s"):
        evaluator.materialize_from_autonomous_geometry(
            geometry,
            v_by_beta=tuple(bad_grid),
            psi_by_beta=psi,
            x_by_beta_sign=x,
            xi=np.asarray((-1.0, 0.0, 1.0)),
            chi_g=np.ones(3),
            pulse_binding="repository_autonomous_source_compatible",
        )

    with pytest.raises(ValueError, match="pulse_binding"):
        evaluator.materialize_from_autonomous_geometry(
            geometry,
            v_by_beta=v,
            psi_by_beta=psi,
            x_by_beta_sign=x,
            xi=np.asarray((-1.0, 0.0, 1.0)),
            chi_g=np.ones(3),
            pulse_binding="source_exact",
        )


def test_receipt_separates_source_formula_from_candidate_quadrature() -> None:
    receipt = KokunoSourceSignedCovarianceMass().receipt()
    assert receipt["schema"] == "kokuno-source-signed-covariance-mass-v1"
    assert receipt["source"]["release_date"] == "2026-09-09"
    assert "h_sigma=D_g*c_i*I0_sigma" in receipt["source_formulas"]["signed_mass"]
    assert receipt["numerical_contract"]["source_hidden_numeric_values_are_not_inferred"] is True
    truth = receipt["truth_boundary"]
    assert truth["source_h_sigma_integral_formula_executable"] is True
    assert truth["candidate_h_sigma_numerically_materializable"] is True
    assert truth["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert truth["actual_source_pulse_samples_recovered"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
