from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_pulse_end_compensator import (
    CURRENT_LAMBDA,
    KokunoPublicPulseEndCompensator,
)


def test_source_geometry_and_correct_row_slopes_are_frozen() -> None:
    comp = KokunoPublicPulseEndCompensator()
    assert comp.lambda_value == 0.05
    assert comp.s1 == pytest.approx(0.45, rel=0.0, abs=1e-15)
    # Important corrected-source row: s2=1/2-2*lambda, not 3/2-lambda.
    assert comp.s2 == pytest.approx(0.40, rel=0.0, abs=1e-15)
    assert comp.y1 == pytest.approx(257.0)
    assert comp.y2 == pytest.approx(259.0)
    assert comp.width == pytest.approx(0.3)
    assert comp.source_formulas["row_slopes"] == "s1=1/2-lambda, s2=1/2-2lambda"


def test_autonomous_bumps_are_nonnegative_compact_and_distinct() -> None:
    comp = KokunoPublicPulseEndCompensator()
    half = 0.5 * comp.width
    for center, beta in ((comp.y1, comp.beta1), (comp.y2, comp.beta2)):
        assert float(beta(center)) == pytest.approx(1.0, abs=1e-15)
        assert float(beta(center - half)) == 0.0
        assert float(beta(center + half)) == 0.0
        assert float(beta(center - half - 1e-6)) == 0.0
        assert float(beta(center + half + 1e-6)) == 0.0
        samples = beta(np.linspace(center - half, center + half, 41))
        assert np.all(samples >= 0.0)
        assert np.max(samples) > 0.99
    assert float(comp.beta1(comp.y2)) == 0.0
    assert float(comp.beta2(comp.y1)) == 0.0


def test_bump_derivative_matches_centered_difference() -> None:
    comp = KokunoPublicPulseEndCompensator()
    for center, beta, beta_prime in (
        (comp.y1, comp.beta1, comp.beta1_prime),
        (comp.y2, comp.beta2, comp.beta2_prime),
    ):
        y = center + 0.037
        h = 2.0e-7
        fd = (float(beta(y + h)) - float(beta(y - h))) / (2.0 * h)
        assert float(beta_prime(y)) == pytest.approx(fd, rel=2e-7, abs=2e-9)


def test_row_scaled_system_is_finite_and_well_resolved() -> None:
    comp = KokunoPublicPulseEndCompensator()
    matrix = comp.matrix_scaled
    main = comp.main_moments_scaled
    assert matrix.shape == (2, 2)
    assert main.shape == (2,)
    assert np.all(np.isfinite(matrix))
    assert np.all(matrix > 0.0)
    assert np.all(np.isfinite(main))
    assert np.all(main > 0.0)
    # Deterministic receipts for the frozen autonomous bump realization.
    assert matrix[0, 0] == pytest.approx(0.18110027, rel=2e-7)
    assert matrix[0, 1] == pytest.approx(0.44543478, rel=2e-7)
    assert matrix[1, 0] == pytest.approx(0.18108658, rel=2e-7)
    assert matrix[1, 1] == pytest.approx(0.40301559, rel=2e-7)
    assert main[0] == pytest.approx(1.6873990e-8, rel=3e-6)
    assert main[1] == pytest.approx(1.9410552e-7, rel=3e-6)
    assert 1.0 < comp.row_condition_number < 100.0


def test_solver_closes_both_public_rows_for_scalar_and_vector_inputs() -> None:
    comp = KokunoPublicPulseEndCompensator()
    scalar = comp.solve(1.01, 0.0, 0.0)
    assert abs(float(scalar.residual_row1)) < 1e-19
    assert abs(float(scalar.residual_row2)) < 1e-19
    assert float(scalar.c1) != 0.0
    assert float(scalar.c2) != 0.0

    amp = np.array([0.97, 1.01, 1.08])
    m0 = np.array([2e-7, 0.0, -3e-7])
    j0 = np.array([-1e-7, 0.0, 4e-7])
    vec = comp.solve(amp, m0, j0)
    assert vec.c1.shape == (3,)
    assert vec.c2.shape == (3,)
    assert np.max(np.abs(vec.residual_row1)) < 1e-18
    assert np.max(np.abs(vec.residual_row2)) < 1e-18


def test_affine_dependence_on_amplitude_is_exact_to_roundoff() -> None:
    comp = KokunoPublicPulseEndCompensator()
    a0, a1, amid = 0.95, 1.15, 1.05
    s0 = comp.solve(a0, 2e-7, -3e-7)
    s1 = comp.solve(a1, 2e-7, -3e-7)
    sm = comp.solve(amid, 2e-7, -3e-7)
    assert float(sm.c1) == pytest.approx(0.5 * (float(s0.c1) + float(s1.c1)), rel=2e-13, abs=1e-18)
    assert float(sm.c2) == pytest.approx(0.5 * (float(s0.c2) + float(s1.c2)), rel=2e-13, abs=1e-18)


def test_eta_jet_matches_independent_centered_parameter_replay() -> None:
    comp = KokunoPublicPulseEndCompensator()
    amp, m0, j0 = 1.01, 1.7e-7, -2.3e-7
    amp_eta, m_eta, j_eta = 3.0e-4, -4.0e-7, 5.0e-7
    sol, c1_eta, c2_eta = comp.solve_eta_jet(
        amp, m0, j0, amp_eta, m_eta, j_eta
    )
    h = 2e-5
    plus = comp.solve(amp + h * amp_eta, m0 + h * m_eta, j0 + h * j_eta)
    minus = comp.solve(amp - h * amp_eta, m0 - h * m_eta, j0 - h * j_eta)
    fd1 = (float(plus.c1) - float(minus.c1)) / (2.0 * h)
    fd2 = (float(plus.c2) - float(minus.c2)) / (2.0 * h)
    assert float(c1_eta) == pytest.approx(fd1, rel=2e-8, abs=2e-12)
    assert float(c2_eta) == pytest.approx(fd2, rel=2e-8, abs=2e-12)
    assert np.isfinite(float(sol.c1)) and np.isfinite(float(sol.c2))


def test_ratio_correction_derivative_is_analytic_and_localized() -> None:
    comp = KokunoPublicPulseEndCompensator()
    sol = comp.solve(1.01, 0.0, 0.0)
    y = comp.y2 + 0.031
    h = 2e-7
    fd = (
        float(comp.ratio_correction(y + h, sol.c1, sol.c2))
        - float(comp.ratio_correction(y - h, sol.c1, sol.c2))
    ) / (2.0 * h)
    analytic = float(comp.ratio_correction_y(y, sol.c1, sol.c2))
    assert analytic == pytest.approx(fd, rel=3e-7, abs=2e-12)
    assert float(comp.ratio_correction(11.0 / CURRENT_LAMBDA, sol.c1, sol.c2)) == 0.0


def test_truth_boundary_keeps_missing_current_J_and_cartesian_composition_closed() -> None:
    comp = KokunoPublicPulseEndCompensator()
    truth = comp.truth_boundary
    assert truth["public_end_compensation_linear_system_materialized"] is True
    assert truth["repository_autonomous_bump_shape_materialized"] is True
    assert truth["source_exact_bump_shape_recovered"] is False
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["current_lineage_J_entry_materialized"] is False
    assert truth["current_cartesian_end_compensation_composed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_configuration_roundtrip_and_mutations_fail_closed() -> None:
    comp = KokunoPublicPulseEndCompensator()
    config = comp.configuration()
    restored = KokunoPublicPulseEndCompensator.from_configuration(config)
    assert restored.semantic_sha256 == comp.semantic_sha256
    np.testing.assert_allclose(restored.matrix_scaled, comp.matrix_scaled, rtol=0.0, atol=0.0)

    bad_source = copy.deepcopy(config)
    bad_source["source"]["commit"] = "paper-exact-hidden-source"
    with pytest.raises(ValueError):
        KokunoPublicPulseEndCompensator.from_configuration(bad_source)

    bad_width = copy.deepcopy(config)
    bad_width["public_geometry"]["bump_width"] = 0.31
    with pytest.raises(ValueError):
        KokunoPublicPulseEndCompensator.from_configuration(bad_width)

    bad_truth = copy.deepcopy(config)
    bad_truth["truth_boundary"]["current_cartesian_end_compensation_composed"] = True
    with pytest.raises(ValueError):
        KokunoPublicPulseEndCompensator.from_configuration(bad_truth)

    bad_lambda = copy.deepcopy(config)
    bad_lambda["lambda_value"] = 0.051
    with pytest.raises(ValueError):
        KokunoPublicPulseEndCompensator.from_configuration(bad_lambda)

    with pytest.raises(ValueError):
        KokunoPublicPulseEndCompensator(lambda_value=0.051)
