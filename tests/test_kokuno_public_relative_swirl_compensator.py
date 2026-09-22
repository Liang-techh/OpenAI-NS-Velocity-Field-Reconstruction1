import copy
import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_relative_swirl_compensator import (
    CURRENT_LAMBDA,
    PARENT_EXACT_HEAD,
    PUBLIC_BUMP_WIDTH,
    SOURCE_BLOB,
    SOURCE_COMMIT,
    KokunoPublicRelativeSwirlCompensator,
)


def test_public_geometry_and_corrected_source_binding():
    comp = KokunoPublicRelativeSwirlCompensator()
    expected_hold = 30.0 * math.log(1.0 / CURRENT_LAMBDA)
    assert comp.hold_length == pytest.approx(expected_hold, abs=1.0e-14)
    assert comp.y1 == pytest.approx(expected_hold - 3.0, abs=1.0e-14)
    assert comp.y2 == pytest.approx(expected_hold - 1.0, abs=1.0e-14)
    assert comp.width == PUBLIC_BUMP_WIDTH == 0.3
    assert comp.angular_slope == pytest.approx(1.0 - CURRENT_LAMBDA)
    assert comp.pressure_slope == pytest.approx(-1.0 - 2.0 * CURRENT_LAMBDA)

    cfg = comp.configuration()
    assert cfg["source_commit"] == SOURCE_COMMIT
    assert cfg["source_blob"] == SOURCE_BLOB
    assert cfg["parent_exact_head"] == PARENT_EXACT_HEAD
    assert cfg["source_release_date"] == "2026-09-09"
    assert "I=XH/(1-lambda)" in cfg["source_formulas"]["angular_target"]
    assert "E^2-E_unedited^2" in cfg["source_formulas"]["pressure_target"]


def test_autonomous_bumps_are_compact_nonnegative_and_have_analytic_derivative():
    comp = KokunoPublicRelativeSwirlCompensator()
    half = 0.5 * comp.width
    probes = np.array(
        [comp.y1 - half - 1.0e-3, comp.y1, comp.y1 + half + 1.0e-3]
    )
    beta = comp.beta1(probes)
    assert beta[0] == 0.0
    assert beta[1] == pytest.approx(1.0)
    assert beta[2] == 0.0
    assert np.all(beta >= 0.0)

    y = comp.y1 + 0.037
    h = 1.0e-6
    fd = float((comp.beta1(y + h) - comp.beta1(y - h)) / (2.0 * h))
    analytic = float(comp.beta1_prime(y))
    assert analytic == pytest.approx(fd, rel=2.0e-8, abs=2.0e-8)

    truth = comp.truth_boundary
    assert truth["repository_autonomous_pointwise_bump_shape_materialized"] is True
    assert truth["source_exact_pointwise_bump_shape_recovered"] is False


def test_scaled_rows_are_nondegenerate_and_deterministic():
    comp = KokunoPublicRelativeSwirlCompensator()
    np.testing.assert_allclose(
        comp.angular_row_scaled,
        np.array([0.18132583737155128, 1.2123254083240895]),
        rtol=0.0,
        atol=3.0e-14,
    )
    np.testing.assert_allclose(
        comp.pressure_linear_row_scaled,
        np.array([0.18142499039700033, 0.02010246194184373]),
        rtol=0.0,
        atol=3.0e-14,
    )
    np.testing.assert_allclose(
        comp.pressure_quadratic_matrix_scaled,
        np.array(
            [[0.14773802507690695, 0.0], [0.0, 0.01636983978873497]]
        ),
        rtol=0.0,
        atol=4.0e-14,
    )
    assert comp.linearized_condition_number == pytest.approx(
        3.4954201700562058, rel=1.0e-13
    )


def test_zero_target_selects_continuous_zero_branch():
    comp = KokunoPublicRelativeSwirlCompensator()
    sol = comp.solve(0.0)
    assert float(sol.c1) == pytest.approx(0.0, abs=1.0e-15)
    assert float(sol.c2) == pytest.approx(0.0, abs=1.0e-15)
    assert float(sol.angular_residual) == pytest.approx(0.0, abs=1.0e-15)
    assert float(sol.pressure_residual) == pytest.approx(0.0, abs=1.0e-15)
    assert float(sol.discriminant) > 0.0


def test_vectorized_solver_hits_exact_angular_and_quadratic_pressure_rows():
    comp = KokunoPublicRelativeSwirlCompensator()
    target = np.array([-2.0e-7, 0.0, 1.0e-7, 1.0e-4])
    sol = comp.solve(target)
    assert sol.c1.shape == target.shape
    assert sol.c2.shape == target.shape
    np.testing.assert_allclose(sol.angular_residual, 0.0, rtol=0.0, atol=3.0e-16)
    np.testing.assert_allclose(sol.pressure_residual, 0.0, rtol=0.0, atol=8.0e-16)
    assert np.all(np.abs(sol.c1) + np.abs(sol.c2) < 1.0)

    # Independent direct algebra replay of the exact pressure row.
    for k in range(target.size):
        c = np.array([sol.c1[k], sol.c2[k]])
        pressure = (
            2.0 * comp.pressure_linear_row_scaled @ c
            + c @ comp.pressure_quadratic_matrix_scaled @ c
        )
        assert pressure == pytest.approx(0.0, abs=8.0e-16)
        assert comp.angular_row_scaled @ c == pytest.approx(target[k], abs=3.0e-16)


def test_target_derivative_and_eta_jet_are_analytic():
    comp = KokunoPublicRelativeSwirlCompensator()
    target = 1.0e-5
    sol = comp.solve(target)
    h = 1.0e-8
    plus = comp.solve(target + h)
    minus = comp.solve(target - h)
    fd1 = float((plus.c1 - minus.c1) / (2.0 * h))
    fd2 = float((plus.c2 - minus.c2) / (2.0 * h))
    assert float(sol.dc1_dtarget) == pytest.approx(fd1, rel=2.0e-8, abs=2.0e-10)
    assert float(sol.dc2_dtarget) == pytest.approx(fd2, rel=2.0e-8, abs=2.0e-10)

    target_eta = -3.25e-6
    _, c1_eta, c2_eta = comp.solve_target_jet(target, target_eta)
    assert float(c1_eta) == pytest.approx(float(sol.dc1_dtarget) * target_eta)
    assert float(c2_eta) == pytest.approx(float(sol.dc2_dtarget) * target_eta)


def test_rI_helper_does_not_invent_entry_state_and_recovers_fixed_point():
    comp = KokunoPublicRelativeSwirlCompensator()
    r_star = 1.0 / (1.0 - CURRENT_LAMBDA)
    assert float(comp.angular_target_from_rI(r_star)) == pytest.approx(0.0, abs=1.0e-14)

    r0 = r_star + 2.0e-6
    y0 = comp.hold_length - 3.5
    alpha = 1.0 - CURRENT_LAMBDA
    r_unedited_end = r_star + (r0 - r_star) * math.exp(
        -alpha * (comp.hold_length - y0)
    )
    expected = math.exp(alpha * 3.0) * (r_star - r_unedited_end)
    assert float(comp.angular_target_from_rI(r0)) == pytest.approx(expected, rel=2e-12)

    with pytest.raises(ValueError):
        comp.angular_target_from_rI(r0, y_entry=comp.y1)


def test_relative_edit_derivative_and_E_positivity():
    comp = KokunoPublicRelativeSwirlCompensator()
    sol = comp.solve(1.0e-4)
    y = comp.y2 + 0.031
    h = 1.0e-6
    fd = float(
        (
            comp.relative_edit(y + h, sol.c1, sol.c2)
            - comp.relative_edit(y - h, sol.c1, sol.c2)
        )
        / (2.0 * h)
    )
    analytic = float(comp.relative_edit_y(y, sol.c1, sol.c2))
    assert analytic == pytest.approx(fd, rel=2.0e-8, abs=2.0e-10)

    e = comp.edited_E(np.array([2.0, 3.0]), np.array([comp.y1, comp.y2]), sol.c1, sol.c2)
    assert np.all(np.isfinite(e))
    assert np.all(e > 0.0)


def test_save_load_semantic_identity_and_provenance_mutation_fail_closed(tmp_path):
    comp = KokunoPublicRelativeSwirlCompensator()
    path = tmp_path / "relative_swirl.json"
    comp.save_configuration(path)
    loaded = KokunoPublicRelativeSwirlCompensator.load_configuration(path)
    assert loaded.semantic_sha256() == comp.semantic_sha256()
    assert loaded.configuration() == comp.configuration()

    payload = json.loads(path.read_text())
    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["source_exact_pointwise_bump_shape_recovered"] = True
    with pytest.raises(ValueError):
        KokunoPublicRelativeSwirlCompensator.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bump_width"] = 0.31
    with pytest.raises(ValueError):
        KokunoPublicRelativeSwirlCompensator.from_configuration(mutated)


def test_truth_boundary_and_public_api_expose_no_residual_tuning_surface():
    comp = KokunoPublicRelativeSwirlCompensator()
    truth = comp.truth_boundary
    assert truth["public_relative_swirl_two_bump_algebra_materialized"] is True
    assert truth["current_lineage_angular_entry_I_materialized"] is False
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["complete_terminal_hold_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "optimizer",
        "heldout",
        "threshold",
        "target_residual",
    }
    for method_name in ("solve", "solve_target_jet", "angular_target_from_rI"):
        params = set(inspect.signature(getattr(comp, method_name)).parameters)
        assert not (params & forbidden)
