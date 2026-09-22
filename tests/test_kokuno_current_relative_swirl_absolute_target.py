import json
import math
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_relative_swirl_absolute_target import (
    ETA_DERIVATIVE_STEP,
    MAX_LOG_PANEL_WIDTH,
    PARENT_EXACT_HEAD,
    QUADRATURE_ORDER,
    KokunoCurrentRelativeSwirlAbsoluteTarget,
)


def test_absolute_target_materializes_missing_base_without_late_O1_subtraction():
    obj = KokunoCurrentRelativeSwirlAbsoluteTarget()
    eta = np.array([-0.2, 0.0, 0.2])
    report = obj.decomposition_with_eta(eta)

    r = report["absolute_r_I_flatten_end"]
    total = report["total_target_current"]
    correction = report["current_PA17_correction_target"]
    base = report["imported_base_target"]

    assert np.all(np.isfinite(r))
    assert np.all(r > 0.0)
    assert np.all(np.isfinite(total))
    assert np.all(np.isfinite(base))
    assert np.allclose(base + correction, total, rtol=2e-14, atol=0.0)
    assert np.max(np.abs(total)) > 0.0

    truth = obj.truth_boundary
    assert truth["absolute_current_r_I_at_flattening_endpoint_materialized"]
    assert truth["current_relative_swirl_total_target_materialized"]
    assert truth["imported_base_absolute_angular_target_materialized"]
    assert not truth["late_hold_O1_subtraction_used"]
    assert not truth["current_cartesian_relative_swirl_composed"]
    assert not truth["heldout_ns_residual_assessed"]
    assert not truth["pde_validated"]
    assert not truth["paper_exact"]
    assert not truth["openai_field_identified"]


def test_public_stable_target_formula_and_current_solver_closure():
    obj = KokunoCurrentRelativeSwirlAbsoluteTarget()
    eta = 0.17
    r_flat = float(obj.absolute_r_flat(eta))
    expected = -(r_flat - obj.r_star) * math.exp(
        -(1.0 - obj.lambda_value) * obj.compensator.y1
    )
    actual = float(obj.total_target(eta))
    assert actual == pytest.approx(expected, rel=2e-15, abs=0.0)

    solution, c1_eta, c2_eta = obj.solve_current_target(eta)
    assert math.isfinite(float(solution.c1))
    assert math.isfinite(float(solution.c2))
    assert math.isfinite(float(c1_eta))
    assert math.isfinite(float(c2_eta))
    assert abs(float(solution.angular_residual)) <= 2e-14
    assert abs(float(solution.pressure_residual)) <= 2e-14


def test_eta_jet_matches_independent_coarser_centered_difference():
    obj = KokunoCurrentRelativeSwirlAbsoluteTarget()
    eta = 0.13
    target, target_eta = obj.target_with_eta(eta)
    h = 4.0e-5
    independent = (
        -float(obj.total_target(eta + 2*h))
        + 8.0*float(obj.total_target(eta + h))
        - 8.0*float(obj.total_target(eta - h))
        + float(obj.total_target(eta - 2*h))
    ) / (12.0*h)
    scale = max(abs(independent), abs(float(target_eta)), 1e-300)
    assert abs(float(target_eta) - independent) / scale <= 5e-5
    assert math.isfinite(float(target))


def test_axis_anchor_is_actual_current_prefix_I_not_zero_or_fixed_point():
    obj = KokunoCurrentRelativeSwirlAbsoluteTarget()
    r0, e0 = obj._prefix_anchor(0.0)
    physical = obj.prefix.physical_prefix_moments_at_Xi(0.0)
    values = obj.prefix.profile_values(obj.prefix.X_i, 0.0)
    E = float(np.asarray(values["E_actual_prefix"]))
    H = math.sqrt(2.0 * obj.prefix.X_i) * E
    direct = float(physical[1]) / (obj.prefix.X_i * H)
    assert r0 == pytest.approx(direct, rel=2e-14, abs=0.0)
    assert e0 == pytest.approx(E, rel=0.0, abs=0.0)
    assert float(physical[1]) != 0.0


def test_configuration_save_load_and_provenance_fail_closed(tmp_path: Path):
    obj = KokunoCurrentRelativeSwirlAbsoluteTarget()
    cfg = obj.configuration()
    assert cfg["parent_exact_head"] == PARENT_EXACT_HEAD
    assert cfg["quadrature_order"] == QUADRATURE_ORDER
    assert cfg["max_log_panel_width"] == MAX_LOG_PANEL_WIDTH
    assert cfg["eta_derivative_step"] == ETA_DERIVATIVE_STEP
    assert cfg["source_commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"

    path = tmp_path / "target.json"
    obj.save_configuration(path)
    loaded = KokunoCurrentRelativeSwirlAbsoluteTarget.load_configuration(path)
    assert loaded.semantic_sha256() == obj.semantic_sha256()

    bad = json.loads(path.read_text())
    bad["max_log_panel_width"] = 0.5
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoCurrentRelativeSwirlAbsoluteTarget.from_configuration(bad)


def test_api_has_no_residual_forcing_or_optimizer_knob():
    import inspect

    sig = inspect.signature(KokunoCurrentRelativeSwirlAbsoluteTarget.total_target)
    forbidden = {
        "residual", "forcing", "pressure", "viscosity", "optimizer",
        "loss", "gain", "threshold", "target_value",
    }
    assert forbidden.isdisjoint(sig.parameters)
