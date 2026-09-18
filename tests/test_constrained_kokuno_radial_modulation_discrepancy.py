import json
import math

import numpy as np

from openai_ns_reconstruction.kokuno_outer_base_schedule import KokunoOuterBaseSchedule
from openai_ns_reconstruction.kokuno_radial_modulation_discrepancy import (
    KokunoRadialModulationDiscrepancy,
)


def test_autonomous_modulation_support_is_power_stage_and_precedes_i1():
    modulation = KokunoRadialModulationDiscrepancy()
    report = modulation.source_interval_report()
    power_low, power_high = report["RF40_power_log_interval"]
    support_low, support_high = report["autonomous_modulation_support_log_X"]
    i1_low, _ = report["source_I1_log_interval"]

    assert power_low < support_low < support_high < power_high
    assert support_high < i1_low
    assert report["gap_to_I1_log_units"] > 0.0
    assert report["support_is_autonomous"] is True
    assert report["source_admissible_loop_reconstructed"] is False


def test_autonomous_loop_has_exact_source_mean_contract():
    modulation = KokunoRadialModulationDiscrepancy()
    left, right = modulation.support_log_interval
    y = 0.5 * (left + right)
    phase = (np.arange(4096, dtype=float) + 0.5) / 4096.0
    loop = modulation.loop_shear(
        np.full_like(phase, y),
        np.full_like(phase, 0.2),
        phase=phase,
    )

    assert abs(float(np.mean(loop["a_L"])) - modulation.base_a) < 2.0e-15
    assert abs(float(np.mean(loop["b_L"]))) < 2.0e-15
    assert float(np.max(np.abs(loop["a_L"] - modulation.base_a))) > 0.0
    assert float(np.max(np.abs(loop["b_L"]))) > 0.0


def test_source_finite_n_shear_identities_and_profile_derivatives():
    modulation = KokunoRadialModulationDiscrepancy()
    left, right = modulation.support_log_interval
    y = 0.5 * (left + right) + 0.173
    eta = 0.23
    values = modulation.profile_values_logX(y, eta)

    assert math.isclose(
        float(values["a_N"]), float(values["a_direct"]), rel_tol=2.0e-13, abs_tol=2.0e-13
    )
    assert math.isclose(
        float(values["b_N"]), float(values["b_direct"]), rel_tol=2.0e-13, abs_tol=2.0e-13
    )

    step = 2.0e-7
    plus = modulation.profile_values_logX(y + step, eta)
    minus = modulation.profile_values_logX(y - step, eta)
    fd_e = (float(plus["E"]) - float(minus["E"])) / (2.0 * step)
    fd_u = (float(plus["U"]) - float(minus["U"])) / (2.0 * step)
    assert math.isclose(fd_e, float(values["D_X_E"]), rel_tol=3.0e-6, abs_tol=1.0e-72)
    assert math.isclose(fd_u, float(values["D_X_U"]), rel_tol=3.0e-6, abs_tol=1.0e-72)


def test_modulation_is_exactly_base_outside_compact_support():
    modulation = KokunoRadialModulationDiscrepancy()
    base = KokunoOuterBaseSchedule(outer_schedule=modulation.outer_schedule)
    left, _ = modulation.support_log_interval
    y = left - 0.4
    eta = -0.31

    modified = modulation.profile_values_logX(y, eta)
    reference = base.profile_values_logX(y, eta)
    for key in ("E", "F", "U", "D_X_E", "D_X_F", "D_X_U", "E_eta", "F_eta", "U_eta"):
        assert np.array_equal(np.asarray(modified[key]), np.asarray(reference[key]))
    assert float(modified["A"]) == 0.0
    assert float(modified["B_over_E0"]) == 0.0


def test_five_moment_discrepancy_is_stable_nonzero_and_routes_to_pa17():
    modulation = KokunoRadialModulationDiscrepancy()
    eta = 0.2
    discrepancy = np.asarray(modulation.normalized_discrepancy(eta), dtype=float)
    target = np.asarray(modulation.repair_target_normalized(eta), dtype=float)

    assert discrepancy.shape == (5,)
    assert np.all(np.isfinite(discrepancy))
    assert np.array_equal(target, -discrepancy)
    # This autonomous loop is deliberately nontrivial; C_p is the dominant
    # finite-N moment channel for the default symmetric diagnostic insertion.
    assert 5.0e-5 < abs(float(discrepancy[4])) < 2.0e-4

    refined = KokunoRadialModulationDiscrepancy(
        quadrature_order=12,
        panels_per_period=5,
    )
    refined_discrepancy = np.asarray(refined.normalized_discrepancy(eta), dtype=float)
    assert np.allclose(discrepancy, refined_discrepancy, rtol=2.0e-7, atol=2.0e-11)

    solution = modulation.solve_repair_at_eta(eta)
    assert solution.success
    assert max(abs(value) for value in solution.coefficients) < modulation._repair.coefficient_limit
    assert solution.max_abs_residual < 1.0e-10
    assert np.allclose(
        np.asarray(solution.achieved_normalized),
        target,
        rtol=0.0,
        atol=1.0e-10,
    )


def test_vectorized_discrepancy_and_frequency_refinement():
    etas = np.asarray([-0.4, 0.0, 0.4])
    default = KokunoRadialModulationDiscrepancy(frequency=16)
    values = default.normalized_discrepancy(etas)
    assert values.shape == (3, 5)
    assert np.all(np.isfinite(values))
    assert np.allclose(values[0], values[2], rtol=2.0e-13, atol=2.0e-15)

    coarse = KokunoRadialModulationDiscrepancy(frequency=8)
    fine = KokunoRadialModulationDiscrepancy(frequency=32)
    cp8 = abs(float(coarse.normalized_discrepancy(0.2)[4]))
    cp16 = abs(float(default.normalized_discrepancy(0.2)[4]))
    cp32 = abs(float(fine.normalized_discrepancy(0.2)[4]))
    assert cp8 > cp16 > cp32 > 0.0


def test_serialization_preserves_truth_boundary_and_tamper_fails_closed(tmp_path):
    modulation = KokunoRadialModulationDiscrepancy()
    path = modulation.save_json(tmp_path / "radial_modulation.json")
    loaded = KokunoRadialModulationDiscrepancy.load_json(path)
    assert loaded.sha256 == modulation.sha256
    assert loaded.to_payload() == modulation.to_payload()

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["source_admissible_loop_reconstructed"] is False
    assert payload["truth_boundary"]["actual_source_admissible_loop_discrepancy_supplied"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False

    payload["truth_boundary"]["source_admissible_loop_reconstructed"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        KokunoRadialModulationDiscrepancy.load_json(path)
    except ValueError as exc:
        assert "hash or content mismatch" in str(exc)
    else:
        raise AssertionError("tampered truth metadata must fail closed")
