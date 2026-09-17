from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_heat_discrepancy_repair import (
    KokunoHeatDiscrepancyRepair,
)
from openai_ns_reconstruction.kokuno_independent_heat_repair import (
    _independent_increment,
    run_audit,
)


def test_independent_heat_repair_structural_preflight_passes() -> None:
    report = run_audit(levels=(513, 1025, 2049), sample_count=4)
    summary = report["summary"]
    assert summary["structural_preflight_passed"] is True
    assert summary["max_finest_public_vs_independent_relative_error"] < 2.0e-5
    assert summary["max_analytic_vs_independent_jacobian_relative_error"] < 2.0e-4
    assert summary["finest_profile_derivative_relative_error"] < 2.0e-4
    assert summary["min_parameter_mutation_abs_detection"] > 1.0e-6
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["fixed_references"]["normalized_ns_residual"] == 1.0e-3
    assert report["fixed_references"]["divergence"] == 1.0e-5


def test_independent_simpson_refines_toward_public_map() -> None:
    repair = KokunoHeatDiscrepancyRepair()
    coefficients = np.asarray([0.003, -0.002, 0.0015])
    f_eta = float(repair.source_f(0.2))
    public = np.asarray(repair.normalized_increment(coefficients, f_eta=f_eta))
    coarse = _independent_increment(repair, coefficients, f_eta, 257)
    medium = _independent_increment(repair, coefficients, f_eta, 513)
    fine = _independent_increment(repair, coefficients, f_eta, 1025)
    coarse_error = np.linalg.norm(coarse - public)
    medium_error = np.linalg.norm(medium - public)
    fine_error = np.linalg.norm(fine - public)
    assert fine_error < medium_error < coarse_error


def test_parameter_mutation_is_detected_without_gate_change() -> None:
    repair = KokunoHeatDiscrepancyRepair()
    coefficients = np.asarray([0.003, -0.002, 0.0015])
    f_eta = float(repair.source_f(-0.31))
    reference = _independent_increment(repair, coefficients, f_eta, 2049)
    mutated = coefficients.copy()
    mutated[0] += 0.002
    mutated_public = np.asarray(repair.normalized_increment(mutated, f_eta=f_eta))
    assert np.linalg.norm(mutated_public - reference) > 1.0e-6
