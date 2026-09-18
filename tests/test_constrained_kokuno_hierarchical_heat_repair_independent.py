from __future__ import annotations

import math

import numpy as np

from openai_ns_reconstruction.kokuno_hierarchical_heat_repair import (
    KokunoHierarchicalHeatRepair,
)
from openai_ns_reconstruction.kokuno_hierarchical_heat_repair_independent import (
    ETA_PROBES,
    RESOLUTIONS,
    build_report,
    independent_moment_tensors,
)


def test_independent_simpson_moment_tensors_refine_without_construction_map() -> None:
    repair = KokunoHierarchicalHeatRepair()
    lam = repair.outer_schedule.lambda_outer
    f_eta = 1.0 / (1.0 + 0.2**2)
    medium = independent_moment_tensors(
        lambda_outer=lam, f_eta=f_eta, resolution=RESOLUTIONS[-2]
    )
    fine = independent_moment_tensors(
        lambda_outer=lam, f_eta=f_eta, resolution=RESOLUTIONS[-1]
    )

    assert medium.linear.shape == (3, 3)
    assert medium.quadratic.shape == (3, 3, 3)
    assert np.all(np.isfinite(medium.linear))
    assert np.all(np.isfinite(medium.quadratic))
    assert np.max(np.abs(fine.linear - medium.linear)) <= 1.0e-10
    assert np.max(np.abs(fine.quadratic - medium.quadratic)) <= 1.0e-10
    # Disjoint support is an independently encoded structural property.
    for row in range(3):
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert fine.quadratic[row, i, j] == 0.0


def test_independent_audit_exposes_float64_moment_precision_barrier() -> None:
    report = build_report()
    summary = report["summary"]
    guards = report["local_guards"]

    # The preregistered independent-I_sub closure guard actually fails.  That
    # scientific rejection is the result to preserve, not a reason to weaken
    # the tolerance after seeing the data.
    assert report["local_audit_completed"] is False
    assert guards["independent_tensor_refinement_max_abs_le_1e-10"] is True
    assert guards["dominant_I_sub_channel_relative_error_le_1e-10"] is False
    assert guards["Cp_independent_mismatch_exceeds_target_by_50_decades"] is True
    assert guards["S_independent_mismatch_exceeds_target_by_50_decades"] is True
    assert summary["max_tensor_refinement_abs"] <= 1.0e-10
    assert summary["max_I_sub_relative_error"] > 1.0e-10
    # The tiny C_p/S channels miss their targets by many decades when the same
    # public coefficient receipt is evaluated through an independent moment
    # tensor.  Thus frozen-discrete closure is not continuous-source closure.
    assert summary["minimum_Cp_log10_relative_error"] >= 50.0
    assert summary["minimum_S_log10_relative_error"] >= 50.0
    assert summary["independent_continuous_moment_closure"] is False
    assert summary["frozen_discrete_map_only"] is True

    assert [case["eta"] for case in report["cases"]] == list(ETA_PROBES)
    assert all(
        len(case["resolution_evaluations"]) == len(RESOLUTIONS)
        for case in report["cases"]
    )


def test_mutation_is_detected_and_truth_boundary_stays_fail_closed() -> None:
    report = build_report()
    assert report["summary"]["minimum_mutated_I_sub_relative_error"] >= 1.0e-6
    assert report["local_guards"]["coefficient_mutation_detected_in_I_sub_ge_1e-6"] is True
    assert report["formal_project_gates"]["normalized_momentum_max"] == 1.0e-3
    assert report["formal_project_gates"]["divergence_max"] == 1.0e-5
    assert report["formal_project_gates"]["assessed_here"] is False
    assert report["formal_project_gates"]["pde_validated"] is False
    assert report["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["truth_boundary"]["continuous_source_moment_compensation_certified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert math.isfinite(report["summary"]["minimum_Cp_log10_relative_error"])
