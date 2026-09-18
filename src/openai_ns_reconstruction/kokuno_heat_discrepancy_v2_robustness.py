"""Independent robustness audit for the corrected Kokuno heat discrepancy.

Agent 1 PR #289 fixes the first-order terminal ``Delta C_p`` factor found by
Agent 4 PR #283.  The production implementation remains Gauss-Legendre /
Gauss-Laguerre based.  This audit deliberately keeps the independent PR #283
small-z series plus adaptive log-space quadrature path and asks a new question:
does the correction remain correct under parameter perturbations and a
three-level production-quadrature ladder, rather than only at the default
checkpoint?

This is a local source-moment validation.  It is not a velocity residual test,
is not directly comparable with the retained ST006 full-domain momentum
benchmark, and cannot set ``pde_validated=true``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_heat_replacement_discrepancy import (
    KokunoHeatReplacementDiscrepancy,
)
from .kokuno_heat_replacement_discrepancy_independent import (
    IndependentHeatDiscrepancyAudit,
)
from .kokuno_heat_replacement_discrepancy_v2 import (
    KokunoHeatReplacementDiscrepancyV2,
)


TASK_ID = "KOKUNO-A4-HEAT-V2-ROBUSTNESS-INDEPENDENT-010"
SEED = 9173061  # provenance only; the fixed grid below is deterministic.
ETAS = np.asarray([-0.8, -0.2, 0.0, 0.2, 0.8], dtype=float)

PARAMETER_CASES = (
    {
        "name": "lower_scale",
        "h": 0.0035,
        "X_tail": 400.0,
        "c_inf": 0.8,
        "rho_o": 0.00015,
    },
    {
        "name": "baseline",
        "h": 0.005,
        "X_tail": 1000.0,
        "c_inf": 1.0,
        "rho_o": 0.00025,
    },
    {
        "name": "upper_scale",
        "h": 0.0075,
        "X_tail": 2500.0,
        "c_inf": 1.2,
        "rho_o": 0.00035,
    },
)

QUADRATURE_LEVELS = (
    {
        "label": "coarse",
        "heat_quadrature_order": 48,
        "transition_quadrature_order": 96,
        "tail_quadrature_order": 48,
    },
    {
        "label": "medium",
        "heat_quadrature_order": 96,
        "transition_quadrature_order": 192,
        "tail_quadrature_order": 96,
    },
    {
        "label": "fine",
        "heat_quadrature_order": 192,
        "transition_quadrature_order": 384,
        "tail_quadrature_order": 192,
    },
)

# Local implementation-consistency guards registered before executing this audit.
# They are not PDE acceptance thresholds.
FINE_COMPONENT_RELATIVE_GUARDS = np.asarray([2.0e-5, 2.0e-6, 2.0e-5])
MEDIUM_TO_FINE_RELATIVE_GUARD = 5.0e-5
TAIL_FIX_RELATIVE_GUARD = 5.0e-5
SYMMETRY_RELATIVE_GUARD = 5.0e-12
MUTATION_CP_RELATIVE_FLOOR = 1.0e-4


def _model_kwargs(case: dict[str, Any], level: dict[str, Any]) -> dict[str, Any]:
    return {
        "h": float(case["h"]),
        "X_tail": float(case["X_tail"]),
        "c_inf": float(case["c_inf"]),
        "rho_o": float(case["rho_o"]),
        "heat_quadrature_order": int(level["heat_quadrature_order"]),
        "transition_quadrature_order": int(level["transition_quadrature_order"]),
        "tail_quadrature_order": int(level["tail_quadrature_order"]),
    }


def _relative_error(value: np.ndarray, reference: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    reference = np.asarray(reference, dtype=float)
    return np.abs(value - reference) / np.maximum(np.abs(reference), 1.0e-300)


def _max_relative_between(value: np.ndarray, reference: np.ndarray) -> float:
    return float(np.max(_relative_error(value, reference)))


def run_audit() -> dict[str, Any]:
    audit = IndependentHeatDiscrepancyAudit(series_order=16, epsrel=1.0e-10)
    case_reports: list[dict[str, Any]] = []

    finest_component_rel = np.zeros(3, dtype=float)
    medium_to_fine_max = 0.0
    tail_fix_rel_max = 0.0
    channel_locality_abs_max = 0.0
    endpoint_abs_max = 0.0
    symmetry_rel_max = 0.0
    mutation_cp_relative_min = float("inf")

    for case in PARAMETER_CASES:
        fine_level = QUADRATURE_LEVELS[-1]
        reference_model = KokunoHeatReplacementDiscrepancyV2(
            **_model_kwargs(case, fine_level)
        )
        independent = np.stack(
            [audit.discrepancy(reference_model, float(eta)) for eta in ETAS], axis=0
        )

        level_reports: list[dict[str, Any]] = []
        public_by_label: dict[str, np.ndarray] = {}
        for level in QUADRATURE_LEVELS:
            model = KokunoHeatReplacementDiscrepancyV2(**_model_kwargs(case, level))
            public = np.asarray(model.discrepancy(ETAS), dtype=float)
            public_by_label[str(level["label"])] = public
            rel = _relative_error(public, independent)
            level_reports.append(
                {
                    "label": level["label"],
                    "orders": {
                        "heat": level["heat_quadrature_order"],
                        "transition": level["transition_quadrature_order"],
                        "tail": level["tail_quadrature_order"],
                    },
                    "max_relative_by_component": np.max(rel, axis=0).tolist(),
                    "max_relative_all": float(np.max(rel)),
                }
            )

        fine = public_by_label["fine"]
        medium = public_by_label["medium"]
        fine_rel = np.max(_relative_error(fine, independent), axis=0)
        finest_component_rel = np.maximum(finest_component_rel, fine_rel)
        case_medium_to_fine = _max_relative_between(medium, fine)
        medium_to_fine_max = max(medium_to_fine_max, case_medium_to_fine)

        old = KokunoHeatReplacementDiscrepancy(**_model_kwargs(case, fine_level))
        before = np.asarray(old.discrepancy(ETAS), dtype=float)
        change = fine - before
        independent_first_order = np.asarray(
            [audit.first_order_tail_cp(reference_model, float(eta)) for eta in ETAS],
            dtype=float,
        )
        expected_added_copy = 0.5 * independent_first_order
        active = np.abs(expected_added_copy) > 1.0e-300
        case_tail_rel = float(
            np.max(
                np.abs(change[active, 0] - expected_added_copy[active])
                / np.abs(expected_added_copy[active])
            )
        )
        tail_fix_rel_max = max(tail_fix_rel_max, case_tail_rel)
        case_channel_locality = float(np.max(np.abs(change[:, 1:])))
        channel_locality_abs_max = max(channel_locality_abs_max, case_channel_locality)

        endpoints = np.asarray(reference_model.discrepancy(np.asarray([-1.0, 1.0])))
        case_endpoint_abs = float(np.max(np.abs(endpoints)))
        endpoint_abs_max = max(endpoint_abs_max, case_endpoint_abs)

        plus = np.asarray(reference_model.discrepancy(np.asarray([0.2, 0.8])))
        minus = np.asarray(reference_model.discrepancy(np.asarray([-0.2, -0.8])))
        case_symmetry = _max_relative_between(plus, minus)
        symmetry_rel_max = max(symmetry_rel_max, case_symmetry)

        # Mutation calibration: remove exactly the newly restored copy from the
        # corrected public result.  This recreates the audited factor-of-two
        # defect without changing any validation reference.
        mutated = fine.copy()
        mutated[:, 0] -= change[:, 0]
        mutation_rel = _relative_error(mutated[:, 0], independent[:, 0])
        active_reference = np.abs(independent[:, 0]) > 1.0e-300
        case_mutation_min = float(np.min(mutation_rel[active_reference]))
        mutation_cp_relative_min = min(mutation_cp_relative_min, case_mutation_min)

        case_reports.append(
            {
                "name": case["name"],
                "parameters": {
                    key: float(case[key]) for key in ("h", "X_tail", "c_inf", "rho_o")
                },
                "etas": ETAS.tolist(),
                "independent_reference": independent.tolist(),
                "levels": level_reports,
                "medium_to_fine_max_relative": case_medium_to_fine,
                "tail_fix_vs_independent_half_first_order_max_relative": case_tail_rel,
                "unchanged_delta_s_i_max_abs": case_channel_locality,
                "endpoint_eta_pm1_max_abs": case_endpoint_abs,
                "eta_even_symmetry_max_relative": case_symmetry,
                "missing-copy_mutation_min_cp_relative_error": case_mutation_min,
            }
        )

    checks = {
        "fine_component_relative_guard": bool(
            np.all(finest_component_rel <= FINE_COMPONENT_RELATIVE_GUARDS)
        ),
        "medium_to_fine_relative_guard": bool(
            medium_to_fine_max <= MEDIUM_TO_FINE_RELATIVE_GUARD
        ),
        "tail_fix_relative_guard": bool(tail_fix_rel_max <= TAIL_FIX_RELATIVE_GUARD),
        "unchanged_channels_exact": bool(channel_locality_abs_max == 0.0),
        "eta_endpoints_exact_zero": bool(endpoint_abs_max == 0.0),
        "eta_even_symmetry_guard": bool(symmetry_rel_max <= SYMMETRY_RELATIVE_GUARD),
        "missing_copy_mutation_detected": bool(
            mutation_cp_relative_min >= MUTATION_CP_RELATIVE_FLOOR
        ),
    }
    passed = bool(all(checks.values()))

    return {
        "task_id": TASK_ID,
        "seed": SEED,
        "operator_independence": {
            "production": "Gauss-Legendre transition + Gauss-Laguerre terminal tail",
            "validation": (
                "small-z heat series + scipy adaptive quadrature in log X and "
                "s=h(y-3) tail coordinate"
            ),
            "training_loss_used": False,
            "candidate_internal_training_tensor_used": False,
        },
        "parameter_cases": case_reports,
        "summary": {
            "finest_max_relative_by_component": finest_component_rel.tolist(),
            "medium_to_fine_max_relative": medium_to_fine_max,
            "tail_fix_vs_independent_half_first_order_max_relative": tail_fix_rel_max,
            "unchanged_delta_s_i_max_abs": channel_locality_abs_max,
            "endpoint_eta_pm1_max_abs": endpoint_abs_max,
            "eta_even_symmetry_max_relative": symmetry_rel_max,
            "missing_copy_mutation_min_cp_relative_error": mutation_cp_relative_min,
        },
        "registered_local_guards": {
            "fine_component_relative": FINE_COMPONENT_RELATIVE_GUARDS.tolist(),
            "medium_to_fine_relative": MEDIUM_TO_FINE_RELATIVE_GUARD,
            "tail_fix_relative": TAIL_FIX_RELATIVE_GUARD,
            "eta_even_symmetry_relative": SYMMETRY_RELATIVE_GUARD,
            "missing_copy_mutation_cp_relative_floor": MUTATION_CP_RELATIVE_FLOOR,
        },
        "checks": checks,
        "structural_preflight_passed": passed,
        "st006_directly_comparable": False,
        "st006_comparability_reason": (
            "this audit validates a local source heat-moment map, not the frozen "
            "4096-point full-domain normalized momentum protocol"
        ),
        "fixed_project_references": {
            "normalized_full_momentum": 1.0e-3,
            "divergence_max": 1.0e-5,
            "changed": False,
        },
        "formal_full_domain_pde_gate_assessed": False,
        "normalized_ns_residual_le_1e-3_claimed": False,
        "pde_validated": False,
    }


def write_report(path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    report = run_audit()
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return destination


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/heat_v2_robustness_report.json",
    )
    args = parser.parse_args()
    destination = write_report(args.output)
    report = json.loads(destination.read_text(encoding="utf-8"))
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    print(f"structural_preflight_passed={report['structural_preflight_passed']}")


if __name__ == "__main__":
    main()
