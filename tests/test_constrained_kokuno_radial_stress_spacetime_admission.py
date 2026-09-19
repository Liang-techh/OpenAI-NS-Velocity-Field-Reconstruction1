from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.kokuno_radial_stress_spacetime_admission import (
    admit_spacetime_radial_stress_report,
)


def _channel(exponent: int, error: float, weighted_moment: float) -> dict[str, float | int]:
    return {
        "exponent": exponent,
        "weighted_moment": weighted_moment,
        "operator_relative_rms": error,
        "operator_relative_max": 1.5 * error,
        "moment_complement_relative": 1.0e-16,
        "edge_relative": 1.0e-15,
        "sign_flip_mutation_relative_rms": 1.999,
        "rhs_rms": 10.0,
        "stress_rms": 2.0,
        "fd8_interior_count": 17,
    }


def _state(index: int, t: float, z: float, scale: float) -> dict[str, object]:
    theta_errors = [1.8e-2 * scale, 1.5e-3 * scale, 1.0e-4 * scale]
    axial_errors = [3.8e-2 * scale, 3.0e-3 * scale, 2.0e-4 * scale]
    levels = []
    for count, theta_error, axial_error in zip((25, 49, 97), theta_errors, axial_errors):
        levels.append(
            {
                "radial_count": count,
                "radial_spacing": 1.0 / count,
                "boundary_margin": 0.0125,
                "stencil_reach": 0.01,
                "velocity_rms": 100.0,
                "self_defect_rms": 1.0e6,
                "theta_e2": _channel(2, theta_error, 3.0),
                "axial_e1": _channel(1, axial_error, -2.0),
            }
        )

    def convergence(name: str) -> dict[str, object]:
        rows = [level[name] for level in levels]
        errors = [float(row["operator_relative_rms"]) for row in rows]
        return {
            "relative_rms_by_radial_count": errors,
            "refinement_ratios": [errors[0] / errors[1], errors[1] / errors[2]],
            "finest_relative_max": float(rows[-1]["operator_relative_max"]),
            "finest_moment_complement_relative": float(rows[-1]["moment_complement_relative"]),
            "finest_edge_relative": float(rows[-1]["edge_relative"]),
            "finest_sign_flip_mutation_relative_rms": float(rows[-1]["sign_flip_mutation_relative_rms"]),
            "finest_weighted_moment": float(rows[-1]["weighted_moment"]),
        }

    return {
        "state_index": index,
        "t": t,
        "z": z,
        "levels": levels,
        "convergence": {
            "theta_e2": convergence("theta_e2"),
            "axial_e1": convergence("axial_e1"),
        },
        "failed_guards": [],
        "state_passed": True,
    }


def _fixture_report() -> dict[str, object]:
    return {
        "schema": "kokuno-a4-radial-stress-spacetime-generalization-v2",
        "task": "KOKUNO-A4-RADIAL-STRESS-SPACETIME-GENERALIZATION-044",
        "provenance": {
            "parent_agent3_pr": 598,
            "parent_agent3_head": "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233",
            "admitted_agent2_pr": 561,
            "admitted_agent2_head": "732800ce4990464b49c8aa32d0dff4580f6684d4",
            "prior_agent4_pr": 595,
            "prior_agent4_head": "19296acad4f84ba05d3b095bb8c01bf5a2c95891",
            "prior_anchor_state": {"t": 0.50, "z": 0.08},
        },
        "protocol_repair": {
            "first_attempt_workflow": 35433365648,
            "first_attempt_scientific_values_emitted": False,
            "first_attempt_failure": "129-node radial spacing .009375 < FD8 four-step reach .01",
            "repair": "use 25/49/97 radial ladder and 16 physical angles; preserve FD8 step, states, metrics and all scientific guards",
        },
        "protocol": {
            "radial_counts": [25, 49, 97],
            "angular_count": 16,
            "new_spacetime_states": [{"t": 0.37, "z": -0.31}, {"t": 0.63, "z": 0.31}],
            "prior_anchor_state": {"t": 0.50, "z": 0.08},
            "cartesian_derivative": "centered FD8, independently implemented",
            "cartesian_fd8_step": 0.0025,
            "radial_derivative": "separate centered FD8",
            "nu": 0.01,
            "guards": {
                "finest_relative_rms_max": 2.0e-2,
                "finest_relative_max_max": 5.0e-2,
                "minimum_refinement_ratio": 2.0,
                "moment_complement_relative_max": 1.0e-10,
                "edge_relative_max": 1.0e-8,
                "sign_flip_mutation_min": 5.0e-1,
                "nontrivial_self_defect_rms_min": 1.0e-8,
            },
        },
        "states": [_state(0, 0.37, -0.31, 1.0), _state(1, 0.63, 0.31, 0.8)],
        "failed_guards": [],
        "spacetime_radial_stress_generalization_passed": True,
        "truth_boundary": {
            "public_black_box_velocity_consumed": True,
            "public_black_box_time_derivative_consumed": True,
            "agent3_self_defect_operator_used": False,
            "agent3_compact_stress_constructor_used": False,
            "agent3_radial_derivative_used": False,
            "pressure_fitted": False,
            "forcing_fitted": False,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "public_velocity_correction_materialized": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": 1.0e-3,
            "final_normalized_divergence_gate": 1.0e-5,
        },
    }


def test_admits_only_spacetime_radial_operator_reuse() -> None:
    receipt = admit_spacetime_radial_stress_report(_fixture_report())
    assert receipt["independent_spacetime_radial_stress_generalization_admitted"] is True
    assert receipt["finite_correction_cycle_radial_operator_reuse_allowed_without_retuning"] is True
    assert receipt["full_composite_radial_stress_execution_allowed_when_actual_defect_available"] is True
    assert receipt["full_same_cycle_composite_requested_stress_materialized"] is False
    assert receipt["signed_mean_inverse_input_ready"] is False
    assert receipt["finite_correction_cycle_run"] is False
    assert receipt["heldout_ns_momentum_residual_assessed"] is False
    assert receipt["pde_validated"] is False


def test_rejects_threshold_laundering() -> None:
    report = deepcopy(_fixture_report())
    report["protocol"]["guards"]["finest_relative_rms_max"] = 0.2
    with pytest.raises(ValueError, match="frozen scientific guard changed"):
        admit_spacetime_radial_stress_report(report)


def test_rejects_scientific_failure_even_if_pass_bit_stays_true() -> None:
    report = deepcopy(_fixture_report())
    fine = report["states"][0]["levels"][-1]["theta_e2"]
    fine["operator_relative_rms"] = 0.03
    report["states"][0]["convergence"]["theta_e2"]["relative_rms_by_radial_count"][-1] = 0.03
    report["states"][0]["convergence"]["theta_e2"]["refinement_ratios"][-1] = (
        report["states"][0]["convergence"]["theta_e2"]["relative_rms_by_radial_count"][-2] / 0.03
    )
    with pytest.raises(ValueError, match="declared failed_guards disagrees"):
        admit_spacetime_radial_stress_report(report)


def test_rejects_truth_boundary_promotion() -> None:
    report = deepcopy(_fixture_report())
    report["truth_boundary"]["public_velocity_correction_materialized"] = True
    with pytest.raises(ValueError, match="cannot promote public_velocity_correction_materialized"):
        admit_spacetime_radial_stress_report(report)


def test_rejects_wrong_upstream_identity() -> None:
    report = deepcopy(_fixture_report())
    report["provenance"]["admitted_agent2_head"] = "0" * 40
    with pytest.raises(ValueError, match="Agent-4 provenance changed"):
        admit_spacetime_radial_stress_report(report)
