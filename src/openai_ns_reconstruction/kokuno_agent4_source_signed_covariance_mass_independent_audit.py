"""Independent quadrature audit of the signed pulse covariance-mass handoff.

This module is deliberately downstream of Agent 2's
``KokunoSourceSignedCovarianceMass`` public interface.  It supplies labelled
nonlinear candidate pulse/cutoff samples and compares only the returned public
values against a separate Gauss--Legendre oracle.  It does not call Agent 2's
trapezoid helper, band-covering helper, signed covariance inverse, complete-curl
implementation, training tensors/loss, pressure, or forcing.

A PASS is local numerical evidence for the displayed candidate pulse-mass law
only.  The source does not publish the actual pulse/background/mode realization,
and no public global oscillatory velocity or full Navier--Stokes candidate exists
on this ancestry.  Therefore the formal PDE gate remains unassessed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)
from .kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)
from .kokuno_source_signed_covariance_mass import KokunoSourceSignedCovarianceMass

SEED = 9173221
SOURCE_B_G = math.sqrt(2.0) - 1.0
H_CASES = (0.0025, 0.0065, 0.0090)
LADDER = (33, 65, 129)
FIXED_FINE_SAMPLES = 2049
GL_ORDER = 256

FROZEN_GUARDS = {
    "pulse_finest_relative_rms": 1.0e-5,
    "pulse_refinement_ratio_min": 3.5,
    "transverse_finest_relative": 1.0e-10,
    "transverse_refinement_ratio_min": 20.0,
    "combined_finest_relative_rms": 1.0e-5,
    "combined_finest_relative_max": 2.0e-5,
    "wrong_psi_power_mutation_relative_rms_min": 0.25,
    "wrong_transverse_prefactor_mutation_relative_min": 0.25,
}


def _relative_rms(actual: np.ndarray, expected: np.ndarray) -> float:
    numerator = float(np.sqrt(np.mean((actual - expected) ** 2)))
    denominator = max(1.0e-300, float(np.sqrt(np.mean(expected**2))))
    return numerator / denominator


def _relative_max(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = np.maximum(1.0e-300, np.abs(expected))
    return float(np.max(np.abs(actual - expected) / scale))


def _profile_phase(case_index: int, beta_index: int) -> float:
    # Deterministic but non-commensurate enough to avoid one repeated profile.
    return 0.03 + 0.047 * case_index + 0.071 * beta_index


def _pulse_profiles(s: np.ndarray, phase: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Smooth positive candidate profiles on normalized pulse coordinate s in [0,1]."""
    psi = (
        0.58
        + 0.20 * np.cos(2.0 * math.pi * (s + phase))
        + 0.08 * np.sin(math.pi * s) ** 2
    )
    x_plus = 1.15 + 0.25 * s + 0.12 * np.sin(2.0 * math.pi * (s + phase))
    x_minus = 1.55 + 0.18 * (1.0 - s) + 0.10 * np.cos(3.0 * math.pi * s - phase)
    return psi, x_plus, x_minus


def _chi_profile(xi: np.ndarray) -> np.ndarray:
    # Compact-looking smooth numerical cutoff profile on the registered [-1,1] sample interval.
    return (1.0 - xi * xi) ** 2 * (
        0.80 + 0.10 * xi + 0.05 * np.cos(math.pi * xi)
    )


def _geometry(h: float) -> dict[str, Any]:
    # This geometry seam was independently audited in Agent-4 PR #504.  Here it
    # is only used to supply the public L_s/r0/beta schedule expected by #513.
    partition = KokunoSourceCompatiblePartitionRealization(ell_min=5).evaluate(
        q=np.asarray(2.0**-5.5),
        D_r_q=np.asarray(0.0),
        D_z_q=np.asarray(0.0),
        slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
    )
    return KokunoAutonomousSignedRectangleGeometry(h=h).instantiate(partition)


def _sample_pulses(
    geometry: dict[str, Any], case_index: int, sample_count: int
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    v_rows: list[np.ndarray] = []
    psi_rows: list[np.ndarray] = []
    x_rows: list[np.ndarray] = []
    s = np.linspace(0.0, 1.0, sample_count)
    for beta_index, L_s in enumerate(geometry["L_s_by_beta"]):
        phase = _profile_phase(case_index, beta_index)
        psi, x_plus, x_minus = _pulse_profiles(s, phase)
        v_rows.append(float(L_s) * s)
        psi_rows.append(psi)
        x_rows.append(np.stack((x_plus, x_minus), axis=-1))
    return tuple(v_rows), tuple(psi_rows), tuple(x_rows)


def _independent_reference(
    geometry: dict[str, Any], case_index: int
) -> dict[str, np.ndarray | float]:
    # Independent high-order quadrature path.  No Agent-2 integration helper is used.
    nodes, weights = leggauss(GL_ORDER)
    s = 0.5 * (nodes + 1.0)
    s_weights = 0.5 * weights

    normalized_i0: list[tuple[float, float]] = []
    wrong_normalized_i0: list[tuple[float, float]] = []
    for beta_index, _ in enumerate(geometry["beta_labels"]):
        phase = _profile_phase(case_index, beta_index)
        psi, x_plus, x_minus = _pulse_profiles(s, phase)
        normalized_i0.append(
            (
                float(np.sum(s_weights * psi * psi * x_plus * x_plus)),
                float(np.sum(s_weights * psi * psi * x_minus * x_minus)),
            )
        )
        # Calibrated mutation: incorrectly use psi instead of psi^2.
        wrong_normalized_i0.append(
            (
                float(np.sum(s_weights * psi * x_plus * x_plus)),
                float(np.sum(s_weights * psi * x_minus * x_minus)),
            )
        )

    chi_squared_integral = float(np.sum(weights * _chi_profile(nodes) ** 2))
    D_g = 0.5 * (1.0 + SOURCE_B_G * SOURCE_B_G) * chi_squared_integral
    r0 = float(geometry["r0"])
    normalized = np.asarray(normalized_i0, dtype=float)
    wrong_normalized = np.asarray(wrong_normalized_i0, dtype=float)
    # Since L_s=2*r0/c_i, h_sigma=D_g*c_i*L_s*Ibar=2*r0*D_g*Ibar.
    h_sigma = 2.0 * r0 * D_g * normalized
    wrong_h_sigma_psi_power = 2.0 * r0 * D_g * wrong_normalized
    wrong_D_g = chi_squared_integral  # mutation: drop source (1+b_g^2)/2 factor
    wrong_h_sigma_prefactor = 2.0 * r0 * wrong_D_g * normalized
    return {
        "normalized_i0": normalized,
        "D_g": D_g,
        "chi_squared_integral": chi_squared_integral,
        "h_sigma": h_sigma,
        "wrong_h_sigma_psi_power": wrong_h_sigma_psi_power,
        "wrong_h_sigma_prefactor": wrong_h_sigma_prefactor,
    }


def _materialize(
    evaluator: KokunoSourceSignedCovarianceMass,
    geometry: dict[str, Any],
    case_index: int,
    pulse_samples: int,
    transverse_samples: int,
) -> dict[str, Any]:
    v, psi, x = _sample_pulses(geometry, case_index, pulse_samples)
    xi = np.linspace(-1.0, 1.0, transverse_samples)
    chi = _chi_profile(xi)
    return evaluator.materialize_from_autonomous_geometry(
        geometry,
        v_by_beta=v,
        psi_by_beta=psi,
        x_by_beta_sign=x,
        xi=xi,
        chi_g=chi,
        pulse_binding="repository_autonomous_source_compatible",
    )


def _refinement_ratios(errors: list[float]) -> list[float]:
    ratios: list[float] = []
    for coarse, fine in zip(errors[:-1], errors[1:]):
        ratios.append(float(coarse / max(1.0e-300, fine)))
    return ratios


def run_audit() -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    # The seed perturbs only harmless case ordering, not thresholds or profile formulas.
    case_order = tuple(int(x) for x in rng.permutation(len(H_CASES)))

    pulse_errors_by_level: list[list[float]] = [[] for _ in LADDER]
    transverse_errors_by_level: list[list[float]] = [[] for _ in LADDER]
    combined_actual_by_level: list[list[np.ndarray]] = [[] for _ in LADDER]
    combined_expected_by_level: list[list[np.ndarray]] = [[] for _ in LADDER]
    wrong_psi_actual: list[np.ndarray] = []
    wrong_prefactor_actual: list[np.ndarray] = []
    correct_reference: list[np.ndarray] = []
    case_receipts: list[dict[str, Any]] = []
    truth_flags_preserved = True

    for case_index in case_order:
        h = H_CASES[case_index]
        geometry = _geometry(h)
        evaluator = KokunoSourceSignedCovarianceMass(h=h)
        reference = _independent_reference(geometry, case_index)
        normalized_i0_ref = np.asarray(reference["normalized_i0"], dtype=float)
        h_ref = np.asarray(reference["h_sigma"], dtype=float)
        D_ref = float(reference["D_g"])

        pulse_case_errors: list[float] = []
        transverse_case_errors: list[float] = []
        combined_case_errors: list[float] = []
        for level_index, sample_count in enumerate(LADDER):
            pulse_out = _materialize(
                evaluator,
                geometry,
                case_index,
                sample_count,
                FIXED_FINE_SAMPLES,
            )
            normalized_i0_public = (
                np.asarray(pulse_out["I0_by_beta_sign"], dtype=float)
                / np.asarray(geometry["L_s_by_beta"], dtype=float)[:, None]
            )
            pulse_error = _relative_rms(normalized_i0_public, normalized_i0_ref)
            pulse_errors_by_level[level_index].append(pulse_error)
            pulse_case_errors.append(pulse_error)

            transverse_out = _materialize(
                evaluator,
                geometry,
                case_index,
                FIXED_FINE_SAMPLES,
                sample_count,
            )
            transverse_error = abs(float(transverse_out["D_g"]) - D_ref) / max(
                1.0e-300, abs(D_ref)
            )
            transverse_errors_by_level[level_index].append(transverse_error)
            transverse_case_errors.append(transverse_error)

            combined_out = _materialize(
                evaluator,
                geometry,
                case_index,
                sample_count,
                sample_count,
            )
            h_public = np.asarray(combined_out["h_sigma_by_beta_sign"], dtype=float)
            combined_actual_by_level[level_index].append(h_public)
            combined_expected_by_level[level_index].append(h_ref)
            combined_case_errors.append(_relative_rms(h_public, h_ref))

        wrong_psi_actual.append(np.asarray(reference["wrong_h_sigma_psi_power"], dtype=float))
        wrong_prefactor_actual.append(np.asarray(reference["wrong_h_sigma_prefactor"], dtype=float))
        correct_reference.append(h_ref)
        case_receipts.append(
            {
                "case_index": case_index,
                "h": h,
                "beta_count": len(geometry["beta_labels"]),
                "r0": float(geometry["r0"]),
                "pulse_relative_rms": pulse_case_errors,
                "transverse_relative": transverse_case_errors,
                "combined_h_sigma_relative_rms": combined_case_errors,
            }
        )

        receipt_truth = evaluator.receipt()["truth_boundary"]
        truth_flags_preserved = truth_flags_preserved and (
            receipt_truth["actual_source_h_sigma_pulse_integrals_bound"] is False
            and receipt_truth["actual_source_pulse_samples_recovered"] is False
            and receipt_truth["actual_positive_order_background_bound"] is False
            and receipt_truth["actual_auxiliary_torus_mode_family_bound"] is False
            and receipt_truth["public_xyz_t_velocity_correction_materialized"] is False
            and receipt_truth["formal_full_domain_pde_gate_assessed"] is False
            and receipt_truth["pde_validated"] is False
        )

    pulse_ladder = [
        float(np.sqrt(np.mean(np.asarray(level, dtype=float) ** 2)))
        for level in pulse_errors_by_level
    ]
    transverse_ladder = [
        float(np.sqrt(np.mean(np.asarray(level, dtype=float) ** 2)))
        for level in transverse_errors_by_level
    ]
    combined_ladder: list[float] = []
    combined_max_ladder: list[float] = []
    for actual_rows, expected_rows in zip(combined_actual_by_level, combined_expected_by_level):
        actual = np.concatenate([row.reshape(-1) for row in actual_rows])
        expected = np.concatenate([row.reshape(-1) for row in expected_rows])
        combined_ladder.append(_relative_rms(actual, expected))
        combined_max_ladder.append(_relative_max(actual, expected))

    pulse_refinement = _refinement_ratios(pulse_ladder)
    transverse_refinement = _refinement_ratios(transverse_ladder)
    combined_refinement = _refinement_ratios(combined_ladder)

    expected_all = np.concatenate([row.reshape(-1) for row in correct_reference])
    wrong_psi_all = np.concatenate([row.reshape(-1) for row in wrong_psi_actual])
    wrong_prefactor_all = np.concatenate([row.reshape(-1) for row in wrong_prefactor_actual])
    wrong_psi_mutation = _relative_rms(wrong_psi_all, expected_all)
    wrong_prefactor_mutation = _relative_rms(wrong_prefactor_all, expected_all)

    local_pass = (
        pulse_ladder[-1] <= FROZEN_GUARDS["pulse_finest_relative_rms"]
        and min(pulse_refinement) >= FROZEN_GUARDS["pulse_refinement_ratio_min"]
        and transverse_ladder[-1] <= FROZEN_GUARDS["transverse_finest_relative"]
        and min(transverse_refinement) >= FROZEN_GUARDS["transverse_refinement_ratio_min"]
        and combined_ladder[-1] <= FROZEN_GUARDS["combined_finest_relative_rms"]
        and combined_max_ladder[-1] <= FROZEN_GUARDS["combined_finest_relative_max"]
        and wrong_psi_mutation
        >= FROZEN_GUARDS["wrong_psi_power_mutation_relative_rms_min"]
        and wrong_prefactor_mutation
        >= FROZEN_GUARDS["wrong_transverse_prefactor_mutation_relative_min"]
        and truth_flags_preserved
    )

    return {
        "schema": "kokuno-agent4-source-signed-covariance-mass-independent-audit-v1",
        "seed": SEED,
        "case_order": list(case_order),
        "h_cases": list(H_CASES),
        "resolution_ladder_samples": list(LADDER),
        "fixed_fine_samples": FIXED_FINE_SAMPLES,
        "independent_gauss_legendre_order": GL_ORDER,
        "oracle": {
            "agent2_trapezoid_helper_reused": False,
            "agent2_band_covering_helper_reused": False,
            "agent2_signed_inverse_reused": False,
            "agent2_complete_curl_reused": False,
            "training_tensor_or_loss_read": False,
            "pressure_or_forcing_fit_used": False,
            "reference_quadrature": "numpy.polynomial.legendre.leggauss",
            "pulse_and_transverse_resolution_varied_separately": True,
        },
        "frozen_guards": dict(FROZEN_GUARDS),
        "metrics": {
            "pulse_relative_rms_ladder": pulse_ladder,
            "pulse_refinement_ratios": pulse_refinement,
            "transverse_relative_ladder": transverse_ladder,
            "transverse_refinement_ratios": transverse_refinement,
            "combined_h_sigma_relative_rms_ladder": combined_ladder,
            "combined_h_sigma_relative_max_ladder": combined_max_ladder,
            "combined_refinement_ratios": combined_refinement,
            "wrong_psi_power_mutation_relative_rms": wrong_psi_mutation,
            "wrong_transverse_prefactor_mutation_relative_rms": wrong_prefactor_mutation,
        },
        "cases": case_receipts,
        "local_structural_preflight_passed": bool(local_pass),
        "physical_coordinate_probe_applicable": False,
        "axis_near_probe_applicable": False,
        "formal_full_domain_pde_gate_assessed": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "truth_boundary": {
            "candidate_quadrature_only": True,
            "source_displayed_mass_formula_audited": True,
            "actual_source_pulse_samples_recovered": False,
            "actual_source_h_sigma_bound": False,
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "public_source_bound_velocity_osc_materialized": False,
            "public_velocity_correction_materialized": False,
            "full_ns_momentum_gate_value": None,
            "full_ns_divergence_gate_value": None,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = run_audit()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.out is None:
        print(text, end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
