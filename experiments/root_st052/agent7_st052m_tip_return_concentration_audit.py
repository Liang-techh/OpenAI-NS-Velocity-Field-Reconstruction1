"""Audit whether the p=9,m=4 tip reshape concentrates unavoidable return flow.

Issue #749 fixes the diagnostic scope. This is a representation-capacity audit only:
it compares the axial radial-response factor of the frozen p=4,m=4 and p=9,m=4
one-channel compact poloidal envelopes. It does not change or select a velocity.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_tip_return_flow_relocation as relocation

TASK_ID = "CR003-ST052M-TIP-RETURN-CONCENTRATION-AUDIT-104"
ISSUE = 749
SOURCE_NONLINEAR_PR = 740
SOURCE_NONLINEAR_HEAD = "21019df839557e64b8d8689f68cdf8fc685a3a36"
SOURCE_RESHAPE_PR = 732
BASELINE_P = 4
RESHAPED_P = 9
M = 4
GRID_SIZES = (25001, 50001, 100001)
STABILITY_TOL = 1.0e-5
SIGNED_BALANCE_TOL = 1.0e-10

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "second_poloidal_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_fit_performed": False,
    "trajectory_replay_performed": False,
    "fresh_714_path_data_used": False,
    "held_out_pde_residual_evaluated": False,
    "pressure_or_force_changed": False,
    "public_image_numeric_target_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def axial_response(z: np.ndarray, p: int) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    bump, bump_q = relocation._asymmetric_axial_bump(z * z, int(p), M)
    return bump + 2.0 * z * z * bump_q


def _metrics(p: int, n: int) -> dict[str, float]:
    z = np.linspace(relocation.AXIAL_Z_INNER, relocation.AXIAL_Z_OUTER, int(n))
    g = axial_response(z, p)
    energy = float(np.trapezoid(g * g, z))
    if not (np.isfinite(energy) and energy > 0.0):
        raise RuntimeError("degenerate axial response energy")
    outward_energy = float(np.trapezoid(np.where(g < 0.0, g * g, 0.0), z))
    dg = np.gradient(g, z, edge_order=2)
    derivative_energy = float(np.trapezoid(dg * dg, z))
    rms = math.sqrt(energy / (z[-1] - z[0]))
    positive = float(np.trapezoid(np.where(g > 0.0, g, 0.0), z))
    negative = float(-np.trapezoid(np.where(g < 0.0, g, 0.0), z))
    q_star, z_star = relocation.physical_sign_change(int(p), M)
    return {
        "grid_points": int(n),
        "physical_q_sign_change": float(q_star),
        "physical_abs_z_sign_change": float(z_star),
        "return_collar_thickness": float(relocation.AXIAL_Z_OUTER - z_star),
        "response_energy": energy,
        "outward_squared_response_fraction": float(outward_energy / energy),
        "normalized_derivative_sharpness": float(math.sqrt(derivative_energy / energy)),
        "normalized_peak_response": float(np.max(np.abs(g)) / rms),
        "positive_signed_response_integral": positive,
        "negative_signed_response_magnitude": negative,
        "signed_integral_balance": float(positive - negative),
    }


def _cosine(n: int) -> float:
    z = np.linspace(relocation.AXIAL_Z_INNER, relocation.AXIAL_Z_OUTER, int(n))
    g0 = axial_response(z, BASELINE_P)
    g1 = axial_response(z, RESHAPED_P)
    e0 = float(np.trapezoid(g0 * g0, z))
    e1 = float(np.trapezoid(g1 * g1, z))
    cross = float(np.trapezoid(g0 * g1, z))
    return float(cross / math.sqrt(e0 * e1))


def evaluate() -> dict[str, Any]:
    selected_p, _, _ = relocation.select_inner_exponent()
    if selected_p != RESHAPED_P or relocation.OUTER_EXPONENT != M:
        raise RuntimeError("frozen #732 envelope identity drifted")

    levels: dict[str, Any] = {}
    for n in GRID_SIZES:
        levels[str(n)] = {
            "baseline_p4": _metrics(BASELINE_P, n),
            "reshaped_p9": _metrics(RESHAPED_P, n),
            "profile_cosine": _cosine(n),
        }

    medium = levels[str(GRID_SIZES[-2])]
    fine = levels[str(GRID_SIZES[-1])]
    stable_fields = (
        "outward_squared_response_fraction",
        "normalized_derivative_sharpness",
        "normalized_peak_response",
    )
    stability = {}
    for field in stable_fields:
        stability[field] = {
            "baseline_abs_change": abs(fine["baseline_p4"][field] - medium["baseline_p4"][field]),
            "reshaped_abs_change": abs(fine["reshaped_p9"][field] - medium["reshaped_p9"][field]),
        }
    stability["profile_cosine_abs_change"] = abs(
        fine["profile_cosine"] - medium["profile_cosine"]
    )
    stable = bool(
        all(
            item["baseline_abs_change"] <= STABILITY_TOL
            and item["reshaped_abs_change"] <= STABILITY_TOL
            for item in stability.values()
            if isinstance(item, dict)
        )
        and stability["profile_cosine_abs_change"] <= STABILITY_TOL
    )
    if not stable:
        raise RuntimeError("axial response metrics are not stable across frozen grid refinement")

    b = fine["baseline_p4"]
    r = fine["reshaped_p9"]
    if abs(b["signed_integral_balance"]) > SIGNED_BALANCE_TOL:
        raise RuntimeError("baseline compact-lobe signed balance failed")
    if abs(r["signed_integral_balance"]) > SIGNED_BALANCE_TOL:
        raise RuntimeError("reshaped compact-lobe signed balance failed")

    collar_thinner = bool(r["return_collar_thickness"] < b["return_collar_thickness"])
    outward_fraction_increased = bool(
        r["outward_squared_response_fraction"] > b["outward_squared_response_fraction"]
    )
    sharpness_increased = bool(
        r["normalized_derivative_sharpness"] > b["normalized_derivative_sharpness"]
    )
    peak_increased = bool(r["normalized_peak_response"] > b["normalized_peak_response"])
    concentrated = bool(collar_thinner and outward_fraction_increased and sharpness_increased)

    return {
        "task_id": TASK_ID,
        "issue": ISSUE,
        "source_nonlinear_pr": {
            "pr": SOURCE_NONLINEAR_PR,
            "head": SOURCE_NONLINEAR_HEAD,
        },
        "source_reshape_pr": SOURCE_RESHAPE_PR,
        "frozen_comparison": {
            "baseline": {"p": BASELINE_P, "m": M},
            "reshape": {"p": RESHAPED_P, "m": M},
            "z_interval": [relocation.AXIAL_Z_INNER, relocation.AXIAL_Z_OUTER],
            "response_factor": "G_p(z)=B_p,4(z^2)+2*z^2*dB_p,4/dq(z^2)",
            "grid_sizes": list(GRID_SIZES),
            "stability_tolerance": STABILITY_TOL,
        },
        "grid_refinement": levels,
        "stability": {**stability, "passes": stable},
        "finest_comparison": {
            "baseline_p4": b,
            "reshaped_p9": r,
            "profile_cosine": fine["profile_cosine"],
            "return_collar_reduction_fraction": float(
                1.0 - r["return_collar_thickness"] / b["return_collar_thickness"]
            ),
            "outward_squared_response_fraction_absolute_change": float(
                r["outward_squared_response_fraction"] - b["outward_squared_response_fraction"]
            ),
            "outward_squared_response_fraction_relative_change": float(
                r["outward_squared_response_fraction"]
                / b["outward_squared_response_fraction"]
                - 1.0
            ),
            "normalized_derivative_sharpness_ratio": float(
                r["normalized_derivative_sharpness"]
                / b["normalized_derivative_sharpness"]
            ),
            "normalized_peak_response_ratio": float(
                r["normalized_peak_response"] / b["normalized_peak_response"]
            ),
        },
        "decision": {
            "geometric_return_collar_thinner": collar_thinner,
            "outward_squared_response_fraction_increased": outward_fraction_increased,
            "normalized_derivative_sharpness_increased": sharpness_increased,
            "normalized_peak_response_increased": peak_increased,
            "return_flow_concentrated_by_reshape": concentrated,
            "geometric_collar_thinning_is_monotonic_free_gain": bool(
                collar_thinner and not (outward_fraction_increased or sharpness_increased)
            ),
            "increase_inner_exponent_again_justified_now": False,
            "second_poloidal_direction_added_now": False,
            "routing": (
                "Do not increase p again on collar-thickness evidence alone. "
                "Use #740 nonlinear outcome first; if it fails, the next capacity experiment "
                "should be one independently controlled compact poloidal return-flow direction "
                "rather than a still sharper one-channel axial envelope."
            ),
            "closer_visualization_delivery_this_increment": False,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/agent7-st052m-tip-return-concentration-audit/report.json"),
    )
    args = parser.parse_args()
    report = evaluate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
