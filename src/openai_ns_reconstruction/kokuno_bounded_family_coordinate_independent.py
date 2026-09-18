"""Independent black-box audit of Agent 2's bounded family coordinates.

This module deliberately does not read Agent 2's public tangent arrays.  It
constructs a fresh multi-band physical by-beta velocity family, calls only the
public bounded-coordinate value map, and reconstructs velocity/covariance
responses from modulated totals.  It is a local interface/numerical preflight,
not a Navier--Stokes residual validation and not source-family evidence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates

SCHEMA = "kokuno-agent4-bounded-family-coordinate-independent-audit-v1"
TASK_ID = "KOKUNO-A4-BOUNDED-FAMILY-COORDINATE-INDEPENDENT-AUDIT-022"
BASE_PR = 400
BASE_HEAD = "d4464acc3e7072f473746df9b02df823cdbcb8c4"
SEED = 9173131
N_RADII = 7
N_ANGLES = 24
RESPONSE_STEPS = (0.04, 0.02, 0.01)
TANGENT_STEP = 0.03
MAX_L1_UPDATE = 0.125

# Frozen before this audit is executed.  These are local coordinate-map guards,
# not replacements for the formal PDE gates.
COVARIANCE_FINE_RELATIVE_RMS_GUARD = 1.0e-2
COVARIANCE_MIN_REFINEMENT_RATIO_GUARD = 1.8
COMMON_QUADRATIC_IDENTITY_RELATIVE_GUARD = 5.0e-13
PERMUTATION_RELATIVE_RMS_GUARD = 5.0e-13
MINIMUM_LABEL_MULTIPLIER_GUARD = 1.0 - MAX_L1_UPDATE - 5.0e-13
LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR = 1.0e-1
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

EXPECTED_PARAMETER_NAMES = ("delta_common", "delta_band")
EXPECTED_PARAMETER_UNIT = "dimensionless_fractional_multiplier_of_Q_scaled_physical_velocity"


def _vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (2,):
        raise ValueError("covariance vectors must have final dimension 2")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _velocity_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("velocity vectors must have final dimension 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _fresh_multiband_family() -> dict[str, Any]:
    """Return a deterministic fresh physical family unrelated to Agent-2 curl code."""
    rng = np.random.default_rng(SEED)
    labels = (
        (5, (0, 0, 0)),
        (5, (1, 0, 0)),
        (6, (0, 0, 0)),
        (6, (1, 0, 0)),
        (7, (0, 0, 0)),
        (7, (1, 0, 0)),
    )
    contrast = np.asarray((-1.0, -1.0, 0.0, 0.0, 1.0, 1.0), dtype=float)
    radii = np.sort(rng.uniform(0.11, 0.47, N_RADII))
    angle_offset = float(rng.uniform(0.1, 0.9))
    angles = 2.0 * np.pi * (np.arange(N_ANGLES, dtype=float) + angle_offset) / N_ANGLES
    phase = 0.23 * (np.arange(6, dtype=float) + 1.0) + rng.uniform(-0.04, 0.04, 6)
    amplitude = 0.35 + 0.06 * (np.arange(6, dtype=float) + 1.0) + rng.uniform(-0.015, 0.015, 6)

    r = radii[:, None, None]
    theta = angles[None, :, None]
    label_index = np.arange(6, dtype=float)[None, None, :]
    g = contrast[None, None, :]
    p = phase[None, None, :]
    a = amplitude[None, None, :]

    radial = a * (1.0 + 0.12 * r) * np.cos(theta + p) + 0.025 * (label_index - 2.0)
    azimuthal = (
        a * (0.62 + 0.11 * g) * (1.0 + 0.05 * r) * np.cos(theta + p + 0.31)
        + 0.035 * (1.0 + g)
        + 0.0 * theta
    )
    axial = (
        a * (0.38 - 0.08 * g) * (1.0 - 0.04 * r) * np.cos(theta + p - 0.22)
        + 0.02 * (label_index - 1.5)
        + 0.0 * theta
    )
    velocity = np.stack((radial, azimuthal, axial), axis=-1)
    return {
        "velocity": velocity,
        "labels": labels,
        "independent_expected_contrast": contrast,
        "radii": radii,
        "angles": angles,
    }


def _public_total(
    contract: KokunoBoundedFamilyCoefficientCoordinates,
    velocity: np.ndarray,
    labels: tuple[Any, ...],
    *,
    delta_common: float = 0.0,
    delta_band: float = 0.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    out = contract.evaluate(
        velocity,
        labels,
        delta_common=delta_common,
        delta_band=delta_band,
    )
    total = np.asarray(out["velocity_physical_cylindrical_total_modulated"], dtype=float)
    return total, out


def _covariance(total: np.ndarray) -> np.ndarray:
    total = np.asarray(total, dtype=float)
    if total.shape != (N_RADII, N_ANGLES, 3):
        raise ValueError("unexpected total-velocity shape")
    return np.stack(
        (
            np.mean(total[..., 0] * total[..., 1], axis=1),
            np.mean(total[..., 0] * total[..., 2], axis=1),
        ),
        axis=-1,
    )


def _product_rule_covariance_response(base: np.ndarray, tangent: np.ndarray) -> np.ndarray:
    return np.stack(
        (
            np.mean(tangent[..., 0] * base[..., 1] + base[..., 0] * tangent[..., 1], axis=1),
            np.mean(tangent[..., 0] * base[..., 2] + base[..., 0] * tangent[..., 2], axis=1),
        ),
        axis=-1,
    )


def _black_box_velocity_tangent(
    contract: KokunoBoundedFamilyCoefficientCoordinates,
    velocity: np.ndarray,
    labels: tuple[Any, ...],
    direction: str,
) -> np.ndarray:
    kwargs_plus = {"delta_common": 0.0, "delta_band": 0.0}
    kwargs_minus = dict(kwargs_plus)
    key = "delta_common" if direction == "common" else "delta_band"
    kwargs_plus[key] = TANGENT_STEP
    kwargs_minus[key] = -TANGENT_STEP
    plus, _ = _public_total(contract, velocity, labels, **kwargs_plus)
    minus, _ = _public_total(contract, velocity, labels, **kwargs_minus)
    return (plus - minus) / (2.0 * TANGENT_STEP)


def _direction_report(
    contract: KokunoBoundedFamilyCoefficientCoordinates,
    velocity: np.ndarray,
    labels: tuple[Any, ...],
    base_total: np.ndarray,
    base_covariance: np.ndarray,
    direction: str,
) -> dict[str, Any]:
    tangent = _black_box_velocity_tangent(contract, velocity, labels, direction)
    reference = _product_rule_covariance_response(base_total, tangent)
    ladder: list[dict[str, float]] = []
    for step in RESPONSE_STEPS:
        kwargs = {"delta_common": 0.0, "delta_band": 0.0}
        kwargs["delta_common" if direction == "common" else "delta_band"] = step
        perturbed, _ = _public_total(contract, velocity, labels, **kwargs)
        response = (_covariance(perturbed) - base_covariance) / step
        relative = _vector_rms(response - reference) / max(_vector_rms(reference), np.finfo(float).tiny)
        ladder.append(
            {
                "step": float(step),
                "response_relative_rms": float(relative),
                "response_vector_rms": _vector_rms(response),
            }
        )
    errors = [item["response_relative_rms"] for item in ladder]
    refinement = [errors[i] / max(errors[i + 1], np.finfo(float).tiny) for i in range(2)]
    return {
        "direction": direction,
        "black_box_velocity_tangent_rms": _velocity_rms(tangent),
        "reference_covariance_response_rms": _vector_rms(reference),
        "reference_covariance_response": reference.tolist(),
        "forward_response_ladder": ladder,
        "refinement_ratios": [float(x) for x in refinement],
    }


def generate_report() -> dict[str, Any]:
    family = _fresh_multiband_family()
    velocity = np.asarray(family["velocity"], dtype=float)
    labels = tuple(family["labels"])
    expected_contrast = np.asarray(family["independent_expected_contrast"], dtype=float)
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=MAX_L1_UPDATE)

    base_total, base_out = _public_total(contract, velocity, labels)
    base_covariance = _covariance(base_total)
    names = tuple(base_out["parameter_names"])
    units = tuple(base_out["parameter_units"])
    public_contrast = np.asarray(base_out["band_contrast"], dtype=float)

    common = _direction_report(contract, velocity, labels, base_total, base_covariance, "common")
    band = _direction_report(contract, velocity, labels, base_total, base_covariance, "band")
    common_reference = np.asarray(common["reference_covariance_response"], dtype=float)
    common_identity_relative = _vector_rms(common_reference - 2.0 * base_covariance) / max(
        2.0 * _vector_rms(base_covariance), np.finfo(float).tiny
    )

    # Correct simultaneous permutation of labels and velocity must be invisible.
    permutation = np.asarray((4, 5, 2, 3, 0, 1), dtype=int)
    permuted_velocity = velocity[..., permutation, :]
    permuted_labels = tuple(labels[i] for i in permutation)
    probe_kwargs = {"delta_common": 0.03, "delta_band": -0.04}
    original_probe, _ = _public_total(contract, velocity, labels, **probe_kwargs)
    permuted_probe, _ = _public_total(contract, permuted_velocity, permuted_labels, **probe_kwargs)
    permutation_relative = _velocity_rms(permuted_probe - original_probe) / max(
        _velocity_rms(original_probe), np.finfo(float).tiny
    )

    # At the aggregate L1 boundary every label must stay positive.
    boundary_updates = (
        (MAX_L1_UPDATE, 0.0),
        (-MAX_L1_UPDATE, 0.0),
        (0.0, MAX_L1_UPDATE),
        (0.0, -MAX_L1_UPDATE),
        (0.5 * MAX_L1_UPDATE, 0.5 * MAX_L1_UPDATE),
        (-0.5 * MAX_L1_UPDATE, -0.5 * MAX_L1_UPDATE),
    )
    boundary_minimum = 1.0
    for d0, d1 in boundary_updates:
        _, out = _public_total(contract, velocity, labels, delta_common=d0, delta_band=d1)
        boundary_minimum = min(boundary_minimum, float(out["minimum_label_multiplier"]))

    # Mutation: keep the public labels fixed while swapping the outer ell-band
    # velocity blocks.  The total/common direction is unchanged, but the band
    # response must move substantially if the oracle is label-sensitive.
    mutated_velocity = velocity[..., permutation, :]
    mutated_base, _ = _public_total(contract, mutated_velocity, labels)
    mutated_band_tangent = _black_box_velocity_tangent(contract, mutated_velocity, labels, "band")
    mutated_band_response = _product_rule_covariance_response(mutated_base, mutated_band_tangent)
    reference_band = np.asarray(band["reference_covariance_response"], dtype=float)
    mutation_relative = _vector_rms(mutated_band_response - reference_band) / max(
        _vector_rms(reference_band), np.finfo(float).tiny
    )

    checks = {
        "parameter_contract_exact": names == EXPECTED_PARAMETER_NAMES
        and units == (EXPECTED_PARAMETER_UNIT, EXPECTED_PARAMETER_UNIT),
        "independent_band_contrast_matches_public": bool(np.array_equal(public_contrast, expected_contrast)),
        "common_quadratic_identity": common_identity_relative <= COMMON_QUADRATIC_IDENTITY_RELATIVE_GUARD,
        "common_fine_covariance_response": common["forward_response_ladder"][-1]["response_relative_rms"]
        <= COVARIANCE_FINE_RELATIVE_RMS_GUARD,
        "band_fine_covariance_response": band["forward_response_ladder"][-1]["response_relative_rms"]
        <= COVARIANCE_FINE_RELATIVE_RMS_GUARD,
        "common_response_refinement": min(common["refinement_ratios"])
        >= COVARIANCE_MIN_REFINEMENT_RATIO_GUARD,
        "band_response_refinement": min(band["refinement_ratios"])
        >= COVARIANCE_MIN_REFINEMENT_RATIO_GUARD,
        "label_permutation_invariance": permutation_relative <= PERMUTATION_RELATIVE_RMS_GUARD,
        "aggregate_boundary_stays_positive": boundary_minimum >= MINIMUM_LABEL_MULTIPLIER_GUARD,
        "label_misalignment_mutation_detected": mutation_relative
        >= LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR,
    }

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_pr": BASE_PR,
        "base_head": BASE_HEAD,
        "seed": SEED,
        "sample_contract": {
            "radial_count": N_RADII,
            "angular_count": N_ANGLES,
            "beta_count": len(labels),
            "beta_labels": [list(label) for label in labels],
            "response_steps": list(RESPONSE_STEPS),
            "black_box_tangent_step": TANGENT_STEP,
            "uses_public_tangent_arrays": False,
            "uses_training_tensor_or_loss": False,
            "uses_pressure_or_forcing_fit": False,
        },
        "guards": {
            "covariance_fine_relative_rms": COVARIANCE_FINE_RELATIVE_RMS_GUARD,
            "covariance_min_refinement_ratio": COVARIANCE_MIN_REFINEMENT_RATIO_GUARD,
            "common_quadratic_identity_relative": COMMON_QUADRATIC_IDENTITY_RELATIVE_GUARD,
            "label_permutation_relative_rms": PERMUTATION_RELATIVE_RMS_GUARD,
            "minimum_label_multiplier": MINIMUM_LABEL_MULTIPLIER_GUARD,
            "label_misalignment_response_relative_floor": LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR,
            "formal_momentum_gate_unchanged": FORMAL_MOMENTUM_GATE,
            "formal_divergence_gate_unchanged": FORMAL_DIVERGENCE_GATE,
        },
        "base_covariance_vector_rms": _vector_rms(base_covariance),
        "common_direction": common,
        "band_direction": band,
        "common_quadratic_identity_relative_error": float(common_identity_relative),
        "label_permutation_relative_velocity_rms": float(permutation_relative),
        "minimum_boundary_label_multiplier": float(boundary_minimum),
        "label_misalignment_band_response_relative_change": float(mutation_relative),
        "checks": checks,
        "structural_preflight_passed": bool(all(checks.values())),
        "truth_boundary": {
            "repository_autonomous_coordinate_map_independently_preflighted": bool(all(checks.values())),
            "actual_source_positive_order_background_bound": False,
            "actual_source_xyz_t_oscillatory_velocity_ready": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/kokuno_agent4/bounded_family_coordinate_independent_audit.json"),
    )
    args = parser.parse_args()
    report = generate_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
