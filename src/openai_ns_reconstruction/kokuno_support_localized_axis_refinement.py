"""Independent small-radius refinement of the support-localized Kokuno curl.

Agent 4 PR #353 used a Cartesian centered-FD2 oracle and retained a scientific
local rejection because the finest absolute divergence maximum was concentrated
in the small-radius probe set.  This follow-up keeps the same local acceptance
numbers but changes the validation derivative to Cartesian centered FD4 and
uses a smaller three-resolution ladder on fresh probes.

The manufactured geometry/public adapter is inherited only to hold the tested
object fixed.  The validation derivative below is not Agent 2's cylindrical
FD4 oracle and is not the previous Agent-4 Cartesian FD2 operator.  This is a
local structural audit, not the formal full-domain Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_support_localized_cartesian_independent import (
    CASES,
    AuditCase,
    _direct_cartesian_potential,
    _public_cartesian_velocity,
    _rms,
)

SCHEMA = "kokuno-agent4-support-localized-axis-refinement-v1"
TASK_ID = "KOKUNO-A4-SUPPORT-LOCALIZED-AXIS-REFINEMENT-018"
BASE_PR = 353
BASE_HEAD = "a9f24388245996eb7b0e2cc66e4847eae2833f3d"
SEED = 9173091
FD4_STEPS = (0.001, 0.0005, 0.00025)
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# Keep the preregistered #353 local guards unchanged.  In particular the
# previous failed divergence maximum is not relaxed in response to that result.
LOCAL_CURL_FINE_REL_GUARD = 2.0e-4
LOCAL_DIVERGENCE_FINE_ABS_GUARD = 2.0e-4
LOCAL_REFINEMENT_RATIO_GUARD = 2.5
MISSING_GRADIENT_RELATIVE_ERROR_FLOOR = 1.0e-3
MISSING_GRADIENT_DIVERGENCE_FLOOR = 1.0e-4


def _fd4_first(
    fn: Callable[[np.ndarray], np.ndarray],
    xyz: np.ndarray,
    axis: int,
    step: float,
) -> np.ndarray:
    """Centered fourth-order first derivative in one Cartesian coordinate."""
    point = np.asarray(xyz, dtype=float)
    m2 = point.copy()
    m1 = point.copy()
    p1 = point.copy()
    p2 = point.copy()
    m2[axis] -= 2.0 * step
    m1[axis] -= step
    p1[axis] += step
    p2[axis] += 2.0 * step
    return (fn(m2) - 8.0 * fn(m1) + 8.0 * fn(p1) - fn(p2)) / (12.0 * step)


def _cartesian_curl_fd4(case: AuditCase, xyz: np.ndarray, step: float) -> np.ndarray:
    jacobian = np.empty((3, 3), dtype=np.complex128)
    for axis in range(3):
        jacobian[:, axis] = _fd4_first(
            lambda point: _direct_cartesian_potential(case, point),
            xyz,
            axis,
            step,
        )
    return np.asarray(
        (
            jacobian[2, 1] - jacobian[1, 2],
            jacobian[0, 2] - jacobian[2, 0],
            jacobian[1, 0] - jacobian[0, 1],
        ),
        dtype=np.complex128,
    )


def _cartesian_divergence_fd4(
    case: AuditCase,
    xyz: np.ndarray,
    step: float,
    *,
    omit_support_gradients: bool = False,
) -> complex:
    diagonal: list[complex] = []
    for axis in range(3):
        derivative = _fd4_first(
            lambda point: _public_cartesian_velocity(
                case, point, omit_support_gradients=omit_support_gradients
            ),
            xyz,
            axis,
            step,
        )
        diagonal.append(complex(derivative[axis]))
    return complex(sum(diagonal))


def _sample_probes() -> list[dict[str, Any]]:
    """Fresh probes, with deliberate extra weight on the rejected small-R seam."""
    rng = np.random.default_rng(SEED)
    probes: list[dict[str, Any]] = []
    for case in CASES:
        for region, count, r_bounds, z_bounds in (
            ("off_grid", 5, (0.56, 1.08), (-0.55, 0.55)),
            ("small_radius", 7, (0.14, 0.28), (-0.32, 0.32)),
        ):
            for _ in range(count):
                radius = float(rng.uniform(*r_bounds))
                theta = float(rng.uniform(0.25, 2.85))
                z_phys = float(rng.uniform(*z_bounds))
                xyz = np.asarray(
                    (radius * np.cos(theta), radius * np.sin(theta), z_phys),
                    dtype=float,
                )
                probes.append({"case": case, "region": region, "xyz": xyz})
    return probes


def _evaluate_step(probes: list[dict[str, Any]], step: float) -> dict[str, Any]:
    curl_errors: list[float] = []
    curl_reference_norms: list[float] = []
    point_relative_errors: list[float] = []
    divergences: list[float] = []
    region_metrics: dict[str, dict[str, list[float]]] = {}
    case_metrics: dict[str, dict[str, list[float]]] = {}

    for probe in probes:
        case = probe["case"]
        xyz = probe["xyz"]
        oracle = _cartesian_curl_fd4(case, xyz, step)
        public = _public_cartesian_velocity(case, xyz)
        difference = float(np.linalg.norm(public - oracle))
        reference = float(np.linalg.norm(oracle))
        divergence = float(abs(_cartesian_divergence_fd4(case, xyz, step)))
        curl_errors.append(difference)
        curl_reference_norms.append(reference)
        point_relative_errors.append(difference / max(reference, 1.0e-300))
        divergences.append(divergence)

        for target, name in (
            (region_metrics, probe["region"]),
            (case_metrics, case.label),
        ):
            bucket = target.setdefault(
                name, {"curl_error": [], "reference": [], "divergence": []}
            )
            bucket["curl_error"].append(difference)
            bucket["reference"].append(reference)
            bucket["divergence"].append(divergence)

    denom = max(_rms(curl_reference_norms), 1.0e-300)
    output: dict[str, Any] = {
        "step": float(step),
        "sample_count": len(probes),
        "curl_error_rms": _rms(curl_errors),
        "curl_reference_rms": _rms(curl_reference_norms),
        "curl_relative_rms": _rms(curl_errors) / denom,
        "curl_point_relative_max": float(max(point_relative_errors)),
        "divergence_max_abs": float(max(divergences)),
        "divergence_rms_abs": _rms(divergences),
    }

    output["regions"] = {}
    for name, bucket in region_metrics.items():
        region_denom = max(_rms(bucket["reference"]), 1.0e-300)
        output["regions"][name] = {
            "sample_count": len(bucket["curl_error"]),
            "curl_relative_rms": _rms(bucket["curl_error"]) / region_denom,
            "divergence_max_abs": float(max(bucket["divergence"])),
            "divergence_rms_abs": _rms(bucket["divergence"]),
        }

    output["parameter_cases"] = {}
    for name, bucket in case_metrics.items():
        case_denom = max(_rms(bucket["reference"]), 1.0e-300)
        output["parameter_cases"][name] = {
            "sample_count": len(bucket["curl_error"]),
            "curl_relative_rms": _rms(bucket["curl_error"]) / case_denom,
            "divergence_max_abs": float(max(bucket["divergence"])),
            "divergence_rms_abs": _rms(bucket["divergence"]),
        }
    return output


def _mutation_metrics(probes: list[dict[str, Any]], step: float) -> dict[str, Any]:
    relative_errors: list[float] = []
    mutated_divergences: list[float] = []
    correct_divergences: list[float] = []
    for probe in probes:
        case = probe["case"]
        xyz = probe["xyz"]
        oracle = _cartesian_curl_fd4(case, xyz, step)
        mutated = _public_cartesian_velocity(case, xyz, omit_support_gradients=True)
        relative_errors.append(
            float(np.linalg.norm(mutated - oracle))
            / max(float(np.linalg.norm(oracle)), 1.0e-300)
        )
        mutated_divergences.append(
            float(
                abs(
                    _cartesian_divergence_fd4(
                        case, xyz, step, omit_support_gradients=True
                    )
                )
            )
        )
        correct_divergences.append(
            float(abs(_cartesian_divergence_fd4(case, xyz, step)))
        )
    return {
        "mutation": "set D_r_eta=D_z_eta=0 while retaining the same nonconstant eta",
        "curl_point_relative_min": float(min(relative_errors)),
        "curl_point_relative_rms": _rms(relative_errors),
        "curl_point_relative_max": float(max(relative_errors)),
        "mutated_divergence_max_abs": float(max(mutated_divergences)),
        "mutated_divergence_rms_abs": _rms(mutated_divergences),
        "correct_divergence_max_abs_same_step": float(max(correct_divergences)),
        "correct_divergence_rms_abs_same_step": _rms(correct_divergences),
    }


def _ratio(coarse: float, fine: float) -> float:
    return float(coarse / max(fine, 1.0e-300))


def build_report() -> dict[str, Any]:
    probes = _sample_probes()
    ladder = [_evaluate_step(probes, step) for step in FD4_STEPS]
    finest = ladder[-1]
    curl_ratios = [
        _ratio(ladder[index]["curl_relative_rms"], ladder[index + 1]["curl_relative_rms"])
        for index in range(len(ladder) - 1)
    ]
    divergence_ratios = [
        _ratio(
            ladder[index]["divergence_rms_abs"],
            ladder[index + 1]["divergence_rms_abs"],
        )
        for index in range(len(ladder) - 1)
    ]
    mutation = _mutation_metrics(probes, FD4_STEPS[-1])

    checks = {
        "finest_cartesian_curl_relative_rms": bool(
            finest["curl_relative_rms"] <= LOCAL_CURL_FINE_REL_GUARD
        ),
        "finest_cartesian_divergence_max_abs": bool(
            finest["divergence_max_abs"] <= LOCAL_DIVERGENCE_FINE_ABS_GUARD
        ),
        "curl_refinement_ratio": bool(
            min(curl_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD
        ),
        "divergence_refinement_ratio": bool(
            min(divergence_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD
        ),
        "missing_support_gradient_curl_detected": bool(
            mutation["curl_point_relative_min"]
            >= MISSING_GRADIENT_RELATIVE_ERROR_FLOOR
        ),
        "missing_support_gradient_divergence_detected": bool(
            mutation["mutated_divergence_max_abs"]
            >= MISSING_GRADIENT_DIVERGENCE_FLOOR
        ),
    }
    structural_pass = bool(all(checks.values()))

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent4_parent_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "agent2_artifact": "KokunoSourceSupportLocalizedCurl",
        },
        "seed": SEED,
        "parameter_cases": [
            {"label": case.label, "epsilon": case.epsilon, "m": case.m}
            for case in CASES
        ],
        "operator": {
            "oracle": (
                "independently recoded Cartesian vector potential differentiated by "
                "centered Cartesian FD4"
            ),
            "public_path": (
                "KokunoSourceSupportLocalizedCurl.localized_mode output converted "
                "from cylindrical to Cartesian"
            ),
            "fd4_steps": list(FD4_STEPS),
            "construction_cylindrical_fd4_oracle_reused": False,
            "previous_agent4_cartesian_fd2_validation_operator_reused": False,
            "manufactured_geometry_held_fixed_from_parent_audit": True,
            "training_tensor_or_loss_read": False,
            "free_forcing_used": False,
        },
        "sampling": {
            "total_points": len(probes),
            "regions": {
                "off_grid": sum(probe["region"] == "off_grid" for probe in probes),
                "small_radius": sum(
                    probe["region"] == "small_radius" for probe in probes
                ),
            },
            "fresh_seed_relative_to_parent": True,
            "note": (
                "small-radius points stay safely off R=0; this local harmonic has no "
                "claimed axis extension"
            ),
        },
        "resolution_ladder": ladder,
        "refinement": {
            "curl_relative_rms_coarse_to_fine_ratios": curl_ratios,
            "divergence_rms_coarse_to_fine_ratios": divergence_ratios,
        },
        "mutation": mutation,
        "local_guards": {
            "thresholds": {
                "finest_curl_relative_rms_max": LOCAL_CURL_FINE_REL_GUARD,
                "finest_divergence_max_abs_max": LOCAL_DIVERGENCE_FINE_ABS_GUARD,
                "minimum_refinement_ratio": LOCAL_REFINEMENT_RATIO_GUARD,
                "missing_gradient_relative_error_min": MISSING_GRADIENT_RELATIVE_ERROR_FLOOR,
                "missing_gradient_divergence_max_abs_min": MISSING_GRADIENT_DIVERGENCE_FLOOR,
            },
            "same_thresholds_as_parent_rejection": True,
            "checks": checks,
            "structural_preflight_passed": structural_pass,
            "scope": (
                "local manufactured small-radius/support-localized curl implementation "
                "consistency only"
            ),
        },
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "reason": (
                "no global matched Kokuno leading pressure, actual-source public xyz,t "
                "oscillatory velocity, or promoted after-correction composite exists"
            ),
        },
        "truth_boundary": {
            "actual_source_partition_bumps_reconstructed": False,
            "actual_source_public_oscillatory_xyz_t_velocity_available": False,
            "global_leading_velocity_pressure_available": False,
            "after_correction_global_velocity_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def write_report(path: str | Path) -> Path:
    target = Path(path)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing audit report: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(build_report(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/support_localized_axis_refinement.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    print(path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
