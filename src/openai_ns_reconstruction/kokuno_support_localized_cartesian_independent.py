"""Independent Cartesian audit of the Kokuno support-localized complete curl.

Agent 2 localizes one source harmonic at the vector-potential level before
applying the cylindrical complete curl.  Its focused construction tests use a
cylindrical fourth-order finite-difference oracle.  This Agent-4 module uses a
different coordinate/operator path: it reconstructs the same *manufactured,
auxiliary-independent* vector potential directly in Cartesian coordinates and
differentiates that Cartesian potential with centered second-order stencils.

This is a structural audit of one local complex harmonic, not a full-domain
Navier--Stokes validation.  There is no global Kokuno leading+oscillatory+
correction velocity/pressure/fixed-forcing composite on this stack, so the
registered momentum <=1e-3 and divergence <=1e-5 project gates remain
unassessed and ``pde_validated`` remains false.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import ceil
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl


SCHEMA = "kokuno-agent4-support-localized-cartesian-audit-v1"
TASK_ID = "KOKUNO-A4-SUPPORT-LOCALIZED-CARTESIAN-AUDIT-017"
BASE_PR = 351
BASE_HEAD = "26112140d8a57873d20c6104ea19f32d8b9afa1d"
SEED = 9173081
FD_STEPS = (0.002, 0.001, 0.0005)
COEFFICIENT_DERIVATIVE_STEP = 4.0e-4
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# These are implementation-consistency guards for this local manufactured
# structural seam.  They are not replacements for the project PDE gates.
LOCAL_CURL_FINE_REL_GUARD = 2.0e-4
LOCAL_DIVERGENCE_FINE_ABS_GUARD = 2.0e-4
LOCAL_REFINEMENT_RATIO_GUARD = 2.5
MISSING_GRADIENT_RELATIVE_ERROR_FLOOR = 1.0e-3
MISSING_GRADIENT_DIVERGENCE_FLOOR = 1.0e-4


@dataclass(frozen=True)
class AuditCase:
    epsilon: float
    m: int
    label: str


CASES = (
    AuditCase(0.2, 1, "eps0p2_m1"),
    AuditCase(0.125, -2, "eps0p125_mminus2"),
    AuditCase(0.3, 2, "eps0p3_m2"),
)


def _fd6_first(fn: Callable[[float], np.ndarray], x: float, h: float) -> np.ndarray:
    """Sixth-order centered first derivative used only for caller derivative data."""
    return (
        -fn(x - 3.0 * h)
        + 9.0 * fn(x - 2.0 * h)
        - 45.0 * fn(x - h)
        + 45.0 * fn(x + h)
        - 9.0 * fn(x + 2.0 * h)
        + fn(x + 3.0 * h)
    ) / (60.0 * h)


def _phase_covector_transverse(
    epsilon: float, R: float, theta: float, Z: float
) -> tuple[float, np.ndarray, np.ndarray]:
    """Manufactured source-normalized phase/covector and transverse amplitude.

    The geometry is deliberately auxiliary independent.  Hence source D_r is
    ordinary partial_R here, while source D_z = epsilon*partial_Z.  The slow
    amplitude is theta independent, as required by the single-fast-phase curl
    contract.
    """
    p = 1.35
    p_z = 0.62
    x0 = 1.08
    pulse = 0.23
    a, b, c = 0.075, -0.052, 0.064
    H = a * R * R + b * R * Z + c * Z * Z
    H_R = 2.0 * a * R + b * Z
    H_Z = b * R + 2.0 * c * Z
    phase = p * theta + p_z * Z / epsilon + x0 * R - pulse * H
    n = np.asarray(
        (
            x0 - pulse * H_R,
            p / R,
            p_z - epsilon * pulse * H_Z,
        ),
        dtype=float,
    )
    raw = np.asarray(
        (
            0.71 + 0.065 * R - 0.018 * Z,
            -0.29 + 0.041 * Z + 0.012 * R,
            0.24 + 0.026 * R * Z - 0.009 * R,
        ),
        dtype=np.complex128,
    )
    n2 = float(np.dot(n, n))
    t = raw - n * (np.dot(n, raw) / n2)
    return float(phase), n, t


def _direct_coefficient(
    epsilon: float, m: int, R: float, theta: float, Z: float
) -> np.ndarray:
    """Recode C_m directly instead of calling the production coefficient helper."""
    _, n, t = _phase_covector_transverse(epsilon, R, theta, Z)
    k_m = ceil(epsilon ** -0.5) * int(m)
    return 1j * np.cross(n, t) / (k_m * float(np.dot(n, n)))


def _direct_coefficient_derivatives(
    epsilon: float, m: int, R: float, theta: float, Z: float
) -> tuple[np.ndarray, np.ndarray]:
    h = COEFFICIENT_DERIVATIVE_STEP
    D_r = _fd6_first(
        lambda rr: _direct_coefficient(epsilon, m, rr, theta, Z), R, h
    )
    D_z = epsilon * _fd6_first(
        lambda zz: _direct_coefficient(epsilon, m, R, theta, zz), Z, h
    )
    return D_r, D_z


def _eta(R: float, Z: float) -> float:
    # Smooth nonconstant support weight, bounded safely inside |eta|<1 on all
    # preregistered probes.  This is manufactured audit data, not a claimed
    # reconstruction of chi_ell or chi_{ell,a}.
    return float(0.54 + 0.085 * R - 0.047 * Z + 0.021 * R * Z)


def _eta_derivatives(epsilon: float, R: float, Z: float) -> tuple[float, float]:
    return float(0.085 + 0.021 * Z), float(epsilon * (-0.047 + 0.021 * R))


def _cyl_to_cart(vector: np.ndarray, theta: float) -> np.ndarray:
    vector = np.asarray(vector, dtype=np.complex128)
    ct = np.cos(theta)
    st = np.sin(theta)
    return np.asarray(
        (
            vector[0] * ct - vector[1] * st,
            vector[0] * st + vector[1] * ct,
            vector[2],
        ),
        dtype=np.complex128,
    )


def _xyz_to_source(epsilon: float, xyz: np.ndarray) -> tuple[float, float, float]:
    x, y, z_phys = np.asarray(xyz, dtype=float)
    R = float(np.hypot(x, y))
    if R <= 0.0:
        raise ValueError("Cartesian audit probes must stay off the cylindrical axis")
    theta = float(np.arctan2(y, x))
    # With Z = epsilon*z_phys, ordinary Cartesian d/dz_phys equals the source
    # normalized D_z = epsilon*d/dZ used by the complete-curl contract.
    Z = float(epsilon * z_phys)
    return R, theta, Z


def _direct_cartesian_potential(case: AuditCase, xyz: np.ndarray) -> np.ndarray:
    R, theta, Z = _xyz_to_source(case.epsilon, xyz)
    phase, _, _ = _phase_covector_transverse(case.epsilon, R, theta, Z)
    C = _direct_coefficient(case.epsilon, case.m, R, theta, Z)
    k_m = ceil(case.epsilon ** -0.5) * int(case.m)
    A_cyl = _eta(R, Z) * C * np.exp(1j * k_m * phase)
    return _cyl_to_cart(A_cyl, theta)


def _public_cartesian_velocity(
    case: AuditCase, xyz: np.ndarray, *, omit_support_gradients: bool = False
) -> np.ndarray:
    R, theta, Z = _xyz_to_source(case.epsilon, xyz)
    phase, n, t = _phase_covector_transverse(case.epsilon, R, theta, Z)
    D_r_C, D_z_C = _direct_coefficient_derivatives(
        case.epsilon, case.m, R, theta, Z
    )
    D_r_eta, D_z_eta = _eta_derivatives(case.epsilon, R, Z)
    if omit_support_gradients:
        D_r_eta = 0.0
        D_z_eta = 0.0
    localizer = KokunoSourceSupportLocalizedCurl(epsilon=case.epsilon, m=case.m)
    velocity_cyl = localizer.localized_mode(
        R,
        phase,
        n,
        t,
        D_r_C,
        D_z_C,
        _eta(R, Z),
        D_r_eta,
        D_z_eta,
    )["velocity"]
    return _cyl_to_cart(velocity_cyl, theta)


def _cartesian_curl_fd2(case: AuditCase, xyz: np.ndarray, step: float) -> np.ndarray:
    jacobian = np.empty((3, 3), dtype=np.complex128)
    for axis in range(3):
        plus = np.asarray(xyz, dtype=float).copy()
        minus = np.asarray(xyz, dtype=float).copy()
        plus[axis] += step
        minus[axis] -= step
        jacobian[:, axis] = (
            _direct_cartesian_potential(case, plus)
            - _direct_cartesian_potential(case, minus)
        ) / (2.0 * step)
    return np.asarray(
        (
            jacobian[2, 1] - jacobian[1, 2],
            jacobian[0, 2] - jacobian[2, 0],
            jacobian[1, 0] - jacobian[0, 1],
        ),
        dtype=np.complex128,
    )


def _cartesian_divergence_fd2(
    case: AuditCase,
    xyz: np.ndarray,
    step: float,
    *,
    omit_support_gradients: bool = False,
) -> complex:
    diagonal = []
    for axis in range(3):
        plus = np.asarray(xyz, dtype=float).copy()
        minus = np.asarray(xyz, dtype=float).copy()
        plus[axis] += step
        minus[axis] -= step
        vp = _public_cartesian_velocity(
            case, plus, omit_support_gradients=omit_support_gradients
        )
        vm = _public_cartesian_velocity(
            case, minus, omit_support_gradients=omit_support_gradients
        )
        diagonal.append((vp[axis] - vm[axis]) / (2.0 * step))
    return complex(sum(diagonal))


def _sample_probes() -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    probes: list[dict[str, Any]] = []
    for case in CASES:
        for region, count, r_bounds, z_bounds in (
            ("off_grid", 6, (0.56, 1.08), (-0.55, 0.55)),
            ("small_radius", 5, (0.18, 0.27), (-0.32, 0.32)),
        ):
            for _ in range(count):
                R = float(rng.uniform(*r_bounds))
                theta = float(rng.uniform(0.31, 2.69))
                z_phys = float(rng.uniform(*z_bounds))
                xyz = np.asarray(
                    (R * np.cos(theta), R * np.sin(theta), z_phys), dtype=float
                )
                probes.append(
                    {
                        "case": case,
                        "region": region,
                        "xyz": xyz,
                    }
                )
    return probes


def _rms(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


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
        oracle = _cartesian_curl_fd2(case, xyz, step)
        public = _public_cartesian_velocity(case, xyz)
        difference = float(np.linalg.norm(public - oracle))
        reference = float(np.linalg.norm(oracle))
        divergence = abs(_cartesian_divergence_fd2(case, xyz, step))
        curl_errors.append(difference)
        curl_reference_norms.append(reference)
        point_relative_errors.append(difference / max(reference, 1.0e-300))
        divergences.append(float(divergence))

        for key, name in ((region_metrics, probe["region"]), (case_metrics, case.label)):
            bucket = key.setdefault(name, {"curl_error": [], "reference": [], "divergence": []})
            bucket["curl_error"].append(difference)
            bucket["reference"].append(reference)
            bucket["divergence"].append(float(divergence))

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
        denom_region = max(_rms(bucket["reference"]), 1.0e-300)
        output["regions"][name] = {
            "sample_count": len(bucket["curl_error"]),
            "curl_relative_rms": _rms(bucket["curl_error"]) / denom_region,
            "divergence_max_abs": float(max(bucket["divergence"])),
            "divergence_rms_abs": _rms(bucket["divergence"]),
        }
    output["parameter_cases"] = {}
    for name, bucket in case_metrics.items():
        denom_case = max(_rms(bucket["reference"]), 1.0e-300)
        output["parameter_cases"][name] = {
            "sample_count": len(bucket["curl_error"]),
            "curl_relative_rms": _rms(bucket["curl_error"]) / denom_case,
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
        oracle = _cartesian_curl_fd2(case, xyz, step)
        mutated = _public_cartesian_velocity(case, xyz, omit_support_gradients=True)
        relative_errors.append(
            float(np.linalg.norm(mutated - oracle)) / max(float(np.linalg.norm(oracle)), 1.0e-300)
        )
        mutated_divergences.append(
            float(abs(_cartesian_divergence_fd2(case, xyz, step, omit_support_gradients=True)))
        )
        correct_divergences.append(
            float(abs(_cartesian_divergence_fd2(case, xyz, step)))
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
    ladder = [_evaluate_step(probes, step) for step in FD_STEPS]
    finest = ladder[-1]
    curl_ratios = [
        _ratio(ladder[index]["curl_relative_rms"], ladder[index + 1]["curl_relative_rms"])
        for index in range(len(ladder) - 1)
    ]
    divergence_ratios = [
        _ratio(ladder[index]["divergence_rms_abs"], ladder[index + 1]["divergence_rms_abs"])
        for index in range(len(ladder) - 1)
    ]
    mutation = _mutation_metrics(probes, FD_STEPS[-1])

    checks = {
        "finest_cartesian_curl_relative_rms": bool(
            finest["curl_relative_rms"] <= LOCAL_CURL_FINE_REL_GUARD
        ),
        "finest_cartesian_divergence_max_abs": bool(
            finest["divergence_max_abs"] <= LOCAL_DIVERGENCE_FINE_ABS_GUARD
        ),
        "curl_refinement_ratio": bool(min(curl_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD),
        "divergence_refinement_ratio": bool(
            min(divergence_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD
        ),
        "missing_support_gradient_curl_detected": bool(
            mutation["curl_point_relative_min"] >= MISSING_GRADIENT_RELATIVE_ERROR_FLOOR
        ),
        "missing_support_gradient_divergence_detected": bool(
            mutation["mutated_divergence_max_abs"] >= MISSING_GRADIENT_DIVERGENCE_FLOOR
        ),
    }
    structural_pass = bool(all(checks.values()))

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent2_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "artifact": "KokunoSourceSupportLocalizedCurl",
        },
        "seed": SEED,
        "parameter_cases": [
            {"label": case.label, "epsilon": case.epsilon, "m": case.m}
            for case in CASES
        ],
        "operator": {
            "oracle": (
                "recode manufactured vector potential in Cartesian coordinates, map "
                "Z=epsilon*z_phys, then centered Cartesian FD2 curl"
            ),
            "public_path": (
                "KokunoSourceSupportLocalizedCurl.localized_mode output only; "
                "converted from cylindrical to Cartesian"
            ),
            "caller_coefficient_derivatives": (
                "separately recoded C_m with sixth-order centered derivatives; "
                "not the Cartesian validation operator"
            ),
            "fd_steps": list(FD_STEPS),
            "coefficient_derivative_step": COEFFICIENT_DERIVATIVE_STEP,
            "construction_cylindrical_fd4_oracle_reused": False,
            "training_tensor_or_loss_read": False,
            "free_forcing_used": False,
        },
        "sampling": {
            "total_points": len(probes),
            "regions": {
                "off_grid": sum(probe["region"] == "off_grid" for probe in probes),
                "small_radius": sum(probe["region"] == "small_radius" for probe in probes),
            },
            "note": (
                "small-radius probes remain safely off R=0 because this local cylindrical "
                "single-harmonic contract does not provide an axis extension"
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
            "checks": checks,
            "structural_preflight_passed": structural_pass,
            "scope": "local manufactured support-localized curl implementation consistency only",
        },
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "reason": (
                "Agent 2 #351 still has no actual-source public xyz,t oscillatory velocity; "
                "Agent 1 segmented leading remains non-global and Agent 3 has no promoted "
                "independent correction cycle"
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
    target.write_text(json.dumps(build_report(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/support_localized_cartesian_audit.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    print(path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
