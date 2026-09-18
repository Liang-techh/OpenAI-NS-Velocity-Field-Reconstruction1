"""Independent Cartesian audit on Kokuno's source-grid/source-wave-shell seam.

Agent 2 PR #358 bound the support-localized complete curl to the corrected
reader's source angular grid ``p=j/k`` and active wave shell ``X>=4/Lambda``;
PR #369 additionally made the slow squared-partition family contract
executable.  Earlier Agent-4 Cartesian checks deliberately used a generic
manufactured phase, including ``p=1.35``, so they did not validate this newer
source-domain wrapper.

This module closes only that seam.  It uses source-grid-compatible labels and
held-out points inside the source wave shell, but reconstructs the same local
vector potential independently in Cartesian coordinates and differentiates it
with centered Cartesian FD4.  The public path is only
``KokunoSourceWaveShellLocalizedCurl.localized_mode``.  The manufactured slow
weight is one member of an exactly squared-partition-embeddable pair
``(cos(alpha), sin(alpha))``; it is *not* claimed to be the source's unrecovered
``chi_ell*chi_{ell,a}``.

This remains a local structural preflight, not the formal full-domain
Navier--Stokes gate.  No global Kokuno leading pressure, source-derived public
oscillatory xyz,t velocity, or promoted after-correction composite exists on
this stack, so ``pde_validated`` remains false.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from math import ceil, sqrt
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from .kokuno_source_wave_shell import KokunoSourceWaveShellLocalizedCurl

SCHEMA = "kokuno-agent4-source-wave-shell-cartesian-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-WAVE-SHELL-CARTESIAN-AUDIT-019"
BASE_PR = 369
BASE_HEAD = "6050ea435536d01600d90a9dc02ebbab27293b94"
SOURCE_WAVE_SHELL_PR = 358
SOURCE_WAVE_SHELL_HEAD = "e24d918bb0cc16b1f1006710183058296fbf760f"
SEED = 9173101
FD4_STEPS = (0.001, 0.0005, 0.00025)
COEFFICIENT_DERIVATIVE_STEP = 4.0e-4
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# Frozen local guards inherited unchanged from Agent 4 PRs #353/#360.
LOCAL_CURL_FINE_REL_GUARD = 2.0e-4
LOCAL_DIVERGENCE_FINE_ABS_GUARD = 2.0e-4
LOCAL_REFINEMENT_RATIO_GUARD = 2.5
MISSING_GRADIENT_RELATIVE_ERROR_FLOOR = 1.0e-3
MISSING_GRADIENT_DIVERGENCE_FLOOR = 1.0e-4


@dataclass(frozen=True)
class SourceAuditCase:
    epsilon: float
    m: int
    j: int
    Lambda: float
    p_z: float
    x0: float
    label: str

    @property
    def k(self) -> int:
        return int(ceil(self.epsilon ** -0.5))

    @property
    def p(self) -> float:
        return self.j / self.k

    @property
    def k_m(self) -> int:
        return self.k * self.m

    @property
    def azimuthal_mode(self) -> int:
        return self.m * self.j

    @property
    def phase_contract(self) -> KokunoOscillatoryPhaseContract:
        return KokunoOscillatoryPhaseContract(
            p=self.p,
            p_z=self.p_z,
            x0=self.x0,
            epsilon=self.epsilon,
            m=self.m,
        )

    @property
    def shell(self) -> KokunoSourceWaveShellLocalizedCurl:
        return KokunoSourceWaveShellLocalizedCurl(Lambda=self.Lambda)


CASES = (
    SourceAuditCase(0.2, 1, 1, 36.0, 0.62, 1.08, "eps0p2_m1_j1_L36"),
    SourceAuditCase(0.125, -2, -2, 64.0, 0.51, 0.94, "eps0p125_mminus2_jminus2_L64"),
    SourceAuditCase(0.3, 2, 1, 49.0, 0.73, 1.17, "eps0p3_m2_j1_L49"),
)


def _fd6_first(fn: Callable[[float], np.ndarray], x: float, h: float) -> np.ndarray:
    """Sixth-order centered derivative used only to supply caller C derivatives."""
    return (
        -fn(x - 3.0 * h)
        + 9.0 * fn(x - 2.0 * h)
        - 45.0 * fn(x - h)
        + 45.0 * fn(x + h)
        - 9.0 * fn(x + 2.0 * h)
        + fn(x + 3.0 * h)
    ) / (60.0 * h)


def _phase_covector_transverse(
    case: SourceAuditCase, R: float, theta: float, Z: float
) -> tuple[float, np.ndarray, np.ndarray]:
    """Independent source-form phase/covector with a slow manufactured H_Phi.

    The source formula is
    ``Phi=p*theta+p_z*Z/epsilon+x0*R-v*H_Phi``.  Here ``H_Phi`` is a smooth
    autonomous test polynomial and the pulse amplitude is fixed.  It can be
    realized by ``G=0, F=H_Phi/p, V=R*F``; no Agent-1 background tensor is read.
    """
    pulse = 0.23
    a, b, c = 0.075, -0.052, 0.064
    H = a * R * R + b * R * Z + c * Z * Z
    H_R = 2.0 * a * R + b * Z
    H_Z = b * R + 2.0 * c * Z
    phase = case.p * theta + case.p_z * Z / case.epsilon + case.x0 * R - pulse * H
    n = np.asarray(
        (
            case.x0 - pulse * H_R,
            case.p / R,
            case.p_z - case.epsilon * pulse * H_Z,
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
    t = raw - n * (np.dot(n, raw) / float(np.dot(n, n)))
    return float(phase), n, t


def _direct_coefficient(
    case: SourceAuditCase, R: float, theta: float, Z: float
) -> np.ndarray:
    """Independently recode C_m=i(n_Phi x t_m)/(k_m |n_Phi|^2)."""
    _, n, t = _phase_covector_transverse(case, R, theta, Z)
    return 1j * np.cross(n, t) / (case.k_m * float(np.dot(n, n)))


def _direct_coefficient_derivatives(
    case: SourceAuditCase, R: float, theta: float, Z: float
) -> tuple[np.ndarray, np.ndarray]:
    h = COEFFICIENT_DERIVATIVE_STEP
    D_r = _fd6_first(lambda rr: _direct_coefficient(case, rr, theta, Z), R, h)
    D_z = case.epsilon * _fd6_first(
        lambda zz: _direct_coefficient(case, R, theta, zz), Z, h
    )
    return D_r, D_z


def _partition_alpha(R: float, Z: float) -> float:
    return float(0.61 + 0.17 * R - 0.09 * Z + 0.03 * R * Z)


def _eta(R: float, Z: float) -> float:
    """One exact member of the manufactured pair (cos(alpha), sin(alpha))."""
    return float(np.cos(_partition_alpha(R, Z)))


def _eta_derivatives(epsilon: float, R: float, Z: float) -> tuple[float, float]:
    alpha = _partition_alpha(R, Z)
    alpha_R = 0.17 + 0.03 * Z
    alpha_Z = -0.09 + 0.03 * R
    return (
        float(-np.sin(alpha) * alpha_R),
        float(epsilon * (-np.sin(alpha) * alpha_Z)),
    )


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


def _xyz_to_source(case: SourceAuditCase, xyz: np.ndarray) -> tuple[float, float, float]:
    x, y, z_phys = np.asarray(xyz, dtype=float)
    R = float(np.hypot(x, y))
    if R <= 0.0:
        raise ValueError("source-wave-shell audit probes must remain away from the axis")
    theta = float(np.arctan2(y, x))
    # Z=epsilon*z_phys makes physical d/dz exactly the source D_z=epsilon*d/dZ.
    Z = float(case.epsilon * z_phys)
    return R, theta, Z


def _direct_cartesian_potential(case: SourceAuditCase, xyz: np.ndarray) -> np.ndarray:
    R, theta, Z = _xyz_to_source(case, xyz)
    phase, _, _ = _phase_covector_transverse(case, R, theta, Z)
    C = _direct_coefficient(case, R, theta, Z)
    A_cyl = _eta(R, Z) * C * np.exp(1j * case.k_m * phase)
    return _cyl_to_cart(A_cyl, theta)


def _public_cartesian_velocity(
    case: SourceAuditCase,
    s_Q: float,
    xyz: np.ndarray,
    *,
    omit_support_gradients: bool = False,
) -> np.ndarray:
    """Use only the public #358 source-wave-shell localized-mode contract."""
    R, theta, Z = _xyz_to_source(case, xyz)
    phase, n, t = _phase_covector_transverse(case, R, theta, Z)
    D_r_C, D_z_C = _direct_coefficient_derivatives(case, R, theta, Z)
    D_r_eta, D_z_eta = _eta_derivatives(case.epsilon, R, Z)
    if omit_support_gradients:
        D_r_eta = 0.0
        D_z_eta = 0.0
    out = case.shell.localized_mode(
        case.phase_contract,
        R,
        s_Q,
        phase,
        n,
        t,
        D_r_C,
        D_z_C,
        _eta(R, Z),
        D_r_eta,
        D_z_eta,
    )
    return _cyl_to_cart(out["velocity"], theta)


def _fd4_first(
    fn: Callable[[np.ndarray], np.ndarray],
    xyz: np.ndarray,
    axis: int,
    step: float,
) -> np.ndarray:
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


def _cartesian_curl_fd4(
    case: SourceAuditCase, xyz: np.ndarray, step: float
) -> np.ndarray:
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
    case: SourceAuditCase,
    s_Q: float,
    xyz: np.ndarray,
    step: float,
    *,
    omit_support_gradients: bool = False,
) -> complex:
    diagonal: list[complex] = []
    for axis in range(3):
        derivative = _fd4_first(
            lambda point: _public_cartesian_velocity(
                case,
                s_Q,
                point,
                omit_support_gradients=omit_support_gradients,
            ),
            xyz,
            axis,
            step,
        )
        diagonal.append(complex(derivative[axis]))
    return complex(sum(diagonal))


def _sample_probes() -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    probes: list[dict[str, Any]] = []
    for case in CASES:
        for region, count, x_ratio_bounds in (
            ("source_shell_interior", 6, (1.8, 4.0)),
            ("near_inner_shell_edge", 6, (1.15, 1.50)),
        ):
            for _ in range(count):
                s_Q = float(rng.uniform(0.58, 1.75))
                X = case.shell.X_a * float(rng.uniform(*x_ratio_bounds))
                R = sqrt(2.0 * s_Q * X)
                theta = float(rng.uniform(0.27, 2.83))
                z_phys = float(rng.uniform(-0.45, 0.45))
                xyz = np.asarray(
                    (R * np.cos(theta), R * np.sin(theta), z_phys), dtype=float
                )
                probes.append(
                    {
                        "case": case,
                        "region": region,
                        "s_Q": s_Q,
                        "xyz": xyz,
                    }
                )
    return probes


def _rms(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _check_fd4_stencils(probes: list[dict[str, Any]]) -> dict[str, float | bool]:
    min_X_ratio = float("inf")
    min_axis_ratio = float("inf")
    for probe in probes:
        case = probe["case"]
        s_Q = probe["s_Q"]
        xyz = np.asarray(probe["xyz"], dtype=float)
        for step in FD4_STEPS:
            for axis in range(3):
                for multiplier in (-2.0, -1.0, 0.0, 1.0, 2.0):
                    point = xyz.copy()
                    point[axis] += multiplier * step
                    R, _, _ = _xyz_to_source(case, point)
                    shell = case.shell.source_shell_coordinates(R, s_Q)
                    X = float(np.asarray(shell["X"]))
                    min_X_ratio = min(min_X_ratio, X / case.shell.X_a)
                    min_axis_ratio = min(
                        min_axis_ratio, R / case.shell.axis_lower_bound
                    )
    return {
        "all_fd4_stencil_points_source_shell_valid": True,
        "minimum_X_over_X_a_across_stencils": min_X_ratio,
        "minimum_R_over_axis_lower_bound_across_stencils": min_axis_ratio,
    }


def _evaluate_step(probes: list[dict[str, Any]], step: float) -> dict[str, Any]:
    curl_errors: list[float] = []
    curl_reference_norms: list[float] = []
    point_relative_errors: list[float] = []
    divergences: list[float] = []
    region_metrics: dict[str, dict[str, list[float]]] = {}
    case_metrics: dict[str, dict[str, list[float]]] = {}

    for probe in probes:
        case = probe["case"]
        s_Q = probe["s_Q"]
        xyz = probe["xyz"]
        oracle = _cartesian_curl_fd4(case, xyz, step)
        public = _public_cartesian_velocity(case, s_Q, xyz)
        difference = float(np.linalg.norm(public - oracle))
        reference = float(np.linalg.norm(oracle))
        divergence = float(abs(_cartesian_divergence_fd4(case, s_Q, xyz, step)))
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
        s_Q = probe["s_Q"]
        xyz = probe["xyz"]
        oracle = _cartesian_curl_fd4(case, xyz, step)
        mutated = _public_cartesian_velocity(
            case, s_Q, xyz, omit_support_gradients=True
        )
        relative_errors.append(
            float(np.linalg.norm(mutated - oracle))
            / max(float(np.linalg.norm(oracle)), 1.0e-300)
        )
        mutated_divergences.append(
            float(
                abs(
                    _cartesian_divergence_fd4(
                        case,
                        s_Q,
                        xyz,
                        step,
                        omit_support_gradients=True,
                    )
                )
            )
        )
        correct_divergences.append(
            float(abs(_cartesian_divergence_fd4(case, s_Q, xyz, step)))
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


def _negative_controls() -> dict[str, bool]:
    off_grid_rejected = True
    below_shell_rejected = True
    for case in CASES:
        try:
            KokunoOscillatoryPhaseContract(
                p=case.p + 0.07,
                p_z=case.p_z,
                x0=case.x0,
                epsilon=case.epsilon,
                m=case.m,
            )
            off_grid_rejected = False
        except ValueError:
            pass
        try:
            R = 0.98 * case.shell.axis_lower_bound
            case.shell.source_shell_coordinates(R, 0.5)
            below_shell_rejected = False
        except ValueError:
            pass
    return {
        "off_source_angular_grid_phase_rejected": off_grid_rejected,
        "below_inner_source_wave_shell_rejected": below_shell_rejected,
    }


def _ratio(coarse: float, fine: float) -> float:
    return float(coarse / max(fine, 1.0e-300))


def build_report() -> dict[str, Any]:
    probes = _sample_probes()
    stencil = _check_fd4_stencils(probes)
    ladder = [_evaluate_step(probes, step) for step in FD4_STEPS]
    finest = ladder[-1]
    curl_ratios = [
        _ratio(
            ladder[index]["curl_relative_rms"],
            ladder[index + 1]["curl_relative_rms"],
        )
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
    negative_controls = _negative_controls()

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
        "all_fd4_stencil_points_source_shell_valid": bool(
            stencil["all_fd4_stencil_points_source_shell_valid"]
        ),
        "off_source_angular_grid_phase_rejected": bool(
            negative_controls["off_source_angular_grid_phase_rejected"]
        ),
        "below_inner_source_wave_shell_rejected": bool(
            negative_controls["below_inner_source_wave_shell_rejected"]
        ),
    }
    structural_pass = bool(all(checks.values()))

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent2_parent_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "source_wave_shell_pr": SOURCE_WAVE_SHELL_PR,
            "source_wave_shell_head": SOURCE_WAVE_SHELL_HEAD,
            "public_artifact": "KokunoSourceWaveShellLocalizedCurl",
        },
        "seed": SEED,
        "parameter_cases": [
            {
                "label": case.label,
                "epsilon": case.epsilon,
                "k": case.k,
                "m": case.m,
                "j": case.j,
                "p": case.p,
                "physical_azimuthal_mode_mj": case.azimuthal_mode,
                "Lambda": case.Lambda,
                "X_a": case.shell.X_a,
                "axis_lower_bound": case.shell.axis_lower_bound,
            }
            for case in CASES
        ],
        "operator": {
            "oracle": (
                "independently recoded source-form phase/covector and localized Cartesian "
                "vector potential, differentiated by centered Cartesian FD4"
            ),
            "public_path": (
                "KokunoSourceWaveShellLocalizedCurl.localized_mode output only, converted "
                "from cylindrical to Cartesian"
            ),
            "caller_coefficient_derivatives": (
                "independently recoded C_m with centered FD6; not the Cartesian "
                "validation derivative"
            ),
            "fd4_steps": list(FD4_STEPS),
            "coefficient_derivative_step": COEFFICIENT_DERIVATIVE_STEP,
            "agent2_cylindrical_validation_oracle_reused": False,
            "previous_agent4_generic_p_phase_reused": False,
            "training_tensor_or_loss_read": False,
            "construction_residual_operator_reused": False,
            "free_forcing_used": False,
        },
        "sampling": {
            "total_points": len(probes),
            "regions": {
                "source_shell_interior": sum(
                    probe["region"] == "source_shell_interior" for probe in probes
                ),
                "near_inner_shell_edge": sum(
                    probe["region"] == "near_inner_shell_edge" for probe in probes
                ),
            },
            "fresh_seed_relative_to_prior_agent4_audits": True,
            "stencil_source_domain": stencil,
        },
        "manufactured_partition_family": {
            "eta_1": "cos(alpha(R,Z))",
            "eta_2": "sin(alpha(R,Z))",
            "sum_eta_beta_squared_exact_algebraically": True,
            "concrete_source_chi_bumps_claimed": False,
        },
        "resolution_ladder": ladder,
        "refinement": {
            "curl_relative_rms_coarse_to_fine_ratios": curl_ratios,
            "divergence_rms_coarse_to_fine_ratios": divergence_ratios,
        },
        "mutation": mutation,
        "negative_controls": negative_controls,
        "local_guards": {
            "thresholds": {
                "finest_curl_relative_rms_max": LOCAL_CURL_FINE_REL_GUARD,
                "finest_divergence_max_abs_max": LOCAL_DIVERGENCE_FINE_ABS_GUARD,
                "minimum_refinement_ratio": LOCAL_REFINEMENT_RATIO_GUARD,
                "missing_gradient_relative_error_min": MISSING_GRADIENT_RELATIVE_ERROR_FLOOR,
                "missing_gradient_divergence_max_abs_min": MISSING_GRADIENT_DIVERGENCE_FLOOR,
            },
            "same_thresholds_as_agent4_prs_353_360": True,
            "checks": checks,
            "source_wave_shell_independent_cartesian_preflight_passed": structural_pass,
            "scope": (
                "local source-grid/source-wave-shell support-localized curl implementation "
                "consistency only"
            ),
        },
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "reason": (
                "global matched Kokuno leading pressure, actual-source public xyz,t "
                "oscillatory velocity, and promoted after-correction composite are absent"
            ),
        },
        "truth_boundary": {
            "actual_source_Lambda_recovered": False,
            "concrete_source_partition_bumps_reconstructed": False,
            "actual_source_background_path_instantiated": False,
            "actual_source_public_oscillatory_xyz_t_velocity_available": False,
            "global_leading_velocity_pressure_available": False,
            "genuinely_independent_second_covariance_column_available": False,
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
        default="artifacts/kokuno_agent4/source_wave_shell_cartesian_audit.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    print(path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
