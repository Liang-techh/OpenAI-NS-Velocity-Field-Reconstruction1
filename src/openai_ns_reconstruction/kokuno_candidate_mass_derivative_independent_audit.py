"""Independent Agent-4 audit of candidate h_sigma derivatives in signed curls.

This validator consumes the public candidate-mass -> phase/complete-curl bridge
from Agent 2, but it does not use the bridge's directional-derivative helper as
an oracle.  It differentiates an explicitly declared nonlinear candidate pulse
model independently, checks the public reference-amplitude derivative with a
separate FD6 value oracle, and reconstructs the derivative-only complete-curl
velocity contribution directly from the displayed coefficient/remainder
formulas.

The audit is local structural validation only.  Background/mode prototypes are
still candidate supplied and no global pressure/forcing or held-out NS gate is
assessed here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)
from .kokuno_candidate_mass_bound_phase_curl import (
    KokunoCandidateMassBoundPhaseCurlFamily,
)
from .kokuno_source_compatible_partition import (
    KokunoSourceCompatiblePartitionRealization,
)
from .kokuno_source_signed_covariance_mass import KokunoSourceSignedCovarianceMass

SEED = 9173231
PULSE_RESOLUTIONS = (33, 65, 129)
FD6_STEPS = (0.12, 0.06, 0.03)
H_VALUES = (0.0025, 0.0065, 0.0090)

# Frozen before execution.  These are local structural tolerances, not NS gates.
GUARDS = {
    "reference_algebra_relative_max": 2.0e-11,
    "fd6_finest_relative_rms": 5.0e-8,
    "fd6_min_refinement_ratio": 10.0,
    "scaled_coefficient_relative_max": 2.0e-11,
    "physical_derivative_delta_relative_rms": 5.0e-10,
    "physical_derivative_delta_relative_max": 2.0e-9,
    "signflip_mutation_relative_rms_min": 1.5,
    "minimum_by_sign_delta_norm": 1.0e-12,
}


def _manual_trapezoid(grid: np.ndarray, values: np.ndarray) -> float:
    widths = np.diff(np.asarray(grid, dtype=float))
    vals = np.asarray(values, dtype=float)
    return float(np.sum(widths * (vals[:-1] + vals[1:]) * 0.5))


def _rel_rms(actual: np.ndarray, expected: np.ndarray) -> float:
    a = np.asarray(actual)
    e = np.asarray(expected)
    den = float(np.sum(np.abs(e) ** 2))
    num = float(np.sum(np.abs(a - e) ** 2))
    return math.sqrt(num / max(den, 1.0e-300))


def _rel_max(actual: np.ndarray, expected: np.ndarray) -> float:
    a = np.asarray(actual)
    e = np.asarray(expected)
    return float(np.max(np.abs(a - e), initial=0.0)) / max(
        float(np.max(np.abs(e), initial=0.0)), 1.0e-300
    )


def _geometry(h: float) -> dict[str, Any]:
    partition = KokunoSourceCompatiblePartitionRealization(ell_min=5).evaluate(
        q=np.asarray(2.0**-5.5),
        D_r_q=np.asarray(0.0),
        D_z_q=np.asarray(0.0),
        slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
    )
    return KokunoAutonomousSignedRectangleGeometry(h=h).instantiate(partition)


def _profiles(v: np.ndarray, L_s: float, s_r: float, s_z: float) -> tuple[np.ndarray, ...]:
    s = np.asarray(v, dtype=float) / float(L_s)
    psi = (
        0.58
        + 0.07 * np.sin(np.pi * s) ** 2
        + s_r * (0.05 + 0.02 * s)
        + s_z * (0.03 * s * (1.0 - s))
    )
    dpsi_r = 0.05 + 0.02 * s
    dpsi_z = 0.03 * s * (1.0 - s)

    x_plus = 1.40 + 0.25 * s + s_r * (0.08 + 0.03 * s) + s_z * 0.04 * s * s
    x_minus = (
        2.00
        + 0.18 * (1.0 - s)
        + s_r * 0.05 * (1.0 - s)
        + s_z * (0.07 + 0.02 * s)
    )
    dx_r = np.stack((0.08 + 0.03 * s, 0.05 * (1.0 - s)), axis=-1)
    dx_z = np.stack((0.04 * s * s, 0.07 + 0.02 * s), axis=-1)
    x = np.stack((x_plus, x_minus), axis=-1)
    return psi, x, dpsi_r, dpsi_z, dx_r, dx_z


def _materialized_mass(
    geometry: dict[str, Any], h: float, resolution: int, s_r: float, s_z: float
) -> dict[str, Any]:
    v_by_beta: list[np.ndarray] = []
    psi_by_beta: list[np.ndarray] = []
    x_by_beta_sign: list[np.ndarray] = []
    for L_s in np.asarray(geometry["L_s_by_beta"], dtype=float):
        v = np.linspace(0.0, float(L_s), int(resolution))
        psi, x, *_ = _profiles(v, float(L_s), s_r, s_z)
        v_by_beta.append(v)
        psi_by_beta.append(psi)
        x_by_beta_sign.append(x)
    return KokunoSourceSignedCovarianceMass(h=h).materialize_from_autonomous_geometry(
        geometry,
        v_by_beta=tuple(v_by_beta),
        psi_by_beta=tuple(psi_by_beta),
        x_by_beta_sign=tuple(x_by_beta_sign),
        xi=np.linspace(-1.0, 1.0, 257),
        chi_g=np.ones(257),
        pulse_binding="candidate_derived_from_executable_background",
    )


def _independent_mass_derivatives(
    geometry: dict[str, Any], resolution: int, s_r: float, s_z: float
) -> tuple[np.ndarray, np.ndarray]:
    """Differentiate the discrete pulse law without Agent-2 derivative helpers."""
    Ls = np.asarray(geometry["L_s_by_beta"], dtype=float)
    r0 = float(geometry["r0"])
    out_r = np.empty((len(Ls), 2), dtype=float)
    out_z = np.empty((len(Ls), 2), dtype=float)
    b_g = math.sqrt(2.0) - 1.0
    # chi_g == 1 on [-1,1], so the source displayed transverse factor is exact.
    D_g = 1.0 + b_g * b_g
    for j, L_s in enumerate(Ls):
        v = np.linspace(0.0, float(L_s), int(resolution))
        psi, x, dpsi_r, dpsi_z, dx_r, dx_z = _profiles(v, float(L_s), s_r, s_z)
        c_i = 2.0 * r0 / float(L_s)
        for sign in range(2):
            df_r = (
                2.0 * psi * dpsi_r * x[:, sign] ** 2
                + 2.0 * psi**2 * x[:, sign] * dx_r[:, sign]
            )
            df_z = (
                2.0 * psi * dpsi_z * x[:, sign] ** 2
                + 2.0 * psi**2 * x[:, sign] * dx_z[:, sign]
            )
            out_r[j, sign] = D_g * c_i * _manual_trapezoid(v, df_r)
            out_z[j, sign] = D_g * c_i * _manual_trapezoid(v, df_z)
    return out_r, out_z


def _phase_inputs(
    bridge: KokunoCandidateMassBoundPhaseCurlFamily,
    geometry: dict[str, Any],
    *,
    R_base: float,
    Z_base: float,
    theta: float,
) -> dict[str, Any]:
    labels = tuple(geometry["beta_labels"])
    n = len(labels)
    L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
    R = R_base + np.linspace(0.0, 0.012, n)
    Z = Z_base + np.linspace(-0.008, 0.008, n)
    u_star = np.linspace(2.0, 2.2, n)
    data: dict[str, Any] = {
        "R": R,
        "Z": Z,
        "theta": float(theta),
        "v": np.stack((0.27 * L_s, 0.73 * L_s), axis=-1),
        "beta_labels": labels,
        "R0": np.full(n, 0.82),
        "F0": np.full(n, 1.2),
        "a": np.full(n, 3.2),
        "b_s": np.full(n, 0.2),
        "u_star": u_star,
        "L_s": L_s,
        "F": np.full(n, 1.1),
        "G": np.full(n, 0.24),
        "F_R": np.full(n, 0.08),
        "G_R": np.full(n, -0.05),
        "F_Z": np.full(n, 0.04),
        "G_Z": np.full(n, 0.03),
    }
    frame = bridge.phase_family.phase_frame(**data)
    n_phi = np.asarray(frame["n_phi"], dtype=float)
    base = np.stack(
        (n_phi[..., 1], -n_phi[..., 0], np.zeros_like(n_phi[..., 0])), axis=-1
    )
    t_plus = (1.0 + 0.17j) * base.astype(np.complex128)
    A_c = -np.asarray(frame["c0_by_beta"], dtype=float) * np.sqrt(1.0 + u_star * u_star)
    eta = np.full(n, 1.0 / np.sqrt(float(n)))
    data.update(
        {
            "t_plus_prototype": t_plus,
            "D_r_C_plus_prototype": np.zeros_like(t_plus),
            "D_z_C_plus_prototype": np.zeros_like(t_plus),
            "eta": eta,
            "D_r_eta": np.zeros(n),
            "D_z_eta": np.zeros(n),
            "A_c": A_c,
            "T_N": -0.40 * A_c,
            "T_K": 0.05 * u_star,
        }
    )
    return data


def _bridge_call(
    bridge: KokunoCandidateMassBoundPhaseCurlFamily,
    mass: dict[str, Any],
    inputs: dict[str, Any],
    D_r_h: np.ndarray,
    D_z_h: np.ndarray,
) -> dict[str, Any]:
    return bridge.physical_family_from_materialized_mass(
        mass,
        D_r_h_by_beta_sign=D_r_h,
        D_z_h_by_beta_sign=D_z_h,
        h_derivative_binding="candidate_derived_from_same_pulse_model",
        **inputs,
    )


def _fd6(values: dict[int, np.ndarray], step: float) -> np.ndarray:
    return (
        -values[-3]
        + 9.0 * values[-2]
        - 45.0 * values[-1]
        + 45.0 * values[1]
        - 9.0 * values[2]
        + values[3]
    ) / (60.0 * step)


def _amplitudes_at(
    bridge: KokunoCandidateMassBoundPhaseCurlFamily,
    geometry: dict[str, Any],
    inputs: dict[str, Any],
    h: float,
    resolution: int,
    s_r: float,
    s_z: float,
) -> np.ndarray:
    mass = _materialized_mass(geometry, h, resolution, s_r, s_z)
    zeros = np.zeros_like(np.asarray(mass["h_sigma_by_beta_sign"], dtype=float))
    return np.asarray(_bridge_call(bridge, mass, inputs, zeros, zeros)["reference_amplitudes"])


def _fd6_direction(
    bridge: KokunoCandidateMassBoundPhaseCurlFamily,
    geometry: dict[str, Any],
    inputs: dict[str, Any],
    h: float,
    resolution: int,
    s_r: float,
    s_z: float,
    step: float,
    direction: str,
) -> np.ndarray:
    values: dict[int, np.ndarray] = {}
    for k in (-3, -2, -1, 1, 2, 3):
        rr = s_r + k * step if direction == "r" else s_r
        zz = s_z + k * step if direction == "z" else s_z
        values[k] = _amplitudes_at(bridge, geometry, inputs, h, resolution, rr, zz)
    return _fd6(values, step)


def _independent_coefficient(
    n_phi: np.ndarray, t0: np.ndarray, epsilon_by_beta: np.ndarray
) -> np.ndarray:
    out = np.empty_like(np.asarray(t0, dtype=np.complex128))
    for j, eps in enumerate(np.asarray(epsilon_by_beta, dtype=float)):
        k = int(math.ceil(float(eps) ** -0.5))
        for sign in range(2):
            n = np.asarray(n_phi[j, sign], dtype=float)
            t = np.asarray(t0[j, sign], dtype=np.complex128)
            out[j, sign] = 1j * np.cross(n, t) / (k * float(np.dot(n, n)))
    return out


def _independent_velocity_delta(
    *,
    inputs: dict[str, Any],
    frame: dict[str, Any],
    epsilon_by_beta: np.ndarray,
    Q_by_beta: np.ndarray,
    A: float,
    delta_DrC: np.ndarray,
    delta_DzC: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    labels = tuple(inputs["beta_labels"])
    eta = np.asarray(inputs["eta"], dtype=float)
    phase = np.asarray(frame["phase"], dtype=float)
    theta = float(inputs["theta"])
    by_sign_cyl = np.empty(delta_DrC.shape, dtype=float)
    for j, eps in enumerate(np.asarray(epsilon_by_beta, dtype=float)):
        k = int(math.ceil(float(eps) ** -0.5))
        qscale = float(Q_by_beta[j]) ** (-float(A))
        for sign in range(2):
            dr = np.asarray(delta_DrC[j, sign], dtype=np.complex128)
            dz = np.asarray(delta_DzC[j, sign], dtype=np.complex128)
            rem = eta[j] * np.asarray(
                (-dz[1], dz[0] - dr[2], dr[1]), dtype=np.complex128
            )
            by_sign_cyl[j, sign] = qscale * 2.0 * np.real(
                rem * np.exp(1j * k * phase[j, sign])
            )
    c = math.cos(theta)
    s = math.sin(theta)
    by_sign_cart = np.empty_like(by_sign_cyl)
    by_sign_cart[..., 0] = by_sign_cyl[..., 0] * c - by_sign_cyl[..., 1] * s
    by_sign_cart[..., 1] = by_sign_cyl[..., 0] * s + by_sign_cyl[..., 1] * c
    by_sign_cart[..., 2] = by_sign_cyl[..., 2]
    by_beta = np.sum(by_sign_cart, axis=-2)
    order = tuple(sorted(range(len(labels)), key=lambda j: (int(labels[j][0]), repr(labels[j][1]))))
    total = np.sum(np.take(by_beta, order, axis=-2), axis=-2)
    return by_sign_cart, total


def run_audit() -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    fd_errors: list[list[float]] = [[], [], []]
    algebra_relmax: list[float] = []
    coeff_relmax: list[float] = []
    velocity_delta_rms: list[float] = []
    velocity_delta_max: list[float] = []
    mutation_rms: list[float] = []
    by_sign_delta_norms: list[float] = []
    total_delta_norms: list[float] = []
    case_records: list[dict[str, Any]] = []
    axis_near = 0
    off_grid = 0

    case_id = 0
    for h in H_VALUES:
        geometry = _geometry(h)
        bridge = KokunoCandidateMassBoundPhaseCurlFamily(h=h)
        for resolution in PULSE_RESOLUTIONS:
            s_r = float(rng.uniform(-0.08, 0.08))
            s_z = float(rng.uniform(-0.08, 0.08))
            if case_id < 3:
                R_base = 0.055 + 0.007 * case_id
                axis_near += 1
            else:
                R_base = float(rng.uniform(0.43, 0.91))
                off_grid += 1
            Z_base = float(rng.uniform(-0.11, 0.11))
            theta = float(rng.uniform(0.19, 1.17))
            inputs = _phase_inputs(
                bridge, geometry, R_base=R_base, Z_base=Z_base, theta=theta
            )
            mass = _materialized_mass(geometry, h, resolution, s_r, s_z)
            D_r_h, D_z_h = _independent_mass_derivatives(
                geometry, resolution, s_r, s_z
            )
            out = _bridge_call(bridge, mass, inputs, D_r_h, D_z_h)
            amps = np.asarray(out["reference_amplitudes"], dtype=float)
            hmass = np.asarray(mass["h_sigma_by_beta_sign"], dtype=float)
            exact_Dr_a = -0.5 * amps * D_r_h / hmass
            exact_Dz_a = -0.5 * amps * D_z_h / hmass
            pub_Dr_a = np.asarray(out["D_r_reference_amplitudes"], dtype=float)
            pub_Dz_a = np.asarray(out["D_z_reference_amplitudes"], dtype=float)
            algebra_relmax.append(
                max(_rel_max(pub_Dr_a, exact_Dr_a), _rel_max(pub_Dz_a, exact_Dz_a))
            )

            for idx, step in enumerate(FD6_STEPS):
                fd_r = _fd6_direction(
                    bridge, geometry, inputs, h, resolution, s_r, s_z, step, "r"
                )
                fd_z = _fd6_direction(
                    bridge, geometry, inputs, h, resolution, s_r, s_z, step, "z"
                )
                actual = np.concatenate((pub_Dr_a.ravel(), pub_Dz_a.ravel()))
                expected = np.concatenate((fd_r.ravel(), fd_z.ravel()))
                fd_errors[idx].append(_rel_rms(actual, expected))

            frame = out["source_phase_frame"]
            n_phi = np.asarray(frame["n_phi"], dtype=float)
            t0 = np.asarray(inputs["t_plus_prototype"], dtype=np.complex128)
            eps = np.asarray(out["epsilon_by_beta"], dtype=float)
            Q = np.asarray(out["Q_by_beta"], dtype=float)
            C0 = _independent_coefficient(n_phi, t0, eps)
            delta_DrC = np.sqrt(eps)[:, None, None] * exact_Dr_a[..., None] * C0
            delta_DzC = np.sqrt(eps)[:, None, None] * exact_Dz_a[..., None] * C0
            pub_DrC = np.asarray(out["scaled_D_r_C_plus_prototype"])
            pub_DzC = np.asarray(out["scaled_D_z_C_plus_prototype"])
            coeff_relmax.append(
                max(_rel_max(pub_DrC, delta_DrC), _rel_max(pub_DzC, delta_DzC))
            )

            zeros = np.zeros_like(D_r_h)
            zero_out = _bridge_call(bridge, mass, inputs, zeros, zeros)
            observed_by_sign = (
                np.asarray(out["candidate_velocity_osc_cartesian_by_beta_sign"], dtype=float)
                - np.asarray(zero_out["candidate_velocity_osc_cartesian_by_beta_sign"], dtype=float)
            )
            observed_total = (
                np.asarray(out["candidate_velocity_osc_cartesian_total"], dtype=float)
                - np.asarray(zero_out["candidate_velocity_osc_cartesian_total"], dtype=float)
            )
            expected_by_sign, expected_total = _independent_velocity_delta(
                inputs=inputs,
                frame=frame,
                epsilon_by_beta=eps,
                Q_by_beta=Q,
                A=float(out["A"]),
                delta_DrC=delta_DrC,
                delta_DzC=delta_DzC,
            )
            velocity_delta_rms.append(_rel_rms(observed_by_sign, expected_by_sign))
            velocity_delta_max.append(_rel_max(observed_by_sign, expected_by_sign))
            by_sign_delta_norms.append(float(np.linalg.norm(expected_by_sign)))
            total_delta_norms.append(float(np.linalg.norm(expected_total)))

            flipped = _bridge_call(bridge, mass, inputs, -D_r_h, -D_z_h)
            flipped_deriv = np.concatenate(
                (
                    np.asarray(flipped["D_r_reference_amplitudes"]).ravel(),
                    np.asarray(flipped["D_z_reference_amplitudes"]).ravel(),
                )
            )
            correct_deriv = np.concatenate((exact_Dr_a.ravel(), exact_Dz_a.ravel()))
            mutation_rms.append(_rel_rms(flipped_deriv, correct_deriv))

            case_records.append(
                {
                    "case": case_id,
                    "h": h,
                    "pulse_resolution": resolution,
                    "slow_coordinate": [s_r, s_z],
                    "R_base": R_base,
                    "Z_base": Z_base,
                    "axis_near": R_base < 0.1,
                    "reference_algebra_relative_max": algebra_relmax[-1],
                    "scaled_coefficient_relative_max": coeff_relmax[-1],
                    "physical_delta_relative_rms": velocity_delta_rms[-1],
                    "physical_delta_relative_max": velocity_delta_max[-1],
                    "signflip_mutation_relative_rms": mutation_rms[-1],
                    "expected_by_sign_delta_norm": by_sign_delta_norms[-1],
                    "expected_total_delta_norm": total_delta_norms[-1],
                }
            )
            case_id += 1

    fd_rms = [float(math.sqrt(np.mean(np.square(level)))) for level in fd_errors]
    refinement = [fd_rms[0] / fd_rms[1], fd_rms[1] / fd_rms[2]]
    metrics = {
        "case_count": len(case_records),
        "axis_near_case_count": axis_near,
        "off_grid_case_count": off_grid,
        "h_values": list(H_VALUES),
        "pulse_resolutions": list(PULSE_RESOLUTIONS),
        "fd6_steps": list(FD6_STEPS),
        "reference_algebra_relative_max": max(algebra_relmax),
        "fd6_reference_amplitude_relative_rms": fd_rms,
        "fd6_refinement_ratios": refinement,
        "scaled_coefficient_relative_max": max(coeff_relmax),
        "physical_derivative_delta_relative_rms": float(
            math.sqrt(np.mean(np.square(velocity_delta_rms)))
        ),
        "physical_derivative_delta_relative_max": max(velocity_delta_max),
        "minimum_by_sign_delta_norm": min(by_sign_delta_norms),
        "minimum_total_delta_norm": min(total_delta_norms),
        "minimum_signflip_mutation_relative_rms": min(mutation_rms),
    }
    checks = {
        "reference_algebra": metrics["reference_algebra_relative_max"]
        <= GUARDS["reference_algebra_relative_max"],
        "fd6_finest": fd_rms[-1] <= GUARDS["fd6_finest_relative_rms"],
        "fd6_refinement": min(refinement) >= GUARDS["fd6_min_refinement_ratio"],
        "scaled_coefficient": metrics["scaled_coefficient_relative_max"]
        <= GUARDS["scaled_coefficient_relative_max"],
        "physical_derivative_delta_rms": metrics["physical_derivative_delta_relative_rms"]
        <= GUARDS["physical_derivative_delta_relative_rms"],
        "physical_derivative_delta_max": metrics["physical_derivative_delta_relative_max"]
        <= GUARDS["physical_derivative_delta_relative_max"],
        "mutation_detected": metrics["minimum_signflip_mutation_relative_rms"]
        >= GUARDS["signflip_mutation_relative_rms_min"],
        "nontrivial_derivative_delta": metrics["minimum_by_sign_delta_norm"]
        >= GUARDS["minimum_by_sign_delta_norm"],
        "coverage": axis_near >= 3 and off_grid >= 3,
    }
    passed = all(checks.values())
    return {
        "schema": "kokuno-agent4-candidate-mass-derivative-independent-audit-v1",
        "seed": SEED,
        "parent_handoff": "Agent-2 #522 candidate h_sigma -> phase/complete-curl bridge",
        "independence": {
            "public_values_only": True,
            "parent_directional_derivative_helper_used_as_oracle": False,
            "independent_discrete_mass_derivative": True,
            "independent_fd6_reference_amplitude_oracle": True,
            "independent_complete_curl_coefficient_formula": True,
            "independent_derivative_only_velocity_delta_formula": True,
            "training_tensor_or_loss_used": False,
            "pressure_or_forcing_fit_used": False,
        },
        "guards": dict(GUARDS),
        "metrics": metrics,
        "checks": checks,
        "local_structural_preflight_passed": passed,
        "truth_boundary": {
            "candidate_mass_derivative_bridge_independently_audited": passed,
            "candidate_background_and_mode_prototypes_still_caller_supplied": True,
            "public_xyz_t_velocity_osc_provider_available": False,
            "leading_plus_oscillatory_plus_correction_candidate_materialized": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "cases": case_records,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["metrics"], sort_keys=True))
    return 0 if report["local_structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
