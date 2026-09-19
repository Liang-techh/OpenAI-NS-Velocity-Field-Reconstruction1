"""Public vector potential for the frozen admitted Kokuno A2 oscillatory field.

The admitted velocity remains the Agent-2 #561 provider
``kokuno_public_z_pullback_velocity.velocity_osc``.  This module does not retune
or rewrite that field.  It exposes the same frozen localized complete-curl
construction one level earlier, as a public-coordinate vector potential whose
ordinary Cartesian curl is intended to reproduce the admitted velocity.

The corrected 2026-09-09 Kokuno reconstruction displays the localized harmonic
potential

    A_beta = eta_beta C_m exp(i k_m Phi),

before taking the complete cylindrical curl.  The frozen repository candidate
already materializes the source phase/frame, signed covariance amplitudes,
localized coefficient ``C_m``, real m=+/-1 pairing, and per-band Q scaling.
For the public pulled-back coordinate chart used by #561 we therefore expose

    A_public,beta,sigma
      = Q_beta^(-A) * 2 Re[
          eta_beta sqrt(epsilon_beta) a_sigma C_{+,beta,sigma}
          exp(i k_beta Phi_beta,sigma)
        ].

The ``Q^(-A)`` here is deliberately a *public-coordinate gauge realization*:
its ordinary curl has the same per-band scaling as the already admitted public
velocity.  It is not relabelled as Kokuno's source ``A_phys`` scaling, it is not
recovered hidden data, and it is not paper-exact.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-PUBLIC-VECTOR-POTENTIAL-042"
SCHEMA = "kokuno-a2-public-oscillatory-vector-potential-v1"
PARENT_AGENT2_PR = 579
PARENT_AGENT2_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
FD4_STEPS = (0.02, 0.01, 0.005)

PotentialEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _rotate_cylindrical_to_cartesian(values: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Rotate (...,r,theta,z) vectors into (...,x,y,z) on sample axes."""
    v = np.asarray(values, dtype=float)
    angle = np.asarray(theta, dtype=float)
    if v.ndim < 2 or v.shape[-1] != 3:
        raise ValueError("values must have trailing vector dimension 3")
    if angle.ndim != 1 or v.shape[0] != angle.size:
        raise ValueError("theta must provide one angle per leading sample")
    extra = (1,) * (v.ndim - 2)
    c = np.cos(angle).reshape((angle.size,) + extra)
    s = np.sin(angle).reshape((angle.size,) + extra)
    return np.stack(
        (
            v[..., 0] * c - v[..., 1] * s,
            v[..., 0] * s + v[..., 1] * c,
            v[..., 2],
        ),
        axis=-1,
    )


def _inside_vector_potential(field: Any, core: dict[str, Any], theta: np.ndarray) -> dict[str, np.ndarray]:
    """Reconstruct the real signed public potential from frozen complete-curl data."""
    frame = core["source_phase_frame"]
    phase = np.asarray(frame["phase"], dtype=float)
    k = np.asarray(frame["k_by_beta"], dtype=float)
    amplitudes = np.asarray(core["reference_amplitudes"], dtype=float)
    C_plus = np.asarray(core["candidate_mode_C_plus_prototype"], dtype=np.complex128)
    Q = np.asarray(core["Q_by_beta"], dtype=float)
    epsilon = np.asarray(core["epsilon_by_beta"], dtype=float)
    exponent = float(core["A"])
    order = tuple(int(j) for j in core["canonical_beta_indices"])
    eta = np.asarray(field._static["geometry"]["eta"], dtype=float)  # frozen A2 realization

    n_sample = theta.size
    n_beta = len(field.beta_labels)
    expected_sign = (n_sample, n_beta, 2)
    if phase.shape != expected_sign or amplitudes.shape != expected_sign:
        raise RuntimeError("frozen phase/amplitude shapes changed")
    if C_plus.shape != expected_sign + (3,):
        raise RuntimeError("frozen localized coefficient shape changed")
    if k.shape != (n_sample, n_beta):
        raise RuntimeError("frozen carrier shape changed")
    if Q.shape != (n_beta,) or epsilon.shape != (n_beta,) or eta.shape != (n_beta,):
        raise RuntimeError("frozen band geometry shape changed")
    if np.any(Q <= 0.0) or np.any(epsilon <= 0.0):
        raise RuntimeError("frozen band scales must stay positive")

    # m=-1 is the exact conjugate reality partner of the m=+1 potential.
    # C_plus already contains the autonomous radial/axial vector-potential support.
    prefactor = (
        eta[None, :, None]
        * np.sqrt(epsilon)[None, :, None]
        * amplitudes
        * Q[None, :, None] ** (-exponent)
    )
    oscillation = np.exp(1j * k[..., None] * phase)
    cylindrical_by_beta_sign = 2.0 * np.real(
        prefactor[..., None] * C_plus * oscillation[..., None]
    )
    cartesian_by_beta_sign = _rotate_cylindrical_to_cartesian(
        cylindrical_by_beta_sign, theta
    )
    cylindrical_by_beta = np.sum(cylindrical_by_beta_sign, axis=-2)
    cartesian_by_beta = np.sum(cartesian_by_beta_sign, axis=-2)
    cylindrical_total = np.sum(np.take(cylindrical_by_beta, order, axis=-2), axis=-2)
    cartesian_total = np.sum(np.take(cartesian_by_beta, order, axis=-2), axis=-2)
    return {
        "vector_potential_cylindrical_by_beta_sign": cylindrical_by_beta_sign,
        "vector_potential_cartesian_by_beta_sign": cartesian_by_beta_sign,
        "vector_potential_cylindrical_by_beta": cylindrical_by_beta,
        "vector_potential_cartesian_by_beta": cartesian_by_beta,
        "vector_potential_cylindrical_total": cylindrical_total,
        "vector_potential_cartesian_total": cartesian_total,
    }


def evaluate_vector_potential_osc(x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
    """Evaluate the frozen public oscillatory vector potential.

    The public support mask is identical to the admitted velocity provider.  No
    numerical differentiation is used here; the potential is reconstructed from
    the already materialized localized complete-curl coefficient and signed
    covariance amplitudes of the frozen candidate.
    """
    field = default_field()
    x, y, z, t = np.broadcast_arrays(
        _finite(x, "x"), _finite(y, "y"), _finite(z, "z"), _finite(t, "t")
    )
    if np.any(t < field.time_min) or np.any(t > field.time_max):
        raise ValueError("t is outside the registered candidate interval")

    sample_shape = x.shape
    n_beta = len(field.beta_labels)
    total = np.zeros(sample_shape + (3,), dtype=float)
    by_beta = np.zeros(sample_shape + (n_beta, 3), dtype=float)
    by_beta_sign = np.zeros(sample_shape + (n_beta, 2, 3), dtype=float)
    R = np.hypot(x, y)
    inside = (
        (R > field.radial_inner)
        & (R < field.radial_outer)
        & (z > field.axial_lower)
        & (z < field.axial_upper)
    )

    indices = np.flatnonzero(inside.ravel())
    if indices.size:
        xf = x.ravel()[indices]
        yf = y.ravel()[indices]
        zf = z.ravel()[indices]
        tf = t.ravel()[indices]
        Rf = np.hypot(xf, yf)
        thetaf = np.arctan2(yf, xf)
        core = field._inside_family(Rf, thetaf, zf, tf)
        potential = _inside_vector_potential(field, core, thetaf)
        total.reshape((-1, 3))[indices] = potential["vector_potential_cartesian_total"]
        by_beta.reshape((-1, n_beta, 3))[indices] = potential[
            "vector_potential_cartesian_by_beta"
        ]
        by_beta_sign.reshape((-1, n_beta, 2, 3))[indices] = potential[
            "vector_potential_cartesian_by_beta_sign"
        ]

    return {
        "vector_potential_cartesian_total": total,
        "vector_potential_cartesian_by_beta": by_beta,
        "vector_potential_cartesian_by_beta_sign": by_beta_sign,
        "support_mask": inside,
        "beta_labels": field.beta_labels,
        "sign_labels": ("sigma_plus", "sigma_minus"),
        "coordinate_contract": (
            "public R=hypot(x,y), theta=atan2(y,x), z; source Z_beta=epsilon_beta*z"
        ),
        "potential_contract": (
            "public-coordinate gauge potential reconstructed from the frozen localized "
            "m=+/-1 complete-curl coefficient; ordinary curl_xyz targets admitted velocity_osc"
        ),
        "source_localized_potential_formula_reused": True,
        "velocity_candidate_changed": False,
        "source_formula_changed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def vector_potential_osc(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return the public-coordinate oscillatory vector potential as ``[...,3]``."""
    return evaluate_vector_potential_osc(x, y, z, t)["vector_potential_cartesian_total"]


def _fd4_axis_derivative(
    evaluator: PotentialEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
    axis: int,
) -> np.ndarray:
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]

    def shifted(multiplier: float) -> np.ndarray:
        args = [value.copy() for value in coords]
        args[axis] = args[axis] + multiplier * step
        return np.asarray(evaluator(args[0], args[1], args[2], t), dtype=float)

    return (
        shifted(-2.0) - 8.0 * shifted(-1.0)
        + 8.0 * shifted(1.0) - shifted(2.0)
    ) / (12.0 * step)


def _fd4_curl(
    evaluator: PotentialEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    d_dx = _fd4_axis_derivative(evaluator, x, y, z, t, step, 0)
    d_dy = _fd4_axis_derivative(evaluator, x, y, z, t, step, 1)
    d_dz = _fd4_axis_derivative(evaluator, x, y, z, t, step, 2)
    return np.stack(
        (
            d_dy[..., 2] - d_dz[..., 1],
            d_dz[..., 0] - d_dx[..., 2],
            d_dx[..., 1] - d_dy[..., 0],
        ),
        axis=-1,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Deterministic support-interior cloud, disjoint from prior A2/A4 audits."""
    radii = (0.48, 0.72, 0.96)
    z_values = (-0.72, 0.58)
    times = (0.34, 0.50, 0.66)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    index = 0
    for time in times:
        for z_value in z_values:
            for radius in radii:
                theta = 0.173 + (index + 1) * golden
                rows.append(
                    (
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        z_value,
                        time,
                    )
                )
                index += 1
    array = np.asarray(rows, dtype=float)
    return array[:, 0], array[:, 1], array[:, 2], array[:, 3]


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def verification_receipt() -> dict[str, Any]:
    """Build a three-resolution public ``curl(A)`` versus ``velocity_osc`` check."""
    x, y, z, t = _verification_cloud()
    target = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    target_rms = _vector_rms(target)
    target_max = float(np.max(np.linalg.norm(target, axis=-1)))
    potential = np.asarray(vector_potential_osc(x, y, z, t), dtype=float)
    potential_rms = _vector_rms(potential)
    rows = []
    for step in FD4_STEPS:
        numerical = _fd4_curl(vector_potential_osc, x, y, z, t, step)
        error = numerical - target
        error_rms = _vector_rms(error)
        error_max = float(np.max(np.linalg.norm(error, axis=-1)))
        rows.append(
            {
                "step": step,
                "absolute_rms_error": error_rms,
                "absolute_max_error": error_max,
                "relative_rms_error": error_rms / max(target_rms, 1.0e-300),
                "relative_max_error": error_max / max(target_max, 1.0e-300),
            }
        )
    refinement = [
        rows[j]["relative_rms_error"] / max(rows[j + 1]["relative_rms_error"], 1.0e-300)
        for j in range(len(rows) - 1)
    ]

    exterior = np.asarray(
        vector_potential_osc(
            np.asarray((0.0, 1.8, 0.8)),
            np.asarray((0.0, 0.0, 0.0)),
            np.asarray((0.0, 0.0, 2.1)),
            np.asarray((0.5, 0.5, 0.5)),
        ),
        dtype=float,
    )
    exterior_max = float(np.max(np.abs(exterior)))

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "velocity_provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
        "vector_potential_provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential:vector_potential_osc"
        ),
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "fd4_steps": list(FD4_STEPS),
        "velocity_rms": target_rms,
        "velocity_sampled_max": target_max,
        "vector_potential_rms": potential_rms,
        "curl_comparison": rows,
        "rms_refinement_ratios": refinement,
        "support_exterior_absolute_max": exterior_max,
        "provenance": {
            "source_formula": "A_beta=eta_beta*C_m*exp(i*k_m*Phi), complete curl before real-pair/Q aggregation",
            "public_realization": (
                "repository-autonomous frozen #561 phase/background/modes/support with public-z pullback; "
                "Q^(-A) potential scaling chosen in public coordinates to match the admitted velocity scaling"
            ),
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "source_formula_changed": False,
            "public_vector_potential_materialized": True,
            "three_resolution_fd4_curl_self_check_performed": True,
            "independent_agent4_vector_potential_audit_required": True,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
