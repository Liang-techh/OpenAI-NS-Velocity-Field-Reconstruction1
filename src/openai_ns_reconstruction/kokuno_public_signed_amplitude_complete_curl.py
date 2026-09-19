"""Typed signed-amplitude complete-curl correction for the frozen Kokuno A2 family.

This module is an Agent-2-owned correction kernel.  It does not recompute the
Agent-3 mean/radial defect chain.  Instead it accepts a provenance-bearing
signed radial amplitude differential and realizes it through the *same*
localized vector-potential / complete-curl machinery as the frozen public A2
oscillatory family.

Source-displayed structure reused here:
- signed covariance axis sigma=+/- with g_sigma=sqrt(epsilon)*a_sigma;
- product rule D(g C)=g DC + (Dg) C inside the complete curl;
- coefficient/support localization before curling;
- m=+/-1 conjugate reality pairing;
- per-band Q^(-A) physical scaling.

Repository realizations, not paper-exact data:
- the public-z pullback inherited from A2 #561;
- a C2 compact quintic spline that lifts a radial delta-a table to arbitrary R;
- the same frozen bounded A2 autonomous time modulation, normalized to equal
  the supplied delta-a profile at ``reference_time``;
- the public-coordinate gauge potential inherited from A2 #616.

The correction velocity is assembled as curl(delta A), rather than by fitting
three Cartesian components independently.  No residual, forcing, pressure,
gain, target, or scientific threshold enters this API.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.interpolate import make_interp_spline

from .kokuno_public_candidate_velocity import _pulse_x
from .kokuno_public_oscillatory_time_derivative import target_log_amplitude_derivative
from .kokuno_public_z_pullback_velocity import default_field
from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl

TASK = "KOKUNO-A2-SIGNED-AMPLITUDE-COMPLETE-CURL-049"
SCHEMA = "kokuno-a2-signed-amplitude-complete-curl-v1"
PARENT_AGENT2_PR = 716
PARENT_AGENT2_HEAD = "b8d0a4965a95ffb024ca28b6649b795f72f00ce3"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
PROFILE_ENDPOINT_RTOL = 1.0e-11
CURL_FD4_STEPS = (0.012, 0.006, 0.003)
TIME_FD6_STEPS = (0.012, 0.006, 0.003)

VectorEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


@dataclass(frozen=True)
class SignedRadialAmplitudeCorrection:
    """Static signed ``delta_a_sigma(R)`` profile with provenance.

    The profile is *not* a free oscillation fit.  It is intended as the typed
    handoff target for Agent 3's ``delta_a=(delta_a_plus,delta_a_minus)`` radial
    output.  Scientific source certification remains an explicit external bit.
    """

    radii: tuple[float, ...]
    delta_a_plus: tuple[float, ...]
    delta_a_minus: tuple[float, ...]
    reference_time: float
    producer_kind: str
    provenance: str
    source_mean_amplitude_differential_certified: bool = False

    def __post_init__(self) -> None:
        field = default_field()
        radii = _finite(self.radii, "radii")
        if radii.ndim != 1 or radii.size < 8:
            raise ValueError("radii must be a one-dimensional profile with at least 8 nodes")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        if not (field.radial_inner < radii[0] < radii[-1] < field.radial_outer):
            raise ValueError("radial correction support must stay strictly inside A2 support")

        plus = _finite(self.delta_a_plus, "delta_a_plus")
        minus = _finite(self.delta_a_minus, "delta_a_minus")
        if plus.shape != radii.shape or minus.shape != radii.shape:
            raise ValueError("signed delta-a profiles must match the radial grid exactly")

        scale = max(1.0, float(np.max(np.abs(np.concatenate((plus, minus))))))
        endpoint_tol = PROFILE_ENDPOINT_RTOL * scale
        endpoints = np.asarray((plus[0], plus[-1], minus[0], minus[-1]), dtype=float)
        if np.max(np.abs(endpoints)) > endpoint_tol:
            raise ValueError(
                "compact correction requires signed delta-a endpoints to vanish "
                "within the registered roundoff tolerance"
            )
        # Make the zero extension literal once the caller has passed the
        # preregistered roundoff closure guard.
        plus = plus.copy()
        minus = minus.copy()
        plus[[0, -1]] = 0.0
        minus[[0, -1]] = 0.0

        reference_time = float(self.reference_time)
        if not math.isfinite(reference_time) or not (
            field.time_min <= reference_time <= field.time_max
        ):
            raise ValueError("reference_time must lie in the frozen A2 candidate interval")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not isinstance(self.source_mean_amplitude_differential_certified, bool):
            raise TypeError("source_mean_amplitude_differential_certified must be bool")

        object.__setattr__(self, "radii", tuple(float(v) for v in radii))
        object.__setattr__(self, "delta_a_plus", tuple(float(v) for v in plus))
        object.__setattr__(self, "delta_a_minus", tuple(float(v) for v in minus))
        object.__setattr__(self, "reference_time", reference_time)

    def stacked(self) -> np.ndarray:
        return np.stack(
            (
                np.asarray(self.delta_a_plus, dtype=float),
                np.asarray(self.delta_a_minus, dtype=float),
            ),
            axis=-1,
        )

    def to_receipt(self) -> dict[str, Any]:
        return {
            "radii": list(self.radii),
            "delta_a_plus": list(self.delta_a_plus),
            "delta_a_minus": list(self.delta_a_minus),
            "reference_time": self.reference_time,
            "producer_kind": self.producer_kind,
            "provenance": self.provenance,
            "source_mean_amplitude_differential_certified": (
                self.source_mean_amplitude_differential_certified
            ),
            "radial_interpolant": (
                "repository C2 compact quintic B-spline; endpoint value/first/second "
                "derivatives match the exterior zero extension"
            ),
            "time_lift": (
                "repository-autonomous frozen A2 signed-amplitude modulation normalized "
                "to one at reference_time"
            ),
        }


def profile_from_delta_a(
    radii: Any,
    delta_a: Any,
    *,
    reference_time: float,
    producer_kind: str,
    provenance: str,
    source_mean_amplitude_differential_certified: bool = False,
) -> SignedRadialAmplitudeCorrection:
    """Build the A2 correction witness from an ``(N,2)`` signed delta-a table."""
    rr = _finite(radii, "radii")
    da = _finite(delta_a, "delta_a")
    if rr.ndim != 1 or da.shape != (rr.size, 2):
        raise ValueError("delta_a must have shape (len(radii),2) ordered plus/minus")
    return SignedRadialAmplitudeCorrection(
        radii=tuple(float(v) for v in rr),
        delta_a_plus=tuple(float(v) for v in da[:, 0]),
        delta_a_minus=tuple(float(v) for v in da[:, 1]),
        reference_time=reference_time,
        producer_kind=producer_kind,
        provenance=provenance,
        source_mean_amplitude_differential_certified=(
            source_mean_amplitude_differential_certified
        ),
    )


def _profile_value_derivative(
    correction: SignedRadialAmplitudeCorrection, radius: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the compact C2 quintic lift and its ordinary radial derivative."""
    r = _finite(radius, "radius")
    nodes = np.asarray(correction.radii, dtype=float)
    samples = correction.stacked()
    values = np.zeros(r.shape + (2,), dtype=float)
    derivatives = np.zeros_like(values)
    inside = (r >= nodes[0]) & (r <= nodes[-1])
    if not np.any(inside):
        return values, derivatives

    for sign in range(2):
        spline = make_interp_spline(
            nodes,
            samples[:, sign],
            k=5,
            bc_type=(
                [(1, 0.0), (2, 0.0)],
                [(1, 0.0), (2, 0.0)],
            ),
        )
        values[..., sign][inside] = spline(r[inside])
        derivatives[..., sign][inside] = spline.derivative(1)(r[inside])
    return values, derivatives


def _target_signed_gap(t: Any) -> np.ndarray:
    """Return gaps proportional to a_plus^2 and a_minus^2 for the frozen A2 law."""
    field = default_field()
    tt = _finite(t, "t")
    if np.any(tt < field.time_min) or np.any(tt > field.time_max):
        raise ValueError("t is outside the frozen A2 candidate interval")
    tau = (tt - field.time_min) / (field.time_max - field.time_min)
    angle = 2.0 * math.pi * tau
    normal = field.normal_target_ratio * (
        1.0 + field.time_modulation * np.sin(angle)
    )
    cross = field.cross_target_ratio * (
        1.0 + field.time_modulation * np.cos(angle)
    )
    gaps = np.stack((normal - cross, normal + cross), axis=-1)
    if np.any(gaps <= 0.0):
        raise RuntimeError("frozen A2 signed-amplitude gap lost positivity")
    return gaps


def _temporal_lift(
    t: np.ndarray, reference_time: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return signed amplitude lift and exact time derivative of that lift."""
    gaps = _target_signed_gap(t)
    ref = _target_signed_gap(np.asarray(reference_time, dtype=float))
    lift = np.sqrt(gaps / ref)
    log_derivative = target_log_amplitude_derivative(t)
    return lift, lift * log_derivative


def _rotate_cylindrical_to_cartesian(
    values: np.ndarray, theta: np.ndarray
) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    angle = np.asarray(theta, dtype=float)
    c = np.cos(angle)[..., None, None]
    s = np.sin(angle)[..., None, None]
    return np.stack(
        (
            v[..., 0] * c - v[..., 1] * s,
            v[..., 0] * s + v[..., 1] * c,
            v[..., 2],
        ),
        axis=-1,
    )


def _inside_correction(
    correction: SignedRadialAmplitudeCorrection,
    R: np.ndarray,
    theta: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> dict[str, np.ndarray]:
    """Assemble correction potential and complete curl on interior samples."""
    field = default_field()
    core = field._inside_family(R, theta, z, t)
    frame = core["source_phase_frame"]
    phase = np.asarray(frame["phase"], dtype=float)
    n_phi = np.asarray(frame["n_phi"], dtype=float)
    k = np.asarray(frame["k_by_beta"], dtype=float)
    C0 = np.asarray(core["candidate_mode_C_plus_prototype"], dtype=np.complex128)
    D_r_C0 = np.asarray(
        core["candidate_mode_D_r_C_plus_prototype"], dtype=np.complex128
    )
    D_z_C0 = np.asarray(
        core["candidate_mode_D_z_C_plus_prototype"], dtype=np.complex128
    )
    Q = np.asarray(core["Q_by_beta"], dtype=float)
    epsilon = np.asarray(core["epsilon_by_beta"], dtype=float)
    exponent = float(core["A"])
    order = tuple(int(j) for j in core["canonical_beta_indices"])

    m = R.size
    n_beta = len(field.beta_labels)
    expected = (m, n_beta, 2)
    if phase.shape != expected or n_phi.shape != expected + (3,):
        raise RuntimeError("frozen phase/frame shape changed")
    if C0.shape != expected + (3,):
        raise RuntimeError("frozen localized coefficient shape changed")
    if Q.shape != (n_beta,) or epsilon.shape != (n_beta,):
        raise RuntimeError("frozen band schedule shape changed")

    radial_envelope = np.asarray(core["candidate_radial_envelope"], dtype=float)
    axial_envelope = np.asarray(core["candidate_axial_envelope"], dtype=float)
    envelope = radial_envelope * axial_envelope
    pulse_factors = _pulse_x(np.asarray((0.25, 0.75)))[
        np.arange(2), np.arange(2)
    ]
    alpha = complex(1.0, field.mode_imaginary_ratio)
    base_direction = np.stack(
        (n_phi[..., 1], -n_phi[..., 0], np.zeros_like(n_phi[..., 0])),
        axis=-1,
    ).astype(np.complex128)
    t0 = (
        envelope[:, None, None, None]
        * alpha
        * pulse_factors[None, None, :, None]
        * base_direction
    )

    profile_value, profile_dr = _profile_value_derivative(correction, R)
    lift, lift_dt = _temporal_lift(t, correction.reference_time)
    amplitude = profile_value * lift
    amplitude_dr = profile_dr * lift
    amplitude_dt = profile_value * lift_dt

    geometry = field._static["geometry"]
    eta0 = np.asarray(geometry["eta"], dtype=float)
    D_r_eta0 = np.asarray(geometry["D_r_eta"], dtype=float)
    D_z_eta0 = np.asarray(geometry["D_z_eta"], dtype=float)
    eta = np.broadcast_to(eta0, (m, n_beta))
    D_r_eta = np.broadcast_to(D_r_eta0, (m, n_beta))
    D_z_eta = np.broadcast_to(D_z_eta0, (m, n_beta))

    cyl_by_beta_sign = np.zeros((m, n_beta, 2, 3), dtype=float)
    A_cyl_by_beta_sign = np.zeros_like(cyl_by_beta_sign)
    for j in range(n_beta):
        root_eps = float(np.sqrt(epsilon[j]))
        q_scale = float(Q[j] ** (-exponent))
        for sign in range(2):
            local_plus = KokunoSourceSupportLocalizedCurl(
                epsilon=float(epsilon[j]), m=1, partition_atol=field.bridge.partition_atol
            )
            n = n_phi[:, j, sign, :]
            g = root_eps * amplitude[:, sign]
            D_r_g = root_eps * amplitude_dr[:, sign]
            t_mode = g[:, None] * t0[:, j, sign, :]
            D_r_C = (
                g[:, None] * D_r_C0[:, j, sign, :]
                + D_r_g[:, None] * C0[:, j, sign, :]
            )
            D_z_C = g[:, None] * D_z_C0[:, j, sign, :]
            common = (R, phase[:, j, sign], n)
            plus = local_plus.localized_mode(
                *common,
                t_mode,
                D_r_C,
                D_z_C,
                eta[:, j],
                D_r_eta[:, j],
                D_z_eta[:, j],
            )
            minus = KokunoSourceSupportLocalizedCurl(
                epsilon=float(epsilon[j]), m=-1, partition_atol=field.bridge.partition_atol
            ).localized_mode(
                *common,
                np.conjugate(t_mode),
                np.conjugate(D_r_C),
                np.conjugate(D_z_C),
                eta[:, j],
                D_r_eta[:, j],
                D_z_eta[:, j],
            )
            pair = np.asarray(plus["velocity"]) + np.asarray(minus["velocity"])
            scale = max(1.0, float(np.max(np.abs(pair.real), initial=0.0)))
            if float(np.max(np.abs(pair.imag), initial=0.0)) > (
                field.bridge.phase_family.complete_curl_family.reality_atol * scale
            ):
                raise RuntimeError("correction m=+/-1 pair failed the frozen reality guard")
            cyl_by_beta_sign[:, j, sign, :] = q_scale * pair.real

            oscillation = np.exp(1j * k[:, j] * phase[:, j, sign])
            A_cyl_by_beta_sign[:, j, sign, :] = q_scale * 2.0 * np.real(
                eta[:, j, None]
                * g[:, None]
                * C0[:, j, sign, :]
                * oscillation[:, None]
            )

    cart_by_beta_sign = _rotate_cylindrical_to_cartesian(
        cyl_by_beta_sign, theta
    )
    A_cart_by_beta_sign = _rotate_cylindrical_to_cartesian(
        A_cyl_by_beta_sign, theta
    )
    cart_by_beta = np.sum(cart_by_beta_sign, axis=-2)
    A_cart_by_beta = np.sum(A_cart_by_beta_sign, axis=-2)
    cart_total = np.sum(np.take(cart_by_beta, order, axis=-2), axis=-2)
    A_cart_total = np.sum(np.take(A_cart_by_beta, order, axis=-2), axis=-2)

    # The only candidate time dependence in this correction is the normalized
    # signed amplitude lift; phase/support/basis remain frozen exactly as in #579.
    log_derivative = target_log_amplitude_derivative(t)
    velocity_dt_by_beta_sign = (
        cart_by_beta_sign * log_derivative[:, None, :, None]
    )
    A_dt_by_beta_sign = (
        A_cart_by_beta_sign * log_derivative[:, None, :, None]
    )
    velocity_dt_by_beta = np.sum(velocity_dt_by_beta_sign, axis=-2)
    A_dt_by_beta = np.sum(A_dt_by_beta_sign, axis=-2)
    velocity_dt_total = np.sum(
        np.take(velocity_dt_by_beta, order, axis=-2), axis=-2
    )
    A_dt_total = np.sum(np.take(A_dt_by_beta, order, axis=-2), axis=-2)

    return {
        "velocity_cartesian_total": cart_total,
        "velocity_cartesian_by_beta": cart_by_beta,
        "velocity_cartesian_by_beta_sign": cart_by_beta_sign,
        "vector_potential_cartesian_total": A_cart_total,
        "vector_potential_cartesian_by_beta": A_cart_by_beta,
        "vector_potential_cartesian_by_beta_sign": A_cart_by_beta_sign,
        "velocity_dt_cartesian_total": velocity_dt_total,
        "velocity_dt_cartesian_by_beta": velocity_dt_by_beta,
        "velocity_dt_cartesian_by_beta_sign": velocity_dt_by_beta_sign,
        "vector_potential_dt_cartesian_total": A_dt_total,
        "vector_potential_dt_cartesian_by_beta": A_dt_by_beta,
        "vector_potential_dt_cartesian_by_beta_sign": A_dt_by_beta_sign,
        "profile_amplitude_by_sign": amplitude,
        "profile_amplitude_dr_by_sign": amplitude_dr,
        "profile_amplitude_dt_by_sign": amplitude_dt,
    }


def evaluate_signed_amplitude_correction(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> dict[str, Any]:
    """Evaluate ``delta A``, ``curl(delta A)`` and their candidate time derivatives."""
    if not isinstance(correction, SignedRadialAmplitudeCorrection):
        raise TypeError("correction must be a SignedRadialAmplitudeCorrection")
    field = default_field()
    x, y, z, t = np.broadcast_arrays(
        _finite(x, "x"), _finite(y, "y"), _finite(z, "z"), _finite(t, "t")
    )
    if np.any(t < field.time_min) or np.any(t > field.time_max):
        raise ValueError("t is outside the frozen A2 candidate interval")
    shape = x.shape
    n_beta = len(field.beta_labels)

    keys = (
        "velocity_cartesian_total",
        "vector_potential_cartesian_total",
        "velocity_dt_cartesian_total",
        "vector_potential_dt_cartesian_total",
    )
    totals = {key: np.zeros(shape + (3,), dtype=float) for key in keys}
    by_sign_keys = (
        "velocity_cartesian_by_beta_sign",
        "vector_potential_cartesian_by_beta_sign",
        "velocity_dt_cartesian_by_beta_sign",
        "vector_potential_dt_cartesian_by_beta_sign",
    )
    by_sign = {
        key: np.zeros(shape + (n_beta, 2, 3), dtype=float)
        for key in by_sign_keys
    }

    R = np.hypot(x, y)
    nodes = np.asarray(correction.radii, dtype=float)
    inside = (
        (R > nodes[0])
        & (R < nodes[-1])
        & (R > field.radial_inner)
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
        interior = _inside_correction(correction, Rf, thetaf, zf, tf)
        for key in keys:
            totals[key].reshape((-1, 3))[indices] = interior[key]
        for key in by_sign_keys:
            by_sign[key].reshape((-1, n_beta, 2, 3))[indices] = interior[key]

    return {
        **totals,
        **by_sign,
        "support_mask": inside,
        "beta_labels": field.beta_labels,
        "sign_labels": ("sigma_plus", "sigma_minus"),
        "coordinate_contract": (
            "public R=hypot(x,y), theta=atan2(y,x), z with per-beta "
            "source Z_beta=epsilon_beta*z"
        ),
        "correction_contract": (
            "shared signed radial delta-a profile -> C2 compact lift -> "
            "localized vector potential -> source complete curl; no free Cartesian fit"
        ),
        "source_complete_curl_product_rule_reused": True,
        "source_mean_amplitude_differential_certified": (
            correction.source_mean_amplitude_differential_certified
        ),
        "agent3_mean_radial_chain_reimplemented": False,
        "public_velocity_correction_materialized": True,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def correction_vector_potential(
    correction: SignedRadialAmplitudeCorrection, x: Any, y: Any, z: Any, t: Any
) -> np.ndarray:
    return evaluate_signed_amplitude_correction(correction, x, y, z, t)[
        "vector_potential_cartesian_total"
    ]


def correction_velocity(
    correction: SignedRadialAmplitudeCorrection, x: Any, y: Any, z: Any, t: Any
) -> np.ndarray:
    return evaluate_signed_amplitude_correction(correction, x, y, z, t)[
        "velocity_cartesian_total"
    ]


def correction_velocity_dt(
    correction: SignedRadialAmplitudeCorrection, x: Any, y: Any, z: Any, t: Any
) -> np.ndarray:
    return evaluate_signed_amplitude_correction(correction, x, y, z, t)[
        "velocity_dt_cartesian_total"
    ]


def _fd4_axis_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
    axis: int,
) -> np.ndarray:
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]

    def shifted(multiplier: float) -> np.ndarray:
        args = [v.copy() for v in coords]
        args[axis] = args[axis] + multiplier * step
        return np.asarray(evaluator(args[0], args[1], args[2], t), dtype=float)

    return (
        shifted(-2.0) - 8.0 * shifted(-1.0)
        + 8.0 * shifted(1.0) - shifted(2.0)
    ) / (12.0 * step)


def _fd4_curl(
    evaluator: VectorEvaluator,
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


def _fd4_divergence(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    d_dx = _fd4_axis_derivative(evaluator, x, y, z, t, step, 0)
    d_dy = _fd4_axis_derivative(evaluator, x, y, z, t, step, 1)
    d_dz = _fd4_axis_derivative(evaluator, x, y, z, t, step, 2)
    return d_dx[..., 0] + d_dy[..., 1] + d_dz[..., 2]


def _fd6_time_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    return (
        -np.asarray(evaluator(x, y, z, t - 3.0 * step), dtype=float)
        + 9.0 * np.asarray(evaluator(x, y, z, t - 2.0 * step), dtype=float)
        - 45.0 * np.asarray(evaluator(x, y, z, t - step), dtype=float)
        + 45.0 * np.asarray(evaluator(x, y, z, t + step), dtype=float)
        - 9.0 * np.asarray(evaluator(x, y, z, t + 2.0 * step), dtype=float)
        + np.asarray(evaluator(x, y, z, t + 3.0 * step), dtype=float)
    ) / (60.0 * step)


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.24, 1.26, 17)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 6
    delta = np.stack(
        (
            0.012 * bump * (1.0 + 0.15 * np.cos(2.0 * math.pi * s)),
            -0.009 * bump * (1.0 - 0.10 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.50,
        producer_kind="analytic-regression-not-agent3-real-candidate",
        provenance=(
            "bounded signed radial regression used only to verify Agent-2 "
            "complete-curl/profile plumbing"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    radii = (0.39, 0.67, 0.97)
    z_values = (-0.61, 0.63)
    times = (0.37, 0.50, 0.63)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    idx = 0
    for time in times:
        for z_value in z_values:
            for radius in radii:
                theta = 0.231 + (idx + 1) * golden
                rows.append(
                    (
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        z_value,
                        time,
                    )
                )
                idx += 1
    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def verification_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()

    potential = lambda xx, yy, zz, tt: correction_vector_potential(
        correction, xx, yy, zz, tt
    )
    velocity = lambda xx, yy, zz, tt: correction_velocity(
        correction, xx, yy, zz, tt
    )
    velocity_dt = lambda xx, yy, zz, tt: correction_velocity_dt(
        correction, xx, yy, zz, tt
    )

    exact = velocity(x, y, z, t)
    exact_dt = velocity_dt(x, y, z, t)
    velocity_rms = _vector_rms(exact)
    velocity_max = float(np.max(np.linalg.norm(exact, axis=-1)))
    dt_rms = _vector_rms(exact_dt)
    dt_max = float(np.max(np.linalg.norm(exact_dt, axis=-1)))

    curl_rows: list[dict[str, float]] = []
    for step in CURL_FD4_STEPS:
        numerical = _fd4_curl(potential, x, y, z, t, step)
        error = numerical - exact
        abs_rms = _vector_rms(error)
        abs_max = float(np.max(np.linalg.norm(error, axis=-1)))
        curl_rows.append(
            {
                "step": step,
                "relative_rms_error": abs_rms / max(velocity_rms, 1.0e-300),
                "relative_max_error": abs_max / max(velocity_max, 1.0e-300),
            }
        )
    curl_refinement = [
        curl_rows[i]["relative_rms_error"]
        / max(curl_rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(curl_rows) - 1)
    ]

    finest = CURL_FD4_STEPS[-1]
    divergence = _fd4_divergence(velocity, x, y, z, t, finest)
    divergence_rms = float(np.sqrt(np.mean(divergence * divergence)))
    divergence_max = float(np.max(np.abs(divergence)))
    divergence_scale = max(velocity_rms / finest, 1.0e-300)

    time_rows: list[dict[str, float]] = []
    for step in TIME_FD6_STEPS:
        numerical = _fd6_time_derivative(velocity, x, y, z, t, step)
        error = numerical - exact_dt
        abs_rms = _vector_rms(error)
        abs_max = float(np.max(np.linalg.norm(error, axis=-1)))
        time_rows.append(
            {
                "step": step,
                "relative_rms_error": abs_rms / max(dt_rms, 1.0e-300),
                "relative_max_error": abs_max / max(dt_max, 1.0e-300),
            }
        )
    time_refinement = [
        time_rows[i]["relative_rms_error"]
        / max(time_rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(time_rows) - 1)
    ]

    exterior = velocity(
        np.asarray((0.0, 0.18, 1.31, 0.75)),
        np.asarray((0.0, 0.0, 0.0, 0.0)),
        np.asarray((0.0, 0.0, 0.0, 2.05)),
        np.asarray((0.5, 0.5, 0.5, 0.5)),
    )
    exterior_max = float(np.max(np.abs(exterior)))

    nodes = np.asarray(correction.radii, dtype=float)
    boundary_values, boundary_dr = _profile_value_derivative(
        correction, np.asarray((nodes[0], nodes[-1]), dtype=float)
    )
    samples = correction.stacked()
    second_boundary = []
    for sign in range(2):
        spline = make_interp_spline(
            nodes,
            samples[:, sign],
            k=5,
            bc_type=(
                [(1, 0.0), (2, 0.0)],
                [(1, 0.0), (2, 0.0)],
            ),
        )
        second_boundary.extend(
            np.asarray(spline.derivative(2)(np.asarray((nodes[0], nodes[-1])))).tolist()
        )
    boundary_c2_max = float(
        max(
            np.max(np.abs(boundary_values)),
            np.max(np.abs(boundary_dr)),
            np.max(np.abs(second_boundary)),
        )
    )

    lift_at_reference, _ = _temporal_lift(
        np.asarray((correction.reference_time,)), correction.reference_time
    )
    reference_lift_error = float(np.max(np.abs(lift_at_reference - 1.0)))

    guards = {
        "nontrivial_velocity": velocity_rms >= 1.0e-8,
        "nontrivial_velocity_dt": dt_rms >= 1.0e-8,
        "curl_finest_relative_rms": curl_rows[-1]["relative_rms_error"] <= 2.0e-5,
        "curl_finest_relative_max": curl_rows[-1]["relative_max_error"] <= 5.0e-5,
        "curl_refinement": min(curl_refinement) >= 8.0,
        "divergence_relative_rms": divergence_rms / divergence_scale <= 2.0e-5,
        "divergence_relative_max": divergence_max / divergence_scale <= 5.0e-5,
        "time_finest_relative_rms": time_rows[-1]["relative_rms_error"] <= 5.0e-8,
        "time_finest_relative_max": time_rows[-1]["relative_max_error"] <= 1.0e-7,
        "time_refinement": min(time_refinement) >= 20.0,
        "support_exterior": exterior_max <= 1.0e-12,
        "compact_c2_boundary": boundary_c2_max <= 1.0e-12,
        "reference_time_lift_identity": reference_lift_error <= 1.0e-14,
    }
    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "sample_count": int(t.size),
        "correction": correction.to_receipt(),
        "curl_fd4_steps": list(CURL_FD4_STEPS),
        "curl_A_vs_complete_curl": curl_rows,
        "curl_rms_refinement_ratios": curl_refinement,
        "velocity_rms": velocity_rms,
        "velocity_sampled_max": velocity_max,
        "divergence_finest_step": finest,
        "divergence_rms": divergence_rms,
        "divergence_max": divergence_max,
        "divergence_relative_rms": divergence_rms / divergence_scale,
        "divergence_relative_max": divergence_max / divergence_scale,
        "time_fd6_steps": list(TIME_FD6_STEPS),
        "time_derivative_comparison": time_rows,
        "time_rms_refinement_ratios": time_refinement,
        "support_exterior_absolute_max": exterior_max,
        "compact_c2_boundary_max_abs": boundary_c2_max,
        "reference_time_lift_error": reference_lift_error,
        "guards": guards,
        "failed_guards": [name for name, passed in guards.items() if not passed],
        "provenance": {
            "source_structure": (
                "signed g=sqrt(epsilon)*a amplitude, complete-curl coefficient product "
                "rule, support-before-curl, conjugate m pair, Q^(-A) scaling"
            ),
            "repository_realization": (
                "public-z pullback, shared-across-beta signed radial profile, C2 quintic "
                "interpolant, frozen autonomous time lift, public gauge potential"
            ),
            "agent3_lane_reimplemented": False,
        },
        "truth_boundary": truth_boundary(),
    }


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(profile_from_delta_a)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "requested_stress",
        "forcing",
        "pressure",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
    }
    return {
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "source_complete_curl_product_rule_reused": True,
        "vector_potential_first_contract": True,
        "signed_radial_delta_a_handoff_executable": True,
        "correction_velocity_and_time_derivative_executable": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_forcing_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "agent3_mean_radial_chain_reimplemented": False,
        "shared_across_beta_profile_is_repository_realization": True,
        "quintic_profile_lift_is_repository_realization": True,
        "autonomous_time_lift_is_repository_realization": True,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "independent_agent4_correction_curl_audit_required": True,
        "real_agent3_delta_a_bound": False,
        "real_full_candidate_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["failed_guards"]:
        raise SystemExit(
            "preregistered signed-amplitude complete-curl guards failed: "
            + ", ".join(receipt["failed_guards"])
        )


if __name__ == "__main__":
    main()
