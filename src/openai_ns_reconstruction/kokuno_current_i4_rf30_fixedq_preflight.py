"""Candidate-bound fixed-Q preflight for the current-I4 RF30 mean-correction seam.

This increment deliberately stops one semantic step before the typed RF30/RF31
backend.  It consumes the exact Agent-3 current-I4 nonlinear backend, which is
itself bound to Agent-2 PR #1080 and its existing PR #960 oscillatory
velocity/Jacobian interface, and makes the corrected-reader fixed-band map
executable:

    R = Q^(-1/2) r,   Z = Q^(-D) z,   T = Q^(-1) (1-t),
    s_Q = q/Q,        X = R^2/(2 s_Q).

The dyadic band is selected canonically from q before any covariance is
observed, with Q <= q < 2Q.  On each physical ring the existing A2 oscillatory
velocity is scaled by Q^A and rotated into cylindrical components.  Its
physical-azimuthal quadratic covariance and the two RF30 derivative observables
needed later,

    d_R <w_r^2>_theta,   d_Z <w_z w_r>_theta,

are computed from the already-frozen A2 Cartesian Jacobian; no curl/Jacobian is
reimplemented here.

Crucially, the corrected Kokuno RF30 bar is normalized Haar average on the
source auxiliary T^2.  The current A2 #1080 field exposes no such auxiliary
T^2 provider.  A physical theta average is not silently relabelled as that Haar
average.  Therefore this module records real candidate observables and the
fixed-Q kinematic map but *does not* authorize an RF30PreMeanState, defect,
correction, or residual claim.  This fail-closed distinction is the purpose of
the increment: it identifies exactly what remains missing before #1126 can
consume a repository-candidate RF30 state.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

import numpy as np

from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_HEAD,
    AGENT2_DIFFERENTIAL_PR,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from .kokuno_oscillatory_quadratic_mean import ANGULAR_ORDERS

TASK = "KOKUNO-A3-CURRENT-I4-RF30-FIXEDQ-PREFLIGHT-126"
SCHEMA = "kokuno-a3-current-i4-rf30-fixedq-preflight-v1"
PARENT_AGENT3_PR = 1150
PARENT_AGENT3_HEAD = "e0b48ccfa964db0caa9d657961d451a21e1d53bf"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FIXED_Q_MAP = (
    "R=Q^(-1/2)r; Z=Q^(-D)z; T=Q^(-1)(1-t); "
    "s_Q=q/Q; X=R^2/(2s_Q); Z=s_Q^D eta; T=s_Q(1-eta^2)"
)
SOURCE_RF30_BAR = "normalized Haar average on source auxiliary T^2"

MAP_RELATIVE_GATE = 5.0e-11
Q_RELATIVE_SPREAD_GATE = 5.0e-13
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4RF30FixedQPreflightError(RuntimeError):
    """Raised when the current-I4 fixed-Q preflight contract is violated."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _finite_scalar(value: Any, label: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise CurrentI4RF30FixedQPreflightError(f"{label} must be finite")
    return out


def _validate_radius(radius: Any) -> np.ndarray:
    out = np.asarray(radius, dtype=float)
    if out.ndim != 1 or out.size < 3:
        raise CurrentI4RF30FixedQPreflightError(
            "radius must be one increasing vector with at least three nodes"
        )
    if not np.all(np.isfinite(out)) or np.any(out <= 0.0):
        raise CurrentI4RF30FixedQPreflightError("radius must be finite and positive")
    if np.any(np.diff(out) <= 0.0):
        raise CurrentI4RF30FixedQPreflightError("radius must be strictly increasing")
    return out


def _relative_difference(new: np.ndarray, old: np.ndarray) -> float:
    n = np.asarray(new, dtype=float)
    o = np.asarray(old, dtype=float)
    denom = float(np.linalg.norm(n.ravel()))
    numer = float(np.linalg.norm((n - o).ravel()))
    if denom == 0.0:
        return 0.0 if numer == 0.0 else math.inf
    return numer / denom


def _canonical_band(q: float) -> tuple[int, float, float]:
    """Choose one deterministic dyadic band before covariance evaluation."""
    if not (math.isfinite(q) and 0.0 < q <= 1.0):
        raise CurrentI4RF30FixedQPreflightError(
            "current fixed-Q preflight requires source-domain 0 < q <= 1"
        )
    ell = int(math.ceil(-math.log2(q)))
    Q = float(math.ldexp(1.0, -ell))
    s_q = q / Q
    tol = 64.0 * np.finfo(float).eps
    if s_q < 1.0 - tol or s_q >= 2.0 + tol:
        raise CurrentI4RF30FixedQPreflightError(
            "canonical dyadic band failed 1 <= q/Q < 2"
        )
    return ell, Q, s_q


@dataclass(frozen=True)
class CurrentI4FixedQObservablePreflight:
    candidate_semantic_sha256: str
    oscillatory_runtime_sha256: str
    leading_semantic_sha256: str
    band_ell: int
    Q: float
    epsilon: float
    A: float
    D: float
    q: float
    s_Q: float
    radius_physical: tuple[float, ...]
    radius_R: tuple[float, ...]
    Z: float
    T: float
    eta: float
    X: tuple[float, ...]
    all_strict_I4: bool
    fixed_q_map_closure_relative_max: float
    q_relative_spread: float
    mean_W_rr_physical_theta: tuple[float, ...]
    mean_W_zr_physical_theta: tuple[float, ...]
    mean_W_thetatheta_physical_theta: tuple[float, ...]
    mean_W_ztheta_physical_theta: tuple[float, ...]
    mean_W_zz_physical_theta: tuple[float, ...]
    dR_mean_W_rr_physical_theta: tuple[float, ...]
    dZ_mean_W_zr_physical_theta: tuple[float, ...]
    angular_orders: tuple[int, ...]
    successive_covariance_relative_differences: tuple[float, ...]
    finest_oscillatory_interior_fraction: float
    finest_covariance_trace_rms: float
    normalized_physical_azimuthal_mean_used: bool
    normalized_source_auxiliary_t2_haar_mean_used: bool
    source_auxiliary_t2_provider_available: bool
    rf30_repository_candidate_state_authorized: bool
    preflight_sha256: str

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "schema": SCHEMA,
            "task": TASK,
            "source_provenance": source_provenance(),
            "truth_boundary": truth_boundary(),
        }


def source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_READER_REPO,
        "commit": SOURCE_READER_HEAD,
        "date": SOURCE_READER_DATE,
        "path": SOURCE_READER_PATH,
        "blob_sha1": SOURCE_READER_BLOB,
        "fixed_q_map": SOURCE_FIXED_Q_MAP,
        "rf30_bar": SOURCE_RF30_BAR,
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _validate_backend(backend: ExactCurrentI4NonlinearBackend) -> Any:
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise CurrentI4RF30FixedQPreflightError("A2 #1080 source blob drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise CurrentI4RF30FixedQPreflightError("A2 #960 differential blob drifted")
    field = backend.composite_field
    leading = getattr(field, "leading_backend", None)
    if leading is None:
        raise CurrentI4RF30FixedQPreflightError("A2 #1080 leading backend is missing")
    for name in ("similarity_coordinates", "velocity"):
        if not callable(getattr(leading, name, None)):
            raise CurrentI4RF30FixedQPreflightError(
                f"A1 #1079 leading backend lacks {name}"
            )
    return leading


def _fixed_q_map(
    backend: ExactCurrentI4NonlinearBackend,
    radius: np.ndarray,
    z: float,
    t: float,
) -> dict[str, object]:
    leading = _validate_backend(backend)
    zeros = np.zeros_like(radius)
    zz = np.full_like(radius, z)
    tt = np.full_like(radius, t)
    coords = leading.similarity_coordinates(radius, zeros, zz, tt)
    q_values = np.asarray(coords["q"], dtype=float)
    eta_values = np.asarray(coords["eta"], dtype=float)
    X_values = np.asarray(coords["X"], dtype=float)
    if q_values.shape != radius.shape or eta_values.shape != radius.shape or X_values.shape != radius.shape:
        raise CurrentI4RF30FixedQPreflightError("A1 similarity coordinates returned unexpected shapes")
    if not all(np.all(np.isfinite(v)) for v in (q_values, eta_values, X_values)):
        raise CurrentI4RF30FixedQPreflightError("A1 similarity coordinates became non-finite")

    q = float(np.mean(q_values))
    q_spread = float(np.max(np.abs(q_values - q), initial=0.0) / max(abs(q), np.finfo(float).tiny))
    if q_spread > Q_RELATIVE_SPREAD_GATE:
        raise CurrentI4RF30FixedQPreflightError(
            "q varied across one fixed-(z,t) radial line"
        )
    ell, Q, s_q = _canonical_band(q)
    A = _finite_scalar(getattr(leading, "A", np.nan), "A")
    D = _finite_scalar(getattr(leading, "D", np.nan), "D")
    h = A - 0.5
    if not (0.0 < h < 0.5 and math.isclose(A + D, 1.0, rel_tol=0.0, abs_tol=5.0e-12)):
        raise CurrentI4RF30FixedQPreflightError("A1 exponents violate A=1/2+h, D=1/2-h")

    R = radius / math.sqrt(Q)
    Z = z / (Q**D)
    T = (1.0 - t) / Q
    X_replay = (R * R) / (2.0 * s_q)
    eta_replay = Z / (s_q**D)
    T_replay = s_q * (1.0 - eta_replay * eta_replay)
    Z_replay = (s_q**D) * eta_replay

    closure = max(
        float(np.max(np.abs(X_replay - X_values), initial=0.0)),
        float(np.max(np.abs(eta_values - eta_replay), initial=0.0)),
        abs(T_replay - T),
        abs(Z_replay - Z),
    )
    scale = max(
        1.0,
        float(np.max(np.abs(X_values), initial=0.0)),
        abs(T),
        abs(Z),
    )
    normalized_closure = closure / scale
    if normalized_closure > MAP_RELATIVE_GATE:
        raise CurrentI4RF30FixedQPreflightError(
            "fixed-Q map disagrees with A1 similarity coordinates"
        )

    x_i4_start = _finite_scalar(getattr(leading, "X_I4_start", np.nan), "X_I4_start")
    x_i4_end = _finite_scalar(getattr(leading, "X_I4_end", np.nan), "X_I4_end")
    if not 0.0 < x_i4_start < x_i4_end:
        raise CurrentI4RF30FixedQPreflightError("A1 I4 interval is malformed")
    all_i4 = bool(np.all((X_values > x_i4_start) & (X_values < x_i4_end)))

    return {
        "leading": leading,
        "ell": ell,
        "Q": Q,
        "s_Q": s_q,
        "A": A,
        "D": D,
        "h": h,
        "R": R,
        "Z": float(Z),
        "T": float(T),
        "eta": float(eta_replay),
        "X": X_values,
        "all_strict_I4": all_i4,
        "fixed_q_map_closure_relative_max": normalized_closure,
        "q_relative_spread": q_spread,
    }


def _covariance_for_order(
    backend: ExactCurrentI4NonlinearBackend,
    radius: np.ndarray,
    z: float,
    t: float,
    Q: float,
    A: float,
    D: float,
    order: int,
) -> tuple[np.ndarray, float]:
    theta = 2.0 * np.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)[None, :]
    s = np.sin(theta)[None, :]
    x = radius[:, None] * c
    y = radius[:, None] * s
    zz = np.full_like(x, z)
    tt = np.full_like(x, t)
    points = np.stack((x, y, zz), axis=-1)
    result = backend.differential_function(points, tt)
    velocity = np.asarray(getattr(result, "velocity"), dtype=float)
    jac = np.asarray(getattr(result, "velocity_jacobian"), dtype=float)
    mask = np.asarray(getattr(result, "interior_mask"), dtype=bool)
    expected_v = (radius.size, order, 3)
    expected_j = (radius.size, order, 3, 3)
    if velocity.shape != expected_v or jac.shape != expected_j or mask.shape != (radius.size, order):
        raise CurrentI4RF30FixedQPreflightError("A2 differential interface returned unexpected shapes")
    if not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(jac)):
        raise CurrentI4RF30FixedQPreflightError("A2 differential interface became non-finite")

    scale = Q**A
    ur = scale * (velocity[..., 0] * c + velocity[..., 1] * s)
    ut = scale * (-velocity[..., 0] * s + velocity[..., 1] * c)
    uz = scale * velocity[..., 2]

    er = np.stack(
        (
            np.broadcast_to(c, (radius.size, order)),
            np.broadcast_to(s, (radius.size, order)),
            np.zeros((radius.size, order), dtype=float),
        ),
        axis=-1,
    )
    radial_directional = np.einsum("...ij,...j->...i", jac, er)
    dR_cart = (Q ** (A + 0.5)) * radial_directional
    dZ_cart = (Q ** (A + D)) * jac[..., :, 2]
    dR_ur = dR_cart[..., 0] * c + dR_cart[..., 1] * s
    dZ_ur = dZ_cart[..., 0] * c + dZ_cart[..., 1] * s
    dZ_uz = dZ_cart[..., 2]

    Wrr = np.mean(ur * ur, axis=1)
    Wzr = np.mean(uz * ur, axis=1)
    Wtt = np.mean(ut * ut, axis=1)
    Wzt = np.mean(uz * ut, axis=1)
    Wzz = np.mean(uz * uz, axis=1)
    dR_Wrr = np.mean(2.0 * ur * dR_ur, axis=1)
    dZ_Wzr = np.mean(dZ_uz * ur + uz * dZ_ur, axis=1)
    values = np.stack((Wrr, Wzr, Wtt, Wzt, Wzz, dR_Wrr, dZ_Wzr), axis=-1)
    if not np.all(np.isfinite(values)):
        raise CurrentI4RF30FixedQPreflightError("physical-angle covariance became non-finite")
    return values, float(np.mean(mask))


def materialize_current_i4_rf30_fixedq_preflight(
    backend: ExactCurrentI4NonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4FixedQObservablePreflight:
    """Materialize real fixed-Q observables without laundering theta mean into RF30 Haar mean."""
    radius_arr = _validate_radius(radius)
    zf = _finite_scalar(z, "z")
    tf = _finite_scalar(t, "t")
    if not (0.0 < tf < 1.0):
        raise CurrentI4RF30FixedQPreflightError("t must satisfy 0 < t < 1")
    fixed = _fixed_q_map(backend, radius_arr, zf, tf)

    samples: list[np.ndarray] = []
    interior_fractions: list[float] = []
    for order in ANGULAR_ORDERS:
        values, fraction = _covariance_for_order(
            backend,
            radius_arr,
            zf,
            tf,
            float(fixed["Q"]),
            float(fixed["A"]),
            float(fixed["D"]),
            int(order),
        )
        samples.append(values)
        interior_fractions.append(fraction)
    convergence = tuple(
        _relative_difference(samples[i], samples[i - 1]) for i in range(1, len(samples))
    )
    finest = samples[-1]
    trace = finest[:, 0] + finest[:, 2] + finest[:, 4]
    trace_rms = float(np.sqrt(np.mean(trace * trace)))

    leading = fixed["leading"]
    leading_semantic = str(getattr(leading, "semantic_sha256", ""))
    if len(leading_semantic) != 64 or any(ch not in "0123456789abcdef" for ch in leading_semantic):
        raise CurrentI4RF30FixedQPreflightError("A1 semantic identity is malformed")

    base_payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_semantic_sha256": backend.composite_semantic_sha256,
        "oscillatory_runtime_sha256": backend.oscillatory_runtime_sha256,
        "leading_semantic_sha256": leading_semantic,
        "band_ell": int(fixed["ell"]),
        "Q": float(fixed["Q"]),
        "A": float(fixed["A"]),
        "D": float(fixed["D"]),
        "q": float(fixed["Q"]) * float(fixed["s_Q"]),
        "s_Q": float(fixed["s_Q"]),
        "radius_physical": radius_arr.tolist(),
        "radius_R": np.asarray(fixed["R"], dtype=float).tolist(),
        "Z": float(fixed["Z"]),
        "T": float(fixed["T"]),
        "eta": float(fixed["eta"]),
        "X": np.asarray(fixed["X"], dtype=float).tolist(),
        "all_strict_I4": bool(fixed["all_strict_I4"]),
        "covariance": finest.tolist(),
        "angular_orders": list(ANGULAR_ORDERS),
        "convergence": list(convergence),
        "source_auxiliary_t2_provider_available": False,
        "rf30_repository_candidate_state_authorized": False,
    }
    preflight_sha = _canonical_sha(base_payload)
    return CurrentI4FixedQObservablePreflight(
        candidate_semantic_sha256=backend.composite_semantic_sha256,
        oscillatory_runtime_sha256=backend.oscillatory_runtime_sha256,
        leading_semantic_sha256=leading_semantic,
        band_ell=int(fixed["ell"]),
        Q=float(fixed["Q"]),
        epsilon=float(fixed["Q"]) ** float(fixed["h"]),
        A=float(fixed["A"]),
        D=float(fixed["D"]),
        q=float(fixed["Q"]) * float(fixed["s_Q"]),
        s_Q=float(fixed["s_Q"]),
        radius_physical=tuple(float(v) for v in radius_arr),
        radius_R=tuple(float(v) for v in np.asarray(fixed["R"], dtype=float)),
        Z=float(fixed["Z"]),
        T=float(fixed["T"]),
        eta=float(fixed["eta"]),
        X=tuple(float(v) for v in np.asarray(fixed["X"], dtype=float)),
        all_strict_I4=bool(fixed["all_strict_I4"]),
        fixed_q_map_closure_relative_max=float(fixed["fixed_q_map_closure_relative_max"]),
        q_relative_spread=float(fixed["q_relative_spread"]),
        mean_W_rr_physical_theta=tuple(float(v) for v in finest[:, 0]),
        mean_W_zr_physical_theta=tuple(float(v) for v in finest[:, 1]),
        mean_W_thetatheta_physical_theta=tuple(float(v) for v in finest[:, 2]),
        mean_W_ztheta_physical_theta=tuple(float(v) for v in finest[:, 3]),
        mean_W_zz_physical_theta=tuple(float(v) for v in finest[:, 4]),
        dR_mean_W_rr_physical_theta=tuple(float(v) for v in finest[:, 5]),
        dZ_mean_W_zr_physical_theta=tuple(float(v) for v in finest[:, 6]),
        angular_orders=tuple(int(v) for v in ANGULAR_ORDERS),
        successive_covariance_relative_differences=convergence,
        finest_oscillatory_interior_fraction=float(interior_fractions[-1]),
        finest_covariance_trace_rms=trace_rms,
        normalized_physical_azimuthal_mean_used=True,
        normalized_source_auxiliary_t2_haar_mean_used=False,
        source_auxiliary_t2_provider_available=False,
        rf30_repository_candidate_state_authorized=False,
        preflight_sha256=preflight_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf30_fixedq_preflight)
    forbidden = {
        "residual", "defect", "forcing", "pressure", "target", "gain", "damping",
        "threshold", "viscosity", "nu", "correction", "haar", "surrogate",
        "amplitude", "phase", "spatial_step", "derivative_step", "angular_order",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "source_reader_head": SOURCE_READER_HEAD,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(signature.parameters)),
        "fixed_q_coordinate_map_materialized": True,
        "canonical_band_selected_before_covariance_observation": True,
        "actual_a2_oscillatory_velocity_consumed": True,
        "actual_a2_oscillatory_jacobian_consumed": True,
        "agent2_curl_or_jacobian_reimplemented": False,
        "physical_azimuthal_covariance_materialized": True,
        "physical_azimuthal_covariance_derivatives_materialized": True,
        "source_auxiliary_t2_provider_available": False,
        "normalized_source_auxiliary_t2_haar_mean_used": False,
        "physical_theta_mean_promoted_to_source_auxiliary_haar": False,
        "rf30_repository_candidate_state_authorized": False,
        "rf30_defect_materialized": False,
        "rf31_five_row_system_materialized_from_this_preflight": False,
        "correction_velocity_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proof_claimed": False,
        "pde_validated": False,
    }
