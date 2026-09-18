"""Execute Kokuno's real homogeneous primary pulse on a supplied source path.

The corrected 2026-09-09 reader separates the zero-data forced inverse from
its *real homogeneous primary pulse*.  For the latter, at harmonic m=1 and
f=0, the source uses the moving transverse plane

    U = [e_r - s_a K_a, N_a],
    M = [[1, 1], [c0*sqrt(1+s(v)^2), -c0*sqrt(1+s(v)^2)]],
    B = U M,

and coordinate data ``z_+(0)=P(0)>0, z_-(0)=0``.  The physical transverse
amplitude is ``t=Bz``.  Its envelope is

    P(v)=exp(int_{L_s/2}^v a_net(|s(w)|) dw),
    s(v)=sigma*(u_*/2 + u_* v/L_s),
    a_net(y)=lambda0/sqrt(1+y^2)
             -lambda0*(1+y^2)/(1+u_*^2)^(3/2).

This module evaluates those identities and reuses the repository's guarded
RK4 projected-pulse stepper with f=0 and the source nonzero entrance datum.
The envelope integral is evaluated by an exact elementary antiderivative;
RK4/path interpolation and end-step projection remain autonomous numerical
approximations and are reported as such.

The caller must still supply the actual source path ``n_Phi, n_Phi', K`` and
the pulse constants.  No positive-order background, discrete source label,
auxiliary-torus evaluation, support partition, D_r/D_z curl sensitivity, or
direct ``velocity(x,y,z,t)`` is invented here.  This is therefore not a
paper-exact or OpenAI-identified field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_pulse_inverse import KokunoProjectedPulseInverseContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-homogeneous-primary-pulse-v1"

_SOURCE_FORMULAS = {
    "pulse_coordinate": "s(v)=sigma*(u_*/2+u_*v/L_s)",
    "net_growth": (
        "a_net(y)=lambda0/sqrt(1+y^2)-"
        "lambda0*(1+y^2)/(1+u_*^2)^(3/2)"
    ),
    "envelope": "P(v)=exp(int_(L_s/2)^v a_net(|s(w)|)dw)",
    "moving_plane": (
        "K_a=n_tan/|n_tan|; N_a=((K_a)_z,-(K_a)_theta); "
        "s_a=n_r/|n_tan|; U=[e_r-s_a*K_a,N_a]"
    ),
    "eigenbasis": (
        "M=[[1,1],[c0*sqrt(1+s^2),-c0*sqrt(1+s^2)]]; B=U*M"
    ),
    "entrance_data": "m=1,f=0: z_+(0)=P(0), z_-(0)=0; t(0)=B(0)z(0)",
    "real_fourier_pair": "for real t*cos(k*Phi), t_+=t_-=t/2",
}

_TRUTH_BOUNDARY = {
    "source_real_homogeneous_primary_pulse_equations_executable": True,
    "source_envelope_exact_antiderivative_executable": True,
    "source_moving_plane_entrance_datum_executable": True,
    "real_cosine_m_plus_one_coefficient_executable": True,
    "source_complete_curl_C_plus_executable": True,
    "source_actual_discrete_label_selected": False,
    "source_actual_background_path_instantiated": False,
    "positive_order_background_corrections_included": False,
    "source_actual_auxiliary_torus_evaluation_instantiated": False,
    "source_actual_support_partition_instantiated": False,
    "source_homogeneous_D_r_C_plus_instantiated": False,
    "source_homogeneous_D_z_C_plus_instantiated": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_real(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


@dataclass(frozen=True)
class KokunoSourceHomogeneousPrimaryPulse:
    """Source homogeneous m=1 pulse on caller-supplied path geometry."""

    epsilon: float = 0.25
    u_star: float = 2.0
    lambda0: float = 1.0
    c0: float = -0.5
    L_s: float = 1.0
    sigma: int = 1
    min_tangential_norm: float = 1.0e-10
    reality_atol: float = 2.0e-11

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        u_star = float(self.u_star)
        lambda0 = float(self.lambda0)
        c0 = float(self.c0)
        L_s = float(self.L_s)
        min_tangential_norm = float(self.min_tangential_norm)
        reality_atol = float(self.reality_atol)
        sigma = int(self.sigma)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        if not np.isfinite(u_star) or not (0.05 <= u_star <= 20.0):
            raise ValueError("u_star must lie in [0.05,20]")
        if not np.isfinite(lambda0) or not (1.0e-8 <= lambda0 <= 100.0):
            raise ValueError("lambda0 must lie in [1e-8,100]")
        if not np.isfinite(c0) or not (-20.0 <= c0 <= -1.0e-8):
            raise ValueError("c0 must be strictly negative and bounded")
        if not np.isfinite(L_s) or not (1.0e-4 <= L_s <= 100.0):
            raise ValueError("L_s must lie in [1e-4,100]")
        if sigma not in (-1, 1) or sigma != self.sigma:
            raise ValueError("sigma must be exactly +1 or -1")
        if not np.isfinite(min_tangential_norm) or not (
            0.0 < min_tangential_norm <= 1.0e-4
        ):
            raise ValueError("min_tangential_norm must lie in (0,1e-4]")
        if not np.isfinite(reality_atol) or not (0.0 < reality_atol <= 1.0e-7):
            raise ValueError("reality_atol must lie in (0,1e-7]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "u_star", u_star)
        object.__setattr__(self, "lambda0", lambda0)
        object.__setattr__(self, "c0", c0)
        object.__setattr__(self, "L_s", L_s)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "min_tangential_norm", min_tangential_norm)
        object.__setattr__(self, "reality_atol", reality_atol)

    @property
    def pulse_solver(self) -> KokunoProjectedPulseInverseContract:
        return KokunoProjectedPulseInverseContract(epsilon=self.epsilon, m=1)

    def s(self, pulse_v: Any) -> np.ndarray:
        v = _finite_real(pulse_v, "pulse_v")
        return self.sigma * self.u_star * (0.5 + v / self.L_s)

    def a_net(self, pulse_v: Any) -> np.ndarray:
        y = np.abs(self.s(pulse_v))
        return self.lambda0 / np.sqrt(1.0 + y * y) - self.lambda0 * (
            1.0 + y * y
        ) / (1.0 + self.u_star**2) ** 1.5

    def log_envelope(self, pulse_v: Any) -> np.ndarray:
        """Return exact ``log P(v)`` from the elementary antiderivative."""
        v = _finite_real(pulse_v, "pulse_v")
        y = np.abs(self.s(v))
        u = self.u_star
        denom = (1.0 + u * u) ** 1.5

        def primitive(value: np.ndarray | float) -> np.ndarray:
            value = np.asarray(value, dtype=float)
            return np.arcsinh(value) - (value + value**3 / 3.0) / denom

        return (self.lambda0 * self.L_s / u) * (primitive(y) - primitive(u))

    def envelope(self, pulse_v: Any) -> np.ndarray:
        value = np.exp(self.log_envelope(pulse_v))
        if not np.all(np.isfinite(value)) or np.any(value <= 0.0):
            raise RuntimeError("source homogeneous envelope became nonfinite/nonpositive")
        return value

    def _validate_path_geometry(
        self, pulse_v: Any, n_phi: Any, n_phi_prime: Any, K: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        v = _finite_real(pulse_v, "pulse_v")
        if v.ndim != 1 or v.size < 3:
            raise ValueError("pulse_v must be a one-dimensional path with at least three nodes")
        if np.any(np.diff(v) <= 0.0):
            raise ValueError("pulse_v nodes must be strictly increasing")
        tol = 5.0e-13 * max(1.0, self.L_s)
        if abs(float(v[0])) > tol or abs(float(v[-1] - self.L_s)) > tol:
            raise ValueError("homogeneous source path must run from v=0 to v=L_s")
        n = _finite_real(n_phi, "n_phi")
        n_prime = _finite_real(n_phi_prime, "n_phi_prime")
        K = _finite_real(K, "K")
        if n.shape != (v.size, 3) or n_prime.shape != (v.size, 3):
            raise ValueError("n_phi and n_phi_prime must have shape (len(pulse_v),3)")
        if K.shape != (v.size, 3, 3):
            raise ValueError("K must have shape (len(pulse_v),3,3)")
        tangential_norm = np.linalg.norm(n[:, 1:], axis=-1)
        if np.any(tangential_norm <= self.min_tangential_norm):
            raise ValueError("n_Phi tangential norm is below the moving-plane guard")
        return v, n, n_prime, K

    def moving_plane(self, pulse_v: Any, n_phi: Any) -> dict[str, np.ndarray]:
        v = _finite_real(pulse_v, "pulse_v")
        n = _finite_real(n_phi, "n_phi")
        if n.shape != v.shape + (3,):
            raise ValueError("n_phi must have shape pulse_v.shape+(3,)")
        n_tan = n[..., 1:]
        tan_norm = np.linalg.norm(n_tan, axis=-1)
        if np.any(tan_norm <= self.min_tangential_norm):
            raise ValueError("n_Phi tangential norm is below the moving-plane guard")
        K_a = n_tan / tan_norm[..., None]
        N_a = np.stack((K_a[..., 1], -K_a[..., 0]), axis=-1)
        s_a = n[..., 0] / tan_norm
        u1 = np.concatenate(
            (np.ones(v.shape + (1,)), -s_a[..., None] * K_a), axis=-1
        )
        u2 = np.concatenate((np.zeros(v.shape + (1,)), N_a), axis=-1)
        s_ref = self.s(v)
        gamma = self.c0 * np.sqrt(1.0 + s_ref * s_ref)
        B_plus = u1 + gamma[..., None] * u2
        B_minus = u1 - gamma[..., None] * u2
        B = np.stack((B_plus, B_minus), axis=-1)
        transverse_defect = np.max(np.abs(np.einsum("...i,...ij->...j", n, B)))
        scale = max(1.0, float(np.max(np.abs(B))))
        if transverse_defect > 2.0e-11 * scale:
            raise RuntimeError("moving-plane basis failed n_Phi-perpendicularity")
        return {
            "K_a": K_a,
            "N_a": N_a,
            "s_a": s_a,
            "s": s_ref,
            "gamma": gamma,
            "B": B,
            "transverse_defect_abs": np.asarray(transverse_defect),
        }

    def plane_coordinates(
        self, pulse_v: Any, n_phi: Any, t_m: Any
    ) -> dict[str, np.ndarray]:
        v = _finite_real(pulse_v, "pulse_v")
        t = np.asarray(t_m, dtype=np.complex128)
        if t.shape != v.shape + (3,) or not np.all(np.isfinite(t.real)) or not np.all(
            np.isfinite(t.imag)
        ):
            raise ValueError("t_m must be a finite complex path with trailing dimension 3")
        plane = self.moving_plane(v, n_phi)
        w0 = t[..., 0]
        w1 = np.sum(plane["N_a"] * t[..., 1:], axis=-1)
        gamma = plane["gamma"]
        z_plus = 0.5 * (w0 + w1 / gamma)
        z_minus = 0.5 * (w0 - w1 / gamma)
        z = np.stack((z_plus, z_minus), axis=-1)
        return {
            **plane,
            "z": z,
            "x": z_plus + z_minus,
            "y": gamma * (z_plus - z_minus),
        }

    def solve_path(
        self, pulse_v: Any, n_phi: Any, n_phi_prime: Any, K: Any
    ) -> dict[str, Any]:
        """Solve the source real homogeneous m=1 pulse on supplied geometry."""
        v, n, n_prime, K = self._validate_path_geometry(pulse_v, n_phi, n_phi_prime, K)
        plane = self.moving_plane(v, n)
        P = self.envelope(v)
        z0 = np.array([P[0], 0.0], dtype=float)
        t0 = plane["B"][0] @ z0
        if abs(float(np.dot(n[0], t0))) > 2.0e-11 * max(1.0, float(np.linalg.norm(t0))):
            raise RuntimeError("source homogeneous entrance datum is not transverse")

        forcing = np.zeros((v.size, 3), dtype=np.complex128)
        solved = self.pulse_solver.solve_path(v, n, n_prime, K, forcing, t0=t0)
        t_complex = np.asarray(solved["t_m"], dtype=np.complex128)
        imag_scale = max(1.0, float(np.max(np.abs(t_complex.real))))
        if float(np.max(np.abs(t_complex.imag))) > self.reality_atol * imag_scale:
            raise RuntimeError("real homogeneous pulse acquired excessive imaginary leakage")
        t_real = t_complex.real
        coords = self.plane_coordinates(v, n, t_real)
        initial_coordinate_error = float(np.max(np.abs(coords["z"][0].real - z0)))
        if initial_coordinate_error > 5.0e-11 * max(1.0, float(P[0])):
            raise RuntimeError("moving-plane entrance coordinates do not recover source z(0)")

        t_plus = 0.5 * t_real.astype(np.complex128)
        C_plus = self.pulse_solver.complete_curl_coefficients(n, t_plus)
        return {
            **solved,
            "t_real_cosine": t_real,
            "t_plus": t_plus,
            "C_plus": C_plus,
            "envelope_P": P,
            "log_envelope_P": self.log_envelope(v),
            "s": plane["s"],
            "moving_plane_B": plane["B"],
            "z_coordinates": coords["z"],
            "x_coordinate": coords["x"],
            "y_coordinate": coords["y"],
            "source_initial_z": z0,
            "source_initial_t": t0,
            "initial_coordinate_error_abs": initial_coordinate_error,
            "forcing_origin": "source homogeneous primary pulse: f=0",
            "fourier_convention": "real t*cos(k*Phi) -> t_plus=t/2 at m=+1",
        }

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "real homogeneous primary pulse and moving-plane entrance data",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "u_star": self.u_star,
                "lambda0": self.lambda0,
                "c0": self.c0,
                "L_s": self.L_s,
                "sigma": self.sigma,
                "min_tangential_norm": self.min_tangential_norm,
                "reality_atol": self.reality_atol,
                "origin": (
                    "pulse constants are caller-declared bounded replay values unless supplied by the "
                    "source label constructor; defaults are not claimed recovered"
                ),
            },
            "numerical_method": {
                "envelope": "exact elementary antiderivative of the source a_net formula",
                "pulse_stepper": "existing RK4 with linear midpoint coefficient interpolation",
                "transverse_stabilization": "existing guarded end-step orthogonal projection",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceHomogeneousPrimaryPulse":
        if not isinstance(payload, dict):
            raise ValueError("homogeneous pulse payload must be an object")
        required = {"schema", "source", "parameters", "numerical_method", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("homogeneous pulse payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported homogeneous pulse schema")
        params = payload["parameters"]
        keys = {
            "epsilon", "u_star", "lambda0", "c0", "L_s", "sigma",
            "min_tangential_norm", "reality_atol", "origin"
        }
        if not isinstance(params, dict) or set(params) != keys:
            raise ValueError("homogeneous pulse parameters changed")
        obj = cls(
            epsilon=params["epsilon"],
            u_star=params["u_star"],
            lambda0=params["lambda0"],
            c0=params["c0"],
            L_s=params["L_s"],
            sigma=params["sigma"],
            min_tangential_norm=params["min_tangential_norm"],
            reality_atol=params["reality_atol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"homogeneous pulse {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("homogeneous pulse sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceHomogeneousPrimaryPulse":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
