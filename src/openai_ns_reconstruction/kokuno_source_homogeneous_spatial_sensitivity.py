"""Spatial sensitivity of Kokuno's real homogeneous primary pulse.

The corrected 2026-09-09 reader gives the real homogeneous primary pulse and,
on a band rectangle, the exact identities ``D_r v=D_z v=0`` for its pulse
coordinate ``v``.  This module keeps those source statements separate from the
repository-derived differentiation used here.

For one caller-selected source-normalized spatial derivative ``xi`` in
``{D_r,D_z}``, the caller supplies ``n_Phi, n_Phi', K`` and their xi
derivatives along the same pulse path.  Because ``xi v=0`` and the pulse
constants are held fixed, the envelope and eigenvalue factor have zero xi
derivative.  Differentiating the moving transverse plane gives the nonzero
entrance sensitivity

    t_xi(0) = B_xi(0) [P(0), 0]^T.

The coupled pulse/sensitivity ODE then uses the already implemented analytic
projected-pulse sensitivity equation.  This closes the nonzero-entrance seam
needed to obtain ``partial_xi t_+`` and ``partial_xi C_+`` for the homogeneous
source pulse without fitting velocity components or inventing forcing.

Actual positive-order background data, auxiliary-torus derivatives, source
labels/support partitions, and direct ``velocity(x,y,z,t)`` remain external.
The result is therefore not paper-exact and is not an identified OpenAI field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_homogeneous_pulse import KokunoSourceHomogeneousPrimaryPulse
from .kokuno_source_pulse_sensitivity import KokunoProjectedPulseSensitivityContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-homogeneous-spatial-sensitivity-v1"

_SOURCE_AND_DERIVED_FORMULAS = {
    "source_spatial_pulse_coordinate": "on a band rectangle: D_r v=D_z v=0, t_* v=1",
    "source_homogeneous_entrance": "z_+(0)=P(0), z_-(0)=0; t(0)=B(0)z(0)",
    "source_real_fourier_pair": "real t*cos(k*Phi) -> t_+=t/2 at m=+1",
    "derived_moving_plane_sensitivity": (
        "B_xi is obtained by differentiating K_a=n_tan/|n_tan|, "
        "N_a=(K_a,z,-K_a,theta), s_a=n_r/|n_tan| with xi(v)=0"
    ),
    "derived_entrance_sensitivity": "t_xi(0)=B_xi(0)*[P(0),0]^T",
    "derived_positive_mode_sensitivity": (
        "t_plus_xi=t_xi/2; C_plus_xi=partial_xi[i*(n_Phi cross t_plus)/(k*|n_Phi|^2)]"
    ),
}

_TRUTH_BOUNDARY = {
    "source_homogeneous_primary_pulse_reused": True,
    "source_spatial_pulse_coordinate_invariance_used": True,
    "derived_nonzero_entrance_sensitivity_executable": True,
    "source_homogeneous_t_plus_sensitivity_executable": True,
    "source_homogeneous_C_plus_sensitivity_executable": True,
    "source_actual_background_path_instantiated": False,
    "positive_order_background_corrections_included": False,
    "source_actual_auxiliary_torus_evaluation_instantiated": False,
    "source_actual_support_partition_instantiated": False,
    "source_actual_D_r_path_instantiated": False,
    "source_actual_D_z_path_instantiated": False,
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
class KokunoSourceHomogeneousSpatialSensitivity:
    """Differentiate the real homogeneous primary pulse in ``D_r`` or ``D_z``."""

    derivative_label: str = "D_z"
    epsilon: float = 0.25
    u_star: float = 2.0
    lambda0: float = 1.0
    c0: float = -0.5
    L_s: float = 1.0
    sigma: int = 1
    min_tangential_norm: float = 1.0e-10
    reality_atol: float = 2.0e-11

    def __post_init__(self) -> None:
        if self.derivative_label not in ("D_r", "D_z"):
            raise ValueError("derivative_label must be exactly 'D_r' or 'D_z'")
        pulse = self.pulse
        object.__setattr__(self, "epsilon", pulse.epsilon)
        object.__setattr__(self, "u_star", pulse.u_star)
        object.__setattr__(self, "lambda0", pulse.lambda0)
        object.__setattr__(self, "c0", pulse.c0)
        object.__setattr__(self, "L_s", pulse.L_s)
        object.__setattr__(self, "sigma", pulse.sigma)
        object.__setattr__(self, "min_tangential_norm", pulse.min_tangential_norm)
        object.__setattr__(self, "reality_atol", pulse.reality_atol)

    @property
    def pulse(self) -> KokunoSourceHomogeneousPrimaryPulse:
        return KokunoSourceHomogeneousPrimaryPulse(
            epsilon=self.epsilon,
            u_star=self.u_star,
            lambda0=self.lambda0,
            c0=self.c0,
            L_s=self.L_s,
            sigma=self.sigma,
            min_tangential_norm=self.min_tangential_norm,
            reality_atol=self.reality_atol,
        )

    @property
    def sensitivity(self) -> KokunoProjectedPulseSensitivityContract:
        return KokunoProjectedPulseSensitivityContract(epsilon=self.epsilon, m=1)

    def moving_plane_sensitivity(
        self, pulse_v: Any, n_phi: Any, n_phi_xi: Any
    ) -> dict[str, np.ndarray]:
        """Differentiate the source moving plane at fixed pulse coordinate ``v``."""
        v = _finite_real(pulse_v, "pulse_v")
        n = _finite_real(n_phi, "n_phi")
        n_xi = _finite_real(n_phi_xi, "n_phi_xi")
        if n.shape != v.shape + (3,) or n_xi.shape != n.shape:
            raise ValueError("n_phi and n_phi_xi must have shape pulse_v.shape+(3,)")

        plane = self.pulse.moving_plane(v, n)
        n_tan = n[..., 1:]
        n_tan_xi = n_xi[..., 1:]
        rho = np.linalg.norm(n_tan, axis=-1)
        if np.any(rho <= self.min_tangential_norm):
            raise ValueError("n_Phi tangential norm is below the moving-plane guard")
        rho_xi = np.sum(n_tan * n_tan_xi, axis=-1) / rho

        K_a = plane["K_a"]
        K_a_xi = n_tan_xi / rho[..., None] - (
            n_tan * rho_xi[..., None] / (rho * rho)[..., None]
        )
        N_a_xi = np.stack((K_a_xi[..., 1], -K_a_xi[..., 0]), axis=-1)
        s_a = plane["s_a"]
        s_a_xi = n_xi[..., 0] / rho - n[..., 0] * rho_xi / (rho * rho)

        u1_tan_xi = -(
            s_a_xi[..., None] * K_a + s_a[..., None] * K_a_xi
        )
        u1_xi = np.concatenate((np.zeros(v.shape + (1,)), u1_tan_xi), axis=-1)
        u2_xi = np.concatenate((np.zeros(v.shape + (1,)), N_a_xi), axis=-1)

        # Source D_r v=D_z v=0 and fixed pulse constants imply gamma_xi=0.
        gamma = plane["gamma"]
        B_plus_xi = u1_xi + gamma[..., None] * u2_xi
        B_minus_xi = u1_xi - gamma[..., None] * u2_xi
        B_xi = np.stack((B_plus_xi, B_minus_xi), axis=-1)

        differentiated_defect = (
            np.einsum("...i,...ij->...j", n_xi, plane["B"])
            + np.einsum("...i,...ij->...j", n, B_xi)
        )
        defect_abs = float(np.max(np.abs(differentiated_defect)))
        scale = max(
            1.0,
            float(np.max(np.abs(n_xi))) * float(np.max(np.abs(plane["B"]))),
            float(np.max(np.abs(n))) * float(np.max(np.abs(B_xi))),
        )
        if defect_abs > 5.0e-11 * scale:
            raise RuntimeError("differentiated moving plane failed transverse identity")

        return {
            **plane,
            "rho_tan": rho,
            "rho_tan_xi": rho_xi,
            "K_a_xi": K_a_xi,
            "N_a_xi": N_a_xi,
            "s_a_xi": s_a_xi,
            "B_xi": B_xi,
            "differentiated_transverse_defect_abs": np.asarray(defect_abs),
        }

    def _validate_derivative_geometry(
        self,
        v: np.ndarray,
        n: np.ndarray,
        n_prime: np.ndarray,
        K: np.ndarray,
        n_xi: Any,
        n_prime_xi: Any,
        K_xi: Any,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_xi = _finite_real(n_xi, "n_phi_xi")
        n_prime_xi = _finite_real(n_prime_xi, "n_phi_prime_xi")
        K_xi = _finite_real(K_xi, "K_xi")
        if n_xi.shape != n.shape:
            raise ValueError("n_phi_xi must have shape (len(pulse_v),3)")
        if n_prime_xi.shape != n_prime.shape:
            raise ValueError("n_phi_prime_xi must have shape (len(pulse_v),3)")
        if K_xi.shape != K.shape:
            raise ValueError("K_xi must have shape (len(pulse_v),3,3)")
        return n_xi, n_prime_xi, K_xi

    def solve_path(
        self,
        pulse_v: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        n_phi_xi: Any,
        n_phi_prime_xi: Any,
        K_xi: Any,
    ) -> dict[str, Any]:
        """Solve the homogeneous pulse and one source-normalized spatial sensitivity."""
        pulse = self.pulse
        v, n, n_prime, K = pulse._validate_path_geometry(pulse_v, n_phi, n_phi_prime, K)
        n_xi, n_prime_xi, K_xi = self._validate_derivative_geometry(
            v, n, n_prime, K, n_phi_xi, n_phi_prime_xi, K_xi
        )
        plane = self.moving_plane_sensitivity(v, n, n_xi)
        P_env = pulse.envelope(v)
        z0 = np.array([P_env[0], 0.0], dtype=float)
        t0 = plane["B"][0] @ z0
        s0 = plane["B_xi"][0] @ z0

        entrance_constraint = float(abs(np.dot(n[0], t0)))
        entrance_sensitivity_constraint = float(
            abs(np.dot(n[0], s0) + np.dot(n_xi[0], t0))
        )
        entrance_scale = max(1.0, float(np.linalg.norm(t0)), float(np.linalg.norm(s0)))
        if entrance_constraint > 3.0e-11 * entrance_scale:
            raise RuntimeError("homogeneous entrance datum is not transverse")
        if entrance_sensitivity_constraint > 5.0e-11 * entrance_scale:
            raise RuntimeError("homogeneous entrance sensitivity violates differentiated constraint")

        forcing = np.zeros((v.size, 3), dtype=np.complex128)
        forcing_xi = np.zeros_like(forcing)
        contract = self.sensitivity
        t = np.zeros((v.size, 3), dtype=np.complex128)
        s = np.zeros_like(t)
        t[0] = t0
        s[0] = s0
        t_correction = np.zeros(v.size, dtype=float)
        s_correction = np.zeros(v.size, dtype=float)

        def pair_rhs(coeffs, y, z):
            nn, np_, KK, ff, nn_x, np_x, KK_x, ff_x = coeffs
            return (
                contract.pulse.projected_rhs(nn, np_, KK, ff, y),
                contract.sensitivity_rhs(
                    nn, np_, KK, ff, nn_x, np_x, KK_x, ff_x, y, z
                ),
            )

        for j in range(v.size - 1):
            step = float(v[j + 1] - v[j])
            left = (
                n[j], n_prime[j], K[j], forcing[j],
                n_xi[j], n_prime_xi[j], K_xi[j], forcing_xi[j],
            )
            right = (
                n[j + 1], n_prime[j + 1], K[j + 1], forcing[j + 1],
                n_xi[j + 1], n_prime_xi[j + 1], K_xi[j + 1], forcing_xi[j + 1],
            )
            mid = tuple(0.5 * (a + b) for a, b in zip(left, right))
            y = t[j]
            z = s[j]
            k1y, k1z = pair_rhs(left, y, z)
            k2y, k2z = pair_rhs(mid, y + 0.5 * step * k1y, z + 0.5 * step * k1z)
            k3y, k3z = pair_rhs(mid, y + 0.5 * step * k2y, z + 0.5 * step * k2z)
            k4y, k4z = pair_rhs(right, y + step * k3y, z + step * k3z)
            raw_t = y + (step / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
            raw_s = z + (step / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z)

            _, _, _, _, projector, projector_xi = contract._projector_data(
                n[j + 1], n_xi[j + 1]
            )
            stable_t = projector @ raw_t
            stable_s = projector @ raw_s + projector_xi @ raw_t
            tc = float(np.linalg.norm(stable_t - raw_t))
            sc = float(np.linalg.norm(stable_s - raw_s))
            if tc > contract.max_projection_correction * max(1.0, float(np.linalg.norm(raw_t))):
                raise ValueError("pulse projection correction exceeded guard; refine path or fix inputs")
            if sc > contract.max_sensitivity_projection_correction * max(
                1.0, float(np.linalg.norm(raw_s))
            ):
                raise ValueError(
                    "sensitivity projection correction exceeded guard; refine path or fix derivatives"
                )
            t[j + 1] = stable_t
            s[j + 1] = stable_s
            t_correction[j + 1] = tc
            s_correction[j + 1] = sc

        pulse_constraint = np.abs(np.sum(n * t, axis=-1))
        sensitivity_constraint = np.abs(
            np.sum(n * s, axis=-1) + np.sum(n_xi * t, axis=-1)
        )
        real_scale = max(1.0, float(np.max(np.abs(t.real))), float(np.max(np.abs(s.real))))
        imag_leak = max(float(np.max(np.abs(t.imag))), float(np.max(np.abs(s.imag))))
        if imag_leak > self.reality_atol * real_scale:
            raise RuntimeError("real homogeneous pulse sensitivity acquired excessive imaginary leakage")

        t_real = t.real
        s_real = s.real
        t_plus = 0.5 * t_real.astype(np.complex128)
        t_plus_xi = 0.5 * s_real.astype(np.complex128)
        C_plus = contract.pulse.complete_curl_coefficients(n, t_plus)
        C_plus_xi = contract.coefficient_sensitivity(n, n_xi, t_plus, t_plus_xi)
        return {
            "path": v.copy(),
            "t_real_cosine": t_real,
            "t_real_cosine_xi": s_real,
            "t_plus": t_plus,
            "t_plus_xi": t_plus_xi,
            "C_plus": C_plus,
            "C_plus_xi": C_plus_xi,
            "moving_plane_B": plane["B"],
            "moving_plane_B_xi": plane["B_xi"],
            "source_initial_z": z0,
            "source_initial_t": t0,
            "source_initial_t_xi": s0,
            "pulse_constraint_defect_abs": pulse_constraint,
            "sensitivity_constraint_defect_abs": sensitivity_constraint,
            "pulse_projection_correction_norm": t_correction,
            "sensitivity_projection_correction_norm": s_correction,
            "max_pulse_constraint_defect_abs": float(np.max(pulse_constraint)),
            "max_sensitivity_constraint_defect_abs": float(np.max(sensitivity_constraint)),
            "max_pulse_projection_correction_norm": float(np.max(t_correction)),
            "max_sensitivity_projection_correction_norm": float(np.max(s_correction)),
            "entrance_constraint_defect_abs": entrance_constraint,
            "entrance_sensitivity_constraint_defect_abs": entrance_sensitivity_constraint,
            "derivative_label": self.derivative_label,
            "forcing_origin": "source homogeneous primary pulse: f=0 and spatial derivative f_xi=0",
            "spatial_coordinate_contract": "source band rectangle has D_r v=D_z v=0",
            "integration_method": (
                "coupled RK4 with linear midpoint coefficient interpolation; differentiated "
                "end-step transverse projection is autonomous numerical stabilization"
            ),
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
                "formula_scope": (
                    "real homogeneous primary pulse plus source band-rectangle D_r v=D_z v=0 identity"
                ),
                "formulas": dict(_SOURCE_AND_DERIVED_FORMULAS),
                "derivation_boundary": (
                    "pulse, entrance data and D_r/D_z pulse-coordinate invariance are sourced; "
                    "moving-plane/nonzero-entrance sensitivities are analytic repository derivations"
                ),
            },
            "parameters": {
                "derivative_label": self.derivative_label,
                "epsilon": self.epsilon,
                "u_star": self.u_star,
                "lambda0": self.lambda0,
                "c0": self.c0,
                "L_s": self.L_s,
                "sigma": self.sigma,
                "min_tangential_norm": self.min_tangential_norm,
                "reality_atol": self.reality_atol,
            },
            "numerical_method": {
                "pulse_stepper": "coupled RK4 with linear midpoint coefficient interpolation",
                "transverse_stabilization": "guarded projection and exact derivative of that numerical projection",
                "origin": "autonomous numerical approximation; not a source hidden discretization",
            },
            "routing": {
                "caller_must_supply": (
                    "n_Phi,n_Phi_prime,K and their selected source-normalized D_r or D_z derivatives"
                ),
                "output": "t_plus, selected spatial derivative of t_plus, C_plus, selected spatial derivative of C_plus",
                "not_yet_output": (
                    "actual positive-order/auxiliary-dependent source path or direct public velocity(x,y,z,t)"
                ),
                "st006_comparison": (
                    "not directly comparable; this is a local source-structure derivative contract, not the ST006 full-domain PDE protocol"
                ),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceHomogeneousSpatialSensitivity":
        if not isinstance(payload, dict):
            raise ValueError("homogeneous spatial sensitivity payload must be an object")
        required = {"schema", "source", "parameters", "numerical_method", "routing", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("homogeneous spatial sensitivity payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported homogeneous spatial sensitivity schema")
        params = payload["parameters"]
        keys = {
            "derivative_label", "epsilon", "u_star", "lambda0", "c0", "L_s",
            "sigma", "min_tangential_norm", "reality_atol"
        }
        if not isinstance(params, dict) or set(params) != keys:
            raise ValueError("homogeneous spatial sensitivity parameters changed")
        obj = cls(**params)
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"homogeneous spatial sensitivity {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("homogeneous spatial sensitivity sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceHomogeneousSpatialSensitivity":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
