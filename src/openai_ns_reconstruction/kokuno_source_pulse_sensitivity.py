"""Parameter sensitivity contract for Kokuno's projected pulse path.

The corrected Kokuno reader supplies the projected nonzero-harmonic pulse ODE
implemented in :mod:`kokuno_source_pulse_inverse`.  This module differentiates
that *already sourced* equation with respect to one caller-declared slow or
transverse parameter ``xi``.  The sensitivity equation itself is an analytic
repository derivation, not an additional formula claimed from the reader.

Writing

    t' = M t + g,
    M = A_Phi - m^2 d I,
    g = -P f,
    P = I - n n^T / |n|^2,

its parameter derivative ``s=partial_xi t`` obeys

    s' = M s + M_xi t + g_xi,
    n dot s + n_xi dot t = 0.

The same derivative gives the source complete-curl coefficient sensitivity

    C_xi = i/k_m * [ (n_xi x t + n x s)/q
                     - (n x t) q_xi/q^2 ],
    q=|n|^2, q_xi=2 n dot n_xi.

When ``xi`` is later instantiated as the reader's source-normalized radial or
axial direction, this is the missing algebraic bridge to ``D_r C_m`` or
``D_z C_m`` required by the executable complete-curl contract.  This module
does not guess those actual source paths, forcing labels, background
corrections, or physical Q scaling.
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
SCHEMA = "kokuno-source-projected-pulse-sensitivity-v1"

_SOURCE_AND_DERIVED_FORMULAS = {
    "source_projected_equation": "t_m'=A_Phi*t_m-m^2*d*t_m-P_perp*f_m",
    "source_projected_matrix": (
        "A_Phi=-K+n_Phi*(n_Phi^T*K-(n_Phi')^T)/|n_Phi|^2"
    ),
    "source_diffusion": "d=epsilon*k^2*|n_Phi|^2",
    "derived_sensitivity_equation": "s_xi'=M*s_xi+M_xi*t_m+g_xi",
    "derived_constraint": "n_Phi dot s_xi+n_Phi_xi dot t_m=0",
    "derived_coefficient_sensitivity": (
        "C_m_xi=i/k_m*((n_xi cross t+n cross s_xi)/q-"
        "(n cross t)*(2*n dot n_xi)/q^2)"
    ),
}

_TRUTH_BOUNDARY = {
    "source_projected_pulse_equation_reused": True,
    "generic_parameter_sensitivity_operator_executable": True,
    "differentiated_transverse_constraint_executable": True,
    "complete_curl_coefficient_sensitivity_executable": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_D_r_path_instantiated": False,
    "source_actual_D_z_path_instantiated": False,
    "source_pulse_parameter_derivatives_reconstructed": False,
    "source_pulse_inverse_reconstructed": False,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "real_conjugate_pair_materialized": False,
    "physical_Q_scaled_velocity_materialized": False,
    "public_velocity_correction_materialized": False,
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


def _finite_complex(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=np.complex128)
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _vector3(value: Any, name: str, *, real: bool = False) -> np.ndarray:
    out = _finite_real(value, name) if real else _finite_complex(value, name)
    if out.ndim == 0 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing dimension 3")
    return out


def _matrix33(value: Any, name: str) -> np.ndarray:
    out = _finite_complex(value, name)
    if out.ndim < 2 or out.shape[-2:] != (3, 3):
        raise ValueError(f"{name} must have trailing shape (3,3)")
    return out


@dataclass(frozen=True)
class KokunoProjectedPulseSensitivityContract:
    """Differentiate a caller-supplied source projected-pulse path in ``xi``."""

    epsilon: float = 0.25
    m: int = 1
    min_covector_norm: float = 1.0e-12
    transverse_rtol: float = 2.0e-10
    max_projection_correction: float = 2.0e-2
    max_sensitivity_projection_correction: float = 5.0e-2
    max_nodes: int = 16385

    def __post_init__(self) -> None:
        pulse = KokunoProjectedPulseInverseContract(
            epsilon=self.epsilon,
            m=self.m,
            min_covector_norm=self.min_covector_norm,
            transverse_rtol=self.transverse_rtol,
            max_projection_correction=self.max_projection_correction,
            max_nodes=self.max_nodes,
        )
        object.__setattr__(self, "epsilon", pulse.epsilon)
        object.__setattr__(self, "m", pulse.m)
        object.__setattr__(self, "min_covector_norm", pulse.min_covector_norm)
        object.__setattr__(self, "transverse_rtol", pulse.transverse_rtol)
        object.__setattr__(self, "max_projection_correction", pulse.max_projection_correction)
        object.__setattr__(self, "max_nodes", pulse.max_nodes)
        guard = float(self.max_sensitivity_projection_correction)
        if not np.isfinite(guard) or not (0.0 < guard <= 0.5):
            raise ValueError("max_sensitivity_projection_correction must lie in (0,0.5]")
        object.__setattr__(self, "max_sensitivity_projection_correction", guard)

    @property
    def pulse(self) -> KokunoProjectedPulseInverseContract:
        return KokunoProjectedPulseInverseContract(
            epsilon=self.epsilon,
            m=self.m,
            min_covector_norm=self.min_covector_norm,
            transverse_rtol=self.transverse_rtol,
            max_projection_correction=self.max_projection_correction,
            max_nodes=self.max_nodes,
        )

    @property
    def k(self) -> int:
        return self.pulse.k

    @property
    def k_m(self) -> int:
        return self.pulse.k_m

    def _projector_data(self, n_phi: Any, n_phi_xi: Any) -> tuple[np.ndarray, ...]:
        n = _vector3(n_phi, "n_phi", real=True)
        n_xi = _vector3(n_phi_xi, "n_phi_xi", real=True)
        shape = np.broadcast_shapes(n.shape[:-1], n_xi.shape[:-1])
        n = np.broadcast_to(n, shape + (3,))
        n_xi = np.broadcast_to(n_xi, shape + (3,))
        q = np.sum(n * n, axis=-1)
        if np.any(q <= self.min_covector_norm**2):
            raise ValueError("n_phi norm is below the declared source-wave guard")
        q_xi = 2.0 * np.sum(n * n_xi, axis=-1)
        identity = np.broadcast_to(np.eye(3), shape + (3, 3))
        nn = n[..., :, None] * n[..., None, :]
        nxi_n = n_xi[..., :, None] * n[..., None, :]
        n_nxi = n[..., :, None] * n_xi[..., None, :]
        P = identity - nn / q[..., None, None]
        P_xi = (
            -(nxi_n + n_nxi) / q[..., None, None]
            + nn * q_xi[..., None, None] / (q * q)[..., None, None]
        )
        return n, n_xi, q, q_xi, P, P_xi

    def differentiated_operators(
        self,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        n_phi_xi: Any,
        n_phi_prime_xi: Any,
        K_xi: Any,
        forcing_xi: Any,
    ) -> dict[str, np.ndarray]:
        """Return ``M,g`` and their exact analytic ``xi`` derivatives.

        The source gives ``A_Phi,d,P``.  Their ``xi`` derivatives here follow
        by ordinary differentiation and are therefore repository-derived rather
        than quoted as an additional source identity.
        """
        n, n_xi, q, q_xi, P, P_xi = self._projector_data(n_phi, n_phi_xi)
        n_prime = _vector3(n_phi_prime, "n_phi_prime", real=True)
        n_prime_xi = _vector3(n_phi_prime_xi, "n_phi_prime_xi", real=True)
        K = _matrix33(K, "K")
        K_xi = _matrix33(K_xi, "K_xi")
        forcing = _vector3(forcing, "forcing")
        forcing_xi = _vector3(forcing_xi, "forcing_xi")
        shape = np.broadcast_shapes(
            n.shape[:-1],
            n_prime.shape[:-1],
            n_prime_xi.shape[:-1],
            K.shape[:-2],
            K_xi.shape[:-2],
            forcing.shape[:-1],
            forcing_xi.shape[:-1],
        )
        n = np.broadcast_to(n, shape + (3,))
        n_xi = np.broadcast_to(n_xi, shape + (3,))
        q = np.broadcast_to(q, shape)
        q_xi = np.broadcast_to(q_xi, shape)
        P = np.broadcast_to(P, shape + (3, 3))
        P_xi = np.broadcast_to(P_xi, shape + (3, 3))
        n_prime = np.broadcast_to(n_prime, shape + (3,))
        n_prime_xi = np.broadcast_to(n_prime_xi, shape + (3,))
        K = np.broadcast_to(K, shape + (3, 3))
        K_xi = np.broadcast_to(K_xi, shape + (3, 3))
        forcing = np.broadcast_to(forcing, shape + (3,))
        forcing_xi = np.broadcast_to(forcing_xi, shape + (3,))

        b = np.einsum("...i,...ij->...j", n, K) - n_prime
        b_xi = (
            np.einsum("...i,...ij->...j", n_xi, K)
            + np.einsum("...i,...ij->...j", n, K_xi)
            - n_prime_xi
        )
        nb = n[..., :, None] * b[..., None, :]
        A = -K + nb / q[..., None, None]
        A_xi = (
            -K_xi
            + (
                n_xi[..., :, None] * b[..., None, :]
                + n[..., :, None] * b_xi[..., None, :]
            )
            / q[..., None, None]
            - nb * q_xi[..., None, None] / (q * q)[..., None, None]
        )
        d = self.epsilon * (self.k**2) * q
        d_xi = self.epsilon * (self.k**2) * q_xi
        identity = np.broadcast_to(np.eye(3), shape + (3, 3))
        M = A - (self.m**2 * d)[..., None, None] * identity
        M_xi = A_xi - (self.m**2 * d_xi)[..., None, None] * identity
        g = -np.einsum("...ij,...j->...i", P, forcing)
        g_xi = -(
            np.einsum("...ij,...j->...i", P_xi, forcing)
            + np.einsum("...ij,...j->...i", P, forcing_xi)
        )
        return {
            "P": P,
            "P_xi": P_xi,
            "A_Phi": A,
            "A_Phi_xi": A_xi,
            "d": d,
            "d_xi": d_xi,
            "M": M,
            "M_xi": M_xi,
            "g": g,
            "g_xi": g_xi,
        }

    def sensitivity_rhs(
        self,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        n_phi_xi: Any,
        n_phi_prime_xi: Any,
        K_xi: Any,
        forcing_xi: Any,
        t_m: Any,
        t_m_xi: Any,
    ) -> np.ndarray:
        """Evaluate ``s'=M s+M_xi t+g_xi`` for ``s=partial_xi t_m``."""
        ops = self.differentiated_operators(
            n_phi,
            n_phi_prime,
            K,
            forcing,
            n_phi_xi,
            n_phi_prime_xi,
            K_xi,
            forcing_xi,
        )
        t = _vector3(t_m, "t_m")
        s = _vector3(t_m_xi, "t_m_xi")
        shape = np.broadcast_shapes(ops["M"].shape[:-2], t.shape[:-1], s.shape[:-1])
        M = np.broadcast_to(ops["M"], shape + (3, 3))
        M_xi = np.broadcast_to(ops["M_xi"], shape + (3, 3))
        g_xi = np.broadcast_to(ops["g_xi"], shape + (3,))
        t = np.broadcast_to(t, shape + (3,))
        s = np.broadcast_to(s, shape + (3,))
        return (
            np.einsum("...ij,...j->...i", M, s)
            + np.einsum("...ij,...j->...i", M_xi, t)
            + g_xi
        )

    def coefficient_sensitivity(
        self, n_phi: Any, n_phi_xi: Any, t_m: Any, t_m_xi: Any
    ) -> np.ndarray:
        """Return analytic ``partial_xi C_m`` for the source curl coefficient."""
        n, n_xi, q, q_xi, _, _ = self._projector_data(n_phi, n_phi_xi)
        t = _vector3(t_m, "t_m")
        s = _vector3(t_m_xi, "t_m_xi")
        shape = np.broadcast_shapes(n.shape[:-1], t.shape[:-1], s.shape[:-1])
        n = np.broadcast_to(n, shape + (3,))
        n_xi = np.broadcast_to(n_xi, shape + (3,))
        q = np.broadcast_to(q, shape)
        q_xi = np.broadcast_to(q_xi, shape)
        t = np.broadcast_to(t, shape + (3,))
        s = np.broadcast_to(s, shape + (3,))
        cross = np.cross(n, t)
        cross_xi = np.cross(n_xi, t) + np.cross(n, s)
        return (1j / self.k_m) * (
            cross_xi / q[..., None]
            - cross * q_xi[..., None] / (q * q)[..., None]
        )

    def _validate_path(
        self,
        path: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        n_phi_xi: Any,
        n_phi_prime_xi: Any,
        K_xi: Any,
        forcing_xi: Any,
    ) -> tuple[np.ndarray, ...]:
        path = _finite_real(path, "path")
        if path.ndim != 1 or path.size < 2 or path.size > self.max_nodes:
            raise ValueError(f"path must be a one-dimensional array with 2..{self.max_nodes} nodes")
        if np.any(np.diff(path) <= 0.0):
            raise ValueError("path nodes must be strictly increasing")
        N = path.size
        n = _vector3(n_phi, "n_phi", real=True)
        n_prime = _vector3(n_phi_prime, "n_phi_prime", real=True)
        K = _matrix33(K, "K")
        forcing = _vector3(forcing, "forcing")
        n_xi = _vector3(n_phi_xi, "n_phi_xi", real=True)
        n_prime_xi = _vector3(n_phi_prime_xi, "n_phi_prime_xi", real=True)
        K_xi = _matrix33(K_xi, "K_xi")
        forcing_xi = _vector3(forcing_xi, "forcing_xi")
        expected = {
            "n_phi": (N, 3),
            "n_phi_prime": (N, 3),
            "K": (N, 3, 3),
            "forcing": (N, 3),
            "n_phi_xi": (N, 3),
            "n_phi_prime_xi": (N, 3),
            "K_xi": (N, 3, 3),
            "forcing_xi": (N, 3),
        }
        actual = {
            "n_phi": n.shape,
            "n_phi_prime": n_prime.shape,
            "K": K.shape,
            "forcing": forcing.shape,
            "n_phi_xi": n_xi.shape,
            "n_phi_prime_xi": n_prime_xi.shape,
            "K_xi": K_xi.shape,
            "forcing_xi": forcing_xi.shape,
        }
        for name, shape in expected.items():
            if actual[name] != shape:
                raise ValueError(f"{name} must have shape {shape}")
        self._projector_data(n, n_xi)
        return path, n, n_prime, K, forcing, n_xi, n_prime_xi, K_xi, forcing_xi

    def solve_path(
        self,
        path: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        n_phi_xi: Any,
        n_phi_prime_xi: Any,
        K_xi: Any,
        forcing_xi: Any,
    ) -> dict[str, np.ndarray | float | str]:
        """Solve the zero-entrance pulse and one analytic parameter sensitivity.

        The RK4 coefficient interpolation matches the parent pulse solver.  At
        each accepted node, ``t`` is stabilized with ``P t_raw`` and the
        sensitivity uses the exact derivative of that numerical projection,
        ``P s_raw + P_xi t_raw``.  Both corrections are reported as autonomous
        numerical stabilization rather than source mathematics.
        """
        data = self._validate_path(
            path,
            n_phi,
            n_phi_prime,
            K,
            forcing,
            n_phi_xi,
            n_phi_prime_xi,
            K_xi,
            forcing_xi,
        )
        path, n, n_prime, K, forcing, n_xi, n_prime_xi, K_xi, forcing_xi = data
        N = path.size
        t = np.zeros((N, 3), dtype=np.complex128)
        s = np.zeros((N, 3), dtype=np.complex128)
        t_correction = np.zeros(N, dtype=float)
        s_correction = np.zeros(N, dtype=float)

        def pair_rhs(coeffs: tuple[np.ndarray, ...], y: np.ndarray, z: np.ndarray):
            nn, np_, KK, ff, nn_x, np_x, KK_x, ff_x = coeffs
            y_rhs = self.pulse.projected_rhs(nn, np_, KK, ff, y)
            z_rhs = self.sensitivity_rhs(
                nn, np_, KK, ff, nn_x, np_x, KK_x, ff_x, y, z
            )
            return y_rhs, z_rhs

        for j in range(N - 1):
            step = float(path[j + 1] - path[j])
            left = (
                n[j],
                n_prime[j],
                K[j],
                forcing[j],
                n_xi[j],
                n_prime_xi[j],
                K_xi[j],
                forcing_xi[j],
            )
            right = (
                n[j + 1],
                n_prime[j + 1],
                K[j + 1],
                forcing[j + 1],
                n_xi[j + 1],
                n_prime_xi[j + 1],
                K_xi[j + 1],
                forcing_xi[j + 1],
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

            _, _, _, _, P, P_xi = self._projector_data(n[j + 1], n_xi[j + 1])
            stable_t = P @ raw_t
            stable_s = P @ raw_s + P_xi @ raw_t
            tc = float(np.linalg.norm(stable_t - raw_t))
            sc = float(np.linalg.norm(stable_s - raw_s))
            if tc > self.max_projection_correction * max(1.0, float(np.linalg.norm(raw_t))):
                raise ValueError("pulse projection correction exceeded guard; refine path or fix inputs")
            if sc > self.max_sensitivity_projection_correction * max(
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
        C = np.zeros_like(t)
        C_xi = np.zeros_like(s)
        active = np.linalg.norm(t, axis=-1) > 0.0
        if np.any(active):
            C[active] = self.pulse.complete_curl_coefficients(n[active], t[active])
            C_xi[active] = self.coefficient_sensitivity(
                n[active], n_xi[active], t[active], s[active]
            )
        return {
            "path": path.copy(),
            "t_m": t,
            "t_m_xi": s,
            "C_m": C,
            "C_m_xi": C_xi,
            "pulse_constraint_defect_abs": pulse_constraint,
            "sensitivity_constraint_defect_abs": sensitivity_constraint,
            "pulse_projection_correction_norm": t_correction,
            "sensitivity_projection_correction_norm": s_correction,
            "max_pulse_constraint_defect_abs": float(np.max(pulse_constraint)),
            "max_sensitivity_constraint_defect_abs": float(np.max(sensitivity_constraint)),
            "max_pulse_projection_correction_norm": float(np.max(t_correction)),
            "max_sensitivity_projection_correction_norm": float(np.max(s_correction)),
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
                    "differentiate the sourced projected pulse equation with respect to one caller-declared parameter"
                ),
                "formulas": dict(_SOURCE_AND_DERIVED_FORMULAS),
                "derivation_boundary": (
                    "projected pulse equation is sourced; xi-sensitivity and C_m derivative are analytic repository derivations"
                ),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "m": self.m,
                "k": self.k,
                "k_m": self.k_m,
                "min_covector_norm": self.min_covector_norm,
                "transverse_rtol": self.transverse_rtol,
                "max_projection_correction": self.max_projection_correction,
                "max_sensitivity_projection_correction": self.max_sensitivity_projection_correction,
                "max_nodes": self.max_nodes,
            },
            "routing": {
                "caller_must_supply": (
                    "n_Phi,n_Phi_prime,K,f_m and their xi derivatives on one increasing source path"
                ),
                "future_binding": (
                    "instantiate xi as source-normalized D_r or D_z only after actual background/pulse path derivatives exist"
                ),
                "output": "t_m, partial_xi t_m, C_m, partial_xi C_m",
                "not_yet_output": "actual D_r C_m/D_z C_m, real conjugate pair, Q-scaled public [u,v,w]",
                "st006_comparison": (
                    "not directly comparable; this is a local manufactured sensitivity/refinement contract, not the ST006 full-domain residual protocol"
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoProjectedPulseSensitivityContract":
        if not isinstance(payload, dict):
            raise ValueError("pulse sensitivity payload must be a JSON object")
        expected = {"schema", "source", "parameters", "routing", "truth_boundary"}
        if set(payload) - {"sha256"} != expected:
            raise ValueError("pulse sensitivity payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported pulse sensitivity schema")
        parameters = payload["parameters"]
        allowed = {
            "epsilon",
            "m",
            "k",
            "k_m",
            "min_covector_norm",
            "transverse_rtol",
            "max_projection_correction",
            "max_sensitivity_projection_correction",
            "max_nodes",
        }
        if not isinstance(parameters, dict) or set(parameters) != allowed:
            raise ValueError("pulse sensitivity parameter keys do not match")
        obj = cls(
            epsilon=parameters["epsilon"],
            m=parameters["m"],
            min_covector_norm=parameters["min_covector_norm"],
            transverse_rtol=parameters["transverse_rtol"],
            max_projection_correction=parameters["max_projection_correction"],
            max_sensitivity_projection_correction=parameters[
                "max_sensitivity_projection_correction"
            ],
            max_nodes=parameters["max_nodes"],
        )
        expected_payload = obj.to_payload()
        if payload["source"] != expected_payload["source"]:
            raise ValueError("pulse sensitivity source/derivation metadata does not match")
        if payload["routing"] != expected_payload["routing"]:
            raise ValueError("pulse sensitivity routing metadata does not match")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("pulse sensitivity truth boundary does not match")
        if parameters["k"] != obj.k or parameters["k_m"] != obj.k_m:
            raise ValueError("pulse sensitivity derived wave numbers do not match")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("pulse sensitivity SHA256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoProjectedPulseSensitivityContract":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
