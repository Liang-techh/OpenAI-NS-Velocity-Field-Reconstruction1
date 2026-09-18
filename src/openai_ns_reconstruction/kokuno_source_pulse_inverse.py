"""Numerical contract for Kokuno's projected nonzero-harmonic pulse equation.

This module implements the *source equation* recorded in the corrected
2026-09-09 Kokuno reconstruction reader, while keeping the discretization and
all still-missing source inputs explicit.

For one nonzero harmonic ``m`` the reader writes

    t_m' + K t_m + m^2 d t_m + i k_m n_Phi pi_m = -f_m,
    n_Phi . t_m = 0,
    d = epsilon k^2 |n_Phi|^2,

with

    pi_m = i/k_m * (n.K.t_m - n'.t_m + n.f_m) / |n|^2,

and therefore the projected transverse equation

    t_m' = A_Phi t_m - m^2 d t_m - proj_(n^perp) f_m,
    A_Phi = -K + n (n^T K - (n')^T) / |n|^2.

The source proves that zero transverse entrance data stay transverse.  The
implementation below solves a *caller-supplied tabulated path* with deterministic
RK4 and linear midpoint interpolation.  End-step orthogonal projection is an
autonomous numerical stabilization, not a source formula; its correction is
reported and guarded.  Consequently this module does not claim to reconstruct
the source's actual forcing labels, positive-order background, pulse family, or
paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import ceil
from pathlib import Path
from typing import Any

import numpy as np

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-projected-pulse-inverse-v1"

_SOURCE_FORMULAS = {
    "full_mode_equation": (
        "t_m'+K*t_m+m^2*d*t_m+i*k_m*n_Phi*pi_m=-f_m; n_Phi dot t_m=0"
    ),
    "diffusion": "d=epsilon*k^2*|n_Phi|^2",
    "pressure": (
        "pi_m=i/k_m*(n_Phi dot K*t_m-n_Phi' dot t_m+n_Phi dot f_m)/|n_Phi|^2"
    ),
    "projected_equation": (
        "t_m'=A_Phi*t_m-m^2*d*t_m-proj_(n_Phi^perp)f_m"
    ),
    "projected_matrix": (
        "A_Phi=-K+n_Phi*(n_Phi^T*K-(n_Phi')^T)/|n_Phi|^2"
    ),
    "constraint_evolution": "(n_Phi dot t_m)'=-m^2*d*(n_Phi dot t_m)",
}

_TRUTH_BOUNDARY = {
    "source_projected_pulse_equation_executable": True,
    "zero_entrance_tabulated_path_solver_executable": True,
    "source_pressure_recovery_executable": True,
    "source_complete_curl_coefficient_bridge_executable": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_pulse_parameter_derivatives_reconstructed": False,
    "source_pulse_inverse_reconstructed": False,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "real_conjugate_pair_materialized": False,
    "physical_Q_scaled_velocity_materialized": False,
    "public_velocity_correction_materialized": False,
    "complete_kokuno_composite_velocity": False,
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


@dataclass(frozen=True)
class KokunoProjectedPulseInverseContract:
    """Deterministic tabulated-path solver for the source projected pulse ODE."""

    epsilon: float = 0.25
    m: int = 1
    min_covector_norm: float = 1.0e-12
    transverse_rtol: float = 2.0e-10
    max_projection_correction: float = 2.0e-2
    max_nodes: int = 16385

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        m = int(self.m)
        min_covector_norm = float(self.min_covector_norm)
        transverse_rtol = float(self.transverse_rtol)
        max_projection_correction = float(self.max_projection_correction)
        max_nodes = int(self.max_nodes)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        if m != self.m or m == 0 or abs(m) > 32:
            raise ValueError("m must be a nonzero integer with |m|<=32")
        if not np.isfinite(min_covector_norm) or not (0.0 < min_covector_norm <= 1.0e-4):
            raise ValueError("min_covector_norm must lie in (0,1e-4]")
        if not np.isfinite(transverse_rtol) or not (0.0 < transverse_rtol <= 1.0e-5):
            raise ValueError("transverse_rtol must lie in (0,1e-5]")
        if not np.isfinite(max_projection_correction) or not (
            0.0 < max_projection_correction <= 0.25
        ):
            raise ValueError("max_projection_correction must lie in (0,0.25]")
        if max_nodes != self.max_nodes or not (2 <= max_nodes <= 100000):
            raise ValueError("max_nodes must be an integer in [2,100000]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "min_covector_norm", min_covector_norm)
        object.__setattr__(self, "transverse_rtol", transverse_rtol)
        object.__setattr__(self, "max_projection_correction", max_projection_correction)
        object.__setattr__(self, "max_nodes", max_nodes)

    @property
    def k(self) -> int:
        return int(ceil(self.epsilon ** -0.5))

    @property
    def k_m(self) -> int:
        return self.k * self.m

    def _guard_n(self, n_phi: Any) -> tuple[np.ndarray, np.ndarray]:
        n = _vector3(n_phi, "n_phi", real=True)
        n2 = np.sum(n * n, axis=-1)
        if np.any(n2 <= self.min_covector_norm**2):
            raise ValueError("n_phi norm is below the declared source-wave guard")
        return n, n2

    def project_transverse(self, n_phi: Any, value: Any) -> np.ndarray:
        """Return the bilinear orthogonal projection onto ``n_Phi^perp``."""
        n, n2 = self._guard_n(n_phi)
        value = _vector3(value, "value")
        shape = np.broadcast_shapes(n.shape[:-1], value.shape[:-1])
        n = np.broadcast_to(n, shape + (3,))
        value = np.broadcast_to(value, shape + (3,))
        n2 = np.broadcast_to(n2, shape)
        normal = np.sum(n * value, axis=-1) / n2
        return value - n * normal[..., None]

    def diffusion(self, n_phi: Any) -> np.ndarray:
        """Return source ``d=epsilon*k^2*|n_Phi|^2``."""
        _, n2 = self._guard_n(n_phi)
        return self.epsilon * (self.k**2) * n2

    def projected_matrix(self, n_phi: Any, n_phi_prime: Any, K: Any) -> np.ndarray:
        """Return source ``A_Phi`` in the projected transverse equation."""
        n, n2 = self._guard_n(n_phi)
        n_prime = _vector3(n_phi_prime, "n_phi_prime", real=True)
        K = _finite_complex(K, "K")
        if K.ndim < 2 or K.shape[-2:] != (3, 3):
            raise ValueError("K must have trailing shape (3,3)")
        shape = np.broadcast_shapes(n.shape[:-1], n_prime.shape[:-1], K.shape[:-2])
        n = np.broadcast_to(n, shape + (3,))
        n_prime = np.broadcast_to(n_prime, shape + (3,))
        K = np.broadcast_to(K, shape + (3, 3))
        n2 = np.broadcast_to(n2, shape)
        nTK = np.einsum("...i,...ij->...j", n, K)
        return -K + n[..., :, None] * (nTK - n_prime)[..., None, :] / n2[..., None, None]

    def projected_rhs(
        self,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        t_m: Any,
    ) -> np.ndarray:
        """Evaluate the exact source projected right-hand side."""
        t_m = _vector3(t_m, "t_m")
        A = self.projected_matrix(n_phi, n_phi_prime, K)
        forcing_perp = self.project_transverse(n_phi, forcing)
        d = self.diffusion(n_phi)
        shape = np.broadcast_shapes(A.shape[:-2], t_m.shape[:-1], forcing_perp.shape[:-1])
        A = np.broadcast_to(A, shape + (3, 3))
        t_m = np.broadcast_to(t_m, shape + (3,))
        forcing_perp = np.broadcast_to(forcing_perp, shape + (3,))
        d = np.broadcast_to(d, shape)
        return (
            np.einsum("...ij,...j->...i", A, t_m)
            - (self.m**2 * d)[..., None] * t_m
            - forcing_perp
        )

    def pressure(
        self,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        t_m: Any,
    ) -> np.ndarray:
        """Recover the source complex pressure amplitude ``pi_m``."""
        n, n2 = self._guard_n(n_phi)
        n_prime = _vector3(n_phi_prime, "n_phi_prime", real=True)
        K = _finite_complex(K, "K")
        forcing = _vector3(forcing, "forcing")
        t_m = _vector3(t_m, "t_m")
        if K.ndim < 2 or K.shape[-2:] != (3, 3):
            raise ValueError("K must have trailing shape (3,3)")
        shape = np.broadcast_shapes(
            n.shape[:-1], n_prime.shape[:-1], K.shape[:-2], forcing.shape[:-1], t_m.shape[:-1]
        )
        n = np.broadcast_to(n, shape + (3,))
        n_prime = np.broadcast_to(n_prime, shape + (3,))
        K = np.broadcast_to(K, shape + (3, 3))
        forcing = np.broadcast_to(forcing, shape + (3,))
        t_m = np.broadcast_to(t_m, shape + (3,))
        n2 = np.broadcast_to(n2, shape)
        Kt = np.einsum("...ij,...j->...i", K, t_m)
        numerator = (
            np.sum(n * Kt, axis=-1)
            - np.sum(n_prime * t_m, axis=-1)
            + np.sum(n * forcing, axis=-1)
        )
        return (1j / self.k_m) * numerator / n2

    def full_equation_residual(
        self,
        t_m_prime: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        t_m: Any,
    ) -> np.ndarray:
        """Evaluate the source full mode equation after pressure recovery."""
        t_m_prime = _vector3(t_m_prime, "t_m_prime")
        n, _ = self._guard_n(n_phi)
        K = _finite_complex(K, "K")
        forcing = _vector3(forcing, "forcing")
        t_m = _vector3(t_m, "t_m")
        pi_m = self.pressure(n, n_phi_prime, K, forcing, t_m)
        d = self.diffusion(n)
        shape = np.broadcast_shapes(
            t_m_prime.shape[:-1], n.shape[:-1], K.shape[:-2], forcing.shape[:-1], t_m.shape[:-1]
        )
        t_m_prime = np.broadcast_to(t_m_prime, shape + (3,))
        n = np.broadcast_to(n, shape + (3,))
        K = np.broadcast_to(K, shape + (3, 3))
        forcing = np.broadcast_to(forcing, shape + (3,))
        t_m = np.broadcast_to(t_m, shape + (3,))
        pi_m = np.broadcast_to(pi_m, shape)
        d = np.broadcast_to(d, shape)
        return (
            t_m_prime
            + np.einsum("...ij,...j->...i", K, t_m)
            + (self.m**2 * d)[..., None] * t_m
            + (1j * self.k_m * pi_m)[..., None] * n
            + forcing
        )

    def _validate_path(
        self,
        path: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        path = _finite_real(path, "path")
        if path.ndim != 1 or path.size < 2 or path.size > self.max_nodes:
            raise ValueError(f"path must be a one-dimensional array with 2..{self.max_nodes} nodes")
        if np.any(np.diff(path) <= 0.0):
            raise ValueError("path nodes must be strictly increasing")
        n, _ = self._guard_n(n_phi)
        n_prime = _vector3(n_phi_prime, "n_phi_prime", real=True)
        K = _finite_complex(K, "K")
        forcing = _vector3(forcing, "forcing")
        N = path.size
        if n.shape != (N, 3) or n_prime.shape != (N, 3):
            raise ValueError("n_phi and n_phi_prime must have shape (len(path),3)")
        if K.shape != (N, 3, 3):
            raise ValueError("K must have shape (len(path),3,3)")
        if forcing.shape != (N, 3):
            raise ValueError("forcing must have shape (len(path),3)")
        return path, n, n_prime, K, forcing

    def solve_path(
        self,
        path: Any,
        n_phi: Any,
        n_phi_prime: Any,
        K: Any,
        forcing: Any,
        *,
        t0: Any | None = None,
    ) -> dict[str, np.ndarray | float | str]:
        """Solve the projected source ODE on a supplied increasing path.

        Coefficients are supplied at path nodes.  RK4 uses their arithmetic
        midpoint values.  The resulting tabulated-coefficient method is
        typically second-order for smooth non-affine coefficient paths because
        the coefficient interpolation, not RK4, is the limiting approximation.

        The source starts the pulse from zero entrance data.  ``t0`` is exposed
        only for manufactured validation and must be transverse when nonzero.
        """
        path, n, n_prime, K, forcing = self._validate_path(
            path, n_phi, n_phi_prime, K, forcing
        )
        N = path.size
        t = np.zeros((N, 3), dtype=np.complex128)
        if t0 is not None:
            initial = _vector3(t0, "t0")
            if initial.shape != (3,):
                raise ValueError("t0 must have shape (3,)")
            initial_perp = self.project_transverse(n[0], initial)
            scale = max(1.0, float(np.linalg.norm(initial)))
            if float(np.linalg.norm(initial_perp - initial)) > self.transverse_rtol * scale:
                raise ValueError("t0 must satisfy n_Phi dot t0=0")
            t[0] = initial

        projection_correction = np.zeros(N, dtype=float)
        for j in range(N - 1):
            step = float(path[j + 1] - path[j])
            n_mid = 0.5 * (n[j] + n[j + 1])
            n_prime_mid = 0.5 * (n_prime[j] + n_prime[j + 1])
            K_mid = 0.5 * (K[j] + K[j + 1])
            forcing_mid = 0.5 * (forcing[j] + forcing[j + 1])

            y = t[j]
            k1 = self.projected_rhs(n[j], n_prime[j], K[j], forcing[j], y)
            k2 = self.projected_rhs(
                n_mid, n_prime_mid, K_mid, forcing_mid, y + 0.5 * step * k1
            )
            k3 = self.projected_rhs(
                n_mid, n_prime_mid, K_mid, forcing_mid, y + 0.5 * step * k2
            )
            k4 = self.projected_rhs(
                n[j + 1], n_prime[j + 1], K[j + 1], forcing[j + 1], y + step * k3
            )
            raw = y + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            stabilized = self.project_transverse(n[j + 1], raw)
            correction = float(np.linalg.norm(stabilized - raw))
            scale = max(1.0, float(np.linalg.norm(raw)))
            if correction > self.max_projection_correction * scale:
                raise ValueError(
                    "transverse projection correction exceeded guard; refine the path or fix source inputs"
                )
            projection_correction[j + 1] = correction
            t[j + 1] = stabilized

        constraint = np.abs(np.sum(n * t, axis=-1))
        pressure = self.pressure(n, n_prime, K, forcing, t)
        return {
            "path": path.copy(),
            "t_m": t,
            "pi_m": pressure,
            "constraint_defect_abs": constraint,
            "projection_correction_norm": projection_correction,
            "max_constraint_defect_abs": float(np.max(constraint)),
            "max_projection_correction_norm": float(np.max(projection_correction)),
            "integration_method": "RK4 with linear midpoint coefficient interpolation and guarded end-step transverse projection",
        }

    def complete_curl_coefficients(self, n_phi: Any, t_m: Any) -> np.ndarray:
        """Feed solved transverse amplitudes into Agent-2's source curl coefficient."""
        from .kokuno_source_complete_curl import KokunoSourceCompleteCurlContract

        curl = KokunoSourceCompleteCurlContract(
            epsilon=self.epsilon,
            m=self.m,
            transverse_rtol=max(self.transverse_rtol, 2.0e-11),
            min_covector_norm=self.min_covector_norm,
        )
        return curl.coefficient(n_phi, t_m)

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "nonzero-harmonic projected pulse equation and pressure recovery",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "m": self.m,
                "k": self.k,
                "k_m": self.k_m,
                "min_covector_norm": self.min_covector_norm,
                "transverse_rtol": self.transverse_rtol,
                "max_projection_correction": self.max_projection_correction,
                "max_nodes": self.max_nodes,
            },
            "numerical_method": {
                "path_representation": "caller-supplied increasing tabulated path",
                "stepper": "classical RK4",
                "midpoint_coefficients": "arithmetic interpolation of adjacent source-data nodes",
                "constraint_stabilization": "orthogonal projection onto n_Phi^perp after each accepted step",
                "origin": "autonomous numerical approximation; not a source hidden parameter or paper-exact discretization",
            },
            "routing": {
                "input_dependency": "caller must supply source-consistent n_Phi,n_Phi_prime,K,f_m path data",
                "output": "complex transverse t_m and pi_m; can feed C_m coefficient contract",
                "not_yet_output": "actual source pulse family, D_r/D_z coefficient sensitivities, real Q-scaled public velocity",
                "st006_comparison": "not directly comparable; this is a local ODE consistency contract, not the ST006 full-domain residual protocol",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoProjectedPulseInverseContract":
        if not isinstance(payload, dict):
            raise ValueError("projected pulse payload must be a JSON object")
        expected = {"schema", "source", "parameters", "numerical_method", "routing", "truth_boundary"}
        if set(payload) - {"sha256"} != expected:
            raise ValueError("projected pulse payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported projected pulse schema")
        if payload["source"] != cls().to_payload()["source"]:
            raise ValueError("projected pulse source provenance does not match pinned source")
        if payload["numerical_method"] != cls().to_payload()["numerical_method"]:
            raise ValueError("projected pulse numerical method metadata does not match")
        if payload["routing"] != cls().to_payload()["routing"]:
            raise ValueError("projected pulse routing metadata does not match")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("projected pulse truth boundary does not match")
        parameters = payload["parameters"]
        if not isinstance(parameters, dict):
            raise ValueError("projected pulse parameters must be an object")
        allowed = {
            "epsilon", "m", "k", "k_m", "min_covector_norm", "transverse_rtol",
            "max_projection_correction", "max_nodes"
        }
        if set(parameters) != allowed:
            raise ValueError("projected pulse parameter keys do not match")
        obj = cls(
            epsilon=parameters["epsilon"],
            m=parameters["m"],
            min_covector_norm=parameters["min_covector_norm"],
            transverse_rtol=parameters["transverse_rtol"],
            max_projection_correction=parameters["max_projection_correction"],
            max_nodes=parameters["max_nodes"],
        )
        if parameters["k"] != obj.k or parameters["k_m"] != obj.k_m:
            raise ValueError("projected pulse derived wave numbers do not match")
        unsigned = obj.to_payload()
        if "sha256" in payload and payload["sha256"] != hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest():
            raise ValueError("projected pulse SHA256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoProjectedPulseInverseContract":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
