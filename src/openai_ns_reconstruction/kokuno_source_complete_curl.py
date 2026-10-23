"""Executable source stage-9 cylindrical complete-curl algebra.

This module implements only the exact single-harmonic curl identities recorded
in KokunoYumeto's corrected 2026-09-09 Navier--Stokes reconstruction reader.
It deliberately stops before reconstructing the imported pulse inverse that
produces the transverse amplitude ``t_m`` and before identifying Agent-1's
leading-only ``V/G`` bridge with the source's positive-order actual background.

For one nonzero harmonic ``m`` with ``k=ceil(epsilon**(-1/2))`` and
``k_m=k*m`` the source uses

    C_m = i (n_Phi x t_m) / (k_m |n_Phi|^2),
    curl_*(C_m exp(i k_m Phi)) = (t_m + r_m) exp(i k_m Phi),

where

    r_m = (-D_z C_theta,
           D_z C_r - D_r C_z,
           (D_r + R^-1) C_theta).

The caller supplies the source covector, a transverse amplitude, and the exact
coefficient derivatives ``D_r C_m`` / ``D_z C_m``.  In the source ``D_r`` may
include an auxiliary-torus derivative, so this module never silently replaces
it by a Cartesian or ordinary ``partial_R`` derivative.  The resulting object
is a complex normalized single harmonic; constructing the real conjugate pair
and applying the physical Q scaling remain separate later steps.
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
SCHEMA = "kokuno-source-complete-curl-contract-v1"

_SOURCE_FORMULAS = {
    "wave_number": "k=ceil(epsilon^(-1/2)); k_m=k*m, m integer nonzero",
    "coefficient": "C_m=i*(n_Phi x t_m)/(k_m*|n_Phi|^2)",
    "transverse": "n_Phi dot t_m=0",
    "complete_curl": "curl_*(C_m*exp(i*k_m*Phi))=(t_m+r_m)*exp(i*k_m*Phi)",
    "remainder": (
        "r_m=(-D_z(C_m)_theta, D_z(C_m)_r-D_r(C_m)_z, "
        "(D_r+R^(-1))(C_m)_theta)"
    ),
    "curl_star": (
        "(R^(-1)*d_theta A_z-D_z A_theta, D_z A_r-D_r A_z, "
        "(D_r+R^(-1))*A_theta-R^(-1)*d_theta A_r)"
    ),
    "physical_potential": "A_phys=Q^(1/2-A)*A_*",
}

_TRUTH_BOUNDARY = {
    "source_single_harmonic_coefficient_executable": True,
    "source_complete_curl_remainder_executable": True,
    "source_cylindrical_frame_terms_retained": True,
    "caller_supplies_transverse_pulse_amplitude": True,
    "caller_supplies_exact_coefficient_D_r_D_z": True,
    "source_pulse_inverse_reconstructed": False,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "real_conjugate_pair_materialized": False,
    "physical_Q_scaled_velocity_materialized": False,
    "public_velocity_correction_materialized": False,
    "complete_kokuno_composite_velocity": False,
    "phase_label_selected": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_real(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite real values")
    return array


def _finite_complex(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.complex128)
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _vector3(value: Any, name: str, *, real: bool = False) -> np.ndarray:
    array = _finite_real(value, name) if real else _finite_complex(value, name)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing vector dimension 3")
    return array


@dataclass(frozen=True)
class KokunoSourceCompleteCurlContract:
    """Exact normalized cylindrical complete-curl contract for one source mode."""

    epsilon: float = 0.25
    m: int = 1
    transverse_rtol: float = 2.0e-11
    min_covector_norm: float = 1.0e-12

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        m = int(self.m)
        transverse_rtol = float(self.transverse_rtol)
        min_covector_norm = float(self.min_covector_norm)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        if m != self.m or m == 0 or abs(m) > 32:
            raise ValueError("m must be a nonzero integer with |m|<=32")
        if not np.isfinite(transverse_rtol) or not (0.0 < transverse_rtol <= 1.0e-6):
            raise ValueError("transverse_rtol must lie in (0,1e-6]")
        if not np.isfinite(min_covector_norm) or not (0.0 < min_covector_norm <= 1.0e-4):
            raise ValueError("min_covector_norm must lie in (0,1e-4]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "transverse_rtol", transverse_rtol)
        object.__setattr__(self, "min_covector_norm", min_covector_norm)

    @property
    def k(self) -> int:
        return int(ceil(self.epsilon ** -0.5))

    @property
    def k_m(self) -> int:
        return self.k * self.m

    @classmethod
    def from_phase_contract(cls, phase_contract: Any, **kwargs: Any) -> "KokunoSourceCompleteCurlContract":
        """Bind ``epsilon,m`` to the executable source phase contract without copying it."""
        if not hasattr(phase_contract, "epsilon") or not hasattr(phase_contract, "m"):
            raise TypeError("phase_contract must expose epsilon and m")
        return cls(epsilon=float(phase_contract.epsilon), m=int(phase_contract.m), **kwargs)

    def _broadcast_n_t(self, n_phi: Any, t_m: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = _vector3(n_phi, "n_phi", real=True)
        t = _vector3(t_m, "t_m")
        shape = np.broadcast_shapes(n.shape[:-1], t.shape[:-1])
        n = np.broadcast_to(n, shape + (3,))
        t = np.broadcast_to(t, shape + (3,))
        n_norm_sq = np.sum(n * n, axis=-1)
        if np.any(n_norm_sq <= self.min_covector_norm**2):
            raise ValueError("n_phi norm is below the declared source-wave guard")
        t_norm = np.sqrt(np.sum(np.abs(t) ** 2, axis=-1))
        if np.any(t_norm <= 0.0):
            raise ValueError("t_m must be nonzero")
        transverse_defect = np.abs(np.sum(n * t, axis=-1))
        scale = np.sqrt(n_norm_sq) * t_norm
        if np.any(transverse_defect > self.transverse_rtol * scale):
            raise ValueError("t_m violates the source transverse condition n_Phi dot t_m=0")
        return n, t, n_norm_sq

    def coefficient(self, n_phi: Any, t_m: Any) -> np.ndarray:
        """Return exact source ``C_m=i(n_Phi x t_m)/(k_m |n_Phi|^2)``."""
        n, t, n_norm_sq = self._broadcast_n_t(n_phi, t_m)
        return 1j * np.cross(n, t) / (self.k_m * n_norm_sq[..., None])

    def phase_derivative_leading(self, n_phi: Any, coefficient: Any) -> np.ndarray:
        """Return the fast exponential contribution ``i*k_m*n_Phi x C_m``."""
        n = _vector3(n_phi, "n_phi", real=True)
        C = _vector3(coefficient, "coefficient")
        shape = np.broadcast_shapes(n.shape[:-1], C.shape[:-1])
        return 1j * self.k_m * np.cross(
            np.broadcast_to(n, shape + (3,)), np.broadcast_to(C, shape + (3,))
        )

    def remainder(self, R: Any, coefficient: Any, D_r_coefficient: Any, D_z_coefficient: Any) -> np.ndarray:
        """Return the exact derivative/frame remainder ``r_m``.

        ``D_r_coefficient`` and ``D_z_coefficient`` are source normalized
        derivatives.  They are deliberately caller-supplied because source
        ``D_r`` may contain an auxiliary-torus derivative and cannot in general
        be replaced by ordinary ``partial_R``.
        """
        R = _finite_real(R, "R")
        C = _vector3(coefficient, "coefficient")
        D_r_C = _vector3(D_r_coefficient, "D_r_coefficient")
        D_z_C = _vector3(D_z_coefficient, "D_z_coefficient")
        shape = np.broadcast_shapes(R.shape, C.shape[:-1], D_r_C.shape[:-1], D_z_C.shape[:-1])
        R = np.broadcast_to(R, shape)
        if np.any(R <= 0.0):
            raise ValueError("source cylindrical complete curl requires R>0")
        C = np.broadcast_to(C, shape + (3,))
        D_r_C = np.broadcast_to(D_r_C, shape + (3,))
        D_z_C = np.broadcast_to(D_z_C, shape + (3,))
        return np.stack(
            (
                -D_z_C[..., 1],
                D_z_C[..., 0] - D_r_C[..., 2],
                D_r_C[..., 1] + C[..., 1] / R,
            ),
            axis=-1,
        )

    def complete_amplitude(
        self,
        R: Any,
        n_phi: Any,
        t_m: Any,
        D_r_coefficient: Any,
        D_z_coefficient: Any,
    ) -> dict[str, np.ndarray]:
        """Return ``C_m``, ``r_m`` and the complete amplitude ``a_m=t_m+r_m``."""
        n, t, _ = self._broadcast_n_t(n_phi, t_m)
        C = self.coefficient(n, t)
        r = self.remainder(R, C, D_r_coefficient, D_z_coefficient)
        shape = np.broadcast_shapes(t.shape[:-1], r.shape[:-1])
        t = np.broadcast_to(t, shape + (3,))
        C = np.broadcast_to(C, shape + (3,))
        r = np.broadcast_to(r, shape + (3,))
        return {"C_m": C, "r_m": r, "a_m": t + r}

    def mode_vector_potential(self, phase: Any, coefficient: Any) -> np.ndarray:
        """Return normalized complex single-harmonic potential ``C_m exp(i k_m Phi)``."""
        phase = _finite_real(phase, "phase")
        C = _vector3(coefficient, "coefficient")
        shape = np.broadcast_shapes(phase.shape, C.shape[:-1])
        phase = np.broadcast_to(phase, shape)
        C = np.broadcast_to(C, shape + (3,))
        return C * np.exp(1j * self.k_m * phase)[..., None]

    def mode_velocity(
        self,
        R: Any,
        phase: Any,
        n_phi: Any,
        t_m: Any,
        D_r_coefficient: Any,
        D_z_coefficient: Any,
    ) -> dict[str, np.ndarray]:
        """Return the exact normalized complex source mode and its constituents."""
        data = self.complete_amplitude(R, n_phi, t_m, D_r_coefficient, D_z_coefficient)
        phase = _finite_real(phase, "phase")
        shape = np.broadcast_shapes(phase.shape, data["a_m"].shape[:-1])
        phase = np.broadcast_to(phase, shape)
        exponential = np.exp(1j * self.k_m * phase)
        C = np.broadcast_to(data["C_m"], shape + (3,))
        r = np.broadcast_to(data["r_m"], shape + (3,))
        a = np.broadcast_to(data["a_m"], shape + (3,))
        return {
            "C_m": C,
            "r_m": r,
            "a_m": a,
            "vector_potential": C * exponential[..., None],
            "velocity": a * exponential[..., None],
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
                "formula_scope": "stage-9 normalized cylindrical single-harmonic complete curl",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "m": self.m,
                "k": self.k,
                "k_m": self.k_m,
                "transverse_rtol": self.transverse_rtol,
                "min_covector_norm": self.min_covector_norm,
                "origin": "bounded numerical guards only; no hidden/source label recovery",
            },
            "caller_contract": {
                "required": [
                    "R>0",
                    "Phi",
                    "n_Phi",
                    "nonzero transverse t_m",
                    "exact source D_r C_m",
                    "exact source D_z C_m",
                ],
                "important": (
                    "D_r may include the source auxiliary-torus derivative; ordinary partial_R "
                    "is valid only for auxiliary-independent manufactured checks"
                ),
                "output": "complex normalized single harmonic; real conjugate pair and physical Q scaling pending",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceCompleteCurlContract":
        if not isinstance(payload, dict):
            raise ValueError("source complete-curl payload must be a JSON object")
        if set(payload) - {"sha256"} != {"schema", "source", "parameters", "caller_contract", "truth_boundary"}:
            raise ValueError("source complete-curl payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source complete-curl schema")
        parameters = payload["parameters"]
        expected_keys = {
            "epsilon", "m", "k", "k_m", "transverse_rtol", "min_covector_norm", "origin"
        }
        if set(parameters) != expected_keys:
            raise ValueError("source complete-curl parameter metadata changed")
        obj = cls(
            epsilon=float(parameters["epsilon"]),
            m=int(parameters["m"]),
            transverse_rtol=float(parameters["transverse_rtol"]),
            min_covector_norm=float(parameters["min_covector_norm"]),
        )
        expected = obj.to_payload()
        if payload["source"] != expected["source"]:
            raise ValueError("source complete-curl provenance metadata changed")
        if payload["caller_contract"] != expected["caller_contract"]:
            raise ValueError("source complete-curl caller contract changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("source complete-curl truth boundary changed")
        if parameters["k"] != obj.k or parameters["k_m"] != obj.k_m:
            raise ValueError("source complete-curl derived wave number changed")
        if parameters["origin"] != expected["parameters"]["origin"]:
            raise ValueError("source complete-curl parameter origin changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("source complete-curl payload SHA mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceCompleteCurlContract":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
