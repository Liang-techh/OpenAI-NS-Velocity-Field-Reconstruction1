"""Materialize Kokuno's real m=±1 complete-curl pair and physical velocity.

The corrected 2026-09-09 reader states that real waves are formed from
conjugate nonzero harmonics.  If ``t_+`` is the m=+1 Fourier coefficient,
then ``t_- = conj(t_+)`` and the complete-curl potentials, remainders and
velocities are conjugate as well.  The normalized real field is therefore

    u_* = u_+ + u_- = 2 Re(u_+).

The same reader gives ``u_phys = Q^{-A} E u_*`` with
``A = 1/2 + h``.  This module executes only those identities and the ordinary
cylindrical-to-Cartesian frame rotation.  It consumes already-computed source
``Phi, n_Phi, t_+, D_r C_+, D_z C_+``; it does not invent the actual supported
pulse forcing, positive-order background, auxiliary-dependent D_r path, or a
direct ``velocity(x,y,z,t)`` map.

Source formulas and the repository's autonomous numerical/input guards are
kept separate in the serialized truth boundary.  In particular this is not a
paper-exact or OpenAI-identified field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_complete_curl import KokunoSourceCompleteCurlContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-real-conjugate-pair-v1"

_SOURCE_FORMULAS = {
    "real_pair": (
        "real waves use conjugate harmonic pairs; for real t*cos(k*Phi), "
        "the Fourier velocity coefficients are t/2 at m=+1 and m=-1"
    ),
    "single_mode": "curl_*(C_m*exp(i*k*m*Phi))=(t_m+r_m)*exp(i*k*m*Phi)",
    "coefficient": "C_m=i*(n_Phi x t_m)/(k*m*|n_Phi|^2)",
    "physical_velocity": "u_phys=Q^(-A)*E(u_*), A=1/2+h",
    "cylindrical_frame": (
        "u_x=u_r*cos(theta)-u_theta*sin(theta); "
        "u_y=u_r*sin(theta)+u_theta*cos(theta); u_z=u_z"
    ),
}

_TRUTH_BOUNDARY = {
    "source_real_conjugate_pair_identity_executable": True,
    "source_physical_Q_velocity_scaling_executable": True,
    "cylindrical_to_cartesian_rotation_executable": True,
    "real_conjugate_pair_materialized_from_caller_mode_data": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_auxiliary_dependent_D_r_path_instantiated": False,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "public_xyz_t_velocity_correction_materialized": False,
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


def _complex_vector3(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=np.complex128)
    if out.ndim == 0 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing vector dimension 3")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoSourceRealConjugatePair:
    """Real m=±1 materializer for already reconstructed source mode data."""

    epsilon: float = 0.25
    h: float = 0.005
    reality_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        h = float(self.h)
        reality_atol = float(self.reality_atol)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        # The completed source construction fixes 0<h<1/100.
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the source range 0<h<1/100")
        if not np.isfinite(reality_atol) or not (0.0 < reality_atol <= 1.0e-8):
            raise ValueError("reality_atol must lie in (0,1e-8]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "reality_atol", reality_atol)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def plus_contract(self) -> KokunoSourceCompleteCurlContract:
        return KokunoSourceCompleteCurlContract(epsilon=self.epsilon, m=1)

    @property
    def minus_contract(self) -> KokunoSourceCompleteCurlContract:
        return KokunoSourceCompleteCurlContract(epsilon=self.epsilon, m=-1)

    def normalized_pair(
        self,
        R: Any,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
    ) -> dict[str, np.ndarray]:
        """Return the normalized real m=±1 complete-curl pair.

        ``t_plus`` is the +1 Fourier coefficient, not the cosine amplitude.
        For a real cosine amplitude ``t`` the source convention is
        ``t_plus=t_minus=t/2``.  The negative-mode coefficient derivatives are
        the complex conjugates of the supplied positive-mode derivatives.
        """
        t_plus = _complex_vector3(t_plus, "t_plus")
        D_r_C_plus = _complex_vector3(D_r_C_plus, "D_r_C_plus")
        D_z_C_plus = _complex_vector3(D_z_C_plus, "D_z_C_plus")

        plus = self.plus_contract.mode_velocity(
            R, phase, n_phi, t_plus, D_r_C_plus, D_z_C_plus
        )
        minus = self.minus_contract.mode_velocity(
            R,
            phase,
            n_phi,
            np.conjugate(t_plus),
            np.conjugate(D_r_C_plus),
            np.conjugate(D_z_C_plus),
        )
        pair_complex = plus["velocity"] + minus["velocity"]
        scale = np.maximum(1.0, np.max(np.abs(pair_complex.real)))
        if np.max(np.abs(pair_complex.imag)) > self.reality_atol * scale:
            raise RuntimeError("conjugate source modes failed the declared reality guard")
        velocity = pair_complex.real
        two_re_plus = 2.0 * plus["velocity"].real
        if not np.allclose(velocity, two_re_plus, rtol=0.0, atol=self.reality_atol * scale):
            raise RuntimeError("m=±1 sum disagrees with 2*Re(m=+1) identity")
        return {
            "velocity_cylindrical": velocity,
            "plus_velocity": plus["velocity"],
            "minus_velocity": minus["velocity"],
            "plus_vector_potential": plus["vector_potential"],
            "minus_vector_potential": minus["vector_potential"],
            "plus_complete_amplitude": plus["a_m"],
            "minus_complete_amplitude": minus["a_m"],
        }

    def physical_pair(
        self,
        Q: Any,
        R: Any,
        theta: Any,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
    ) -> dict[str, np.ndarray]:
        """Return source-scaled physical cylindrical and Cartesian velocity.

        This is a pointwise materializer from source mode data, not yet a
        public ``velocity(x,y,z,t)`` because the actual source forcing,
        background and label-to-physical coordinate assembly remain upstream.
        """
        Q = _finite_real(Q, "Q")
        theta = _finite_real(theta, "theta")
        if np.any(Q <= 0.0) or np.any(Q > 1.0):
            raise ValueError("Q must lie in (0,1] for the source dyadic bands")

        pair = self.normalized_pair(R, phase, n_phi, t_plus, D_r_C_plus, D_z_C_plus)
        u_star = pair["velocity_cylindrical"]
        shape = np.broadcast_shapes(Q.shape, theta.shape, u_star.shape[:-1])
        Q = np.broadcast_to(Q, shape)
        theta = np.broadcast_to(theta, shape)
        u_star = np.broadcast_to(u_star, shape + (3,))

        q_scale = Q ** (-self.A)
        cyl = q_scale[..., None] * u_star
        c = np.cos(theta)
        s = np.sin(theta)
        cart = np.stack(
            (
                cyl[..., 0] * c - cyl[..., 1] * s,
                cyl[..., 0] * s + cyl[..., 1] * c,
                cyl[..., 2],
            ),
            axis=-1,
        )
        if not np.all(np.isfinite(cart)):
            raise RuntimeError("physical real-pair evaluation produced a nonfinite result")
        return {
            **pair,
            "Q": Q,
            "A": np.full(shape, self.A, dtype=float),
            "physical_scale": q_scale,
            "velocity_physical_cylindrical": cyl,
            "velocity_physical_cartesian": cart,
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
                "formula_scope": "stage-9 real conjugate pair plus physical Q scaling",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "h": self.h,
                "A": self.A,
                "reality_atol": self.reality_atol,
                "origin": "source h-range plus a bounded autonomous floating-point reality guard",
            },
            "caller_contract": {
                "required": [
                    "Q in (0,1]",
                    "R>0",
                    "theta",
                    "Phi",
                    "n_Phi",
                    "m=+1 Fourier coefficient t_plus",
                    "source-normalized D_r C_plus",
                    "source-normalized D_z C_plus",
                ],
                "negative_mode": "constructed by exact complex conjugation with m=-1",
                "output": "real physical velocity from supplied source-mode data; not yet velocity(x,y,z,t)",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceRealConjugatePair":
        if not isinstance(payload, dict):
            raise ValueError("real-pair payload must be an object")
        required = {"schema", "source", "parameters", "caller_contract", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("real-pair payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported real-pair schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "epsilon", "h", "A", "reality_atol", "origin"
        }:
            raise ValueError("real-pair parameters changed")
        obj = cls(
            epsilon=params["epsilon"],
            h=params["h"],
            reality_atol=params["reality_atol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"real-pair {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("real-pair sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceRealConjugatePair":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
