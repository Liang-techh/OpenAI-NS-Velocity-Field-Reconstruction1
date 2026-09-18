"""Localize Kokuno source harmonics at vector-potential / complete-curl level.

The corrected 2026-09-09 reconstruction assembles slow labels with the squared
partition

    eta_beta = chi_ell(q) * chi_{ell,a},     sum_beta eta_beta**2 = 1.

A source oscillatory wave must be localized before taking its curl.  For one
harmonic with complete-curl coefficient ``C_m`` this means

    A_beta = eta_beta * C_m * exp(i*k_m*Phi).

Because the slow support weight is independent of the fast angular phase, the
source-normalized product derivatives are exactly

    D_r(eta C_m) = eta D_r C_m + (D_r eta) C_m,
    D_z(eta C_m) = eta D_z C_m + (D_z eta) C_m.

Feeding those product derivatives into the full cylindrical curl retains the
cutoff-gradient terms.  In particular this module never forms a divergence-free
velocity first and then multiplies that velocity by a cutoff.

The concrete partition functions ``chi_ell`` and ``chi_{ell,a}``, and their
source-normalized derivatives, are deliberately caller supplied: the retained
corrected-reader statement specifies the squared-partition structure but does
not supply a unique numerical bump implementation here.  No autonomous cutoff
is therefore mislabeled as Kokuno exact.
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
SCHEMA = "kokuno-source-support-localized-curl-v1"

_SOURCE_FORMULAS = {
    "squared_partition": (
        "eta_beta=chi_ell(q)*chi_{ell,a}; sum_beta eta_beta^2=1, "
        "with duplicated signs retained as distinct rectangle labels"
    ),
    "localized_potential": "A_beta=eta_beta*C_m*exp(i*k_m*Phi)",
    "radial_product_rule": "D_r(eta*C_m)=eta*D_r(C_m)+(D_r eta)*C_m",
    "axial_product_rule": "D_z(eta*C_m)=eta*D_z(C_m)+(D_z eta)*C_m",
    "complete_curl": "u_beta=curl_*(A_beta), retaining cylindrical frame terms",
}

_TRUTH_BOUNDARY = {
    "source_squared_partition_structure_executable": True,
    "vector_potential_level_localization_executable": True,
    "support_gradient_complete_curl_terms_retained": True,
    "caller_supplies_partition_weight_eta": True,
    "caller_supplies_source_normalized_D_r_eta_D_z_eta": True,
    "concrete_source_partition_bumps_reconstructed": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_auxiliary_dependent_D_r_path_instantiated": False,
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


def _vector3(value: Any, name: str, *, real: bool = False) -> np.ndarray:
    dtype = float if real else np.complex128
    out = np.asarray(value, dtype=dtype)
    if out.ndim == 0 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing vector dimension 3")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoSourceSupportLocalizedCurl:
    """Apply one caller-supplied source slow support weight before curling."""

    epsilon: float = 0.25
    m: int = 1
    partition_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        m = int(self.m)
        partition_atol = float(self.partition_atol)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        if m != self.m or m == 0 or abs(m) > 32:
            raise ValueError("m must be a nonzero integer with |m|<=32")
        if not np.isfinite(partition_atol) or not (0.0 < partition_atol <= 1.0e-8):
            raise ValueError("partition_atol must lie in (0,1e-8]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "partition_atol", partition_atol)

    @property
    def complete_curl(self) -> KokunoSourceCompleteCurlContract:
        return KokunoSourceCompleteCurlContract(epsilon=self.epsilon, m=self.m)

    def localized_coefficient_data(
        self,
        n_phi: Any,
        t_m: Any,
        D_r_C_m: Any,
        D_z_C_m: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
    ) -> dict[str, np.ndarray]:
        """Return exact product-rule data for ``eta*C_m``.

        ``D_r_eta`` and ``D_z_eta`` must use the same normalized source
        derivatives as ``D_r_C_m`` and ``D_z_C_m``.  This is essential when
        ``D_r`` contains the source auxiliary-torus directional derivative.
        """
        contract = self.complete_curl
        n = _vector3(n_phi, "n_phi", real=True)
        t = _vector3(t_m, "t_m")
        C = contract.coefficient(n, t)
        D_r_C = _vector3(D_r_C_m, "D_r_C_m")
        D_z_C = _vector3(D_z_C_m, "D_z_C_m")
        eta = _finite_real(eta, "eta")
        D_r_eta = _finite_real(D_r_eta, "D_r_eta")
        D_z_eta = _finite_real(D_z_eta, "D_z_eta")

        shape = np.broadcast_shapes(
            C.shape[:-1], D_r_C.shape[:-1], D_z_C.shape[:-1],
            eta.shape, D_r_eta.shape, D_z_eta.shape,
        )
        C = np.broadcast_to(C, shape + (3,))
        t = np.broadcast_to(t, shape + (3,))
        D_r_C = np.broadcast_to(D_r_C, shape + (3,))
        D_z_C = np.broadcast_to(D_z_C, shape + (3,))
        eta = np.broadcast_to(eta, shape)
        D_r_eta = np.broadcast_to(D_r_eta, shape)
        D_z_eta = np.broadcast_to(D_z_eta, shape)

        # From sum_beta eta_beta^2=1 every real partition factor has |eta|<=1.
        # The small tolerance is only a floating-point guard.
        if np.any(np.abs(eta) > 1.0 + self.partition_atol):
            raise ValueError("eta is incompatible with a real squared partition: |eta| must be <=1")

        C_local = eta[..., None] * C
        D_r_C_local = eta[..., None] * D_r_C + D_r_eta[..., None] * C
        D_z_C_local = eta[..., None] * D_z_C + D_z_eta[..., None] * C
        t_local = eta[..., None] * t
        return {
            "C_m": C,
            "t_m": t,
            "eta": eta,
            "D_r_eta": D_r_eta,
            "D_z_eta": D_z_eta,
            "C_local": C_local,
            "t_local": t_local,
            "D_r_C_local": D_r_C_local,
            "D_z_C_local": D_z_C_local,
        }

    def localized_mode(
        self,
        R: Any,
        phase: Any,
        n_phi: Any,
        t_m: Any,
        D_r_C_m: Any,
        D_z_C_m: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
    ) -> dict[str, np.ndarray]:
        """Return ``curl_*(eta*C_m*exp(i*k_m*Phi))`` exactly algebraically.

        The leading fast-phase term is computed from the localized coefficient
        instead of calling ``mode_velocity`` with ``eta*t_m``.  This permits
        exact support zeros while retaining the base mode's transversality
        validation.
        """
        contract = self.complete_curl
        data = self.localized_coefficient_data(
            n_phi, t_m, D_r_C_m, D_z_C_m, eta, D_r_eta, D_z_eta
        )
        R = _finite_real(R, "R")
        phase = _finite_real(phase, "phase")
        leading = contract.phase_derivative_leading(n_phi, data["C_local"])
        remainder = contract.remainder(
            R, data["C_local"], data["D_r_C_local"], data["D_z_C_local"]
        )
        shape = np.broadcast_shapes(
            phase.shape, data["C_local"].shape[:-1], leading.shape[:-1], remainder.shape[:-1]
        )
        phase = np.broadcast_to(phase, shape)
        C_local = np.broadcast_to(data["C_local"], shape + (3,))
        leading = np.broadcast_to(leading, shape + (3,))
        remainder = np.broadcast_to(remainder, shape + (3,))
        exponential = np.exp(1j * contract.k_m * phase)
        complete_amplitude = leading + remainder
        return {
            **data,
            "leading_local": leading,
            "remainder_local": remainder,
            "a_local": complete_amplitude,
            "vector_potential": C_local * exponential[..., None],
            "velocity": complete_amplitude * exponential[..., None],
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
                "formula_scope": "slow squared-partition support applied before stage-9 complete curl",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "epsilon": self.epsilon,
                "m": self.m,
                "partition_atol": self.partition_atol,
                "origin": "source support algebra plus a bounded floating-point partition guard",
            },
            "caller_contract": {
                "required": [
                    "source n_Phi and nonzero transverse t_m",
                    "source-normalized D_r C_m and D_z C_m",
                    "real eta_beta=chi_ell*chi_{ell,a}",
                    "source-normalized D_r eta_beta and D_z eta_beta",
                ],
                "important": (
                    "localization is applied to the vector potential coefficient; multiplying the "
                    "already-curled velocity by eta would omit cutoff-gradient curl terms"
                ),
                "partition_realization": (
                    "caller supplied; this module does not invent chi_ell or chi_{ell,a}"
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceSupportLocalizedCurl":
        if not isinstance(payload, dict):
            raise ValueError("support-localized curl payload must be an object")
        required = {"schema", "source", "parameters", "caller_contract", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("support-localized curl payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported support-localized curl schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "epsilon", "m", "partition_atol", "origin"
        }:
            raise ValueError("support-localized curl parameters changed")
        obj = cls(
            epsilon=params["epsilon"],
            m=params["m"],
            partition_atol=params["partition_atol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"support-localized curl {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("support-localized curl sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceSupportLocalizedCurl":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
