"""Compose caller-supplied Kokuno slow labels into a real localized curl family.

This module is deliberately below the public ``velocity(x,y,z,t)`` boundary.
It combines three source-facing contracts that are already executable in this
repository:

* the slow squared partition ``sum_beta eta_beta**2 = 1``;
* vector-potential localization before the complete cylindrical curl; and
* real ``m=+/-1`` harmonic reconstruction with physical ``Q**(-A)`` scaling.

For every slow label ``beta=(ell,a)`` the caller supplies the +1 source-mode
coefficient data and the corresponding partition values/derivatives.  The
module forms

    A_{beta,+} = eta_beta C_{beta,+} exp(i k Phi_beta),
    A_{beta,-} = conjugate(A_{beta,+}),

through the existing support-localized complete-curl implementation, then sums
the real conjugate pairs over beta.  The source rectangle sign ``sigma`` is not
the same object as the Fourier harmonic sign ``m=+/-1`` and must not be inserted
as an extra partition axis.

Concrete ``chi_ell`` / ``chi_{ell,a}`` bumps, actual source labels, actual
positive-order background data and the map from physical ``(x,y,z,t)`` to the
supplied mode data remain upstream.  Consequently this is a supplied-data
composition contract, not a paper-exact Kokuno field and not the missing public
second covariance column.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_source_partition_family import KokunoSourcePartitionFamilyContract
from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-localized-real-pair-family-v1"

_SOURCE_FORMULAS = {
    "slow_partition": "eta_beta=chi_ell(q)*chi_{ell,a}; sum_beta eta_beta^2=1",
    "localized_potential": "A_{beta,m}=eta_beta*C_{beta,m}*exp(i*k*m*Phi_beta)",
    "real_harmonics": "real waves are reconstructed from conjugate m=+/-1 harmonics",
    "physical_scaling": "u_phys=Q^(-A)*E(u_*), A=1/2+h",
    "slow_family_assembly": "source oscillatory field is assembled by summing localized label contributions",
}

_AUTONOMOUS_COMPOSITION = {
    "family_reality_guard": "check each supplied-data +/- harmonic sum against 2*Re(m=+1)",
    "partition_first": "validate the whole beta family before evaluating any localized curl",
    "per_label_Q_scaling": "apply caller-supplied Q_beta^(-A) before summing labels",
}

_TRUTH_BOUNDARY = {
    "source_partition_family_guard_reused": True,
    "source_vector_potential_localization_reused": True,
    "source_real_conjugate_harmonic_reconstruction_executable": True,
    "supplied_mode_family_composition_executable": True,
    "per_label_physical_Q_scaling_executable": True,
    "concrete_source_partition_bumps_reconstructed": False,
    "source_actual_partition_labels_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_auxiliary_dependent_D_r_path_instantiated": False,
    "positive_order_background_corrections_included": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "genuinely_independent_second_covariance_column_ready": False,
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
    if out.ndim < 2 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have a beta axis followed by vector dimension 3")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _real_vector3(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 2 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have a beta axis followed by vector dimension 3")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoSourceLocalizedRealPairFamily:
    """Real supplied-mode family assembled after whole-partition validation."""

    epsilon: float = 0.25
    h: float = 0.005
    partition_atol: float = 2.0e-12
    derivative_rtol: float = 2.0e-10
    reality_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        h = float(self.h)
        partition_atol = float(self.partition_atol)
        derivative_rtol = float(self.derivative_rtol)
        reality_atol = float(self.reality_atol)
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the source range 0<h<1/100")
        if not np.isfinite(partition_atol) or not (0.0 < partition_atol <= 1.0e-8):
            raise ValueError("partition_atol must lie in (0,1e-8]")
        if not np.isfinite(derivative_rtol) or not (0.0 < derivative_rtol <= 1.0e-6):
            raise ValueError("derivative_rtol must lie in (0,1e-6]")
        if not np.isfinite(reality_atol) or not (0.0 < reality_atol <= 1.0e-8):
            raise ValueError("reality_atol must lie in (0,1e-8]")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "partition_atol", partition_atol)
        object.__setattr__(self, "derivative_rtol", derivative_rtol)
        object.__setattr__(self, "reality_atol", reality_atol)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def partition_contract(self) -> KokunoSourcePartitionFamilyContract:
        return KokunoSourcePartitionFamilyContract(
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
        )

    @property
    def plus_localizer(self) -> KokunoSourceSupportLocalizedCurl:
        return KokunoSourceSupportLocalizedCurl(
            epsilon=self.epsilon, m=1, partition_atol=self.partition_atol
        )

    @property
    def minus_localizer(self) -> KokunoSourceSupportLocalizedCurl:
        return KokunoSourceSupportLocalizedCurl(
            epsilon=self.epsilon, m=-1, partition_atol=self.partition_atol
        )

    def _validated_inputs(
        self,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        family = self.partition_contract.validate_family(
            eta, D_r_eta, D_z_eta, beta_labels
        )
        eta_array = family["eta"]
        phase_array = _finite_real(phase, "phase")
        n_array = _real_vector3(n_phi, "n_phi")
        t_array = _complex_vector3(t_plus, "t_plus")
        Dr_array = _complex_vector3(D_r_C_plus, "D_r_C_plus")
        Dz_array = _complex_vector3(D_z_C_plus, "D_z_C_plus")

        expected_scalar = eta_array.shape
        expected_vector = expected_scalar + (3,)
        if phase_array.shape != expected_scalar:
            raise ValueError("phase must have exactly the eta shape, including the beta axis")
        for name, value in (
            ("n_phi", n_array),
            ("t_plus", t_array),
            ("D_r_C_plus", Dr_array),
            ("D_z_C_plus", Dz_array),
        ):
            if value.shape != expected_vector:
                raise ValueError(f"{name} must have shape eta.shape+(3,)")
        return {
            **family,
            "phase": phase_array,
            "n_phi": n_array,
            "t_plus": t_array,
            "D_r_C_plus": Dr_array,
            "D_z_C_plus": Dz_array,
        }

    def normalized_family(
        self,
        R: Any,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        """Return the real normalized localized family summed over beta.

        ``R`` may be scalar, sample-shaped, or beta-shaped.  It is broadcast to
        the partition family because source labels may use different normalized
        radial coordinates even at one physical sample.
        """
        data = self._validated_inputs(
            phase, n_phi, t_plus, D_r_C_plus, D_z_C_plus,
            eta, D_r_eta, D_z_eta, beta_labels,
        )
        family_shape = data["eta"].shape
        R_array = _finite_real(R, "R")
        try:
            R_family = np.broadcast_to(R_array, family_shape)
        except ValueError:
            # A sample-only radius naturally needs the final beta singleton.
            if R_array.shape == family_shape[:-1]:
                R_family = np.broadcast_to(R_array[..., None], family_shape)
            else:
                raise ValueError("R must broadcast to the eta family shape") from None

        plus = self.plus_localizer.localized_mode(
            R_family,
            data["phase"],
            data["n_phi"],
            data["t_plus"],
            data["D_r_C_plus"],
            data["D_z_C_plus"],
            data["eta"],
            data["D_r_eta"],
            data["D_z_eta"],
        )
        minus = self.minus_localizer.localized_mode(
            R_family,
            data["phase"],
            data["n_phi"],
            np.conjugate(data["t_plus"]),
            np.conjugate(data["D_r_C_plus"]),
            np.conjugate(data["D_z_C_plus"]),
            data["eta"],
            data["D_r_eta"],
            data["D_z_eta"],
        )
        pair_complex = plus["velocity"] + minus["velocity"]
        real_scale = np.maximum(1.0, np.max(np.abs(pair_complex.real)))
        if np.max(np.abs(pair_complex.imag)) > self.reality_atol * real_scale:
            raise RuntimeError("localized conjugate family failed the declared reality guard")
        pair_by_beta = pair_complex.real
        two_re_plus = 2.0 * plus["velocity"].real
        if not np.allclose(
            pair_by_beta, two_re_plus, rtol=0.0, atol=self.reality_atol * real_scale
        ):
            raise RuntimeError("localized m=+/-1 family disagrees with 2*Re(m=+1)")
        total = np.sum(pair_by_beta, axis=-2)
        return {
            "beta_labels": data["beta_labels"],
            "partition_error": data["partition_error"],
            "partition_radial_closure": data["radial_closure"],
            "partition_axial_closure": data["axial_closure"],
            "velocity_cylindrical_by_beta": pair_by_beta,
            "velocity_cylindrical_total": total,
            "plus_velocity_by_beta": plus["velocity"],
            "minus_velocity_by_beta": minus["velocity"],
            "plus_vector_potential_by_beta": plus["vector_potential"],
            "minus_vector_potential_by_beta": minus["vector_potential"],
        }

    def physical_family(
        self,
        Q: Any,
        R: Any,
        theta: Any,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        """Apply per-label Q scaling and return supplied-data physical velocity.

        ``Q`` is allowed to vary across beta because overlapping dyadic labels
        can carry different source scales.  ``theta`` is a physical sample
        angle and therefore has no beta axis.  This function still is not a
        public ``velocity(x,y,z,t)`` map: every source mode quantity is supplied
        explicitly by the caller.
        """
        result = self.normalized_family(
            R, phase, n_phi, t_plus, D_r_C_plus, D_z_C_plus,
            eta, D_r_eta, D_z_eta, beta_labels,
        )
        pair_by_beta = result["velocity_cylindrical_by_beta"]
        family_shape = pair_by_beta.shape[:-1]
        sample_shape = family_shape[:-1]

        Q_array = _finite_real(Q, "Q")
        try:
            Q_family = np.broadcast_to(Q_array, family_shape)
        except ValueError:
            if Q_array.shape == sample_shape:
                Q_family = np.broadcast_to(Q_array[..., None], family_shape)
            else:
                raise ValueError("Q must broadcast to the beta family shape") from None
        if np.any(Q_family <= 0.0) or np.any(Q_family > 1.0):
            raise ValueError("every Q_beta must lie in (0,1]")

        theta_array = _finite_real(theta, "theta")
        try:
            theta_sample = np.broadcast_to(theta_array, sample_shape)
        except ValueError:
            raise ValueError("theta must broadcast to sample axes only, without a beta axis") from None

        q_scale = Q_family ** (-self.A)
        physical_by_beta = q_scale[..., None] * pair_by_beta
        cyl_total = np.sum(physical_by_beta, axis=-2)
        c = np.cos(theta_sample)
        s = np.sin(theta_sample)
        cart_total = np.stack(
            (
                cyl_total[..., 0] * c - cyl_total[..., 1] * s,
                cyl_total[..., 0] * s + cyl_total[..., 1] * c,
                cyl_total[..., 2],
            ),
            axis=-1,
        )
        if not np.all(np.isfinite(cart_total)):
            raise RuntimeError("localized real-pair family produced a nonfinite physical velocity")
        return {
            **result,
            "Q_by_beta": Q_family,
            "A": np.full(sample_shape, self.A, dtype=float),
            "physical_scale_by_beta": q_scale,
            "velocity_physical_cylindrical_by_beta": physical_by_beta,
            "velocity_physical_cylindrical_total": cyl_total,
            "velocity_physical_cartesian_total": cart_total,
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
                "formula_scope": "slow-label localized complete curls assembled into supplied-data real harmonic family",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "autonomous_composition": dict(_AUTONOMOUS_COMPOSITION),
            "parameters": {
                "epsilon": self.epsilon,
                "h": self.h,
                "A": self.A,
                "partition_atol": self.partition_atol,
                "derivative_rtol": self.derivative_rtol,
                "reality_atol": self.reality_atol,
                "origin": "source ranges plus bounded floating-point guards; no source bump/background is invented",
            },
            "caller_contract": {
                "required": [
                    "one partition entry per beta=(ell,a), never duplicated by rectangle sigma",
                    "eta_beta and source-normalized D_r eta_beta, D_z eta_beta",
                    "per-beta Phi, n_Phi, t_plus, D_r C_plus, D_z C_plus",
                    "per-beta source Q and positive R plus a physical sample theta",
                ],
                "negative_harmonic": "constructed by exact complex conjugation at m=-1",
                "output": "real supplied-mode family; not a public velocity(x,y,z,t)",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceLocalizedRealPairFamily":
        if not isinstance(payload, dict):
            raise ValueError("localized real-pair family payload must be an object")
        required = {
            "schema", "source", "autonomous_composition", "parameters",
            "caller_contract", "truth_boundary",
        }
        if set(payload) - {"sha256"} != required:
            raise ValueError("localized real-pair family payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported localized real-pair family schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "epsilon", "h", "A", "partition_atol", "derivative_rtol",
            "reality_atol", "origin",
        }:
            raise ValueError("localized real-pair family parameters changed")
        obj = cls(
            epsilon=params["epsilon"],
            h=params["h"],
            partition_atol=params["partition_atol"],
            derivative_rtol=params["derivative_rtol"],
            reality_atol=params["reality_atol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"localized real-pair family {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("localized real-pair family sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceLocalizedRealPairFamily":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
