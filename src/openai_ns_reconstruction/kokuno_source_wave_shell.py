"""Fail-closed source wave-shell adapter for Kokuno localized complete curls.

The corrected 2026-09-09 reconstruction fixes two structural facts before the
oscillatory complete curl is used:

* the angular phase label lies on ``p=j/k`` with ``j`` a nonzero integer and
  ``k=ceil(epsilon**(-1/2))``; hence one harmonic has integer physical
  azimuthal index ``m*j``;
* on a dyadic chart ``X=R**2/(2*s_Q)``, with ``s_Q=q/Q`` in ``[1/2,2]``.
  The active annular wave region begins at ``X_a=4/Lambda``.  Therefore every
  active source wave point obeys ``R>=sqrt(X_a)`` and is separated from the
  cylindrical axis.

This adapter binds those source-domain contracts to the support-localized curl
implemented in :mod:`kokuno_source_support_localized_curl`.  It deliberately
does not invent the source's existentially chosen numerical ``Lambda``: the
caller supplies it, and the serialized truth boundary records that the actual
source value has not been recovered.  Likewise, this is only a lower annulus
/ axis-safety guard; it does not reconstruct the outer annulus edge, concrete
partition bumps, positive-order background, or a public ``velocity(x,y,z,t)``.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import sqrt
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-wave-shell-localized-curl-v1"

_SOURCE_FORMULAS = {
    "angular_grid": "p=j/k, j integer nonzero, k=ceil(epsilon^(-1/2))",
    "azimuthal_index": "k_m*p=(k*m)*(j/k)=m*j in Z\\{0}",
    "chart_profile_map": "X=R^2/(2*s_Q), s_Q=q/Q, 1/2<=s_Q<=2",
    "inner_annulus_edge": "X_a=4/Lambda",
    "axis_separation": "X>=X_a and s_Q>=1/2 imply R>=sqrt(X_a)>0",
    "localized_potential": "A_beta=eta_beta*C_m*exp(i*k_m*Phi)",
}

_TRUTH_BOUNDARY = {
    "source_angular_label_grid_enforced": True,
    "source_integer_azimuthal_index_executable": True,
    "source_inner_annulus_guard_executable": True,
    "source_axis_separation_guard_executable": True,
    "phase_covector_theta_consistency_guard_executable": True,
    "source_actual_Lambda_recovered": False,
    "source_outer_annulus_edge_reconstructed": False,
    "concrete_source_partition_bumps_reconstructed": False,
    "source_actual_background_path_instantiated": False,
    "positive_order_background_corrections_included": False,
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


def _vector3_real(value: Any, name: str) -> np.ndarray:
    out = _finite_real(value, name)
    if out.ndim == 0 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing vector dimension 3")
    return out


@dataclass(frozen=True)
class KokunoSourceWaveShellLocalizedCurl:
    """Bind source periodicity and the away-from-axis annulus to one curl lane.

    ``Lambda`` is the source large parameter entering ``X_a=4/Lambda``.  The
    corrected reader fixes that relation but chooses ``Lambda`` only after
    preceding constants are fixed; this class therefore requires a caller
    value and does not label any default/replay number as source exact.
    """

    Lambda: float
    shell_rtol: float = 2.0e-12
    covector_rtol: float = 2.0e-11

    def __post_init__(self) -> None:
        Lambda = float(self.Lambda)
        shell_rtol = float(self.shell_rtol)
        covector_rtol = float(self.covector_rtol)
        if not np.isfinite(Lambda) or not (1.0 <= Lambda <= 1.0e12):
            raise ValueError("Lambda must be finite in [1,1e12]")
        if not np.isfinite(shell_rtol) or not (0.0 < shell_rtol <= 1.0e-8):
            raise ValueError("shell_rtol must lie in (0,1e-8]")
        if not np.isfinite(covector_rtol) or not (0.0 < covector_rtol <= 1.0e-7):
            raise ValueError("covector_rtol must lie in (0,1e-7]")
        object.__setattr__(self, "Lambda", Lambda)
        object.__setattr__(self, "shell_rtol", shell_rtol)
        object.__setattr__(self, "covector_rtol", covector_rtol)

    @property
    def X_a(self) -> float:
        return 4.0 / self.Lambda

    @property
    def axis_lower_bound(self) -> float:
        return sqrt(self.X_a)

    @staticmethod
    def _check_phase_contract(phase_contract: Any) -> KokunoOscillatoryPhaseContract:
        if not isinstance(phase_contract, KokunoOscillatoryPhaseContract):
            raise TypeError("phase_contract must be KokunoOscillatoryPhaseContract")
        # The phase contract itself enforces p=j/k.  Recheck the derived
        # physical index here so the localized-curl seam fails closed if that
        # upstream contract changes in a future schema.
        index = phase_contract.m * phase_contract.j
        if index == 0:
            raise ValueError("source physical azimuthal index m*j must be nonzero")
        expected = phase_contract.k_m * phase_contract.p
        if not np.isclose(expected, index, rtol=0.0, atol=2.0e-12):
            raise ValueError("phase contract no longer satisfies k_m*p=m*j")
        return phase_contract

    @staticmethod
    def azimuthal_mode_number(phase_contract: Any) -> int:
        phase_contract = KokunoSourceWaveShellLocalizedCurl._check_phase_contract(
            phase_contract
        )
        return int(phase_contract.m * phase_contract.j)

    def source_shell_coordinates(self, R: Any, s_Q: Any) -> dict[str, np.ndarray]:
        """Return ``X=R^2/(2*s_Q)`` after the source lower-shell guards.

        This checks only the inner active-annulus edge.  The source outer edge
        is intentionally not guessed in this Agent-2 increment.
        """
        R, s_Q = np.broadcast_arrays(_finite_real(R, "R"), _finite_real(s_Q, "s_Q"))
        if np.any(R <= 0.0):
            raise ValueError("source wave shell requires R>0")
        if np.any((s_Q < 0.5) | (s_Q > 2.0)):
            raise ValueError("source dyadic chart requires 1/2<=s_Q<=2")
        X = R * R / (2.0 * s_Q)
        if np.any(X < self.X_a * (1.0 - self.shell_rtol)):
            raise ValueError("point lies below source active annulus X_a=4/Lambda")
        # Redundant algebraic guard, retained because this is the singularity
        # boundary actually needed by cylindrical complete-curl evaluation.
        if np.any(R < self.axis_lower_bound * (1.0 - self.shell_rtol)):
            raise ValueError("point violates source away-from-axis wave-shell bound")
        return {
            "R": np.asarray(R, dtype=float),
            "s_Q": np.asarray(s_Q, dtype=float),
            "X": np.asarray(X, dtype=float),
        }

    def _check_covector_theta(
        self, phase_contract: KokunoOscillatoryPhaseContract, R: np.ndarray, n_phi: Any
    ) -> np.ndarray:
        n = _vector3_real(n_phi, "n_phi")
        shape = np.broadcast_shapes(R.shape, n.shape[:-1])
        Rb = np.broadcast_to(R, shape)
        nb = np.broadcast_to(n, shape + (3,))
        expected = phase_contract.p / Rb
        scale = np.maximum(1.0, np.abs(expected))
        if np.any(np.abs(nb[..., 1] - expected) > self.covector_rtol * scale):
            raise ValueError("n_phi theta component is inconsistent with source p/R")
        return nb

    def localized_mode(
        self,
        phase_contract: Any,
        R: Any,
        s_Q: Any,
        phase: Any,
        n_phi: Any,
        t_m: Any,
        D_r_C_m: Any,
        D_z_C_m: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
    ) -> dict[str, np.ndarray | int | float]:
        """Evaluate #351 only after source periodicity and shell validation."""
        phase_contract = self._check_phase_contract(phase_contract)
        shell = self.source_shell_coordinates(R, s_Q)
        n = self._check_covector_theta(phase_contract, shell["R"], n_phi)
        localizer = KokunoSourceSupportLocalizedCurl(
            epsilon=phase_contract.epsilon,
            m=phase_contract.m,
        )
        out = localizer.localized_mode(
            shell["R"],
            phase,
            n,
            t_m,
            D_r_C_m,
            D_z_C_m,
            eta,
            D_r_eta,
            D_z_eta,
        )
        return {
            **out,
            "X": shell["X"],
            "s_Q": shell["s_Q"],
            "X_a": self.X_a,
            "axis_lower_bound": self.axis_lower_bound,
            "azimuthal_mode_number": self.azimuthal_mode_number(phase_contract),
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
                "formula_scope": "source angular-grid plus inner active-annulus guard before support-localized complete curl",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "Lambda": self.Lambda,
                "X_a": self.X_a,
                "axis_lower_bound": self.axis_lower_bound,
                "shell_rtol": self.shell_rtol,
                "covector_rtol": self.covector_rtol,
                "origin": (
                    "source X_a=4/Lambda relation with caller-supplied Lambda; "
                    "bounded floating-point guards are autonomous"
                ),
            },
            "caller_contract": {
                "required": [
                    "source-grid KokunoOscillatoryPhaseContract",
                    "R and s_Q=q/Q with 1/2<=s_Q<=2",
                    "n_Phi with theta component p/R",
                    "the #351 t_m, D_r C_m, D_z C_m, eta, D_r eta, D_z eta data",
                ],
                "lower_shell_only": (
                    "X>=X_a is enforced; the unreconstructed outer annulus and concrete support bumps remain caller/upstream work"
                ),
                "periodicity": "one fast harmonic has integer Cartesian azimuthal mode m*j",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceWaveShellLocalizedCurl":
        if not isinstance(payload, dict):
            raise ValueError("source wave-shell payload must be an object")
        required = {"schema", "source", "parameters", "caller_contract", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("source wave-shell payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source wave-shell schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "Lambda", "X_a", "axis_lower_bound", "shell_rtol", "covector_rtol", "origin"
        }:
            raise ValueError("source wave-shell parameters changed")
        obj = cls(
            Lambda=params["Lambda"],
            shell_rtol=params["shell_rtol"],
            covector_rtol=params["covector_rtol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"source wave-shell {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("source wave-shell sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceWaveShellLocalizedCurl":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
