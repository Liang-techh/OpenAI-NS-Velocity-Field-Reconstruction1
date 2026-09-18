"""Fail-closed structural contract for Kokuno slow squared partitions.

The corrected 2026-09-09 reconstruction uses slow labels ``beta=(ell,a)`` with

    eta_beta = chi_ell(q) * chi_{ell,a},
    sum_beta eta_beta**2 = 1.

The duplicated signs belong to rectangle labels ``gamma=(ell,a,sigma)`` but
are *not* counted twice in the squared-partition sum.  This module makes those
source bookkeeping facts executable without inventing concrete numerical
``chi_ell`` / ``chi_{ell,a}`` bumps.

If the caller also supplies source-normalized derivatives, differentiating the
squared partition gives

    sum_beta eta_beta * D_r eta_beta = 0,
    sum_beta eta_beta * D_z eta_beta = 0.

Those two derivative identities are elementary algebraic consequences of the
source partition, not separately quoted source formulas.  They are useful
fail-closed guards before the per-label product-rule complete curl in
``kokuno_source_support_localized_curl``.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-partition-family-v1"

_SOURCE_FORMULAS = {
    "slow_label": "beta=(ell,a)",
    "rectangle_label": "gamma=(ell,a,sigma), with duplicated signs assigned after the beta partition",
    "partition_weight": "eta_beta=chi_ell(q)*chi_{ell,a}",
    "squared_partition": "sum_beta eta_beta^2=1",
    "dyadic_overlap": "overlapping dyadic supports imply |ell-ell'|<=2",
    "product_grid_mesh": "chi_{ell,a} comes from a product-grid squared partition of mesh S_*^(-3)",
}

_AUTONOMOUS_DERIVATIONS = {
    "radial_closure": "D_r(sum_beta eta_beta^2)=2*sum_beta eta_beta*D_r eta_beta=0",
    "axial_closure": "D_z(sum_beta eta_beta^2)=2*sum_beta eta_beta*D_z eta_beta=0",
}

_TRUTH_BOUNDARY = {
    "source_squared_partition_family_structure_executable": True,
    "source_beta_vs_gamma_counting_guard_executable": True,
    "derived_partition_D_r_closure_guard_executable": True,
    "derived_partition_D_z_closure_guard_executable": True,
    "concrete_source_partition_bumps_reconstructed": False,
    "source_actual_S_star_recovered": False,
    "source_actual_partition_labels_instantiated": False,
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


def _normalize_beta_labels(labels: Sequence[Any], n_beta: int) -> tuple[tuple[Any, Any], ...]:
    if isinstance(labels, (str, bytes)):
        raise TypeError("beta_labels must be a sequence of beta=(ell,a) labels")
    raw = list(labels)
    if len(raw) != n_beta:
        raise ValueError("beta_labels length must match the partition axis")
    normalized: list[tuple[Any, Any]] = []
    keys: list[str] = []
    for label in raw:
        if isinstance(label, (str, bytes)):
            raise ValueError("each partition label must be beta=(ell,a), not a string")
        try:
            parts = tuple(label)
        except TypeError as exc:
            raise ValueError("each partition label must be beta=(ell,a)") from exc
        if len(parts) != 2:
            raise ValueError(
                "partition sum is over beta=(ell,a); do not duplicate gamma=(ell,a,sigma) signs"
            )
        ell, a = parts
        if isinstance(ell, (bool, np.bool_)) or not isinstance(ell, (int, np.integer)):
            raise ValueError("beta ell must be an integer dyadic-band label")
        item = (int(ell), a)
        # repr is intentionally only a runtime uniqueness key; source labels
        # are not serialized or claimed recovered by this contract.
        key = repr(item)
        if key in keys:
            raise ValueError("beta_labels must be unique; duplicated signs are not partition entries")
        keys.append(key)
        normalized.append(item)
    return tuple(normalized)


@dataclass(frozen=True)
class KokunoSourcePartitionFamilyContract:
    """Validate a caller-supplied slow squared-partition family.

    The final array axis is the source ``beta=(ell,a)`` partition axis.  All
    preceding axes are arbitrary batch/sample axes.  ``D_r_eta`` and
    ``D_z_eta`` must use the same source-normalized derivative convention as
    the complete-curl lane.
    """

    partition_atol: float = 2.0e-12
    derivative_rtol: float = 2.0e-10

    def __post_init__(self) -> None:
        partition_atol = float(self.partition_atol)
        derivative_rtol = float(self.derivative_rtol)
        if not np.isfinite(partition_atol) or not (0.0 < partition_atol <= 1.0e-8):
            raise ValueError("partition_atol must lie in (0,1e-8]")
        if not np.isfinite(derivative_rtol) or not (0.0 < derivative_rtol <= 1.0e-6):
            raise ValueError("derivative_rtol must lie in (0,1e-6]")
        object.__setattr__(self, "partition_atol", partition_atol)
        object.__setattr__(self, "derivative_rtol", derivative_rtol)

    def validate_family(
        self,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        """Validate source partition normalization and its first derivatives."""
        eta = _finite_real(eta, "eta")
        D_r_eta = _finite_real(D_r_eta, "D_r_eta")
        D_z_eta = _finite_real(D_z_eta, "D_z_eta")
        if eta.ndim < 1 or eta.shape[-1] < 1:
            raise ValueError("eta must have a nonempty final beta partition axis")
        if D_r_eta.shape != eta.shape or D_z_eta.shape != eta.shape:
            raise ValueError("D_r_eta and D_z_eta must have exactly the eta shape")

        labels = _normalize_beta_labels(beta_labels, eta.shape[-1])
        squared_sum = np.sum(eta * eta, axis=-1)
        partition_error = squared_sum - 1.0
        if np.any(np.abs(partition_error) > self.partition_atol):
            raise ValueError("caller weights fail source squared partition sum_beta eta_beta^2=1")

        radial_closure = np.sum(eta * D_r_eta, axis=-1)
        axial_closure = np.sum(eta * D_z_eta, axis=-1)
        eta_norm = np.linalg.norm(eta, axis=-1)
        radial_scale = np.maximum(1.0, eta_norm * np.linalg.norm(D_r_eta, axis=-1))
        axial_scale = np.maximum(1.0, eta_norm * np.linalg.norm(D_z_eta, axis=-1))
        if np.any(np.abs(radial_closure) > self.derivative_rtol * radial_scale):
            raise ValueError("D_r_eta fails differentiated squared-partition closure")
        if np.any(np.abs(axial_closure) > self.derivative_rtol * axial_scale):
            raise ValueError("D_z_eta fails differentiated squared-partition closure")

        return {
            "eta": eta,
            "D_r_eta": D_r_eta,
            "D_z_eta": D_z_eta,
            "beta_labels": labels,
            "squared_sum": squared_sum,
            "partition_error": partition_error,
            "radial_closure": radial_closure,
            "axial_closure": axial_closure,
            "radial_scaled_error": np.abs(radial_closure) / radial_scale,
            "axial_scaled_error": np.abs(axial_closure) / axial_scale,
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
                "formula_scope": "slow beta squared-partition bookkeeping before support-localized complete curls",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "autonomous_derivations": dict(_AUTONOMOUS_DERIVATIONS),
            "parameters": {
                "partition_atol": self.partition_atol,
                "derivative_rtol": self.derivative_rtol,
                "origin": "bounded floating-point guards only; no numerical source cutoff is invented",
            },
            "caller_contract": {
                "required": [
                    "real eta_beta values indexed once by beta=(ell,a)",
                    "source-normalized D_r eta_beta and D_z eta_beta",
                    "unique beta labels; duplicated +/- rectangle signs are not partition entries",
                ],
                "partition_axis": "last array axis",
                "important": (
                    "this contract validates caller-supplied partition data only; concrete chi_ell, "
                    "chi_{ell,a}, S_*, actual labels, and positive-order background remain upstream"
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourcePartitionFamilyContract":
        if not isinstance(payload, dict):
            raise ValueError("source partition-family payload must be an object")
        required = {
            "schema", "source", "autonomous_derivations", "parameters",
            "caller_contract", "truth_boundary",
        }
        if set(payload) - {"sha256"} != required:
            raise ValueError("source partition-family payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source partition-family schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "partition_atol", "derivative_rtol", "origin"
        }:
            raise ValueError("source partition-family parameters changed")
        obj = cls(
            partition_atol=params["partition_atol"],
            derivative_rtol=params["derivative_rtol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"source partition-family {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("source partition-family sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourcePartitionFamilyContract":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
