"""Low-dimensional bounded coefficient coordinates for Kokuno-derived label families.

This module sits downstream of :mod:`kokuno_source_localized_real_pair_family`.
It does not reconstruct another pulse, cutoff, background, or mean correction.
Instead it gives an already-Q-scaled physical ``by_beta`` velocity family an
explicit two-coordinate, dimensionless coefficient normalization that Agent 3
can differentiate and budget without introducing one free parameter per label.

For physical label velocities ``w_beta`` it uses the repository-autonomous map

    m_beta(delta_0, delta_1) = 1 + delta_0 + delta_1 g_beta,

where ``g_beta`` is a deterministic band contrast in ``[-1,1]`` obtained only
from the source slow-band label ``ell``.  The admissible set is

    |delta_0| + |delta_1| <= B <= 1/2,

so every label multiplier remains positive.  Consequently the exact tangent
columns at the unmodified family are

    d_{delta_0} u = sum_beta w_beta,
    d_{delta_1} u = sum_beta g_beta w_beta.

These coordinates are autonomous engineering choices, not Kokuno's unpublished
primary amplitudes or the source signed differential inverse.  In particular,
``source_coefficient_unit_mapping_available`` remains false.  The useful claim
is narrower: the repository now has explicit, bounded, dimensionless units for
any supplied physical Agent-2 family, with a fixed parameter count independent
of the number of active labels.
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
SCHEMA = "kokuno-bounded-family-coefficient-coordinates-v1"

_SOURCE_SCOPE = {
    "slow_label": "beta=(ell,a)",
    "physical_family_input": "consume already-localized real-pair label velocities after per-label Q^(-A) scaling",
    "source_primary_amplitudes": "not recovered by this module",
    "source_signed_inverse": "not reimplemented by this module",
}

_AUTONOMOUS_PARAMETERIZATION = {
    "coordinates": "m_beta=1+delta_common+delta_band*g_beta",
    "band_contrast": "g_beta is affine in ell over the active ell range and lies in [-1,1]",
    "bound": "|delta_common|+|delta_band|<=max_l1_update<=1/2",
    "units": "both deltas are dimensionless fractional multipliers of already-Q-scaled physical by-beta velocity",
    "anti_budget_laundering": "parameter count is fixed at two and does not grow with beta count",
    "phase_policy": "no phase-copy coordinate is introduced; phase copies are not used as a substitute for an independent covariance direction",
}

_TRUTH_BOUNDARY = {
    "repository_coefficient_unit_mapping_available": True,
    "low_dimensional_bounded_parameterization_available": True,
    "parameter_count_independent_of_label_count": True,
    "source_coefficient_unit_mapping_available": False,
    "source_primary_amplitudes_recovered": False,
    "source_signed_differential_inverse_instantiated": False,
    "actual_positive_order_background_bound": False,
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


def _normalize_labels(beta_labels: Sequence[Any], n_beta: int) -> tuple[tuple[Any, ...], ...]:
    if len(beta_labels) != n_beta:
        raise ValueError("beta_labels length must equal the by-beta velocity axis")
    normalized: list[tuple[Any, ...]] = []
    seen: set[str] = set()
    for label in beta_labels:
        if not isinstance(label, (tuple, list)) or len(label) != 2:
            raise ValueError("every beta label must be the slow label pair (ell,a)")
        ell = label[0]
        if isinstance(ell, bool) or not isinstance(ell, (int, np.integer)):
            raise ValueError("beta label ell must be an integer")
        item = (int(ell), label[1])
        marker = repr(item)
        if marker in seen:
            raise ValueError("beta labels must be unique; per-label duplication can launder coefficient budgets")
        seen.add(marker)
        normalized.append(item)
    return tuple(normalized)


def _rotate_cylindrical(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c = np.cos(theta)
    s = np.sin(theta)
    return np.stack(
        (
            v[..., 0] * c - v[..., 1] * s,
            v[..., 0] * s + v[..., 1] * c,
            v[..., 2],
        ),
        axis=-1,
    )


@dataclass(frozen=True)
class KokunoBoundedFamilyCoefficientCoordinates:
    """Two bounded dimensionless coordinates on a physical by-beta family."""

    max_l1_update: float = 0.125
    consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        bound = float(self.max_l1_update)
        atol = float(self.consistency_atol)
        if not np.isfinite(bound) or not (0.0 < bound <= 0.5):
            raise ValueError("max_l1_update must lie in (0,0.5]")
        if not np.isfinite(atol) or not (0.0 < atol <= 1.0e-8):
            raise ValueError("consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "max_l1_update", bound)
        object.__setattr__(self, "consistency_atol", atol)

    @staticmethod
    def parameter_names() -> tuple[str, str]:
        return ("delta_common", "delta_band")

    @staticmethod
    def parameter_units() -> tuple[str, str]:
        unit = "dimensionless_fractional_multiplier_of_Q_scaled_physical_velocity"
        return (unit, unit)

    def contrast_weights(self, beta_labels: Sequence[Any], n_beta: int) -> dict[str, Any]:
        labels = _normalize_labels(beta_labels, n_beta)
        ell = np.asarray([label[0] for label in labels], dtype=float)
        ell_min = float(np.min(ell))
        ell_max = float(np.max(ell))
        if ell_max == ell_min:
            weights = np.zeros(n_beta, dtype=float)
            active = False
            center = ell_min
            half_span = 0.0
        else:
            center = 0.5 * (ell_min + ell_max)
            half_span = 0.5 * (ell_max - ell_min)
            weights = (ell - center) / half_span
            active = True
        if np.max(np.abs(weights), initial=0.0) > 1.0 + 4.0 * np.finfo(float).eps:
            raise RuntimeError("internal band contrast escaped [-1,1]")
        return {
            "beta_labels": labels,
            "ell": ell,
            "ell_center": center,
            "ell_half_span": half_span,
            "band_contrast": weights,
            "band_contrast_active": active,
        }

    def _validate_deltas(self, delta_common: float, delta_band: float, contrast_active: bool) -> tuple[float, float]:
        d0 = float(delta_common)
        d1 = float(delta_band)
        if not np.isfinite(d0) or not np.isfinite(d1):
            raise ValueError("coefficient deltas must be finite")
        if abs(d0) + abs(d1) > self.max_l1_update + 8.0 * np.finfo(float).eps:
            raise ValueError("coefficient update exceeds the declared aggregate L1 bound")
        if not contrast_active and d1 != 0.0:
            raise ValueError("delta_band is unavailable when the active family contains only one ell band")
        return d0, d1

    def evaluate(
        self,
        velocity_physical_cylindrical_by_beta: Any,
        beta_labels: Sequence[Any],
        *,
        delta_common: float = 0.0,
        delta_band: float = 0.0,
        theta: Any | None = None,
    ) -> dict[str, Any]:
        """Apply bounded coordinates and return exact physical tangent columns.

        The input must already include localization, real m=+/-1 reconstruction,
        and per-label Q^(-A) scaling.  Its shape is ``sample_shape+(beta,3)``.
        """
        velocity = _finite_real(
            velocity_physical_cylindrical_by_beta,
            "velocity_physical_cylindrical_by_beta",
        )
        if velocity.ndim < 2 or velocity.shape[-1] != 3:
            raise ValueError("by-beta velocity must have shape sample_shape+(beta,3)")
        n_beta = int(velocity.shape[-2])
        if n_beta < 1:
            raise ValueError("by-beta velocity must contain at least one beta label")
        labels = self.contrast_weights(beta_labels, n_beta)
        d0, d1 = self._validate_deltas(delta_common, delta_band, labels["band_contrast_active"])
        g = labels["band_contrast"]

        multipliers = 1.0 + d0 + d1 * g
        if np.any(multipliers <= 0.0):
            raise RuntimeError("bounded coefficient contract produced a nonpositive label multiplier")
        modulated_by_beta = velocity * multipliers.reshape((1,) * (velocity.ndim - 2) + (n_beta, 1))
        base_total = np.sum(velocity, axis=-2)
        modulated_total = np.sum(modulated_by_beta, axis=-2)
        common_tangent = base_total
        band_tangent = np.sum(
            velocity * g.reshape((1,) * (velocity.ndim - 2) + (n_beta, 1)), axis=-2
        )
        tangents = np.stack((common_tangent, band_tangent), axis=-2)

        out: dict[str, Any] = {
            **labels,
            "parameter_names": self.parameter_names(),
            "parameter_units": self.parameter_units(),
            "parameter_count": 2,
            "max_l1_update": self.max_l1_update,
            "delta_common": d0,
            "delta_band": d1,
            "aggregate_l1_update": abs(d0) + abs(d1),
            "label_multipliers": multipliers,
            "minimum_label_multiplier": float(np.min(multipliers)),
            "velocity_physical_cylindrical_by_beta_base": velocity,
            "velocity_physical_cylindrical_total_base": base_total,
            "velocity_physical_cylindrical_by_beta_modulated": modulated_by_beta,
            "velocity_physical_cylindrical_total_modulated": modulated_total,
            "velocity_parameter_tangent_cylindrical": tangents,
            "repository_coefficient_unit_mapping_available": True,
            "source_coefficient_unit_mapping_available": False,
        }

        if theta is not None:
            sample_shape = base_total.shape[:-1]
            theta_array = _finite_real(theta, "theta")
            try:
                theta_sample = np.broadcast_to(theta_array, sample_shape)
            except ValueError:
                raise ValueError("theta must broadcast to sample_shape") from None
            out["velocity_physical_cartesian_total_base"] = _rotate_cylindrical(base_total, theta_sample)
            out["velocity_physical_cartesian_total_modulated"] = _rotate_cylindrical(modulated_total, theta_sample)
            tangent_theta = np.broadcast_to(theta_sample[..., None], sample_shape + (2,))
            out["velocity_parameter_tangent_cartesian"] = _rotate_cylindrical(tangents, tangent_theta)
        return out

    def consume_physical_family(
        self,
        physical_family_result: dict[str, Any],
        *,
        delta_common: float = 0.0,
        delta_band: float = 0.0,
        theta: Any | None = None,
    ) -> dict[str, Any]:
        """Consume the exact result dictionary emitted by Agent-2 #380."""
        required = (
            "beta_labels",
            "velocity_physical_cylindrical_by_beta",
            "velocity_physical_cylindrical_total",
        )
        missing = [key for key in required if key not in physical_family_result]
        if missing:
            raise ValueError(f"physical family result is missing required keys: {missing}")
        velocity = _finite_real(
            physical_family_result["velocity_physical_cylindrical_by_beta"],
            "velocity_physical_cylindrical_by_beta",
        )
        declared_total = _finite_real(
            physical_family_result["velocity_physical_cylindrical_total"],
            "velocity_physical_cylindrical_total",
        )
        computed_total = np.sum(velocity, axis=-2)
        scale = max(1.0, float(np.max(np.abs(computed_total), initial=0.0)))
        if declared_total.shape != computed_total.shape or not np.allclose(
            declared_total, computed_total, rtol=0.0, atol=self.consistency_atol * scale
        ):
            raise ValueError("physical family declared total does not equal the by-beta sum")
        out = self.evaluate(
            velocity,
            physical_family_result["beta_labels"],
            delta_common=delta_common,
            delta_band=delta_band,
            theta=theta,
        )
        out["consumed_agent2_physical_family_contract"] = True
        return out

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "scope": dict(_SOURCE_SCOPE),
            },
            "autonomous_parameterization": dict(_AUTONOMOUS_PARAMETERIZATION),
            "parameters": {
                "max_l1_update": self.max_l1_update,
                "consistency_atol": self.consistency_atol,
                "parameter_names": list(self.parameter_names()),
                "parameter_units": list(self.parameter_units()),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        return output

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoBoundedFamilyCoefficientCoordinates":
        if payload.get("schema") != SCHEMA:
            raise ValueError("unexpected coefficient-coordinate schema")
        source = payload.get("source")
        if not isinstance(source, dict):
            raise ValueError("source metadata is missing")
        expected_source = {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
            "scope": dict(_SOURCE_SCOPE),
        }
        if source != expected_source:
            raise ValueError("source provenance does not match the pinned corrected reconstruction")
        if payload.get("autonomous_parameterization") != _AUTONOMOUS_PARAMETERIZATION:
            raise ValueError("autonomous_parameterization metadata is inconsistent")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth_boundary metadata is inconsistent")
        parameters = payload.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError("parameters metadata is missing")
        if tuple(parameters.get("parameter_names", ())) != cls.parameter_names():
            raise ValueError("parameter_names metadata is inconsistent")
        if tuple(parameters.get("parameter_units", ())) != cls.parameter_units():
            raise ValueError("parameter_units metadata is inconsistent")
        return cls(
            max_l1_update=float(parameters["max_l1_update"]),
            consistency_atol=float(parameters["consistency_atol"]),
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoBoundedFamilyCoefficientCoordinates":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        expected_sha = payload.pop("sha256", None)
        obj = cls.from_payload(payload)
        if expected_sha is not None and expected_sha != obj.sha256:
            raise ValueError("serialized coefficient-coordinate SHA256 does not match payload")
        return obj
