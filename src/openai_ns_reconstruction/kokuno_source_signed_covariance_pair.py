"""Source-facing signed auxiliary-rectangle covariance pair for Kokuno waves.

The corrected 2026-09-09 reconstruction has two *rectangle signs* ``sigma=+/-``
inside each slow box.  This sign is distinct from the Fourier harmonic sign
``m=+/-1`` used to make one real complete curl.  The two sigma fields have
disjoint auxiliary-torus rectangles and provide the two covariance columns
needed by the signed stress construction.

At the reference-column level the source gives

    H_ref = [[-A_c h_+, -A_c h_-],
             [-u_* h_+, +u_* h_-]],

with ``det(H_ref) = -2 A_c u_* h_+ h_-``.  For a target
``T=(T_N,T_K)``, the positive squared amplitudes ``y=(a_+^2,a_-^2)`` solve
``H_ref y = T``.  Equivalently,

    h_+ y_+ = 1/2 (-T_N/A_c - T_K/u_*),
    h_- y_- = 1/2 (-T_N/A_c + T_K/u_*).

This module makes that source seam executable and differentiable.  It does not
claim that the repository has recovered the actual source background values,
actual pulse integrals ``h_sigma``, actual rectangle labels, or the final
complete-curl physical velocity.  Those remain upstream/downstream bindings.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-signed-covariance-pair-v1"

_SOURCE_FORMULAS = {
    "covariance": "C(w)=<<w_r w_tan>_theta>_Y, component order (rtheta,rz)",
    "rectangle_sign_fields": "b_sigma=chi_g(xi_g) psi(v) t_sigma^h cos(k Phi_sigma), sigma=+/- on disjoint auxiliary rectangles",
    "reference_columns": "H_ref=[[-A_c h_+,-A_c h_-],[-u_* h_+,+u_* h_-]]",
    "determinant": "det(H_ref)=-2 A_c u_* h_+ h_-",
    "positive_inverse": "h_+ y_+=(-T_N/A_c-T_K/u_*)/2; h_- y_-=(-T_N/A_c+T_K/u_*)/2",
    "primary_wave": "W_0=sqrt(epsilon) sum_{sigma=+/-} a_sigma b_sigma, a_sigma=sqrt(y_sigma)",
    "primary_covariance": "C(W_0)=epsilon H (a_+^2,a_-^2)^T=epsilon T_{0,*}",
}

_TRUTH_BOUNDARY = {
    "source_rectangle_sign_axis_identified": True,
    "rectangle_sign_sigma_distinct_from_fourier_harmonic_m": True,
    "source_reference_two_column_covariance_map_executable": True,
    "source_reference_positive_amplitude_inverse_executable": True,
    "source_reference_amplitude_directional_derivative_executable": True,
    "source_primary_tangent_wave_assembly_executable": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_rectangle_labels_instantiated": False,
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


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _broadcast_named(**values: Any) -> dict[str, np.ndarray]:
    arrays = {name: _finite(value, name) for name, value in values.items()}
    try:
        shape = np.broadcast_shapes(*(arr.shape for arr in arrays.values()))
    except ValueError as exc:
        raise ValueError("source covariance inputs must be broadcast-compatible") from exc
    return {name: np.broadcast_to(arr, shape) for name, arr in arrays.items()}


@dataclass(frozen=True)
class KokunoSourceSignedCovariancePair:
    """Reference two-sign covariance inverse and amplitude differential contract."""

    consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        atol = float(self.consistency_atol)
        if not np.isfinite(atol) or not (0.0 < atol <= 1.0e-8):
            raise ValueError("consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "consistency_atol", atol)

    @staticmethod
    def sign_labels() -> tuple[str, str]:
        return ("sigma_plus", "sigma_minus")

    @staticmethod
    def source_input_units() -> dict[str, str]:
        return {
            "A_c": "source_reference_covariance_normal_coordinate",
            "u_star": "source_reference_covariance_tangential_coordinate",
            "h_plus": "source_positive_pulse_covariance_mass",
            "h_minus": "source_positive_pulse_covariance_mass",
            "T_N": "source_target_normal_stress_coordinate",
            "T_K": "source_target_tangential_stress_coordinate",
            "y_sigma": "squared_multiplier_of_source_sigma_prototype",
            "a_sigma": "multiplier_of_source_sigma_prototype",
        }

    def _validated_base(
        self,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
        T_N: Any,
        T_K: Any,
    ) -> dict[str, np.ndarray]:
        data = _broadcast_named(
            A_c=A_c,
            u_star=u_star,
            h_plus=h_plus,
            h_minus=h_minus,
            T_N=T_N,
            T_K=T_K,
        )
        for name in ("A_c", "u_star", "h_plus", "h_minus"):
            if np.any(data[name] <= 0.0):
                raise ValueError(f"{name} must be strictly positive on the source open shell")
        return data

    def reference_column_matrix(
        self,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
    ) -> np.ndarray:
        data = _broadcast_named(A_c=A_c, u_star=u_star, h_plus=h_plus, h_minus=h_minus)
        for name in ("A_c", "u_star", "h_plus", "h_minus"):
            if np.any(data[name] <= 0.0):
                raise ValueError(f"{name} must be strictly positive on the source open shell")
        return np.stack(
            (
                np.stack((-data["A_c"] * data["h_plus"], -data["A_c"] * data["h_minus"]), axis=-1),
                np.stack((-data["u_star"] * data["h_plus"], +data["u_star"] * data["h_minus"]), axis=-1),
            ),
            axis=-2,
        )

    def solve_reference_amplitudes(
        self,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
        T_N: Any,
        T_K: Any,
        *,
        direction_gap_eta: float | None = None,
    ) -> dict[str, Any]:
        """Solve the source reference 2x2 inverse and require positive amplitudes.

        ``direction_gap_eta`` is optional.  When supplied it verifies the source
        comparison inequality ``|q| <= (1-eta) p`` for
        ``p=-T_N/A_c`` and ``q=T_K/u_star``.  The value itself is caller-bound;
        this module does not invent the source compact-set margin.
        """
        data = self._validated_base(A_c, u_star, h_plus, h_minus, T_N, T_K)
        p = -data["T_N"] / data["A_c"]
        q = data["T_K"] / data["u_star"]
        if direction_gap_eta is not None:
            eta = float(direction_gap_eta)
            if not np.isfinite(eta) or not (0.0 < eta < 1.0):
                raise ValueError("direction_gap_eta must lie strictly in (0,1)")
            scale = np.maximum(1.0, np.abs(p))
            if np.any(np.abs(q) > (1.0 - eta) * p + self.consistency_atol * scale):
                raise ValueError("target violates the declared strict source direction gap")
        else:
            eta = None

        y_plus = 0.5 * (p - q) / data["h_plus"]
        y_minus = 0.5 * (p + q) / data["h_minus"]
        y = np.stack((y_plus, y_minus), axis=-1)
        if np.any(y <= 0.0):
            raise ValueError("target is outside the positive two-sign source covariance cone")
        amplitudes = np.sqrt(y)

        H = self.reference_column_matrix(
            data["A_c"], data["u_star"], data["h_plus"], data["h_minus"]
        )
        target = np.stack((data["T_N"], data["T_K"]), axis=-1)
        reconstructed = np.einsum("...ij,...j->...i", H, y)
        scale = np.maximum(1.0, np.max(np.abs(target), axis=-1, keepdims=True))
        if not np.all(np.abs(reconstructed - target) <= self.consistency_atol * scale):
            raise RuntimeError("analytic source reference inverse failed reconstruction check")

        determinant = np.linalg.det(H)
        determinant_formula = -2.0 * data["A_c"] * data["u_star"] * data["h_plus"] * data["h_minus"]
        det_scale = np.maximum(1.0, np.abs(determinant_formula))
        if not np.all(np.abs(determinant - determinant_formula) <= self.consistency_atol * det_scale):
            raise RuntimeError("reference covariance determinant disagrees with source identity")
        singular_values = np.linalg.svd(H, compute_uv=False)
        ratio = singular_values[..., 1] / singular_values[..., 0]

        return {
            **data,
            "sign_labels": self.sign_labels(),
            "p": p,
            "q": q,
            "direction_gap_eta": eta,
            "reference_column_matrix": H,
            "reference_column_determinant": determinant,
            "squared_amplitudes": y,
            "amplitudes": amplitudes,
            "target": target,
            "reconstructed_target": reconstructed,
            "reference_column_singular_values": singular_values,
            "reference_column_min_singular_ratio": float(np.min(ratio)),
            "reference_covariance_rank_two": bool(np.all(singular_values[..., 1] > 0.0)),
            "positive_reference_coefficients": True,
            "sigma_axis_is_not_fourier_harmonic_axis": True,
            "actual_source_mode_family_bound": False,
            "genuinely_independent_second_covariance_column_ready": False,
        }

    def directional_derivative(
        self,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
        T_N: Any,
        T_K: Any,
        *,
        dA_c: Any = 0.0,
        du_star: Any = 0.0,
        dh_plus: Any = 0.0,
        dh_minus: Any = 0.0,
        dT_N: Any = 0.0,
        dT_K: Any = 0.0,
    ) -> dict[str, Any]:
        """Differentiate ``H y=T`` exactly along one slow-coordinate direction."""
        base = self.solve_reference_amplitudes(A_c, u_star, h_plus, h_minus, T_N, T_K)
        tangent = _broadcast_named(
            dA_c=dA_c,
            du_star=du_star,
            dh_plus=dh_plus,
            dh_minus=dh_minus,
            dT_N=dT_N,
            dT_K=dT_K,
            anchor=np.asarray(base["A_c"]),
        )
        A = np.asarray(base["A_c"])
        u = np.asarray(base["u_star"])
        hp = np.asarray(base["h_plus"])
        hm = np.asarray(base["h_minus"])
        dA = tangent["dA_c"]
        du = tangent["du_star"]
        dhp = tangent["dh_plus"]
        dhm = tangent["dh_minus"]
        dH = np.stack(
            (
                np.stack((-(dA * hp + A * dhp), -(dA * hm + A * dhm)), axis=-1),
                np.stack((-(du * hp + u * dhp), +(du * hm + u * dhm)), axis=-1),
            ),
            axis=-2,
        )
        dT = np.stack((tangent["dT_N"], tangent["dT_K"]), axis=-1)
        y = np.asarray(base["squared_amplitudes"])
        rhs = dT - np.einsum("...ij,...j->...i", dH, y)
        dy = np.linalg.solve(np.asarray(base["reference_column_matrix"]), rhs)
        a = np.asarray(base["amplitudes"])
        da = 0.5 * dy / a
        identity_error = np.einsum("...ij,...j->...i", dH, y) + np.einsum(
            "...ij,...j->...i", np.asarray(base["reference_column_matrix"]), dy
        ) - dT
        scale = np.maximum(1.0, np.max(np.abs(dT), axis=-1, keepdims=True))
        if not np.all(np.abs(identity_error) <= self.consistency_atol * scale):
            raise RuntimeError("differentiated source covariance identity failed")
        return {
            "base": base,
            "d_reference_column_matrix": dH,
            "d_target": dT,
            "d_squared_amplitudes": dy,
            "d_amplitudes": da,
            "differentiated_identity_error": identity_error,
        }

    def assemble_primary_tangent_wave(
        self,
        sign_prototypes: Any,
        amplitudes: Any,
        epsilon: Any,
    ) -> np.ndarray:
        """Evaluate only the source primary ``W_0`` formula, before complete curl.

        ``sign_prototypes`` has shape ``sample_shape+(2,3)`` ordered
        ``(sigma_plus,sigma_minus)``.  This helper is intentionally not called a
        divergence-free complete curl; support-gradient curl remainders belong
        to the existing Agent-2 complete-curl path.
        """
        b = _finite(sign_prototypes, "sign_prototypes")
        a = _finite(amplitudes, "amplitudes")
        eps = _finite(epsilon, "epsilon")
        if b.ndim < 2 or b.shape[-2:] != (2, 3):
            raise ValueError("sign_prototypes must have shape sample_shape+(2,3)")
        sample_shape = b.shape[:-2]
        try:
            a = np.broadcast_to(a, sample_shape + (2,))
            eps = np.broadcast_to(eps, sample_shape)
        except ValueError as exc:
            raise ValueError("amplitudes/epsilon must broadcast to sign_prototypes sample axes") from exc
        if np.any(a <= 0.0) or np.any(eps <= 0.0):
            raise ValueError("amplitudes and epsilon must be strictly positive")
        return np.sqrt(eps)[..., None] * np.sum(a[..., :, None] * b, axis=-2)

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {"consistency_atol": self.consistency_atol},
            "source_input_units": self.source_input_units(),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceSignedCovariancePair":
        if payload.get("schema") != SCHEMA:
            raise ValueError("unexpected signed-covariance-pair schema")
        expected_source = {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
            "formulas": dict(_SOURCE_FORMULAS),
        }
        if payload.get("source") != expected_source:
            raise ValueError("source provenance/formulas do not match the pinned corrected reconstruction")
        if payload.get("source_input_units") != cls.source_input_units():
            raise ValueError("source_input_units metadata is inconsistent")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth_boundary metadata is inconsistent")
        parameters = payload.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError("parameters are missing")
        obj = cls(consistency_atol=float(parameters["consistency_atol"]))
        declared = payload.get("sha256")
        if declared is not None and declared != obj.sha256:
            raise ValueError("signed covariance pair SHA does not match payload")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceSignedCovariancePair":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("signed covariance pair JSON must contain an object")
        return cls.from_payload(payload)
