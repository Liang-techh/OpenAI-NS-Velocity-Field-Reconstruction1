"""Carry Kokuno signed covariance amplitudes through localized complete curls.

The corrected 2026-09-09 reader uses auxiliary-rectangle sign ``sigma=+/-``
inside each slow label ``beta=(ell,a)``.  This sign is the covariance axis and
is distinct from the Fourier reality pair ``m=+/-1``.  This module composes the
existing signed reference inverse, per-band schedule, slow partition and
support-localized complete-curl contracts.  Caller-supplied source-normalized
background/mode data remain explicit and are not relabeled as recovered data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_source_band_covering import KokunoSourceBandCovering
from .kokuno_source_partition_family import KokunoSourcePartitionFamilyContract
from .kokuno_source_signed_covariance_pair import KokunoSourceSignedCovariancePair
from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-signed-complete-curl-family-v1"

_SOURCE_FORMULAS = {
    "rectangle_sign_axis": "gamma=(ell,a,sigma), sigma=+/-; sigma is distinct from Fourier m=+/-1",
    "reference_inverse": "H_ref=[[-A_c h_+,-A_c h_-],[-u_* h_+,+u_* h_-]]; H_ref y=T; y_sigma=a_sigma^2",
    "primary_scale": "W_0=sqrt(epsilon_ell) sum_sigma a_sigma b_sigma",
    "band_scales": "Q_ell=2^(-ell); epsilon_ell=Q_ell^h",
    "localized_potential": "A_{beta,sigma,m}=eta_beta*C_{beta,sigma,m}*exp(i*k_{ell,m}*Phi_{beta,sigma})",
    "physical_scaling": "u_phys=Q_ell^(-A) u_*; A=1/2+h",
}
_DERIVED_PRODUCT_RULES = {
    "amplitude": "g_sigma=sqrt(epsilon_ell)*a_sigma; Dg_sigma=sqrt(epsilon_ell)*Da_sigma",
    "coefficient": "D(g_sigma*C_sigma)=g_sigma*D(C_sigma)+Dg_sigma*C_sigma",
    "support": "D(eta_beta*g_sigma*C_sigma)=eta_beta*D(g_sigma*C_sigma)+(D eta_beta)*g_sigma*C_sigma",
}
_TRUTH_BOUNDARY = {
    "source_signed_reference_inverse_reused": True,
    "source_per_band_scale_schedule_reused": True,
    "signed_amplitude_gradient_complete_curl_terms_retained": True,
    "slow_support_gradient_complete_curl_terms_retained": True,
    "sign_resolved_real_fourier_pairs_executable": True,
    "q_scaled_sign_by_beta_total_physical_family_executable": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_partition_labels_instantiated": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "agent4_physical_covariance_audit_required": True,
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


def _signed_vector(value: Any, name: str, *, real: bool) -> np.ndarray:
    out = np.asarray(value, dtype=float if real else np.complex128)
    if out.ndim < 3 or out.shape[-2:] != (2, 3):
        raise ValueError(f"{name} must have shape sample_shape+(beta,2,3)")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _canonical_beta_indices(labels: Sequence[tuple[Any, Any]]) -> tuple[int, ...]:
    return tuple(sorted(range(len(labels)), key=lambda j: (int(labels[j][0]), repr(labels[j][1]))))


def _rotate_signed(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c = np.cos(theta)[..., None, None]
    s = np.sin(theta)[..., None, None]
    return np.stack((v[..., 0] * c - v[..., 1] * s,
                     v[..., 0] * s + v[..., 1] * c,
                     v[..., 2]), axis=-1)


@dataclass(frozen=True)
class KokunoSourceSignedCompleteCurlFamily:
    """Sigma-resolved physical complete-curl family from supplied source data."""

    h: float = 0.005
    partition_atol: float = 2.0e-12
    derivative_rtol: float = 2.0e-10
    reality_atol: float = 2.0e-12
    covariance_consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        h = float(self.h)
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        for name, value, upper in (
            ("partition_atol", self.partition_atol, 1.0e-8),
            ("derivative_rtol", self.derivative_rtol, 1.0e-6),
            ("reality_atol", self.reality_atol, 1.0e-8),
            ("covariance_consistency_atol", self.covariance_consistency_atol, 1.0e-8),
        ):
            val = float(value)
            if not np.isfinite(val) or not (0.0 < val <= upper):
                raise ValueError(f"{name} must lie in (0,{upper:g}]")
            object.__setattr__(self, name, val)
        object.__setattr__(self, "h", h)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def partition_contract(self) -> KokunoSourcePartitionFamilyContract:
        return KokunoSourcePartitionFamilyContract(self.partition_atol, self.derivative_rtol)

    @property
    def covariance_pair(self) -> KokunoSourceSignedCovariancePair:
        return KokunoSourceSignedCovariancePair(self.covariance_consistency_atol)

    @staticmethod
    def _beta_scalar(value: Any, name: str, shape: tuple[int, ...]) -> np.ndarray:
        value = _finite_real(value, name)
        try:
            return np.broadcast_to(value, shape)
        except ValueError:
            raise ValueError(f"{name} must broadcast to sample_shape+(beta,)") from None

    def physical_family(
        self,
        R: Any,
        theta: Any,
        phase: Any,
        n_phi: Any,
        t_plus_prototype: Any,
        D_r_C_plus_prototype: Any,
        D_z_C_plus_prototype: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
        *,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
        T_N: Any,
        T_K: Any,
        D_r_A_c: Any = 0.0,
        D_r_u_star: Any = 0.0,
        D_r_h_plus: Any = 0.0,
        D_r_h_minus: Any = 0.0,
        D_r_T_N: Any = 0.0,
        D_r_T_K: Any = 0.0,
        D_z_A_c: Any = 0.0,
        D_z_u_star: Any = 0.0,
        D_z_h_plus: Any = 0.0,
        D_z_h_minus: Any = 0.0,
        D_z_T_N: Any = 0.0,
        D_z_T_K: Any = 0.0,
        direction_gap_eta: float | None = None,
    ) -> dict[str, Any]:
        part = self.partition_contract.validate_family(eta, D_r_eta, D_z_eta, beta_labels)
        beta_shape = part["eta"].shape
        sample_shape = beta_shape[:-1]
        sign_shape = beta_shape + (2,)
        sign_vector_shape = beta_shape + (2, 3)

        phase = _finite_real(phase, "phase")
        n_phi = _signed_vector(n_phi, "n_phi", real=True)
        t_plus_prototype = _signed_vector(t_plus_prototype, "t_plus_prototype", real=False)
        D_r_C_plus_prototype = _signed_vector(D_r_C_plus_prototype, "D_r_C_plus_prototype", real=False)
        D_z_C_plus_prototype = _signed_vector(D_z_C_plus_prototype, "D_z_C_plus_prototype", real=False)
        if phase.shape != sign_shape:
            raise ValueError("phase must have shape eta.shape+(2,)")
        for name, value in (("n_phi", n_phi), ("t_plus_prototype", t_plus_prototype),
                            ("D_r_C_plus_prototype", D_r_C_plus_prototype),
                            ("D_z_C_plus_prototype", D_z_C_plus_prototype)):
            if value.shape != sign_vector_shape:
                raise ValueError(f"{name} must have shape eta.shape+(2,3)")

        R = _finite_real(R, "R")
        try:
            R_beta = np.broadcast_to(R, beta_shape)
        except ValueError:
            if R.shape == sample_shape:
                R_beta = np.broadcast_to(R[..., None], beta_shape)
            else:
                raise ValueError("R must broadcast to sample_shape+(beta,)") from None
        if np.any(R_beta <= 0.0):
            raise ValueError("source cylindrical complete curl requires R>0")
        theta = _finite_real(theta, "theta")
        try:
            theta = np.broadcast_to(theta, sample_shape)
        except ValueError:
            raise ValueError("theta must broadcast to sample axes only") from None

        scalar_names = (
            "A_c", "u_star", "h_plus", "h_minus", "T_N", "T_K",
            "D_r_A_c", "D_r_u_star", "D_r_h_plus", "D_r_h_minus", "D_r_T_N", "D_r_T_K",
            "D_z_A_c", "D_z_u_star", "D_z_h_plus", "D_z_h_minus", "D_z_T_N", "D_z_T_K",
        )
        scalar_raw = (
            A_c, u_star, h_plus, h_minus, T_N, T_K,
            D_r_A_c, D_r_u_star, D_r_h_plus, D_r_h_minus, D_r_T_N, D_r_T_K,
            D_z_A_c, D_z_u_star, D_z_h_plus, D_z_h_minus, D_z_T_N, D_z_T_K,
        )
        s = {name: self._beta_scalar(value, name, beta_shape)
             for name, value in zip(scalar_names, scalar_raw)}

        pair = self.covariance_pair
        base = pair.solve_reference_amplitudes(s["A_c"], s["u_star"], s["h_plus"], s["h_minus"],
                                               s["T_N"], s["T_K"], direction_gap_eta=direction_gap_eta)
        radial = pair.directional_derivative(
            s["A_c"], s["u_star"], s["h_plus"], s["h_minus"], s["T_N"], s["T_K"],
            dA_c=s["D_r_A_c"], du_star=s["D_r_u_star"], dh_plus=s["D_r_h_plus"],
            dh_minus=s["D_r_h_minus"], dT_N=s["D_r_T_N"], dT_K=s["D_r_T_K"])
        axial = pair.directional_derivative(
            s["A_c"], s["u_star"], s["h_plus"], s["h_minus"], s["T_N"], s["T_K"],
            dA_c=s["D_z_A_c"], du_star=s["D_z_u_star"], dh_plus=s["D_z_h_plus"],
            dh_minus=s["D_z_h_minus"], dT_N=s["D_z_T_N"], dT_K=s["D_z_T_K"])
        amplitudes = np.asarray(base["amplitudes"])
        D_r_amplitudes = np.asarray(radial["d_amplitudes"])
        D_z_amplitudes = np.asarray(axial["d_amplitudes"])

        ell_by_beta = tuple(int(label[0]) for label in part["beta_labels"])
        schedules = tuple(KokunoSourceBandCovering(ell, self.h) for ell in ell_by_beta)
        signed_cyl = []
        scaled_t = []
        scaled_DrC = []
        scaled_DzC = []
        for j, schedule in enumerate(schedules):
            sign_fields = []
            sign_t = []
            sign_DrC = []
            sign_DzC = []
            eps_root = float(np.sqrt(schedule.epsilon))
            for sigma_index in range(2):
                plus_localizer = KokunoSourceSupportLocalizedCurl(
                    schedule.epsilon, 1, self.partition_atol)
                n = n_phi[..., j, sigma_index, :]
                t_proto = t_plus_prototype[..., j, sigma_index, :]
                Dr_proto = D_r_C_plus_prototype[..., j, sigma_index, :]
                Dz_proto = D_z_C_plus_prototype[..., j, sigma_index, :]
                C_proto = plus_localizer.complete_curl.coefficient(n, t_proto)
                g = eps_root * amplitudes[..., j, sigma_index]
                D_r_g = eps_root * D_r_amplitudes[..., j, sigma_index]
                D_z_g = eps_root * D_z_amplitudes[..., j, sigma_index]
                t_scaled = g[..., None] * t_proto
                DrC_scaled = g[..., None] * Dr_proto + D_r_g[..., None] * C_proto
                DzC_scaled = g[..., None] * Dz_proto + D_z_g[..., None] * C_proto
                plus = plus_localizer.localized_mode(
                    R_beta[..., j], phase[..., j, sigma_index], n, t_scaled, DrC_scaled, DzC_scaled,
                    part["eta"][..., j], part["D_r_eta"][..., j], part["D_z_eta"][..., j])
                minus = KokunoSourceSupportLocalizedCurl(
                    schedule.epsilon, -1, self.partition_atol).localized_mode(
                    R_beta[..., j], phase[..., j, sigma_index], n, np.conjugate(t_scaled),
                    np.conjugate(DrC_scaled), np.conjugate(DzC_scaled), part["eta"][..., j],
                    part["D_r_eta"][..., j], part["D_z_eta"][..., j])
                pair_complex = plus["velocity"] + minus["velocity"]
                scale = max(1.0, float(np.max(np.abs(pair_complex.real), initial=0.0)))
                if float(np.max(np.abs(pair_complex.imag), initial=0.0)) > self.reality_atol * scale:
                    raise RuntimeError("signed m=+/-1 pair failed the declared reality guard")
                real_pair = pair_complex.real
                if not np.allclose(real_pair, 2.0 * plus["velocity"].real, rtol=0.0,
                                   atol=self.reality_atol * scale):
                    raise RuntimeError("signed real pair disagrees with 2*Re(m=+1)")
                sign_fields.append(schedule.Q ** (-self.A) * real_pair)
                sign_t.append(t_scaled)
                sign_DrC.append(DrC_scaled)
                sign_DzC.append(DzC_scaled)
            signed_cyl.append(np.stack(sign_fields, axis=-2))
            scaled_t.append(np.stack(sign_t, axis=-2))
            scaled_DrC.append(np.stack(sign_DrC, axis=-2))
            scaled_DzC.append(np.stack(sign_DzC, axis=-2))

        cyl_by_beta_sign = np.stack(signed_cyl, axis=-3)
        cart_by_beta_sign = _rotate_signed(cyl_by_beta_sign, theta)
        cyl_by_beta = np.sum(cyl_by_beta_sign, axis=-2)
        cart_by_beta = np.sum(cart_by_beta_sign, axis=-2)
        canonical = _canonical_beta_indices(part["beta_labels"])
        cyl_total = np.sum(np.take(cyl_by_beta, canonical, axis=-2), axis=-2)
        cart_total = np.sum(np.take(cart_by_beta, canonical, axis=-2), axis=-2)
        if not np.all(np.isfinite(cart_total)):
            raise RuntimeError("signed complete-curl family produced nonfinite physical velocity")

        return {
            "beta_labels": part["beta_labels"],
            "sigma_labels": pair.sign_labels(),
            "ell_by_beta": ell_by_beta,
            "canonical_beta_indices": canonical,
            "Q_by_beta": np.asarray([item.Q for item in schedules], dtype=float),
            "epsilon_by_beta": np.asarray([item.epsilon for item in schedules], dtype=float),
            "A": self.A,
            "partition_error": part["partition_error"],
            "partition_radial_closure": part["radial_closure"],
            "partition_axial_closure": part["axial_closure"],
            "reference_squared_amplitudes": np.asarray(base["squared_amplitudes"]),
            "reference_amplitudes": amplitudes,
            "D_r_reference_amplitudes": D_r_amplitudes,
            "D_z_reference_amplitudes": D_z_amplitudes,
            "scaled_t_plus_prototype": np.stack(scaled_t, axis=-3),
            "scaled_D_r_C_plus_prototype": np.stack(scaled_DrC, axis=-3),
            "scaled_D_z_C_plus_prototype": np.stack(scaled_DzC, axis=-3),
            "velocity_physical_cylindrical_by_beta_sign": cyl_by_beta_sign,
            "velocity_physical_cartesian_by_beta_sign": cart_by_beta_sign,
            "velocity_physical_cylindrical_by_beta": cyl_by_beta,
            "velocity_physical_cartesian_by_beta": cart_by_beta,
            "velocity_physical_cylindrical_total": cyl_total,
            "velocity_physical_cartesian_total": cart_total,
            "reference_covariance_rank_two": bool(base["reference_covariance_rank_two"]),
            "physical_complete_curl_covariance_rank_two_assessed": False,
            "genuinely_independent_second_covariance_column_ready": False,
        }

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {"repository": SOURCE_REPOSITORY, "commit": SOURCE_COMMIT, "path": SOURCE_PATH,
                       "corrected_release": CORRECTED_RELEASE,
                       "corrected_release_date": CORRECTED_RELEASE_DATE,
                       "formula_scope": "signed covariance amplitudes through localized complete curls",
                       "source_formulas": dict(_SOURCE_FORMULAS)},
            "derived_product_rules": dict(_DERIVED_PRODUCT_RULES),
            "parameters": {"h": self.h, "partition_atol": self.partition_atol,
                           "derivative_rtol": self.derivative_rtol, "reality_atol": self.reality_atol,
                           "covariance_consistency_atol": self.covariance_consistency_atol,
                           "origin": "source formulas plus bounded floating-point consistency guards"},
            "coefficient_units": {"a_sigma": "source reference multiplier of one signed prototype",
                                  "sqrt_epsilon_a_sigma": "normalized primary-wave multiplier before Q^(-A)",
                                  "physical_velocity": "Q_ell^(-A) times normalized complete-curl velocity",
                                  "important": "source h_sigma/background/rectangle numerics remain caller supplied"},
            "caller_contract": {"axis": "R must stay strictly positive; no unsourced axis regularization",
                                "signs": "sigma=+/- is covariance/rectangle; m=+/-1 is Fourier reality",
                                "partition": "one eta entry per beta=(ell,a), never duplicated over sigma"},
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceSignedCompleteCurlFamily":
        if not isinstance(payload, dict):
            raise ValueError("signed complete-curl family payload must be an object")
        required = {"schema", "source", "derived_product_rules", "parameters", "coefficient_units",
                    "caller_contract", "truth_boundary"}
        if set(payload) - {"sha256"} != required or payload.get("schema") != SCHEMA:
            raise ValueError("signed complete-curl family payload keys or schema do not match")
        params = payload["parameters"]
        obj = cls(params["h"], params["partition_atol"], params["derivative_rtol"],
                  params["reality_atol"], params["covariance_consistency_atol"])
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"signed complete-curl family {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("signed complete-curl family sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceSignedCompleteCurlFamily":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
