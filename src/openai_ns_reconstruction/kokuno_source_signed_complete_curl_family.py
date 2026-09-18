"""Carry Kokuno signed covariance amplitudes through localized complete curls.

The corrected 2026-09-09 reader uses auxiliary-rectangle sign ``sigma=+/-``
inside each slow label ``beta=(ell,a)``. This is the covariance axis, distinct
from the Fourier reality pair ``m=+/-1``. The actual signed rectangles,
positive-order background, h_sigma data and mode labels remain caller supplied.
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
    "sign_axis": "gamma=(ell,a,sigma), sigma=+/-; sigma is distinct from Fourier m=+/-1",
    "reference_inverse": "H_ref y=T, y_sigma=a_sigma^2",
    "physical_covariance_units": "epsilon_ell H_ref y=Delta C, so reference T=Delta C/epsilon_ell",
    "primary_scale": "W_0=sqrt(epsilon_ell) sum_sigma a_sigma b_sigma",
    "band_scales": "Q_ell=2^(-ell); epsilon_ell=Q_ell^h",
    "physical_scaling": "u_phys=Q_ell^(-A) u_*, A=1/2+h",
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


def _canonical_json(x: dict[str, Any]) -> str:
    return json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _real(x: Any, name: str) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain only finite real values")
    return a


def _signed_vector(x: Any, name: str, real: bool) -> np.ndarray:
    a = np.asarray(x, dtype=float if real else np.complex128)
    if a.ndim < 3 or a.shape[-2:] != (2, 3):
        raise ValueError(f"{name} must have shape sample_shape+(beta,2,3)")
    if not np.all(np.isfinite(a.real)) or not np.all(np.isfinite(a.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return a


def _canonical_indices(labels: Sequence[tuple[Any, Any]]) -> tuple[int, ...]:
    return tuple(sorted(range(len(labels)), key=lambda j: (int(labels[j][0]), repr(labels[j][1]))))


def _rotate(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c = np.cos(theta)[..., None, None]
    s = np.sin(theta)[..., None, None]
    return np.stack((v[..., 0] * c - v[..., 1] * s,
                     v[..., 0] * s + v[..., 1] * c,
                     v[..., 2]), axis=-1)


@dataclass(frozen=True)
class KokunoSourceSignedCompleteCurlFamily:
    h: float = 0.005
    partition_atol: float = 2e-12
    derivative_rtol: float = 2e-10
    reality_atol: float = 2e-12
    covariance_consistency_atol: float = 2e-12

    def __post_init__(self) -> None:
        h = float(self.h)
        if not np.isfinite(h) or not (0 < h < 1e-2):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        for name, value, upper in (
            ("partition_atol", self.partition_atol, 1e-8),
            ("derivative_rtol", self.derivative_rtol, 1e-6),
            ("reality_atol", self.reality_atol, 1e-8),
            ("covariance_consistency_atol", self.covariance_consistency_atol, 1e-8),
        ):
            value = float(value)
            if not np.isfinite(value) or not (0 < value <= upper):
                raise ValueError(f"{name} must lie in (0,{upper:g}]")
            object.__setattr__(self, name, value)
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
    def _beta_scalar(x: Any, name: str, shape: tuple[int, ...]) -> np.ndarray:
        try:
            return np.broadcast_to(_real(x, name), shape)
        except ValueError:
            raise ValueError(f"{name} must broadcast to sample_shape+(beta,)") from None

    @staticmethod
    def _scalar_directional_amplitudes(
        pair: KokunoSourceSignedCovariancePair,
        s: dict[str, np.ndarray],
        prefix: str,
        shape: tuple[int, ...],
    ) -> np.ndarray:
        """Use the audited #472 scalar derivative pointwise over batch axes."""
        out = np.empty(shape + (2,), dtype=float)
        for index in np.ndindex(shape):
            r = pair.directional_derivative(
                float(s["A_c"][index]), float(s["u_star"][index]),
                float(s["h_plus"][index]), float(s["h_minus"][index]),
                float(s["T_N"][index]), float(s["T_K"][index]),
                dA_c=float(s[f"{prefix}_A_c"][index]),
                du_star=float(s[f"{prefix}_u_star"][index]),
                dh_plus=float(s[f"{prefix}_h_plus"][index]),
                dh_minus=float(s[f"{prefix}_h_minus"][index]),
                dT_N=float(s[f"{prefix}_T_N"][index]),
                dT_K=float(s[f"{prefix}_T_K"][index]),
            )
            out[index] = np.asarray(r["d_amplitudes"], dtype=float)
        return out

    def physical_family(
        self, R: Any, theta: Any, phase: Any, n_phi: Any,
        t_plus_prototype: Any, D_r_C_plus_prototype: Any, D_z_C_plus_prototype: Any,
        eta: Any, D_r_eta: Any, D_z_eta: Any, beta_labels: Sequence[Any], *,
        A_c: Any, u_star: Any, h_plus: Any, h_minus: Any, T_N: Any, T_K: Any,
        D_r_A_c: Any = 0, D_r_u_star: Any = 0, D_r_h_plus: Any = 0,
        D_r_h_minus: Any = 0, D_r_T_N: Any = 0, D_r_T_K: Any = 0,
        D_z_A_c: Any = 0, D_z_u_star: Any = 0, D_z_h_plus: Any = 0,
        D_z_h_minus: Any = 0, D_z_T_N: Any = 0, D_z_T_K: Any = 0,
        direction_gap_eta: float | None = None,
    ) -> dict[str, Any]:
        part = self.partition_contract.validate_family(eta, D_r_eta, D_z_eta, beta_labels)
        beta_shape = part["eta"].shape
        sample_shape = beta_shape[:-1]
        sign_shape = beta_shape + (2,)
        vector_shape = beta_shape + (2, 3)

        phase = _real(phase, "phase")
        n_phi = _signed_vector(n_phi, "n_phi", True)
        t_plus_prototype = _signed_vector(t_plus_prototype, "t_plus_prototype", False)
        Dr_proto_all = _signed_vector(D_r_C_plus_prototype, "D_r_C_plus_prototype", False)
        Dz_proto_all = _signed_vector(D_z_C_plus_prototype, "D_z_C_plus_prototype", False)
        if phase.shape != sign_shape:
            raise ValueError("phase must have shape eta.shape+(2,)")
        for name, value in (("n_phi", n_phi), ("t_plus_prototype", t_plus_prototype),
                            ("D_r_C_plus_prototype", Dr_proto_all),
                            ("D_z_C_plus_prototype", Dz_proto_all)):
            if value.shape != vector_shape:
                raise ValueError(f"{name} must have shape eta.shape+(2,3)")

        R = _real(R, "R")
        try:
            R_beta = np.broadcast_to(R, beta_shape)
        except ValueError:
            if R.shape != sample_shape:
                raise ValueError("R must broadcast to sample_shape+(beta,)") from None
            R_beta = np.broadcast_to(R[..., None], beta_shape)
        if np.any(R_beta <= 0):
            raise ValueError("source cylindrical complete curl requires R>0")
        try:
            theta = np.broadcast_to(_real(theta, "theta"), sample_shape)
        except ValueError:
            raise ValueError("theta must broadcast to sample axes only") from None

        names = ("A_c","u_star","h_plus","h_minus","T_N","T_K",
                 "D_r_A_c","D_r_u_star","D_r_h_plus","D_r_h_minus","D_r_T_N","D_r_T_K",
                 "D_z_A_c","D_z_u_star","D_z_h_plus","D_z_h_minus","D_z_T_N","D_z_T_K")
        raw = (A_c,u_star,h_plus,h_minus,T_N,T_K,
               D_r_A_c,D_r_u_star,D_r_h_plus,D_r_h_minus,D_r_T_N,D_r_T_K,
               D_z_A_c,D_z_u_star,D_z_h_plus,D_z_h_minus,D_z_T_N,D_z_T_K)
        s = {name: self._beta_scalar(value, name, beta_shape) for name, value in zip(names, raw)}

        pair = self.covariance_pair
        base = pair.solve_reference_amplitudes(
            s["A_c"], s["u_star"], s["h_plus"], s["h_minus"], s["T_N"], s["T_K"],
            direction_gap_eta=direction_gap_eta)
        amplitudes = np.asarray(base["amplitudes"], dtype=float)
        Dr_a = self._scalar_directional_amplitudes(pair, s, "D_r", beta_shape)
        Dz_a = self._scalar_directional_amplitudes(pair, s, "D_z", beta_shape)

        ell_by_beta = tuple(int(label[0]) for label in part["beta_labels"])
        schedules = tuple(KokunoSourceBandCovering(ell, self.h) for ell in ell_by_beta)
        signed_cyl, scaled_t, scaled_DrC, scaled_DzC = [], [], [], []
        for j, schedule in enumerate(schedules):
            sign_fields, sign_t, sign_DrC, sign_DzC = [], [], [], []
            root_eps = float(np.sqrt(schedule.epsilon))
            for sigma_index in range(2):
                plus_loc = KokunoSourceSupportLocalizedCurl(
                    epsilon=schedule.epsilon, m=1, partition_atol=self.partition_atol)
                n = n_phi[..., j, sigma_index, :]
                t0 = t_plus_prototype[..., j, sigma_index, :]
                Dr0 = Dr_proto_all[..., j, sigma_index, :]
                Dz0 = Dz_proto_all[..., j, sigma_index, :]
                C0 = plus_loc.complete_curl.coefficient(n, t0)
                g = root_eps * amplitudes[..., j, sigma_index]
                Drg = root_eps * Dr_a[..., j, sigma_index]
                Dzg = root_eps * Dz_a[..., j, sigma_index]
                t = g[..., None] * t0
                DrC = g[..., None] * Dr0 + Drg[..., None] * C0
                DzC = g[..., None] * Dz0 + Dzg[..., None] * C0
                common = (R_beta[..., j], phase[..., j, sigma_index], n)
                plus = plus_loc.localized_mode(
                    *common, t, DrC, DzC, part["eta"][..., j],
                    part["D_r_eta"][..., j], part["D_z_eta"][..., j])
                minus = KokunoSourceSupportLocalizedCurl(
                    epsilon=schedule.epsilon, m=-1, partition_atol=self.partition_atol
                ).localized_mode(
                    *common, np.conjugate(t), np.conjugate(DrC), np.conjugate(DzC),
                    part["eta"][..., j], part["D_r_eta"][..., j], part["D_z_eta"][..., j])
                z = plus["velocity"] + minus["velocity"]
                scale = max(1.0, float(np.max(np.abs(z.real), initial=0.0)))
                if float(np.max(np.abs(z.imag), initial=0.0)) > self.reality_atol * scale:
                    raise RuntimeError("signed m=+/-1 pair failed the declared reality guard")
                if not np.allclose(z.real, 2 * plus["velocity"].real, rtol=0,
                                   atol=self.reality_atol * scale):
                    raise RuntimeError("signed real pair disagrees with 2*Re(m=+1)")
                sign_fields.append(schedule.Q ** (-self.A) * z.real)
                sign_t.append(t); sign_DrC.append(DrC); sign_DzC.append(DzC)
            signed_cyl.append(np.stack(sign_fields, axis=-2))
            scaled_t.append(np.stack(sign_t, axis=-2))
            scaled_DrC.append(np.stack(sign_DrC, axis=-2))
            scaled_DzC.append(np.stack(sign_DzC, axis=-2))

        cyl_bs = np.stack(signed_cyl, axis=-3)
        cart_bs = _rotate(cyl_bs, theta)
        cyl_b = np.sum(cyl_bs, axis=-2)
        cart_b = np.sum(cart_bs, axis=-2)
        order = _canonical_indices(part["beta_labels"])
        cyl_total = np.sum(np.take(cyl_b, order, axis=-2), axis=-2)
        cart_total = np.sum(np.take(cart_b, order, axis=-2), axis=-2)
        return {
            "beta_labels": part["beta_labels"], "sigma_labels": pair.sign_labels(),
            "ell_by_beta": ell_by_beta, "canonical_beta_indices": order,
            "Q_by_beta": np.asarray([q.Q for q in schedules]),
            "epsilon_by_beta": np.asarray([q.epsilon for q in schedules]), "A": self.A,
            "partition_error": part["partition_error"],
            "reference_squared_amplitudes": np.asarray(base["squared_amplitudes"]),
            "reference_amplitudes": amplitudes,
            "D_r_reference_amplitudes": Dr_a, "D_z_reference_amplitudes": Dz_a,
            "scaled_t_plus_prototype": np.stack(scaled_t, axis=-3),
            "scaled_D_r_C_plus_prototype": np.stack(scaled_DrC, axis=-3),
            "scaled_D_z_C_plus_prototype": np.stack(scaled_DzC, axis=-3),
            "velocity_physical_cylindrical_by_beta_sign": cyl_bs,
            "velocity_physical_cartesian_by_beta_sign": cart_bs,
            "velocity_physical_cylindrical_by_beta": cyl_b,
            "velocity_physical_cartesian_by_beta": cart_b,
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
                       "source_formulas": dict(_SOURCE_FORMULAS)},
            "derived_product_rules": dict(_DERIVED_PRODUCT_RULES),
            "parameters": {"h": self.h, "partition_atol": self.partition_atol,
                           "derivative_rtol": self.derivative_rtol,
                           "reality_atol": self.reality_atol,
                           "covariance_consistency_atol": self.covariance_consistency_atol},
            "caller_contract": {
                "reference_target_units": "T_N,T_K are H_ref y=T units; physical Delta C requires T=Delta C/epsilon_ell",
                "axis": "R>0; no unsourced axis regularization",
                "partition": "one beta=(ell,a) entry, never duplicated over sigma",
                "actual_source_inputs": "background, h_sigma, signed rectangles/modes and partition labels remain caller supplied",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode()).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload(); payload["sha256"] = self.sha256
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return path

    @classmethod
    def from_payload(cls, p: dict[str, Any]) -> "KokunoSourceSignedCompleteCurlFamily":
        if not isinstance(p, dict):
            raise ValueError("signed complete-curl family payload must be an object")
        params = p.get("parameters", {})
        obj = cls(params.get("h"), params.get("partition_atol"), params.get("derivative_rtol"),
                  params.get("reality_atol"), params.get("covariance_consistency_atol"))
        expected = obj.to_payload()
        if set(p) - {"sha256"} != set(expected):
            raise ValueError("signed complete-curl family payload keys or schema do not match")
        for key in expected:
            if p[key] != expected[key]:
                raise ValueError(f"signed complete-curl family {key} metadata changed")
        if "sha256" in p and p["sha256"] != obj.sha256:
            raise ValueError("signed complete-curl family sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceSignedCompleteCurlFamily":
        return cls.from_payload(json.loads(Path(path).read_text()))
