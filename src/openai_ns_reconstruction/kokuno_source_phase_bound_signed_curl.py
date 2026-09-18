"""Bind Kokuno's displayed phase/frame formulas to the signed complete-curl family.

The corrected 2026-09-09 reader specifies, for each slow box beta=(ell,a)
and auxiliary-rectangle sign sigma=+/- ,

    k = ceil(epsilon**(-1/2)),
    s(v) = sigma*(u_*/2 + u_* v/L_s),
    x0 = sigma B_s u_*/2,
    B_s**2 = lambda0 / (epsilon*k**2*(1+u_***2)**(3/2)),
    (tilde_p/R0, p_z)
        = B_s*(K - sigma*u_* g0/(L_s*|g0|**2)),
    Phi = p theta + p_z Z/epsilon + x0 R - v H_Phi,
    H_Phi = p F + p_z G,
    n_Phi = (x0-v(H_Phi)_R, p/R, p_z-epsilon*v(H_Phi)_Z),

with p=(kp)/k and kp a nearest nonzero integer to k*tilde_p.

This module makes that phase/gradient seam executable and passes the derived
``Phi`` and ``n_Phi`` directly into the already audited signed complete-curl
family. The positive-order background values F,G and their derivatives, the
homogeneous pulse prototypes, slow partition, rectangle pulse coordinate v,
and source h_sigma integrals are still caller supplied. Therefore this is a
source-formula-bound adapter, not a recovered paper-exact oscillatory field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_source_band_covering import KokunoSourceBandCovering
from .kokuno_source_signed_complete_curl_family import KokunoSourceSignedCompleteCurlFamily

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-phase-bound-signed-complete-curl-v1"

_SOURCE_FORMULAS = {
    "cone_frame": (
        "t_s=-b_s/a; N=-(1,t_s)/sqrt(1+t_s^2); "
        "K=(t_s,-1)/sqrt(1+t_s^2); |g0|=a*F0*sqrt(1+t_s^2)"
    ),
    "lambda_c0": (
        "lambda0^2=2*a*F0^2*(1-2/v_s); "
        "c0=lambda0/(2*F0*N_theta)<0; v_s=a*(1+t_s^2)"
    ),
    "carrier": "k=ceil(epsilon^(-1/2)); 1<=epsilon*k^2<=4",
    "B_s": "B_s^2=lambda0/(epsilon*k^2*(1+u_star^2)^(3/2))",
    "signed_path": "s(v)=sigma*(u_star/2+u_star*v/L_s); x0=sigma*B_s*u_star/2",
    "phase_tilt": (
        "(tilde_p/R0,p_z)=B_s*(K-sigma*u_star*g0/(L_s*|g0|^2)); "
        "p=(kp)/k with kp nearest nonzero integer to k*tilde_p"
    ),
    "phase": "Phi=p*theta+p_z*Z/epsilon+x0*R-v*(p*F+p_z*G)",
    "phase_gradient": (
        "n_Phi=(x0-v*(p*F_R+p_z*G_R), p/R, "
        "p_z-epsilon*v*(p*F_Z+p_z*G_Z))"
    ),
}

_DERIVED_NUMERICS = {
    "rounding": (
        "repository deterministic nearest-nonzero-integer tie policy; "
        "away-from-zero on half-integer ties, +1 only for an exact zero tie"
    ),
    "bridge": (
        "derived Phi and n_Phi are passed unchanged into "
        "KokunoSourceSignedCompleteCurlFamily.physical_family"
    ),
}

_TRUTH_BOUNDARY = {
    "source_phase_formula_bound": True,
    "source_phase_gradient_formula_bound": True,
    "source_carrier_integer_rule_executable": True,
    "signed_phase_is_not_caller_free": True,
    "q_scaled_signed_complete_curl_bridge_executable": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_partition_labels_instantiated": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "genuinely_independent_second_covariance_column_ready": False,
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


def _nearest_nonzero_integer(x: np.ndarray) -> np.ndarray:
    """Deterministic nearest nonzero integer with an explicit tie rule."""
    x = _real(x, "k*tilde_p")
    nearest = np.where(x >= 0.0, np.floor(x + 0.5), np.ceil(x - 0.5)).astype(np.int64)
    zero = nearest == 0
    if np.any(zero):
        nearest = nearest.copy()
        nearest[zero] = np.where(x[zero] < 0.0, -1, 1)
    return nearest


@dataclass(frozen=True)
class KokunoSourcePhaseBoundSignedCurlFamily:
    h: float = 0.005
    partition_atol: float = 2e-12
    derivative_rtol: float = 2e-10
    reality_atol: float = 2e-12
    covariance_consistency_atol: float = 2e-12

    def __post_init__(self) -> None:
        family = KokunoSourceSignedCompleteCurlFamily(
            h=self.h,
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
            reality_atol=self.reality_atol,
            covariance_consistency_atol=self.covariance_consistency_atol,
        )
        object.__setattr__(self, "h", family.h)
        object.__setattr__(self, "partition_atol", family.partition_atol)
        object.__setattr__(self, "derivative_rtol", family.derivative_rtol)
        object.__setattr__(self, "reality_atol", family.reality_atol)
        object.__setattr__(
            self, "covariance_consistency_atol", family.covariance_consistency_atol
        )

    @property
    def complete_curl_family(self) -> KokunoSourceSignedCompleteCurlFamily:
        return KokunoSourceSignedCompleteCurlFamily(
            h=self.h,
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
            reality_atol=self.reality_atol,
            covariance_consistency_atol=self.covariance_consistency_atol,
        )

    @staticmethod
    def _to_beta(x: Any, name: str, beta_shape: tuple[int, ...]) -> np.ndarray:
        try:
            return np.broadcast_to(_real(x, name), beta_shape)
        except ValueError:
            raise ValueError(f"{name} must broadcast to sample_shape+(beta,)") from None

    @staticmethod
    def _to_beta_sign(x: Any, name: str, sign_shape: tuple[int, ...]) -> np.ndarray:
        try:
            return np.broadcast_to(_real(x, name), sign_shape)
        except ValueError:
            raise ValueError(f"{name} must broadcast to sample_shape+(beta,2)") from None

    def phase_frame(
        self,
        *,
        R: Any,
        Z: Any,
        theta: Any,
        v: Any,
        beta_labels: Sequence[Any],
        R0: Any,
        F0: Any,
        a: Any,
        b_s: Any,
        u_star: Any,
        L_s: Any,
        F: Any,
        G: Any,
        F_R: Any,
        G_R: Any,
        F_Z: Any,
        G_Z: Any,
    ) -> dict[str, Any]:
        labels = tuple(tuple(label) for label in beta_labels)
        if not labels:
            raise ValueError("beta_labels must not be empty")
        n_beta = len(labels)

        R_raw = _real(R, "R")
        if R_raw.ndim == 0:
            sample_shape: tuple[int, ...] = ()
            R_beta = np.broadcast_to(R_raw, (n_beta,))
        elif R_raw.shape[-1:] == (n_beta,):
            sample_shape = R_raw.shape[:-1]
            R_beta = R_raw
        else:
            sample_shape = R_raw.shape
            R_beta = np.broadcast_to(R_raw[..., None], sample_shape + (n_beta,))
        beta_shape = sample_shape + (n_beta,)
        sign_shape = beta_shape + (2,)

        Z_beta = self._to_beta(Z, "Z", beta_shape)
        R0_beta = self._to_beta(R0, "R0", beta_shape)
        F0_beta = self._to_beta(F0, "F0", beta_shape)
        a_beta = self._to_beta(a, "a", beta_shape)
        bs_beta = self._to_beta(b_s, "b_s", beta_shape)
        us_beta = self._to_beta(u_star, "u_star", beta_shape)
        Ls_beta = self._to_beta(L_s, "L_s", beta_shape)
        F_beta = self._to_beta(F, "F", beta_shape)
        G_beta = self._to_beta(G, "G", beta_shape)
        FR_beta = self._to_beta(F_R, "F_R", beta_shape)
        GR_beta = self._to_beta(G_R, "G_R", beta_shape)
        FZ_beta = self._to_beta(F_Z, "F_Z", beta_shape)
        GZ_beta = self._to_beta(G_Z, "G_Z", beta_shape)
        v_sign = self._to_beta_sign(v, "v", sign_shape)

        try:
            theta_sample = np.broadcast_to(_real(theta, "theta"), sample_shape)
        except ValueError:
            raise ValueError("theta must broadcast to sample axes only") from None
        theta_beta_sign = np.broadcast_to(theta_sample[..., None, None], sign_shape)

        if np.any(R_beta <= 0.0):
            raise ValueError("source phase frame requires R>0")
        if np.any(R0_beta <= 0.0):
            raise ValueError("R0 must be strictly positive")
        if np.any(F0_beta <= 0.0) or np.any(a_beta <= 0.0):
            raise ValueError("source cone frame requires F0>0 and a>0")
        if np.any(us_beta <= 0.0) or np.any(Ls_beta <= 0.0):
            raise ValueError("source phase frame requires u_star>0 and L_s>0")
        if np.any(v_sign < 0.0) or np.any(v_sign > Ls_beta[..., None]):
            raise ValueError("pulse coordinate v must satisfy 0<=v<=L_s")

        t_s = -bs_beta / a_beta
        denom = np.sqrt(1.0 + t_s * t_s)
        N = np.stack((-1.0 / denom, -t_s / denom), axis=-1)
        K = np.stack((t_s / denom, -1.0 / denom), axis=-1)
        v_s = a_beta * (1.0 + t_s * t_s)
        if np.any(v_s <= 2.0):
            raise ValueError("source phase cone requires v_s=a*(1+t_s^2)>2")

        g0_abs = a_beta * F0_beta * denom
        g0 = g0_abs[..., None] * N
        lambda_sq = 2.0 * a_beta * F0_beta * F0_beta * (1.0 - 2.0 / v_s)
        if np.any(lambda_sq <= 0.0):
            raise ValueError("source lambda0^2 must be positive")
        lambda0 = np.sqrt(lambda_sq)
        c0 = lambda0 / (2.0 * F0_beta * N[..., 0])
        if np.any(c0 >= 0.0):
            raise RuntimeError("source cone orientation must give c0<0")

        phase = np.empty(sign_shape, dtype=float)
        n_phi = np.empty(sign_shape + (3,), dtype=float)
        k_by_beta = np.empty(beta_shape, dtype=np.int64)
        kp = np.empty(sign_shape, dtype=np.int64)
        p = np.empty(sign_shape, dtype=float)
        p_z = np.empty(sign_shape, dtype=float)
        tilde_p = np.empty(sign_shape, dtype=float)
        B_s = np.empty(beta_shape, dtype=float)
        x0 = np.empty(sign_shape, dtype=float)
        signed_s = np.empty(sign_shape, dtype=float)
        rounding_error = np.empty(sign_shape, dtype=float)

        sigma_values = (1.0, -1.0)
        for j, label in enumerate(labels):
            ell = int(label[0])
            schedule = KokunoSourceBandCovering(ell, self.h)
            k = int(np.ceil(schedule.epsilon ** -0.5))
            if k < 1:
                raise RuntimeError("source carrier integer k must be positive")
            k_by_beta[..., j] = k
            Bs = np.sqrt(
                lambda0[..., j]
                / (
                    schedule.epsilon
                    * (k * k)
                    * (1.0 + us_beta[..., j] ** 2) ** 1.5
                )
            )
            B_s[..., j] = Bs

            for sidx, sigma in enumerate(sigma_values):
                sv = sigma * (
                    us_beta[..., j] / 2.0
                    + us_beta[..., j] * v_sign[..., j, sidx] / Ls_beta[..., j]
                )
                signed_s[..., j, sidx] = sv
                xx0 = sigma * Bs * us_beta[..., j] / 2.0
                x0[..., j, sidx] = xx0

                tilt = (
                    K[..., j, :]
                    - sigma
                    * us_beta[..., j, None]
                    * g0[..., j, :]
                    / (Ls_beta[..., j, None] * g0_abs[..., j, None] ** 2)
                )
                pair = Bs[..., None] * tilt
                tp = R0_beta[..., j] * pair[..., 0]
                pz = pair[..., 1]
                carrier = _nearest_nonzero_integer(k * tp)
                pp = carrier.astype(float) / float(k)

                tilde_p[..., j, sidx] = tp
                kp[..., j, sidx] = carrier
                p[..., j, sidx] = pp
                p_z[..., j, sidx] = pz
                rounding_error[..., j, sidx] = np.abs(pp - tp)

                H = pp * F_beta[..., j] + pz * G_beta[..., j]
                H_R = pp * FR_beta[..., j] + pz * GR_beta[..., j]
                H_Z = pp * FZ_beta[..., j] + pz * GZ_beta[..., j]
                vv = v_sign[..., j, sidx]
                phase[..., j, sidx] = (
                    pp * theta_beta_sign[..., j, sidx]
                    + pz * Z_beta[..., j] / schedule.epsilon
                    + xx0 * R_beta[..., j]
                    - vv * H
                )
                n_phi[..., j, sidx, 0] = xx0 - vv * H_R
                n_phi[..., j, sidx, 1] = pp / R_beta[..., j]
                n_phi[..., j, sidx, 2] = pz - schedule.epsilon * vv * H_Z

                if np.any(rounding_error[..., j, sidx] > 1.0 / k + 1e-15):
                    raise RuntimeError("nearest-nonzero carrier violates source |p-tilde_p|<=1/k")

        tangent_norm = np.linalg.norm(n_phi[..., 1:], axis=-1)
        if np.any(tangent_norm <= 0.0):
            raise ValueError("source phase tangential gradient must remain nonzero")

        return {
            "beta_labels": labels,
            "sigma_labels": self.complete_curl_family.covariance_pair.sign_labels(),
            "phase": phase,
            "n_phi": n_phi,
            "k_by_beta": k_by_beta,
            "kp_by_beta_sign": kp,
            "p_by_beta_sign": p,
            "p_z_by_beta_sign": p_z,
            "tilde_p_by_beta_sign": tilde_p,
            "rounding_error_by_beta_sign": rounding_error,
            "B_s_by_beta": B_s,
            "s_by_beta_sign": signed_s,
            "x0_by_beta_sign": x0,
            "N_by_beta": N,
            "K_by_beta": K,
            "g0_by_beta": g0,
            "g0_abs_by_beta": g0_abs,
            "lambda0_by_beta": lambda0,
            "c0_by_beta": c0,
            "v_s_by_beta": v_s,
            "source_phase_formula_bound": True,
            "background_values_caller_supplied": True,
        }

    def physical_family(
        self,
        *,
        R: Any,
        Z: Any,
        theta: Any,
        v: Any,
        beta_labels: Sequence[Any],
        R0: Any,
        F0: Any,
        a: Any,
        b_s: Any,
        u_star: Any,
        L_s: Any,
        F: Any,
        G: Any,
        F_R: Any,
        G_R: Any,
        F_Z: Any,
        G_Z: Any,
        t_plus_prototype: Any,
        D_r_C_plus_prototype: Any,
        D_z_C_plus_prototype: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        A_c: Any,
        h_plus: Any,
        h_minus: Any,
        T_N: Any,
        T_K: Any,
        D_r_A_c: Any = 0,
        D_r_u_star: Any = 0,
        D_r_h_plus: Any = 0,
        D_r_h_minus: Any = 0,
        D_r_T_N: Any = 0,
        D_r_T_K: Any = 0,
        D_z_A_c: Any = 0,
        D_z_u_star: Any = 0,
        D_z_h_plus: Any = 0,
        D_z_h_minus: Any = 0,
        D_z_T_N: Any = 0,
        D_z_T_K: Any = 0,
        direction_gap_eta: float | None = None,
    ) -> dict[str, Any]:
        frame = self.phase_frame(
            R=R, Z=Z, theta=theta, v=v, beta_labels=beta_labels,
            R0=R0, F0=F0, a=a, b_s=b_s, u_star=u_star, L_s=L_s,
            F=F, G=G, F_R=F_R, G_R=G_R, F_Z=F_Z, G_Z=G_Z,
        )
        result = self.complete_curl_family.physical_family(
            R=R,
            theta=theta,
            phase=frame["phase"],
            n_phi=frame["n_phi"],
            t_plus_prototype=t_plus_prototype,
            D_r_C_plus_prototype=D_r_C_plus_prototype,
            D_z_C_plus_prototype=D_z_C_plus_prototype,
            eta=eta,
            D_r_eta=D_r_eta,
            D_z_eta=D_z_eta,
            beta_labels=beta_labels,
            A_c=A_c,
            u_star=u_star,
            h_plus=h_plus,
            h_minus=h_minus,
            T_N=T_N,
            T_K=T_K,
            D_r_A_c=D_r_A_c,
            D_r_u_star=D_r_u_star,
            D_r_h_plus=D_r_h_plus,
            D_r_h_minus=D_r_h_minus,
            D_r_T_N=D_r_T_N,
            D_r_T_K=D_r_T_K,
            D_z_A_c=D_z_A_c,
            D_z_u_star=D_z_u_star,
            D_z_h_plus=D_z_h_plus,
            D_z_h_minus=D_z_h_minus,
            D_z_T_N=D_z_T_N,
            D_z_T_K=D_z_T_K,
            direction_gap_eta=direction_gap_eta,
        )
        return {
            **result,
            "source_phase_frame": frame,
            "source_phase_formula_bound": True,
            "source_phase_gradient_formula_bound": True,
            "actual_positive_order_background_bound": False,
            "actual_signed_auxiliary_rectangles_bound": False,
            "public_xyz_t_velocity_correction_materialized": False,
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
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "derived_numerics": dict(_DERIVED_NUMERICS),
            "parameters": {
                "h": self.h,
                "partition_atol": self.partition_atol,
                "derivative_rtol": self.derivative_rtol,
                "reality_atol": self.reality_atol,
                "covariance_consistency_atol": self.covariance_consistency_atol,
            },
            "caller_contract": {
                "background": (
                    "F,G,F_R,G_R,F_Z,G_Z and representative F0,a,b_s,R0 "
                    "remain caller supplied; this adapter does not reconstruct Agent-1 background"
                ),
                "pulse_coordinate": (
                    "v and L_s remain caller supplied rectangle/pulse data; "
                    "actual source rectangles are not claimed recovered"
                ),
                "pulse_prototype": (
                    "t_plus and prototype D_r C_plus/D_z C_plus remain caller supplied"
                ),
                "partition": "eta_beta and derivatives remain caller supplied",
                "units": (
                    "R,Z are source wave-chart coordinates; returned velocity is the "
                    "Q_ell^(-A)-scaled physical Cartesian/cylindrical correction from #482"
                ),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode()).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return path

    @classmethod
    def from_payload(cls, p: dict[str, Any]) -> "KokunoSourcePhaseBoundSignedCurlFamily":
        if not isinstance(p, dict):
            raise ValueError("source phase-bound payload must be an object")
        params = p.get("parameters", {})
        obj = cls(
            params.get("h"),
            params.get("partition_atol"),
            params.get("derivative_rtol"),
            params.get("reality_atol"),
            params.get("covariance_consistency_atol"),
        )
        expected = obj.to_payload()
        if set(p) - {"sha256"} != set(expected):
            raise ValueError("source phase-bound payload keys or schema do not match")
        for key in expected:
            if p[key] != expected[key]:
                raise ValueError(f"source phase-bound {key} metadata changed")
        if "sha256" in p and p["sha256"] != obj.sha256:
            raise ValueError("source phase-bound sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourcePhaseBoundSignedCurlFamily":
        return cls.from_payload(json.loads(Path(path).read_text()))
