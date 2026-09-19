"""Public, deterministic Kokuno-structured candidate oscillatory velocity.

This module binds the already-audited Agent-2 source-formula pipeline to one
fully specified repository-autonomous numerical realization and exposes a
callable ``velocity_osc(x,y,z,t)``.  It does *not* recover Kokuno's unpublished
positive-order background, pulse data, auxiliary mode family, or partition
realization.  Those choices remain explicitly candidate-only.

Source-displayed parts reused here are the per-band schedule, signed covariance
inverse, phase/frame equations, vector-potential complete curl, real m=+/-1
pairing, and Q_ell physical scaling.  The frozen numerical background,
radial/axial compact envelopes, pulse/cutoff samples, transverse mode prototype,
slow partition realization and time modulation are repository choices.  No
forcing or pressure is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_autonomous_signed_rectangle_geometry import KokunoAutonomousSignedRectangleGeometry
from .kokuno_candidate_mass_bound_phase_curl import KokunoCandidateMassBoundPhaseCurlFamily
from .kokuno_source_compatible_partition import KokunoSourceCompatiblePartitionRealization
from .kokuno_source_signed_covariance_mass import KokunoSourceSignedCovarianceMass

SCHEMA = "kokuno-public-candidate-oscillatory-velocity-v2"
PARENT_HEAD = "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e"
PROJECT_AXIAL_LIMIT = 2.0

_SOURCE_BOUND = {
    "band_schedule": "Q_ell=2^(-ell), epsilon_ell=Q_ell^h and source covering schedule",
    "signed_covariance": "H_ref y=T with sigma=+/- and y_sigma=a_sigma^2",
    "phase_frame": "corrected-reader signed phase/frame and carrier construction",
    "curl": "u_osc is assembled from localized vector potentials through the complete cylindrical curl",
    "reality": "m=+1 and m=-1 are conjugately paired before Cartesian conversion",
    "physical_scaling": "each beta contribution is multiplied by Q_ell^(-(1/2+h))",
}

_AUTONOMOUS = {
    "slow_partition": "KokunoSourceCompatiblePartitionRealization evaluated at frozen q=2^-5.5 and slow=(0,0,0)",
    "signed_rectangles": "deterministic rational source-compatible witness, not recovered source rectangles",
    "background": "constant candidate F,G,F0,a,b_s,u_star,R0 values with zero slow spatial derivatives",
    "pulse_mass": "analytic positive signed pulse samples plus compact numerical cutoff, integrated by the existing mass adapter",
    "mode": "signed transverse prototype t=alpha*x_sigma(v/L_s)*E_R(R)*E_Z(Z)*(n_theta,-n_r,0)",
    "mode_derivative": "analytic ordinary radial and axial derivatives of C=i(n x t)/(k|n|^2); auxiliary-independent candidate route",
    "support": "repository C7 cos(pi*s/2)^8 radial and axial envelopes, applied at vector-potential coefficient level before complete curl",
    "time": "bounded autonomous modulation of covariance targets only; spatial derivative inputs remain unchanged",
}

_TRUTH = {
    "public_xyz_t_velocity_correction_materialized": True,
    "candidate_autonomous_background_bound": True,
    "candidate_autonomous_signed_mode_bound": True,
    "candidate_autonomous_partition_bound": True,
    "candidate_project_axial_support_bound": True,
    "vector_potential_complete_curl_path_used": True,
    "axis_safe_by_compact_radial_support": True,
    "project_axial_support_applied_before_curl": True,
    "source_formula_chain_preserved": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_partition_labels_instantiated": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _real(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _pulse_x(s: np.ndarray) -> np.ndarray:
    """Positive autonomous signed pulse prototype on normalized v/L_s in [0,1]."""
    s = np.asarray(s, dtype=float)
    c = np.cos(2.0 * math.pi * (s - 0.5))
    return np.stack((2.0 + 0.15 * c, 3.0 - 0.12 * c), axis=-1)


def _compact_cos8_envelope(
    coordinate: np.ndarray,
    center: float,
    halfwidth: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a C7 compact cos^8 envelope and ordinary coordinate derivative."""
    s = (np.asarray(coordinate, dtype=float) - center) / halfwidth
    value = np.zeros_like(s)
    derivative = np.zeros_like(s)
    mask = np.abs(s) < 1.0
    if np.any(mask):
        q = 0.5 * math.pi * s[mask]
        c = np.cos(q)
        value[mask] = c**8
        derivative[mask] = -4.0 * math.pi * np.sin(q) * c**7 / halfwidth
    return value, derivative


def _radial_envelope(R: np.ndarray, center: float, halfwidth: float) -> tuple[np.ndarray, np.ndarray]:
    """Backward-compatible radial envelope helper used by focused tests."""
    return _compact_cos8_envelope(R, center, halfwidth)


@dataclass(frozen=True)
class KokunoPublicCandidateOscillatoryVelocity:
    """One frozen, low-dimensional, provenance-labelled oscillatory candidate.

    Coordinates are repository wave-chart coordinates: ``R=hypot(x,y)``,
    ``theta=atan2(y,x)``, and ``Z=z``.  They are dimensionless here.  The
    construction is identically zero outside its radial annulus and outside the
    registered project axial support, and therefore never evaluates the
    cylindrical source formulas on the axis or beyond the compact support.
    """

    h: float = 0.005
    radial_center: float = 0.75
    radial_halfwidth: float = 0.60
    axial_center: float = 0.0
    axial_halfwidth: float = 2.0
    time_min: float = 0.25
    time_max: float = 0.75
    normal_target_ratio: float = 0.40
    cross_target_ratio: float = 0.05
    time_modulation: float = 0.08
    mode_imaginary_ratio: float = 0.17
    pulse_samples: int = 129
    transverse_samples: int = 257

    def __post_init__(self) -> None:
        values = {
            "h": self.h,
            "radial_center": self.radial_center,
            "radial_halfwidth": self.radial_halfwidth,
            "axial_center": self.axial_center,
            "axial_halfwidth": self.axial_halfwidth,
            "time_min": self.time_min,
            "time_max": self.time_max,
            "normal_target_ratio": self.normal_target_ratio,
            "cross_target_ratio": self.cross_target_ratio,
            "time_modulation": self.time_modulation,
            "mode_imaginary_ratio": self.mode_imaginary_ratio,
        }
        for name, value in values.items():
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        if not (0.0 < self.h < 0.01):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        if not (self.radial_center > self.radial_halfwidth > 0.0):
            raise ValueError("radial support must stay strictly away from the cylindrical axis")
        if not self.axial_halfwidth > 0.0:
            raise ValueError("axial_halfwidth must be positive")
        if abs(self.axial_center) + self.axial_halfwidth > PROJECT_AXIAL_LIMIT:
            raise ValueError("axial vector-potential support must stay within registered |z|<2 support")
        if not self.time_min < self.time_max:
            raise ValueError("time_min must be less than time_max")
        if not (0.10 <= self.normal_target_ratio <= 0.80):
            raise ValueError("normal_target_ratio must stay in [0.10,0.80]")
        if not (0.0 < self.cross_target_ratio <= 0.20):
            raise ValueError("cross_target_ratio must stay in (0,0.20]")
        if not (0.0 <= self.time_modulation <= 0.20):
            raise ValueError("time_modulation must stay in [0,0.20]")
        if self.cross_target_ratio * (1.0 + self.time_modulation) >= self.normal_target_ratio * (1.0 - self.time_modulation):
            raise ValueError("signed covariance target must retain a strict positive amplitude gap")
        if abs(self.mode_imaginary_ratio) > 0.5:
            raise ValueError("mode_imaginary_ratio must stay in [-0.5,0.5]")
        for name, value, minimum in (
            ("pulse_samples", self.pulse_samples, 33),
            ("transverse_samples", self.transverse_samples, 65),
        ):
            if isinstance(value, bool) or int(value) != value or int(value) < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        object.__setattr__(self, "pulse_samples", int(self.pulse_samples))
        object.__setattr__(self, "transverse_samples", int(self.transverse_samples))

    @property
    def radial_inner(self) -> float:
        return self.radial_center - self.radial_halfwidth

    @property
    def radial_outer(self) -> float:
        return self.radial_center + self.radial_halfwidth

    @property
    def axial_lower(self) -> float:
        return self.axial_center - self.axial_halfwidth

    @property
    def axial_upper(self) -> float:
        return self.axial_center + self.axial_halfwidth

    @cached_property
    def _static(self) -> dict[str, Any]:
        partition = KokunoSourceCompatiblePartitionRealization(ell_min=5).evaluate(
            q=np.asarray(2.0**-5.5),
            D_r_q=np.asarray(0.0),
            D_z_q=np.asarray(0.0),
            slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
            D_r_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
            D_z_slow_coordinates=np.asarray((0.0, 0.0, 0.0)),
        )
        geometry = KokunoAutonomousSignedRectangleGeometry(h=self.h).instantiate(partition)
        v_by_beta: list[np.ndarray] = []
        psi_by_beta: list[np.ndarray] = []
        x_by_beta_sign: list[np.ndarray] = []
        for L_s in np.asarray(geometry["L_s_by_beta"], dtype=float):
            s = np.linspace(0.0, 1.0, self.pulse_samples)
            v_by_beta.append(float(L_s) * s)
            psi_by_beta.append(np.sin(math.pi * s) ** 4)
            x_by_beta_sign.append(_pulse_x(s))
        xi = np.linspace(-1.0, 1.0, self.transverse_samples)
        chi_g = np.cos(0.5 * math.pi * xi) ** 8
        mass = KokunoSourceSignedCovarianceMass(h=self.h).materialize_from_autonomous_geometry(
            geometry,
            v_by_beta=tuple(v_by_beta),
            psi_by_beta=tuple(psi_by_beta),
            x_by_beta_sign=tuple(x_by_beta_sign),
            xi=xi,
            chi_g=chi_g,
            pulse_binding="repository_autonomous_source_compatible",
        )
        return {"partition": partition, "geometry": geometry, "mass": mass}

    @property
    def beta_labels(self) -> tuple[tuple[Any, ...], ...]:
        return tuple(self._static["geometry"]["beta_labels"])

    @property
    def bridge(self) -> KokunoCandidateMassBoundPhaseCurlFamily:
        return KokunoCandidateMassBoundPhaseCurlFamily(h=self.h)

    def _background(self) -> dict[str, Any]:
        n = len(self.beta_labels)
        return {
            "R0": np.full(n, 0.82),
            "F0": np.full(n, 1.20),
            "a": np.full(n, 3.20),
            "b_s": np.full(n, 0.20),
            "u_star": np.full(n, 2.10),
            "F": np.full(n, 1.10),
            "G": np.full(n, 0.24),
            "F_R": np.zeros(n),
            "G_R": np.zeros(n),
            "F_Z": np.zeros(n),
            "G_Z": np.zeros(n),
        }

    def _inside_family(self, R: np.ndarray, theta: np.ndarray, Z: np.ndarray, t: np.ndarray) -> dict[str, Any]:
        """Evaluate points strictly inside both compact vector-potential supports."""
        labels = self.beta_labels
        n_beta = len(labels)
        geometry = self._static["geometry"]
        mass = self._static["mass"]
        background = self._background()
        m = R.size
        R_beta = np.broadcast_to(R[:, None], (m, n_beta))
        Z_beta = np.broadcast_to(Z[:, None], (m, n_beta))
        L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
        v = np.stack((0.25 * L_s, 0.75 * L_s), axis=-1)
        frame_inputs = {
            "R": R_beta,
            "Z": Z_beta,
            "theta": theta,
            "v": v,
            "beta_labels": labels,
            "L_s": L_s,
            **background,
        }
        frame = self.bridge.phase_family.phase_frame(**frame_inputs)
        n_phi = np.asarray(frame["n_phi"], dtype=float)

        radial_envelope, D_r_radial_envelope = _compact_cos8_envelope(
            R, self.radial_center, self.radial_halfwidth
        )
        axial_envelope, D_z_axial_envelope = _compact_cos8_envelope(
            Z, self.axial_center, self.axial_halfwidth
        )
        envelope = radial_envelope * axial_envelope
        D_r_envelope = D_r_radial_envelope * axial_envelope
        D_z_envelope = radial_envelope * D_z_axial_envelope

        pulse_factors = _pulse_x(np.asarray((0.25, 0.75)))[np.arange(2), np.arange(2)]
        alpha = complex(1.0, self.mode_imaginary_ratio)
        base_direction = np.stack(
            (n_phi[..., 1], -n_phi[..., 0], np.zeros_like(n_phi[..., 0])), axis=-1
        ).astype(np.complex128)
        base_t = alpha * pulse_factors[None, None, :, None] * base_direction

        D_r_n = np.zeros_like(n_phi)
        D_r_n[..., 1] = -n_phi[..., 1] / R[:, None, None]
        D_r_base_direction = np.stack(
            (D_r_n[..., 1], -D_r_n[..., 0], np.zeros_like(D_r_n[..., 0])), axis=-1
        ).astype(np.complex128)
        D_r_base_t = alpha * pulse_factors[None, None, :, None] * D_r_base_direction

        n_norm_sq = np.sum(n_phi * n_phi, axis=-1)
        cross = np.cross(n_phi, base_t)
        D_r_cross = np.cross(D_r_n, base_t) + np.cross(n_phi, D_r_base_t)
        D_r_norm_sq = 2.0 * np.sum(n_phi * D_r_n, axis=-1)
        k = np.asarray(frame["k_by_beta"], dtype=float)
        denom = k[..., None, None] * n_norm_sq[..., None]
        C_base = 1j * cross / denom
        D_r_C_base = 1j * (
            D_r_cross / denom
            - cross * D_r_norm_sq[..., None] / (k[..., None, None] * n_norm_sq[..., None] ** 2)
        )

        t_plus = envelope[:, None, None, None] * base_t
        D_r_C = (
            D_r_envelope[:, None, None, None] * C_base
            + envelope[:, None, None, None] * D_r_C_base
        )
        # For this frozen constant background, the unlocalized C_base has no
        # ordinary Z dependence.  The axial vector-potential localization is
        # therefore the complete ordinary Z coefficient derivative.
        D_z_C = D_z_envelope[:, None, None, None] * C_base

        eta0 = np.asarray(geometry["eta"], dtype=float)
        D_r_eta0 = np.asarray(geometry["D_r_eta"], dtype=float)
        D_z_eta0 = np.asarray(geometry["D_z_eta"], dtype=float)
        eta = np.broadcast_to(eta0, (m, n_beta))
        D_r_eta = np.broadcast_to(D_r_eta0, (m, n_beta))
        D_z_eta = np.broadcast_to(D_z_eta0, (m, n_beta))

        u_star = np.asarray(background["u_star"], dtype=float)
        A_c = -np.asarray(frame["c0_by_beta"], dtype=float) * np.sqrt(1.0 + u_star[None, :] ** 2)
        tau = (t - self.time_min) / (self.time_max - self.time_min)
        normal_ratio = self.normal_target_ratio * (1.0 + self.time_modulation * np.sin(2.0 * math.pi * tau))
        cross_ratio = self.cross_target_ratio * (1.0 + self.time_modulation * np.cos(2.0 * math.pi * tau))
        T_N = -normal_ratio[:, None] * A_c
        T_K = cross_ratio[:, None] * u_star[None, :]

        zeros_h = np.zeros((n_beta, 2), dtype=float)
        out = self.bridge.physical_family_from_materialized_mass(
            mass,
            D_r_h_by_beta_sign=zeros_h,
            D_z_h_by_beta_sign=zeros_h,
            h_derivative_binding="explicit_frozen_local_zero",
            **frame_inputs,
            t_plus_prototype=t_plus,
            D_r_C_plus_prototype=D_r_C,
            D_z_C_plus_prototype=D_z_C,
            eta=eta,
            D_r_eta=D_r_eta,
            D_z_eta=D_z_eta,
            A_c=A_c,
            T_N=T_N,
            T_K=T_K,
            D_r_A_c=np.zeros_like(A_c),
            D_r_u_star=np.zeros_like(A_c),
            D_r_T_N=np.zeros_like(A_c),
            D_r_T_K=np.zeros_like(A_c),
            D_z_A_c=np.zeros_like(A_c),
            D_z_u_star=np.zeros_like(A_c),
            D_z_T_N=np.zeros_like(A_c),
            D_z_T_K=np.zeros_like(A_c),
        )
        return {
            **out,
            "candidate_mode_C_plus_prototype": envelope[:, None, None, None] * C_base,
            "candidate_mode_D_r_C_plus_prototype": D_r_C,
            "candidate_mode_D_z_C_plus_prototype": D_z_C,
            "candidate_radial_envelope": radial_envelope,
            "candidate_radial_envelope_D_r": D_r_radial_envelope,
            "candidate_axial_envelope": axial_envelope,
            "candidate_axial_envelope_D_z": D_z_axial_envelope,
            "candidate_product_envelope": envelope,
            "candidate_product_envelope_D_r": D_r_envelope,
            "candidate_product_envelope_D_z": D_z_envelope,
        }

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
        """Evaluate total/by-beta/by-sign correction on broadcastable inputs.

        The public result is exactly zero outside the radial annulus or outside
        the registered axial support.  Both cutoffs are applied at the
        vector-potential coefficient level, before the complete curl; the
        corresponding coefficient derivatives are carried analytically.
        """
        x, y, z, t = np.broadcast_arrays(
            _real(x, "x"), _real(y, "y"), _real(z, "z"), _real(t, "t")
        )
        if np.any(t < self.time_min) or np.any(t > self.time_max):
            raise ValueError("t is outside the registered candidate interval")
        sample_shape = x.shape
        n_beta = len(self.beta_labels)
        total = np.zeros(sample_shape + (3,), dtype=float)
        by_beta = np.zeros(sample_shape + (n_beta, 3), dtype=float)
        by_beta_sign = np.zeros(sample_shape + (n_beta, 2, 3), dtype=float)
        R = np.hypot(x, y)
        inside = (
            (R > self.radial_inner)
            & (R < self.radial_outer)
            & (z > self.axial_lower)
            & (z < self.axial_upper)
        )

        indices = np.flatnonzero(inside.ravel())
        if indices.size:
            xf = x.ravel()[indices]
            yf = y.ravel()[indices]
            zf = z.ravel()[indices]
            tf = t.ravel()[indices]
            Rf = np.hypot(xf, yf)
            thetaf = np.arctan2(yf, xf)
            core = self._inside_family(Rf, thetaf, zf, tf)
            total.reshape((-1, 3))[indices] = np.asarray(core["candidate_velocity_osc_cartesian_total"], dtype=float)
            by_beta.reshape((-1, n_beta, 3))[indices] = np.asarray(core["candidate_velocity_osc_cartesian_by_beta"], dtype=float)
            by_beta_sign.reshape((-1, n_beta, 2, 3))[indices] = np.asarray(core["candidate_velocity_osc_cartesian_by_beta_sign"], dtype=float)

        return {
            "velocity_cartesian_total": total,
            "velocity_cartesian_by_beta": by_beta,
            "velocity_cartesian_by_beta_sign": by_beta_sign,
            "support_mask": inside,
            "beta_labels": self.beta_labels,
            "sign_labels": ("sigma_plus", "sigma_minus"),
            "coordinate_contract": "dimensionless repository wave chart: R=hypot(x,y), theta=atan2(y,x), Z=z",
            "time_interval": (self.time_min, self.time_max),
            "radial_support": (self.radial_inner, self.radial_outer),
            "axial_support": (self.axial_lower, self.axial_upper),
            "public_xyz_t_velocity_correction_materialized": True,
            "pde_validated": False,
            "paper_exact": False,
        }

    def velocity_osc(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the public candidate oscillatory correction as ``[...,3]``."""
        return self.evaluate(x, y, z, t)["velocity_cartesian_total"]

    def to_payload(self) -> dict[str, Any]:
        geometry = self._static["geometry"]
        mass = self._static["mass"]
        payload = {
            "schema": SCHEMA,
            "parent_agent2_head": PARENT_HEAD,
            "source_bound": dict(_SOURCE_BOUND),
            "autonomous_choices": dict(_AUTONOMOUS),
            "parameters": {
                "h": self.h,
                "radial_center": self.radial_center,
                "radial_halfwidth": self.radial_halfwidth,
                "axial_center": self.axial_center,
                "axial_halfwidth": self.axial_halfwidth,
                "time_min": self.time_min,
                "time_max": self.time_max,
                "normal_target_ratio": self.normal_target_ratio,
                "cross_target_ratio": self.cross_target_ratio,
                "time_modulation": self.time_modulation,
                "mode_imaginary_ratio": self.mode_imaginary_ratio,
                "pulse_samples": self.pulse_samples,
                "transverse_samples": self.transverse_samples,
            },
            "frozen_realization": {
                "beta_labels": [[int(label[0]), [int(v) for v in label[1]]] for label in self.beta_labels],
                "eta": np.asarray(geometry["eta"], dtype=float).tolist(),
                "L_s_by_beta": np.asarray(geometry["L_s_by_beta"], dtype=float).tolist(),
                "center_index_by_beta_sign": np.asarray(geometry["center_index_by_beta_sign"], dtype=int).tolist(),
                "h_sigma_by_beta_sign": np.asarray(mass["h_sigma_by_beta_sign"], dtype=float).tolist(),
                "background": {
                    "R0": 0.82,
                    "F0": 1.20,
                    "a": 3.20,
                    "b_s": 0.20,
                    "u_star": 2.10,
                    "F": 1.10,
                    "G": 0.24,
                    "F_R": 0.0,
                    "G_R": 0.0,
                    "F_Z": 0.0,
                    "G_Z": 0.0,
                },
                "pulse_sign_fractions": [0.25, 0.75],
                "mode_complex_factor": [1.0, self.mode_imaginary_ratio],
                "support_contract": {
                    "radial": [self.radial_inner, self.radial_outer],
                    "axial": [self.axial_lower, self.axial_upper],
                    "project_axial_limit": PROJECT_AXIAL_LIMIT,
                    "applied_at": "vector_potential_coefficient_before_complete_curl",
                },
            },
            "truth_boundary": dict(_TRUTH),
        }
        return payload

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return path

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPublicCandidateOscillatoryVelocity":
        if not isinstance(payload, dict):
            raise ValueError("public oscillatory candidate payload must be an object")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("public oscillatory candidate payload is missing parameters")
        obj = cls(**params)
        expected = obj.to_payload()
        body = {k: v for k, v in payload.items() if k != "sha256"}
        if body != expected:
            raise ValueError("public oscillatory candidate payload/provenance changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("public oscillatory candidate sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPublicCandidateOscillatoryVelocity":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


@lru_cache(maxsize=1)
def default_field() -> KokunoPublicCandidateOscillatoryVelocity:
    return KokunoPublicCandidateOscillatoryVelocity()


def velocity_osc(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Public default candidate correction: ``velocity_osc(x,y,z,t)->[...,3]``."""
    return default_field().velocity_osc(x, y, z, t)
