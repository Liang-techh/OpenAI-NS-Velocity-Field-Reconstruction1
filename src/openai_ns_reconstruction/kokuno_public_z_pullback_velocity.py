"""Public-coordinate pullback for the carrier-resolved Kokuno complete-curl candidate.

The corrected 2026-09-09 reconstruction uses the normalized axial operator

    D_z = epsilon * partial_Z

inside ``n_Phi = grad_* Phi`` and the stage-9 complete cylindrical curl.  The
previous public candidate identified the source chart coordinate directly with
the public coordinate (``Z=z``).  That is convenient for evaluating the source
formulas, but it does not make the source-normalized curl an ordinary Cartesian
curl with respect to public ``z``.

This candidate applies the minimal bandwise pullback

    Z_beta = epsilon_beta * z,

so that, for every beta band separately,

    partial_z = epsilon_beta * partial_{Z_beta} = D_z.

No source phase, carrier, covariance, complete-curl, Q/epsilon scaling, support,
or autonomous background formula is changed.  The radial/axial compact support
still lives in public coordinates and is still applied to the vector-potential
coefficient before curling.  This is a repository coordinate realization, not
recovered Kokuno hidden data and not a paper-exact field.
"""
from __future__ import annotations

from functools import lru_cache
import math
from typing import Any

import numpy as np

from .kokuno_public_candidate_velocity import _compact_cos8_envelope, _pulse_x
from .kokuno_public_carrier_resolved_velocity import (
    KokunoCarrierResolvedCandidateOscillatoryVelocity,
)

SCHEMA = "kokuno-public-z-pullback-oscillatory-velocity-v1"
PARENT_AGENT2_HEAD = "92852046e6e1ec0a5889b53f989df3ea5d79cb3b"


class KokunoPublicZPullbackCandidateOscillatoryVelocity(
    KokunoCarrierResolvedCandidateOscillatoryVelocity
):
    """Carrier-resolved candidate with source ``D_z`` pulled back to public z."""

    def _source_Z_by_beta(self, public_z: Any) -> np.ndarray:
        z = np.asarray(public_z, dtype=float)
        if z.ndim != 1 or not np.all(np.isfinite(z)):
            raise ValueError("public_z must be one finite sample vector")
        epsilon = np.asarray(self._static["geometry"]["epsilon_by_beta"], dtype=float)
        if epsilon.shape != (len(self.beta_labels),) or np.any(epsilon <= 0.0):
            raise RuntimeError("autonomous geometry lost one positive epsilon per beta")
        return z[:, None] * epsilon[None, :]

    def _inside_family(
        self,
        R: np.ndarray,
        theta: np.ndarray,
        Z: np.ndarray,
        t: np.ndarray,
    ) -> dict[str, Any]:
        """Evaluate interior points after the bandwise source/public z pullback.

        ``Z`` on entry is the public coordinate z.  The source phase sees
        ``Z_beta=epsilon_beta*z``.  Axial support and its derivative remain
        functions of public z, so the coefficient derivative supplied to the
        normalized source curl is exactly the ordinary public-z derivative.
        """
        labels = self.beta_labels
        n_beta = len(labels)
        geometry = self._static["geometry"]
        mass = self._static["mass"]
        background = self._background()
        m = R.size
        R_beta = np.broadcast_to(R[:, None], (m, n_beta))
        source_Z_beta = self._source_Z_by_beta(Z)
        L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
        epsilon = np.asarray(geometry["epsilon_by_beta"], dtype=float)
        v = np.stack((0.25 * L_s, 0.75 * L_s), axis=-1)
        frame_inputs = {
            "R": R_beta,
            "Z": source_Z_beta,
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
        # Because source_Z_beta=epsilon_beta*z and source D_z=epsilon_beta*d/dZ,
        # source D_z of a public-z envelope equals its ordinary d/dz derivative.
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
            - cross
            * D_r_norm_sq[..., None]
            / (k[..., None, None] * n_norm_sq[..., None] ** 2)
        )

        t_plus = envelope[:, None, None, None] * base_t
        D_r_C = (
            D_r_envelope[:, None, None, None] * C_base
            + envelope[:, None, None, None] * D_r_C_base
        )
        # The frozen background makes unlocalized C_base independent of source Z.
        # After the pullback D_z(source)=partial_z(public), so this is the complete
        # public-z coefficient derivative for the autonomous candidate.
        D_z_C = D_z_envelope[:, None, None, None] * C_base

        eta0 = np.asarray(geometry["eta"], dtype=float)
        D_r_eta0 = np.asarray(geometry["D_r_eta"], dtype=float)
        D_z_eta0 = np.asarray(geometry["D_z_eta"], dtype=float)
        eta = np.broadcast_to(eta0, (m, n_beta))
        D_r_eta = np.broadcast_to(D_r_eta0, (m, n_beta))
        D_z_eta = np.broadcast_to(D_z_eta0, (m, n_beta))

        u_star = np.asarray(background["u_star"], dtype=float)
        A_c = -np.asarray(frame["c0_by_beta"], dtype=float) * np.sqrt(
            1.0 + u_star[None, :] ** 2
        )
        tau = (t - self.time_min) / (self.time_max - self.time_min)
        normal_ratio = self.normal_target_ratio * (
            1.0 + self.time_modulation * np.sin(2.0 * math.pi * tau)
        )
        cross_ratio = self.cross_target_ratio * (
            1.0 + self.time_modulation * np.cos(2.0 * math.pi * tau)
        )
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
            "candidate_source_Z_by_beta": source_Z_beta,
            "candidate_epsilon_by_beta": epsilon,
            "public_source_D_z_pullback_applied": True,
        }

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
        out = super().evaluate(x, y, z, t)
        out["coordinate_contract"] = (
            "public R=hypot(x,y), theta=atan2(y,x), z; per beta source "
            "Z_beta=epsilon_beta*z, so public partial_z equals source "
            "D_z=epsilon_beta*partial_Z_beta"
        )
        out["public_source_D_z_pullback_applied"] = True
        return out

    def carrier_resolution_receipt(self) -> dict[str, Any]:
        """Record carrier scales in the pulled-back public coordinate."""
        labels = self.beta_labels
        n_beta = len(labels)
        geometry = self._static["geometry"]
        L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
        epsilon = np.asarray(geometry["epsilon_by_beta"], dtype=float)
        radius = 0.42
        background = self._background()
        v = np.stack((0.25 * L_s, 0.75 * L_s), axis=-1)
        frame = self.bridge.phase_family.phase_frame(
            R=np.full(n_beta, radius),
            Z=np.zeros(n_beta),
            theta=np.asarray(0.0),
            v=v,
            beta_labels=labels,
            L_s=L_s,
            **background,
        )
        k = np.asarray(frame["k_by_beta"], dtype=float)
        kp = np.asarray(frame["kp_by_beta_sign"], dtype=float)
        pz = np.asarray(frame["p_z_by_beta_sign"], dtype=float)
        x0 = np.asarray(frame["x0_by_beta_sign"], dtype=float)
        radial = k[:, None] * x0
        tangential = kp / radius
        # Source phase has k*p_z*Z/epsilon; Z=epsilon*z cancels epsilon.
        axial = k[:, None] * pz
        norm = np.sqrt(radial * radial + tangential * tangential + axial * axial)
        steps = (0.02, 0.01, 0.005)
        increments = {
            f"{step:.6g}": {
                "max_local_phase_increment": float(step * np.max(norm)),
                "max_tangential_phase_increment": float(step * np.max(np.abs(tangential))),
                "max_axial_phase_increment": float(step * np.max(np.abs(axial))),
                "max_radial_phase_increment": float(step * np.max(np.abs(radial))),
            }
            for step in steps
        }
        return {
            "schema": "kokuno-public-z-pullback-carrier-diagnostic-v1",
            "background_binding": "inherits frozen bounded carrier-resolved background from Agent-2 #551",
            "coordinate_pullback": "Z_beta=epsilon_beta*z; partial_z=D_z(source)",
            "diagnostic": {
                "reference_radius": radius,
                "max_abs_kp": int(np.max(np.abs(kp))),
                "max_local_wavenumber": float(np.max(norm)),
                "minimum_local_wavelength": float(2.0 * math.pi / np.max(norm)),
                "phase_increments_by_step": increments,
                "epsilon_by_beta": epsilon.tolist(),
            },
            "source_formulas_changed": False,
            "autonomous_background_changed_from_parent": False,
            "vector_potential_complete_curl_path_changed": False,
            "public_coordinate_pullback_changed": True,
            "independent_agent4_reaudit_required": True,
            "pde_validated": False,
            "paper_exact": False,
        }

    def to_payload(self) -> dict[str, Any]:
        payload = super().to_payload()
        payload["schema"] = SCHEMA
        payload["parent_agent2_head"] = PARENT_AGENT2_HEAD
        frozen = dict(payload["frozen_realization"])
        frozen["public_z_pullback"] = self.carrier_resolution_receipt()
        payload["frozen_realization"] = frozen
        autonomous = dict(payload["autonomous_choices"])
        autonomous["public_coordinate_pullback"] = (
            "per-beta Z_beta=epsilon_beta*z chosen so the source normalized "
            "D_z becomes ordinary public partial_z; not recovered hidden source data"
        )
        payload["autonomous_choices"] = autonomous
        truth = dict(payload["truth_boundary"])
        truth.update(
            {
                "public_source_D_z_pullback_applied": True,
                "public_coordinate_curl_consistency_targeted": True,
                "independent_fd4_divergence_reaudit_required": True,
                "independent_fd4_divergence_passed": False,
                "formal_full_domain_pde_gate_assessed": False,
                "pde_validated": False,
                "paper_exact": False,
            }
        )
        payload["truth_boundary"] = truth
        return payload


@lru_cache(maxsize=1)
def default_field() -> KokunoPublicZPullbackCandidateOscillatoryVelocity:
    return KokunoPublicZPullbackCandidateOscillatoryVelocity()


def velocity_osc(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return the public-z-pullback oscillatory candidate as ``[...,3]``."""
    return default_field().velocity_osc(x, y, z, t)
