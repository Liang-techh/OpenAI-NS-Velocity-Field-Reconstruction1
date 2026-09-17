"""Compact radial-stress reconstruction driven by a real phase-mean defect.

This module implements the finite-stage radial inverse used by the corrected
KokunoYumeto reconstruction, while keeping a strict boundary between public
source structure and repository-local engineering choices.

For e in {1, 2}, the source defines

    M_e F = integral r**e F(r) dr,
    P_e F = F - b_e M_e F,
    T_e F = -r**(-e) integral_0**r s**e P_e F(s) ds,

with integral r**e b_e dr = 1.  Consequently

    (d/dr + e/r) T_e F = -F + b_e M_e F,

and the equal-moment subtraction is what makes the stress compactly supported.

The current Agent-2 oscillatory field is an autonomous Cartesian box-supported
surrogate, not the source annular hierarchy.  For executable repository use we
therefore first take a cylindrical ring average of Agent-3's *measured* public
phase-mean NS defect, then apply a fixed C4 annular gate before the source radial
inverse.  This gate is an explicit autonomous adapter; it is not a hidden fit and
is not claimed to be paper-exact.  No free forcing is introduced and no residual
reduction is claimed by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.interpolate import CubicSpline

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_mean_defect import PhaseMeanDefectContract, cylindrical_theta_component


SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_FORMULA = "R6-R9 / C12 compact radial moment-complement inverse"


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _compact_c4_annulus_window(radii: np.ndarray, r_inner: float, r_outer: float) -> np.ndarray:
    """Fixed C4 scalar gate ``(1-s^2)^5_+`` on one radial annulus."""
    radii = np.asarray(radii, dtype=float)
    if not np.isfinite([r_inner, r_outer]).all() or not 0.0 < r_inner < r_outer:
        raise ValueError("annulus bounds must satisfy 0 < r_inner < r_outer")
    s = 2.0 * (radii - r_inner) / (r_outer - r_inner) - 1.0
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    return base**5


@dataclass(frozen=True)
class RadialMeanDefectProfile:
    """One cylindrical mean-defect channel sampled on a fixed annulus.

    ``raw_values`` are obtained directly from the real phase-mean defect after
    a uniform physical-angle ring average.  ``values`` are those same values
    multiplied by the fixed annular adapter window.  The distinction is kept so
    downstream reports cannot relabel the gated source as the whole mean defect.
    """

    radii: np.ndarray
    raw_values: np.ndarray
    values: np.ndarray
    exponent: int
    component: str
    time: float
    z: float
    angular_count: int
    annulus: tuple[float, float]
    source_kind: str = "real_phase_mean_defect_ring_average"

    def __post_init__(self) -> None:
        radii = np.asarray(self.radii, dtype=float)
        raw = np.asarray(self.raw_values, dtype=float)
        values = np.asarray(self.values, dtype=float)
        if radii.ndim != 1 or len(radii) < 9 or raw.shape != radii.shape or values.shape != radii.shape:
            raise ValueError("radii/raw_values/values must be matching 1-D arrays with at least 9 nodes")
        if not np.isfinite(radii).all() or not np.isfinite(raw).all() or not np.isfinite(values).all():
            raise ValueError("profile arrays must be finite")
        if np.any(np.diff(radii) <= 0.0) or radii[0] <= 0.0:
            raise ValueError("radii must be strictly increasing and bounded away from zero")
        if self.exponent not in (1, 2):
            raise ValueError("source radial weights are e=1 (axial) or e=2 (angular)")
        if self.component not in ("theta", "z"):
            raise ValueError("component must be 'theta' or 'z'")
        if (self.exponent, self.component) not in ((2, "theta"), (1, "z")):
            raise ValueError("source pairing is e=2/theta and e=1/z")
        if self.angular_count < 8 or self.angular_count % 2:
            raise ValueError("angular_count must be an even integer >= 8")
        if self.source_kind != "real_phase_mean_defect_ring_average":
            raise ValueError("scientific profile must identify the real phase-mean defect source")
        if not np.isclose(radii[0], self.annulus[0]) or not np.isclose(radii[-1], self.annulus[1]):
            raise ValueError("radii must include both declared annulus endpoints")
        # The adapter gate is exactly zero at both boundaries.  Fail closed if a
        # caller tries to pass an ungated source while claiming compact support.
        scale = max(1.0, float(np.max(np.abs(values))))
        if abs(values[0]) > 1.0e-13 * scale or abs(values[-1]) > 1.0e-13 * scale:
            raise ValueError("gated radial source must vanish at both annulus endpoints")

        object.__setattr__(self, "radii", radii.copy())
        object.__setattr__(self, "raw_values", raw.copy())
        object.__setattr__(self, "values", values.copy())
        self.radii.setflags(write=False)
        self.raw_values.setflags(write=False)
        self.values.setflags(write=False)

    @property
    def gate_capture_rms_ratio(self) -> float:
        return _rms(self.values) / max(_rms(self.raw_values), np.finfo(float).tiny)


@dataclass(frozen=True)
class CompactRadialStressInverse:
    """Executable moment-complement radial inverse for one sampled channel.

    A cubic interpolant is used only as a numerical representation of the fixed
    sampled radial source.  The scalar moment and the normalizing moment of the
    autonomous bump are integrated from those interpolants.  The stress itself
    is then evaluated from their antiderivatives, preserving the equal-moment
    cancellation at the outer edge to floating precision.
    """

    profile: RadialMeanDefectProfile
    bump_power: int = 5

    def __post_init__(self) -> None:
        if self.bump_power < 3:
            raise ValueError("bump_power must be >= 3")

    @property
    def r_inner(self) -> float:
        return float(self.profile.radii[0])

    @property
    def r_outer(self) -> float:
        return float(self.profile.radii[-1])

    def _splines(self) -> tuple[CubicSpline, CubicSpline, float, float]:
        r = self.profile.radii
        e = self.profile.exponent
        weighted_source = r**e * self.profile.values
        bump = _compact_c4_annulus_window(r, self.r_inner, self.r_outer) ** (self.bump_power / 5.0)
        weighted_bump = r**e * bump
        source_spline = CubicSpline(r, weighted_source)
        bump_spline = CubicSpline(r, weighted_bump)
        moment = float(source_spline.integrate(self.r_inner, self.r_outer))
        bump_moment = float(bump_spline.integrate(self.r_inner, self.r_outer))
        if not np.isfinite([moment, bump_moment]).all() or bump_moment <= 0.0:
            raise ValueError("radial source/bump moments must be finite and bump moment positive")
        return source_spline, bump_spline, moment, bump_moment

    @property
    def moment(self) -> float:
        return self._splines()[2]

    @property
    def bump_moment_before_normalization(self) -> float:
        return self._splines()[3]

    def source(self, radius) -> np.ndarray:
        """Interpolated gated ``F_e`` with exact zero extension."""
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = (q > self.r_inner) & (q < self.r_outer)
        if np.any(inside):
            source_spline, _, _, _ = self._splines()
            qq = q[inside]
            out[inside] = source_spline(qq) / qq**self.profile.exponent
        return out

    def normalized_bump(self, radius) -> np.ndarray:
        """Interpolated ``b_e`` normalized so ``integral r^e b_e dr = 1``."""
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = (q > self.r_inner) & (q < self.r_outer)
        if np.any(inside):
            _, bump_spline, _, bump_moment = self._splines()
            qq = q[inside]
            out[inside] = bump_spline(qq) / (bump_moment * qq**self.profile.exponent)
        return out

    def moment_complement(self, radius) -> np.ndarray:
        q = np.asarray(radius, dtype=float)
        return self.source(q) - self.normalized_bump(q) * self.moment

    def stress(self, radius) -> np.ndarray:
        """Evaluate ``sigma_e=T_e F`` with exact zero extension."""
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = (q > self.r_inner) & (q < self.r_outer)
        if not np.any(inside):
            return out

        source_spline, bump_spline, moment, bump_moment = self._splines()
        qq = q[inside]
        integrals = np.array(
            [
                source_spline.integrate(self.r_inner, float(value))
                - (moment / bump_moment) * bump_spline.integrate(self.r_inner, float(value))
                for value in qq
            ],
            dtype=float,
        )
        out[inside] = -integrals / qq**self.profile.exponent
        return out

    def metrics(self, *, query_count: int = 257, fd_step: float = 1.0e-5) -> dict[str, float]:
        if query_count < 33:
            raise ValueError("query_count must be >= 33")
        width = self.r_outer - self.r_inner
        if not 0.0 < fd_step < width / 100.0:
            raise ValueError("fd_step must be positive and small relative to the annulus")

        # Stay several FD steps away from the zero-extension edge so the identity
        # audit probes the reconstructed interior rather than a piecewise stencil.
        margin = max(10.0 * fd_step, 0.002 * width)
        query = np.linspace(self.r_inner + margin, self.r_outer - margin, query_count)
        sigma = self.stress(query)
        derivative = (self.stress(query + fd_step) - self.stress(query - fd_step)) / (2.0 * fd_step)
        lhs = derivative + self.profile.exponent * sigma / query
        rhs = -self.source(query) + self.normalized_bump(query) * self.moment
        error = lhs - rhs

        source_spline, bump_spline, moment, bump_moment = self._splines()
        complement_moment = float(
            source_spline.integrate(self.r_inner, self.r_outer)
            - (moment / bump_moment) * bump_spline.integrate(self.r_inner, self.r_outer)
        )
        raw_tail_without_subtraction = abs(moment) / self.r_outer**self.profile.exponent
        return {
            "source_rms": _rms(self.source(query)),
            "source_max_abs": float(np.max(np.abs(self.source(query)))),
            "source_weighted_moment": moment,
            "normalized_bump_weighted_moment": 1.0,
            "moment_complement_weighted_moment_abs": abs(complement_moment),
            "stress_rms": _rms(sigma),
            "stress_max_abs": float(np.max(np.abs(sigma))),
            "inner_edge_stress_abs": abs(float(self.stress(np.array([self.r_inner]))[0])),
            "outer_edge_stress_abs": abs(float(self.stress(np.array([self.r_outer]))[0])),
            "outside_stress_max_abs": float(
                np.max(
                    np.abs(
                        self.stress(
                            np.array(
                                [self.r_inner - 0.1 * width, self.r_outer + 0.1 * width],
                                dtype=float,
                            )
                        )
                    )
                )
            ),
            "raw_inverse_outer_tail_without_moment_subtraction": raw_tail_without_subtraction,
            "radial_identity_fd_rms": _rms(error),
            "radial_identity_fd_max_abs": float(np.max(np.abs(error))),
            "gate_capture_rms_ratio": self.profile.gate_capture_rms_ratio,
        }


def sample_real_phase_mean_radial_profiles(
    base_velocity: Callable[[np.ndarray, float], np.ndarray],
    correction: KokunoCompleteCurlCorrection,
    contract: PhaseMeanDefectContract,
    *,
    time: float = 0.5,
    z: float = 0.0,
    r_inner: float = 0.12,
    r_outer: float = 0.72,
    radial_count: int = 65,
    angular_count: int = 16,
) -> tuple[RadialMeanDefectProfile, RadialMeanDefectProfile, dict[str, float]]:
    """Project the *real* phase-mean defect onto radial theta/z channels.

    The ring average and C4 annular gate are repository-local adapters.  They are
    fixed independently of the measured defect and do not alter the upstream
    candidate/correction.  The returned profiles are the only scientific input
    accepted by the current report path; no hand-entered residual array is used.
    """
    if not callable(base_velocity):
        raise TypeError("base_velocity must be callable")
    if not isinstance(correction, KokunoCompleteCurlCorrection):
        raise TypeError("correction must be KokunoCompleteCurlCorrection")
    if not isinstance(contract, PhaseMeanDefectContract):
        raise TypeError("contract must be PhaseMeanDefectContract")
    center = np.asarray(correction.center, dtype=float)
    if np.linalg.norm(center[:2]) > 1.0e-12:
        raise ValueError("current cylindrical adapter requires Agent-2 correction centered on the physical axis")
    if radial_count < 17 or angular_count < 8 or angular_count % 2:
        raise ValueError("radial_count >= 17 and even angular_count >= 8 are required")
    if not 0.0 < r_inner < r_outer:
        raise ValueError("invalid radial annulus")
    if r_outer >= min(correction.half_widths[0], correction.half_widths[1]):
        raise ValueError("adapter annulus must stay inside the current correction's transverse box widths")
    if not abs(z - center[2]) < correction.half_widths[2]:
        raise ValueError("z slice must lie inside the current correction support")

    radii = np.linspace(r_inner, r_outer, radial_count)
    angles = 2.0 * np.pi * np.arange(angular_count, dtype=float) / angular_count
    rr, aa = np.meshgrid(radii, angles, indexing="ij")
    points = np.stack((rr * np.cos(aa), rr * np.sin(aa), np.full_like(rr, z)), axis=-1)
    flat_points = points.reshape(-1, 3)

    sample = contract.evaluate(base_velocity, correction, flat_points, float(time))
    increment = sample.mean_defect_increment
    theta_values = cylindrical_theta_component(increment, flat_points).reshape(radial_count, angular_count)
    z_values = increment[:, 2].reshape(radial_count, angular_count)
    raw_theta = np.mean(theta_values, axis=1)
    raw_z = np.mean(z_values, axis=1)

    window = _compact_c4_annulus_window(radii, r_inner, r_outer)
    theta_profile = RadialMeanDefectProfile(
        radii=radii,
        raw_values=raw_theta,
        values=raw_theta * window,
        exponent=2,
        component="theta",
        time=float(time),
        z=float(z),
        angular_count=angular_count,
        annulus=(float(r_inner), float(r_outer)),
    )
    axial_profile = RadialMeanDefectProfile(
        radii=radii,
        raw_values=raw_z,
        values=raw_z * window,
        exponent=1,
        component="z",
        time=float(time),
        z=float(z),
        angular_count=angular_count,
        annulus=(float(r_inner), float(r_outer)),
    )
    projection_receipt = {
        "full_mean_defect_increment_rms": _rms(increment),
        "ring_theta_raw_rms": _rms(raw_theta),
        "ring_z_raw_rms": _rms(raw_z),
        "theta_gate_capture_rms_ratio": theta_profile.gate_capture_rms_ratio,
        "z_gate_capture_rms_ratio": axial_profile.gate_capture_rms_ratio,
        "phase_count": float(contract.phase_count),
        "angular_count": float(angular_count),
    }
    return theta_profile, axial_profile, projection_receipt


def generate_current_candidate_radial_stress_report(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    output: str | Path = "artifacts/kokuno_agent3/radial_stress_report.json",
    time: float = 0.5,
    radial_count: int = 65,
    angular_count: int = 16,
) -> dict:
    """Apply the radial inverse to measured mean-defect channels of the real artifact."""
    from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate

    candidate_path = Path(candidate_path)
    candidate = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    correction = KokunoCompleteCurlCorrection()
    contract = PhaseMeanDefectContract(phase_count=8)
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
        candidate.at_points,
        correction,
        contract,
        time=time,
        z=0.0,
        r_inner=0.12,
        r_outer=0.72,
        radial_count=radial_count,
        angular_count=angular_count,
    )

    theta_inverse = CompactRadialStressInverse(theta_profile)
    axial_inverse = CompactRadialStressInverse(axial_profile)
    report = {
        "task": "KOKUNO-A3-RADIAL-STRESS-INVERSE-002",
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "formula": SOURCE_FORMULA,
            "source_pairing": {"theta": "e=2", "z": "e=1"},
            "source_requirement": "equal weighted moment subtraction is required for compact radial stress support",
        },
        "candidate_path": str(candidate_path),
        "candidate_sha256": candidate.sha256,
        "correction": correction.metadata(),
        "real_defect_contract": {
            "nu": contract.nu,
            "derivative_step": contract.step,
            "phase_count": contract.phase_count,
            "time": float(time),
            "z": 0.0,
        },
        "engineering_adapter": {
            "radial_annulus": [0.12, 0.72],
            "radial_count": radial_count,
            "angular_count": angular_count,
            "ring_average": "uniform physical cylindrical angle",
            "annular_gate": "fixed C4 (1-s^2)^5 window chosen before defect metrics",
            "moving_q_scale_used": False,
            "paper_exact": False,
        },
        "projection": projection,
        "theta_e2": theta_inverse.metrics(),
        "axial_e1": axial_inverse.metrics(),
        "truth_boundary": {
            "consumes_real_candidate_artifact": True,
            "consumes_real_phase_mean_defect": True,
            "surrogate_defect_used": False,
            "radial_stress_inverse_applied_to_measured_defect_projection": True,
            "mean_velocity_correction_applied": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Reconstruct compact radial stresses from the real phase-mean defect")
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--output", default="artifacts/kokuno_agent3/radial_stress_report.json")
    parser.add_argument("--time", type=float, default=0.5)
    parser.add_argument("--radial-count", type=int, default=65)
    parser.add_argument("--angular-count", type=int, default=16)
    args = parser.parse_args()
    report = generate_current_candidate_radial_stress_report(
        candidate_path=args.candidate,
        output=args.output,
        time=args.time,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
    )
    print(json.dumps({"candidate_sha256": report["candidate_sha256"], "projection": report["projection"], "theta_e2": report["theta_e2"], "axial_e1": report["axial_e1"]}, indent=2))


if __name__ == "__main__":
    main()
