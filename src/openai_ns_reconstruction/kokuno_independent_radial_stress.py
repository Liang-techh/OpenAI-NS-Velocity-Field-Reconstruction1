"""Independent numerical preflight for Agent-3's Kokuno radial-stress inverse.

This module is validation-only.  It deliberately does not reuse Agent-3's
``CubicSpline`` antiderivative path to reconstruct the radial inverse.  Instead
it consumes only the sampled ``RadialMeanDefectProfile`` and the public
``stress(radius)`` evaluator, builds an independent PCHIP reconstruction of the
weighted source/bump, and differentiates the public stress with a fourth-order
finite-difference stencil on held-out off-grid radii.

The full registered Navier--Stokes gate remains unassessed here.  Agent 1 still
does not expose a complete 3-D Kokuno leading velocity, and Agent 3 has not yet
applied a finite mean-velocity correction cycle.  Therefore the fixed 1e-3
normalized full-momentum threshold is recorded but cannot be promoted by this
structural preflight.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    RadialMeanDefectProfile,
    sample_real_phase_mean_radial_profiles,
)


HELD_OUT_SEED = 9172761
RADIAL_COUNTS = (33, 49, 65)
ANGULAR_COUNT = 16
REGISTERED_FULL_NS_THRESHOLD = 1.0e-3
MUTATION_SCALE = 0.05


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _window(radius: np.ndarray, r_inner: float, r_outer: float, bump_power: int) -> np.ndarray:
    """Independent copy of the fixed compact scalar shape from its public formula."""
    radius = np.asarray(radius, dtype=float)
    s = 2.0 * (radius - r_inner) / (r_outer - r_inner) - 1.0
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    return base**bump_power


class IndependentPchipRadialInverse:
    """Independent PCHIP/quadrature reconstruction from profile samples only."""

    def __init__(self, profile: RadialMeanDefectProfile, *, bump_power: int = 5):
        if not isinstance(profile, RadialMeanDefectProfile):
            raise TypeError("profile must be RadialMeanDefectProfile")
        if bump_power < 3:
            raise ValueError("bump_power must be >= 3")
        self.profile = profile
        self.bump_power = int(bump_power)
        r = profile.radii
        e = profile.exponent
        weighted_source = r**e * profile.values
        bump = _window(r, float(r[0]), float(r[-1]), self.bump_power)
        weighted_bump = r**e * bump
        self._source = PchipInterpolator(r, weighted_source, extrapolate=False)
        self._bump = PchipInterpolator(r, weighted_bump, extrapolate=False)
        self._source_anti = self._source.antiderivative()
        self._bump_anti = self._bump.antiderivative()
        self.r_inner = float(r[0])
        self.r_outer = float(r[-1])
        self.moment = float(self._source_anti(self.r_outer) - self._source_anti(self.r_inner))
        self.bump_moment = float(self._bump_anti(self.r_outer) - self._bump_anti(self.r_inner))
        if not np.isfinite([self.moment, self.bump_moment]).all() or self.bump_moment <= 0.0:
            raise ValueError("independent radial moments must be finite and bump moment positive")

    def _inside(self, radius: np.ndarray) -> np.ndarray:
        return (radius > self.r_inner) & (radius < self.r_outer)

    def source(self, radius) -> np.ndarray:
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = self._inside(q)
        if np.any(inside):
            qq = q[inside]
            out[inside] = self._source(qq) / qq**self.profile.exponent
        return out

    def normalized_bump(self, radius) -> np.ndarray:
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = self._inside(q)
        if np.any(inside):
            qq = q[inside]
            out[inside] = self._bump(qq) / (self.bump_moment * qq**self.profile.exponent)
        return out

    def stress(self, radius) -> np.ndarray:
        q = np.asarray(radius, dtype=float)
        if not np.isfinite(q).all():
            raise ValueError("radius must be finite")
        out = np.zeros_like(q, dtype=float)
        inside = self._inside(q)
        if np.any(inside):
            qq = q[inside]
            source_int = self._source_anti(qq) - self._source_anti(self.r_inner)
            bump_int = self._bump_anti(qq) - self._bump_anti(self.r_inner)
            complement_int = source_int - (self.moment / self.bump_moment) * bump_int
            out[inside] = -complement_int / qq**self.profile.exponent
        return out


def _fd4_public_stress(reference: CompactRadialStressInverse, radii: np.ndarray, step: float) -> np.ndarray:
    radii = np.asarray(radii, dtype=float)
    return (
        reference.stress(radii - 2.0 * step)
        - 8.0 * reference.stress(radii - step)
        + 8.0 * reference.stress(radii + step)
        - reference.stress(radii + 2.0 * step)
    ) / (12.0 * step)


def _held_out_radii(profile: RadialMeanDefectProfile, *, seed: int, count: int = 257) -> np.ndarray:
    rng = np.random.default_rng(seed)
    width = float(profile.radii[-1] - profile.radii[0])
    margin = 0.025 * width
    values = rng.uniform(float(profile.radii[0] + margin), float(profile.radii[-1] - margin), count)
    values.sort()
    return values


def validate_profile_independently(
    profile: RadialMeanDefectProfile,
    *,
    seed: int = HELD_OUT_SEED,
    held_out_count: int = 257,
) -> dict[str, float]:
    """Cross-check one Agent-3 profile with independent interpolation/differentiation."""
    reference = CompactRadialStressInverse(profile)
    independent = IndependentPchipRadialInverse(profile)
    radii = _held_out_radii(profile, seed=seed, count=held_out_count)
    width = independent.r_outer - independent.r_inner
    # Keep the FD stencil well inside the exact-zero extension and much smaller
    # than the sample spacing.  This is deliberately unrelated to Agent-3's 1e-5
    # centered-difference audit.
    step = width / (64.0 * (len(profile.radii) - 1))

    public_stress = reference.stress(radii)
    independent_stress = independent.stress(radii)
    stress_error = public_stress - independent_stress
    stress_scale = max(_rms(public_stress), np.finfo(float).tiny)

    derivative = _fd4_public_stress(reference, radii, step)
    lhs = derivative + profile.exponent * public_stress / radii
    rhs = -independent.source(radii) + independent.normalized_bump(radii) * independent.moment
    identity_error = lhs - rhs

    mutated_stress = (1.0 + MUTATION_SCALE) * public_stress
    mutated_derivative = (1.0 + MUTATION_SCALE) * derivative
    mutated_lhs = mutated_derivative + profile.exponent * mutated_stress / radii
    mutation_error = mutated_lhs - rhs

    outside = np.array(
        [
            independent.r_inner - 0.05 * width,
            independent.r_inner,
            independent.r_outer,
            independent.r_outer + 0.05 * width,
        ],
        dtype=float,
    )
    outside_public = reference.stress(outside)
    outside_independent = independent.stress(outside)

    raw_tail = abs(independent.moment) / independent.r_outer**profile.exponent
    moment_relative_difference = abs(reference.moment - independent.moment) / max(abs(reference.moment), 1.0e-30)

    return {
        "radial_count": float(len(profile.radii)),
        "held_out_count": float(held_out_count),
        "fd4_step": float(step),
        "public_stress_rms": _rms(public_stress),
        "independent_stress_rms": _rms(independent_stress),
        "stress_disagreement_rms": _rms(stress_error),
        "stress_disagreement_max_abs": float(np.max(np.abs(stress_error))),
        "stress_disagreement_normalized_rms": _rms(stress_error) / stress_scale,
        "identity_rms": _rms(identity_error),
        "identity_max_abs": float(np.max(np.abs(identity_error))),
        "mutation_identity_rms": _rms(mutation_error),
        "mutation_to_baseline_identity_ratio": _rms(mutation_error) / max(_rms(identity_error), 1.0e-30),
        "reference_vs_independent_moment_relative_difference": float(moment_relative_difference),
        "raw_inverse_outer_tail_without_moment_subtraction": float(raw_tail),
        "public_support_edge_outside_max_abs": float(np.max(np.abs(outside_public))),
        "independent_support_edge_outside_max_abs": float(np.max(np.abs(outside_independent))),
        "gate_capture_rms_ratio": float(profile.gate_capture_rms_ratio),
    }


def _observed_order(coarse: float, fine: float, coarse_n: int, fine_n: int) -> float | None:
    if coarse <= 0.0 or fine <= 0.0:
        return None
    hc = 1.0 / (coarse_n - 1)
    hf = 1.0 / (fine_n - 1)
    return float(np.log(coarse / fine) / np.log(hc / hf))


def run_real_candidate_preflight(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    output: str | Path = "artifacts/kokuno_independent_validation/radial_stress_preflight.json",
    radial_counts: tuple[int, ...] = RADIAL_COUNTS,
    angular_count: int = ANGULAR_COUNT,
    time: float = 0.5,
) -> dict[str, object]:
    """Run three-resolution independent validation on the current real artifact."""
    from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate

    if len(radial_counts) < 3 or any(n < 17 for n in radial_counts):
        raise ValueError("at least three radial resolutions >= 17 are required")
    if tuple(sorted(radial_counts)) != tuple(radial_counts) or len(set(radial_counts)) != len(radial_counts):
        raise ValueError("radial_counts must be strictly increasing")

    candidate_path = Path(candidate_path)
    candidate = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    correction = KokunoCompleteCurlCorrection()
    contract = PhaseMeanDefectContract(phase_count=8)

    levels: list[dict[str, object]] = []
    for index, radial_count in enumerate(radial_counts):
        theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
            candidate.at_points,
            correction,
            contract,
            time=float(time),
            radial_count=int(radial_count),
            angular_count=int(angular_count),
        )
        theta = validate_profile_independently(theta_profile, seed=HELD_OUT_SEED + 17 * index)
        axial = validate_profile_independently(axial_profile, seed=HELD_OUT_SEED + 1000 + 17 * index)
        levels.append(
            {
                "radial_count": int(radial_count),
                "projection": projection,
                "theta_e2": theta,
                "axial_e1": axial,
            }
        )

    theta_disagreement = [float(level["theta_e2"]["stress_disagreement_normalized_rms"]) for level in levels]
    theta_identity = [float(level["theta_e2"]["identity_rms"]) for level in levels]
    theta_orders = [
        _observed_order(theta_disagreement[i], theta_disagreement[i + 1], radial_counts[i], radial_counts[i + 1])
        for i in range(len(radial_counts) - 1)
    ]
    identity_orders = [
        _observed_order(theta_identity[i], theta_identity[i + 1], radial_counts[i], radial_counts[i + 1])
        for i in range(len(radial_counts) - 1)
    ]

    finest = levels[-1]["theta_e2"]
    structural_preflight_passed = bool(
        finest["stress_disagreement_normalized_rms"] < 5.0e-3
        and finest["identity_rms"] < 5.0e-6
        and finest["public_support_edge_outside_max_abs"] == 0.0
        and finest["independent_support_edge_outside_max_abs"] == 0.0
        and finest["mutation_to_baseline_identity_ratio"] > 20.0
        and finest["raw_inverse_outer_tail_without_moment_subtraction"] > 1.0e-8
    )

    report: dict[str, object] = {
        "task_id": "KOKUNO-VAL-RADIAL-STRESS-PREFLIGHT-002",
        "agent3_exact_head": "cac3f4de3de2ca47364423e9a1b70407dd819526",
        "candidate_path": str(candidate_path),
        "candidate_sha256": candidate.sha256,
        "time": float(time),
        "held_out_seed": HELD_OUT_SEED,
        "radial_counts": list(radial_counts),
        "angular_count": int(angular_count),
        "operator": "independent PCHIP weighted antiderivative plus fourth-order finite difference on public stress",
        "levels": levels,
        "theta_normalized_stress_disagreement_orders": theta_orders,
        "theta_identity_orders": identity_orders,
        "structural_preflight_passed": structural_preflight_passed,
        "registered_full_ns_gate": {
            "threshold": REGISTERED_FULL_NS_THRESHOLD,
            "assessed": False,
            "reason": (
                "No complete 3-D Kokuno leading velocity and no applied finite mean-velocity correction cycle with "
                "compatible pressure/fixed-or-restricted forcing are present in this ancestry."
            ),
        },
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently validate the real Kokuno radial-stress primitive")
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--output", default="artifacts/kokuno_independent_validation/radial_stress_preflight.json")
    parser.add_argument("--time", type=float, default=0.5)
    parser.add_argument("--angular-count", type=int, default=16)
    args = parser.parse_args()
    report = run_real_candidate_preflight(
        candidate_path=args.candidate,
        output=args.output,
        time=args.time,
        angular_count=args.angular_count,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
