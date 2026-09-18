"""Fail-closed gain bookkeeping for Kokuno Agent-3 finite correction cycles.

The corrected 2026-09-09 Kokuno reader proves *exponent-class* gains for its
source four-step cycle.  Those gains are not raw residual contraction ratios.
This module makes that distinction executable and then replays Agent 3's most
recent real one-column finite cycle through a small routing guard.

The guard deliberately does not construct an oscillatory column, fit a force or
pressure, alter the registered PDE thresholds, or reinterpret a local
phase-mean defect as the repository's full-domain normalized NS residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_signed_amplitude_curl_cycle import (
    generate_signed_amplitude_curl_cycle_report,
)


TASK = "KOKUNO-A3-CORRECTION-CYCLE-GAIN-GUARD-009"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Four-step correction-cycle gain and invariant audit"
SOURCE_KAPPA_S = 1.0e-5
SOURCE_REQUIRED_EXPONENT_GAIN = 0.1
REGISTERED_DIVERGENCE_MAX = 1.0e-5

# Repository-wide comparison reference surfaced by coordination PR #290.  These
# numbers are only used if a caller explicitly declares a directly comparable
# full-domain residual receipt.  The local Agent-3 phase-mean report generated
# by this module is deliberately marked non-comparable.
ST006_SHA256 = "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3"
ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_VOLUME_L2 = 0.10758432876230622
ST006_VALIDATION_SEED = 9_172_801
ST006_VALIDATION_POINTS = 4096
ST006_FINEST_SPATIAL_STEP = 0.005


def source_exponent_gain_margins(kappa_s: float = SOURCE_KAPPA_S) -> dict[str, Any]:
    """Return the source's displayed exponent gains without calling them residual ratios."""
    kappa_s = float(kappa_s)
    if not np.isfinite(kappa_s) or not 0.0 < kappa_s < 0.01:
        raise ValueError("kappa_s must be finite and lie in (0,0.01)")
    wave = min(0.4, 0.4 - kappa_s, 0.5 - 2.0 * kappa_s, 0.5 - 2.0 * kappa_s)
    tangential_mean = min(0.17, 1.0 - 4.0 * kappa_s, 1.0 - 4.0 * kappa_s)
    defect = 0.9 - 4.0 * kappa_s
    gains = {
        "wave_residual_exponent_gain": float(wave),
        "full_tangential_mean_exponent_gain": float(tangential_mean),
        "defect_exponent_gain": float(defect),
    }
    return {
        "kappa_s": kappa_s,
        "required_exponent_gain": SOURCE_REQUIRED_EXPONENT_GAIN,
        **gains,
        "all_source_displayed_gains_exceed_required": bool(
            min(gains.values()) > SOURCE_REQUIRED_EXPONENT_GAIN
        ),
        "interpretation": (
            "These are source asymptotic exponent-class gains after the source four-step cycle, "
            "not measured residual_before/residual_after contraction ratios."
        ),
    }


def empirical_exponent_gain(before: float, after: float, q_scale: float) -> float:
    """Convert a measured contraction to an exponent gain only when a Q scale is meaningful."""
    before = float(before)
    after = float(after)
    q_scale = float(q_scale)
    if not np.isfinite([before, after, q_scale]).all():
        raise ValueError("before, after and q_scale must be finite")
    if before <= 0.0 or after <= 0.0:
        raise ValueError("before and after must be positive")
    if q_scale <= 1.0:
        raise ValueError("q_scale must be greater than one")
    return float(-np.log(after / before) / np.log(q_scale))


@dataclass(frozen=True)
class ComparableFullDomainResidual:
    """Explicit full-domain receipt used only for an ST006-like comparison."""

    momentum_max: float
    volume_l2: float
    seed: int
    point_count: int
    finest_spatial_step: float

    def __post_init__(self) -> None:
        numeric = np.asarray(
            [self.momentum_max, self.volume_l2, self.finest_spatial_step], dtype=float
        )
        if not np.isfinite(numeric).all() or np.any(numeric <= 0.0):
            raise ValueError("full-domain residual metrics and step must be positive and finite")
        if int(self.seed) < 0 or int(self.point_count) <= 0:
            raise ValueError("seed must be nonnegative and point_count positive")


def _require_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def evaluate_correction_cycle_gain_guard(
    cycle_report: Mapping[str, Any],
    *,
    source_class_comparable: bool = False,
    source_q_scale: float | None = None,
    full_domain_receipt: ComparableFullDomainResidual | None = None,
) -> dict[str, Any]:
    """Reject divergent cycles while keeping source and empirical gain notions separate."""
    report = _require_mapping(cycle_report, "cycle_report")
    inputs = _require_mapping(report.get("inputs"), "cycle_report.inputs")
    held_in = _require_mapping(report.get("held_in"), "cycle_report.held_in")
    held_out = _require_mapping(report.get("held_out"), "cycle_report.held_out")
    divergence = _require_mapping(
        report.get("independent_signed_correction_divergence"),
        "cycle_report.independent_signed_correction_divergence",
    )
    truth = _require_mapping(report.get("truth_boundary"), "cycle_report.truth_boundary")

    if bool(inputs.get("surrogate_defect_used", True)):
        raise ValueError("cycle guard refuses a surrogate defect receipt")
    if not bool(truth.get("real_candidate_defect_consumed", False)):
        raise ValueError("cycle guard requires a real candidate defect receipt")
    if not bool(truth.get("finite_correction_cycle_run", False)):
        raise ValueError("cycle guard requires an actually executed finite correction cycle")

    held_in_ratio = float(held_in["rms_ratio"])
    held_out_ratio = float(held_out["aggregate_mean_defect_rms_ratio"])
    theta_ratio = float(held_out["aggregate_theta_rms_ratio"])
    pressure_ratio = float(held_out["pressure_inclusive_raw_phase_mean_operator_rms_ratio"])
    correction_rms = float(held_out["signed_correction_rms"])
    divergence_max = float(divergence["max"])
    values = np.asarray(
        [held_in_ratio, held_out_ratio, theta_ratio, pressure_ratio, correction_rms, divergence_max],
        dtype=float,
    )
    if not np.isfinite(values).all() or np.any(values < 0.0):
        raise ValueError("cycle metrics must be finite and nonnegative")

    checks = {
        "held_in_total_mean_defect_improves": bool(held_in_ratio < 1.0),
        "held_out_total_mean_defect_improves": bool(held_out_ratio < 1.0),
        "held_out_theta_mean_defect_nonworsening": bool(theta_ratio <= 1.0),
        "pressure_inclusive_raw_operator_nonworsening_0p1pct": bool(pressure_ratio <= 1.001),
        "independent_correction_divergence_max_le_1e5": bool(
            divergence_max <= REGISTERED_DIVERGENCE_MAX
        ),
        "correction_nontrivial": bool(correction_rms > 1.0e-12),
    }
    accepted = bool(all(checks.values()))

    source_gains = source_exponent_gain_margins()
    if source_class_comparable:
        if source_q_scale is None:
            raise ValueError("source_q_scale is required when source_class_comparable=true")
        held_out_empirical_gain = empirical_exponent_gain(1.0, held_out_ratio, source_q_scale)
        theta_empirical_gain = empirical_exponent_gain(1.0, theta_ratio, source_q_scale)
        empirical_source_comparison: dict[str, Any] | None = {
            "q_scale": float(source_q_scale),
            "held_out_total_empirical_exponent_gain": held_out_empirical_gain,
            "held_out_theta_empirical_exponent_gain": theta_empirical_gain,
            "held_out_total_meets_source_required_gain": bool(
                held_out_empirical_gain >= SOURCE_REQUIRED_EXPONENT_GAIN
            ),
            "held_out_theta_meets_source_required_gain": bool(
                theta_empirical_gain >= SOURCE_REQUIRED_EXPONENT_GAIN
            ),
            "warning": (
                "Even a positive empirical exponent gain does not replace the repository full-domain PDE gate."
            ),
        }
    else:
        if source_q_scale is not None:
            raise ValueError("source_q_scale must be omitted when source_class_comparable=false")
        empirical_source_comparison = None

    st006_reference = {
        "candidate_sha256": ST006_SHA256,
        "momentum_max": ST006_MOMENTUM_MAX,
        "volume_l2": ST006_VOLUME_L2,
        "validation_seed": ST006_VALIDATION_SEED,
        "validation_points": ST006_VALIDATION_POINTS,
        "finest_spatial_step": ST006_FINEST_SPATIAL_STEP,
        "pde_validated": False,
    }
    if full_domain_receipt is None:
        st006_comparison: dict[str, Any] = {
            "directly_comparable": False,
            "repository_level_pde_improvement_claimed": False,
            "reason": (
                "Agent-3 cycle metrics are core-local phase-mean defect diagnostics, not the ST006 "
                "4096-point full-domain normalized momentum protocol."
            ),
        }
    else:
        st006_comparison = {
            "directly_comparable": True,
            "momentum_max_ratio_to_st006": float(
                full_domain_receipt.momentum_max / ST006_MOMENTUM_MAX
            ),
            "volume_l2_ratio_to_st006": float(
                full_domain_receipt.volume_l2 / ST006_VOLUME_L2
            ),
            "beats_st006_on_both_reported_momentum_metrics": bool(
                full_domain_receipt.momentum_max < ST006_MOMENTUM_MAX
                and full_domain_receipt.volume_l2 < ST006_VOLUME_L2
            ),
            "repository_level_pde_improvement_claimed": False,
            "reason": (
                "Beating ST006 would be a benchmark improvement only; the fixed 1e-3 PDE gate and "
                "independent acceptance protocol would still be required."
            ),
        }

    return {
        "source_exponent_gain_audit": source_gains,
        "source_class_comparable_to_this_cycle": bool(source_class_comparable),
        "empirical_source_exponent_comparison": empirical_source_comparison,
        "measured_cycle": {
            "held_in_total_mean_defect_ratio": held_in_ratio,
            "held_out_total_mean_defect_ratio": held_out_ratio,
            "held_out_theta_mean_defect_ratio": theta_ratio,
            "pressure_inclusive_raw_operator_ratio": pressure_ratio,
            "signed_correction_rms": correction_rms,
            "independent_correction_divergence_max": divergence_max,
        },
        "finite_cycle_checks": checks,
        "accepted_for_next_cycle": accepted,
        "st006_reference": st006_reference,
        "st006_comparison": st006_comparison,
    }


def generate_correction_cycle_gain_guard_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/correction_cycle_gain_guard_report.json",
    cycle_output: str | Path = "artifacts/kokuno_agent3/correction_cycle_gain_guard_cycle_replay.json",
) -> dict[str, Any]:
    """Replay the real one-column cycle and bind a fail-closed gain/routing receipt."""
    cycle_report = generate_signed_amplitude_curl_cycle_report(output=cycle_output)
    guard = evaluate_correction_cycle_gain_guard(
        cycle_report,
        source_class_comparable=False,
        full_domain_receipt=None,
    )
    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "section": SOURCE_SECTION,
            "scope": (
                "source exponent-class bookkeeping is recorded as a divergence guard only; "
                "it is not relabeled as a measured residual contraction"
            ),
        },
        "real_cycle_replay": {
            "task": cycle_report.get("task"),
            "inputs": cycle_report.get("inputs"),
            "held_in": cycle_report.get("held_in"),
            "held_out": cycle_report.get("held_out"),
            "independent_signed_correction_divergence": cycle_report.get(
                "independent_signed_correction_divergence"
            ),
            "finite_cycle_decision": cycle_report.get("finite_cycle_decision"),
            "truth_boundary": cycle_report.get("truth_boundary"),
        },
        "gain_guard": guard,
        "routing": {
            "rerun_rejected_one_column_cycle": False,
            "widen_one_column_amplitude_or_fit_degree": False,
            "wait_for_genuinely_independent_second_public_covariance_column": True,
            "second_column_must_first_pass_missing_target_rank_coverage_contract": True,
            "then_rerun_real_held_in_held_out_cycle": True,
            "reason": (
                "The existing public one-column correction worsens the real held-in, held-out, and theta "
                "mean-defect metrics. Positive source exponent margins do not override measured divergence."
            ),
        },
        "truth_boundary": {
            "surrogate_defect_used": False,
            "real_candidate_defect_consumed": True,
            "finite_correction_cycle_replayed": True,
            "source_exponent_gain_replayed": True,
            "source_exponent_gain_used_as_raw_residual_ratio": False,
            "source_class_comparability_established": False,
            "st006_direct_comparability_established": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/correction_cycle_gain_guard_report.json",
    )
    parser.add_argument(
        "--cycle-output",
        default="artifacts/kokuno_agent3/correction_cycle_gain_guard_cycle_replay.json",
    )
    args = parser.parse_args()
    report = generate_correction_cycle_gain_guard_report(
        output=args.output,
        cycle_output=args.cycle_output,
    )
    guard = report["gain_guard"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "source_wave_exponent_gain": guard["source_exponent_gain_audit"][
                    "wave_residual_exponent_gain"
                ],
                "source_tangential_mean_exponent_gain": guard[
                    "source_exponent_gain_audit"
                ]["full_tangential_mean_exponent_gain"],
                "source_defect_exponent_gain": guard["source_exponent_gain_audit"][
                    "defect_exponent_gain"
                ],
                "held_out_total_mean_defect_ratio": guard["measured_cycle"][
                    "held_out_total_mean_defect_ratio"
                ],
                "held_out_theta_mean_defect_ratio": guard["measured_cycle"][
                    "held_out_theta_mean_defect_ratio"
                ],
                "accepted_for_next_cycle": guard["accepted_for_next_cycle"],
                "st006_directly_comparable": guard["st006_comparison"]["directly_comparable"],
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
