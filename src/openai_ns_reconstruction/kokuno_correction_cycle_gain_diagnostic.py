"""Source-structured correction-cycle gain gate for Kokuno Agent 3.

The corrected 2026-09-09 reconstruction fixes the stage arithmetic

    kappa_s = 1e-5,
    sigma_j = 1/5 + j/10,
    B_j = 1/2 + sigma_j,
    C_j = 1 + sigma_j,

and records, after the complete mean/radial/temporal/compatibility cycle,
source-side exponent gains including 0.39999 for the wave residual,
0.17 for the full tangential mean, and 0.89996 for compatibility defects.
Those source exponents are not numerical Navier--Stokes residual ratios.

This module freezes that public arithmetic and composes it with the existing
Agent-3 finite-cycle ledger.  The public admission path accepts raw states and
corrections, reruns the actual momentum defect through #645, and rejects any
empirically divergent intermediate correction before reporting the source
bookkeeping.  Callers cannot supply a residual, defect, source exponent,
gain threshold, B/C/sigma value, or kappa_s.

Passing this gate is therefore only a consistency diagnostic: it is neither
Kokuno's theorem, nor an independent source proof, nor evidence that the real
full candidate reaches the repository's final 1e-3 momentum gate.
"""

from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from .kokuno_finite_correction_cycle_contract import CorrectionFieldProvider
from .kokuno_finite_correction_cycle_ledger import (
    CycleStateProviders,
    FiniteCorrectionCycleLedgerReport,
    _analytic_cycle,
    _regression_points,
    admit_finite_correction_cycle,
)

TASK = "KOKUNO-A3-SOURCE-CORRECTION-CYCLE-GAIN-062"
SCHEMA = "kokuno-a3-source-correction-cycle-gain-v1"
PARENT_AGENT3_PR = 753
PARENT_AGENT3_HEAD = "af6a18dfdff6a2e80c311f085a8c1d8e321290bf"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

_KAPPA_S = Fraction(1, 100_000)
_REQUIRED_STAGE_GAIN = Fraction(1, 10)
_MIN_SIGNED_BUDGET = Fraction(68, 100)
_RETAINED_MEAN_BUDGET = Fraction(9, 10)
_RETAINED_RADIAL_BUDGET = Fraction(19, 10)

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


def _float(value: Fraction) -> float:
    return float(value.numerator) / float(value.denominator)


@dataclass(frozen=True)
class SourceCorrectionStageGain:
    """Frozen source exponent arithmetic for one finite correction stage."""

    stage_index: int
    kappa_s: float
    sigma_j: float
    B_j: float
    C_j: float
    H_j: float
    covariance_error_gains: tuple[float, float, float]
    covariance_min_gain: float
    wave_residual_gain: float
    tangential_mean_gain: float
    compatibility_defect_gain: float
    H_minus_B_gain: float
    signed_increment_budget: float
    temporal_increment_budget: float
    retained_mean_budget: float
    retained_radial_budget: float
    next_B: float
    next_C: float
    required_stage_gain: float
    source_bound_passed: bool

    def to_receipt(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SourceStructuredCycleGainReport:
    """Actual-defect ledger plus frozen source exponent bookkeeping."""

    cycle_id: str
    step_count: int
    source_stages: tuple[SourceCorrectionStageGain, ...]
    ledger: FiniteCorrectionCycleLedgerReport
    minimum_source_wave_gain: float
    minimum_source_mean_gain: float
    minimum_source_defect_gain: float
    source_gain_bounds_passed: bool
    empirical_actual_defect_cycle_passed: bool
    admitted: bool

    def to_receipt(self) -> dict[str, object]:
        return {
            "cycle_id": self.cycle_id,
            "step_count": self.step_count,
            "source_stages": [item.to_receipt() for item in self.source_stages],
            "ledger": self.ledger.to_receipt(),
            "minimum_source_wave_gain": self.minimum_source_wave_gain,
            "minimum_source_mean_gain": self.minimum_source_mean_gain,
            "minimum_source_defect_gain": self.minimum_source_defect_gain,
            "source_gain_bounds_passed": self.source_gain_bounds_passed,
            "empirical_actual_defect_cycle_passed": self.empirical_actual_defect_cycle_passed,
            "admitted": self.admitted,
        }


def source_correction_stage_gain(stage_index: int) -> SourceCorrectionStageGain:
    """Evaluate the corrected reader's fixed exponent arithmetic for stage j.

    The caller supplies only the integer stage label.  Every exponent, source
    constant, and required 1/10 gain is frozen here from the public reader.
    """

    if isinstance(stage_index, bool) or not isinstance(stage_index, int):
        raise TypeError("stage_index must be an integer")
    if stage_index < 0:
        raise ValueError("stage_index must be non-negative")

    j = Fraction(stage_index, 1)
    sigma = Fraction(1, 5) + j / 10
    B = Fraction(1, 2) + sigma
    C = 1 + sigma
    H = C - 2 * _KAPPA_S

    covariance = (
        Fraction(18, 100) - 2 * _KAPPA_S,
        Fraction(1, 2) - 3 * _KAPPA_S,
        sigma - 3 * _KAPPA_S,
    )
    covariance_min = min(covariance)

    # Stage-(ii)/(iii) closure arithmetic in the corrected reader.
    wave_gain = min(
        Fraction(1, 2) - 4 * _KAPPA_S,
        Fraction(4, 10) - _KAPPA_S,
        Fraction(1, 2) - 2 * _KAPPA_S,
        B - 3 * _KAPPA_S,
    )
    tangential_mean_gain = min(Fraction(17, 100), 1 - 4 * _KAPPA_S)
    compatibility_defect_gain = Fraction(9, 10) - 4 * _KAPPA_S
    h_minus_b = H - B
    signed_increment = B - _KAPPA_S
    temporal_increment = H

    # Preserve the stronger displayed source inequalities, not merely > 0.1.
    passed = bool(
        B >= Fraction(7, 10)
        and covariance_min > Fraction(17, 100)
        and wave_gain >= Fraction(4, 10) - _KAPPA_S
        and h_minus_b > _REQUIRED_STAGE_GAIN
        and tangential_mean_gain > _REQUIRED_STAGE_GAIN
        and compatibility_defect_gain > _REQUIRED_STAGE_GAIN
        and signed_increment > _MIN_SIGNED_BUDGET
        and temporal_increment > _RETAINED_MEAN_BUDGET
        and _RETAINED_RADIAL_BUDGET >= Fraction(19, 10)
    )

    return SourceCorrectionStageGain(
        stage_index=stage_index,
        kappa_s=_float(_KAPPA_S),
        sigma_j=_float(sigma),
        B_j=_float(B),
        C_j=_float(C),
        H_j=_float(H),
        covariance_error_gains=tuple(_float(value) for value in covariance),
        covariance_min_gain=_float(covariance_min),
        wave_residual_gain=_float(wave_gain),
        tangential_mean_gain=_float(tangential_mean_gain),
        compatibility_defect_gain=_float(compatibility_defect_gain),
        H_minus_B_gain=_float(h_minus_b),
        signed_increment_budget=_float(signed_increment),
        temporal_increment_budget=_float(temporal_increment),
        retained_mean_budget=_float(_RETAINED_MEAN_BUDGET),
        retained_radial_budget=_float(_RETAINED_RADIAL_BUDGET),
        next_B=_float(B + Fraction(1, 10)),
        next_C=_float(C + Fraction(1, 10)),
        required_stage_gain=_float(_REQUIRED_STAGE_GAIN),
        source_bound_passed=passed,
    )


def admit_source_structured_finite_cycle(
    states: Sequence[CycleStateProviders],
    corrections: Sequence[CorrectionFieldProvider],
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> SourceStructuredCycleGainReport:
    """Admit a raw finite cycle only if actual defect and source gains both close.

    The existing finite-cycle ledger owns the empirical divergence rejection:
    every adjacent actual NS defect is recomputed from raw providers on one
    disjoint held-in/held-out split.  Only after it admits the cycle do we attach
    the source exponent diagnostic.

    The source construction starts at j=0.  Requiring the first state to carry
    cycle_index 0 prevents a caller from relabeling the same cycle as a later,
    formally easier stage.
    """

    state_tuple = tuple(states)
    if not state_tuple:
        raise ValueError("states must be nonempty")
    if state_tuple[0].identity.cycle_index != 0:
        raise ValueError("source-structured correction cycle must start at stage index 0")

    ledger = admit_finite_correction_cycle(
        state_tuple,
        corrections,
        held_in_points,
        held_out_points,
        update_check_points,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    source_stages = tuple(
        source_correction_stage_gain(state.identity.cycle_index)
        for state in state_tuple[:-1]
    )
    source_passed = bool(
        len(source_stages) == ledger.step_count
        and all(item.source_bound_passed for item in source_stages)
    )
    empirical_passed = bool(
        ledger.every_step_actual_defect_gain_guard_passed
        and ledger.empirical_cycle_contraction_passed
    )
    if not source_passed:
        raise ValueError("corrected-source correction-cycle exponent gain bound failed")
    if not empirical_passed:
        # Defensive: admit_finite_correction_cycle already fail-closes this.
        raise ValueError("actual-defect finite correction cycle failed contraction guard")

    return SourceStructuredCycleGainReport(
        cycle_id=ledger.cycle_id,
        step_count=ledger.step_count,
        source_stages=source_stages,
        ledger=ledger,
        minimum_source_wave_gain=min(item.wave_residual_gain for item in source_stages),
        minimum_source_mean_gain=min(item.tangential_mean_gain for item in source_stages),
        minimum_source_defect_gain=min(
            item.compatibility_defect_gain for item in source_stages
        ),
        source_gain_bounds_passed=source_passed,
        empirical_actual_defect_cycle_passed=empirical_passed,
        admitted=True,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(admit_source_structured_finite_cycle)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "delta_y",
        "delta_a",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "required_gain",
        "kappa_s",
        "sigma_j",
        "B",
        "C",
        "damping_grid",
    }
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "source_kappa_s_frozen": _float(_KAPPA_S),
        "source_required_stage_gain_frozen": _float(_REQUIRED_STAGE_GAIN),
        "source_stage_B_C_derived_from_cycle_index": True,
        "source_exponent_gain_arithmetic_executable": True,
        "source_gain_composed_with_actual_defect_ledger": True,
        "each_step_actual_defect_recomputed_from_raw_providers": True,
        "divergent_empirical_step_rejected": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_source_gain_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "formal_kokuno_correction_cycle_theorem_claimed": False,
        "source_reader_is_final_independent_validation": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def deterministic_receipt() -> dict[str, object]:
    """Run a non-f=R mechanics regression through the combined gate."""

    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.25))
    report = admit_source_structured_finite_cycle(
        states,
        corrections,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )

    bad_states, bad_corrections = _analytic_cycle((1.0, 0.5, 0.75))
    empirical_expansion_rejected = False
    try:
        admit_source_structured_finite_cycle(
            bad_states,
            bad_corrections,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )
    except ValueError:
        empirical_expansion_rejected = True

    negative_stage_rejected = False
    try:
        source_correction_stage_gain(-1)
    except ValueError:
        negative_stage_rejected = True

    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "source_reader_repository": SOURCE_READER_REPO,
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_date": SOURCE_READER_DATE,
            "source_formula_summary": (
                "kappa_s=1e-5; sigma_j=1/5+j/10; B_j=1/2+sigma_j; "
                "C_j=1+sigma_j; B_{j+1}=B_j+0.1; C_{j+1}=C_j+0.1"
            ),
        },
        "analytic_mechanics_regression": {
            "field_family": "u_s=(s*t*y,0,0), p=0, fixed f=0",
            "scales": [1.0, 0.5, 0.25],
            "report": report.to_receipt(),
            "empirical_expansion_mutation_rejected": empirical_expansion_rejected,
            "negative_source_stage_rejected": negative_stage_rejected,
        },
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
