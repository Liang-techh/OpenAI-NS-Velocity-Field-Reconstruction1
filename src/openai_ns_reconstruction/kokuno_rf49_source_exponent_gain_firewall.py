"""Source-specific Kokuno RF49 exponent-gain firewall.

This module stacks on Agent-3 PR #1142.  The corrected 2026-09-09
reconstruction records two different notions that must not be conflated:

* asymptotic *class-exponent* gain in the finite-stage correction cycle; and
* a numerical norm ratio measured from one concrete post-update defect.

For the source value ``kappa_s = 1e-5`` the retained cycle arithmetic gives a
wave-residual exponent gain of at least ``0.39999``, a full tangential-mean
exponent gain of ``0.17``, and a defect exponent gain of ``0.89996``.  These
all exceed the source's required ``0.1`` exponent increment.  They are not
multiplicative residual contraction factors.

The public entry point rematerializes #1142's actual RF44--RF49 pre/post defect
receipt.  In addition to checking the source exponent arithmetic, it applies a
separate repository-autonomous engineering guard: a correction whose measured
RF44 defect 2-norm grows is rejected before handoff.  The non-expansion guard
is deliberately *not* presented as a Kokuno theorem or as a project PDE
acceptance threshold.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Protocol

from .kokuno_rf44_rf49_postupdate_recompute import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    RF44RF49Backend,
    RF44RF49ContractError,
    RF44RF49RecomputeReceipt,
    _RF44MechanicsBackend,
    materialize_rf44_rf49_postupdate_recompute,
)

TASK = "KOKUNO-A3-RF49-SOURCE-EXPONENT-GAIN-FIREWALL-125"
SCHEMA = "kokuno-a3-rf49-source-exponent-gain-firewall-v1"
PARENT_AGENT3_PR = 1142
PARENT_AGENT3_HEAD = "d4ff4357b1619cc13eeacfc88d453f73976fc17c"
PARENT_SOURCE_BLOB = "f6dac0be5cf88f39ca94a6066cdbf1f9270293fb"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_KAPPA_S = 1.0e-5
SOURCE_REQUIRED_EXPONENT_GAIN = 1.0e-1
SOURCE_RETAINED_TANGENTIAL_MEAN_GAIN = 1.7e-1
AUTONOMOUS_DEFECT_NORM_NONEXPANSION_LIMIT = 1.0
AUTONOMOUS_RATIO_ROUNDOFF = 1.0e-12


class RF49SourceGainError(RF44RF49ContractError):
    """Raised when source-gain provenance is not admissible."""


@dataclass(frozen=True)
class RF49SourceGainContext:
    """Candidate-bound applicability witness; no gain value is caller supplied."""

    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    evidence_kind: str
    stage_index: int
    rf49_class_hypotheses_verified: bool
    stage_exponent_bookkeeping_hypotheses_verified: bool
    stage_index_frozen_before_defect_evaluation: bool
    heldout_samples_used_to_select_stage_or_gain: bool = False
    caller_supplied_gain_used: bool = False
    residual_as_forcing_shortcut_used: bool = False


class RF49SourceGainBackend(RF44RF49Backend, Protocol):
    def rf49_source_gain_context(
        self, candidate: Any, recompute: RF44RF49RecomputeReceipt
    ) -> RF49SourceGainContext: ...


@dataclass(frozen=True)
class RF49SourceExponentArithmetic:
    stage_index: int
    sigma_j: float
    B_j: float
    C_j: float
    kappa_s: float
    wave_residual_exponent_gain: float
    tangential_mean_raw_margin_018_minus_2k: float
    tangential_mean_raw_margin_sigma_minus_3k: float
    tangential_mean_raw_margin_1_minus_4k: float
    retained_tangential_mean_exponent_gain: float
    defect_exponent_gain: float
    source_required_exponent_gain: float
    source_exponent_gain_certified: bool


@dataclass(frozen=True)
class RF49SourceGainFirewallReceipt:
    candidate_id: str
    candidate_sha256: str
    evidence_kind: str
    source_chart_id: str
    source_chart_sha256: str
    parent_recompute_candidate_state_sha256_before: str
    parent_recompute_candidate_state_sha256_after: str
    parent_rf47_closure_max_relative: float
    source_arithmetic: RF49SourceExponentArithmetic
    observed_rf44_defect_norm_ratio_after_over_before: float
    observed_defect_nonexpansion: bool
    autonomous_defect_norm_nonexpansion_limit: float
    autonomous_divergence_guard_passed: bool
    rejection_reasons: tuple[str, ...]
    source_specific_mean_handoff_accepted: bool
    repository_candidate_gain_evidence: bool
    current_i4_candidate_gain_evidence: bool
    complete_ns_cycle_gain: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, object]:
        return {**asdict(self), "truth_boundary": truth_boundary()}


def _source_arithmetic(stage_index: int) -> RF49SourceExponentArithmetic:
    if isinstance(stage_index, bool) or not isinstance(stage_index, int) or stage_index < 0:
        raise RF49SourceGainError("RF49 stage_index must be a nonnegative integer")

    k = SOURCE_KAPPA_S
    sigma = 1.0 / 5.0 + stage_index / 10.0
    B = 0.5 + sigma
    C = 1.0 + sigma
    wave = 0.4 - k
    mean_1 = 0.18 - 2.0 * k
    mean_2 = sigma - 3.0 * k
    mean_3 = 1.0 - 4.0 * k
    mean_retained = SOURCE_RETAINED_TANGENTIAL_MEAN_GAIN
    defect = 0.9 - 4.0 * k
    required = SOURCE_REQUIRED_EXPONENT_GAIN

    # The reader retains 0.17 as the full tangential-mean gain after showing
    # all displayed component margins are strictly larger.  Preserve that
    # conservative distinction instead of silently replacing 0.17 by 0.17998.
    certified = bool(
        wave > required
        and mean_retained > required
        and defect > required
        and min(mean_1, mean_2, mean_3) > mean_retained
    )
    return RF49SourceExponentArithmetic(
        stage_index=stage_index,
        sigma_j=float(sigma),
        B_j=float(B),
        C_j=float(C),
        kappa_s=k,
        wave_residual_exponent_gain=float(wave),
        tangential_mean_raw_margin_018_minus_2k=float(mean_1),
        tangential_mean_raw_margin_sigma_minus_3k=float(mean_2),
        tangential_mean_raw_margin_1_minus_4k=float(mean_3),
        retained_tangential_mean_exponent_gain=mean_retained,
        defect_exponent_gain=float(defect),
        source_required_exponent_gain=required,
        source_exponent_gain_certified=certified,
    )


def _check_context(
    context: RF49SourceGainContext, recompute: RF44RF49RecomputeReceipt
) -> None:
    if not isinstance(context, RF49SourceGainContext):
        raise RF49SourceGainError("backend must return RF49SourceGainContext")
    expected = (
        recompute.candidate.candidate_id,
        recompute.candidate.candidate_sha256,
        recompute.source_chart_id,
        recompute.source_chart_sha256,
        recompute.candidate.evidence_kind,
    )
    actual = (
        context.source_candidate_id,
        context.source_candidate_sha256,
        context.source_chart_id,
        context.source_chart_sha256,
        context.evidence_kind,
    )
    if actual != expected:
        raise RF49SourceGainError("RF49 candidate/chart/evidence identity mismatch")
    if context.evidence_kind not in {"mechanics-only", "repository-candidate"}:
        raise RF49SourceGainError("RF49 evidence_kind is invalid")
    if context.heldout_samples_used_to_select_stage_or_gain:
        raise RF49SourceGainError("RF49 stage/gain selection used held-out samples")
    if context.caller_supplied_gain_used:
        raise RF49SourceGainError("RF49 gain must be recomputed, not caller supplied")
    if context.residual_as_forcing_shortcut_used:
        raise RF49SourceGainError("RF49 provenance used forbidden residual-as-forcing shortcut")
    if not all(
        (
            context.rf49_class_hypotheses_verified,
            context.stage_exponent_bookkeeping_hypotheses_verified,
            context.stage_index_frozen_before_defect_evaluation,
        )
    ):
        raise RF49SourceGainError("RF49 source exponent applicability provenance is incomplete")


def assess_rf49_source_exponent_gain_firewall(
    backend: RF49SourceGainBackend, candidate: Any
) -> RF49SourceGainFirewallReceipt:
    """Recompute #1142, certify source exponent arithmetic, and gate handoff.

    Contract/provenance failures raise.  A measured local defect-norm increase
    is instead retained as a negative receipt with handoff rejected, so the
    failed numerical correction cannot be hidden by exception-only reporting.
    """

    if not callable(getattr(backend, "rf49_source_gain_context", None)):
        raise RF49SourceGainError("backend missing rf49_source_gain_context")

    recompute = materialize_rf44_rf49_postupdate_recompute(backend, candidate)
    context = backend.rf49_source_gain_context(candidate, recompute)
    _check_context(context, recompute)
    arithmetic = _source_arithmetic(context.stage_index)
    if not arithmetic.source_exponent_gain_certified:
        raise RF49SourceGainError("RF49 source exponent arithmetic did not clear required gain")

    ratio = float(recompute.defect_gain_ratio_2)
    if not math.isfinite(ratio) or ratio < 0.0:
        raise RF49SourceGainError("RF49 observed defect norm ratio is invalid")
    nonexpanding = ratio <= (
        AUTONOMOUS_DEFECT_NORM_NONEXPANSION_LIMIT + AUTONOMOUS_RATIO_ROUNDOFF
    )
    reasons: list[str] = []
    if not nonexpanding:
        reasons.append("observed_rf44_defect_norm_increased")

    repository_evidence = bool(
        recompute.repository_candidate_recompute_evidence
        and context.evidence_kind == "repository-candidate"
    )
    accepted = bool(arithmetic.source_exponent_gain_certified and nonexpanding)

    return RF49SourceGainFirewallReceipt(
        candidate_id=recompute.candidate.candidate_id,
        candidate_sha256=recompute.candidate.candidate_sha256,
        evidence_kind=context.evidence_kind,
        source_chart_id=recompute.source_chart_id,
        source_chart_sha256=recompute.source_chart_sha256,
        parent_recompute_candidate_state_sha256_before=recompute.pre_defect.state_sha256,
        parent_recompute_candidate_state_sha256_after=recompute.post_defect.state_sha256,
        parent_rf47_closure_max_relative=float(recompute.rf47_closure_max_relative),
        source_arithmetic=arithmetic,
        observed_rf44_defect_norm_ratio_after_over_before=ratio,
        observed_defect_nonexpansion=nonexpanding,
        autonomous_defect_norm_nonexpansion_limit=AUTONOMOUS_DEFECT_NORM_NONEXPANSION_LIMIT,
        autonomous_divergence_guard_passed=nonexpanding,
        rejection_reasons=tuple(reasons),
        source_specific_mean_handoff_accepted=accepted,
        repository_candidate_gain_evidence=repository_evidence,
        current_i4_candidate_gain_evidence=False,
        complete_ns_cycle_gain=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
    )


def truth_boundary() -> dict[str, object]:
    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_source_blob": PARENT_SOURCE_BLOB,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "source_reader_tex_path": SOURCE_READER_TEX_PATH,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "source_kappa_s": SOURCE_KAPPA_S,
        "source_required_exponent_gain": SOURCE_REQUIRED_EXPONENT_GAIN,
        "source_exponent_gain_firewall_materialized": True,
        "source_exponent_gain_is_multiplicative_norm_contraction": False,
        "observed_rf44_defect_norm_nonexpansion_guard_materialized": True,
        "observed_nonexpansion_guard_role": "repository-autonomous engineering rejection only",
        "new_project_scientific_acceptance_threshold_introduced": False,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_candidate_gain_evidence": False,
        "cartesian_correction_velocity_materialized": False,
        "complete_ns_defect": False,
        "complete_ns_cycle_gain": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "finite_stage_small_residual_is_blowup_proof": False,
    }


class _RF49MechanicsBackend(_RF44MechanicsBackend):
    """Mechanics-only source-arithmetic replay; never candidate evidence."""

    def rf49_source_gain_context(
        self, candidate: Any, recompute: RF44RF49RecomputeReceipt
    ) -> RF49SourceGainContext:
        return RF49SourceGainContext(
            source_candidate_id=recompute.candidate.candidate_id,
            source_candidate_sha256=recompute.candidate.candidate_sha256,
            source_chart_id=recompute.source_chart_id,
            source_chart_sha256=recompute.source_chart_sha256,
            evidence_kind=recompute.candidate.evidence_kind,
            stage_index=0,
            rf49_class_hypotheses_verified=True,
            stage_exponent_bookkeeping_hypotheses_verified=True,
            stage_index_frozen_before_defect_evaluation=True,
        )


def build_mechanics_report() -> dict[str, object]:
    receipt = assess_rf49_source_exponent_gain_firewall(_RF49MechanicsBackend(), None)
    return {
        "task": TASK,
        "schema": SCHEMA,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "receipt": receipt.to_dict(),
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(build_mechanics_report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
