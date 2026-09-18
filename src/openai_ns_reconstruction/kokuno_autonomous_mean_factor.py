"""Deterministic autonomous finite-head mean factor with a hard provenance boundary.

The pinned formal Kokuno/OpenAI theorem uses a Mathlib ``ContDiffBump`` whose
numeric realization is not exported into this Python repository.  Agent-3 PR
#503 therefore correctly refuses to invent theorem-machine ``missingWeight``
values.

For repository engineering we still need a deterministic finite-head factor to
exercise downstream correction machinery.  This module supplies one, but gives
it a deliberately different name and provenance:

* ``autonomous_profile`` is an explicit C-infinity bump supported in ``(-1,1)``;
* integer translates on ``log2(q)`` are normalized in l2, so the squared masks
  sum to one up to floating roundoff;
* ``autonomous_missing_weight`` sums those squared masks over the *same integer
  head* ``m=-1,...,N-1`` used by the pinned formal definition.

This is an autonomous repository partition, **not** the theorem-machine bump,
not a replay of Mathlib's nonconstructive choice, and not evidence that the two
missing-weight values agree.  In particular this module never turns the factor
into a physical ``Delta C`` without an independently materialized same-cycle
``requestedStress``.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_actual_mean_factor_dependency import (
    BUMP_NUMERIC_DEPENDENCY,
    MISSING_WEIGHT_DEFINITION,
    PINNED_MATHLIB_COMMIT,
    PINNED_MATHLIB_REPOSITORY,
    missing_weight_indices,
    physical_scale_from_coordinate_q,
)

TASK = "KOKUNO-A3-AUTONOMOUS-MEAN-FACTOR-034"
PROFILE_VERSION = "agent3-cinf-logdyadic-l2-v1"
PROFILE_FORMULA = (
    "b(s)=exp(1-1/(1-s^2)) for |s|<1 else 0; "
    "raw_m(q)=b(log2(q)+m); mask_m=raw_m/sqrt(sum_j raw_j^2)"
)
PROFILE_SUPPORT = "|log2(q)+m|<1"


def _finite_positive(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be strictly positive and finite")
    return result


def autonomous_profile(s: float) -> float:
    """Explicit peak-normalized C-infinity compact bump on ``(-1,1)``."""

    value = float(s)
    if not math.isfinite(value):
        raise ValueError("profile coordinate must be finite")
    if abs(value) >= 1.0:
        return 0.0
    # exp(1) makes b(0)=1 while preserving compact support and smoothness.
    return math.exp(1.0 - 1.0 / (1.0 - value * value))


def autonomous_dyadic_masks(q: float) -> dict[int, float]:
    """Return all nonzero normalized log-dyadic translate weights.

    The compact support means at most two integer translates are nonzero.  We
    nevertheless enumerate a padded deterministic window to make endpoint
    behavior explicit rather than relying on an inferred pair of indices.
    """

    scale = _finite_positive(q, "q")
    log_q = math.log2(scale)
    center = -log_q
    lo = math.floor(center) - 2
    hi = math.ceil(center) + 2
    raw = {
        m: autonomous_profile(log_q + m)
        for m in range(lo, hi + 1)
    }
    raw = {m: value for m, value in raw.items() if value > 0.0}
    norm_sq = math.fsum(value * value for value in raw.values())
    if not math.isfinite(norm_sq) or norm_sq <= 0.0:
        raise ArithmeticError("autonomous translate family failed to cover q")
    norm = math.sqrt(norm_sq)
    masks = {m: value / norm for m, value in raw.items()}
    closure = math.fsum(value * value for value in masks.values())
    if abs(closure - 1.0) > 128.0 * math.ulp(1.0):
        raise ArithmeticError("autonomous squared partition lost normalization")
    return masks


def autonomous_missing_weight(prepared_n: int, q: float) -> float:
    """Sum autonomous squared masks on the source-exact finite integer head."""

    head = set(missing_weight_indices(prepared_n))
    masks = autonomous_dyadic_masks(q)
    value = math.fsum(weight * weight for m, weight in masks.items() if m in head)
    if not math.isfinite(value) or not -1.0e-15 <= value <= 1.0 + 1.0e-15:
        raise ArithmeticError("autonomous missing weight left [0,1]")
    return min(1.0, max(0.0, value))


def _profile_digest() -> str:
    payload = json.dumps(
        {
            "version": PROFILE_VERSION,
            "formula": PROFILE_FORMULA,
            "support": PROFILE_SUPPORT,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class KokunoAutonomousMeanFactor:
    """Typed autonomous factor that cannot masquerade as theorem ``missingWeight``."""

    prepared_n: int
    band: int
    coordinate_q: float

    def evaluate(self) -> dict[str, Any]:
        q_phys = physical_scale_from_coordinate_q(self.band, self.coordinate_q)
        masks = autonomous_dyadic_masks(q_phys)
        head = missing_weight_indices(self.prepared_n)
        weight = autonomous_missing_weight(self.prepared_n, q_phys)
        closure = math.fsum(value * value for value in masks.values())
        return {
            "task": TASK,
            "profile_version": PROFILE_VERSION,
            "profile_sha256": _profile_digest(),
            "profile_formula": PROFILE_FORMULA,
            "profile_support": PROFILE_SUPPORT,
            "prepared_N": int(self.prepared_n),
            "band": int(self.band),
            "coordinate_q": float(self.coordinate_q),
            "physical_q": q_phys,
            "finite_head_indices": head,
            "active_mask_indices": tuple(sorted(masks)),
            "active_mask_values": {str(m): masks[m] for m in sorted(masks)},
            "squared_partition_sum": closure,
            "squared_partition_closure_error": abs(closure - 1.0),
            "autonomous_missing_weight": weight,
            "formal_missing_weight_definition": MISSING_WEIGHT_DEFINITION,
            "formal_bump_dependency": BUMP_NUMERIC_DEPENDENCY,
            "formal_mathlib_repository": PINNED_MATHLIB_REPOSITORY,
            "formal_mathlib_commit": PINNED_MATHLIB_COMMIT,
            "profile_is_repository_autonomous": True,
            "formal_theorem_machine_bump_identity_claimed": False,
            "formal_missing_weight_equality_claimed": False,
            "theorem_missing_weight_replaced": False,
            "theorem_missing_weight_materialized": False,
            "requested_stress_actual_state_values_materialized": False,
            "finite_head_mean_debt_materialized": False,
            "surrogate_defect_used": False,
            "real_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
            "finite_correction_cycle_rerun_allowed": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "blockers": (
                "autonomous factor is not the nonconstructive formal theorem-machine missingWeight",
                "same-cycle physical requestedStress is not materialized by this factor module",
            ),
            "scope": (
                "deterministic autonomous finite-head factor for repository engineering only; "
                "it may exercise downstream code but cannot certify the formal mean debt"
            ),
        }


def write_receipt(
    output: str | Path,
    *,
    prepared_n: int,
    band: int,
    coordinate_q: float,
) -> dict[str, Any]:
    report = KokunoAutonomousMeanFactor(
        prepared_n=prepared_n,
        band=band,
        coordinate_q=coordinate_q,
    ).evaluate()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/autonomous_mean_factor_receipt.json",
    )
    parser.add_argument("--prepared-n", type=int, default=5)
    parser.add_argument("--band", type=int, default=5)
    parser.add_argument("--coordinate-q", type=float, default=1.3)
    args = parser.parse_args()
    report = write_receipt(
        args.output,
        prepared_n=args.prepared_n,
        band=args.band,
        coordinate_q=args.coordinate_q,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
