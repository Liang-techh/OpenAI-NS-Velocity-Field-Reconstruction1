"""Independent interpretation of multi-step held-out PDE validation data.

This module does not evaluate the training loss and does not differentiate the
candidate. It consumes already-produced independent validation rows and asks a
narrow CR009 question: does refinement materially move the observed momentum
error toward the preregistered threshold, or is there a nonzero error plateau?
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping


_METRICS = (
    "residual_sampled_max",
    "residual_L2_estimate",
    "divergence_sampled_max",
    "divergence_L2_estimate",
)


def _finite_nonnegative(value: Any, name: str) -> float:
    value = float(value)
    if not (value >= 0.0 and value < float("inf")):
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _positive(value: Any, name: str) -> float:
    value = float(value)
    if not (value > 0.0 and value < float("inf")):
        raise ValueError(f"{name} must be finite and positive")
    return value


def audit_refinement(
    validation: Mapping[str, Any],
    thresholds: Mapping[str, Any],
    *,
    plateau_relative_tolerance: float = 0.10,
    threshold_margin: float = 100.0,
    divergence_contraction_limit: float = 0.50,
) -> dict[str, Any]:
    """Audit the three-or-more-level refinement already present in validation.

    The two finest levels are used only for a conservative plateau/contraction
    diagnostic. A momentum metric is called plateau-like when it is still at
    least ``threshold_margin`` times above the preregistered threshold while the
    two finest levels differ by at most ``plateau_relative_tolerance``. This is
    evidence of a numerical error floor, not a rigorous h->0 limit.

    Divergence is tracked separately because the candidate representation can be
    analytically divergence-free while finite differences converge to that
    identity. Strong fine/medium contraction therefore argues against treating
    a finite-step divergence failure as the same phenomenon as the momentum
    plateau.
    """
    plateau_relative_tolerance = _positive(
        plateau_relative_tolerance, "plateau_relative_tolerance"
    )
    threshold_margin = _positive(threshold_margin, "threshold_margin")
    divergence_contraction_limit = _positive(
        divergence_contraction_limit, "divergence_contraction_limit"
    )
    if plateau_relative_tolerance >= 1.0 or divergence_contraction_limit >= 1.0:
        raise ValueError("relative tolerances must be < 1")

    rows = validation.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("validation rows are required")

    required_thresholds = {
        "residual_sampled_max": _positive(
            thresholds.get("pde_residual_max"), "pde_residual_max"
        ),
        "residual_L2_estimate": _positive(
            thresholds.get("pde_residual_L2"), "pde_residual_L2"
        ),
        "divergence_sampled_max": _positive(
            thresholds.get("divergence_max"), "divergence_max"
        ),
        "divergence_L2_estimate": _positive(
            thresholds.get("divergence_L2"), "divergence_L2"
        ),
    }

    by_time: dict[float, dict[float, dict[str, float]]] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            raise ValueError("each validation row must be an object")
        time = float(row["time"])
        step = _positive(row["step"], "step")
        metrics = {
            name: _finite_nonnegative(row[name], name) for name in _METRICS
        }
        bucket = by_time.setdefault(time, {})
        if step in bucket:
            raise ValueError(f"duplicate time/step row: {time}/{step}")
        bucket[step] = metrics

    all_steps = sorted({step for bucket in by_time.values() for step in bucket})
    if len(all_steps) < 3:
        raise ValueError("at least three derivative steps are required")
    expected = set(all_steps)
    if any(set(bucket) != expected for bucket in by_time.values()):
        raise ValueError("every validation time must contain the same refinement levels")

    finest = all_steps[0]
    medium = all_steps[1]
    rows_out = []
    plateau_flags = []
    divergence_flags = []
    for time in sorted(by_time):
        fine = by_time[time][finest]
        med = by_time[time][medium]

        residual_rel = {}
        threshold_multiples = {}
        residual_plateau = True
        for name in ("residual_sampled_max", "residual_L2_estimate"):
            scale = max(fine[name], med[name], 1e-300)
            relative_change = abs(fine[name] - med[name]) / scale
            multiple = fine[name] / required_thresholds[name]
            residual_rel[name] = relative_change
            threshold_multiples[name] = multiple
            residual_plateau = residual_plateau and (
                relative_change <= plateau_relative_tolerance
                and multiple >= threshold_margin
            )

        divergence_ratios = {}
        divergence_contracting = True
        for name in ("divergence_sampled_max", "divergence_L2_estimate"):
            ratio = fine[name] / max(med[name], 1e-300)
            divergence_ratios[name] = ratio
            divergence_contracting = divergence_contracting and (
                ratio <= divergence_contraction_limit
            )

        plateau_flags.append(residual_plateau)
        divergence_flags.append(divergence_contracting)
        rows_out.append(
            {
                "time": time,
                "medium_step": medium,
                "finest_step": finest,
                "finest": fine,
                "momentum_relative_change_medium_to_finest": residual_rel,
                "momentum_threshold_multiples_at_finest": threshold_multiples,
                "momentum_plateau_like": bool(residual_plateau),
                "divergence_fine_over_medium": divergence_ratios,
                "divergence_still_contracting": bool(divergence_contracting),
            }
        )

    all_plateau = bool(all(plateau_flags))
    all_divergence_contracting = bool(all(divergence_flags))
    if all_plateau and all_divergence_contracting:
        classification = "momentum_mismatch_not_explained_by_current_fd_refinement"
    elif all_plateau:
        classification = "momentum_plateau_evidence_with_mixed_divergence_refinement"
    else:
        classification = "refinement_interpretation_unresolved"

    fine_residual_max = [row["finest"]["residual_sampled_max"] for row in rows_out]
    fine_residual_l2 = [row["finest"]["residual_L2_estimate"] for row in rows_out]
    return {
        "schema_version": 1,
        "source_candidate": validation.get("candidate"),
        "validation_seed": validation.get("seed"),
        "validation_points": validation.get("points"),
        "refinement_steps": all_steps,
        "medium_step": medium,
        "finest_step": finest,
        "plateau_relative_tolerance": plateau_relative_tolerance,
        "threshold_margin": threshold_margin,
        "divergence_contraction_limit": divergence_contraction_limit,
        "momentum_plateau_all_times": all_plateau,
        "divergence_contracting_all_times": all_divergence_contracting,
        "classification": classification,
        "finest_residual_sampled_max_range": [
            min(fine_residual_max),
            max(fine_residual_max),
        ],
        "finest_residual_L2_range": [min(fine_residual_l2), max(fine_residual_l2)],
        "rows": rows_out,
        "scope": (
            "post-processing of independent held-out validation rows only; "
            "does not reuse training loss, change thresholds, prove an h->0 bound, "
            "or establish PDE acceptance"
        ),
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Audit CR009 finite-difference refinement plateau"
    )
    parser.add_argument("--validation", required=True)
    parser.add_argument("--config", default="configs/constraints.json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    validation = json.loads(Path(args.validation).read_text())
    config = json.loads(Path(args.config).read_text())
    result = audit_refinement(validation, config["validation"]["thresholds"])
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(result["classification"])


if __name__ == "__main__":
    main()
