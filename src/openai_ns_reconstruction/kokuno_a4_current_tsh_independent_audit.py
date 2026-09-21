"""Independent Agent-4 audit of the current-lineage numerical PA.10 ``T_sh`` certificate.

This module deliberately validates the public candidate-side certificate through a
numerical path that is distinct from Agent 1's construction path.  The production
certificate uses Chebyshev-Lobatto construction nodes plus midpoint holdouts to
form a padded numerical ``B_0`` envelope.  Here, after save/reload, Agent 4 reads
only the public Xi ``ell_i(eta)`` surface and samples it on nested *uniform* grids
plus a fresh seeded off-grid set.  The public PA.10 algebra is then reimplemented
locally rather than calling Agent-1 helpers.

Passing this audit means only that the current numerical certificate is internally
consistent and that its selected numerical envelope covers this independent probe
set.  Finite sampling is not the source analytic C^0 theorem.  A correctly detected
geometry failure is valid negative evidence, not an audit failure and not a reason
to retune the outer scale.  This module never promotes a global Cartesian velocity,
matched pressure/forcing, complete NS residual, or ``pde_validated``.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_pa16_current_tsh_certificate import KokunoPA16CurrentTshCertificate


SCHEMA = "kokuno-a4-current-tsh-independent-audit-v1"
SEED = 9_173_641
ETA_INTERVAL = (-1.0, 1.0)
GRID_SIZES = (129, 257, 513)
OFFGRID_COUNT = 257
MAX_ENVELOPE_STABILITY = 2.0e-3
ENVELOPE_COVERAGE_ABS_SLACK = 2.0e-12
FORMULA_REL_TOL = 2.0e-13
GEOMETRY_ABS_TOL = 2.0e-11
NONTRIVIAL_ELL_FLOOR = 1.0e-12

# Re-coded public constants.  Do not import these from the production certificate:
# the point of this seam is to independently replay the displayed source algebra.
INDEPENDENT_SIGMA_PRIME_SUP = 8.0
INDEPENDENT_LOG_F_SUP = math.log(2.0)
INDEPENDENT_LOG_X_RESTORE_START = -8.0

# Final repository admission gates are frozen here only as truth-boundary metadata;
# this scoped source-coordinate audit does not assess them.
FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5


def _scale_relative_error(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(1.0, abs(float(a)), abs(float(b)))


def _uniform_eta(count: int) -> np.ndarray:
    return np.linspace(ETA_INTERVAL[0], ETA_INTERVAL[1], int(count), dtype=float)


def _offgrid_eta() -> np.ndarray:
    """Return a deterministic stratified random set disjoint from interval endpoints."""
    rng = np.random.default_rng(SEED)
    bins = np.arange(OFFGRID_COUNT, dtype=float)
    unit = (bins + rng.random(OFFGRID_COUNT)) / float(OFFGRID_COUNT)
    eta = ETA_INTERVAL[0] + (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * unit
    # Stratification already keeps all samples strictly inside each bin; sort only
    # for deterministic human-readable worst-witness diagnostics.
    return np.sort(eta.astype(float))


def _coerce_ell(values: Any, eta: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != eta.shape:
        raise ValueError(f"ell_i shape mismatch: {array.shape} != {eta.shape}")
    if np.any(~np.isfinite(array)):
        raise ValueError("ell_i contains non-finite values")
    return array


def independent_ell_observation(
    ell_evaluator: Callable[[np.ndarray], Any],
) -> dict[str, Any]:
    """Sample public ``ell_i`` on the frozen implementation-distinct protocol."""
    levels: list[dict[str, Any]] = []
    level_values: list[np.ndarray] = []
    level_eta: list[np.ndarray] = []
    for count in GRID_SIZES:
        eta = _uniform_eta(count)
        values = _coerce_ell(ell_evaluator(eta), eta)
        index = int(np.argmax(np.abs(values)))
        levels.append(
            {
                "eta_points": int(count),
                "max_abs_ell_i": float(abs(values[index])),
                "worst_eta": float(eta[index]),
                "worst_ell_i": float(values[index]),
                "rms_ell_i": float(np.sqrt(np.mean(values * values))),
            }
        )
        level_values.append(values)
        level_eta.append(eta)

    off_eta = _offgrid_eta()
    off_values = _coerce_ell(ell_evaluator(off_eta), off_eta)
    off_index = int(np.argmax(np.abs(off_values)))
    off_max = float(abs(off_values[off_index]))
    fine_max = float(levels[-1]["max_abs_ell_i"])
    medium_max = float(levels[-2]["max_abs_ell_i"])
    scale = max(fine_max, NONTRIVIAL_ELL_FLOOR)
    stability = abs(fine_max - medium_max) / scale

    nested = True
    for coarse_eta, fine_eta in zip(level_eta[:-1], level_eta[1:]):
        if not np.array_equal(coarse_eta, fine_eta[::2]):
            nested = False
            break
    monotone = bool(
        levels[0]["max_abs_ell_i"] <= levels[1]["max_abs_ell_i"] + 1.0e-15
        and levels[1]["max_abs_ell_i"] <= levels[2]["max_abs_ell_i"] + 1.0e-15
    )

    if off_max > fine_max:
        combined_max = off_max
        combined_eta = float(off_eta[off_index])
        combined_value = float(off_values[off_index])
        combined_source = "offgrid"
    else:
        combined_max = fine_max
        combined_eta = float(levels[-1]["worst_eta"])
        combined_value = float(levels[-1]["worst_ell_i"])
        combined_source = "fine_uniform"

    return {
        "levels": levels,
        "offgrid": {
            "seed": SEED,
            "eta_points": OFFGRID_COUNT,
            "max_abs_ell_i": off_max,
            "worst_eta": float(off_eta[off_index]),
            "worst_ell_i": float(off_values[off_index]),
            "rms_ell_i": float(np.sqrt(np.mean(off_values * off_values))),
        },
        "nested_uniform_grids_exact": nested,
        "nested_max_monotone": monotone,
        "medium_to_fine_scale_normalized_max_change": float(stability),
        "combined_observed_max_abs_ell_i": float(combined_max),
        "combined_worst_eta": combined_eta,
        "combined_worst_ell_i": combined_value,
        "combined_worst_source": combined_source,
    }


def _assess_report_from_observation(
    production_report: Mapping[str, Any],
    observation: Mapping[str, Any],
    *,
    X_i: float,
    log_X_R: float,
) -> dict[str, Any]:
    """Compare a production report with independently sampled/public algebra."""
    selected = production_report.get("selected_B0")
    geometry = production_report.get("geometry")
    if not isinstance(selected, Mapping) or not isinstance(geometry, Mapping):
        raise ValueError("production report is missing selected_B0/geometry mappings")

    B0 = float(selected["selected_numerical_B0_envelope"])
    observed_max = float(observation["combined_observed_max_abs_ell_i"])

    independent_lower = 160.0 * (B0 + INDEPENDENT_LOG_F_SUP)
    independent_selected_T = math.nextafter(independent_lower, math.inf)
    independent_log_x_i = math.log(float(X_i)) - float(log_X_R)
    independent_max_T = INDEPENDENT_LOG_X_RESTORE_START - independent_log_x_i
    independent_log_x_sep = independent_log_x_i + independent_selected_T
    independent_margin = independent_max_T - independent_selected_T
    independent_feasible = bool(independent_selected_T < independent_max_T)

    prod_lower = float(geometry["T_sh_lower_bound_from_selected_numerical_B0"])
    prod_selected_T = float(geometry["selected_T_sh"])
    prod_log_X_R = float(geometry["log_X_R"])
    prod_log_x_i = float(geometry["log_x_i"])
    prod_max_T = float(geometry["max_T_sh_for_current_outer_schedule"])
    prod_log_x_sep = float(geometry["selected_log_x_sep"])
    prod_margin = float(geometry["geometry_margin"])
    prod_feasible = bool(geometry["separation_geometry_feasible"])
    prod_route_ready = bool(production_report["selected_route_ready"])
    upstream_holdout_covered = bool(selected["disjoint_holdout_within_envelope"])

    checks = {
        "finite_independent_observation": bool(
            math.isfinite(observed_max)
            and all(
                math.isfinite(float(level["max_abs_ell_i"]))
                for level in observation["levels"]
            )
        ),
        "independent_ell_nontrivial": bool(observed_max > NONTRIVIAL_ELL_FLOOR),
        "nested_uniform_grids_exact": bool(observation["nested_uniform_grids_exact"]),
        "nested_max_monotone": bool(observation["nested_max_monotone"]),
        "independent_max_stable": bool(
            float(observation["medium_to_fine_scale_normalized_max_change"])
            <= MAX_ENVELOPE_STABILITY
        ),
        "selected_B0_covers_independent_probes": bool(
            B0 + ENVELOPE_COVERAGE_ABS_SLACK >= observed_max
        ),
        "production_holdout_reports_covered": upstream_holdout_covered,
        "T_sh_lower_bound_formula_matches": bool(
            _scale_relative_error(prod_lower, independent_lower) <= FORMULA_REL_TOL
        ),
        "selected_T_sh_matches_nextafter": bool(
            _scale_relative_error(prod_selected_T, independent_selected_T)
            <= FORMULA_REL_TOL
        ),
        "log_X_R_identity_matches": bool(
            abs(prod_log_X_R - float(log_X_R)) <= GEOMETRY_ABS_TOL
        ),
        "log_x_i_algebra_matches": bool(
            abs(prod_log_x_i - independent_log_x_i) <= GEOMETRY_ABS_TOL
        ),
        "max_T_sh_geometry_matches": bool(
            abs(prod_max_T - independent_max_T) <= GEOMETRY_ABS_TOL
        ),
        "selected_log_x_sep_matches": bool(
            abs(prod_log_x_sep - independent_log_x_sep) <= GEOMETRY_ABS_TOL
        ),
        "geometry_margin_matches": bool(
            abs(prod_margin - independent_margin) <= GEOMETRY_ABS_TOL
        ),
        "feasibility_boolean_matches": bool(prod_feasible == independent_feasible),
        "route_ready_boolean_matches": bool(
            prod_route_ready == (upstream_holdout_covered and independent_feasible)
        ),
    }

    audit_pass = bool(all(checks.values()))
    return {
        "production": {
            "selected_numerical_B0_envelope": B0,
            "T_sh_lower_bound": prod_lower,
            "selected_T_sh": prod_selected_T,
            "log_X_R": prod_log_X_R,
            "log_x_i": prod_log_x_i,
            "max_T_sh_for_current_outer_schedule": prod_max_T,
            "selected_log_x_sep": prod_log_x_sep,
            "geometry_margin": prod_margin,
            "separation_geometry_feasible": prod_feasible,
            "selected_route_ready": prod_route_ready,
        },
        "independent": {
            "observed_max_abs_ell_i": observed_max,
            "B0_minus_observed_max": float(B0 - observed_max),
            "T_sh_lower_bound": float(independent_lower),
            "selected_T_sh": float(independent_selected_T),
            "log_X_R": float(log_X_R),
            "log_x_i": float(independent_log_x_i),
            "max_T_sh_for_current_outer_schedule": float(independent_max_T),
            "selected_log_x_sep": float(independent_log_x_sep),
            "geometry_margin": float(independent_margin),
            "separation_geometry_feasible": independent_feasible,
        },
        "checks": checks,
        "audit_pass": audit_pass,
        # Important: this is separately reported and is not required to be True
        # for an audit of a correctly negative certificate to pass.
        "candidate_side_route_feasible": independent_feasible,
    }


def _public_ell(certificate: KokunoPA16CurrentTshCertificate, eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    values = certificate.moments.bridge.handoff_at_Xi(eta)["ell_i"]
    return _coerce_ell(values, eta)


def _mutation_firewall(
    certificate: KokunoPA16CurrentTshCertificate,
    production_report: Mapping[str, Any],
    observation: Mapping[str, Any],
) -> dict[str, bool]:
    X_i = float(certificate.moments.X_i)
    log_X_R = float(certificate.moments.outer_schedule.log_X_R)

    understated = deepcopy(dict(production_report))
    understated["selected_B0"] = dict(understated["selected_B0"])
    understated["selected_B0"]["selected_numerical_B0_envelope"] = (
        0.95 * float(understated["selected_B0"]["selected_numerical_B0_envelope"])
    )
    b0_rejected = not _assess_report_from_observation(
        understated, observation, X_i=X_i, log_X_R=log_X_R
    )["audit_pass"]

    formula_drift = deepcopy(dict(production_report))
    formula_drift["geometry"] = dict(formula_drift["geometry"])
    formula_drift["geometry"]["T_sh_lower_bound_from_selected_numerical_B0"] = (
        float(formula_drift["geometry"]["T_sh_lower_bound_from_selected_numerical_B0"])
        + 1.0e-3
    )
    formula_rejected = not _assess_report_from_observation(
        formula_drift, observation, X_i=X_i, log_X_R=log_X_R
    )["audit_pass"]

    scale_drift = deepcopy(dict(production_report))
    scale_drift["geometry"] = dict(scale_drift["geometry"])
    scale_drift["geometry"]["log_X_R"] = float(scale_drift["geometry"]["log_X_R"]) + 1.0e-3
    scale_rejected = not _assess_report_from_observation(
        scale_drift, observation, X_i=X_i, log_X_R=log_X_R
    )["audit_pass"]

    flag_drift = deepcopy(dict(production_report))
    flag_drift["geometry"] = dict(flag_drift["geometry"])
    flag_drift["geometry"]["separation_geometry_feasible"] = not bool(
        flag_drift["geometry"]["separation_geometry_feasible"]
    )
    flag_rejected = not _assess_report_from_observation(
        flag_drift, observation, X_i=X_i, log_X_R=log_X_R
    )["audit_pass"]

    amplified_observation = independent_ell_observation(
        lambda eta: 1.05 * _public_ell(certificate, eta)
    )
    ell_rejected = not _assess_report_from_observation(
        production_report, amplified_observation, X_i=X_i, log_X_R=log_X_R
    )["audit_pass"]

    return {
        "understated_B0_rejected": bool(b0_rejected),
        "T_sh_formula_drift_rejected": bool(formula_rejected),
        "outer_scale_log_X_R_drift_rejected": bool(scale_rejected),
        "feasibility_boolean_flip_rejected": bool(flag_rejected),
        "public_ell_amplitude_mutation_rejected": bool(ell_rejected),
    }


def materialize_current_tsh_independent_audit() -> dict[str, Any]:
    """Save/reload #953's public certificate and materialize the A4 receipt."""
    original = KokunoPA16CurrentTshCertificate()
    with TemporaryDirectory(prefix="kokuno-a4-tsh-") as directory:
        path = Path(directory) / "candidate_tsh_certificate.json"
        original.save_configuration(path)
        rebound = KokunoPA16CurrentTshCertificate.load_configuration(path)

    replay_identity = bool(original.semantic_sha256 == rebound.semantic_sha256)
    production_report = rebound.report()
    observation = independent_ell_observation(lambda eta: _public_ell(rebound, eta))
    assessment = _assess_report_from_observation(
        production_report,
        observation,
        X_i=float(rebound.moments.X_i),
        log_X_R=float(rebound.moments.outer_schedule.log_X_R),
    )
    mutations = _mutation_firewall(rebound, production_report, observation)

    truth = production_report.get("truth_boundary", {})
    truth_fail_closed = bool(
        truth.get("source_B0_analytic_bound_proved") is False
        and truth.get("source_T_sh_lower_bound_verified") is False
        and truth.get("inner_to_outer_join_completed") is False
        and truth.get("outer_global_leading_velocity_materialized") is False
        and truth.get("matched_global_pressure_materialized") is False
        and truth.get("restricted_forcing_materialized") is False
        and truth.get("heldout_ns_residual_assessed") is False
        and truth.get("pde_validated") is False
    )

    audit_pass = bool(
        replay_identity
        and assessment["audit_pass"]
        and all(mutations.values())
        and truth_fail_closed
    )
    return {
        "schema": SCHEMA,
        "upstream_certificate_semantic_sha256": original.semantic_sha256,
        "reloaded_certificate_semantic_sha256": rebound.semantic_sha256,
        "save_reload_semantic_identity_exact": replay_identity,
        "frozen_protocol": {
            "seed": SEED,
            "eta_interval": list(ETA_INTERVAL),
            "uniform_grid_sizes": list(GRID_SIZES),
            "offgrid_count": OFFGRID_COUNT,
            "medium_to_fine_max_stability_gate": MAX_ENVELOPE_STABILITY,
            "envelope_coverage_abs_slack": ENVELOPE_COVERAGE_ABS_SLACK,
            "formula_relative_tolerance": FORMULA_REL_TOL,
            "geometry_absolute_tolerance": GEOMETRY_ABS_TOL,
            "nontrivial_ell_floor": NONTRIVIAL_ELL_FLOOR,
            "final_momentum_gate_unassessed": FINAL_MOMENTUM_GATE,
            "final_divergence_gate_unassessed": FINAL_DIVERGENCE_GATE,
        },
        "independent_ell_observation": observation,
        "certificate_assessment": assessment,
        "mutation_firewall": mutations,
        "truth_boundary_fail_closed": truth_fail_closed,
        "candidate_side_route_feasible": assessment["candidate_side_route_feasible"],
        "source_B0_analytic_bound_proved": False,
        "source_T_sh_lower_bound_verified": False,
        "leading_only_ns_residual_assessed": False,
        "leading_plus_oscillatory_ns_residual_assessed": False,
        "after_correction_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "audit_pass": audit_pass,
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not bool(receipt.get("audit_pass")):
        failed = []
        assessment = receipt.get("certificate_assessment", {})
        for name, passed in assessment.get("checks", {}).items():
            if not passed:
                failed.append(name)
        for name, passed in receipt.get("mutation_firewall", {}).items():
            if not passed:
                failed.append(name)
        if not receipt.get("save_reload_semantic_identity_exact"):
            failed.append("save_reload_semantic_identity_exact")
        if not receipt.get("truth_boundary_fail_closed"):
            failed.append("truth_boundary_fail_closed")
        raise RuntimeError("current T_sh independent audit failed: " + ", ".join(failed))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args(argv)

    receipt = materialize_current_tsh_independent_audit()
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.enforce:
        enforce_receipt(receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
