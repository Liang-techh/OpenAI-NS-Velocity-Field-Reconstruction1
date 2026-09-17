"""Govern comparative-evidence scope for the supported Eq45 quartic Phi(1,0) trial.

This is a claim-boundary audit. It changes no velocity, force, threshold, or
candidate-selection result. The key distinction is that a favorable fixed
holdout and a seed-sensitive fresh-seed ordering are diagnostics, not the
preregistered CR001 full-momentum acceptance gate.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

EXPECTED_SCHEMA = "eq45_phi10_quartic_comparative_evidence_gate_v1"
EXPECTED_SOURCE_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
EXPECTED_CANDIDATE_SHA256 = (
    "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"
)
EXPECTED_PARENT_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)
EXPECTED_FRESH_SEEDS = (914531, 914547, 914563)


def _load(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _same_float(left: Any, right: Any, *, atol: float = 1e-12) -> bool:
    try:
        return abs(float(left) - float(right)) <= atol
    except (TypeError, ValueError):
        return False


def _fractional_change(new: float, old: float) -> float:
    _require(old != 0.0, "comparative RMS baseline must be nonzero")
    return (new - old) / abs(old)


def audit_quartic_comparative_evidence_gate(
    *,
    contract_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    delivery_state_path: str | Path | None = None,
) -> dict[str, Any]:
    """Fail closed if quartic comparative evidence is promoted outside its scope."""
    root = Path(__file__).resolve().parents[2]
    contract = _load(
        contract_path
        or root / "configs/eq45_phi10_quartic_comparative_evidence_gate.json"
    )
    constraints = _load(constraints_path or root / "configs/constraints.json")
    delivery = _load(delivery_state_path or root / "configs/delivery_state_contract.json")

    _require(contract.get("schema") == EXPECTED_SCHEMA, "governance schema drifted")

    canonical_vocab = set(delivery["classification_vocabulary"])
    contract_vocab = set(contract.get("source_vocabulary", []))
    _require(
        canonical_vocab == EXPECTED_SOURCE_VOCABULARY,
        "delivery-state source vocabulary drifted",
    )
    _require(
        contract_vocab == EXPECTED_SOURCE_VOCABULARY,
        "quartic governance source vocabulary drifted",
    )
    classifications = contract.get("source_classification", {})
    _require(
        set(classifications.values()) <= EXPECTED_SOURCE_VOCABULARY,
        "unknown source classification",
    )
    expected_classifications = {
        "callable_3d_time_varying_velocity_delivery": "user_requirement",
        "quartic_temporal_schedule_and_nullspace_balance": "autonomous_design",
        "fixed_holdout_pressure_free_vorticity_diagnostic": "autonomous_design",
        "fresh_seed_pressure_free_vorticity_generalization": "autonomous_design",
        "target_free_late_morphology_diagnostic": "autonomous_design",
        "public_visual_correspondence": "pending_unknown",
        "paper_exact_or_hidden_openai_identity": "pending_unknown",
    }
    _require(
        classifications == expected_classifications,
        "quartic evidence source classification drifted",
    )

    binding = contract["cr001_binding"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    _require(_same_float(binding["nu"], constraints["nu"]), "nu drifted")
    _require(binding["physical_domain"] == domain["physical"], "physical domain drifted")
    _require(binding["evaluation_box"] == domain["evaluation_box"], "evaluation box drifted")
    _require(binding["support"] == domain["support"], "support contract drifted")
    _require(binding["time_interval"] == domain["time_interval"], "time interval drifted")
    _require(binding["forcing_mode"] == forcing["mode"], "forcing mode drifted")
    _require(
        binding["forcing_parameters"] == sorted(forcing["parameters"]),
        "restricted forcing parameter set drifted",
    )
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "restricted forcing no-free-force guard drifted",
    )
    _require(
        int(binding["held_out_points"]) == int(validation["held_out_points"]),
        "held-out sample count drifted",
    )
    _require(
        [float(x) for x in binding["derivative_steps"]]
        == [float(x) for x in validation["derivative_steps"]],
        "derivative ladder drifted",
    )
    threshold_pairs = (
        ("divergence_max_threshold", "divergence_max"),
        ("divergence_L2_threshold", "divergence_L2"),
        ("pde_residual_max_threshold", "pde_residual_max"),
        ("pde_residual_L2_threshold", "pde_residual_L2"),
    )
    for local_key, registered_key in threshold_pairs:
        _require(
            _same_float(binding[local_key], thresholds[registered_key]),
            f"{registered_key} threshold drifted",
        )

    candidate = contract["candidate"]
    _require(
        candidate["family"]
        == "Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate",
        "quartic candidate family drifted",
    )
    _require(
        candidate["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256,
        "quartic candidate identity drifted",
    )
    _require(
        candidate["supported_parent_sha256"] == EXPECTED_PARENT_SHA256,
        "supported parent identity drifted",
    )
    _require(_same_float(candidate["early_delta"], -1.4), "early delta drifted")
    _require(
        _same_float(candidate["nullspace_coefficient"], -0.2831460674157303),
        "target-free quartic balance coefficient drifted",
    )
    _require(
        candidate["nullspace_selection"] == "target_free_derivative_balance",
        "quartic nullspace selection must remain target-free",
    )
    _require(
        candidate["pde_objective_used_to_choose_nullspace_coefficient"] is False,
        "PDE objective may not be relabeled as the quartic balance selector",
    )
    _require(
        candidate["public_image_fitted"] is False,
        "quartic schedule may not be relabeled as public-image fitted",
    )

    evidence = contract["evidence"]
    fixed = evidence["fixed_holdout"]
    fresh = evidence["fresh_seed_generalization"]
    _require(fixed["pr"] == 169, "fixed-holdout provenance drifted")
    _require(fresh["pr"] == 171, "fresh-seed provenance drifted")
    _require(
        fixed["ancestry_status"] == "current_ancestry_evidence"
        and fresh["ancestry_status"] == "current_ancestry_evidence",
        "current ancestry evidence status drifted",
    )
    _require(
        fixed["formal_pde_gate_assessed"] is False
        and fresh["formal_pde_gate_assessed"] is False,
        "pressure-free diagnostics cannot be promoted to the formal PDE gate",
    )
    _require(
        _same_float(fixed["finest_spatial_step"], validation["derivative_steps"][-1]),
        "fixed-holdout finest derivative level drifted",
    )
    _require(
        _same_float(fresh["finest_spatial_step"], validation["derivative_steps"][-1]),
        "fresh-seed finest derivative level drifted",
    )

    fixed_cubic = float(fixed["cubic_zero_rms"])
    fixed_quartic = float(fixed["quartic_zero_rms"])
    fixed_change = _fractional_change(fixed_quartic, fixed_cubic)
    _require(
        fixed_change < 0.0,
        "fixed disjoint holdout must retain the observed quartic-better ordering",
    )

    seeds = tuple(int(x) for x in fresh["seeds"])
    _require(seeds == EXPECTED_FRESH_SEEDS, "fresh held-out seeds drifted")
    rows = fresh["rows"]
    _require(len(rows) == len(EXPECTED_FRESH_SEEDS), "fresh-seed row count drifted")
    row_seeds = tuple(int(row["seed"]) for row in rows)
    _require(row_seeds == EXPECTED_FRESH_SEEDS, "fresh-seed rows drifted")

    q_vs_c = [
        _fractional_change(float(row["quartic_zero_rms"]), float(row["cubic_zero_rms"]))
        for row in rows
    ]
    q_vs_s = [
        _fractional_change(float(row["quartic_zero_rms"]), float(row["static_zero_rms"]))
        for row in rows
    ]
    quartic_better_cubic_count = sum(change < 0.0 for change in q_vs_c)
    quartic_better_static_count = sum(change < 0.0 for change in q_vs_s)
    mean_q_vs_c = sum(q_vs_c) / len(q_vs_c)
    mean_q_vs_s = sum(q_vs_s) / len(q_vs_s)
    _require(
        quartic_better_cubic_count == 1 and any(change > 0.0 for change in q_vs_c),
        "fresh-seed quartic-vs-cubic ordering must remain seed-sensitive",
    )
    _require(
        mean_q_vs_c > 0.0,
        "fresh-seed mean must retain the observed quartic-vs-cubic reversal",
    )
    _require(
        quartic_better_static_count == len(rows) and mean_q_vs_s < 0.0,
        "fresh-seed quartic-vs-static observation drifted",
    )
    _require(
        fresh["quartic_vs_cubic_sign_consistent_across_derivative_levels"] is True,
        "cross-level sign-stability observation drifted",
    )

    for key in ("late_morphology", "delivery_capsule"):
        sibling = evidence[key]
        _require(
            sibling["ancestry_status"] == "open_unconsumed_sibling_evidence",
            f"{key} must remain unconsumed sibling evidence in this ancestry",
        )
    _require(
        evidence["late_morphology"]["may_establish_public_visual_correspondence"] is False,
        "target-free morphology cannot establish public visual correspondence",
    )
    _require(
        evidence["delivery_capsule"]["may_establish_pde_validation"] is False
        and evidence["delivery_capsule"]["may_establish_public_visual_correspondence"] is False,
        "delivery packaging cannot establish scientific acceptance",
    )

    conclusions = contract["governed_conclusions"]
    _require(
        conclusions["quartic_vs_cubic_pde_ordering"] == "unresolved_seed_sensitive",
        "seed-sensitive quartic-vs-cubic PDE ordering must remain unresolved",
    )
    false_keys = (
        "single_fixed_holdout_may_establish_pde_superiority",
        "fresh_seed_pressure_free_diagnostic_may_establish_pde_validated",
        "target_free_morphology_may_establish_public_visual_correspondence",
        "green_ci_or_delivery_packaging_may_establish_scientific_acceptance",
        "pde_failure_or_pending_may_block_velocity_export",
    )
    for key in false_keys:
        _require(conclusions[key] is False, f"{key} must remain false")
    _require(
        conclusions["visual_candidate_selection_may_proceed_while_pde_unvalidated"] is True,
        "PDE status may not become a visualization-candidate delivery blocker",
    )

    states = contract["truth_states"]
    expected_states = {
        "velocity_export_ready": True,
        "visualization_candidate_only": True,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    _require(states == expected_states, "quartic truth-state boundary drifted")

    required_forbidden = {
        "fixed_holdout_quartic_better_than_cubic=>formal_pde_superiority",
        "fresh_seed_quartic_better_than_static=>pde_validated",
        "target_free_late_morphology=>public_visual_correspondence",
        "green_ci_or_replayable_capsule=>scientific_acceptance",
        "pde_failed_or_pending=>velocity_export_not_allowed",
        "visual_similarity=>paper_exact_or_openai_field_identified",
    }
    _require(
        set(contract["forbidden_inferences"]) == required_forbidden,
        "forbidden-inference set drifted",
    )

    return {
        "candidate_sha256": candidate["candidate_sha256"],
        "fixed_holdout_quartic_vs_cubic_fractional_change": fixed_change,
        "fresh_seed_quartic_vs_cubic_fractional_changes": q_vs_c,
        "fresh_seed_quartic_vs_cubic_mean_fractional_change": mean_q_vs_c,
        "fresh_seed_quartic_vs_static_mean_fractional_change": mean_q_vs_s,
        "fresh_seed_quartic_better_than_cubic_count": quartic_better_cubic_count,
        "fresh_seed_count": len(rows),
        "quartic_vs_cubic_pde_ordering": conclusions["quartic_vs_cubic_pde_ordering"],
        "velocity_export_ready": states["velocity_export_ready"],
        "pde_validated": states["pde_validated"],
        "visual_correspondence_verified": states["visual_correspondence_verified"],
    }


def main() -> int:
    print(json.dumps(audit_quartic_comparative_evidence_gate(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
