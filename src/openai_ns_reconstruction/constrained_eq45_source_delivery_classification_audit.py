"""CR002 audit for Eq45 source-vs-delivery classification semantics.

The public Eq. (4.1)/(4.5) formulas and the repository's callable supported
Eq45 candidate intentionally have different evidence scopes.  The source
contract still uses an older four-label vocabulary, while current delivery
governance uses ``user_requirement/public_source_fact/autonomous_design/
pending_unknown``.  This audit maps those labels fail-closed and prevents two
opposite errors:

* unresolved public-source numerical profiles must not be read as "no callable
  velocity can be delivered"; and
* an autonomous callable candidate must not be promoted to recovered/paper-
  exact/OpenAI numerical source data merely because it uses the sourced Eq45
  backbone.

No scientific parameter or velocity value is changed here.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


CONTRACT_FILENAME = "eq45_source_delivery_classification_bridge.json"
SOURCE_FILENAME = "eq45_source_contract.json"
DELIVERY_FILENAME = "velocity_delivery_contract.json"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise AssertionError(f"expected JSON object: {path}")
    return obj


def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} drifted: {actual!r} != {expected!r}")


def _source_item_by_prefix(source: dict[str, Any], prefix: str) -> dict[str, Any]:
    matches = [
        row
        for row in source.get("source_classification", [])
        if isinstance(row, dict) and str(row.get("item", "")).startswith(prefix)
    ]
    if len(matches) != 1:
        raise AssertionError(
            f"expected exactly one Eq45 source-classification item for prefix {prefix!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


def _delivery_item_by_prefix(delivery: dict[str, Any], prefix: str) -> dict[str, Any]:
    matches = [
        row
        for row in delivery.get("source_classification", [])
        if isinstance(row, dict) and str(row.get("item", "")).startswith(prefix)
    ]
    if len(matches) != 1:
        raise AssertionError(
            f"expected exactly one delivery source-classification item for prefix {prefix!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


def audit_eq45_source_delivery_classification(
    repo_root: str | Path | None = None,
) -> dict[str, Any]:
    """Fail closed on Eq45 provenance, delivery-state, or CR001 claim drift."""
    root = Path(repo_root) if repo_root is not None else _repo_root()
    contract = _load_json(root / "configs" / CONTRACT_FILENAME)
    source = _load_json(root / "configs" / SOURCE_FILENAME)
    delivery = _load_json(root / "configs" / DELIVERY_FILENAME)
    status = _load_json(root / "project_status.json")
    constraints = _load_json(root / "configs" / "constraints.json")

    _assert_equal(contract.get("schema_version"), 1, "bridge schema version")
    _assert_equal(
        contract.get("task_id"),
        "CR002-EQ45-SOURCE-DELIVERY-CLASSIFICATION-BRIDGE",
        "bridge task id",
    )

    canonical_vocab = set(contract["canonical_classification_vocabulary"])
    _assert_equal(
        canonical_vocab,
        {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "canonical four-class vocabulary",
    )
    legacy_map = contract["legacy_eq45_source_label_map"]
    _assert_equal(
        legacy_map,
        {
            "user_required": "user_requirement",
            "public_source": "public_source_fact",
            "autonomous": "autonomous_design",
            "pending": "pending_unknown",
        },
        "legacy Eq45 source-label map",
    )

    if source.get("task_id") != "CR001-EQ45-SOURCE-GOV-002":
        raise AssertionError("unexpected Eq45 public-source contract identity")
    if source["velocity_contract"].get("leading_field_only") is not True:
        raise AssertionError("Eq45 source formula must remain scoped to the leading field")

    for expected in contract["eq45_source_items"]:
        row = _source_item_by_prefix(source, expected["item_prefix"])
        _assert_equal(
            row.get("classification"),
            expected["legacy_classification"],
            f"legacy classification for {expected['item_prefix']}",
        )
        mapped = legacy_map.get(row["classification"])
        _assert_equal(
            mapped,
            expected["canonical_classification"],
            f"canonical classification for {expected['item_prefix']}",
        )
        if mapped not in canonical_vocab:
            raise AssertionError("mapped Eq45 source classification left canonical vocabulary")

    source_claims = source["claim_status"]
    if source_claims.get("coordinate_formula_publicly_sourced") is not True:
        raise AssertionError("public Eq45 coordinate source fact regressed")
    if source_claims.get("leading_velocity_formula_publicly_sourced") is not True:
        raise AssertionError("public Eq45 leading-velocity source fact regressed")
    for key in (
        "final_profiles_numerically_identified",
        "full_constructed_field_recovered",
        "openai_correspondence_verified",
        "pde_validated",
        "paper_exact",
    ):
        if source_claims.get(key) is not False:
            raise AssertionError(f"Eq45 public-source claim promoted without evidence: {key}")

    bridge = contract["source_to_delivery_bridge"]
    _assert_equal(
        bridge["public_source_formula_available"],
        source_claims["leading_velocity_formula_publicly_sourced"],
        "public-source formula availability",
    )
    _assert_equal(
        bridge["public_source_formula_is_leading_field_only"],
        source["velocity_contract"]["leading_field_only"],
        "leading-field-only scope",
    )
    _assert_equal(
        bridge["public_source_final_profiles_numerically_identified"],
        source_claims["final_profiles_numerically_identified"],
        "public-source final-profile identification state",
    )
    _assert_equal(
        bridge["public_source_full_constructed_field_recovered"],
        source_claims["full_constructed_field_recovered"],
        "public-source full-field recovery state",
    )

    primary = delivery["primary_deliverable"]
    _assert_equal(
        primary["candidate_family"], bridge["canonical_candidate_family"], "candidate family"
    )
    _assert_equal(
        primary["candidate_sha256"], bridge["canonical_candidate_sha256"], "candidate sha256"
    )
    _assert_equal(primary["api"], bridge["canonical_velocity_api"], "canonical velocity API")
    if primary.get("canonical") is not True:
        raise AssertionError("Eq45 supported delivery must remain the canonical callable delivery")

    delivery_candidate_row = _delivery_item_by_prefix(
        delivery, "Eq45 support-connected candidate identity, API, export format, and delivery wrapper"
    )
    _assert_equal(
        delivery_candidate_row.get("classification"),
        bridge["canonical_candidate_numeric_instantiation_class"],
        "canonical candidate provenance class",
    )
    _assert_equal(
        bridge["canonical_candidate_numeric_instantiation_class"],
        "autonomous_design",
        "canonical numerical instantiation class",
    )

    hidden_row = _delivery_item_by_prefix(
        delivery, "exact OpenAI numerical velocity samples, hidden coefficients"
    )
    _assert_equal(hidden_row.get("classification"), "pending_unknown", "hidden OpenAI data class")

    claims = delivery["claim_status"]
    _assert_equal(
        claims.get("velocity_export_ready"),
        bridge["canonical_velocity_export_ready"],
        "canonical velocity export readiness",
    )
    if bridge["callable_delivery_requires_public_source_profile_identification"] is not False:
        raise AssertionError("unresolved public-source profiles must not block callable delivery")
    if bridge["callable_delivery_requires_pde_validation"] is not False:
        raise AssertionError("PDE validation must not be turned into a callable-delivery prerequisite")
    if bridge["source_formula_alone_identifies_canonical_numeric_candidate"] is not False:
        raise AssertionError("source formula alone must not identify the autonomous numerical candidate")
    if bridge["canonical_candidate_may_be_called_exact_openai_field"] is not False:
        raise AssertionError("canonical autonomous candidate must not be called the exact OpenAI field")

    boundary = contract["independent_truth_states"]
    _assert_equal(boundary, claims, "bridge versus delivery independent truth states")
    if boundary["velocity_export_ready"] is not True:
        raise AssertionError("callable/exportable canonical delivery readiness regressed")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if boundary.get(key) is not False:
            raise AssertionError(f"scientific/exact truth state promoted without evidence: {key}")

    _assert_equal(status["candidate_family"], primary["candidate_family"], "project candidate family")
    _assert_equal(status["candidate_sha256"], primary["candidate_sha256"], "project candidate sha256")
    _assert_equal(status["velocity_api"], primary["api"], "project velocity API")
    for key, expected in boundary.items():
        if key not in status["states"]:
            raise AssertionError(f"project status is missing independent truth state {key}")
        _assert_equal(status["states"][key], expected, f"project truth state {key}")

    frozen = contract["canonical_cr001"]
    _assert_equal(constraints["nu"], frozen["nu"], "nu")
    _assert_equal(constraints["domain"]["physical"], frozen["physical_domain"], "physical domain")
    _assert_equal(constraints["domain"]["evaluation_box"], frozen["evaluation_box"], "evaluation box")
    _assert_equal(constraints["domain"]["support"], frozen["support"], "support")
    _assert_equal(constraints["domain"]["time_interval"], frozen["time_interval"], "time interval")
    _assert_equal(constraints["forcing"]["mode"], frozen["force_mode"], "forcing mode")
    _assert_equal(constraints["forcing"]["parameters"], frozen["force_bounds"], "forcing bounds")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("residual-defined/free-force prohibition drifted")
    _assert_equal(
        constraints["nontriviality"]["reference_energy"], frozen["reference_energy"], "reference energy"
    )
    _assert_equal(
        constraints["nontriviality"]["reference_energy_abs_tolerance"],
        frozen["reference_energy_abs_tolerance"],
        "reference energy tolerance",
    )
    if "reject collapsed candidates" not in constraints["nontriviality"]["enforcement"]:
        raise AssertionError("amplitude-collapse rejection drifted")
    validation = constraints["validation"]
    _assert_equal(validation["seed"], frozen["validation_seed"], "validation seed")
    _assert_equal(validation["held_out_points"], frozen["held_out_points"], "held-out points")
    _assert_equal(validation["times"], frozen["validation_times"], "validation times")
    _assert_equal(validation["derivative_steps"], frozen["derivative_steps"], "derivative ladder")
    _assert_equal(
        validation["quadrature_orders_per_axis"],
        frozen["quadrature_orders_per_axis"],
        "quadrature ladder",
    )
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _assert_equal(thresholds[key], frozen[key], key)

    return {
        "task_id": contract["task_id"],
        "legacy_source_labels_mapped_to_canonical_four_classes": True,
        "public_source_and_autonomous_numeric_instantiation_kept_distinct": True,
        "callable_delivery_kept_independent_from_pde_and_exactness": True,
        "cr001_checked_unchanged": True,
        "truth_boundary": boundary,
    }
