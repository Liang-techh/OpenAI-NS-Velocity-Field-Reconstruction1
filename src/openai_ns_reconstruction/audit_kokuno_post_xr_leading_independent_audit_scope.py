"""Fail-closed CR002 audit for the post-X_R leading independent-audit boundary.

This module does not evaluate or alter a velocity field. It verifies that the
machine-readable governance contract remains consistent with the exact A5 #990
registration, the parent A5 #984 through-X_R leading audit registration, the
canonical CR001 contract, and the independent Eq45 delivery state.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping

CONFIG_REL = Path("configs/kokuno_post_xr_leading_independent_audit_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
PROJECT_STATUS_REL = Path("project_status.json")
A5_990_REL = Path("src/openai_ns_reconstruction/kokuno_a5_current_exterior_xr_composite_ingest_contract.py")
A5_984_REL = Path("src/openai_ns_reconstruction/kokuno_a5_current_cartesian_exterior_xr_ingest_contract.py")
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    payload = f"blob {len(data)}\0".encode("ascii") + data
    return hashlib.sha1(payload).hexdigest()


def _load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load module spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _eq(errors: list[str], label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        errors.append(f"{label}: expected {expected!r}, got {actual!r}")


def _false(errors: list[str], label: str, value: Any) -> None:
    if value is not False:
        errors.append(f"{label}: must remain false")


def _true(errors: list[str], label: str, value: Any) -> None:
    if value is not True:
        errors.append(f"{label}: must remain true")


def audit_scope(root: Path | None = None, payload: Mapping[str, Any] | None = None) -> list[str]:
    root = _repo_root() if root is None else Path(root)
    contract = _read_json(root / CONFIG_REL) if payload is None else dict(payload)
    constraints = _read_json(root / CONSTRAINTS_REL)
    status = _read_json(root / PROJECT_STATUS_REL)
    a5_990 = _load_module(root / A5_990_REL, "_cr002_a5_990")
    a5_984 = _load_module(root / A5_984_REL, "_cr002_a5_984")

    errors: list[str] = []

    _eq(errors, "schema", contract.get("schema"), "cr002-kokuno-post-xr-leading-independent-audit-scope-v1")
    _eq(errors, "task_id", contract.get("task_id"), "CR002-KOKUNO-POST-XR-LEADING-INDEPENDENT-AUDIT-SCOPE-101")
    base = contract.get("exact_base", {})
    _eq(errors, "exact_base.pr", base.get("pr"), 990)
    _eq(errors, "exact_base.head", base.get("head"), "92def17ccab4192a7bd270a4f4d637138e953266")

    # Bind the current A5 registration instead of trusting duplicated prose.
    upstream = contract.get("upstream_scope", {})
    a1_cfg = upstream.get("agent1_986", {})
    a1_live = a5_990.AGENT1_POST_XR_FIRST_TURN_SIBLING
    for key in ("pr", "head", "source_blob", "materialized_domain"):
        _eq(errors, f"agent1_986.{key}", a1_cfg.get(key), a1_live.get(key))
    _eq(errors, "agent1_986.pr", a1_cfg.get("pr"), 986)
    _eq(errors, "agent1_986.head", a1_cfg.get("head"), "23c00b98526e187ff04d432745300a084ec859f2")
    _true(errors, "A5 #990 records A1 #986 post-X_R leading sibling", a5_990.TRUTH_BOUNDARY.get("agent1_986_leading_only_post_xr_first_turn_materialized_as_sibling"))
    _false(errors, "A1 #986 consumed by A2 #987", a1_live.get("consumed_by_agent2_987"))
    _false(errors, "A1 #986 consumed by A3 #988", a1_live.get("consumed_by_agent3_988"))
    _false(errors, "A1 #986 covered by A4 #989", a1_live.get("covered_by_agent4_989"))

    a4_989_cfg = upstream.get("agent4_989", {})
    a4_989_live = a5_990.AGENT4_XR_COMPOSITE_AUDIT
    _eq(errors, "agent4_989.pr", a4_989_cfg.get("pr"), a4_989_live.get("pr"))
    _eq(errors, "agent4_989.head", a4_989_cfg.get("head"), a4_989_live.get("head"))
    _false(errors, "A4 #989 post-X_R coverage", a4_989_live.get("covers_agent1_986_post_XR_first_turn"))
    _false(errors, "A4 #989 scientific admission", a4_989_live.get("scientific_admission"))

    a4_983_cfg = upstream.get("agent4_983", {})
    a4_983_live = a5_984.AGENT4_EXTERIOR_TO_XR_AUDIT
    _eq(errors, "agent4_983.pr", a4_983_cfg.get("pr"), a4_983_live.get("pr"))
    _eq(errors, "agent4_983.head", a4_983_cfg.get("head"), a4_983_live.get("head"))
    _eq(errors, "agent4_983 audited leading PR", a4_983_cfg.get("audited_leading_pr"), a4_983_live.get("audited_agent1_pr"))
    _eq(errors, "agent4_983 audited domain", a4_983_cfg.get("audited_domain"), "through X_R")
    _false(errors, "A4 #983 scientific admission", a4_983_live.get("scientific_admission"))

    locked = contract.get("machine_locked_distinctions", {})
    _true(errors, "post-X_R A1 materialization", locked.get("a1_986_post_xr_leading_callable_materialized"))
    _true(errors, "A1 production/source-coordinate checks exist", locked.get("a1_986_has_focused_production_and_source_coordinate_checks"))
    for key in (
        "post_xr_leading_independent_cartesian_public_velocity_audit_available",
        "through_xr_a4_audit_may_be_relabelled_as_post_xr_audit",
        "source_coordinate_handoff_or_derivative_replay_may_be_relabelled_as_independent_cartesian_post_xr_validation",
        "post_xr_leading_materialization_implies_post_xr_composite_materialization",
        "post_xr_leading_materialization_implies_global_project_domain_totality",
        "post_xr_leading_materialization_implies_pde_validation",
    ):
        _false(errors, key, locked.get(key))

    # A5 #990 must continue to stop the integrated composite at X_R.
    _false(errors, "A5 #990 post-X_R composite", a5_990.TRUTH_BOUNDARY.get("current_leading_plus_oscillatory_velocity_post_xr_materialized"))
    _false(errors, "A5 #990 independent post-X_R composite audit", a5_990.TRUTH_BOUNDARY.get("independent_post_xr_composite_audit_available"))
    _false(errors, "A5 #990 global leading", a5_990.TRUTH_BOUNDARY.get("outer_global_leading_velocity_materialized"))
    _false(errors, "A5 #990 complete candidate API", a5_990.TRUTH_BOUNDARY.get("complete_candidate_api_ready"))
    _false(errors, "A5 #990 PDE validation", a5_990.READINESS.get("pde_validated"))

    # Rebind the immutable CR001 contract, including nontriviality and validation gates.
    snapshot = contract.get("cr001_snapshot", {})
    _eq(errors, "constraints blob", _git_blob_sha1(root / CONSTRAINTS_REL), EXPECTED_CONSTRAINTS_BLOB)
    _eq(errors, "snapshot constraints blob", snapshot.get("constraints_blob_sha1"), EXPECTED_CONSTRAINTS_BLOB)
    _eq(errors, "nu", snapshot.get("nu"), constraints.get("nu"))
    domain = constraints.get("domain", {})
    _eq(errors, "physical domain", snapshot.get("physical_domain"), domain.get("physical"))
    _eq(errors, "evaluation box", snapshot.get("evaluation_box"), domain.get("evaluation_box"))
    _eq(errors, "support", snapshot.get("support"), domain.get("support"))
    _eq(errors, "time interval", snapshot.get("time_interval"), domain.get("time_interval"))
    forcing = constraints.get("forcing", {})
    _eq(errors, "forcing mode", snapshot.get("forcing_mode"), forcing.get("mode"))
    if "No residual-dependent" not in str(forcing.get("restriction", "")):
        errors.append("canonical forcing restriction no longer forbids residual-dependent free forcing")
    nontriviality = constraints.get("nontriviality", {})
    _eq(errors, "reference energy", snapshot.get("reference_energy"), nontriviality.get("reference_energy"))
    _eq(errors, "reference energy tolerance", snapshot.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"))
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    _eq(errors, "validation seed", snapshot.get("validation_seed"), validation.get("seed"))
    _eq(errors, "held-out points", snapshot.get("held_out_points"), validation.get("held_out_points"))
    _eq(errors, "derivative steps", snapshot.get("derivative_steps"), validation.get("derivative_steps"))
    _eq(errors, "quadrature ladder", snapshot.get("quadrature_orders_per_axis"), validation.get("quadrature_orders_per_axis"))
    _eq(errors, "momentum max", snapshot.get("momentum_max"), thresholds.get("pde_residual_max"))
    _eq(errors, "momentum L2", snapshot.get("momentum_l2"), thresholds.get("pde_residual_L2"))
    _eq(errors, "divergence max", snapshot.get("divergence_max"), thresholds.get("divergence_max"))
    _eq(errors, "divergence L2", snapshot.get("divergence_l2"), thresholds.get("divergence_L2"))
    _true(errors, "residual-defined forcing forbidden", snapshot.get("residual_defined_free_forcing_forbidden"))
    _true(errors, "candidate collapse forbidden", snapshot.get("candidate_collapse_forbidden"))
    _true(errors, "post-hoc threshold relaxation forbidden", snapshot.get("post_hoc_threshold_relaxation_forbidden"))

    # Kokuno incompleteness cannot downgrade the separate canonical Eq45 delivery.
    eq45 = contract.get("canonical_delivery_independence", {})
    _eq(errors, "canonical candidate family", eq45.get("candidate_family"), status.get("candidate_family"))
    states = status.get("states", {})
    for key in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _eq(errors, f"canonical state {key}", eq45.get(key), states.get(key))
    _true(errors, "canonical Eq45 velocity export", eq45.get("velocity_export_ready"))
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _false(errors, f"canonical {key}", eq45.get(key))

    provenance = contract.get("four_way_provenance", {})
    _eq(errors, "four-way provenance buckets", set(provenance), {"user_requirements", "public_source_facts", "autonomous_repository_facts", "pending_or_unknown"})
    if any("exact OpenAI field" in str(item) and "not" not in str(item).lower() for item in provenance.get("public_source_facts", [])):
        errors.append("public-source bucket improperly upgrades exact OpenAI-field identity")

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    args = parser.parse_args(argv)
    errors = audit_scope(args.root)
    print(json.dumps({"ok": not errors, "errors": errors}, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
