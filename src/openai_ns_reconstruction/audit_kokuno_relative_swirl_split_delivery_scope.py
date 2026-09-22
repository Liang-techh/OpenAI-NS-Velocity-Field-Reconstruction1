"""Fail-closed CR002 audit for the Kokuno relative-swirl split delivery seam.

This module changes no candidate mathematics.  Exact A1 #1171 deliberately
preserves a mathematically nonzero, sub-binary64-relative correction as separate
base and delta Cartesian channels.  That is useful executable evidence, but it
is not yet the user-required single ``velocity(x,y,z,t)->[u,v,w]`` total field.

The audit prevents evidence laundering across that representation boundary.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


CONTRACT_PATH = Path("configs/kokuno_relative_swirl_split_delivery_scope.json")
UPSTREAM_PATH = Path(
    "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_relative_swirl_split.py"
)
CONSTRAINTS_PATH = Path("configs/constraints.json")
PROJECT_STATUS_PATH = Path("project_status.json")
EXPECTED_UPSTREAM_BLOB = "c6d7d93ee33e5551513bdc817b2f118e4c544298"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_PROJECT_STATUS_BLOB = "f8e05e9d25ab6cc83a7221fc431ab42319d38ff0"
EXPECTED_A1_HEAD = "657818dd83e119e6091bd2490d3804c04c3ef723"
EXPECTED_INTEGRATION_HEAD = "6a293b3c870d70b8d8aece1b12d68685c9b8b010"


class ScopeAuditError(RuntimeError):
    """Raised when the registered CR002 representation scope drifts."""


def _git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ScopeAuditError(f"{path} must decode to a JSON object")
    return payload


def validate_contract_payload(payload: Mapping[str, Any]) -> None:
    if payload.get("schema_version") != 1:
        raise ScopeAuditError("unexpected split-delivery scope schema version")
    if payload.get("scope_id") != "cr002-kokuno-relative-swirl-split-delivery-scope-143":
        raise ScopeAuditError("unexpected split-delivery scope id")
    if payload.get("role") != "representation_and_delivery_contract_only":
        raise ScopeAuditError("scope role drifted into candidate/scientific work")

    integration = payload.get("active_integration_reference")
    if not isinstance(integration, Mapping):
        raise ScopeAuditError("active integration reference is missing")
    if dict(integration) != {
        "branch": "codex/cr001-constraints",
        "head": EXPECTED_INTEGRATION_HEAD,
    }:
        raise ScopeAuditError("active constrained integration reference drifted")

    dependency = payload.get("exact_dependency")
    if not isinstance(dependency, Mapping):
        raise ScopeAuditError("exact A1 dependency is missing")
    expected_dependency = {
        "pr": 1171,
        "head": EXPECTED_A1_HEAD,
        "source_path": str(UPSTREAM_PATH),
        "source_blob": EXPECTED_UPSTREAM_BLOB,
    }
    if dict(dependency) != expected_dependency:
        raise ScopeAuditError("exact A1 #1171 dependency identity drifted")

    cr001 = payload.get("canonical_cr001")
    if not isinstance(cr001, Mapping):
        raise ScopeAuditError("canonical CR001 contract is missing")
    expected_cr001 = {
        "path": str(CONSTRAINTS_PATH),
        "blob": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "free_residual_force_forbidden": True,
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max": 0.001,
        "momentum_L2": 0.001,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    if dict(cr001) != expected_cr001:
        raise ScopeAuditError("canonical CR001 values or gates drifted")

    eq45 = payload.get("canonical_eq45_delivery_reference")
    if not isinstance(eq45, Mapping):
        raise ScopeAuditError("canonical Eq45 delivery reference is missing")
    expected_eq45 = {
        "project_status_path": str(PROJECT_STATUS_PATH),
        "project_status_blob_on_exact_dependency": EXPECTED_PROJECT_STATUS_BLOB,
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    if dict(eq45) != expected_eq45:
        raise ScopeAuditError("independent canonical Eq45 delivery state drifted")

    classes = payload.get("provenance_classes")
    if not isinstance(classes, Mapping):
        raise ScopeAuditError("provenance classes are missing")
    if set(classes) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        raise ScopeAuditError("four-way provenance taxonomy drifted")
    if any(not isinstance(classes[key], list) or not classes[key] for key in classes):
        raise ScopeAuditError("every provenance class must remain nonempty")

    state = payload.get("materialization_state")
    if not isinstance(state, Mapping):
        raise ScopeAuditError("materialization state is missing")
    expected_state = {
        "relative_swirl_profile_delta_channel_materialized": True,
        "relative_swirl_cartesian_delta_channel_materialized": True,
        "split_channel_semantic_identity_materialized": True,
        "split_configuration_save_load_materialized": True,
        "binary64_total_relative_swirl_sum_is_resolved": False,
        "single_array_total_velocity_materialized": False,
        "current_cartesian_relative_swirl_composed": False,
        "kokuno_unified_global_cartesian_velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "heldout_complete_ns_residual_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    if dict(state) != expected_state:
        raise ScopeAuditError("split-channel delivery boundary was promoted or drifted")

    rules = payload.get("evidence_transfer_rules")
    if not isinstance(rules, Mapping) or not rules:
        raise ScopeAuditError("evidence-transfer rules are missing")
    if set(rules.values()) != {True}:
        raise ScopeAuditError("every evidence-transfer firewall must remain enabled")

    future = payload.get("future_promotion_requirements")
    if not isinstance(future, list) or len(future) < 6:
        raise ScopeAuditError("future composition requirements were weakened")

    witness = payload.get("autonomous_mechanics_witness")
    if not isinstance(witness, Mapping):
        raise ScopeAuditError("autonomous mechanics witness metadata is missing")
    if witness.get("role") != "autonomous_mechanics_only":
        raise ScopeAuditError("mechanics witness was laundered into source/candidate evidence")
    for key in ("not_source_evidence", "not_openai_numerical_evidence", "not_pde_evidence"):
        if witness.get(key) is not True:
            raise ScopeAuditError(f"mechanics witness firewall {key} was relaxed")

    forbidden = payload.get("forbidden_promotions")
    if not isinstance(forbidden, list) or len(forbidden) < 8:
        raise ScopeAuditError("forbidden-promotion list was weakened")


def _literal_assignment(module: ast.Module, name: str) -> Any:
    for node in module.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise ScopeAuditError(f"required upstream assignment {name} not found")


def _find_class(module: ast.Module, class_name: str) -> ast.ClassDef:
    for node in module.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return node
    raise ScopeAuditError(f"required upstream class {class_name} not found")


def _method_map(class_node: ast.ClassDef) -> dict[str, ast.FunctionDef]:
    return {
        child.name: child
        for child in class_node.body
        if isinstance(child, ast.FunctionDef)
    }


def audit_upstream_source(source_text: str) -> dict[str, Any]:
    module = ast.parse(source_text)
    truth = _literal_assignment(module, "_TRUTH_UPDATES")
    numerical = _literal_assignment(module, "_NUMERICAL_REALIZATION")
    if not isinstance(truth, dict) or not isinstance(numerical, dict):
        raise ScopeAuditError("upstream truth/numerical realization must remain literal dictionaries")

    expected_truth = {
        "relative_swirl_profile_delta_channel_materialized": True,
        "relative_swirl_cartesian_delta_channel_materialized": True,
        "binary64_total_relative_swirl_sum_is_resolved": False,
        "current_cartesian_relative_swirl_composed": False,
        "source_terminal_hold_after_eta_flattening_materialized": False,
        "source_exact_pointwise_bump_shape_recovered": False,
        "source_hidden_parameters_recovered": False,
        "source_exterior_heat_replacement_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "unified_global_cartesian_velocity_export_ready": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    for key, expected in expected_truth.items():
        if truth.get(key) is not expected:
            raise ScopeAuditError(f"upstream truth boundary promoted/drifted at {key}")

    representation = str(numerical.get("representation", ""))
    if "separate float64 channels" not in representation:
        raise ScopeAuditError("upstream split-channel representation declaration drifted")
    if "do not round" not in representation:
        raise ScopeAuditError("upstream anti-rounding representation warning drifted")

    class_node = _find_class(module, "KokunoPA16CurrentCartesianRelativeSwirlSplit")
    methods = _method_map(class_node)
    for required in (
        "velocity_split",
        "configuration",
        "save_configuration",
        "from_configuration",
        "load_configuration",
    ):
        if required not in methods:
            raise ScopeAuditError(f"required split representation method {required} disappeared")
    if "velocity" in methods:
        raise ScopeAuditError(
            "#1171 unexpectedly defines a single-array velocity method; governance identity must be refreshed"
        )

    split_doc = ast.get_docstring(methods["velocity_split"]) or ""
    lowered = split_doc.lower()
    if "not added in float64" not in lowered or "correction" not in lowered:
        raise ScopeAuditError("velocity_split no longer declares the binary64 composition boundary")

    return {
        "profile_delta_channel_materialized": True,
        "cartesian_delta_channel_materialized": True,
        "velocity_split_method_present": True,
        "single_array_velocity_method_defined_on_split_class": False,
        "binary64_total_sum_resolved": False,
        "current_cartesian_relative_swirl_composed": False,
    }


def autonomous_binary64_split_witness() -> dict[str, Any]:
    """Show why a nonzero split delta is not an ordinary-float total delta.

    The numbers are synthetic representation mechanics only.  They are not
    Kokuno/OpenAI values and carry no PDE meaning.
    """

    base = 1.0
    delta = math.ldexp(1.0, -60)
    if delta == 0.0:
        raise AssertionError("synthetic split delta unexpectedly vanished")
    ordinary_total = base + delta
    if ordinary_total != base:
        raise AssertionError("witness delta is not below binary64 spacing at base=1")

    encoded = json.dumps({"base": base, "delta": delta}, sort_keys=True)
    decoded = json.loads(encoded)
    if decoded["base"] != base or decoded["delta"] != delta:
        raise AssertionError("separate-channel JSON round trip did not preserve the witness channels")

    return {
        "role": "autonomous_mechanics_only",
        "base": base,
        "delta": delta,
        "delta_nonzero": True,
        "ordinary_binary64_total": ordinary_total,
        "ordinary_total_equals_base": True,
        "separate_channel_roundtrip_preserves_delta": True,
        "source_evidence": False,
        "openai_numerical_evidence": False,
        "pde_evidence": False,
    }


def _audit_constraints(constraints: Mapping[str, Any]) -> None:
    if constraints.get("nu") != 0.01:
        raise ScopeAuditError("CR001 viscosity drifted")
    domain = constraints.get("domain")
    if not isinstance(domain, Mapping):
        raise ScopeAuditError("CR001 domain contract is missing")
    if domain.get("physical") != "R^3":
        raise ScopeAuditError("CR001 physical domain drifted")
    if domain.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]]:
        raise ScopeAuditError("CR001 evaluation box drifted")
    if domain.get("support") != "r < 2 and abs(z) < 2":
        raise ScopeAuditError("CR001 support drifted")
    if domain.get("time_interval") != [0.25, 0.75]:
        raise ScopeAuditError("CR001 time interval drifted")

    forcing = constraints.get("forcing")
    if not isinstance(forcing, Mapping):
        raise ScopeAuditError("CR001 forcing contract is missing")
    if forcing.get("mode") != "restricted_two_parameter_family":
        raise ScopeAuditError("CR001 restricted forcing family drifted")
    if "No residual-dependent basis or pointwise free force" not in str(forcing.get("restriction", "")):
        raise ScopeAuditError("CR001 no-free-residual-force restriction drifted")

    nontrivial = constraints.get("nontriviality")
    if not isinstance(nontrivial, Mapping):
        raise ScopeAuditError("CR001 nontriviality contract is missing")
    if nontrivial.get("reference_energy") != 1.0:
        raise ScopeAuditError("CR001 reference energy drifted")
    if nontrivial.get("reference_energy_abs_tolerance") != 0.001:
        raise ScopeAuditError("CR001 energy tolerance drifted")

    validation = constraints.get("validation")
    if not isinstance(validation, Mapping):
        raise ScopeAuditError("CR001 validation contract is missing")
    if validation.get("seed") != 914027 or validation.get("held_out_points") != 4096:
        raise ScopeAuditError("CR001 held-out validation identity drifted")
    if validation.get("derivative_steps") != [0.02, 0.01, 0.005]:
        raise ScopeAuditError("CR001 derivative ladder drifted")
    if validation.get("quadrature_orders_per_axis") != [24, 48, 96]:
        raise ScopeAuditError("CR001 quadrature ladder drifted")
    thresholds = validation.get("thresholds")
    if not isinstance(thresholds, Mapping):
        raise ScopeAuditError("CR001 validation thresholds are missing")
    expected = {
        "pde_residual_max": 0.001,
        "pde_residual_L2": 0.001,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    for key, value in expected.items():
        if thresholds.get(key) != value:
            raise ScopeAuditError(f"CR001 threshold drifted at {key}")


def _audit_eq45_project_status(payload: Mapping[str, Any]) -> None:
    if payload.get("velocity_api") != "openai_ns_reconstruction.eq45_supported_delivery:velocity":
        raise ScopeAuditError("canonical Eq45 velocity API drifted")
    states = payload.get("states")
    if not isinstance(states, Mapping):
        raise ScopeAuditError("canonical Eq45 states are missing")
    expected = {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, value in expected.items():
        if states.get(key) is not value:
            raise ScopeAuditError(f"independent canonical Eq45 state drifted at {key}")


def audit_repository(root: str | Path) -> dict[str, Any]:
    root_path = Path(root)
    contract = _read_json(root_path / CONTRACT_PATH)
    validate_contract_payload(contract)

    upstream_bytes = (root_path / UPSTREAM_PATH).read_bytes()
    upstream_blob = _git_blob_sha1(upstream_bytes)
    if upstream_blob != EXPECTED_UPSTREAM_BLOB:
        raise ScopeAuditError(
            "exact A1 #1171 source blob drifted; a changed split representation requires new governance"
        )
    upstream = audit_upstream_source(upstream_bytes.decode("utf-8"))

    constraints_bytes = (root_path / CONSTRAINTS_PATH).read_bytes()
    constraints_blob = _git_blob_sha1(constraints_bytes)
    if constraints_blob != EXPECTED_CONSTRAINTS_BLOB:
        raise ScopeAuditError("canonical CR001 constraints blob drifted")
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    _audit_constraints(constraints)

    status_bytes = (root_path / PROJECT_STATUS_PATH).read_bytes()
    status_blob = _git_blob_sha1(status_bytes)
    if status_blob != EXPECTED_PROJECT_STATUS_BLOB:
        raise ScopeAuditError("exact-parent project_status blob drifted")
    status = json.loads(status_bytes.decode("utf-8"))
    _audit_eq45_project_status(status)

    witness = autonomous_binary64_split_witness()
    return {
        "scope_id": contract["scope_id"],
        "exact_dependency_head": EXPECTED_A1_HEAD,
        "exact_dependency_source_blob": upstream_blob,
        "canonical_constraints_blob": constraints_blob,
        "project_status_blob": status_blob,
        "upstream_scope": upstream,
        "autonomous_mechanics_witness": witness,
        "scientific_state_promoted": False,
        "candidate_bytes_changed_by_this_audit": False,
        "cr001_threshold_changed_by_this_audit": False,
    }


__all__ = [
    "ScopeAuditError",
    "audit_repository",
    "audit_upstream_source",
    "autonomous_binary64_split_witness",
    "validate_contract_payload",
]
