"""CR002 fail-closed audit for current l=-h bridge target-totality scope.

Representation/evidence governance only: no candidate mathematics, CR001 gate,
forcing, or scientific state is changed here.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT = Path("configs/kokuno_current_bridge_target_totality_scope.json")
SOURCE_MODULE = Path("src/openai_ns_reconstruction/kokuno_current_exterior_lminus_h_bridge.py")
CONSTRAINTS = Path("configs/constraints.json")
PROJECT_STATUS = Path("project_status.json")

_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
_CURRENT_TRUE = (
    "current_q_s_release2_endpoint_materialized",
    "current_l_minus_h_q_s_transport_materialized",
    "current_eta_resolved_pointwise_target_time_diagnostic_materialized",
)
_DOWNSTREAM_FALSE = (
    "full_eta_interval_common_scalar_bridge_length_established",
    "current_l_minus_h_matching_bridge_materialized",
    "eta_dependent_cartesian_matching_boundary_materialized",
    "current_cartesian_terminal_multiplier_composed",
    "outer_global_leading_velocity_materialized",
    "unified_global_cartesian_velocity_export_ready",
    "matched_global_pressure_materialized",
    "restricted_forcing_materialized",
    "heldout_ns_residual_assessed",
    "same_protocol_comparable_to_st006",
    "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved",
)


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AssertionError(f"{path} must contain a JSON object")
    return value


def _blob(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _literal(path: Path, name: str) -> Any:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise AssertionError(f"missing literal assignment {name} in {path}")


def _methods(path: Path, class_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return {child.name for child in node.body if isinstance(child, ast.FunctionDef)}
    raise AssertionError(f"missing class {class_name}")


def _strings(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {
        node.value for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }


def _eq(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def _mechanics_full_eta_totality_witness(contract: Mapping[str, Any]) -> dict[str, Any]:
    """Finite-hit versus continuum-totality logic; never source/PDE evidence."""
    w = contract["mechanics_witness"]
    q_p, margin, depth = float(w["q_p"]), float(w["margin"]), float(w["depth"])
    center, width, k = float(w["center"]), float(w["width"]), float(w["k"])
    samples = [float(v) for v in w["finite_eta_samples"]]
    hidden = float(w["hidden_eta_probe"])
    if not (q_p > 0 and margin > 0 and depth > margin and width > 0 and k > 0 and samples):
        raise AssertionError("mechanics witness must define a target-totality gap")

    def q_entry(eta: float) -> float:
        return q_p + margin - depth * math.exp(-((eta - center) / width) ** 2)

    q_samples = [q_entry(eta) for eta in samples]
    q_hidden = q_entry(hidden)
    sample_margins = [value - q_p for value in q_samples]
    hidden_margin = q_hidden - q_p
    if min(sample_margins) <= 0:
        raise AssertionError("finite mechanics probes must all start above Q_p")
    if hidden_margin >= 0:
        raise AssertionError("hidden mechanics probe must expose a no-hit region")
    hit_times = [math.log(value / q_p) / k for value in q_samples]
    if not all(math.isfinite(value) and value >= 0 for value in hit_times):
        raise AssertionError("finite mechanics probes must have nonnegative target hits")
    hidden_formal = math.log(q_hidden / q_p) / k
    if hidden_formal >= 0:
        raise AssertionError("hidden formal hit must lie outside y>=0")
    return {
        "finite_eta_samples": samples,
        "finite_q_entry": q_samples,
        "finite_entry_margin_min": min(sample_margins),
        "finite_target_hit_times": hit_times,
        "hidden_eta_probe": hidden,
        "hidden_q_entry": q_hidden,
        "hidden_entry_margin": hidden_margin,
        "hidden_formal_hit_time": hidden_formal,
        "full_eta_target_totality_established": False,
    }


def audit(root: str | Path | None = None) -> dict[str, Any]:
    repo = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    contract = _load_json(repo / CONTRACT)
    constraints = _load_json(repo / CONSTRAINTS)
    status = _load_json(repo / PROJECT_STATUS)
    bindings = contract["bindings"]

    _eq(contract["schema_version"], 1, "contract schema")
    _eq(set(contract["source_classification"]), _CLASSES, "four-class provenance")
    for name, entries in contract["source_classification"].items():
        if not isinstance(entries, list) or not entries or not all(isinstance(x, str) and x.strip() for x in entries):
            raise AssertionError(f"source_classification.{name} must be nonempty strings")

    _eq(_blob(repo / CONSTRAINTS), bindings["canonical_constraints_blob_sha1"], "CR001 blob")
    _eq(_blob(repo / PROJECT_STATUS), bindings["project_status_blob_sha1"], "status blob")
    _eq(_blob(repo / SOURCE_MODULE), bindings["dependency_module_blob_sha1"], "#1233 module blob")

    source_truth = _literal(repo / SOURCE_MODULE, "_TRUTH_UPDATES")
    if not isinstance(source_truth, dict):
        raise AssertionError("#1233 truth table must remain a literal dict")
    for key in _CURRENT_TRUE:
        _eq(source_truth.get(key), True, f"#1233 truth {key}")
    for key in _DOWNSTREAM_FALSE:
        _eq(source_truth.get(key), False, f"#1233 truth {key}")

    methods = _methods(repo / SOURCE_MODULE, "KokunoCurrentExteriorLMinusHBridge")
    for name in ("state", "pointwise_target_time", "pointwise_target_report"):
        if name not in methods:
            raise AssertionError(f"#1233 missing pointwise interface {name}")
    if "velocity" in methods:
        raise AssertionError("#1233 unexpectedly composes Cartesian velocity")
    literals = _strings(repo / SOURCE_MODULE)
    for marker in (
        "finite_probe_only",
        "full_eta_interval_common_scalar_bridge_length_established",
        "current_l_minus_h_matching_bridge_materialized",
        "eta_dependent_cartesian_matching_boundary_materialized",
    ):
        if marker not in literals:
            raise AssertionError(f"#1233 finite-probe marker missing: {marker}")
    expansions = _literal(repo / SOURCE_MODULE, "_ROOT_EXPANSIONS")
    bisections = _literal(repo / SOURCE_MODULE, "_ROOT_BISECTIONS")
    _eq(expansions, 12, "#1233 bracket expansions")
    _eq(bisections, 96, "#1233 bisections")

    truth = contract["truth_state"]
    for key in (
        "current_l_minus_h_q_s_transport_materialized",
        "eta_resolved_pointwise_target_time_diagnostic_materialized",
        "finite_probe_target_time_diagnostic_interface_materialized",
    ):
        _eq(truth.get(key), True, f"contract truth {key}")
    for key in (
        "full_eta_entry_above_qp_established",
        "full_eta_target_hit_totality_established",
        "full_eta_target_root_uniqueness_established",
        "full_eta_target_transversality_established",
        "smooth_eta_target_time_map_materialized",
        "full_eta_common_scalar_bridge_length_established",
        "eta_dependent_cartesian_matching_boundary_materialized",
        "current_cartesian_terminal_multiplier_composed",
        "kokuno_unified_global_velocity_export_ready",
        "visual_correspondence_verified", "pde_validated", "paper_exact",
        "openai_field_identified", "blowup_proved",
    ):
        _eq(truth.get(key), False, f"contract truth {key}")
    for key, value in contract["evidence_transfer"].items():
        _eq(value, False, f"forbidden transfer {key}")
    if len(contract["promotion_prerequisites"]) < 6:
        raise AssertionError("target-totality promotion gates were dropped")
    _eq(contract["mechanics_witness"]["role"], "autonomous_representation_mechanics_only", "witness role")
    mechanics = _mechanics_full_eta_totality_witness(contract)

    canonical = contract["cr001_guardrails"]
    _eq(constraints["nu"], canonical["nu"], "nu")
    _eq(constraints["domain"]["physical"], canonical["physical_domain"], "physical domain")
    _eq(constraints["domain"]["evaluation_box"], canonical["evaluation_box"], "evaluation box")
    _eq(constraints["domain"]["support"], canonical["support"], "support")
    _eq(constraints["domain"]["time_interval"], canonical["time_interval"], "time interval")
    _eq(constraints["forcing"]["mode"], canonical["forcing_mode"], "forcing mode")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("residual-defined/free pointwise force prohibition drifted")
    _eq(constraints["nontriviality"]["reference_energy"], canonical["reference_energy"], "reference energy")
    _eq(constraints["nontriviality"]["reference_energy_abs_tolerance"], canonical["reference_energy_abs_tolerance"], "energy tolerance")
    validation = constraints["validation"]
    _eq(validation["seed"], canonical["validation_seed"], "validation seed")
    _eq(validation["held_out_points"], canonical["held_out_points"], "heldout")
    _eq(validation["derivative_steps"], canonical["derivative_steps"], "FD ladder")
    _eq(validation["quadrature_orders_per_axis"], canonical["quadrature_orders_per_axis"], "quadrature ladder")
    gates = validation["thresholds"]
    _eq(gates["pde_residual_max"], canonical["momentum_max_threshold"], "momentum max")
    _eq(gates["pde_residual_L2"], canonical["momentum_l2_threshold"], "momentum L2")
    _eq(gates["divergence_max"], canonical["divergence_max_threshold"], "divergence max")
    _eq(gates["divergence_L2"], canonical["divergence_l2_threshold"], "divergence L2")

    states = status["states"]
    independent = contract["canonical_independence"]
    _eq(states["velocity_export_ready"], independent["canonical_eq45_velocity_export_ready_expected"], "Eq45 export")
    _eq(states["visual_correspondence_verified"], independent["canonical_visual_correspondence_verified_expected"], "visual state")
    _eq(states["pde_validated"], independent["canonical_pde_validated_expected"], "PDE state")
    _eq(states["paper_exact"], independent["canonical_paper_exact_expected"], "paper exactness")
    _eq(states["openai_field_identified"], independent["canonical_openai_field_identified_expected"], "OpenAI identity")

    return {
        "contract_id": contract["contract_id"],
        "dependency_pr": bindings["dependency_pr"],
        "dependency_exact_head": bindings["dependency_exact_head"],
        "dependency_module_blob_sha1": bindings["dependency_module_blob_sha1"],
        "pointwise_diagnostic_mechanics": {"root_expansions": expansions, "root_bisections": bisections},
        "mechanics_witness": mechanics,
        "truth_state": truth,
        "canonical_states": {
            "velocity_export_ready": states["velocity_export_ready"],
            "visual_correspondence_verified": states["visual_correspondence_verified"],
            "pde_validated": states["pde_validated"],
            "paper_exact": states["paper_exact"],
            "openai_field_identified": states["openai_field_identified"],
        },
    }
