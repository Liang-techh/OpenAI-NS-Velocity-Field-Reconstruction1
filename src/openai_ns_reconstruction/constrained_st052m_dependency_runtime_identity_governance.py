"""CR002 governance for the ST052-M external execution-environment identity.

The integrated ST052-M whole-child runtime now authenticates the exact historical
Git checkout and isolates/verifies the tracked historical Python module graph.
That is a real reproducibility improvement, but it is a narrower statement than
identifying the complete executable numerical environment: the dedicated replay
workflow still selects Python 3.11 and installs NumPy/SciPy/SymPy from version
*ranges*, while neither the source-runtime identity nor the whole-child identity
binds the resolved interpreter/dependency/build identities.

This module is deliberately governance-only.  It changes no velocity field,
pressure, forcing, candidate bytes, sampling rule, norm, or scientific gate.
"""
from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
from typing import Any

from . import st052_linear_temporal_capsule as capsule
from . import st052_source_runtime_identity as source_identity

CONTRACT_RELATIVE_PATH = Path("configs/st052m_dependency_runtime_identity_contract.json")
CONSTRAINTS_RELATIVE_PATH = Path("configs/constraints.json")
WORKFLOW_RELATIVE_PATH = Path(".github/workflows/agent9-st052-linear-temporal-capsule.yml")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def load_contract(path: str | Path | None = None) -> dict[str, Any]:
    target = Path(path) if path is not None else _repo_root() / CONTRACT_RELATIVE_PATH
    return _read_json(target)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_false(mapping: dict[str, Any], key: str, *, where: str) -> None:
    _require(mapping.get(key) is False, f"{where}.{key} must remain exactly false")


def _require_true(mapping: dict[str, Any], key: str, *, where: str) -> None:
    _require(mapping.get(key) is True, f"{where}.{key} must remain exactly true")


def _audit_cr001(snapshot: dict[str, Any], constraints: dict[str, Any]) -> None:
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    expected = {
        "nu": constraints["nu"],
        "physical_domain": domain["physical"],
        "evaluation_box": domain["evaluation_box"],
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "forcing_mode": forcing["mode"],
        "forcing_a_bounds": forcing["parameters"]["a"],
        "forcing_c_bounds": forcing["parameters"]["c"],
        "reference_energy": nontriviality["reference_energy"],
        "reference_energy_abs_tolerance": nontriviality["reference_energy_abs_tolerance"],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "validation_times": validation["times"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "divergence_max": thresholds["divergence_max"],
        "divergence_L2": thresholds["divergence_L2"],
        "pde_residual_max": thresholds["pde_residual_max"],
        "pde_residual_L2": thresholds["pde_residual_L2"],
    }
    _require(snapshot == expected, "CR001 snapshot drifted from configs/constraints.json")

    _require(constraints["nu"] == 0.01, "CR001 viscosity changed")
    _require(domain["physical"] == "R^3", "CR001 physical domain changed")
    _require(domain["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "CR001 box changed")
    _require(domain["support"] == "r < 2 and abs(z) < 2", "CR001 support changed")
    _require(domain["time_interval"] == [0.25, 0.75], "CR001 time interval changed")
    _require(forcing["mode"] == "restricted_two_parameter_family", "forcing family changed")
    _require(forcing["parameters"] == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "forcing bounds changed")
    _require(nontriviality["reference_energy"] == 1.0, "reference energy changed")
    _require(nontriviality["reference_energy_abs_tolerance"] == 0.001, "energy tolerance changed")
    _require(validation["seed"] == 914027, "validation seed changed")
    _require(validation["held_out_points"] == 4096, "held-out point count changed")
    _require(validation["times"] == [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75], "validation times changed")
    _require(validation["derivative_steps"] == [0.02, 0.01, 0.005], "derivative ladder changed")
    _require(validation["quadrature_orders_per_axis"] == [24, 48, 96], "quadrature ladder changed")
    _require(thresholds["divergence_max"] == 1e-5 and thresholds["divergence_L2"] == 1e-5, "divergence gates changed")
    _require(thresholds["pde_residual_max"] == 1e-3 and thresholds["pde_residual_L2"] == 1e-3, "momentum gates changed")


def audit_contract(contract: dict[str, Any] | None = None) -> dict[str, Any]:
    """Fail closed on dependency-environment identity laundering.

    Positive integrated facts are checked from the live ST052 runtime.  The
    unresolved dependency-environment statements are checked both against the
    machine-readable contract and against the exact workflow/identity payloads.
    """
    root = _repo_root()
    contract = deepcopy(load_contract() if contract is None else contract)
    _require(contract.get("schema") == "st052m-dependency-runtime-identity-governance/v1", "unexpected contract schema")
    _require(contract.get("task_id") == "CR002", "unexpected task id")

    positives = contract.get("integrated_positive_facts")
    environment = contract.get("execution_environment_scope")
    states = contract.get("states")
    snapshot = contract.get("cr001_snapshot")
    classification = contract.get("source_classification")
    _require(isinstance(positives, dict), "integrated_positive_facts missing")
    _require(isinstance(environment, dict), "execution_environment_scope missing")
    _require(isinstance(states, dict), "states missing")
    _require(isinstance(snapshot, dict), "cr001_snapshot missing")
    _require(isinstance(classification, dict), "source_classification missing")
    _require(set(classification) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}, "source classification categories drifted")

    # Preserve the new CR-A9-070 gains instead of repeating the superseded
    # transitive-cache finding from Agent-6 PR #657.
    source_payload = source_identity.source_runtime_identity_payload()
    source_truth = source_identity.TRUTH_BOUNDARY
    whole_truth = capsule.TRUTH_BOUNDARY
    _require(source_payload.get("source_head") == source_identity.SOURCE_HEAD, "source runtime HEAD identity missing")
    _require(source_payload.get("source_tree") == source_identity.SOURCE_TREE, "source runtime tree identity missing")
    guard = source_payload.get("module_import_guard")
    _require(isinstance(guard, dict), "source runtime module-import guard missing")
    _require(guard.get("strategy") == "evict_then_verify_tracked_experiments_modules", "historical module-cache isolation strategy changed")
    _require(guard.get("loaded_module_origin_and_blob_verified") is True, "loaded historical module origin/blob verification missing")
    _require(source_truth.get("exact_source_runtime_identity_closed") is True, "integrated source runtime identity closure lost")
    _require(source_truth.get("historical_module_cache_isolated") is True, "integrated historical module-cache isolation lost")
    _require(whole_truth.get("whole_child_bundle_materialized") is True, "whole-child bundle materialization lost")
    _require(whole_truth.get("whole_child_save_load_ready_with_exact_source_runtime") is True, "authenticated-source save/load readiness lost")
    _require(whole_truth.get("exact_source_runtime_identity_closed") is True, "whole-child source runtime binding lost")
    _require(whole_truth.get("historical_module_cache_isolated") is True, "whole-child module-cache isolation lost")

    for key in (
        "whole_child_bundle_materialized",
        "whole_child_save_load_ready_with_authenticated_exact_source_runtime",
        "exact_source_git_head_tree_authenticated",
        "historical_source_module_graph_cache_isolated",
        "loaded_historical_module_origin_and_git_blob_verified",
        "exact_ci_execution_context_receipt_valid",
    ):
        _require_true(positives, key, where="integrated_positive_facts")

    # The source-code identity payload and whole-child identity bind the exact
    # historical source runtime identity, but no resolved external interpreter /
    # numerical-library identity is part of either checksum at this checkpoint.
    serialized_source = json.dumps(source_payload, sort_keys=True).lower()
    forbidden_identity_markers = (
        "python_version",
        "python_implementation",
        "numpy_version",
        "scipy_version",
        "sympy_version",
        "dependency_lock",
        "dependency_build",
        "blas",
    )
    _require(not any(marker in serialized_source for marker in forbidden_identity_markers), "source runtime identity unexpectedly started binding external dependency environment; re-audit required")

    whole_identity_parameters = set(inspect.signature(capsule._identity_payload).parameters)
    _require("exact_source_runtime_identity_sha256" in whole_identity_parameters, "whole-child identity lost source runtime binding")
    _require(not any("python" in name or "numpy" in name or "scipy" in name or "sympy" in name or "dependency" in name or "blas" in name for name in whole_identity_parameters), "whole-child identity now binds external execution environment; re-audit required")

    workflow = (root / WORKFLOW_RELATIVE_PATH).read_text(encoding="utf-8")
    _require("python-version: '3.11'" in workflow, "dedicated ST052 workflow Python selector changed; re-audit required")
    _require("'numpy>=1.24,<3'" in workflow, "dedicated ST052 NumPy range changed; re-audit required")
    _require("'scipy>=1.10,<2'" in workflow, "dedicated ST052 SciPy range changed; re-audit required")
    _require("'sympy>=1.12,<2'" in workflow, "dedicated ST052 SymPy range changed; re-audit required")

    _require(environment.get("workflow_python_selector") == "3.11", "contract Python selector drifted")
    _require(environment.get("workflow_numpy_requirement") == "numpy>=1.24,<3", "contract NumPy range drifted")
    _require(environment.get("workflow_scipy_requirement") == "scipy>=1.10,<2", "contract SciPy range drifted")
    _require(environment.get("workflow_sympy_requirement") == "sympy>=1.12,<2", "contract SymPy range drifted")
    for key in (
        "python_micro_version_bound_in_whole_candidate_identity",
        "python_implementation_build_bound_in_whole_candidate_identity",
        "numpy_resolved_version_bound_in_whole_candidate_identity",
        "scipy_resolved_version_bound_in_whole_candidate_identity",
        "sympy_resolved_version_bound_in_whole_candidate_identity",
        "binary_wheel_or_blas_build_identity_bound",
        "external_dependency_runtime_identity_bound",
        "complete_executable_environment_identity_closed",
        "exact_ci_receipt_portable_across_dependency_versions",
    ):
        _require_false(environment, key, where="execution_environment_scope")

    # Callable/save-load scope remains independent from visualization/PDE/exact
    # claims.  The current ST052 child is not promoted until its complete
    # executable environment is bound and replayed.
    for key in (
        "standalone_package_parent_runtime_ready",
        "standalone_reproducible_velocity_identity_ready",
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_false(states, key, where="states")
    _require(whole_truth.get("velocity_export_ready") is False, "live ST052 export readiness changed; re-audit required")
    _require(whole_truth.get("pde_validated") is False, "live ST052 PDE truth changed; re-audit required")
    _require(whole_truth.get("visual_correspondence_verified") is False, "live ST052 visual truth changed; re-audit required")

    constraints = _read_json(root / CONSTRAINTS_RELATIVE_PATH)
    _audit_cr001(snapshot, constraints)

    return {
        "passed": True,
        "schema": contract["schema"],
        "source_runtime_identity_sha256": source_identity.source_runtime_identity_sha256(),
        "source_runtime_code_identity_closed": True,
        "historical_module_cache_isolated": True,
        "whole_child_save_load_ready_with_authenticated_exact_source_runtime": True,
        "workflow_python_selector": environment["workflow_python_selector"],
        "workflow_dependency_requirements": {
            "numpy": environment["workflow_numpy_requirement"],
            "scipy": environment["workflow_scipy_requirement"],
            "sympy": environment["workflow_sympy_requirement"],
        },
        "external_dependency_runtime_identity_bound": False,
        "complete_executable_environment_identity_closed": False,
        "standalone_reproducible_velocity_identity_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


__all__ = ["audit_contract", "load_contract"]
