"""Fail-closed CR002 governance for ST052-M exact-source import-cache scope.

This auditor does not change the candidate or runtime. It records a narrow
delivery-integrity fact exposed after the whole-child capsule was integrated:
the historical parent loader currently uses a bare ``replay_st052`` import
after prepending a source directory. In a long-lived Python process, bare-name
imports can be satisfied from ``sys.modules`` before the newly prepended path
is consulted.

The clean exact-checkout PR #632 receipt remains valid for that execution
context. This module only prevents that receipt from being promoted into a
stronger arbitrary-process runtime-identity claim before direct and transitive
module-cache isolation is actually implemented and replayed.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path
from typing import Any

from . import st052_linear_temporal_capsule as capsule

SCHEMA = "st052m-source-runtime-import-cache-governance/v1"
BASE_HEAD = "458ace6f456cba8200048d2a63d102b6596fbc43"
SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "configs"
    / "st052m_source_runtime_import_cache_contract.json"
)

_FALSE_STATES = (
    "direct_replay_import_cache_isolation_verified",
    "transitive_source_runtime_import_cache_isolation_verified",
    "runtime_identity_closed",
    "standalone_package_parent_runtime_ready",
    "standalone_reproducible_velocity_identity_ready",
    "velocity_export_ready",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)

_EXPECTED_CR001 = {
    "nu": 0.01,
    "physical_domain": "R^3",
    "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
    "support": "r < 2 and abs(z) < 2",
    "time_interval": [0.25, 0.75],
    "forcing_mode": "restricted_two_parameter_family",
    "forcing_parameter_bounds": {"a": [0.0, 10.0], "c": [0.0, 10.0]},
    "reference_energy": 1.0,
    "reference_energy_abs_tolerance": 0.001,
    "validation_seed": 914027,
    "held_out_points": 4096,
    "validation_times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
    "derivative_steps": [0.02, 0.01, 0.005],
    "quadrature_orders_per_axis": [24, 48, 96],
    "divergence_max": 1e-5,
    "divergence_L2": 1e-5,
    "pde_residual_max": 0.001,
    "pde_residual_L2": 0.001,
    "free_residual_defined_forcing_allowed": False,
    "amplitude_collapse_success_allowed": False,
}


def load_contract(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("ST052 import-cache contract must be a JSON object")
    return obj


def observe_loader() -> dict[str, bool]:
    """Observe only the direct import mechanics of the current live loader."""
    source = inspect.getsource(capsule._load_exact_source_parent)
    compact = "".join(source.split())
    return {
        "uses_sys_path_prepend": "sys.path.insert(0,inserted)" in compact,
        "uses_bare_replay_import": (
            'importlib.import_module("replay_st052")' in source
            or "importlib.import_module('replay_st052')" in source
        ),
        "references_sys_modules": "sys.modules" in source,
        "checks_imported_module_file_origin": "__file__" in source,
    }


def _append(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def audit_contract(contract: dict[str, Any]) -> dict[str, Any]:
    """Audit the contract against live loader mechanics and frozen CR001 facts."""
    errors: list[str] = []
    observation = observe_loader()

    _append(errors, contract.get("schema") == SCHEMA, "schema drift")
    base = contract.get("integration_base")
    _append(errors, isinstance(base, dict), "integration_base missing")
    if isinstance(base, dict):
        _append(
            errors,
            base.get("branch") == "codex/cr001-constraints",
            "integration branch drift",
        )
        _append(errors, base.get("head") == BASE_HEAD, "integration base head drift")

    source = contract.get("source_runtime")
    _append(errors, isinstance(source, dict), "source_runtime missing")
    if isinstance(source, dict):
        _append(errors, source.get("source_head") == SOURCE_HEAD, "source head drift")
        _append(
            errors,
            source.get("direct_module_name") == "replay_st052",
            "direct module drift",
        )
        _append(
            errors,
            source.get("direct_expected_path")
            == "experiments/root_st052/replay_st052.py",
            "direct replay path drift",
        )
        _append(
            errors,
            source.get("preexisting_direct_module_cache_rejected") is False,
            "direct import-cache rejection cannot be promoted before implementation",
        )
        _append(
            errors,
            source.get("post_import_direct_module_file_origin_checked") is False,
            "direct module-origin verification cannot be promoted before implementation",
        )
        _append(
            errors,
            source.get("transitive_module_cache_isolation_verified") is False,
            "transitive module-cache isolation cannot be promoted before implementation",
        )

    _append(
        errors,
        observation["uses_sys_path_prepend"],
        "live loader no longer matches registered sys.path-prepend mechanism",
    )
    _append(
        errors,
        observation["uses_bare_replay_import"],
        "live loader no longer matches registered bare replay import",
    )
    _append(
        errors,
        observation["references_sys_modules"] is False,
        "loader now references sys.modules; contract must be deliberately revised after auditing the new guard",
    )
    _append(
        errors,
        observation["checks_imported_module_file_origin"] is False,
        "loader now checks module origin; contract must be deliberately revised after auditing the new guard",
    )

    states = contract.get("allowed_states")
    _append(errors, isinstance(states, dict), "allowed_states missing")
    if isinstance(states, dict):
        _append(
            errors,
            states.get("whole_child_bundle_materialized") is True,
            "materialized bundle state lost",
        )
        _append(
            errors,
            states.get("whole_child_save_load_ready_with_exact_source_runtime") is True,
            "exact-source-runtime save/load receipt state lost",
        )
        _append(
            errors,
            states.get("exact_ci_execution_context_evidence_valid") is True,
            "clean exact-CI receipt must not be discarded by this narrower audit",
        )
        for key in _FALSE_STATES:
            _append(
                errors,
                states.get(key) is False,
                f"premature truth-state promotion: {key}",
            )

    evidence = contract.get("evidence_scope")
    _append(errors, isinstance(evidence, dict), "evidence_scope missing")
    if isinstance(evidence, dict):
        _append(
            errors,
            evidence.get("merged_source_pr_exact_head_dedicated_run") == 35442084024,
            "dedicated receipt drift",
        )
        _append(
            errors,
            evidence.get("merged_source_pr_exact_head_standard_run") == 35442084021,
            "standard receipt drift",
        )
        _append(
            errors,
            evidence.get("exact_source_checkout_in_dedicated_workflow") is True,
            "exact checkout receipt drift",
        )
        _append(
            errors,
            evidence.get("exact_source_parity_gate") == 5e-12,
            "parity gate drift",
        )

    classification = contract.get("classification")
    _append(errors, isinstance(classification, dict), "classification missing")
    if isinstance(classification, dict):
        _append(
            errors,
            set(classification)
            == {
                "user_requirement",
                "public_source_fact",
                "autonomous_design",
                "pending_unknown",
            },
            "source-classification keys drift",
        )
        autonomous = classification.get("autonomous_design")
        _append(
            errors,
            isinstance(autonomous, list)
            and any("sys.path" in str(item) for item in autonomous),
            "Python runtime bridge must remain classified as autonomous design",
        )
        pending = classification.get("pending_unknown")
        _append(
            errors,
            isinstance(pending, list)
            and any("sys.modules" in str(item) for item in pending),
            "import-cache closure must remain pending until implemented",
        )

    cr001 = contract.get("cr001_nonmutation")
    _append(errors, cr001 == _EXPECTED_CR001, "CR001 scientific contract drift")

    return {
        "schema": SCHEMA,
        "passed": not errors,
        "errors": errors,
        "loader_observation": observation,
        "truth_boundary": {
            "clean_exact_ci_receipt_retained": True,
            "arbitrary_process_import_cache_isolation_claimed": False,
            "runtime_identity_closed": False,
            "velocity_export_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def audit(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    return audit_contract(load_contract(path))


def mutated(contract: dict[str, Any], *keys: str, value: Any) -> dict[str, Any]:
    """Small helper used by regressions to exercise fail-closed promotions."""
    obj = copy.deepcopy(contract)
    target: Any = obj
    for key in keys[:-1]:
        target = target[key]
    target[keys[-1]] = value
    return obj


def main() -> int:
    report = audit()
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
