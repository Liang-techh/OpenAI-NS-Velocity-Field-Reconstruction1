"""CR002 audit for ST052-M transitive historical-runtime import identity.

PR #649 materially improves the whole-child loader: it authenticates the exact
historical Git checkout and evicts the direct ``replay_st052`` module before the
bare-name import.  That is useful delivery evidence and is retained here.

A narrower Python-process seam remains.  The pinned historical replay imports
``minimax_exchange`` and ``spacetime`` by bare name, and the minimax path imports
``aligned_continuation``.  The #649 loader does not evict or origin-check those
transitive modules.  In a long-lived process, an already-populated
``sys.modules`` entry can therefore be reused even though the checkout itself is
cryptographically authenticated.

This module changes no velocity, candidate, pressure, forcing, optimizer,
sampling rule, norm, scientific threshold, or canonical API.  It only prevents
checkout identity/direct-module cache isolation from being silently promoted to
whole executable-runtime identity closure.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path
from typing import Any

from . import st052_linear_temporal_capsule as capsule
from . import st052_source_runtime_identity as runtime_identity

SCHEMA = "st052m-transitive-runtime-graph-governance/v1"
ACTIVE_INTEGRATION_HEAD = "458ace6f456cba8200048d2a63d102b6596fbc43"
AUDITED_UPSTREAM_HEAD = "f319b57aaaa8672ddbbaf833933b94396d3ccaad"
SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "configs"
    / "st052m_transitive_runtime_graph_contract.json"
)
TRANSITIVE_MODULES = ("minimax_exchange", "spacetime", "aligned_continuation")

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

_FALSE_STATES = (
    "transitive_runtime_import_cache_isolation_verified",
    "whole_executable_runtime_identity_closed",
    "standalone_reproducible_velocity_identity_ready",
    "experimental_st052_velocity_export_ready",
    "visualization_ready",
    "visual_correspondence_verified",
    "held_out_temporal_child_pde_residual_evaluated",
    "pde_validated",
    "source_correspondence_verified",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def load_contract(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("ST052 transitive-runtime contract must be a JSON object")
    return obj


def observe_upstream_loader() -> dict[str, Any]:
    """Observe the #649 loader mechanics that this version is scoped to."""
    source = inspect.getsource(capsule._load_exact_source_parent)
    identity_source = inspect.getsource(runtime_identity.authenticate_source_runtime)
    compact = "".join(source.split())
    transitive_mentions = {
        name: {
            "cache_eviction_mentioned": (
                f'sys.modules.pop("{name}"' in source
                or f"sys.modules.pop('{name}'" in source
            ),
            "module_name_mentioned": name in source,
        }
        for name in TRANSITIVE_MODULES
    }
    return {
        "authenticates_source_before_import": "authenticate_source_runtime(source_root)" in compact,
        "evicts_direct_replay_module": (
            'sys.modules.pop("replay_st052",None)' in compact
            or "sys.modules.pop('replay_st052',None)" in compact
        ),
        "uses_bare_replay_import": (
            'importlib.import_module("replay_st052")' in source
            or "importlib.import_module('replay_st052')" in source
        ),
        "identity_guard_checks_git_head": "rev-parse" in identity_source and "HEAD" in identity_source,
        "identity_guard_checks_git_tree": "HEAD^{tree}" in identity_source,
        "transitive_modules": transitive_mentions,
        "transitive_cache_eviction_present": any(
            item["cache_eviction_mentioned"] for item in transitive_mentions.values()
        ),
        "transitive_origin_check_present": any(
            name in source and "__file__" in source for name in TRANSITIVE_MODULES
        ),
    }


def _require(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def audit_contract(contract: dict[str, Any]) -> dict[str, Any]:
    errors: list[str] = []
    observed = observe_upstream_loader()

    _require(errors, contract.get("schema") == SCHEMA, "schema drift")
    snapshot = contract.get("snapshot")
    _require(errors, isinstance(snapshot, dict), "snapshot missing")
    if isinstance(snapshot, dict):
        _require(
            errors,
            snapshot.get("active_integration_branch") == "codex/cr001-constraints",
            "active integration branch drift",
        )
        _require(
            errors,
            snapshot.get("active_integration_head") == ACTIVE_INTEGRATION_HEAD,
            "active integration head drift",
        )
        _require(errors, snapshot.get("audited_upstream_pr") == 649, "audited PR drift")
        _require(
            errors,
            snapshot.get("audited_upstream_head") == AUDITED_UPSTREAM_HEAD,
            "audited upstream head drift",
        )

    ci = contract.get("upstream_ci")
    _require(errors, isinstance(ci, dict), "upstream_ci missing")
    if isinstance(ci, dict):
        _require(errors, ci.get("dedicated_run") == 35447803624, "dedicated run drift")
        _require(errors, ci.get("dedicated_conclusion") == "success", "dedicated CI not frozen success")
        _require(errors, ci.get("standard_run") == 35447803651, "standard run drift")
        _require(errors, ci.get("standard_conclusion") == "success", "standard CI not frozen success")

    _require(
        errors,
        runtime_identity.SOURCE_HEAD == SOURCE_HEAD,
        "exact historical source head drift",
    )
    _require(
        errors,
        observed["authenticates_source_before_import"],
        "#649 exact-source authentication is no longer observed before import",
    )
    _require(
        errors,
        observed["evicts_direct_replay_module"],
        "#649 direct replay-module cache eviction is no longer observed",
    )
    _require(
        errors,
        observed["uses_bare_replay_import"],
        "registered direct replay import mechanism changed",
    )
    _require(
        errors,
        observed["identity_guard_checks_git_head"],
        "source identity guard no longer checks Git HEAD",
    )
    _require(
        errors,
        observed["identity_guard_checks_git_tree"],
        "source identity guard no longer checks Git tree",
    )
    _require(
        errors,
        observed["transitive_cache_eviction_present"] is False,
        "transitive cache eviction now exists; v1 governance must be deliberately re-audited",
    )
    _require(
        errors,
        observed["transitive_origin_check_present"] is False,
        "transitive origin checking now exists; v1 governance must be deliberately re-audited",
    )

    graph = contract.get("runtime_graph_scope")
    _require(errors, isinstance(graph, dict), "runtime_graph_scope missing")
    if isinstance(graph, dict):
        _require(errors, graph.get("direct_module") == "replay_st052", "direct module drift")
        _require(errors, graph.get("direct_module_cache_eviction_verified") is True, "direct cache guard lost")
        _require(
            errors,
            tuple(graph.get("transitive_modules_observed", [])) == TRANSITIVE_MODULES,
            "transitive module registry drift",
        )
        _require(
            errors,
            graph.get("transitive_module_cache_eviction_verified") is False,
            "transitive cache isolation promoted without implementation",
        )
        _require(
            errors,
            graph.get("transitive_module_origin_checks_verified") is False,
            "transitive origin verification promoted without implementation",
        )
        _require(
            errors,
            graph.get("transitive_runtime_import_graph_identity_closed") is False,
            "whole transitive runtime graph promoted without implementation",
        )

    positive = contract.get("upstream_positive_facts")
    _require(errors, isinstance(positive, dict), "upstream_positive_facts missing")
    if isinstance(positive, dict):
        for key in (
            "whole_child_bundle_materialized",
            "whole_child_save_load_ready_with_authenticated_exact_source_runtime",
            "exact_source_git_head_tree_and_clean_worktree_authenticated",
            "direct_replay_module_cache_evicted_before_import",
            "direct_replay_entrypoint_covered_by_authenticated_source_tree",
            "exact_head_parity_receipt_valid",
        ):
            _require(errors, positive.get(key) is True, f"positive upstream fact lost: {key}")
        _require(
            errors,
            positive.get("standalone_package_parent_runtime_ready") is False,
            "standalone parent runtime promoted without package port",
        )

    states = contract.get("allowed_states")
    _require(errors, isinstance(states, dict), "allowed_states missing")
    if isinstance(states, dict):
        _require(errors, states.get("exact_source_checkout_identity_closed") is True, "checkout identity state lost")
        _require(errors, states.get("direct_replay_import_cache_isolation_verified") is True, "direct cache-isolation state lost")
        for key in _FALSE_STATES:
            _require(errors, states.get(key) is False, f"premature truth-state promotion: {key}")

    classes = contract.get("source_classification")
    _require(errors, isinstance(classes, dict), "source_classification missing")
    if isinstance(classes, dict):
        _require(
            errors,
            set(classes) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
            "source-classification key drift",
        )
        pending = classes.get("pending_unknown")
        _require(
            errors,
            isinstance(pending, list)
            and any("transitive" in str(item).lower() for item in pending),
            "transitive runtime identity must remain pending",
        )

    _require(
        errors,
        contract.get("cr001_nonmutation") == _EXPECTED_CR001,
        "CR001 scientific contract drift",
    )
    mutation_scope = contract.get("mutation_scope")
    _require(errors, isinstance(mutation_scope, dict), "mutation_scope missing")
    if isinstance(mutation_scope, dict):
        for key, value in mutation_scope.items():
            _require(errors, value is False, f"unexpected mutation allowed: {key}")

    return {
        "schema": SCHEMA,
        "passed": not errors,
        "errors": errors,
        "observation": observed,
        "truth_boundary": {
            "exact_source_checkout_identity_closed": True,
            "direct_replay_import_cache_isolation_verified": True,
            "transitive_runtime_import_cache_isolation_verified": False,
            "whole_executable_runtime_identity_closed": False,
            "experimental_st052_velocity_export_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def audit(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    return audit_contract(load_contract(path))


def mutated(contract: dict[str, Any], *keys: str, value: Any) -> dict[str, Any]:
    result = copy.deepcopy(contract)
    target: Any = result
    for key in keys[:-1]:
        target = target[key]
    target[keys[-1]] = value
    return result


def main() -> int:
    report = audit()
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
