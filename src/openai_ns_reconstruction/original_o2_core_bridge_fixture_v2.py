"""Synthetic-only v2 fixture for the static O2-to-core bridge.

Every mapping and source snippet in this module is made locally.  The fixture
does not open repository paths, import a producer, decode a receipt, or call a
source function.  It supplies the five arguments expected by the v2 leaf and
keeps each guard independent so a failed guard retains its raw diagnostic.
"""

from __future__ import annotations

import copy
import traceback
from collections.abc import Callable, Mapping
from typing import Any


SCHEMA = "original_o2_core_bridge_fixture_v2"
LEAF_SCHEMA = "original_o2_core_bridge_v2"


def _o2_packet() -> dict[str, Any]:
    """Return a fresh synthetic O2 evidence packet with the real pin spelling."""

    return {
        "source_definition_binding": {
            "passed": True,
            "synthetic_frame_id": "cobalt-orbit-frame-v2",
        },
        "delta_branch_proof": {
            "passed": True,
            "selected_branch": "exp(-4logPstar-30)",
            "no_rounded_branch_comparison": True,
            "branch_witness": "unrounded-lantern-v2",
        },
        "input_hashes": {
            "lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame.py": "dc2378e7cf690c879b48cb0ba0bbccb1fcf00d26b563269ca447fcc03f5c8a2c",
        },
    }


def _core_packet() -> dict[str, Any]:
    """Return a fresh synthetic core identity-gate packet."""

    return {
        "all_passed": True,
        "checks": {
            "actual_oracle": {
                "passed": True,
                "exact_parameter_and_native_P0_identity": {
                    "passed": True,
                    "exact_graph_logP_delta_logLambda_epsilon_recipes_verified": True,
                    "actual_core_encloses_same_exact_recipes": {
                        "logP": True,
                        "logLambda": True,
                        "delta": True,
                        "Lambda": True,
                        "epsilon": True,
                    },
                    "interval_overlap_not_used_as_a_function_identity": True,
                    "synthetic_gate_marker": "violet-circuit-v2",
                },
            },
        },
        "input_hashes": {
            "lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame.py": "dc2378e7cf690c879b48cb0ba0bbccb1fcf00d26b563269ca447fcc03f5c8a2c",
            "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
        },
    }


def _lineage_packet() -> dict[str, Any]:
    """Return a deliberately nested synthetic lineage packet."""

    return {
        "receipt_wrapper": "nested-coral-envelope-v2",
        "result": {
            "result": {
                "leaf": {
                    "schema": "coral_lineage_probe_v2",
                    "status": "proved",
                    "passed": True,
                    "comparisons": {"overlap_used_as_identity": False},
                    "provenance_scope": {
                        "cross_run_object_identity_claimed": False,
                        "live_endpoint_basis_identity": False,
                        "native_corrected_RQy_identity": False,
                        "source_function_point_value_identity": False,
                        "interval_overlap_promoted_to_identity": False,
                    },
                    "lineage_marker": "nested-coral-leaf-v2",
                }
            }
        },
    }


def _core_source_text() -> str:
    """Return a new source spelling of the four immutable recipe constants."""

    return '''
class SapphireKernelIdentity:
    def parameter_identity(self, witness_token):
        g = self.operator
        c = g.constant
        p = self.parameter_frame
        expected_logP = g.add(g.unary("exp", c(40)), c(11))
        expected_delta = g.unary("exp", g.neg(g.add(g.mul(c(4), p["logP"]), c(30))))
        if p["logP"] != expected_logP or p["delta"] != expected_delta:
            raise RuntimeError("coral recipe witness rejected")
'''


def _bridge_source_text() -> str:
    """Return a bridge source with several local aliases and new identifiers."""

    return '''
class AzureRelayScaffold:
    def __init__(self, pressure_packet, transfer_packet):
        boundary_alias = pressure_packet
        handoff_alias = transfer_packet
        core_transfer = handoff_alias
        star_route = core_transfer
        delta_value = exp(boundary_alias.parameter_bounds.log_delta)
        logp_value = star_route.logPstar
        self.core = core_container(delta=delta_value, logP=logp_value)
'''


def _cross_object_bridge_text() -> str:
    """Return a data-flow-equivalent source using two distinct packet inputs."""

    return '''
class TwinOrbitRelay:
    def __init__(self, pressure_packet_left, transfer_packet_right):
        boundary_alias_left = pressure_packet_left
        handoff_alias_right = transfer_packet_right
        core_transfer = handoff_alias_right
        star_route = core_transfer
        delta_value = exp(boundary_alias_left.parameter_bounds.log_delta)
        logp_value = star_route.logPstar
        self.core = core_container(delta=delta_value, logP=logp_value)
'''


def _failed(result: Any, statuses: str | tuple[str, ...]) -> bool:
    expected = (statuses,) if isinstance(statuses, str) else statuses
    return (
        isinstance(result, Mapping)
        and result.get("passed") is False
        and result.get("status") in expected
    )


def _diagnostic(result: Any) -> dict[str, Any]:
    """Keep a JSON-shaped, lossless-enough view of a verifier return value."""

    if isinstance(result, Mapping):
        return {
            "returned": True,
            "passed": result.get("passed"),
            "status": result.get("status"),
            "errors": result.get("errors"),
            "bridge_status": result.get("bridge_source_AST", {}).get("status")
            if isinstance(result.get("bridge_source_AST"), Mapping)
            else None,
            "core_status": result.get("core_source_AST", {}).get("status")
            if isinstance(result.get("core_source_AST"), Mapping)
            else None,
        }
    return {"returned": False, "value": result}


def _mapping(value: Any) -> Mapping[str, Any]:
    """Treat malformed nested verifier fields as empty evidence, not exceptions."""

    return value if isinstance(value, Mapping) else {}


def _run_guard(
    name: str,
    operation: Callable[[], Any],
    checks: dict[str, bool],
    diagnostics: dict[str, dict[str, Any]],
) -> Any:
    """Run one guard in isolation and preserve exception details verbatim."""

    try:
        value = operation()
    except Exception as exc:
        checks[name] = False
        diagnostics[name] = {
            "returned": False,
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
            "traceback": traceback.format_exc(),
        }
        return None
    checks[name] = bool(value)
    diagnostics[name] = _diagnostic(value) if isinstance(value, Mapping) else {
        "returned": True,
        "value": value,
    }
    return value


def run_fixture(verify: Callable[..., Mapping[str, Any]]) -> dict[str, Any]:
    """Run positive, fail-closed, alias, and nested-lineage synthetic guards."""

    o2 = _o2_packet()
    core = _core_packet()
    lineage = _lineage_packet()
    core_source = _core_source_text()
    bridge_source = _bridge_source_text()
    checks: dict[str, bool] = {}
    diagnostics: dict[str, dict[str, Any]] = {}

    accepted = _run_guard(
        "acceptance_recipe_status",
        lambda: verify(o2, core, lineage, core_source, bridge_source),
        checks,
        diagnostics,
    )
    if isinstance(accepted, Mapping):
        checks["v2_schema_identifier"] = accepted.get("schema") == LEAF_SCHEMA
        diagnostics["v2_schema_identifier"] = {
            "returned": True,
            "observed_schema": accepted.get("schema"),
            "expected_schema": LEAF_SCHEMA,
        }
        checks["recipe_pair_is_certified"] = (
            _mapping(accepted.get("proof_scope")).get("expression_recipes_certified_and_same")
            is True
        )
        diagnostics["recipe_pair_is_certified"] = _diagnostic(accepted)
        checks["pin_pair_is_provenance_only"] = (
            _mapping(accepted.get("source_pin_consistency")).get("same_source_pin") is True
            and _mapping(accepted.get("source_pin_consistency")).get("source_pin_is_provenance_only")
            is True
        )
        diagnostics["pin_pair_is_provenance_only"] = _diagnostic(accepted)
        checks["nested_leaf_status_survives"] = (
            _mapping(accepted.get("lineage_receipt")).get("proved_status_preserved") is True
            and _mapping(accepted.get("lineage_receipt")).get("path") == "result.result.leaf"
        )
        diagnostics["nested_leaf_status_survives"] = _diagnostic(accepted)
        checks["bridge_flow_is_certified"] = (
            _mapping(accepted.get("bridge_source_AST")).get("actual_data_flow_proved") is True
        )
        diagnostics["bridge_flow_is_certified"] = _diagnostic(accepted)
        checks["live_identity_remains_unresolved"] = (
            _mapping(accepted.get("proof_scope")).get("live_parameter_object_identity") is False
        )
        diagnostics["live_identity_remains_unresolved"] = _diagnostic(accepted)
        checks["legacy_scope_flags_are_retained"] = (
            accepted.get("old_212_terms_and_all_errors_unchanged") is True
            and accepted.get("old_errors_removed") is False
        )
        diagnostics["legacy_scope_flags_are_retained"] = _diagnostic(accepted)
        checks["static_execution_flags_are_false"] = (
            accepted.get("source_execution") is False
            and accepted.get("real_receipt_opened_here") is False
            and accepted.get("producer_imported") is False
            and accepted.get("graph_traversed") is False
        )
        diagnostics["static_execution_flags_are_false"] = _diagnostic(accepted)
    else:
        for name in (
            "v2_schema_identifier",
            "recipe_pair_is_certified",
            "pin_pair_is_provenance_only",
            "nested_leaf_status_survives",
            "bridge_flow_is_certified",
            "live_identity_remains_unresolved",
            "legacy_scope_flags_are_retained",
            "static_execution_flags_are_false",
        ):
            checks[name] = False
            diagnostics[name] = {
                "returned": False,
                "reason": "positive verifier call raised or returned a non-mapping",
            }

    wrong_branch = copy.deepcopy(o2)
    wrong_branch["delta_branch_proof"]["selected_branch"] = "exp(-4logPstar-29)"
    _run_guard(
        "branch_variant_closes_as_mismatch",
        lambda: _failed(
            verify(wrong_branch, core, lineage, core_source, bridge_source),
            "mismatch",
        ),
        checks,
        diagnostics,
    )

    missing_gate = copy.deepcopy(core)
    del missing_gate["checks"]["actual_oracle"][
        "exact_parameter_and_native_P0_identity"
    ]
    _run_guard(
        "gate_omission_closes_as_incomplete",
        lambda: _failed(
            verify(o2, missing_gate, lineage, core_source, bridge_source),
            "incomplete",
        ),
        checks,
        diagnostics,
    )

    tampered_core = core_source.replace("c(11)", "c(13)")
    _run_guard(
        "core_recipe_mutation_is_rejected",
        lambda: (
            lambda result: _failed(result, "mismatch")
            and result.get("core_source_AST", {}).get("expected_logP_is_exp_40_plus_11")
            is False
        )(verify(o2, core, lineage, tampered_core, bridge_source)),
        checks,
        diagnostics,
    )

    disguised_bridge = bridge_source.replace(
        "boundary_alias.parameter_bounds.log_delta",
        "boundary_alias.unmapped.log_delta",
    )
    _run_guard(
        "bridge_field_disguise_stays_unresolved",
        lambda: (
            lambda result: _failed(result, ("incomplete", "mismatch"))
            and result.get("bridge_source_AST", {}).get("status") == "unresolved"
        )(verify(o2, core, lineage, core_source, disguised_bridge)),
        checks,
        diagnostics,
    )

    _run_guard(
        "separate_packet_aliases_do_not_claim_identity",
        lambda: (
            lambda result: (
                result.get("status") == "proved"
                and result.get("proof_scope", {}).get("live_parameter_object_identity") is False
                and result.get("provenance_scope", {}).get("cross_source_live_parameter_identity")
                is False
            )
        )(verify(o2, core, lineage, core_source, _cross_object_bridge_text())),
        checks,
        diagnostics,
    )

    nested_lineage_probe = {
        "result": {"result": {"leaf": lineage["result"]["result"]["leaf"]}},
        "envelope_marker": "deep-coral-probe-v2",
    }
    _run_guard(
        "deep_lineage_path_is_retained",
        lambda: (
            lambda result: (
                result.get("status") == "proved"
                and result.get("lineage_receipt", {}).get("path") == "result.result.leaf"
            )
        )(verify(o2, core, nested_lineage_probe, core_source, bridge_source)),
        checks,
        diagnostics,
    )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(checks.values()),
        "synthetic_only": True,
        "real_sources_read": False,
        "implementation_imported": False,
        "checks": checks,
        "diagnostics": diagnostics,
        "source_execution": False,
        "receipt_decoded_here": False,
        "producer_imported": False,
        "source_imported": False,
        "graph_traversed": False,
    }


run = run_fixture


__all__ = ["SCHEMA", "run_fixture", "run"]
