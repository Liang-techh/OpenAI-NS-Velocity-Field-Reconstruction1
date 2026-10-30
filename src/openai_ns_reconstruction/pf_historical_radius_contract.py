"""Source-only contract for the unresolved historical P/F gap radius.

This module consumes only the two already-published JSON result reports for
the deferred-role and saved-rewrite leaves.  It authenticates each report's
raw bytes and each embedded result's raw bytes before parsing direct fields.
It never imports or runs an old leaf, opens a gzip archive, follows graph
edges, or interprets an opaque node id as a numerical value.

The contract records the source-defined historical template
``R_F = exp(logRp_F + xi / mu_F)`` and the selected original gap endpoint
roles.  It deliberately does not claim that the template is the rewritten
historical F radius, the current G radius, or the physical P radius.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "pf_historical_radius_contract_v1"

DEFERRED_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_historical_frame_gap_deferred_roles.json"
)
DEFERRED_REPORT_SHA256 = (
    "caca47450da722fa290d5cb8d29d845e2c747b2d81b311317966ab38a2e33e25"
)
DEFERRED_ARTIFACT = "tmp/historical_frame_gap_deferred_roles_v1/result.json"
DEFERRED_RESULT_SHA256 = (
    "3ba2e27cb66459ca75592de04bf6359fa3295733cf87aa1b02f4efcf62d30ec6"
)

LOOKUP_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_historical_frame_gap_saved_rewrite_lookup.json"
)
LOOKUP_REPORT_SHA256 = (
    "f9e6eff7457afbe591323268655be122455de98b14bf320a235c80137be75b00"
)
LOOKUP_ARTIFACT = "tmp/historical_frame_gap_saved_rewrite_lookup_v1/result.json"
LOOKUP_RESULT_SHA256 = (
    "44da94d04c52e020799692049e8a6c58f1585a1ea4f794a0146ba933c71cfc4d"
)

# These are source-definition anchors, not executable imports.  They document
# the source structure used by the contract and authenticate the bytes before
# the contract may be admitted.  No source graph is parsed or traversed.
STATIC_SOURCE_SHA256 = {
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py":
        "1400066d2da5373a95f6cd30edad19359256211d0358b454203a92fac53c5ced",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py":
        "8cc7442487dbbd3054e29211ea799a04d27c43fac3f060f74f5fadd929737594",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py":
        "3a41add18b7f9bfb63ab3c484bea82c321c514aeddd59f46a5c41ca2c9edaad2",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py":
        "369e6404d885dd4548e4f6451f2d5a5ed4160c0abb79a4aca2849520ede03541",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py":
        "38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py":
        "cb7bc1bc1801e04c198f40458dd88246434e74318f500260e8a9043c49398674",
}

STATIC_SOURCE_ANCHORS = {
    "historical_logRp_owner": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py",
        "lines": "61-73",
        "definition": "original_logRp = phase_rc + log(r_minus) - 2 - 60*log(mu)",
    },
    "gamma_imports_frame_logRp": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py",
        "lines": "212-232",
        "definition": "Gamma.schedule is built from the imported frame original_logRp",
    },
    "own_field_radius": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py",
        "lines": "44-142,156-159",
        "definition": "y=xi/mu; logR=logRp+y; R=exp(logR)",
    },
    "gap_endpoint_builder": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py",
        "lines": "27-63,66-78",
        "definition": "field R supplies chart measures and finite endpoint tensors",
    },
    "high_endpoint_owner": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py",
        "lines": "159-170",
        "definition": "C3 retains per-chart and global endpoint bindings",
    },
    "current_radius_metadata_is_separate": {
        "path": "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py",
        "lines": "84-113",
        "definition": "current pulse_coordinates radius is rebuilt from current geometry",
    },
}

ENDPOINT_FIELDS = (
    "R_q",
    "R_squared_q",
    "R_q_squared",
    "R_squared_q_squared",
    "R_U_squared",
    "R_squared_U_squared",
    "q_plus_q_y",
    "R_q_y",
)
SELECTED_ENDPOINT_FIELDS = ("q_plus_q_y", "R_q_y")


class HistoricalRadiusContractInputError(ValueError):
    """A required authenticated direct-field input is absent or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise HistoricalRadiusContractInputError(message)


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _relative(path: Path) -> str:
    return path.as_posix()


def _read_authenticated_bytes(root: Path, relative: Path, expected: str) -> bytes:
    path = root / relative
    _require(path.is_file(), f"missing required input: {_relative(relative)}")
    raw = path.read_bytes()
    _require(
        _sha256(raw) == expected,
        f"raw SHA-256 mismatch for {_relative(relative)}",
    )
    return raw


def _read_json_report(root: Path, relative: Path, expected: str) -> Mapping[str, Any]:
    raw = _read_authenticated_bytes(root, relative, expected)
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise HistoricalRadiusContractInputError(
            f"invalid JSON report: {_relative(relative)}"
        ) from exc
    _require(isinstance(value, Mapping), f"report is not an object: {_relative(relative)}")
    return value


def _read_embedded_result(
    root: Path,
    report_path: Path,
    report_sha256: str,
    artifact_key: str,
    result_sha256: str,
    expected_schema: str,
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    report = _read_json_report(root, report_path, report_sha256)
    artifacts = report.get("artifacts")
    _require(isinstance(artifacts, Mapping), f"missing artifacts in {_relative(report_path)}")
    artifact = artifacts.get(artifact_key)
    _require(isinstance(artifact, Mapping), f"missing result artifact: {artifact_key}")
    encoded = artifact.get("base64")
    _require(isinstance(encoded, str), f"result artifact has no base64: {artifact_key}")
    try:
        result_raw = base64.b64decode(encoded, validate=True)
    except (ValueError, TypeError) as exc:
        raise HistoricalRadiusContractInputError(
            f"invalid base64 result artifact: {artifact_key}"
        ) from exc
    _require(
        _sha256(result_raw) == result_sha256,
        f"embedded result SHA-256 mismatch: {artifact_key}",
    )
    try:
        packet = json.loads(result_raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise HistoricalRadiusContractInputError(
            f"invalid embedded result JSON: {artifact_key}"
        ) from exc
    _require(isinstance(packet, Mapping), f"embedded result is not an object: {artifact_key}")
    _require(packet.get("status") == "returned", f"result is not returned: {artifact_key}")
    result = packet.get("result")
    _require(isinstance(result, Mapping), f"returned result has no object: {artifact_key}")
    _require(result.get("schema") == expected_schema, f"unexpected result schema: {artifact_key}")
    return result, {
        "report_path": _relative(report_path),
        "report_raw_sha256": report_sha256,
        "result_artifact": artifact_key,
        "result_raw_sha256": result_sha256,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": False,
    }


def _authenticate_static_sources(root: Path) -> dict[str, str]:
    """Hash only the pinned source files; never parse/import them."""

    authenticated: dict[str, str] = {}
    for relative, expected in STATIC_SOURCE_SHA256.items():
        raw = _read_authenticated_bytes(root, Path(relative), expected)
        authenticated[relative] = _sha256(raw)
    return authenticated


def _opaque_ref(value: Any, source_role: str) -> dict[str, Any]:
    if type(value) is int:
        return {
            "status": "saved_opaque_ref",
            "source_role": source_role,
            "source_node_id": value,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
        }
    return {
        "status": "unknown",
        "source_role": source_role,
        "reason": "missing_or_noninteger_saved_reference",
        "source_value_selected": False,
    }


def _row_item(item: Any, source_role: str, row_order: int) -> dict[str, Any]:
    if type(item) is int:
        return dict(_opaque_ref(item, f"{source_role}[{row_order}]"), row_order=row_order)
    if isinstance(item, Mapping):
        original = item.get("original_ref", item.get("source_node_id"))
        record = dict(_opaque_ref(original, f"{source_role}[{row_order}]"))
        record["row_order"] = row_order
        if "status" in item:
            record["saved_status"] = item.get("status")
        if "rewritten_ref" in item:
            record["rewritten_ref"] = _opaque_ref(
                item.get("rewritten_ref"), f"{source_role}[{row_order}].rewritten_ref"
            )
        return record
    return {
        "status": "unknown",
        "source_role": f"{source_role}[{row_order}]",
        "row_order": row_order,
        "reason": "malformed_saved_row",
        "source_value_selected": False,
    }


def _rows(value: Any, source_role: str, expected_length: int = 4) -> dict[str, Any]:
    if not isinstance(value, list):
        return {
            "status": "unknown",
            "source_role": source_role,
            "ordinary_Z_orders": list(range(expected_length)),
            "reason": "saved_row_list_missing",
            "source_value_selected": False,
        }
    return {
        "status": "direct_saved_rows" if len(value) >= expected_length else "unknown",
        "source_role": source_role,
        "ordinary_Z_orders": list(range(expected_length)),
        "rows": [_row_item(item, source_role, index) for index, item in enumerate(value[:expected_length])],
        "row_count": len(value),
        "source_value_selected": False,
        "reason": None if len(value) >= expected_length else "fewer_than_four_saved_rows",
    }


def _original_ref(row: Mapping[str, Any]) -> int | None:
    value = row.get("source_node_id")
    return value if type(value) is int else None


def _row_originals(packet: Mapping[str, Any]) -> list[int | None]:
    rows = packet.get("rows")
    if not isinstance(rows, list):
        return []
    return [_original_ref(row) for row in rows if isinstance(row, Mapping)]


def _direct_metadata(value: Any, source_role: str) -> dict[str, Any]:
    if value is None:
        return {
            "status": "unknown",
            "source_role": source_role,
            "reason": "saved_metadata_missing",
            "source_value_selected": False,
        }
    return {
        "status": "direct_saved_metadata",
        "source_role": source_role,
        "value": value,
        "source_value_selected": False,
    }


def _cross_check_selected_rows(
    deferred: Mapping[str, Any], lookup: Mapping[str, Any]
) -> dict[str, Any]:
    endpoint_roles = deferred.get("original_gap_endpoint_roles")
    lookup_roles = lookup.get("roles")
    _require(isinstance(endpoint_roles, Mapping), "deferred endpoint roles are missing")
    _require(isinstance(lookup_roles, Mapping), "lookup roles are missing")
    checks: dict[str, Any] = {}
    for side in ("lower", "upper"):
        endpoint = endpoint_roles.get(side)
        _require(isinstance(endpoint, Mapping), f"deferred {side} endpoint is missing")
        for field in SELECTED_ENDPOINT_FIELDS:
            name = f"{side}_{field}"
            deferred_rows = _rows(endpoint.get(field), f"deferred.original_gap_endpoint_roles.{side}.{field}")
            lookup_rows = _rows(lookup_roles.get(name), f"lookup.roles.{name}")
            deferred_ids = _row_originals(deferred_rows)
            lookup_ids = _row_originals(lookup_rows)
            checks[name] = {
                "deferred_rows": deferred_rows,
                "lookup_rows": lookup_rows,
                "direct_original_refs_match": (
                    len(deferred_ids) == 4
                    and all(type(value) is int for value in deferred_ids)
                    and deferred_ids == lookup_ids
                ),
                "source_values_selected": False,
            }
    return checks


def _endpoint_tensor(
    deferred: Mapping[str, Any], selected_checks: Mapping[str, Any]
) -> dict[str, Any]:
    endpoints = deferred.get("original_gap_endpoint_roles")
    result: dict[str, Any] = {}
    for side in ("lower", "upper"):
        side_result: dict[str, Any] = {}
        source = endpoints.get(side) if isinstance(endpoints, Mapping) else None
        for field in ENDPOINT_FIELDS:
            key = f"{side}_{field}"
            if field in SELECTED_ENDPOINT_FIELDS and key in selected_checks:
                side_result[field] = selected_checks[key]
            else:
                side_result[field] = {
                    "status": "unknown",
                    "source_role": f"original_gap_endpoint_roles.{side}.{field}",
                    "reason": "full_endpoint_tensor_not_in_deferred_result",
                    "source_value_selected": False,
                }
        side_result["source_side"] = "original_gap_endpoint_roles." + side
        side_result["source_value_selected"] = False
        result[side] = side_result
    return result


def _radius_definition(
    deferred: Mapping[str, Any], lookup: Mapping[str, Any]
) -> dict[str, Any]:
    template = deferred.get("original_own_gap_field_template")
    _require(isinstance(template, Mapping), "deferred original gap field template is missing")
    lookup_roles = lookup.get("roles")
    _require(isinstance(lookup_roles, Mapping), "lookup roles are missing")
    saved_radius = lookup_roles.get("field_radius")
    radius_ref = _opaque_ref(template.get("radius"), "original_own_gap_field_template.radius")
    rewritten = {
        "status": "unknown",
        "source_role": "lookup.roles.field_radius",
        "source_value_selected": False,
    }
    if isinstance(saved_radius, Mapping):
        rewritten = dict(saved_radius)
        rewritten["source_value_selected"] = False
        rewritten["node_id_is_not_a_numeric_value"] = True
    return {
        "formula": "candidate_R_F(Z,xi)=exp(logRp_F(Z)+xi/mu_F(Z))",
        "source_template": {
            "status": "direct_saved_template_role",
            "owner": "original_own_gap_field_template",
            "radius": radius_ref,
            "query": _opaque_ref(template.get("query"), "original_own_gap_field_template.query"),
            "log_radius": _opaque_ref(
                template.get("log_radius"), "original_own_gap_field_template.log_radius"
            ),
            "domain": _direct_metadata(
                template.get("domain"), "original_own_gap_field_template.domain"
            ),
            "ordinary_Z_rows_are_source_refs": template.get(
                "ordinary_Z_rows_are_source_refs"
            ),
            "source_value_selected": False,
        },
        "saved_rewrite_boundary": rewritten,
        "historical_F_owner_bound": False,
        "historical_F_function_identity_proved": False,
        "old_template_is_current_G_owner": False,
        "old_template_is_F_rewritten_owner": False,
        "owner_separation_only": True,
        "function_inequality_claimed": False,
        "radius_Z_formula": (
            "R_F_Z=R_F*(logRp_F_Z+xi_Z/mu_F-"
            "xi*mu_F_Z/mu_F^2)"
        ),
        "radius_Z_proved": False,
        "mu_Z": {
            "status": "unknown",
            "source_role": "historical_F.mu_Z",
            "source_value_selected": False,
        },
        "logRp_Z": {
            "status": "unknown",
            "source_role": "historical_F.logRp_Z",
            "source_value_selected": False,
        },
        "moving_endpoint": {
            "lower": {
                "coordinate_formula": "xi_lower=11",
                "coordinate_Z_formula": "xi_lower_Z (source-owned; not assumed zero)",
                "source_value_selected": False,
            },
            "upper": {
                "coordinate_formula": "xi_upper=13-4*mu_F",
                "coordinate_Z_formula": "xi_upper_Z=-4*mu_F_Z",
                "log_radius_formula": "logRp_F+13/mu_F-4",
                "log_radius_Z_formula": "logRp_F_Z-13*mu_F_Z/mu_F^2",
                "moving_endpoint_chain_rule_proved": False,
                "source_value_selected": False,
            },
        },
    }


def _negative_fixture() -> dict[str, Any]:
    """Return unexecuted negative controls for callers/tests to inspect."""

    return {
        "schema": "pf_historical_radius_contract_negative_fixture_v1",
        "executed": False,
        "cases": {
            "same_opaque_radius_id_does_not_prove_function_equality": {
                "old_ref": 33982,
                "rewritten_ref": 33982,
                "must_be_rejected_as_F_binding": True,
            },
            "missing_endpoint_field_stays_unknown": {
                "field": "R_squared_q",
                "must_not_be_synthesized": True,
            },
            "caller_zero_mu_Z_is_not_source_proof": {
                "mu_Z": 0,
                "must_not_admit_radius_Z_independence": True,
            },
            "current_G_radius_is_not_F_radius": {
                "owner": "G",
                "must_not_be_used_for_F": True,
            },
        },
    }


def build_historical_radius_contract(root: str | Path = ".") -> dict[str, Any]:
    """Build an unevaluated historical F radius/endpoint source contract."""

    repo_root = Path(root)
    deferred, deferred_input = _read_embedded_result(
        repo_root,
        DEFERRED_REPORT,
        DEFERRED_REPORT_SHA256,
        DEFERRED_ARTIFACT,
        DEFERRED_RESULT_SHA256,
        "historical_frame_gap_deferred_roles_v1",
    )
    lookup, lookup_input = _read_embedded_result(
        repo_root,
        LOOKUP_REPORT,
        LOOKUP_REPORT_SHA256,
        LOOKUP_ARTIFACT,
        LOOKUP_RESULT_SHA256,
        "historical_frame_gap_saved_rewrite_lookup_v1",
    )
    _require(
        deferred.get("status") == "deferred_source_roles_not_numeric",
        "deferred result is not the expected source-only status",
    )
    _require(
        lookup.get("status") == "saved_rewrite_map_lookups_not_numerical_values",
        "lookup result is not the expected direct-map status",
    )
    selected = _cross_check_selected_rows(deferred, lookup)
    radius = _radius_definition(deferred, lookup)
    tensor = _endpoint_tensor(deferred, selected)
    template = deferred.get("original_own_gap_field_template")
    endpoint_roles = deferred.get("original_gap_endpoint_roles")
    inherited = {
        "independent_reference_offset_preserved": deferred.get(
            "independent_reference_offset_preserved", "unknown"
        ),
        "physical_correction_roots": deferred.get(
            "actual_physical_correction_roots_preserved", "unknown"
        ),
        "old_errors_removed": False,
        "error_transport_certified": False,
    }
    result = {
        "schema": SCHEMA,
        "status": "historical_F_radius_contract_source_only_unbound",
        "source_only": True,
        "source_values_computed": False,
        "graph_traversed": False,
        "old_leaf_called": False,
        "gzip_decoded": False,
        "source_static_anchors": STATIC_SOURCE_ANCHORS,
        "source_static_sha256": _authenticate_static_sources(repo_root),
        "required_inputs": {
            "deferred_roles": deferred_input,
            "saved_rewrite_lookup": lookup_input,
        },
        "historical_F_radius": radius,
        "original_gap_endpoint_roles": {
            "source_archive": (
                endpoint_roles.get("archive")
                if isinstance(endpoint_roles, Mapping)
                else "unknown"
            ),
            "source_pointer": (
                endpoint_roles.get("pointer")
                if isinstance(endpoint_roles, Mapping)
                else "unknown"
            ),
            "gap_domain": _direct_metadata(
                template.get("domain") if isinstance(template, Mapping) else None,
                "original_own_gap_field_template.domain",
            ),
            "lower": tensor["lower"],
            "upper": tensor["upper"],
            "eight_field_names": list(ENDPOINT_FIELDS),
            "selected_fields": list(SELECTED_ENDPOINT_FIELDS),
            "source_value_selected": False,
        },
        "source_gates": {
            "historical_F_radius_owner_bound": False,
            "historical_F_endpoint_owner_bound": False,
            "full_eight_field_tensor_bound": False,
            "capture_free_inlet_substitution": False,
            "ordinary_Z0_to_Z3_transport": False,
            "common_mu_delta_Pstar": False,
            "common_coordinate_measure_orientation": False,
            "radius_Z_independence": False,
            "moving_endpoint_chain_rule": False,
            "same_source_endpoint_seams": False,
            "local_to_global_m15_identity": False,
            "P_equals_F": False,
            "G_equals_F": False,
            "P_equals_G": False,
            "physical_scale_and_error_transport": False,
            "numeric_use_admissible": False,
        },
        "owner_distinctions": {
            "old_template_vs_current_G": {
                "status": "owner_distinct_not_function_inequality",
                "old_template_is_current_G": False,
                "function_inequality_proved": False,
            },
            "old_template_vs_historical_F": {
                "status": "template_requires_rewrite_and_source_binding",
                "old_template_is_F_rewritten_owner": False,
                "function_identity_proved": False,
            },
            "current_G_radius": {
                "status": "not_consumed",
                "owner": "G.current_geometry_radius",
                "usable_as_F": False,
                "G_equals_F_proved": False,
            },
            "physical_P_radius": {
                "status": "not_consumed",
                "owner": "P.current_physical_source",
                "usable_as_F": False,
                "P_equals_F_proved": False,
            },
        },
        "P_minus_F": {
            "formal_decomposition": "(P-G)+(G-F)",
            "P_minus_G_evaluated": False,
            "G_minus_F_evaluated": False,
            "either_term_zero_proved": False,
        },
        "inherited_error_gates": inherited,
        "negative_fixture": _negative_fixture(),
        "admissible_for_numeric_use": False,
        "full_goal_complete": False,
    }
    return result


def recover(root: str | Path = ".") -> dict[str, Any]:
    """Alias for the source-only contract builder."""

    return build_historical_radius_contract(root)


def synthetic_negative_fixture() -> dict[str, Any]:
    """Expose, but do not execute, the small negative fixture."""

    return _negative_fixture()


__all__ = [
    "DEFERRED_ARTIFACT",
    "DEFERRED_REPORT",
    "DEFERRED_REPORT_SHA256",
    "DEFERRED_RESULT_SHA256",
    "ENDPOINT_FIELDS",
    "HistoricalRadiusContractInputError",
    "LOOKUP_ARTIFACT",
    "LOOKUP_REPORT",
    "LOOKUP_REPORT_SHA256",
    "LOOKUP_RESULT_SHA256",
    "SCHEMA",
    "SELECTED_ENDPOINT_FIELDS",
    "STATIC_SOURCE_ANCHORS",
    "STATIC_SOURCE_SHA256",
    "build_historical_radius_contract",
    "recover",
    "synthetic_negative_fixture",
]
