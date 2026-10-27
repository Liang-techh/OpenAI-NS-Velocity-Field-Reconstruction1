"""Metadata-only lifecycle workflow for the native W parameter interface.

This manager freezes source provenance and coordinates one future plain-dict
``derive(source, diagnostic=None)`` call. It contains no W mathematics,
provider, checker, reducer, normalizer, numerical engine, or ancestor import.
"""

from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

REPO = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
WORK = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")
COMPONENT = "native_W_parameter_interface"
INPUT_SCHEMA = "native_W_parameter_interface_input_v1"
RESULT_SCHEMA = "original_native_W_parameter_interface"
RESULT_SCHEMA_ALIASES = (RESULT_SCHEMA,)
ENTERED_SCHEMA = "native_W_parameter_interface_entered_v1"
DIAGNOSTIC_SCHEMA = "native_W_parameter_interface_diagnostic_v1"
PUBLICATION_SCHEMA = "original_native_W_parameter_interface_publication_v1"
SEAL_SCHEMA = "native_W_parameter_interface_final_seal_v1"

KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_W_parameter_interface.py")
RUNNER_RELATIVE = Path("scripts/run_native_W_parameter_interface.py")
REPORT_RELATIVE = Path("experiments/root_st073/lei_ren_part1_paper_compliant_current_native_W_parameter_interface.json")
MILESTONE_RELATIVE = Path("docs/CURRENT_NATIVE_W_PARAMETER_INTERFACE_2026_10_27.md")
HANDOFF_RELATIVES = (Path("docs/AGENT_TASKS.md"), Path("docs/CROSS_DEVICE_HANDOFF.md"), Path("docs/CURRENT_CHECKPOINT.md"), Path("docs/NEXT_TASKS.md"))
SCOPED_RELATIVES = (Path(".gitattributes"), KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, MILESTONE_RELATIVE, *HANDOFF_RELATIVES)

CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (SOURCE_CHECKPOINT, ENTERED_CHECKPOINT, RESULT_CHECKPOINT, DIAGNOSTIC_CHECKPOINT)
INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_27.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_27.json"

PARENT_REPORT_DEFAULT = REPO / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_O3_transition_W_sector_correlation.json"
PARENT_REPORT_SHA256 = "c6542479047de497772d31459ec3d7893fa34dfc469626dff4590cb782ffe98d"
PARENT_RESULT_DEFAULT = REPO / "tmp/compact_O3_transition_W_sector_correlation/result.json"
PARENT_RESULT_SHA256 = "b2a67a90c45417f3e5ca3f6aa8d2a6cad0391a39347d6fccb7d0ba336f20bf84"
ACCEPTED_BASELINE_REPORT = REPO / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_O3_transition_P2_pressure_gap.json"
ACCEPTED_BASELINE_REPORT_SHA256 = "2956baf8c82250ffd789714ceb1ee80170bd69e57f3c84d1d6b3b8889a1d1466"
PARENT_COMMIT = "4707fc262d51d33677dc5d53fb6f1bb1e6fb53af"
PARENT_PARENT_COMMIT = "fec736e2a0ae6b55a59f7ab22f0e2037d7a5e4f3"
CURRENT_MAIN_BRIEF_SHA256 = "e1bfcaa334513e6869b53d6c4b4eeb2e8872b527cfb36183d4ce4f13ca18916b"

PARENT_SOURCE_PIN_COUNT = 2247
PARENT_CHECKPOINT_COUNT = 256
PARENT_FAILED_COUNT = 183
FAILED_COUNT_IF_NEW_FAILURE = 190
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
SOURCE_FACTS_DEFAULT = WORK / "wsector_source_trace_2026_10_27_v2.json"
SOURCE_FACTS_SHA256 = "eb0af8bce8ab9af4bb72f8bc8a3eb5e03f5a2a6db9688ab30f37104084003dcf"
SOURCE_FACTS_V1_SHA256 = "7d11e08d8b52ef2eab5ff88b7b49e8087cd296b47c59a77d933d662b36475e56"
PARAMETER_ALIAS_TRACE_DEFAULT = WORK / "wsector_parameter_alias_trace_2026_10_27.json"
PARAMETER_ALIAS_TRACE_SHA256 = "b483e1bd0657fc2ad223f4153044875191b5beae7e2028e99f7a08b9d05ba112"
PARAMETER_ALIAS_TRACE_BYTES = 23979
DESIGN_REVIEW_DEFAULT = WORK / "wsector_design_review_2026_10_27_v2.json"
DESIGN_REVIEW_SHA256 = "a7797c06cebf04176d3c4cc5e714e87d9f0669c9f278fd3043cd951c24c90d4c"
PARENT_REVIEW_DEFAULT = WORK / "review_compact_O3_transition_P2_pressure_gap_saved_2026_10_27_v2_receipt.json"
PARENT_REVIEW_SHA256 = "31dfc0eb524bcc767a4d86766c719aeddc840497b1b9b32f50a64b2855ad876d"
KERNEL_DRAFT_DEFAULT = WORK / "native_W_parameter_interface_draft_2026_10_27_v3.py"
KERNEL_DRAFT_SHA256 = "683c9ff5f279acd52b3052bba62c4cf2001c83680c713c6343efd25cfa6d536c"
KERNEL_API_DEFAULT = WORK / "native_W_parameter_interface_api_2026_10_27_v3.json"
KERNEL_API_SHA256 = "ab9823f7f50aa747a3842559438b25fa189f824eee551069c942c8591906a250"
KERNEL_DRAFT_AST_SHA256 = "c0702cfb72e6f4f35094d7cac7d10d1b3036212120944c227db9f3561eae5ba8"
SOURCE_INTERFACE_INPUT_PARENT_SHA256 = "958ba5154ea7a532f33d2e48cafeb12209bf94e376b9a0cc5dafae83bfb676a7"

EXPECTED_COUNTS = {"source_binding_calls": 1, "native_admissions": 2, "rejections": 2, "template_signatures": 1, "new_mathematical_proofs": 0, "new_algebra_family_calls": 0, "numeric_calls": 0, "source_provider_calls": 0, "old_math_calls": 0}
COUNT_KEYS = tuple(EXPECTED_COUNTS)
EXPECTED_FLAGS = {"native_W_parameter_source_handles_bound": True, "native_W_parameter_enclosure_parameters_separated": True, "native_W_parameter_uninstantiated_template_signature_defined": True}
FLAG_KEYS = tuple(EXPECTED_FLAGS)
CALLBACK_EVENT_COUNTS = {
    "source_validation": 1,
    "native_handle_admission": 2,
    "source_reference_rejection": 2,
    "template_signature": 1,
}
CALLBACK_EVENT_COUNT = sum(CALLBACK_EVENT_COUNTS.values())
CALLBACK_EVENT_COUNT_IF_NONE = 0
RETAINED_MAP_NAMES = ("preserved_error_cap_metadata", "retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance")
KERNEL_CONTRACT = {
    "derive": "derive(source, diagnostic=None)",
    "component": COMPONENT,
    "input_schema": INPUT_SCHEMA,
    "result_schema": RESULT_SCHEMA,
    "metadata_only": True,
    "source_only": True,
    "plain_json_result": True,
    "callback_rows": dict(CALLBACK_EVENT_COUNTS),
    "diagnostic_callback_calls": CALLBACK_EVENT_COUNT,
    "payload_arrays_are_not_counters": True,
    "no_math_or_provider_calls": True,
}
NONCLAIM_FLAGS = (
    "full_goal_complete", "physical_layer_completed", "numeric_values_computed",
    "numeric_containment_claimed", "global_closure_proved", "time_order_recursion_claimed",
    "temporal_n_recursion_claimed", "pointwise_parameter_identity_proved",
    "radius_equality_proved", "phase_equality_proved", "cap_error_reduction_claimed",
    "error_cap_reduced", "old_error_terms_removed", "old_errors_removed",
)
REQUIRED_BRIEF_KEYS = (
    "component", "current_component", "execution_rules", "dirty_paths", "protected_paths",
    "authority", "next_task", "read_first", "current_exact_source_result",
    "current_document_sha256", "latest_review", "latest_QA", "latest_seal",
)


def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    return sha_bytes(Path(path).read_bytes())


def encode(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8") + b"\n"


def clone(value: Any) -> Any:
    return json.loads(encode(value).decode("utf-8"))


def brief_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: brief_safe(child) for key, child in value.items() if key != "raw_bytes_base64"}
    if isinstance(value, list):
        return [brief_safe(child) for child in value]
    return value


def portable(value: Any) -> Any:
    if isinstance(value, dict):
        return {key.replace("\\", "/") if isinstance(key, str) else key: portable(child) for key, child in value.items()}
    if isinstance(value, list):
        return [portable(child) for child in value]
    if isinstance(value, str) and ("\\" in value or (len(value) > 2 and value[1] == ":")):
        return value.replace("\\", "/")
    return value


def wire_path(path: Path | str) -> str:
    return Path(path).as_posix()


def write_new(path: Path, value: Any) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_replace(path: Path, value: Any) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("wb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def read_bytes(path: Path, expected: str | None = None) -> bytes:
    raw = Path(path).read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), f"SHA-256 changed: {path}")
    return raw


def read_json(path: Path, expected: str | None = None) -> tuple[bytes, dict]:
    raw = read_bytes(path, expected)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"JSON object required: {path}")
    return raw, value


def valid_sha(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def valid_revision(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 40:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def resolve_path(value: str | Path, base: Path = WORK) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    candidates = (base / candidate, REPO / candidate, REPO / "experiments/root_st073" / candidate, WORK / candidate)
    for item in candidates:
        if item.exists():
            return item.resolve()
    return candidates[0].resolve()


def auth_map(value: Any, label: str, expected_count: int | None = None) -> dict:
    require(isinstance(value, dict), f"{label} must be an object")
    result = {str(key): str(digest).lower() for key, digest in value.items()}
    require(all(valid_sha(digest) for digest in result.values()), f"{label} has malformed SHA")
    if expected_count is not None:
        require(len(result) == expected_count, f"{label} count changed")
    return result


def merge_hashes(*maps: dict) -> dict:
    result = {}
    for mapping in maps:
        for path, digest in mapping.items():
            require(result.get(path) in (None, digest), f"Hash collision changed: {path}")
            result[path] = digest
    return result


def artifact_entry(path: Path, digest: str | None = None, embed: bool = False, outside_science_pins: bool = False) -> dict:
    path = Path(path)
    raw = read_bytes(path, digest)
    entry = {"path": wire_path(path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_bytes_exact": True, "raw_readonly": True}
    if outside_science_pins:
        entry["outside_science_pins"] = True
    if embed:
        entry["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
    return entry


def bounded_parent_view(report: dict, report_sha: str = PARENT_REPORT_SHA256) -> dict:
    keys = ("schema", "component", "attempt_status", "passed", "metadata_only", "source_only", "parent_commit", "parent_head", "source_pin_count", "parent_checkpoint_count", "checkpoint_count", "preserved_failed_artifact_count", "parent_protected_file_count", "parent_dirty_path_count", "expected_parent_false_flag_count", "actual_counts", "new_flags", "full_goal_complete", "parent_false_flags", "parent_false_flag_values", "accepted_progress_physical_baseline_preserved")
    view = portable({key: report[key] for key in keys if key in report})
    view["parent_report_is_authenticated_view"] = True
    view["complete_parent_report_raw_sha256"] = report_sha
    view["omitted_fields_saved_in_complete_parent_report"] = sorted(set(report).difference(keys))
    require(len(encode(view)) < 100_000, "Bounded parent view exceeds 100 KB")
    return view


def _check_hash_map_paths(mapping: dict, base: Path = REPO) -> None:
    for name, digest in mapping.items():
        read_bytes(resolve_path(name, base), digest)


def load_accepted_baseline() -> dict:
    raw, report = read_json(ACCEPTED_BASELINE_REPORT, ACCEPTED_BASELINE_REPORT_SHA256)
    require(report.get("schema") == "original_compact_O3_transition_P2_pressure_gap_publication_v1", "Accepted baseline schema changed")
    require(report.get("component") == "compact_O3_transition_P2_pressure_gap", "Accepted baseline component changed")
    require(report.get("passed") is True and report.get("metadata_only") is True and report.get("source_only") is True, "Accepted baseline is not passed metadata-only")
    require(report.get("parent_commit") == "b7687e8f187ab5e54552d9a3d28201d69a38476c" and report.get("parent_head") == "b7687e8f187ab5e54552d9a3d28201d69a38476c", "Accepted baseline parent commit changed")
    require(report.get("source_pin_count") == 2234 and report.get("parent_checkpoint_count") == 248 and report.get("checkpoint_count") == 252, "Accepted baseline counts changed")
    require(report.get("preserved_failed_artifact_count") == 176 and report.get("parent_protected_file_count") == 16 and report.get("parent_dirty_path_count") == 2, "Accepted baseline ledger counts changed")
    source_pins = auth_map(report.get("source_pins"), "Accepted baseline source pins", 2234)
    parent_checkpoints = auth_map(report.get("parent_checkpoint_hashes"), "Accepted baseline parent checkpoints", 248)
    current_checkpoints = auth_map(report.get("checkpoint_hashes"), "Accepted baseline current checkpoints", 4)
    _check_hash_map_paths(merge_hashes(parent_checkpoints, current_checkpoints))
    protected = auth_map(report.get("protected_file_hashes"), "Accepted baseline protected files", 16)
    dirty = auth_map(report.get("legacy_dirty_path_hashes"), "Accepted baseline dirty paths", 2)
    failed = auth_map(report.get("preserved_failed_artifact_hashes_union"), "Accepted baseline failed artifacts", 176)
    _check_hash_map_paths(source_pins)
    parent_flags = report.get("parent_false_flags")
    parent_values = report.get("parent_false_flag_values")
    require(isinstance(parent_flags, list) and len(parent_flags) == PARENT_FALSE_FLAG_COUNT and len(set(parent_flags)) == PARENT_FALSE_FLAG_COUNT, "Accepted baseline false-flag ledger changed")
    require(isinstance(parent_values, dict) and parent_values == {key: False for key in parent_flags}, "Accepted baseline false values changed")
    for name in RETAINED_MAP_NAMES:
        require(isinstance(report.get(name), dict), f"Accepted retained map missing: {name}")
    result_path = resolve_path(report.get("result_checkpoint", {}).get("path", ""), REPO)
    result_raw, result = read_json(result_path, report.get("full_result_raw_sha256"))
    require(report.get("full_result_raw_sha256") == SOURCE_INTERFACE_INPUT_PARENT_SHA256, "Accepted baseline result SHA changed")
    require(report.get("full_result") == result and result.get("passed") is True, "Accepted baseline full result changed")
    return {
        "report": report,
        "report_raw": raw,
        "report_meta": {"path": wire_path(ACCEPTED_BASELINE_REPORT), "sha256": ACCEPTED_BASELINE_REPORT_SHA256, "bytes": len(raw), "schema": report.get("schema"), "component": report.get("component")},
        "report_view": bounded_parent_view(report, ACCEPTED_BASELINE_REPORT_SHA256),
        "result": result,
        "result_raw": result_raw,
        "result_meta": {"path": wire_path(result_path), "sha256": SOURCE_INTERFACE_INPUT_PARENT_SHA256, "bytes": len(result_raw), "schema": result.get("schema"), "component": result.get("component")},
        "source_pins": source_pins,
        "parent_checkpoints": parent_checkpoints,
        "current_checkpoints": current_checkpoints,
        "protected": protected,
        "dirty": dirty,
        "failed": failed,
        "false_flags": list(parent_flags),
        "false_flag_values": clone(parent_values),
        "retained": {name: clone(report[name]) for name in RETAINED_MAP_NAMES},
    }


def load_parent(args: argparse.Namespace | None = None) -> dict:
    args = args or argparse.Namespace()
    accepted = load_accepted_baseline()
    report_path = resolve_path(getattr(args, "parent_report", None) or PARENT_REPORT_DEFAULT, WORK)
    report_sha = (getattr(args, "parent_report_sha256", None) or PARENT_REPORT_SHA256).lower()
    require(report_sha == PARENT_REPORT_SHA256, "Failed parent report SHA is not authenticated")
    raw, report = read_json(report_path, report_sha)
    require(report.get("schema") == "original_compact_O3_transition_W_sector_correlation_publication_v1", "Failed parent report schema changed")
    require(report.get("component") == "compact_O3_transition_W_sector_correlation", "Failed parent report component changed")
    require(report.get("passed") is False and report.get("attempt_status") == "failed", "Failed parent status changed")
    require(report.get("metadata_only") is True and report.get("source_only") is True, "Failed parent publication markers changed")
    require(report.get("parent_commit") == PARENT_PARENT_COMMIT and report.get("parent_head") == PARENT_PARENT_COMMIT, "Failed parent commit anchor changed")
    require(report.get("parent_report_sha256") == ACCEPTED_BASELINE_REPORT_SHA256, "Failed parent accepted baseline pointer changed")
    require(report.get("parent_source_pin_count") == 2234 and report.get("parent_checkpoint_count") == 252 and report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Failed parent source/checkpoint counts changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT and report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Failed parent checkpoint/failure counts changed")
    require(report.get("parent_protected_file_count") == PARENT_PROTECTED_COUNT and report.get("parent_dirty_path_count") == PARENT_DIRTY_COUNT, "Failed parent protected/dirty counts changed")
    source_pins = auth_map(report.get("source_pins"), "Failed parent source pins", PARENT_SOURCE_PIN_COUNT)
    current_checkpoints = auth_map(report.get("checkpoint_hashes"), "Failed parent current checkpoints", 4)
    _check_hash_map_paths(source_pins)
    _check_hash_map_paths(current_checkpoints)
    protected = auth_map(report.get("protected_file_hashes"), "Failed parent protected files", PARENT_PROTECTED_COUNT)
    dirty = auth_map(report.get("legacy_dirty_path_hashes"), "Failed parent dirty paths", PARENT_DIRTY_COUNT)
    failed = auth_map(report.get("preserved_failed_artifact_hashes_union"), "Failed parent artifacts", PARENT_FAILED_COUNT)
    parent_flags = report.get("parent_false_flags")
    parent_values = report.get("parent_false_flag_values")
    require(isinstance(parent_flags, list) and len(parent_flags) == PARENT_FALSE_FLAG_COUNT and len(set(parent_flags)) == PARENT_FALSE_FLAG_COUNT, "Failed parent false-flag ledger changed")
    require(isinstance(parent_values, dict) and parent_values == {key: False for key in parent_flags}, "Failed parent false values changed")
    retained = {}
    for name in RETAINED_MAP_NAMES:
        require(isinstance(report.get(name), dict), f"Failed parent retained map missing: {name}")
        retained[name] = clone(report[name])
    result_path = resolve_path(PARENT_RESULT_DEFAULT, REPO)
    result_raw, result = read_json(result_path, PARENT_RESULT_SHA256)
    require(report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256 and report.get("full_result") == result, "Failed parent result is not exact")
    require(result.get("schema") == "original_compact_O3_transition_W_sector_correlation" and result.get("component") == "compact_O3_transition_W_sector_correlation", "Failed parent result schema changed")
    require(result.get("passed") is False and result.get("metadata_only") is False and result.get("source_only") is True, "Failed parent result markers changed")
    require(report.get("parent_result_raw_sha256") == SOURCE_INTERFACE_INPUT_PARENT_SHA256, "Failed parent accepted-result pointer changed")
    return {
        "report": report,
        "report_raw": raw,
        "report_meta": {"path": wire_path(report_path), "sha256": report_sha, "bytes": len(raw), "schema": report.get("schema"), "component": report.get("component")},
        "report_view": bounded_parent_view(report, PARENT_REPORT_SHA256),
        "result": result,
        "result_raw": result_raw,
        "result_meta": {"path": wire_path(result_path), "sha256": PARENT_RESULT_SHA256, "bytes": len(result_raw), "schema": result.get("schema"), "component": result.get("component")},
        "source_pins": source_pins,
        "parent_checkpoints": accepted["parent_checkpoints"],
        "current_checkpoints": current_checkpoints,
        "protected": protected,
        "dirty": dirty,
        "failed": failed,
        "false_flags": list(parent_flags),
        "false_flag_values": clone(parent_values),
        "retained": retained,
        "accepted": accepted,
    }


def load_source_facts(path: Path, expected_sha: str) -> tuple[bytes, dict]:
    raw, facts = read_json(path, expected_sha)
    require(facts.get("packet_type") == "wsector_source_trace_v2", "Wrong source trace packet type")
    require(facts.get("status") == "static_source_ast_trace_only", "Source trace is not static-only")
    require(facts.get("v1_artifact", {}).get("sha256") == SOURCE_FACTS_V1_SHA256, "V1 source trace SHA changed")
    raw_inputs = facts.get("raw_inputs")
    require(isinstance(raw_inputs, dict), "Source trace raw-input metadata is missing")
    require(raw_inputs.get("parent_report", {}).get("sha256") == ACCEPTED_BASELINE_REPORT_SHA256, "Source trace parent report pin changed")
    require(raw_inputs.get("saved_result", {}).get("sha256") == SOURCE_INTERFACE_INPUT_PARENT_SHA256, "Source trace saved result pin changed")
    scope = facts.get("scope")
    v2_scope = facts.get("v2_scope")
    require(isinstance(scope, dict) and scope.get("source_bytes_and_AST_only") is True, "Source trace is not bytes/AST-only")
    require(scope.get("project_modules_executed") is False and scope.get("project_modules_imported") is False, "Source trace imported or executed project code")
    require(scope.get("math_helper_provider_checker_normalizer_reducer_quadrature_recurrence_calls") == 0, "Source trace contains consumer calls")
    require(isinstance(v2_scope, dict) and v2_scope.get("new_source_consumer_is_metadata_only") is True, "V2 source consumer is not metadata-only")
    require(v2_scope.get("project_modules_executed") is False and v2_scope.get("project_modules_imported") is False, "V2 source trace imported or executed project code")
    require(v2_scope.get("formula_equivalence_checked") is False and v2_scope.get("ordinary_Z_equivalence_evaluated") is False, "V2 trace performed mathematical comparison")
    anchors = facts.get("v2_source_anchors")
    require(isinstance(anchors, dict), "V2 source anchors missing")
    coordinate = anchors.get("native_added_transport_and_coordinate_bindings", {}).get("coordinate_Z_binding", {})
    require(coordinate.get("targets") == ["z"] and coordinate.get("line_span") == [35, 35], "Source coordinate anchor changed")
    functions = facts.get("source_functions")
    require(isinstance(functions, dict), "Source function IR missing")
    required = {"native_generic_functions", "native_generic_roots_handoff", "native_outer_call_route", "native_outer_transport_geometry", "native_quotient_cell_port", "rebuilt_C3_recover_inputs", "rebuilt_C3_source_quotients_port", "rebuilt_compile_recovery_ast_anchor"}
    require(required.issubset(functions), "Source function IR lacks required native or rebuilt port")
    for key, value in functions.items():
        require(isinstance(value, dict) and valid_sha(value.get("source_sha256")), f"Source function hash missing: {key}")
    return raw, facts


def load_alias_trace(path: Path = PARAMETER_ALIAS_TRACE_DEFAULT) -> tuple[bytes, dict, dict]:
    raw, trace = read_json(path, PARAMETER_ALIAS_TRACE_SHA256)
    require(len(raw) == PARAMETER_ALIAS_TRACE_BYTES, "Parameter alias trace byte count changed")
    require(trace.get("schema") == "wsector_parameter_alias_trace_static_v1", "Wrong parameter alias trace schema")
    require(trace.get("runtime_model_metadata") is None, "Parameter alias trace has runtime metadata")
    audit = trace.get("static_audit")
    require(isinstance(audit, dict), "Parameter alias static audit missing")
    require(audit.get("project_modules_imported") is False and audit.get("project_modules_called") is False, "Parameter alias trace imported or called project code")
    require(audit.get("math_or_helper_consumer_calls") == 0 and audit.get("numeric_values_computed") is False, "Parameter alias trace contains consumer work")
    status = trace.get("alias_and_theorem_status")
    require(isinstance(status, dict), "Parameter alias status missing")
    for key in ("native_to_rebuilt_same_source_object_Pstar", "native_to_rebuilt_same_source_object_delta", "native_to_rebuilt_same_source_object_radius", "explicit_saved_alias_theorem_for_delta_S", "explicit_saved_alias_theorem_for_delta_delta"):
        require(status.get(key) is False, f"Parameter alias status changed: {key}")
    require(trace.get("next_action_readiness", {}).get("no_math_identity_claim") is True, "Parameter alias trace makes an identity claim")
    pointer = artifact_entry(path, PARAMETER_ALIAS_TRACE_SHA256, embed=False, outside_science_pins=True)
    pointer.update({"static_only": True, "no_cross_side_parameter_identity_claim": True, "no_zero_or_error_cap_claim": True})
    return raw, trace, pointer


def load_design_review(path: Path = DESIGN_REVIEW_DEFAULT) -> tuple[bytes, dict]:
    raw, review = read_json(path, DESIGN_REVIEW_SHA256)
    require(review.get("schema") == "wsector_design_review.v2_addendum", "Wrong design review schema")
    require(review.get("status") == "static_advisory_only_no_execution", "Design review is not advisory-only")
    state = review.get("receipt_state", {})
    require(state.get("kernel_executed") is False and state.get("proof_executed") is False, "Design review records execution")
    chart = review.get("held_fixed_chart_contract", {})
    require(chart.get("held_fixed_scalars") == ["Rn", "Rc", "Sn", "Sc", "dn", "dc"], "Design review held-fixed ports changed")
    return raw, review


def load_review(path: Path = PARENT_REVIEW_DEFAULT) -> tuple[bytes, dict]:
    raw, review = read_json(path, PARENT_REVIEW_SHA256)
    require(review.get("schema") == "compact_O3_transition_P2_pressure_gap_saved_review_v2", "Wrong saved review schema")
    require(review.get("review_kind") == "saved_metadata_only" and review.get("independent_review") is True, "Saved review is not independent metadata")
    require(review.get("passed") is True and review.get("source_only") is True, "Saved review did not pass source-only checks")
    meta = {"path": wire_path(path), "sha256": sha_bytes(raw), "bytes": len(raw), "schema": review.get("schema"), "review_kind": review.get("review_kind"), "independent_review": True, "outside_science_pins": True}
    return raw, meta


def load_kernel_api(path: Path = KERNEL_API_DEFAULT) -> tuple[bytes, dict, dict]:
    raw, api = read_json(path, KERNEL_API_SHA256)
    require(api.get("schema") == "native_W_parameter_interface_api_v3", "Wrong kernel API metadata schema")
    require(api.get("component") == COMPONENT, "Kernel API component changed")
    input_contract = api.get("input_contract", {})
    require(input_contract.get("schema") == INPUT_SCHEMA and input_contract.get("signature") == "derive(source: Mapping, diagnostic: Callable | None = None) -> dict", "Kernel API input contract changed")
    require(input_contract.get("alias_trace", {}).get("raw_sha256") == PARAMETER_ALIAS_TRACE_SHA256, "Kernel API alias pin changed")
    require(input_contract.get("source_trace", {}).get("raw_sha256") == SOURCE_FACTS_SHA256, "Kernel API source pin changed")
    result_contract = api.get("result_contract", {})
    require(result_contract.get("schema") == RESULT_SCHEMA and result_contract.get("component") == COMPONENT and result_contract.get("metadata_only_on_failure") is True and result_contract.get("source_only_on_failure") is True, "Kernel API result contract changed")
    require(result_contract.get("new_flags") == list(FLAG_KEYS), "Kernel API flag contract changed")
    require(api.get("actual_count_keys") == list(COUNT_KEYS), "Kernel API count keys changed")
    require(api.get("success_count_contract") == EXPECTED_COUNTS, "Kernel API count values changed")
    diagnostic = input_contract.get("diagnostic", {})
    require(diagnostic.get("success_callback_event_count") == CALLBACK_EVENT_COUNT, "Kernel API callback count changed")
    meta = {"path": wire_path(path), "sha256": sha_bytes(raw), "bytes": len(raw), "schema": api.get("schema"), "status": api.get("artifact", {}).get("status"), "outside_science_pins": True}
    return raw, api, meta


def parent_data_packet(parent: dict) -> dict:
    return {
        "failed_parent_report_view": clone(parent["report_view"]),
        "accepted_baseline_report_view": clone(parent["accepted"]["report_view"]),
        "failed_parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "accepted_baseline_result_raw_sha256": SOURCE_INTERFACE_INPUT_PARENT_SHA256,
        "false_flag_ledger": clone(parent["false_flag_values"]),
        "retained_maps": clone(parent["retained"]),
        "preserve_failed_math_without_promotion": True,
        "opaque_to_native_interface": True,
    }


def build_input_spec(parent: dict, facts_path: Path, facts: dict, facts_raw: bytes, alias_path: Path, alias: dict, alias_raw: bytes, design_review_path: Path, design_review: dict, design_review_raw: bytes, review_meta: dict, api_meta: dict) -> tuple[bytes, dict]:
    facts_artifact = artifact_entry(facts_path, SOURCE_FACTS_SHA256, embed=True, outside_science_pins=True)
    alias_artifact = artifact_entry(alias_path, PARAMETER_ALIAS_TRACE_SHA256, embed=True, outside_science_pins=True)
    design_artifact = artifact_entry(design_review_path, DESIGN_REVIEW_SHA256, embed=True, outside_science_pins=True)
    api_artifact = artifact_entry(KERNEL_API_DEFAULT, KERNEL_API_SHA256, embed=True, outside_science_pins=True)
    draft_artifact = artifact_entry(KERNEL_DRAFT_DEFAULT, KERNEL_DRAFT_SHA256, embed=True, outside_science_pins=True)
    value = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "source": {
            "alias_trace": clone(alias),
            "alias_trace_raw_sha256": PARAMETER_ALIAS_TRACE_SHA256,
            "source_trace": clone(facts),
            "source_trace_raw_sha256": SOURCE_FACTS_SHA256,
            "parent_data": parent_data_packet(parent),
        },
        "alias_trace": clone(alias),
        "alias_trace_raw_sha256": PARAMETER_ALIAS_TRACE_SHA256,
        "source_trace": clone(facts),
        "source_trace_raw_sha256": SOURCE_FACTS_SHA256,
        "parent_data": parent_data_packet(parent),
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "accepted_parent_report_raw_sha256": ACCEPTED_BASELINE_REPORT_SHA256,
        "accepted_parent_result_raw_sha256": SOURCE_INTERFACE_INPUT_PARENT_SHA256,
        "parent_result": clone(parent["result"]),
        "parent_result_data_only": True,
        "parent_math_failure_preserved": True,
        "parent_flags_not_promoted": True,
        "parent_report_view": clone(parent["report_view"]),
        "parent_false_flags": list(parent["false_flags"]),
        "parent_false_flag_values": clone(parent["false_flag_values"]),
        "preserved_error_cap_metadata": clone(parent["retained"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(parent["retained"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(parent["retained"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(parent["retained"]["retained_cap_provenance"]),
        "sector_source_facts_artifact": facts_artifact,
        "parameter_alias_trace_artifact": alias_artifact,
        "design_review_artifact": design_artifact,
        "saved_parent_review": clone(review_meta),
        "review": clone(review_meta),
        "kernel_api_metadata": clone(api_meta),
        "kernel_api_artifact": api_artifact,
        "external_draft": draft_artifact,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "kernel_contract": clone(KERNEL_CONTRACT),
        "claims": {"native_handles_bound": True, "enclosure_parameters_separated": True, "template_uninstantiated": True, "histories_deferred": True, "native_rebuilt_equality": False, "physical_completion": False, "zero_or_error_cap_reduction": False},
        "history_provider_status": "deferred",
        "runtime_model_metadata": None,
        "math_reexecuted": False,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "metadata_only_assembly": True,
    }
    raw = write_new(INPUT_SPEC, value)
    value["raw_sha256"] = sha_bytes(raw)
    return raw, value


def build_source(parent: dict, input_raw: bytes, facts_path: Path, facts: dict, facts_raw: bytes, alias_path: Path, alias: dict, alias_raw: bytes, design_review_path: Path, design_review: dict, design_review_raw: bytes, review_meta: dict, api_meta: dict, kernel_path: Path, runner_path: Path) -> dict:
    kernel_raw = read_bytes(kernel_path)
    runner_raw = read_bytes(runner_path)
    v1_path_value = facts.get("v1_artifact", {}).get("path")
    require(v1_path_value, "V1 source artifact path missing")
    v1_path = resolve_path(v1_path_value, WORK)
    v1_raw = read_bytes(v1_path, SOURCE_FACTS_V1_SHA256)
    additions = {
        wire_path(INPUT_SPEC): sha_bytes(input_raw),
        wire_path(kernel_path): sha_bytes(kernel_raw),
        wire_path(runner_path): sha_bytes(runner_raw),
    }
    source_union = merge_hashes(parent["source_pins"], additions)
    facts_artifact = artifact_entry(facts_path, SOURCE_FACTS_SHA256, embed=True, outside_science_pins=True)
    alias_artifact = artifact_entry(alias_path, PARAMETER_ALIAS_TRACE_SHA256, embed=True, outside_science_pins=True)
    design_artifact = artifact_entry(design_review_path, DESIGN_REVIEW_SHA256, embed=True, outside_science_pins=True)
    api_artifact = artifact_entry(KERNEL_API_DEFAULT, KERNEL_API_SHA256, embed=True, outside_science_pins=True)
    draft_artifact = artifact_entry(KERNEL_DRAFT_DEFAULT, KERNEL_DRAFT_SHA256, embed=True, outside_science_pins=True)
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_parent_commit": PARENT_PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "accepted_parent_report_raw_sha256": ACCEPTED_BASELINE_REPORT_SHA256,
        "parent_report_view": clone(parent["report_view"]),
        "parent_result": clone(parent["result"]),
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "accepted_parent_result_raw_sha256": SOURCE_INTERFACE_INPUT_PARENT_SHA256,
        "parent_result_data_only": True,
        "parent_math_failure_preserved": True,
        "parent_flags_not_promoted": True,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "parent_checkpoint_hashes": merge_hashes(parent["parent_checkpoints"], parent["current_checkpoints"]),
        "parent_source_pins": clone(parent["source_pins"]),
        "parent_protected_hashes": clone(parent["protected"]),
        "parent_dirty_hashes": clone(parent["dirty"]),
        "parent_failed_hashes": clone(parent["failed"]),
        "parent_false_flags": list(parent["false_flags"]),
        "parent_false_flag_values": clone(parent["false_flag_values"]),
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "preserved_error_cap_metadata": clone(parent["retained"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(parent["retained"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(parent["retained"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(parent["retained"]["retained_cap_provenance"]),
        "sector_source_facts": clone(facts),
        "sector_source_facts_raw_sha256": SOURCE_FACTS_SHA256,
        "sector_source_facts_path": wire_path(facts_path),
        "sector_source_facts_artifact": facts_artifact,
        "sector_source_facts_v1_artifact": artifact_entry(v1_path, SOURCE_FACTS_V1_SHA256, embed=True, outside_science_pins=True),
        "alias_trace": clone(alias),
        "alias_trace_raw_sha256": PARAMETER_ALIAS_TRACE_SHA256,
        "parameter_alias_trace_artifact": alias_artifact,
        "source_trace": clone(facts),
        "source_trace_raw_sha256": SOURCE_FACTS_SHA256,
        "design_review": clone(design_review),
        "design_review_artifact": design_artifact,
        "saved_parent_review": clone(review_meta),
        "review": clone(review_meta),
        "kernel_api_metadata": clone(api_meta),
        "parent_data": parent_data_packet(parent),
        "input_spec_metadata": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw), "bytes": len(input_raw), "raw_readonly": True},
        "input_raw_sha256": sha_bytes(input_raw),
        "kernel_path": wire_path(kernel_path),
        "kernel_sha256": sha_bytes(kernel_raw),
        "runner_path": wire_path(runner_path),
        "runner_sha256": sha_bytes(runner_raw),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "kernel_contract": clone(KERNEL_CONTRACT),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "raw_sha_pins": source_union,
        "outside_science_artifacts": [facts_artifact, alias_artifact, design_artifact, api_artifact, draft_artifact, artifact_entry(v1_path, SOURCE_FACTS_V1_SHA256, embed=False, outside_science_pins=True)],
        "source_packet_defined": True,
        "source_packet_definition_accepted": True,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "native_rebuilt_parameter_equality": False,
        "histories_deferred": True,
        "runtime_model_metadata": None,
        "math_reexecuted": False,
        "old_math_reexecuted": False,
        "no_runner_execution": True,
    }


def load_source() -> tuple[bytes, dict]:
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "Source packet identity changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    require(source.get("parent_report_raw_sha256") == PARENT_REPORT_SHA256 and source.get("parent_result_raw_sha256") == PARENT_RESULT_SHA256, "Source failed-parent pins changed")
    require(source.get("accepted_parent_report_raw_sha256") == ACCEPTED_BASELINE_REPORT_SHA256 and source.get("accepted_parent_result_raw_sha256") == SOURCE_INTERFACE_INPUT_PARENT_SHA256, "Source accepted-baseline pins changed")
    require(source.get("sector_source_facts_raw_sha256") == SOURCE_FACTS_SHA256 and source.get("alias_trace_raw_sha256") == PARAMETER_ALIAS_TRACE_SHA256, "Source trace pins changed")
    require(source.get("kernel_contract") == KERNEL_CONTRACT, "Kernel contract changed after source freeze")
    require(source.get("metadata_only") is True and source.get("source_only") is True and source.get("physical_layer_completed") is False, "Source layer marker changed")
    flags = source.get("parent_false_flags")
    values = source.get("parent_false_flag_values")
    require(isinstance(flags, list) and len(flags) == PARENT_FALSE_FLAG_COUNT and len(set(flags)) == PARENT_FALSE_FLAG_COUNT, "Source false-flag ledger changed")
    require(isinstance(values, dict) and values == {key: False for key in flags}, "Source false values changed")
    for name in RETAINED_MAP_NAMES:
        require(isinstance(source.get(name), dict), f"Source retained map missing: {name}")
    facts_artifact = source.get("sector_source_facts_artifact")
    require(isinstance(facts_artifact, dict) and facts_artifact.get("path"), "Nested source facts artifact missing")
    facts_path = resolve_path(facts_artifact["path"], WORK)
    facts_raw, facts = read_json(facts_path, SOURCE_FACTS_SHA256)
    require(source.get("sector_source_facts") == facts and source.get("source_trace") == facts, "Frozen source facts changed")
    alias_artifact = source.get("parameter_alias_trace_artifact")
    require(isinstance(alias_artifact, dict) and alias_artifact.get("path"), "Nested alias artifact missing")
    alias_path = resolve_path(alias_artifact["path"], WORK)
    alias_raw, alias = read_json(alias_path, PARAMETER_ALIAS_TRACE_SHA256)
    require(source.get("alias_trace") == alias, "Frozen alias trace changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name, WORK), digest)
    return raw, source


def append_diagnostic(event: dict) -> bytes:
    require(isinstance(event, dict), "Diagnostic event must be an object")
    raw = encode(event)
    DIAGNOSTIC_CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    with DIAGNOSTIC_CHECKPOINT.open("ab") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


CALLBACK_NAMES = {
    "source_validation": "source_validation",
    "native_handle_admission": "native_handle_admission",
    "source_reference_rejection": "source_reference_rejection",
    "template_signature": "template_signature",
}


def diagnostic_summary() -> dict:
    lines = DIAGNOSTIC_CHECKPOINT.read_bytes().splitlines()
    events = []
    for line in lines:
        require(line, "Empty diagnostic line")
        event = json.loads(line.decode("utf-8"))
        require(isinstance(event, dict), "Diagnostic event is not an object")
        events.append(event)
    require(len(events) >= 2, "Diagnostic manager prefix is incomplete")
    require(events[0].get("event") == "manager_entry_started", "Diagnostic prefix changed")
    require(events[1].get("event") in ("manager_kernel_imported", "manager_kernel_or_derive_failed"), "Diagnostic import prefix changed")
    callbacks = events[2:]
    counts = {}
    typed = {name: 0 for name in CALLBACK_EVENT_COUNTS}
    bare = True
    for event in callbacks:
        if "manager_event" in event:
            bare = False
        name = event.get("event_kind", event.get("event"))
        counts[name] = counts.get(name, 0) + 1
        category = CALLBACK_NAMES.get(name)
        if category is not None:
            typed[category] += 1
    return {"event_count": len(events), "manager_prefix_event_count": 2, "callback_event_count": len(callbacks), "callback_event_counts": counts, "typed_callback_event_counts": typed, "callback_events_are_bare": bare, "manager_kernel_imported": events[1].get("event") == "manager_kernel_imported", "callback_event_raw_rows": [clone(event) for event in callbacks]}


def import_kernel(path: Path):
    spec = importlib.util.spec_from_file_location("_native_W_parameter_interface_kernel", path)
    require(spec is not None and spec.loader is not None, "Kernel module spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "Kernel derive(source, diagnostic=None) is missing")
    return module


def failure_result(error: BaseException, derive_entry_calls: int, derive_calls: int, kernel_import_calls: int) -> dict:
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": True,
        "source_only": True,
        "status": "failed",
        "actual_counts": {name: 0 for name in COUNT_KEYS},
        "new_flags": {name: False for name in FLAG_KEYS},
        "payload": {"status": "manager_failure_prefix", "source_validation_events": [], "native_source_handle_records": [], "rebuilt_original_records": [], "template_signature_records": [], "native_handle_admission_events": [], "source_reference_rejection_events": []},
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_import_calls == 1,
        "error_type": type(error).__name__,
        "error_message": str(error),
        "runtime_model_metadata": None,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "cap_error_reduction_claimed": False,
        "error_cap_reduced": False,
        "old_error_terms_removed": False,
        "old_errors_removed": False,
        "failure": {"kind": "manager_import_or_derive_failure", "error_type": type(error).__name__, "error_message": str(error)},
    }


def validate_result(result: dict, source: dict, derive_entry_calls: int | None = None, derive_calls: int | None = None, diagnostic: dict | None = None) -> dict:
    violations: list[str] = []
    if not isinstance(result, dict):
        return {"passed": False, "protocol_violations": ["result is not an object"], "actual_counts": None, "new_flags": {}, "metadata_only": True, "source_only": True}
    passed = result.get("passed")
    if not isinstance(passed, bool):
        violations.append("passed is not an explicit Boolean")
    if result.get("schema") not in RESULT_SCHEMA_ALIASES:
        violations.append("result schema changed")
    if result.get("component") != COMPONENT:
        violations.append("result component changed")
    if result.get("metadata_only") is not True:
        violations.append("metadata_only marker is not true")
    if result.get("source_only") is not True:
        violations.append("source_only marker is not true")
    if result.get("physical_layer_completed") is True or result.get("full_goal_complete") is True:
        violations.append("physical or full-goal completion was claimed")
    if passed is True and result.get("failure") is not None:
        violations.append("successful result carries a failure payload")
    counts = result.get("actual_counts")
    if passed is True:
        if counts != EXPECTED_COUNTS:
            violations.append("actual_counts do not exactly match the nine-key contract")
    elif not isinstance(counts, dict):
        violations.append("failed result lacks partial actual_counts")
    else:
        if set(counts).difference(COUNT_KEYS):
            violations.append("failed result has unknown actual_counts keys")
        for key, value in counts.items():
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                violations.append(f"failed result count is not nonnegative: {key}")
    flags = result.get("new_flags")
    if passed is True:
        if flags != EXPECTED_FLAGS:
            violations.append("new_flags do not exactly match the three-key contract")
    elif not isinstance(flags, dict) or any(not isinstance(value, bool) for value in flags.values()):
        violations.append("failed result flags are not Boolean metadata")
    payload = result.get("payload")
    if passed is True:
        if not isinstance(payload, dict):
            violations.append("successful result payload is missing")
        else:
            expected_lengths = {"native_source_handle_records": 2, "rebuilt_original_records": 3, "template_signature_records": 1, "native_handle_admission_events": 2, "source_reference_rejection_events": 2}
            for key, length in expected_lengths.items():
                rows = payload.get(key)
                if not isinstance(rows, list) or len(rows) != length or not all(isinstance(row, dict) for row in rows):
                    violations.append(f"payload row namespace changed: {key}")
    elif payload is not None and not isinstance(payload, dict):
        violations.append("failed payload is not an object")
    for key in ("parent_false_flags", "inherited_false_flags"):
        if key in result and result[key] != source.get("parent_false_flags"):
            violations.append(f"inherited false-flag list changed: {key}")
    for key in ("parent_false_flag_values", "inherited_false_flag_values"):
        if key in result and result[key] != source.get("parent_false_flag_values"):
            violations.append(f"inherited false-flag values changed: {key}")
    for key in RETAINED_MAP_NAMES:
        if key in result and result[key] != source.get(key):
            violations.append(f"retained parent map changed: {key}")
    for key in NONCLAIM_FLAGS:
        if key in result and result[key] is not False:
            violations.append(f"forbidden claim is true: {key}")
    claims = result.get("claims")
    if isinstance(claims, dict):
        for key in ("native_rebuilt_equality", "physical_completion", "zero_or_error_cap_reduction", "history_instantiated"):
            if key in claims and claims[key] is not False:
                violations.append(f"forbidden claim is true: claims.{key}")
        if claims.get("histories_deferred") is not True and "histories_deferred" in claims:
            violations.append("history deferral marker changed")
    if derive_entry_calls is not None and derive_entry_calls not in (0, 1):
        violations.append("derive_entry_calls is not a measured zero/one")
    if derive_calls is not None and derive_calls not in (0, 1):
        violations.append("derive_calls is not a measured zero/one")
    if passed is True and (derive_entry_calls != 1 or derive_calls != 1):
        violations.append("successful derive lifecycle did not run exactly once")
    if diagnostic is not None:
        callback_count = diagnostic.get("callback_event_count")
        typed = diagnostic.get("typed_callback_event_counts", {})
        if passed is True:
            if callback_count != CALLBACK_EVENT_COUNT:
                violations.append("successful diagnostic callback count is not six")
            if typed != CALLBACK_EVENT_COUNTS:
                violations.append("successful diagnostic event categories changed")
            if diagnostic.get("callback_events_are_bare") is not True:
                violations.append("kernel callback rows were wrapped by manager")
        elif callback_count is None or not isinstance(callback_count, int) or callback_count < 0 or callback_count > CALLBACK_EVENT_COUNT:
            violations.append("failed diagnostic callback count is not an honest partial count")
    return {"passed": not violations and passed is True, "protocol_violations": violations, "actual_counts": clone(counts) if isinstance(counts, dict) else None, "new_flags": clone(flags) if isinstance(flags, dict) else {}, "derive_entry_calls": derive_entry_calls, "derive_calls": derive_calls, "diagnostic": clone(diagnostic) if isinstance(diagnostic, dict) else None, "metadata_only": True, "source_only": True, "physical_layer_completed": False, "expected_counts": clone(EXPECTED_COUNTS), "expected_new_flags": clone(EXPECTED_FLAGS)}


def update_main_brief(patch: dict) -> tuple[bytes, dict]:
    require(isinstance(patch, dict), "Main BRIEF patch must be an object")
    raw = read_bytes(MAIN_BRIEF)
    current = json.loads(raw.decode("utf-8"))
    require(isinstance(current, dict), "Main BRIEF must be an object")
    safe_patch = brief_safe(portable(clone(patch)))
    if "latest_authoritative_state" in patch:
        authority = brief_safe(clone(patch["latest_authoritative_state"]))
        require(isinstance(authority, dict), "Latest authoritative state must be an object")
        current["latest_authoritative_state"] = authority
        current["authority"] = authority
        if authority.get("verified_commit"):
            current["verified_commit"] = authority["verified_commit"]
    for key, value in safe_patch.items():
        if key != "latest_authoritative_state":
            current[key] = value
    if patch.get("current_component") or patch.get("component"):
        current["latest_component"] = patch.get("current_component") or patch.get("component")
        current["latest_component_pointer"] = safe_patch
    current["current_exact_source_result"] = safe_patch
    milestone_path = REPO / MILESTONE_RELATIVE
    if milestone_path.exists() and patch.get("phase", "").startswith(("published", "resumed", "sealed")):
        current["latest_milestone"] = {
            "path": wire_path(milestone_path),
            "report": wire_path(REPO / REPORT_RELATIVE),
            "phase": patch.get("phase"),
            "report_sha256": patch.get("report_sha256"),
            "milestone_sha256": sha_path(milestone_path),
        }
    compact = brief_safe(current)
    compact_raw = encode(compact)
    if len(compact_raw) >= 100_000:
        archive = WORK / f"RESUME_BRIEF_ARCHIVE_BEFORE_{COMPONENT}_2026_10_27.json"
        if not archive.exists():
            write_new(archive, raw)
        pointer = {"path": wire_path(archive), "sha256": sha_bytes(raw), "bytes": len(raw), "preimage_exact": True}
        for key in ("previous_in_progress", "in_progress"):
            if key in compact:
                compact[key] = {"preimage_pointer": pointer, "selectors_compacted": True}
        compact = brief_safe(compact)
        compact_raw = encode(compact)
    require(len(compact_raw) < 100_000, "Main BRIEF exceeds 100 KB after compaction")
    write_replace(MAIN_BRIEF, compact)
    return compact_raw, compact


def checkpoint_artifacts() -> list[dict]:
    return [artifact_entry(path, embed=True) for path in CHECKPOINTS if path.exists()]


def code_artifacts(kernel_path: Path, runner_path: Path) -> list[dict]:
    return [artifact_entry(kernel_path, embed=True), artifact_entry(runner_path, embed=True)]


def failed_attempt_artifacts(code2: list[dict], input_artifact: dict, current_checkpoints: list[dict]) -> dict:
    additions = {entry["path"]: entry["sha256"] for entry in code2}
    additions[input_artifact["path"]] = input_artifact["sha256"]
    additions.update({entry["path"]: entry["sha256"] for entry in current_checkpoints})
    require(len(additions) == 7, "A failed interface attempt must add exactly seven artifacts")
    return additions


def alias_trace_pointer(alias_artifact: dict) -> dict:
    return {"path": alias_artifact["path"], "sha256": alias_artifact["sha256"], "bytes": alias_artifact["bytes"], "outside_science_pins": True, "static_only": True, "no_cross_side_parameter_identity_claim": True, "no_zero_or_error_cap_claim": True}


def next_task_metadata(passed: bool, alias_artifact: dict | None = None) -> dict:
    if not passed:
        return {"phase": "cached_native_W_parameter_interface_prefix_inspection", "action": "Inspect only the saved native parameter interface failure prefix and immutable result; do not reenter or replay the failed W-sector family.", "source": "saved_native_W_parameter_interface_result_and_diagnostic_prefix", "reentry": "forbidden", "physical_layer_completed": False, "claims": "no native/rebuilt equality, zero, error-cap reduction, closure, temporal recursion, or physical completion"}
    require(isinstance(alias_artifact, dict), "Passed next task requires alias metadata")
    return {"phase": "source_native_history_provider_for_rebuilt_W", "action": "Define accepted source-native history and radius interfaces for the uninstantiated rebuilt W template; keep native Pstar/delta and enclosure-origin records distinct.", "source": "saved_native_W_parameter_interface_and_static_alias_trace", "parameter_alias_trace_pointer": alias_trace_pointer(alias_artifact), "retain_unresolved": ["source-native history providers", "rebuilt own radius", "native/rebuilt Pstar relation", "native/rebuilt delta relation"], "reentry": "forbidden", "physical_layer_completed": False, "claims": "no source value computation, native/rebuilt equality, zero, error-cap deletion, closure, temporal recursion, or physical completion"}


def report_for_publish(source: dict, source_raw: bytes, result: dict, result_raw: bytes, validation: dict, code2: list[dict], input_artifact: dict, review_artifact: dict) -> dict:
    current = checkpoint_artifacts()
    current_hashes = {entry["path"]: entry["sha256"] for entry in current}
    if validation["passed"]:
        failed_hashes = clone(source["parent_failed_hashes"])
    else:
        failed_hashes = merge_hashes(source["parent_failed_hashes"], failed_attempt_artifacts(code2, input_artifact, current))
    source_union = merge_hashes(source["source_union_hashes"], current_hashes)
    helper = artifact_entry(Path(__file__).resolve(), embed=True, outside_science_pins=True)
    alias_artifact = source["parameter_alias_trace_artifact"]
    facts_artifact = source["sector_source_facts_artifact"]
    facts_v1_artifact = source["sector_source_facts_v1_artifact"]
    design_artifact = source["design_review_artifact"]
    api_artifact = artifact_entry(KERNEL_API_DEFAULT, KERNEL_API_SHA256, embed=True, outside_science_pins=True)
    draft_artifact = artifact_entry(KERNEL_DRAFT_DEFAULT, KERNEL_DRAFT_SHA256, embed=True, outside_science_pins=True)
    next_task = next_task_metadata(validation["passed"], alias_artifact if validation["passed"] else None)
    component_result = result.get("component_result", result)
    lossless = [*code2, input_artifact, *current, helper, review_artifact]
    require(len(lossless) == 9, "Lossless publication set must be code2 plus input, four checkpoints, helper, and review")
    return {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "passed": validation["passed"],
        "explicit_result_passed": result.get("passed"),
        "result_status": result.get("status"),
        "metadata_only": True,
        "publication_metadata_only": True,
        "source_only": True,
        "result_metadata_only": result.get("metadata_only"),
        "result_source_only": result.get("source_only"),
        "binding_only": True,
        "source_binding_only": True,
        "component_result_metadata_only": component_result.get("metadata_only") if isinstance(component_result, dict) else None,
        "component_result_source_only": component_result.get("source_only") if isinstance(component_result, dict) else None,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "cap_error_reduction_claimed": False,
        "error_cap_reduced": False,
        "old_error_terms_removed": False,
        "old_errors_removed": False,
        "parent_math_failure_preserved": True,
        "parent_flags_not_promoted": True,
        "accepted_progress_physical_baseline_preserved": True,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_parent_commit": PARENT_PARENT_COMMIT,
        "parent_report": source["parent_report"],
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "accepted_parent_report_raw_sha256": ACCEPTED_BASELINE_REPORT_SHA256,
        "accepted_parent_result_raw_sha256": SOURCE_INTERFACE_INPUT_PARENT_SHA256,
        "parent_result_data_only": True,
        "parent_report_view": source["parent_report_view"],
        "parent_false_flags": source["parent_false_flags"],
        "parent_false_flag_values": source["parent_false_flag_values"],
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "preserved_error_cap_metadata": source["preserved_error_cap_metadata"],
        "retained_construction_data": source["retained_construction_data"],
        "retained_integral_term_inventory": source["retained_integral_term_inventory"],
        "retained_cap_provenance": source["retained_cap_provenance"],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "checkpoint_count": PARENT_CHECKPOINT_COUNT + 4,
        "checkpoint_paths": list(current_hashes),
        "checkpoint_hashes": current_hashes,
        "source_pin_count": len(source_union),
        "source_union_count": len(source_union),
        "source_pins": source_union,
        "source_union_hashes": source_union,
        "raw_sha_pins": source_union,
        "preserved_failed_artifact_hashes_union": failed_hashes,
        "preserved_failed_artifact_count": len(failed_hashes),
        "new_failed_artifact_count": len(failed_hashes) - PARENT_FAILED_COUNT,
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "protected_file_hashes": source["parent_protected_hashes"],
        "legacy_dirty_path_hashes": source["parent_dirty_hashes"],
        "sector_source_facts": source["sector_source_facts"],
        "sector_source_facts_raw_sha256": SOURCE_FACTS_SHA256,
        "sector_source_facts_artifact": facts_artifact,
        "sector_source_facts_v1_artifact": facts_v1_artifact,
        "source_trace": source["source_trace"],
        "source_trace_raw_sha256": SOURCE_FACTS_SHA256,
        "parameter_alias_trace_artifact": alias_artifact,
        "parameter_alias_trace_pointer": alias_trace_pointer(alias_artifact),
        "design_review_artifact": design_artifact,
        "kernel_api_artifact": api_artifact,
        "external_draft": draft_artifact,
        "input_artifact": input_artifact,
        "input_raw_sha256": input_artifact["sha256"],
        "code2": code2,
        "current_code_hashes": {entry["path"]: entry["sha256"] for entry in code2},
        "metadata_helper_artifact": helper,
        "full_result": clone(result),
        "component_result": clone(component_result),
        "component_result_exact": True,
        "full_result_raw_bytes": len(result_raw),
        "full_result_raw_sha256": sha_bytes(result_raw),
        "raw_bytes_lossless": True,
        "raw_result_preserved": True,
        "result_checkpoint": {"path": wire_path(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw)},
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_path(ENTERED_CHECKPOINT),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "actual_counts": validation["actual_counts"],
        "new_flags": validation["new_flags"],
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "protocol_violations": validation["protocol_violations"],
        "durable_result_validation": validation,
        "execution_metadata": clone(result.get("execution_metadata", {})),
        "kernel_import_calls": result.get("kernel_import_calls"),
        "kernel_imported": result.get("kernel_imported"),
        "derive_entry_calls": result.get("derive_entry_calls"),
        "derive_calls": result.get("derive_calls"),
        "runner_executed": False,
        "no_runner_execution": True,
        "lossless_artifacts": lossless,
        "lossless_artifact_count": len(lossless),
        "independent_review": review_artifact,
        "next_executable_task": next_task,
        "runtime_model_metadata": None,
        "exact_nine_scope": True,
        "scoped_paths": [path.as_posix() for path in SCOPED_RELATIVES],
        "input_spec_schema": INPUT_SCHEMA,
    }


def git_output(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def validate_scope_attributes() -> None:
    attrs_path = REPO / ".gitattributes"
    text = attrs_path.read_text(encoding="utf-8") if attrs_path.exists() else ""
    lines = set(text.splitlines())
    for relative in (KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, MILESTONE_RELATIVE):
        require(f"/{relative.as_posix()} -text" in lines, f"Missing -text scope attribute: {relative}")


def verify_published_commit(report: dict) -> str:
    head = git_output("rev-parse", "HEAD")
    require(valid_revision(head), "Local HEAD is not a 40-hex revision")
    require(git_output("rev-parse", "HEAD^").lower() == PARENT_COMMIT, "Published parent changed")
    branch = git_output("branch", "--show-current")
    require(branch, "Published branch is unnamed")
    remote = git_output("ls-remote", "origin", "refs/heads/" + branch).split()
    require(remote and valid_revision(remote[0]) and remote[0].lower() == head.lower(), "Remote branch is not synchronized")
    changed = set(git_output("diff", "--name-only", PARENT_COMMIT + "..HEAD").splitlines())
    require(changed == {path.as_posix() for path in SCOPED_RELATIVES}, "Committed scope is not exactly nine paths")
    validate_scope_attributes()
    for relative in SCOPED_RELATIVES:
        local = REPO / relative
        require(local.exists(), "Published scoped path is missing: " + relative.as_posix())
        require(int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix())) == local.stat().st_size, "Published blob size changed")
        require(git_output("rev-parse", "HEAD:" + relative.as_posix()).lower() == git_output("hash-object", wire_path(local)).lower(), "Published blob bytes changed")
    require(report.get("parent_commit") == PARENT_COMMIT, "Report parent commit changed")
    return head


def publish(args: argparse.Namespace) -> dict:
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT)
    _, run_receipt = read_json(RUN_RECEIPT)
    validation = validate_result(result, source, run_receipt.get("derive_entry_calls"), run_receipt.get("derive_calls"), run_receipt.get("diagnostic_summary"))
    kernel_path = resolve_path(source["kernel_path"], REPO)
    runner_path = resolve_path(source["runner_path"], REPO)
    code2 = code_artifacts(kernel_path, runner_path)
    input_artifact = artifact_entry(resolve_path(source["input_spec_metadata"]["path"], WORK), source["input_spec_metadata"]["sha256"], embed=True)
    review_artifact = artifact_entry(resolve_path(source["saved_parent_review"]["path"], WORK), source["saved_parent_review"]["sha256"], embed=True, outside_science_pins=True)
    report = report_for_publish(source, source_raw, result, result_raw, validation, code2, input_artifact, review_artifact)
    report_path = REPO / REPORT_RELATIVE
    report_raw = write_replace(report_path, report)
    milestone_text = (
        "# Current native W parameter interface (2026-10-27)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Status: `{report['attempt_status']}`; actual counts: `{report['actual_counts']}`.\n\n"
        "This publication is a source-only metadata binding. It records exact native Pstar and delta handles, keeps rebuilt enclosure-origin records distinct, and defines an uninstantiated template signature with history providers deferred. It does not claim native/rebuilt equality, W evaluation, physical completion, numeric values, zero or error-cap reduction, closure, or temporal recursion.\n\n"
        "The saved result is immutable and must not be reentered.\n"
    ).encode("utf-8")
    milestone_path = REPO / MILESTONE_RELATIVE
    milestone_raw = write_replace(milestone_path, milestone_text)
    prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; do not reenter.\n"
        f"Next phase: `{report['next_executable_task']['phase']}`.\n"
        f"Next action: {report['next_executable_task']['action']}\n"
        f"Static parameter-alias trace: `{report['parameter_alias_trace_pointer']['path']}` SHA-256 `{report['parameter_alias_trace_pointer']['sha256']}`; no cross-side identity or zero/error-cap claim.\n"
    ).encode("utf-8")
    preimages = []
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(old)[:12]}.md"
        if not preimage.exists():
            write_new(preimage, old)
        preimages.append({"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(old)})
        write_replace(path, prefix + old)
    attrs_path = REPO / ".gitattributes"
    attrs_text = attrs_path.read_text(encoding="utf-8") if attrs_path.exists() else ""
    for relative in (KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, MILESTONE_RELATIVE):
        line = f"/{relative.as_posix()} -text"
        if line not in attrs_text.splitlines():
            attrs_text += ("\n" if attrs_text and not attrs_text.endswith("\n") else "") + line + "\n"
    write_replace(attrs_path, attrs_text.encode("utf-8"))
    brief_raw, brief = update_main_brief({
        "phase": "published_passed" if report["passed"] else "published_failed",
        "current_component": COMPONENT,
        "report": wire_path(report_path),
        "report_sha256": sha_bytes(report_raw),
        "milestone": wire_path(milestone_path),
        "milestone_sha256": sha_bytes(milestone_raw),
        "attempt_status": report["attempt_status"],
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "latest_review": review_artifact,
        "latest_postrun_saved_review": review_artifact,
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
        "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
        "current_document_sha256": sha_bytes(milestone_raw),
        "current_doc_sha256": sha_bytes(milestone_raw),
        "next_task": report.get("next_executable_task"),
        "next_executable_task": report.get("next_executable_task"),
        "read_first": [wire_path(REPO / "AGENTS.md"), wire_path(REPO / MILESTONE_RELATIVE), wire_path(REPO / REPORT_RELATIVE)],
    })
    receipt = {"schema": COMPONENT + "_publication_receipt_v1", "component": COMPONENT, "report": wire_path(report_path), "report_sha256": sha_bytes(report_raw), "attempt_status": report["attempt_status"], "source_checkpoint_sha256": sha_bytes(source_raw), "result_checkpoint_sha256": sha_bytes(result_raw), "handoff_preimages": preimages, "exact_nine_scope": True, "runtime_model_metadata": None, "main_resume_brief": {"sha256": sha_bytes(brief_raw), "value": brief}}
    receipt_raw = write_replace(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}


def resume(args: argparse.Namespace) -> dict:
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    require(report.get("full_result") == result, "Published result changed")
    stored = report.get("durable_result_validation")
    require(isinstance(stored, dict), "Durable validation is missing")
    validation = validate_result(result, source, stored.get("derive_entry_calls"), stored.get("derive_calls"), stored.get("diagnostic"))
    require(validation == stored, "Durable validation changed")
    head = verify_published_commit(report)
    milestone_sha256 = sha_path(REPO / MILESTONE_RELATIVE)
    brief_raw, brief = update_main_brief({
        "phase": "resumed_published",
        "current_component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "head": head,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "read_first": [wire_path(REPO / "AGENTS.md"), wire_path(REPO / MILESTONE_RELATIVE), wire_path(REPO / REPORT_RELATIVE)],
        "latest_review": report.get("independent_review"),
        "latest_postrun_saved_review": report.get("independent_review"),
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
        "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
        "current_document_sha256": milestone_sha256,
        "current_doc_sha256": milestone_sha256,
        "next_task": report.get("next_executable_task"),
        "next_executable_task": report.get("next_executable_task"),
    })
    receipt = {"schema": COMPONENT + "_resume_receipt_v1", "component": COMPONENT, "report": wire_path(REPO / REPORT_RELATIVE), "report_sha256": sha_bytes(report_raw), "source_checkpoint_sha256": sha_bytes(source_raw), "result_checkpoint_sha256": sha_bytes(result_raw), "head": head, "main_resume_brief_sha256": sha_bytes(brief_raw), "main_resume_brief": brief, "do_not_reenter": True}
    receipt_raw = write_replace(RESUME_RECEIPT, receipt)
    return {"resume": receipt, "resume_sha256": sha_bytes(receipt_raw)}


def seal(args: argparse.Namespace) -> dict:
    review_path = resolve_path(getattr(args, "review", None) or PARENT_REVIEW_DEFAULT, WORK)
    review_raw, review_meta = load_review(review_path)
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    require(report.get("full_result") == result, "Published result changed before seal")
    stored = report.get("durable_result_validation")
    require(isinstance(stored, dict), "Durable validation is missing before seal")
    validation = validate_result(result, source, stored.get("derive_entry_calls"), stored.get("derive_calls"), stored.get("diagnostic"))
    require(validation == stored, "Durable validation changed before seal")
    head = verify_published_commit(report)
    milestone_sha256 = sha_path(REPO / MILESTONE_RELATIVE)
    authority = {
        "schema": COMPONENT + "_authoritative_state_v1",
        "component": COMPONENT,
        "commit": head,
        "head": head,
        "verified_commit": head,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "passed": report.get("passed") is True,
        "attempt_status": report.get("attempt_status"),
        "actual_counts": clone(report.get("actual_counts")),
        "new_flags": clone(report.get("new_flags")),
        "full_goal_complete": False,
        "physical_layer_completed": False,
        "metadata_only": report.get("metadata_only") is True,
        "source_only": report.get("source_only") is True,
        "accepted_interface_only": report.get("passed") is True,
        "old_failure_math_accepted": False,
        "parent_math_failure_preserved": True,
        "do_not_reenter": True,
        "review": review_meta,
        "next_remaining": report.get("next_executable_task"),
        "read_first": [wire_path(REPO / "AGENTS.md"), wire_path(REPO / MILESTONE_RELATIVE), wire_path(REPO / REPORT_RELATIVE)],
    }
    brief_raw, brief = update_main_brief({
        "phase": "sealed",
        "current_component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "seal": wire_path(SEAL_RECEIPT),
        "head": head,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "latest_authoritative_state": authority,
        "latest_review": review_meta,
        "latest_postrun_saved_review": review_meta,
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "current_document_sha256": milestone_sha256,
        "current_doc_sha256": milestone_sha256,
        "next_task": report.get("next_executable_task"),
        "next_executable_task": report.get("next_executable_task"),
    })
    receipt = {"schema": SEAL_SCHEMA, "component": COMPONENT, "report": wire_path(REPO / REPORT_RELATIVE), "report_sha256": sha_bytes(report_raw), "source_checkpoint_sha256": sha_bytes(source_raw), "result_checkpoint_sha256": sha_bytes(result_raw), "review": review_meta, "head": head, "authority": authority, "main_resume_brief_sha256": sha_bytes(brief_raw), "main_resume_brief": brief, "do_not_reenter": True}
    receipt_raw = write_replace(SEAL_RECEIPT, receipt)
    return {"seal": receipt, "seal_sha256": sha_bytes(receipt_raw)}


def prepare(args: argparse.Namespace) -> dict:
    require(not SOURCE_CHECKPOINT.exists(), "Source packet already exists; do not reenter")
    require(not INPUT_SPEC.exists(), "Input specification already exists; do not rewrite")
    require(not BRIEF.exists(), "Component BRIEF already exists; do not rewrite")
    require(all(not path.exists() for path in CHECKPOINTS), "Interface checkpoints already exist; do not reenter")
    require(sha_path(MAIN_BRIEF) == CURRENT_MAIN_BRIEF_SHA256, "Main BRIEF preimage is not the authenticated current document")
    parent = load_parent(args)
    facts_path = resolve_path(getattr(args, "source_facts", None) or SOURCE_FACTS_DEFAULT, WORK)
    facts_sha = (getattr(args, "source_facts_sha256", None) or SOURCE_FACTS_SHA256).lower()
    require(facts_sha == SOURCE_FACTS_SHA256, "Source trace SHA is not authenticated")
    facts_raw, facts = load_source_facts(facts_path, facts_sha)
    alias_path = resolve_path(getattr(args, "alias_trace", None) or PARAMETER_ALIAS_TRACE_DEFAULT, WORK)
    alias_raw, alias, alias_pointer = load_alias_trace(alias_path)
    design_path = resolve_path(getattr(args, "design_review", None) or DESIGN_REVIEW_DEFAULT, WORK)
    design_raw, design_review = load_design_review(design_path)
    review_path = resolve_path(getattr(args, "review", None) or PARENT_REVIEW_DEFAULT, WORK)
    review_raw, review_meta = load_review(review_path)
    api_path = resolve_path(getattr(args, "kernel_api", None) or KERNEL_API_DEFAULT, WORK)
    api_raw, api, api_meta = load_kernel_api(api_path)
    kernel_path = resolve_path(getattr(args, "kernel_path", None) or REPO / KERNEL_RELATIVE, REPO)
    runner_path = resolve_path(getattr(args, "runner_path", None) or REPO / RUNNER_RELATIVE, REPO)
    require(kernel_path.resolve() == (REPO / KERNEL_RELATIVE).resolve(), "Kernel path escaped scope")
    require(runner_path.resolve() == (REPO / RUNNER_RELATIVE).resolve(), "Runner path escaped scope")
    require(kernel_path.exists() and runner_path.exists(), "Kernel and runner must exist before prepare")
    require(kernel_path.resolve() != Path(__file__).resolve() and runner_path.resolve() != Path(__file__).resolve(), "Metadata helper cannot be code2")
    if getattr(args, "kernel_sha256", None):
        require(sha_path(kernel_path) == args.kernel_sha256.lower(), "Kernel SHA does not match")
    if getattr(args, "runner_sha256", None):
        require(sha_path(runner_path) == args.runner_sha256.lower(), "Runner SHA does not match")
    input_raw, input_value = build_input_spec(parent, facts_path, facts, facts_raw, alias_path, alias, alias_raw, design_path, design_review, design_raw, review_meta, api_meta)
    source = build_source(parent, input_raw, facts_path, facts, facts_raw, alias_path, alias, alias_raw, design_path, design_review, design_raw, review_meta, api_meta, kernel_path, runner_path)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    diagnostic_raw = write_new(DIAGNOSTIC_CHECKPOINT, b"")
    brief = {
        "schema": COMPONENT + "_brief_v1",
        "component": COMPONENT,
        "phase": "prepared",
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "source_facts": source["sector_source_facts_artifact"],
        "alias_trace": alias_pointer,
        "kernel_api": api_meta,
        "kernel": {"path": wire_path(kernel_path), "sha256": source["kernel_sha256"]},
        "runner": {"path": wire_path(runner_path), "sha256": source["runner_sha256"]},
        "metadata_only": True,
        "source_only": True,
        "no_runner_execution": True,
        "do_not_reenter": True,
    }
    brief_raw = write_new(BRIEF, brief)
    pending = {
        "schema": COMPONENT + "_pending_authority_v1",
        "component": COMPONENT,
        "current_component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "phase": "prepared_unentered",
        "passed": None,
        "actual_counts": None,
        "new_flags": None,
        "full_goal_complete": False,
        "metadata_only": True,
        "source_only": True,
        "do_not_reenter": True,
    }
    main_raw, main_value = update_main_brief({
        "phase": "prepared_unentered",
        "current_component": COMPONENT,
        "component": COMPONENT,
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "source_union_count": source["source_union_count"],
        "pending_authority": pending,
        "latest_authoritative_state": {"commit": PARENT_COMMIT, "head": PARENT_COMMIT, "verified_commit": PARENT_COMMIT, "passed": None, "actual_counts": None, "new_flags": None, "full_goal_complete": False},
        "execution_rules": {"metadata_only": True, "source_only": True, "no_math_reentry": True, "no_runner_execution": True},
        "dirty_paths": parent["dirty"],
        "protected_paths": parent["protected"],
        "next_task": "native_W_parameter_interface_run_once_after_entry",
        "read_first": [wire_path(SOURCE_CHECKPOINT), wire_path(INPUT_SPEC), wire_path(ENTERED_CHECKPOINT)],
        "metadata_only": True,
        "source_only": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
    })
    receipt = {
        "schema": COMPONENT + "_prepare_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_raw),
        "diagnostic_checkpoint_sha256": sha_bytes(diagnostic_raw),
        "brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief_sha256": sha_bytes(main_raw),
        "main_resume_brief": main_value,
        "pending_authority": pending,
        "source_union_count": source["source_union_count"],
        "parent_commit": PARENT_COMMIT,
        "runtime_model_metadata": None,
        "metadata_only": True,
        "source_only": True,
        "no_runner_execution": True,
    }
    receipt_raw = write_new(PREPARE_RECEIPT, receipt)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw)}


def run(args: argparse.Namespace) -> dict:
    source_raw, source = load_source()
    require(not ENTERED_CHECKPOINT.exists(), "Entered marker already exists; do not reenter")
    require(not RESULT_CHECKPOINT.exists(), "Result already exists; do not rerun")
    require(DIAGNOSTIC_CHECKPOINT.exists(), "Diagnostic checkpoint is missing")
    kernel_path = resolve_path(source["kernel_path"], REPO)
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "entered_at": datetime.now(timezone.utc).isoformat(),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "kernel_path": source["kernel_path"],
        "kernel_sha256": source["kernel_sha256"],
        "derive_contract": KERNEL_CONTRACT["derive"],
        "derive_entry_calls": 0,
        "derive_calls": 0,
        "kernel_import_calls": 0,
        "metadata_only": True,
        "source_only": True,
        "runner_executed": False,
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    derive_entry_calls = 0
    derive_calls = 0
    kernel_import_calls = 0
    kernel_imported = False
    append_diagnostic({"event": "manager_entry_started", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw)})
    try:
        module = import_kernel(kernel_path)
        kernel_import_calls = 1
        kernel_imported = True
        append_diagnostic({"event": "manager_kernel_imported", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "kernel_path": source["kernel_path"], "kernel_sha256": source["kernel_sha256"]})

        def diagnostic(event: dict) -> None:
            require(isinstance(event, dict), "Kernel diagnostic callback payload must be a dict")
            append_diagnostic(event)

        derive_entry_calls = 1
        derive_calls = 1
        component_result = module.derive(source, diagnostic)
        require(isinstance(component_result, dict), "Kernel derive result must be a plain JSON object")
    except BaseException as error:
        if kernel_imported:
            append_diagnostic({"event": "manager_kernel_or_derive_failed", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "error_type": type(error).__name__})
        else:
            append_diagnostic({"event": "manager_kernel_or_derive_failed", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "error_type": type(error).__name__})
        component_result = failure_result(error, derive_entry_calls, derive_calls, kernel_import_calls)
    durable_result = clone(component_result)
    durable_result["component_result"] = clone(component_result)
    durable_result["component_result_exact"] = True
    durable_result["execution_metadata"] = {"kernel_import_calls": kernel_import_calls, "kernel_imported": kernel_imported, "derive_entry_calls": derive_entry_calls, "derive_calls": derive_calls, "runner_executed": False, "math_reexecuted": False, "old_math_reexecuted": False}
    durable_result["kernel_import_calls"] = kernel_import_calls
    durable_result["kernel_imported"] = kernel_imported
    durable_result["derive_entry_calls"] = derive_entry_calls
    durable_result["derive_calls"] = derive_calls
    durable_result["runner_executed"] = False
    result_raw = write_new(RESULT_CHECKPOINT, durable_result)
    diagnostic_summary_value = diagnostic_summary()
    validation = validate_result(durable_result, source, derive_entry_calls, derive_calls, diagnostic_summary_value)
    receipt = {
        "schema": COMPONENT + "_run_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_path(ENTERED_CHECKPOINT),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "diagnostic_summary": diagnostic_summary_value,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "validation": validation,
        "raw_result_saved_before_validation": True,
        "component_result_preserved": True,
        "rerun": False,
        "runner_executed": False,
        "runtime_model_metadata": durable_result.get("runtime_model_metadata"),
    }
    receipt_raw = write_new(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("action", choices=("prepare", "run", "publish", "resume", "seal"))
    value.add_argument("--parent-report")
    value.add_argument("--parent-report-sha256")
    value.add_argument("--source-facts")
    value.add_argument("--source-facts-sha256")
    value.add_argument("--alias-trace")
    value.add_argument("--design-review")
    value.add_argument("--review")
    value.add_argument("--review-sha256")
    value.add_argument("--kernel-api")
    value.add_argument("--kernel-path")
    value.add_argument("--runner-path")
    value.add_argument("--kernel-sha256")
    value.add_argument("--runner-sha256")
    return value


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    if args.action == "prepare":
        value = prepare(args)
    elif args.action == "run":
        value = run(args)
    elif args.action == "publish":
        value = publish(args)
    elif args.action == "resume":
        value = resume(args)
    else:
        value = seal(args)
    print(json.dumps(value, ensure_ascii=True, sort_keys=True))


if __name__ == "__main__":
    main()
