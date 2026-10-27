"""Metadata-only draft for the next native history/radius definition bundle.

This file is deliberately an external WORK manager.  Importing it does not
import project modules, decode protected archives, call a source provider,
evaluate a jet, run a kernel, or perform mathematics.  The lifecycle entry
points perform only the explicitly requested metadata workflow once the root
agent has supplied and frozen the static trace packet, saved parameter-handle
packet, and provider contract.

The draft records the workflow shape needed by that later consumer:

* nine source definition quantities (``V``, ``m``, ``e``, ``P0``, ``p``,
  ``E``, ``h``, ``k``, and ``pressure``) across background and modified
  branches, plus one native radius definition, each carrying a function/jet
  identifier and frame;
* five history channels are data records, not counters;
* native and rebuilt ``Pstar``/``delta``/radius records stay separate;
* all four input checkpoints, two code artifacts, the report, document, and
  four handoff artifacts are represented losslessly;
* the current saved review is supplied at publication/seal time, never copied
  from the parent review constant.

No result is produced by importing this module.  The provider is loaded only
by the explicit ``run`` lifecycle action after durable source and entry
checkpoints have been written; ``prepare`` and all validators remain
metadata-only.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import base64
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


REPO = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
WORK = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")
COMPONENT = "native_history_radius_bundle"

# These are authenticated predecessor facts supplied by the current brief.
# They are provenance anchors only; they are not provider inputs or computed
# source values.
PARENT_COMPONENT = "native_W_parameter_interface"
# The new bundle's parent is the published 23322... state.  The parent
# report itself records 4707... as its internal parent lineage; keep that
# audit field distinct and never relabel the new bundle's parent.
PARENT_COMMIT = "23322d9553753394d0c1be6e085c5d54fcce2c5f"
PARENT_REPORT_INTERNAL_PARENT_COMMIT = "4707fc262d51d33677dc5d53fb6f1bb1e6fb53af"
CURRENT_HEAD_COMMIT = "23322d9553753394d0c1be6e085c5d54fcce2c5f"
PARENT_REPORT = REPO / "experiments/root_st073/lei_ren_part1_paper_compliant_current_native_W_parameter_interface.json"
PARENT_REPORT_SHA256 = "b817a0f0fca623bd40ec8ab074b1e7c851440f5c2b59a934c2ef5e2b3dce9c26"
PARENT_RESULT = REPO / "tmp/native_W_parameter_interface/result.json"
PARENT_RESULT_SHA256 = "7d45286f9a4c7a619f2edf06f822ae709dcf4e2175ca001bae7bd7455d7d3877"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "1cb0ae67e0c6bfe2d8ddc26a8f77757db19d0903bcdf7fc0e62093698923c460"

PARENT_SOURCE_PIN_COUNT = 2254
PARENT_CHECKPOINT_COUNT = 260
PARENT_FAILED_ARTIFACT_COUNT = 183
PARENT_PROTECTED_FILE_COUNT = 16
PARENT_DIRTY_PATH_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

# The four additional checkpoint names and the failure ledger increment are a
# planning fact from the current brief.  They are not acceptance counts until
# the root API contract and current review have been frozen.
PLANNED_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
PLANNED_FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_ARTIFACT_COUNT + 7

HISTORY_TRACE_V1 = WORK / "native_history_bundle_source_trace_2026_10_27.json"
HISTORY_TRACE_V1_SHA256 = "4c74e9ce767ed5d4e6d896d5bd411d5b137e8aa534058a8275cda95769e131ae"
HISTORY_TRACE = WORK / "native_history_bundle_source_trace_2026_10_27_v4.json"
HISTORY_TRACE_SHA256 = "4192b97d6d4135b3d5d3f3ff69c86982baa1b895b06412a188664a189420c114"
RADIUS_TRACE_V1 = WORK / "native_radius_frame_source_trace_2026_10_27.json"
RADIUS_TRACE_V1_SHA256 = "08b0fdbdef53b890e08fdb2c694be7e11ba51f62fd99cef66215e2676dd792ac"
RADIUS_TRACE = WORK / "native_radius_frame_source_trace_2026_10_27_v2.json"
RADIUS_TRACE_SHA256 = "1ecae6e29c77d95ae21545ff20bbfc6db69f3e5207e3fdf5cace1db023f400ee"
SAVED_PARAMETER_RESULT_SHA256 = PARENT_RESULT_SHA256

CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (SOURCE_CHECKPOINT, ENTERED_CHECKPOINT, RESULT_CHECKPOINT, DIAGNOSTIC_CHECKPOINT)
INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_27.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_27.json"
REPORT_RELATIVE = Path("experiments/root_st073/lei_ren_part1_paper_compliant_current_native_history_radius_bundle.json")
MILESTONE_RELATIVE = Path("docs/CURRENT_NATIVE_HISTORY_RADIUS_BUNDLE_2026_10_27.md")
HANDOFF_RELATIVES = (
    Path("docs/AGENT_TASKS.md"),
    Path("docs/CROSS_DEVICE_HANDOFF.md"),
    Path("docs/CURRENT_CHECKPOINT.md"),
    Path("docs/NEXT_TASKS.md"),
)
DEFAULT_KERNEL_PATH = REPO / "src/openai_ns_reconstruction/native_history_radius_bundle.py"
DEFAULT_RUNNER_PATH = REPO / "scripts/run_native_history_radius_bundle.py"
EXPECTED_KERNEL_SHA256 = "d35f2478c202f0017dd5051ed3a2115f0d31fbfe9c60e52b710ac2c3be72c49a"
EXPECTED_KERNEL_API_SHA256 = "68c882439e3415fa1e1df6f453ecac211a4b9999e0f1468a5ae5513f75a4ce8c"
DIAGNOSTIC_SCHEMA = COMPONENT + "_diagnostic_v1"
ENTERED_SCHEMA = COMPONENT + "_entered_v1"

INPUT_SCHEMA = "native_history_radius_bundle_input_v1"
TRACE_SCHEMA = "native_functionjet_row_convention_and_provider_refs_v2"
TRACE_SCHEMA_V1 = "native_history_bundle_source_trace_v1"
RADIUS_TRACE_SCHEMA = "native_radius_frame_source_trace_v1"
RESULT_SCHEMA = "original_native_history_radius_bundle"
REPORT_SCHEMA = "original_native_history_radius_bundle_publication_v1"
SEAL_SCHEMA = "native_history_radius_bundle_final_seal_v1"

ROOT_API_STATUS = "ready_for_root_freeze"
ROOT_API_CANDIDATE_SIGNATURE = "derive(request: dict, callback: Callable | None = None) -> dict"

DEFINITION_KEYS = ("V", "m", "e", "P0", "p", "E", "h", "k", "pressure")
DEFINITION_SLOT_KEYS = (*DEFINITION_KEYS, "radius")
HISTORY_KEYS = ("m", "h", "k", "e", "p")
BRANCH_KEYS = ("background", "modified")
PARAMETER_SIDE_KEYS = ("native", "rebuilt")
PARAMETER_KEYS = ("Pstar", "delta", "radius")
ORDERS = (0, 1, 2, 3, 4)
NATIVE_RADIUS_NODE = 241787
NATIVE_COORDINATE_NODE = 241536
NATIVE_OFFSET_NODE = 241542
NATIVE_OFFSET_LITERAL = "17 + exp(40) + x"
NATIVE_DOMAIN = (0, 1)
P2_BACKGROUND_ROOT_IDS = (242304, 242308, 242314, 242321)
P2_MODIFIED_ROOT_IDS = (249570, 249574, 249580, 249587, 249595)

EXPECTED_DEFINITION_RECORD_COUNT = len(DEFINITION_KEYS) * len(BRANCH_KEYS) * len(ORDERS)
EXPECTED_RADIUS_RECORD_COUNT = 1
EXPECTED_BRANCH_BUNDLE_COUNT = 2
EXPECTED_REJECTION_COUNT = 3
REJECTION_KEYS = (
    "mixed_branch_bundle",
    "rebuilt_radius_frame",
    "duplicate_factorial_request",
)

EXPECTED_COUNTS = {
    "definition_provider_entries": 1,
    "native_branches_served": 2,
    "field_definition_queries": 90,
    "radius_definition_queries": 1,
    "native_bundle_manifests": 2,
    "mixed_branch_rejections": 1,
    "rebuilt_frame_rejections": 1,
    "duplicate_normalization_rejections": 1,
    "new_mathematical_proofs": 0,
    "new_algebra_family_calls": 0,
    "numeric_calls": 0,
    "source_value_provider_calls": 0,
    "old_math_calls": 0,
}
COUNT_KEYS = tuple(EXPECTED_COUNTS)
PROVISIONAL_COUNT_KEYS = COUNT_KEYS
PROVISIONAL_COUNTS = EXPECTED_COUNTS
EXPECTED_FLAGS = {
    "native_history_radius_definitions_served": True,
    "native_branch_frame_separation_enforced": True,
    "native_ordinary_Z_projection_preserved": True,
}
FLAG_KEYS = tuple(EXPECTED_FLAGS)
CALLBACK_EVENT_COUNTS = {
    "source_validation": 1,
    "definition_provider_entry": 1,
    "field_definition_row": EXPECTED_DEFINITION_RECORD_COUNT,
    "radius_definition_row": EXPECTED_RADIUS_RECORD_COUNT,
    "native_bundle_manifest": EXPECTED_BRANCH_BUNDLE_COUNT,
    "context_rejection": EXPECTED_REJECTION_COUNT,
}
CALLBACK_EVENT_COUNT = sum(CALLBACK_EVENT_COUNTS.values())
MANAGER_PREFIX_EVENT_COUNT = 2
JOURNAL_EVENT_COUNT = MANAGER_PREFIX_EVENT_COUNT + CALLBACK_EVENT_COUNT
FINAL_COUNTS = EXPECTED_COUNTS
FINAL_FLAGS = EXPECTED_FLAGS

# Nine report-scope attributes, two code artifacts, the publication document,
# and four handoff artifacts are intentionally named here so the eventual
# report cannot silently flatten or omit them.
REPORT_SCOPE_ATTRIBUTES = (
    "static_trace_packet",
    "saved_parameter_handles",
    "definition_records",
    "history_channel_records",
    "parameter_separation",
    "provider_contract",
    "provisional_counts",
    "ledger_preservation",
    "current_review_and_handoffs",
)
CODE2_ARTIFACT_KEYS = ("kernel", "runner")
HANDOFF_KEYS = ("agent_tasks", "cross_device_handoff", "current_checkpoint", "next_tasks")
INPUT_CHECKPOINT_KEYS = ("source", "entered", "result", "diagnostics")

NONCLAIM_FLAGS = (
    "source_value_computation",
    "native_rebuilt_Pstar_equality",
    "native_rebuilt_delta_equality",
    "native_rebuilt_radius_equality",
    "zero_claim",
    "error_cap_deletion",
    "closure",
    "temporal_recursion",
    "physical_completion",
    "numeric_containment",
    "old_math_reentry",
    "rebuilt_template_instantiation",
)


def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def clone(value: Any) -> Any:
    return deepcopy(value)


def encode(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8") + b"\n"


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def wire_path(path: str | Path) -> str:
    return Path(path).as_posix()


def resolve_path(value: str | Path, base: Path = WORK) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (base / candidate).resolve()


def read_bytes(path: str | Path, expected_sha256: str | None = None) -> bytes:
    raw = Path(path).read_bytes()
    if expected_sha256 is not None:
        require(sha_bytes(raw) == expected_sha256.lower(), f"SHA-256 changed: {path}")
    return raw


def sha_path(path: str | Path) -> str:
    return sha_bytes(Path(path).read_bytes())


def portable(value: Any) -> Any:
    """Return a JSON-safe metadata value without importing project helpers."""

    return json.loads(json.dumps(value, ensure_ascii=True, allow_nan=False))


def read_json(path: str | Path, expected_sha256: str | None = None) -> tuple[bytes, dict[str, Any]]:
    raw = read_bytes(path, expected_sha256)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"JSON object required: {path}")
    return raw, value


def write_new(path: str | Path, value: Any) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_replace(path: str | Path, value: Any) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("wb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def append_jsonl(path: str | Path, value: Mapping[str, Any]) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = encode(dict(value))
    with path.open("ab") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def brief_safe(value: Any) -> Any:
    """Strip only large callback/base64 rows when making a BRIEF view."""

    if isinstance(value, dict):
        return {
            key: brief_safe(child)
            for key, child in value.items()
            if key not in {"raw_bytes_base64", "callback_event_raw_rows"}
        }
    if isinstance(value, list):
        return [brief_safe(child) for child in value]
    return value


def artifact_entry(
    path: str | Path,
    *,
    expected_sha256: str | None = None,
    embed: bool = True,
    outside_science_pins: bool = True,
) -> dict[str, Any]:
    """Return a lossless artifact record without interpreting its contents."""

    path = Path(path)
    raw = read_bytes(path, expected_sha256)
    entry: dict[str, Any] = {
        "path": wire_path(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
    }
    if outside_science_pins:
        entry["outside_science_pins"] = True
    if embed:
        entry["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
    return entry


def load_current_brief(path: str | Path = MAIN_BRIEF) -> dict[str, Any]:
    """Authenticate and return the current brief without writing it."""

    raw, brief = read_json(path, MAIN_BRIEF_SHA256)
    require(len(raw) == 83915, "Current BRIEF byte count changed")
    rules = brief.get("execution_rules")
    require(isinstance(rules, list) and len(rules) == 7, "Current execution_rules are not the restored seven-item list")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "Goal status changed")
    require(brief.get("metadata_only") is True and brief.get("source_only") is True, "Current brief is not metadata/source-only")
    return {
        "path": wire_path(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "value": clone(brief),
        "execution_rules": clone(rules),
        "execution_rules_preserved_exact": True,
    }


def build_brief_preimage_record(current_brief: Mapping[str, Any]) -> dict[str, Any]:
    """Capture the exact sub-100 KB BRIEF preimage before any update.

    The full bytes are retained here.  ``brief_safe`` is reserved for a
    compact decoded view and is never used for this archive record.
    """

    path = current_brief.get("path")
    expected = current_brief.get("sha256")
    require(path and expected, "Current BRIEF pointer is incomplete")
    raw = read_bytes(path, expected)
    require(len(raw) < 100_000, "Current BRIEF preimage exceeds the 100 KB limit")
    return {
        "path": wire_path(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "preserved_exact_bytes": True,
        "raw_bytes_base64": base64.b64encode(raw).decode("ascii"),
        "decoded_view_policy": "brief_safe may strip only raw_bytes_base64 and callback_event_raw_rows from a separate BRIEF view",
    }


def load_parent_authority(
    report_path: str | Path = PARENT_REPORT,
    result_path: str | Path = PARENT_RESULT,
) -> dict[str, Any]:
    """Load predecessor metadata and all retained ledgers, never bare result data."""

    report_raw, report = read_json(report_path, PARENT_REPORT_SHA256)
    result_raw, result = read_json(result_path, PARENT_RESULT_SHA256)
    require(report.get("component") == PARENT_COMPONENT, "Parent component changed")
    require(report.get("passed") is True and report.get("metadata_only") is True and report.get("source_only") is True, "Parent authority markers changed")
    require(report.get("parent_math_failure_preserved") is True, "Parent failed math family was not preserved")
    require(report.get("parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT and report.get("parent_head") == PARENT_REPORT_INTERNAL_PARENT_COMMIT, "Parent report internal lineage changed")
    require(result.get("component") == PARENT_COMPONENT and result.get("passed") is True, "Parent result identity changed")
    require(report.get("full_result") == result or report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256, "Parent full result is not pinned")

    source_pins = report.get("source_pins")
    require(isinstance(source_pins, dict) and len(source_pins) == PARENT_SOURCE_PIN_COUNT, "Parent source pins changed")
    require(all(valid_sha(digest) for digest in source_pins.values()), "Parent source pin digest malformed")
    failed_hashes = report.get("preserved_failed_artifact_hashes_union")
    protected_hashes = report.get("protected_file_hashes")
    dirty_hashes = report.get("legacy_dirty_path_hashes")
    require(isinstance(failed_hashes, dict) and len(failed_hashes) == PARENT_FAILED_ARTIFACT_COUNT, "Parent failed-artifact ledger changed")
    require(isinstance(protected_hashes, dict) and len(protected_hashes) == PARENT_PROTECTED_FILE_COUNT, "Parent protected-file ledger changed")
    require(isinstance(dirty_hashes, dict) and len(dirty_hashes) == PARENT_DIRTY_PATH_COUNT, "Parent dirty-path ledger changed")

    actual = report.get("actual_counts")
    require(isinstance(actual, dict), "Parent actual counts missing")
    required_counts = {
        "source_binding_calls": 1,
        "native_admissions": 2,
        "rejections": 2,
        "template_signatures": 1,
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
    }
    require(actual == required_counts, "Parent actual counts changed")

    retained_maps: dict[str, Any] = {}
    for key in ("preserved_error_cap_metadata", "retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance"):
        require(isinstance(report.get(key), dict), f"Parent retained map missing: {key}")
        retained_maps[key] = clone(report[key])

    false_flags = report.get("parent_false_flags")
    false_values = report.get("parent_false_flag_values")
    require(isinstance(false_flags, list) and len(false_flags) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag ledger changed")
    require(isinstance(false_values, dict) and false_values == {key: False for key in false_flags}, "Parent false-flag values changed")

    return {
        "report": {
            "path": wire_path(report_path),
            "sha256": sha_bytes(report_raw),
            "bytes": len(report_raw),
            "schema": report.get("schema"),
            "component": report.get("component"),
        },
        "result": {
            "path": wire_path(result_path),
            "sha256": sha_bytes(result_raw),
            "bytes": len(result_raw),
            "schema": result.get("schema"),
            "component": result.get("component"),
        },
        "report_sha256": sha_bytes(report_raw),
        "result_sha256": sha_bytes(result_raw),
        "commit": PARENT_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "current_head_commit": CURRENT_HEAD_COMMIT,
        "actual_counts": clone(actual),
        "source_pins": clone(source_pins),
        "failed_artifact_hashes": clone(failed_hashes),
        "protected_file_hashes": clone(protected_hashes),
        "dirty_path_hashes": clone(dirty_hashes),
        "source_pin_count": report.get("source_pin_count", 2254),
        "checkpoint_count": report.get("checkpoint_count", PARENT_CHECKPOINT_COUNT),
        "preserved_failed_artifact_count": report.get("preserved_failed_artifact_count", PARENT_FAILED_ARTIFACT_COUNT),
        "parent_protected_file_count": report.get("parent_protected_file_count", PARENT_PROTECTED_FILE_COUNT),
        "parent_dirty_path_count": report.get("parent_dirty_path_count", PARENT_DIRTY_PATH_COUNT),
        "expected_parent_false_flag_count": report.get("expected_parent_false_flag_count", PARENT_FALSE_FLAG_COUNT),
        "parent_false_flags": clone(false_flags),
        "parent_false_flag_values": clone(false_values),
        "parent_math_failure_preserved": True,
        "retained_maps": retained_maps,
        "report_view": {
            "schema": report.get("schema"),
            "component": report.get("component"),
            "passed": report.get("passed"),
            "metadata_only": report.get("metadata_only"),
            "source_only": report.get("source_only"),
            "actual_counts": clone(actual),
            "full_goal_complete": report.get("full_goal_complete"),
            "parent_report_is_authenticated_view": True,
            "complete_parent_report_raw_sha256": PARENT_REPORT_SHA256,
        },
    }


def _definition_identifier(record: Mapping[str, Any]) -> Any:
    for key in ("function_id", "jet_definition_id", "node_id", "definition_id"):
        if key in record:
            return record[key]
    return None


def validate_definition_record(record: Mapping[str, Any], *, side: str, quantity: str) -> dict[str, Any]:
    """Validate an identifier/frame record without reading or evaluating its source."""

    require(isinstance(record, Mapping), f"Definition record is not an object: {side}.{quantity}")
    identifier = _definition_identifier(record)
    require(identifier is not None, f"Definition identifier missing: {side}.{quantity}")
    require(record.get("frame") is not None, f"Definition frame missing: {side}.{quantity}")
    require(record.get("source_ref") or record.get("source_path"), f"Definition source reference missing: {side}.{quantity}")
    if "source_sha256" in record:
        require(isinstance(record["source_sha256"], str) and len(record["source_sha256"]) == 64, f"Definition source SHA malformed: {side}.{quantity}")
    if "value" in record:
        raise AssertionError(f"Definition record contains a source value: {side}.{quantity}")
    return {
        "side": side,
        "quantity": quantity,
        "identifier": identifier,
        "frame": clone(record["frame"]),
        "source_ref": record.get("source_ref", record.get("source_path")),
        "source_sha256": record.get("source_sha256"),
        "jet_order": record.get("jet_order"),
        "record": clone(dict(record)),
    }


def _native_frame(order: int, branch: str) -> dict[str, Any]:
    return {
        "chart": "O3_transition",
        "branch": branch,
        "coordinate_node": NATIVE_COORDINATE_NODE,
        "coordinate_symbol": "true_outer_O3_transition",
        "ordinary_Z_order": order,
        "domain": list(NATIVE_DOMAIN),
        "offset_node": NATIVE_OFFSET_NODE,
        "offset_literal": NATIVE_OFFSET_LITERAL,
        "physical_measure": "1",
        "physical_log_radius_frame": True,
        "ordinary_Z_is_separate_from_chart_coordinate": True,
    }


def _native_node_ids(packet: Mapping[str, Any], branch: str, quantity: str) -> list[int]:
    if packet.get("schema") == TRACE_SCHEMA:
        selectors = packet["executable_native_definition_provider_plan"]["branch_field_selectors"][branch]
        nodes = selectors.get(quantity)
        require(isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes)), f"Native v2 node rows missing: {branch}.{quantity}")
        require(len(nodes) >= len(ORDERS), f"Native v2 node rows too short: {branch}.{quantity}")
        return [int(nodes[order]) for order in ORDERS]
    handles = packet["native_field_function_handles"][branch]["fields"]
    if quantity in {"E", "h", "k", "pressure"}:
        record = handles["auxiliary_required_by_generic_source"][quantity]
    elif quantity == "P0":
        projection = packet["P0_ordinary_projection"]["ordinary_rows"]
        by_order = {int(item["order"]): int(item["node"]) for item in projection}
        return [by_order[order] for order in ORDERS]
    else:
        record = handles[quantity]
    nodes = record.get("node_ids") if isinstance(record, Mapping) else record
    require(isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes)), f"Native node rows missing: {branch}.{quantity}")
    require(len(nodes) >= len(ORDERS), f"Native node rows too short: {branch}.{quantity}")
    return [int(nodes[order]) for order in ORDERS]


def _native_definition_record(packet: Mapping[str, Any], branch: str, quantity: str, order: int) -> dict[str, Any]:
    nodes = _native_node_ids(packet, branch, quantity)
    node = nodes[order]
    trace_sha256 = HISTORY_TRACE_SHA256 if packet.get("schema") == TRACE_SCHEMA else HISTORY_TRACE_V1_SHA256
    record: dict[str, Any] = {
        "side": branch,
        "branch": branch,
        "quantity": quantity,
        "order": order,
        "node_id": node,
        "definition_kind": "saved_FunctionJet_row_selector",
        "frame": _native_frame(order, branch),
        "source_ref": f"native_field_function_handles.{branch}.{quantity}.node_ids[{order}]",
        "trace_packet_sha256": trace_sha256,
        "value_materialized": False,
        "source_value_present": False,
        "ordinary_Z_factorial_applied_once": False,
    }
    if quantity == "P0":
        if packet.get("schema") == TRACE_SCHEMA:
            guard = packet["P0_factorial_once_guard"]["saved_ordinary_nodes"][f"order_{order}"]
            record.update({
                "source_ref": f"P0_factorial_once_guard.saved_ordinary_nodes.order_{order}",
                "factorial": int(guard["factorial_already_applied"]),
                "ordinary_Z_factorial_applied_once": True,
                "factorial_source_projection": packet["P0_factorial_once_guard"]["source_projection"].get("order_4_literal"),
                "duplicate_factorial_request_rejected": True,
            })
        else:
            projection = packet["P0_ordinary_projection"]["ordinary_rows"][order]
            record.update({
                "source_ref": f"P0_ordinary_projection.ordinary_rows[{order}]",
                "factorial": int(projection["factorial"]),
                "ordinary_Z_factorial_applied_once": True,
                "factorial_source_projection": projection.get("source_projection"),
                "duplicate_factorial_request_rejected": True,
            })
    if quantity == "V":
        record.update({
            "source_literal": "V=0",
            "zero_function_handle": True,
            "zero_identity_claimed": False,
        })
    return record


def extract_native_definition_bundle(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the 90 native definition references from the static packet."""

    require(packet.get("schema") in {TRACE_SCHEMA, TRACE_SCHEMA_V1}, "Native history trace schema changed")
    definitions: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for branch in BRANCH_KEYS:
        definitions[branch] = {}
        for quantity in DEFINITION_KEYS:
            definitions[branch][quantity] = [
                _native_definition_record(packet, branch, quantity, order)
                for order in ORDERS
            ]

    chart = packet.get("native_chart_signature", {})
    frame = chart.get("frame", {})
    if packet.get("schema") == TRACE_SCHEMA:
        raw_radius = packet.get("native_radius_raw_record", {})
        chart_name = raw_radius.get("chart", "O3_transition")
        radius_node = raw_radius.get("radius_node", NATIVE_RADIUS_NODE)
        offset_node = raw_radius.get("offset_node", NATIVE_OFFSET_NODE)
        radius_literal = raw_radius.get("radius_source_literal", "R=Rm*exp(offset)")
        radius_owner = raw_radius.get("radius_owner", "self.Rm=self.leading.inlet['radius']")
        radius_domain = list(NATIVE_DOMAIN)
    else:
        chart_name = chart.get("chart")
        radius_node = frame.get("radius_node", NATIVE_RADIUS_NODE)
        offset_node = frame.get("offset_node", NATIVE_OFFSET_NODE)
        radius_literal = frame.get("radius_literal", "R=Rm*exp(offset)")
        radius_owner = "self.Rm=self.leading.inlet['radius']"
        radius_domain = clone(chart.get("domain", {}).get("semantic", list(NATIVE_DOMAIN)))
    radius_record = {
        "side": "native",
        "branch": "shared_native_chart",
        "quantity": "radius",
        "order": None,
        "node_id": int(radius_node),
        "definition_kind": "native_radius_function_definition",
        "frame": {
            "chart": chart_name,
            "coordinate_node": chart.get("coordinate", {}).get("node", NATIVE_COORDINATE_NODE) if chart else NATIVE_COORDINATE_NODE,
            "domain": radius_domain,
            "offset_node": offset_node,
            "offset_literal": NATIVE_OFFSET_LITERAL,
            "physical_measure": frame.get("physical_measure_node", 1) if frame else "1",
            "radius_literal": radius_literal,
        },
        "source_ref": "native_radius_raw_record.radius_node" if packet.get("schema") == TRACE_SCHEMA else "native_chart_signature.frame.radius_node",
        "radius_owner": radius_owner,
        "trace_packet_sha256": HISTORY_TRACE_SHA256 if packet.get("schema") == TRACE_SCHEMA else HISTORY_TRACE_V1_SHA256,
        "value_materialized": False,
        "source_value_present": False,
        "rebuilt_radius_alias_allowed": False,
    }
    normalization = packet.get("history_normalization_and_Rbase", {})
    p2 = packet.get("P2_branch_and_inputs", {})
    if packet.get("schema") == TRACE_SCHEMA:
        history_rows: dict[str, Any] = {}
        transport_rules: dict[str, Any] = {}
        p2_routes = {
            "background_generic": list(P2_BACKGROUND_ROOT_IDS),
            "modified_generic": list(P2_MODIFIED_ROOT_IDS),
            "mixed_branch_bundle_rejected": True,
        }
        factorial_guard = {
            "P0_factorial_applied_once": packet.get("P0_factorial_once_guard", {}).get("double_factorial_guard") is True,
            "duplicate_factorial_request_rejected": packet.get("P0_factorial_once_guard", {}).get("double_factorial_guard") is True,
        }
    else:
        history_rows = {
            "O2_buffer_incoming": clone(normalization.get("O2_buffer_endpoint_history_nodes", {})),
            "background": clone(normalization.get("background_history_rows", {})),
            "finite_N": clone(normalization.get("finite_N_history_rows", {})),
        }
        transport_rules = clone(normalization.get("transport_rules", {}))
        p2_routes = {
            "background_generic": clone(p2.get("selected_p2_handles", list(P2_BACKGROUND_ROOT_IDS))),
            "modified_generic": clone(p2.get("separate_modified_generic_p2_handles", list(P2_MODIFIED_ROOT_IDS))),
            "mixed_branch_bundle_rejected": True,
        }
        factorial_guard = {
            "P0_factorial_applied_once": packet.get("P0_ordinary_projection", {}).get("factorial_applied_once") is True,
            "duplicate_factorial_request_rejected": packet.get("P0_ordinary_projection", {}).get("do_not_multiply_saved_ordinary_P0_nodes_again") is True,
        }
    return {
        "schema": packet.get("schema"),
        "status": packet.get("status"),
        "definition_records": definitions,
        "radius_definition": radius_record,
        "definition_record_count": EXPECTED_DEFINITION_RECORD_COUNT,
        "radius_definition_count": EXPECTED_RADIUS_RECORD_COUNT,
        "branch_bundle_count": EXPECTED_BRANCH_BUNDLE_COUNT,
        "history_channels": list(HISTORY_KEYS),
        "history_rows": history_rows,
        "transport_rules": transport_rules,
        "P2_routes": p2_routes,
        "mixed_branch_bundle_rejected": True,
        "factorial_guard": factorial_guard,
        "provider_definition_plan": clone(packet.get("frozen_provider_definition_plan", {})),
        "raw_packet": clone(dict(packet)),
    }


def validate_native_history_trace(packet: Mapping[str, Any]) -> dict[str, Any]:
    require(packet.get("schema") in {TRACE_SCHEMA, TRACE_SCHEMA_V1}, "Native history trace schema changed")
    scope = packet.get("scope", {})
    if packet.get("schema") == TRACE_SCHEMA:
        require(packet.get("packet_type") in {"native_functionjet_row_convention_and_provider_refs_v2", "native_history_bundle_source_trace_v2_addendum"}, "Native v2 packet type changed")
        base_packet = packet.get("base_packet", {})
        base_path = base_packet.get("path")
        require(base_path in {wire_path(HISTORY_TRACE_V1), f"work/{HISTORY_TRACE_V1.name}"} and base_packet.get("sha256") == HISTORY_TRACE_V1_SHA256 and base_packet.get("preserved_unchanged") is True, "Native v2 base packet pin changed")
        require(scope.get("source_imported") is False and scope.get("source_executed") is False, "Native history source was imported or executed")
        require(scope.get("graph_evaluated") is False and scope.get("numeric_or_integral_evaluation") is False, "Native history graph or values were evaluated")
        require(scope.get("math_identity_check") is False and scope.get("new_provider_called") is False, "Native history consumer was entered")
        plan = packet.get("executable_native_definition_provider_plan", {})
        require(plan.get("serves_branches") == list(BRANCH_KEYS), "Native branches changed")
        require(plan.get("explicit_branch_argument_required") is True and plan.get("branch_mixing_forbidden") is True, "Branch guards changed")
        require(plan.get("all_nine_generic_inputs") == list(DEFINITION_KEYS), "Native quantity list changed")
        require(plan.get("ordinary_orders") == list(ORDERS), "Native ordinary orders changed")
        reference_count = plan.get("reference_count", {})
        require(reference_count.get("per_branch") == 45 and reference_count.get("branches") == 2 and reference_count.get("total") == EXPECTED_DEFINITION_RECORD_COUNT, "Native definition count changed")
        require(reference_count.get("radius_references") == EXPECTED_RADIUS_RECORD_COUNT, "Native radius count changed")
        raw_radius = packet.get("native_radius_raw_record", {})
        require(raw_radius.get("radius_node") == NATIVE_RADIUS_NODE, "Native radius node changed")
        require(raw_radius.get("radius_source_literal") == "R=Rm*exp(offset)", "Native radius literal changed")
        require(raw_radius.get("radius_owner") == "self.Rm=self.leading.inlet['radius']", "Native radius owner changed")
        require(raw_radius.get("radius_alias_status", {}).get("native_to_rebuilt") is False, "Native radius alias changed")
        zero = packet.get("native_zero_source", {})
        # v2 keeps the graph node as a structured source record.  Preserve
        # that record and its explicit non-padding marker; do not reduce it to
        # a bare node id or infer a value.
        graph_node = zero.get("graph_node", {})
        require(isinstance(graph_node, Mapping), "Native zero graph-node record changed")
        require(graph_node.get("node") == 0, "Native zero graph node changed")
        require(graph_node.get("raw_record", {}).get("operation") == "exact_rational", "Native zero source operation changed")
        require(graph_node.get("raw_record", {}).get("numerator") == 0 and graph_node.get("raw_record", {}).get("denominator") == 1, "Native zero rational record changed")
        require(graph_node.get("not_padding") is True and graph_node.get("not_numeric_evaluation") is True, "Native zero non-padding marker changed")
        guard = packet.get("P0_factorial_once_guard", {})
        require(guard.get("double_factorial_guard") is True, "P0 duplicate-factorial guard changed")
    else:
        require(scope.get("static_only") is True, "Native history trace is not static-only")
        require(scope.get("source_imported") is False and scope.get("source_executed") is False, "Native history source was imported or executed")
        require(scope.get("graph_evaluated") is False and scope.get("numeric_or_integral_evaluation") is False, "Native history graph or values were evaluated")
        require(scope.get("provider_executed") is False and scope.get("source_consumer_executed") is False, "Native history provider was executed")
        chart = packet.get("native_chart_signature", {})
        require(chart.get("chart") == "O3_transition", "Native chart changed")
        require(chart.get("coordinate", {}).get("node") == NATIVE_COORDINATE_NODE, "Native coordinate node changed")
        require(chart.get("domain", {}).get("semantic") == ["0", "1"], "Native chart domain changed")
        frame = chart.get("frame", {})
        require(frame.get("radius_node") == NATIVE_RADIUS_NODE, "Native radius node changed")
        require(frame.get("physical_measure_node") == 1, "Native physical measure metadata changed")
        p2 = packet.get("P2_branch_and_inputs", {})
        require(tuple(p2.get("selected_p2_handles", ())) == P2_BACKGROUND_ROOT_IDS, "Background P2 route changed")
        require(tuple(p2.get("separate_modified_generic_p2_handles", ())) == P2_MODIFIED_ROOT_IDS, "Modified P2 route changed")
        require(p2.get("route_distinction"), "P2 branch distinction missing")
        finite = packet.get("finite_N_literal_source", {})
        require(finite.get("N_materialized") is False and finite.get("J_node18_materialized", False) is False, "Finite-N literal was materialized")
    bundle = extract_native_definition_bundle(packet)
    require(bundle["definition_record_count"] == EXPECTED_DEFINITION_RECORD_COUNT, "Native definition count changed")
    require(bundle["radius_definition_count"] == EXPECTED_RADIUS_RECORD_COUNT, "Native radius count changed")
    bundle["rejection_contract"] = validate_branch_rejection_contract(bundle)
    return bundle


def validate_radius_trace(packet: Mapping[str, Any]) -> dict[str, Any]:
    require(packet.get("schema") == RADIUS_TRACE_SCHEMA, "Native radius trace schema changed")
    boundary = packet.get("execution_boundary", {})
    require(boundary.get("source_only") is True and boundary.get("metadata_only") is True, "Radius trace markers changed")
    require(boundary.get("runtime_model_metadata") is None, "Runtime model metadata is not null")
    for key in ("imports_of_project_modules", "provider_or_helper_calls", "old_math_calls", "new_math_or_proof_calls", "numeric_calls", "new_algebra_family_calls", "source_provider_calls"):
        require(boundary.get(key) == 0, f"Radius trace execution count changed: {key}")
    require(boundary.get("raw_gzip_decoded") is False, "Protected archive was decoded")
    frame = packet.get("native_chart_frame", {})
    require(frame.get("chart") == "O3_transition", "Radius chart changed")
    require(frame.get("coordinate_node") == NATIVE_COORDINATE_NODE, "Radius coordinate changed")
    require(frame.get("O3_transition_domain") == list(NATIVE_DOMAIN), "Radius domain changed")
    require(frame.get("O3_transition_log_offset") == NATIVE_OFFSET_LITERAL, "Radius offset changed")
    require(frame.get("O3_transition_physical_measure") == "1", "Radius physical measure metadata changed")
    require(frame.get("radius_definition") == "R_native = self.Rm * exp(offset)", "Radius definition changed")
    require(frame.get("native_Rm_owner") == "self.Rm = self.leading.inlet['radius']", "Native radius owner changed")
    separation = packet.get("rebuilt_own_radius_separation", {})
    require(separation.get("native_alias_substitution_allowed") is False, "Rebuilt radius alias was allowed")
    template = packet.get("current_template_instantiation", {})
    require(template.get("can_instantiate_uninstantiated_template_that_asks_rebuiltOwnRadius_now") is False, "Rebuilt radius template was instantiated")
    return {
        "schema": packet.get("schema"),
        "radius_definition": {
            "node_id": NATIVE_RADIUS_NODE,
            "frame": clone(frame),
            "value_materialized": False,
            "source_value_present": False,
            "rebuilt_alias_allowed": False,
        },
        "P2_routes": clone(packet.get("native_p2_selectors", {})),
        "rebuilt_radius_separation": clone(separation),
        "template_instantiation": clone(template),
        "raw_packet": clone(dict(packet)),
    }


def validate_static_trace(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one of the two concrete static packets without evaluation."""

    require(isinstance(packet, Mapping), "Static trace packet must be an object")
    schema = packet.get("schema")
    if schema in {TRACE_SCHEMA, TRACE_SCHEMA_V1}:
        return validate_native_history_trace(packet)
    if schema == RADIUS_TRACE_SCHEMA:
        return validate_radius_trace(packet)
    raise AssertionError(f"Unexpected static trace schema: {schema}")


def load_frozen_static_packets(
    history_path: str | Path = HISTORY_TRACE,
    radius_path: str | Path = RADIUS_TRACE,
) -> dict[str, Any]:
    """Read and validate the two frozen static packets without executing code."""

    history_raw, history = read_json(history_path, HISTORY_TRACE_SHA256)
    radius_raw, radius = read_json(radius_path, RADIUS_TRACE_SHA256)
    history_view = validate_native_history_trace(history)
    radius_view = validate_radius_trace(radius)
    return {
        "history": history_view,
        "radius_frame": radius_view,
        "artifacts": {
            "history": {
                "path": wire_path(history_path),
                "sha256": sha_bytes(history_raw),
                "bytes": len(history_raw),
                "raw_bytes_exact": True,
                "raw_bytes_base64": base64.b64encode(history_raw).decode("ascii"),
            },
            "radius_frame": {
                "path": wire_path(radius_path),
                "sha256": sha_bytes(radius_raw),
                "bytes": len(radius_raw),
                "raw_bytes_exact": True,
                "raw_bytes_base64": base64.b64encode(radius_raw).decode("ascii"),
            },
        },
        "packet_count": 2,
        "source_values_evaluated": False,
        "provider_executed": False,
    }


def validate_saved_parameter_handles(handles: Mapping[str, Any]) -> dict[str, Any]:
    """Keep native and rebuilt handles separate and data-only."""

    require(isinstance(handles, Mapping), "Saved parameter handles must be an object")
    normalized: dict[str, dict[str, Any]] = {}
    for side in PARAMETER_SIDE_KEYS:
        side_handles = handles.get(side)
        require(isinstance(side_handles, Mapping), f"Saved parameter side missing: {side}")
        normalized[side] = {}
        for quantity in PARAMETER_KEYS:
            record = side_handles.get(quantity)
            require(isinstance(record, Mapping), f"Saved parameter handle missing: {side}.{quantity}")
            require(record.get("data_only") is True or record.get("source_value_present") is False, f"Parameter handle is not data-only: {side}.{quantity}")
            require("value" not in record, f"Parameter handle contains a source value: {side}.{quantity}")
            normalized[side][quantity] = clone(dict(record))
    return {
        "native": normalized["native"],
        "rebuilt": normalized["rebuilt"],
        "same_source_object_claimed": False,
        "native_rebuilt_equality_claimed": False,
        "separate_origin_records": True,
    }


def build_definition_request(
    history_trace: Mapping[str, Any],
    saved_parameter_handles: Mapping[str, Any],
    radius_trace: Mapping[str, Any] | None = None,
    *,
    history_v1_trace: Mapping[str, Any] | None = None,
    history_v1_artifact: Mapping[str, Any] | None = None,
    history_v2_artifact: Mapping[str, Any] | None = None,
    parameter_result_literal: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the root-owned provider request without invoking it."""

    normalized_history = validate_static_trace(history_trace)
    require(normalized_history.get("schema") == TRACE_SCHEMA, "History trace packet is required")
    require(isinstance(history_v1_trace, Mapping), "The unchanged v1 history packet is required")
    require(history_v1_trace.get("packet_type") == "native_history_bundle_source_trace", "History v1 packet type changed")
    require(history_v1_trace.get("schema") == TRACE_SCHEMA_V1, "History v1 packet schema changed")
    normalized_history_v1 = validate_static_trace(history_v1_trace)
    require(radius_trace is not None, "Radius trace packet is required")
    normalized_radius = validate_static_trace(radius_trace)
    require(normalized_radius.get("schema") == RADIUS_TRACE_SCHEMA, "Radius trace packet is required")
    normalized_handles = validate_saved_parameter_handles(saved_parameter_handles)
    v1_artifact = clone(dict(history_v1_artifact or {
        "path": wire_path(HISTORY_TRACE_V1),
        "sha256": HISTORY_TRACE_V1_SHA256,
        "raw_bytes_exact": True,
    }))
    v2_artifact = clone(dict(history_v2_artifact or {
        "path": wire_path(HISTORY_TRACE),
        "sha256": HISTORY_TRACE_SHA256,
        "raw_bytes_exact": True,
    }))
    require(v1_artifact.get("sha256") == HISTORY_TRACE_V1_SHA256 and v2_artifact.get("sha256") == HISTORY_TRACE_SHA256, "History raw pins changed")
    def auth_pin(entry: Mapping[str, Any]) -> dict[str, Any]:
        return {
            key: clone(entry[key])
            for key in ("path", "sha256", "bytes", "schema", "packet_type", "raw_bytes_exact")
            if key in entry
        }
    if parameter_result_literal is None:
        _, parameter_result_literal = read_json(PARENT_RESULT, SAVED_PARAMETER_RESULT_SHA256)
    require(isinstance(parameter_result_literal, Mapping), "Saved parameter result literal is required")
    history_wrapper = {
        "literal": clone(dict(history_trace)),
        "base_literal": clone(dict(history_v1_trace)),
        "raw_sha256": HISTORY_TRACE_SHA256,
    }
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "contract_status": ROOT_API_STATUS,
        "candidate_signature": ROOT_API_CANDIDATE_SIGNATURE,
        "provider_api_sha256": EXPECTED_KERNEL_API_SHA256,
        "metadata_only": True,
        "source_only": True,
        "history_trace": history_wrapper,
        "history_trace_base": {
            "raw_sha256": HISTORY_TRACE_V1_SHA256,
            "packet_type": "native_history_bundle_source_trace",
            "schema": TRACE_SCHEMA_V1,
            "literal_embedded_in_history_trace": True,
        },
        "radius_trace": {
            "literal": clone(dict(radius_trace)),
            "raw_sha256": RADIUS_TRACE_SHA256,
        },
        "parameter_result": {
            "literal": clone(dict(parameter_result_literal)),
            "raw_sha256": SAVED_PARAMETER_RESULT_SHA256,
        },
        "static_trace_packets": {
            "history": {
                "wrapper_ref": "history_trace",
                "schema": TRACE_SCHEMA,
                "v1_schema": TRACE_SCHEMA_V1,
                "v2_schema": TRACE_SCHEMA,
                "raw_sha256": HISTORY_TRACE_SHA256,
                "base_raw_sha256": HISTORY_TRACE_V1_SHA256,
                "raw_bytes_pins_separate": True,
            },
            "radius_frame": {"wrapper_ref": "radius_trace", "schema": RADIUS_TRACE_SCHEMA, "raw_sha256": RADIUS_TRACE_SHA256},
        },
        "static_trace_packet_count": 2,
        "source_trace_artifacts": {
            "history_v1": auth_pin(v1_artifact),
            "history_v2": auth_pin(v2_artifact),
            "radius_frame_original": {
                "path": wire_path(RADIUS_TRACE_V1),
                "sha256": RADIUS_TRACE_V1_SHA256,
                "schema": RADIUS_TRACE_SCHEMA,
                "raw_bytes_exact": True,
                "preserved_original_before_corrected_trace": True,
            },
            "radius_frame": {"path": wire_path(RADIUS_TRACE), "sha256": RADIUS_TRACE_SHA256, "raw_bytes_exact": True},
        },
        "saved_parameter_handles": normalized_handles,
        "saved_parameter_result_pointer": {
            "path": wire_path(PARENT_RESULT),
            "sha256": SAVED_PARAMETER_RESULT_SHA256,
            "data_only": True,
            "source_values_consumed": False,
        },
        "definition_slots": list(DEFINITION_SLOT_KEYS),
        "served_definition_quantity_order": list(DEFINITION_KEYS),
        "served_orders": list(ORDERS),
        "fields_served": EXPECTED_DEFINITION_RECORD_COUNT,
        "native_radius_definition_records": EXPECTED_RADIUS_RECORD_COUNT,
        "native_branch_bundles": EXPECTED_BRANCH_BUNDLE_COUNT,
        "required_type_context_rejections": EXPECTED_REJECTION_COUNT,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_flags": clone(EXPECTED_FLAGS),
        "history_channels": list(HISTORY_KEYS),
        "branch_order": list(BRANCH_KEYS),
        "parameter_order": list(PARAMETER_KEYS),
        "native_radius_node": NATIVE_RADIUS_NODE,
        "native_radius_frame": {
            "coordinate_node": NATIVE_COORDINATE_NODE,
            "domain": list(NATIVE_DOMAIN),
            "offset_literal": NATIVE_OFFSET_LITERAL,
            "physical_measure": "1",
        },
        "P2_branch_routes": {
            "background_generic": list(P2_BACKGROUND_ROOT_IDS),
            "modified_generic": list(P2_MODIFIED_ROOT_IDS),
            "must_not_blend": True,
        },
        "provider_policy": {
            "definition_provider_entries": 1,
            "source_value_provider_calls": 0,
            "count_contract_final": True,
            "required_branch_getter": "get_native_branch(bundle, branch)",
            "required_branch_set": list(BRANCH_KEYS),
            "mixed_branch_bundle_rejection": True,
            "rebuilt_radius_frame_rejection": True,
            "duplicate_factorial_request_rejection": True,
            "actual_rejection_counts": {
                "mixed_branch_rejections": 1,
                "rebuilt_frame_rejections": 1,
                "duplicate_normalization_rejections": 1,
            },
            "kernel_call_only_once": True,
            "no_fix_or_replay": True,
            "no_math_or_numeric_calls": True,
            "no_source_values": True,
            "no_ancestor_builder_query": True,
        },
        "claims": {
            "native_rebuilt_Pstar_equality": False,
            "native_rebuilt_delta_equality": False,
            "native_rebuilt_radius_equality": False,
            "history_value_computation": False,
            "radius_value_computation": False,
            "closure": False,
            "physical_completion": False,
            "rebuilt_template_instantiation": False,
        },
    }


def get_native_branch(bundle: Mapping[str, Any], branch: str) -> dict[str, Any]:
    """Return exactly one native branch; mixed bundles are rejected."""

    require(isinstance(branch, str), "Native branch getter requires a string branch")
    require(branch in BRANCH_KEYS, f"Unknown native branch: {branch}")
    selected = bundle.get("selected_branch")
    require(selected not in {"mixed", "background+modified", "both"}, "Mixed native branch bundle is rejected")
    definitions = bundle.get("definition_records")
    require(isinstance(definitions, Mapping), "Native definition records missing")
    selected_records = definitions.get(branch)
    require(isinstance(selected_records, Mapping), f"Native branch missing: {branch}")
    return clone(dict(selected_records))


def reject_rebuilt_radius_frame(radius_record: Mapping[str, Any]) -> None:
    """Reject a rebuilt radius frame from the native branch bundle."""

    require(isinstance(radius_record, Mapping), "Radius record must be an object")
    side = radius_record.get("side")
    owner = str(radius_record.get("radius_owner", ""))
    require(side != "rebuilt" and "Rm_factor" not in owner, "Rebuilt radius frame is rejected")
    require(radius_record.get("rebuilt_alias_allowed") is not True, "Native radius alias substitution is rejected")


def reject_duplicate_factorial_request(record: Mapping[str, Any]) -> None:
    """Reject a second factorial on an already projected ordinary P0 row."""

    require(isinstance(record, Mapping), "Factorial request must be an object")
    if record.get("quantity") == "P0":
        require(record.get("ordinary_Z_factorial_applied_once") is not True or record.get("request_factorial") is not True, "Duplicate P0 factorial request is rejected")


def validate_branch_rejection_contract(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Describe the three required type/context rejections without invoking code."""

    radius = bundle.get("radius_definition", {})
    reject_rebuilt_radius_frame(radius)
    p0 = bundle.get("definition_records", {}).get("background", {}).get("P0", [{}])[4]
    reject_duplicate_factorial_request(p0)
    # A mixed bundle is invalid by construction; this marker is deliberately
    # checked without constructing a second branch or provider object.
    require(bundle.get("mixed_branch_bundle_rejected") is not False, "Mixed branch rejection guard missing")
    return {
        "rejection_keys": list(REJECTION_KEYS),
        "required_branch_getter": "get_native_branch(bundle, branch)",
        "native_p2_routes": {
            "background_generic": list(P2_BACKGROUND_ROOT_IDS),
            "modified_generic": list(P2_MODIFIED_ROOT_IDS),
            "must_not_blend": True,
        },
        "type_context_rejections": EXPECTED_REJECTION_COUNT,
        "mixed_branch_bundle": True,
        "rebuilt_radius_frame": True,
        "duplicate_factorial_request": True,
        "provider_calls": 0,
    }


def build_checkpoint_plan() -> dict[str, Any]:
    """Describe the four checkpoint slots without creating any files."""

    root = REPO / "tmp" / COMPONENT
    paths = {
        "source": root / "source.json",
        "entered": root / "entered.json",
        "result": root / "result.json",
        "diagnostics": root / "diagnostics.jsonl",
    }
    return {
        "input_checkpoint_keys": list(INPUT_CHECKPOINT_KEYS),
        "paths": {key: wire_path(path) for key, path in paths.items()},
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "planned_checkpoint_count": PLANNED_CHECKPOINT_COUNT,
        "planned_count_confirmed": False,
        "new_failure_count_if_failure": PLANNED_FAILED_COUNT_IF_NEW_FAILURE,
        "new_failure_count_confirmed": False,
        "write_order": ["source", "entered", "result", "diagnostics"],
        "result_saved_before_validation": True,
        "entered_fsync_before_consumer": True,
    }


def build_ledger_preservation(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Carry the complete parent ledger and all four retained maps forward."""

    require(parent.get("commit") == PARENT_COMMIT, "Parent ledger commit changed")
    retained = parent.get("retained_maps")
    require(isinstance(retained, Mapping) and set(retained) == {
        "preserved_error_cap_metadata",
        "retained_construction_data",
        "retained_integral_term_inventory",
        "retained_cap_provenance",
    }, "All four retained maps are required")
    return {
        "parent_report": clone(parent["report"]),
        "parent_result": clone(parent["result"]),
        "parent_commit": PARENT_COMMIT,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "preserved_failed_artifact_count": PARENT_FAILED_ARTIFACT_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_FILE_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_PATH_COUNT,
        "parent_math_failure_preserved": True,
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "parent_false_flags": clone(parent["parent_false_flags"]),
        "parent_false_flag_values": clone(parent["parent_false_flag_values"]),
        "retained_maps": clone(dict(retained)),
        "caps_retained_without_reduction": True,
        "failed_artifacts_retained_without_promotion": True,
        "parent_result_is_not_a_replacement_for_full_report": True,
    }


def load_current_review(review_path: str | Path, review_sha256: str) -> dict[str, Any]:
    """Authenticate the review supplied to publish/seal at the current call."""

    require(review_path and review_sha256, "Current saved review path and SHA are required")
    raw, review = read_json(review_path, review_sha256)
    require(review.get("passed") is True, "Current saved review did not pass")
    require(review.get("independent_review") is True, "Current saved review is not independent")
    require(review.get("source_binding_accepted") is True, "Current saved review did not accept source binding")
    require(review.get("mathematical_proof_accepted") is False, "Current saved review accepted a mathematical proof")
    return {
        "path": wire_path(review_path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "schema": review.get("schema"),
        "review_kind": review.get("review_kind"),
        "independent_review": True,
        "source_binding_accepted": True,
        "mathematical_proof_accepted": False,
        "outside_science_pins": True,
        "raw_bytes_base64": base64.b64encode(raw).decode("ascii"),
    }


def load_saved_parameter_handles(path: str | Path, expected_sha256: str) -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    """Load only data-only parameter owner records for the provider input."""

    raw, packet = read_json(path, expected_sha256)
    handles = packet.get("saved_parameter_handles", packet)
    normalized = validate_saved_parameter_handles(handles)
    artifact = {
        "path": wire_path(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "schema": packet.get("schema"),
        "data_only": True,
        "source_values_consumed": False,
        "raw_bytes_exact": True,
        "raw_bytes_base64": base64.b64encode(raw).decode("ascii"),
    }
    return raw, normalized, artifact


def merge_hash_maps(*maps: Mapping[str, str]) -> dict[str, str]:
    merged: dict[str, str] = {}
    for mapping in maps:
        for path, digest in mapping.items():
            require(valid_sha(digest), f"Malformed source pin: {path}")
            require(path not in merged or merged[path] == digest, f"Source pin collision: {path}")
            merged[path] = digest
    return merged


def checkpoint_artifact_entries() -> list[dict[str, Any]]:
    return [artifact_entry(path, embed=True, outside_science_pins=True) for path in CHECKPOINTS if path.exists()]


def read_diagnostic_events() -> list[dict[str, Any]]:
    if not DIAGNOSTIC_CHECKPOINT.exists():
        return []
    events: list[dict[str, Any]] = []
    for line in DIAGNOSTIC_CHECKPOINT.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "Diagnostic event must be an object")
            events.append(value)
    return events


def diagnostic_summary() -> dict[str, Any]:
    events = read_diagnostic_events()
    require(len(events) >= MANAGER_PREFIX_EVENT_COUNT, "Diagnostic manager prefix is incomplete")
    require(events[0].get("event") == "manager_entry_started", "Diagnostic entry prefix changed")
    require(events[1].get("event") in {"manager_kernel_imported", "manager_kernel_or_derive_failed"}, "Diagnostic import prefix changed")
    callbacks = events[MANAGER_PREFIX_EVENT_COUNT:]
    counts: dict[str, int] = {}
    for event in callbacks:
        name = event.get("event")
        require(name in CALLBACK_EVENT_COUNTS, f"Unexpected callback event: {name}")
        counts[name] = counts.get(name, 0) + 1
    return {
        "event_count": len(events),
        "manager_prefix_event_count": MANAGER_PREFIX_EVENT_COUNT,
        "callback_event_count": len(callbacks),
        "callback_event_counts": counts,
        "callback_event_raw_rows": clone(callbacks),
        "journal_event_count_expected": JOURNAL_EVENT_COUNT,
    }


def append_diagnostic(event: Mapping[str, Any]) -> bytes:
    row = dict(event)
    require(isinstance(row.get("event"), str), "Diagnostic event name is required")
    return append_jsonl(DIAGNOSTIC_CHECKPOINT, row)


def callback_prefix_from_journal() -> tuple[dict[str, int], list[dict[str, Any]]]:
    """Recover callback progress without replaying or resetting a failed entry."""

    events = read_diagnostic_events()
    if len(events) <= MANAGER_PREFIX_EVENT_COUNT:
        return {}, []
    rows = events[MANAGER_PREFIX_EVENT_COUNT:]
    counts: dict[str, int] = {}
    for row in rows:
        event = row.get("event")
        if event in CALLBACK_EVENT_COUNTS:
            counts[event] = counts.get(event, 0) + 1
    return counts, clone(rows)


def build_report_skeleton(
    *,
    parent: Mapping[str, Any],
    current_brief: Mapping[str, Any],
    input_spec: Mapping[str, Any],
    trace_artifact: Mapping[str, Any],
    parameter_handles_artifact: Mapping[str, Any],
    code_artifacts: Mapping[str, Any],
    document_artifact: Mapping[str, Any] | None,
    handoff_artifacts: Mapping[str, Any],
    current_review: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Create an unaccepted metadata report; no source consumer is entered."""

    require(set(REPORT_SCOPE_ATTRIBUTES) == {
        "static_trace_packet",
        "saved_parameter_handles",
        "definition_records",
        "history_channel_records",
        "parameter_separation",
        "provider_contract",
        "provisional_counts",
        "ledger_preservation",
        "current_review_and_handoffs",
    }, "Report scope changed")
    require(set(code_artifacts) == set(CODE2_ARTIFACT_KEYS), "Exactly two code artifacts are required")
    require(set(handoff_artifacts) == set(HANDOFF_KEYS), "Exactly four handoff artifacts are required")

    return {
        "schema": REPORT_SCHEMA,
        "component": COMPONENT,
        "attempt_status": ROOT_API_STATUS,
        "passed": None,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "parent": clone(dict(parent)),
        "current_brief": {
            "path": current_brief.get("path"),
            "sha256": current_brief.get("sha256"),
            "bytes": current_brief.get("bytes"),
            "execution_rules": clone(current_brief.get("execution_rules")),
            "execution_rules_preserved_exact": True,
        },
        "current_brief_preimage": build_brief_preimage_record(current_brief),
        "report_scope_attributes": list(REPORT_SCOPE_ATTRIBUTES),
        "input_spec": clone(dict(input_spec)),
        "static_trace_packet_count": 2,
        "served_definition_record_count": EXPECTED_DEFINITION_RECORD_COUNT,
        "native_radius_definition_count": EXPECTED_RADIUS_RECORD_COUNT,
        "native_branch_bundle_count": EXPECTED_BRANCH_BUNDLE_COUNT,
        "required_type_context_rejection_count": EXPECTED_REJECTION_COUNT,
        "trace_artifact": clone(dict(trace_artifact)),
        "saved_parameter_handles_artifact": clone(dict(parameter_handles_artifact)),
        "code_artifacts": clone(dict(code_artifacts)),
        "document_artifact": clone(document_artifact) if document_artifact is not None else None,
        "handoff_artifacts": clone(dict(handoff_artifacts)),
        "input_checkpoint_plan": build_checkpoint_plan(),
        "ledger_preservation": build_ledger_preservation(parent),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "count_contract_final": True,
        "current_review": clone(dict(current_review)) if current_review is not None else None,
        "current_review_required_at_publish_and_seal": True,
        "new_flags": None,
        "provider_contract": {
            "status": ROOT_API_STATUS,
            "candidate_signature": ROOT_API_CANDIDATE_SIGNATURE,
            "definition_provider_call_count_finalized": False,
            "math_or_numeric_calls_finalized": False,
        },
        "raw_publication_lossless": True,
        "do_not_reenter_parent_or_failed_kernels": True,
    }


def build_brief_update(current_brief: Mapping[str, Any], *, report_pointer: Mapping[str, Any]) -> dict[str, Any]:
    """Return a proposed update while preserving execution_rules byte-for-byte."""

    rules = current_brief.get("execution_rules")
    require(isinstance(rules, list) and len(rules) == 7, "Execution rules must remain the restored seven-item list")
    return {
        "phase": ROOT_API_STATUS,
        "component": COMPONENT,
        "current_component": COMPONENT,
        "report": clone(dict(report_pointer)),
        "execution_rules": clone(rules),
        "execution_rules_preserved_exact": True,
        "metadata_only": True,
        "source_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_task": "freeze_final_v2_trace_and_source_guards_then_enter_once",
    }


def update_main_brief(patch: Mapping[str, Any], *, expected_sha256: str | None = None) -> tuple[bytes, dict[str, Any]]:
    """Write a compact main BRIEF update while preserving execution_rules."""

    raw, current = read_json(MAIN_BRIEF, expected_sha256)
    if expected_sha256 is not None:
        require(expected_sha256 == MAIN_BRIEF_SHA256 and len(raw) == 83915, "Main BRIEF preimage changed before first update")
    else:
        require(len(raw) < 100_000, "Current compact BRIEF exceeds the 100 KB limit")
    rules = current.get("execution_rules")
    require(isinstance(rules, list) and len(rules) == 7, "Main BRIEF execution_rules changed")
    safe_patch = brief_safe(portable(clone(dict(patch))))
    for key, value in safe_patch.items():
        if key == "execution_rules":
            require(value == rules, "Execution rules update is not byte-value preserving")
        else:
            current[key] = value
    current["execution_rules"] = clone(rules)
    current["execution_rules_preserved_exact"] = True
    compact = brief_safe(current)
    compact_raw = encode(compact)
    if len(compact_raw) >= 100_000:
        archive = WORK / f"RESUME_BRIEF_ARCHIVE_BEFORE_{COMPONENT}_2026_10_27.json"
        if not archive.exists():
            write_new(archive, raw)
        pointer = {"path": wire_path(archive), "sha256": sha_bytes(raw), "bytes": len(raw), "preimage_exact": True}
        compact["brief_preimage"] = pointer
        compact = brief_safe(compact)
        compact_raw = encode(compact)
    require(len(compact_raw) < 100_000, "Main BRIEF exceeds 100 KB")
    write_replace(MAIN_BRIEF, compact)
    return compact_raw, compact


def validate_provider_result(result: Mapping[str, Any], *, request: Mapping[str, Any]) -> dict[str, Any]:
    """Validate a future result shape without assuming its provider implementation."""

    require(isinstance(result, Mapping), "Provider result must be a plain object")
    require(result.get("schema") == RESULT_SCHEMA, "Provider result schema changed")
    require(result.get("component") == COMPONENT, "Provider result component changed")
    require(result.get("metadata_only") is True and result.get("source_only") is True, "Provider result markers changed")
    require(result.get("physical_layer_completed") is False, "Provider result claims physical completion")
    require(result.get("runtime_model_metadata") is None, "Runtime model metadata must be null")
    counts = result.get("actual_counts")
    require(isinstance(counts, Mapping), "Provider actual counts missing")
    require(dict(counts) == EXPECTED_COUNTS, "Provider actual counts changed")
    require(not any(key in result for key in ("source_values", "numeric_values", "evaluated_values")), "Provider result contains source values")
    flags = result.get("new_flags")
    require(isinstance(flags, Mapping) and dict(flags) == EXPECTED_FLAGS, "Provider new flags changed")
    claims = result.get("claims", {})
    require(isinstance(claims, Mapping), "Provider claims missing")
    for flag in NONCLAIM_FLAGS:
        require(claims.get(flag, False) is False, f"Provider result made forbidden claim: {flag}")
    return clone(dict(result))


def code2_artifacts(kernel_path: str | Path, runner_path: str | Path, *, kernel_sha256: str | None = None, runner_sha256: str | None = None) -> dict[str, dict[str, Any]]:
    kernel = Path(kernel_path).resolve()
    runner = Path(runner_path).resolve()
    require(kernel.is_relative_to(REPO.resolve()) and runner.is_relative_to(REPO.resolve()), "Code2 path escaped repository")
    require(kernel.exists() and runner.exists(), "Kernel and runner must exist before prepare")
    return {
        "kernel": artifact_entry(kernel, expected_sha256=kernel_sha256, embed=True, outside_science_pins=True),
        "runner": artifact_entry(runner, expected_sha256=runner_sha256, embed=True, outside_science_pins=True),
    }


def build_source_packet(
    *,
    parent: Mapping[str, Any],
    current_brief: Mapping[str, Any],
    request: Mapping[str, Any],
    input_raw: bytes,
    parameter_handles_artifact: Mapping[str, Any],
    code_artifacts: Mapping[str, Mapping[str, Any]],
    trace_artifacts: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    history_v2_artifact = trace_artifacts.get("history_v2", trace_artifacts.get("history"))
    require(isinstance(history_v2_artifact, Mapping), "v2 history artifact is required")
    additions: dict[str, str] = {
        entry["path"]: entry["sha256"]
        for entry in trace_artifacts.values()
        if isinstance(entry, Mapping) and entry.get("path") and entry.get("sha256")
    }
    additions.update({
        parameter_handles_artifact["path"]: parameter_handles_artifact["sha256"],
        code_artifacts["kernel"]["path"]: code_artifacts["kernel"]["sha256"],
        code_artifacts["runner"]["path"]: code_artifacts["runner"]["sha256"],
    })
    source_union = merge_hash_maps(parent["source_pins"], additions)
    return {
        "schema": COMPONENT + "_source_v1",
        "component": COMPONENT,
        "parent_component": PARENT_COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "current_head_commit": CURRENT_HEAD_COMMIT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "parent_report": clone(parent["report"]),
        "parent_result": clone(parent["result"]),
        "parent_report_view": clone(parent["report_view"]),
        "parent_false_flags": clone(parent["parent_false_flags"]),
        "parent_false_flag_values": clone(parent["parent_false_flag_values"]),
        "retained_maps": clone(parent["retained_maps"]),
        "parent_source_pins": clone(parent["source_pins"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "current_brief": clone({key: current_brief[key] for key in ("path", "sha256", "bytes", "execution_rules")}),
        "current_brief_preimage": build_brief_preimage_record(current_brief),
        "input_spec": {
            "schema": INPUT_SCHEMA,
            "sha256": sha_bytes(input_raw),
            "bytes": len(input_raw),
            "raw_bytes_exact": True,
            "raw_bytes_base64": base64.b64encode(input_raw).decode("ascii"),
        },
        "input": clone(dict(request)),
        "static_trace_artifacts": clone(dict(trace_artifacts)),
        "history_trace_wrapper": {
            "v1_path": trace_artifacts.get("history_v1", {}).get("path"),
            "v1_sha256": trace_artifacts.get("history_v1", {}).get("sha256"),
            "v2_path": history_v2_artifact.get("path"),
            "v2_sha256": history_v2_artifact.get("sha256"),
            "parsed_v1_and_v2_embedded_in_input": True,
            "raw_bytes_pins_separate": True,
        },
        "saved_parameter_handles_artifact": clone(dict(parameter_handles_artifact)),
        "code2": clone(dict(code_artifacts)),
        "provider_api_sha256": EXPECTED_KERNEL_API_SHA256,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "callback_event_counts": clone(CALLBACK_EVENT_COUNTS),
        "callback_event_count": CALLBACK_EVENT_COUNT,
        "manager_prefix_event_count": MANAGER_PREFIX_EVENT_COUNT,
        "journal_event_count": JOURNAL_EVENT_COUNT,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "no_runner_execution": True,
        "no_math_or_provider_execution_during_prepare": True,
        "do_not_reenter_parent_or_failed_kernels": True,
    }


def load_source() -> tuple[bytes, dict[str, Any]]:
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == COMPONENT + "_source_v1" and source.get("component") == COMPONENT, "Source checkpoint identity changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    require(source.get("current_head_commit") == CURRENT_HEAD_COMMIT, "Source current head changed")
    require(source.get("metadata_only") is True and source.get("source_only") is True, "Source markers changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    require(source.get("input_spec", {}).get("schema") == INPUT_SCHEMA, "Input schema changed")
    require(source.get("expected_counts") == EXPECTED_COUNTS, "Source count contract changed")
    require(source.get("expected_new_flags") == EXPECTED_FLAGS, "Source flag contract changed")
    require(source.get("provider_api_sha256") == EXPECTED_KERNEL_API_SHA256, "Provider API pin changed")
    require(PREPARE_RECEIPT.exists(), "Prepare receipt is missing")
    _, prepare_receipt = read_json(PREPARE_RECEIPT)
    require(prepare_receipt.get("source_checkpoint_sha256") == sha_bytes(raw), "Prepare receipt is not bound to source checkpoint")
    require(BRIEF.exists(), "Component pending BRIEF is missing")
    require(prepare_receipt.get("brief_sha256") == sha_path(BRIEF), "Pending component BRIEF changed")
    require(INPUT_SPEC.exists(), "Frozen input specification is missing")
    input_raw, input_value = read_json(INPUT_SPEC, source["input_spec"]["sha256"])
    require(input_value == source.get("input"), "Frozen input payload changed")
    require(source["input_spec"].get("bytes") == len(input_raw), "Frozen input byte count changed")
    encoded_input = source["input_spec"].get("raw_bytes_base64")
    require(isinstance(encoded_input, str) and base64.b64decode(encoded_input.encode("ascii")) == input_raw, "Frozen input raw bytes changed")
    _, history_v1_literal = read_json(HISTORY_TRACE_V1, HISTORY_TRACE_V1_SHA256)
    _, history_v2_literal = read_json(HISTORY_TRACE, HISTORY_TRACE_SHA256)
    _, radius_literal = read_json(RADIUS_TRACE, RADIUS_TRACE_SHA256)
    _, parameter_result_literal = read_json(PARENT_RESULT, SAVED_PARAMETER_RESULT_SHA256)
    history_wrapper = input_value.get("history_trace", {})
    require(history_wrapper.get("literal") == history_v2_literal, "Frozen v2 history literal changed")
    require(history_wrapper.get("base_literal") == history_v1_literal, "Frozen v1 history literal changed")
    require(history_wrapper.get("raw_sha256") == HISTORY_TRACE_SHA256, "Frozen v2 history SHA changed")
    require(input_value.get("history_trace_base", {}).get("raw_sha256") == HISTORY_TRACE_V1_SHA256, "Frozen v1 history SHA changed")
    require(input_value.get("radius_trace", {}).get("literal") == radius_literal, "Frozen radius literal changed")
    require(input_value.get("radius_trace", {}).get("raw_sha256") == RADIUS_TRACE_SHA256, "Frozen radius SHA changed")
    require(input_value.get("parameter_result", {}).get("literal") == parameter_result_literal, "Frozen parameter result literal changed")
    require(input_value.get("parameter_result", {}).get("raw_sha256") == SAVED_PARAMETER_RESULT_SHA256, "Frozen parameter result SHA changed")
    current_brief = source.get("current_brief", {})
    require(current_brief.get("path") == wire_path(MAIN_BRIEF), "Source BRIEF path changed")
    require(current_brief.get("sha256") == MAIN_BRIEF_SHA256, "Source BRIEF preimage SHA changed")
    code2 = source.get("code2", {})
    require(set(code2) == {"kernel", "runner"}, "Source code2 artifact set changed")
    for key in ("kernel", "runner"):
        entry = code2[key]
        code_path = resolve_path(entry["path"], REPO)
        require(wire_path(code_path) == entry["path"], f"Source {key} path changed")
        require(sha_path(code_path) == entry["sha256"], f"Source {key} bytes changed")
    require(len(DIAGNOSTIC_CHECKPOINT.read_bytes()) == 0, "Diagnostic journal is not empty before entry")
    return raw, source


def failure_result(error: BaseException, *, kernel_import_calls: int, kernel_imported: bool, derive_entry_calls: int, derive_calls: int) -> dict[str, Any]:
    callback_counts, callback_rows = callback_prefix_from_journal()
    partial_counts = {key: 0 for key in COUNT_KEYS}
    partial_counts["definition_provider_entries"] = callback_counts.get("definition_provider_entry", 0)
    partial_counts["field_definition_queries"] = callback_counts.get("field_definition_row", 0)
    partial_counts["radius_definition_queries"] = callback_counts.get("radius_definition_row", 0)
    partial_counts["native_bundle_manifests"] = callback_counts.get("native_bundle_manifest", 0)
    partial_counts["native_branches_served"] = partial_counts["native_bundle_manifests"]
    rejection_rows = [row for row in callback_rows if row.get("event") == "context_rejection"]
    rejection_keys = {str(row.get("rejection_key", row.get("kind", ""))) for row in rejection_rows}
    partial_counts["mixed_branch_rejections"] = int("mixed_branch_bundle" in rejection_keys or "mixed_branch_rejection" in rejection_keys)
    partial_counts["rebuilt_frame_rejections"] = int("rebuilt_radius_frame" in rejection_keys or "rebuilt_frame_rejection" in rejection_keys)
    partial_counts["duplicate_normalization_rejections"] = int("duplicate_factorial_request" in rejection_keys or "duplicate_normalization_rejection" in rejection_keys)
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": partial_counts,
        "new_flags": {key: False for key in FLAG_KEYS},
        "claims": {key: False for key in NONCLAIM_FLAGS},
        "failure": {"error_type": type(error).__name__, "message": str(error)},
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "runner_executed": False,
        "runtime_model_metadata": None,
        "callback_event_counts": callback_counts,
        "callback_event_raw_rows": callback_rows,
        "partial_result_from_durable_journal": True,
    }


def validate_result(
    result: Mapping[str, Any],
    source: Mapping[str, Any],
    *,
    derive_entry_calls: int,
    derive_calls: int,
    diagnostic: Mapping[str, Any],
) -> dict[str, Any]:
    violations: list[str] = []
    passed = result.get("passed") is True
    if passed:
        try:
            validate_provider_result(result, request=source["input"])
        except AssertionError as error:
            violations.append(str(error))
    if derive_entry_calls != 1:
        violations.append("derive entry call count changed")
    if derive_calls != 1:
        violations.append("derive call count changed")
    if diagnostic.get("event_count") != JOURNAL_EVENT_COUNT:
        violations.append("diagnostic journal count changed")
    if diagnostic.get("callback_event_counts") != CALLBACK_EVENT_COUNTS:
        violations.append("diagnostic callback counts changed")
    if passed and result.get("runner_executed") is not False:
        violations.append("runner execution marker changed")
    return {
        "passed": not violations and passed,
        "protocol_violations": violations,
        "actual_counts": clone(result.get("actual_counts")),
        "new_flags": clone(result.get("new_flags")),
        "metadata_only": result.get("metadata_only") is True,
        "source_only": result.get("source_only") is True,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "diagnostic": clone(dict(diagnostic)),
        "raw_result_saved_before_validation": True,
    }


def import_provider(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("_native_history_radius_bundle_kernel", path)
    require(spec is not None and spec.loader is not None, "Provider module specification failed")
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    derive = getattr(module, "derive", None)
    require(callable(derive), "Provider derive(request, callback) is missing")
    return module


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    require(all(not path.exists() for path in CHECKPOINTS), "Component checkpoints already exist; do not reenter")
    require(not INPUT_SPEC.exists() and not BRIEF.exists() and not PREPARE_RECEIPT.exists(), "Component preparation artifacts already exist")
    current_brief = load_current_brief()
    parent = load_parent_authority()
    _, parameter_result_literal = read_json(PARENT_RESULT, SAVED_PARAMETER_RESULT_SHA256)
    history_v1_raw, history_v1_packet = read_json(HISTORY_TRACE_V1, HISTORY_TRACE_V1_SHA256)
    history_raw, history_packet = read_json(HISTORY_TRACE, HISTORY_TRACE_SHA256)
    radius_v1_raw, radius_v1_packet = read_json(RADIUS_TRACE_V1, RADIUS_TRACE_V1_SHA256)
    radius_raw, radius_packet = read_json(RADIUS_TRACE, RADIUS_TRACE_SHA256)
    validate_native_history_trace(history_v1_packet)
    history_view = validate_native_history_trace(history_packet)
    validate_radius_trace(radius_v1_packet)
    radius_view = validate_radius_trace(radius_packet)
    require(history_view.get("definition_record_count") == EXPECTED_DEFINITION_RECORD_COUNT, "History definition count changed")
    require(radius_view.get("radius_definition", {}).get("node_id") == NATIVE_RADIUS_NODE, "Radius definition changed")
    parameter_path = resolve_path(args.parameter_handles, WORK)
    parameter_raw, parameter_handles, parameter_artifact = load_saved_parameter_handles(parameter_path, args.parameter_handles_sha256)
    history_v1_artifact = {
        "path": wire_path(HISTORY_TRACE_V1),
        "sha256": sha_bytes(history_v1_raw),
        "bytes": len(history_v1_raw),
        "schema": TRACE_SCHEMA_V1,
        "packet_type": "native_history_bundle_source_trace",
        "raw_bytes_exact": True,
        "raw_bytes_base64": base64.b64encode(history_v1_raw).decode("ascii"),
    }
    history_v2_artifact = {
        "path": wire_path(HISTORY_TRACE),
        "sha256": sha_bytes(history_raw),
        "bytes": len(history_raw),
        "schema": TRACE_SCHEMA,
        "raw_bytes_exact": True,
        "raw_bytes_base64": base64.b64encode(history_raw).decode("ascii"),
    }
    radius_v1_artifact = {
        "path": wire_path(RADIUS_TRACE_V1),
        "sha256": sha_bytes(radius_v1_raw),
        "bytes": len(radius_v1_raw),
        "schema": RADIUS_TRACE_SCHEMA,
        "raw_bytes_exact": True,
        "raw_bytes_base64": base64.b64encode(radius_v1_raw).decode("ascii"),
        "preserved_original_before_corrected_trace": True,
    }
    radius_v2_artifact = {
        "path": wire_path(RADIUS_TRACE),
        "sha256": sha_bytes(radius_raw),
        "bytes": len(radius_raw),
        "schema": RADIUS_TRACE_SCHEMA,
        "raw_bytes_exact": True,
        "raw_bytes_base64": base64.b64encode(radius_raw).decode("ascii"),
    }
    request = build_definition_request(
        history_packet,
        parameter_handles,
        radius_packet,
        history_v1_trace=history_v1_packet,
        history_v1_artifact=history_v1_artifact,
        history_v2_artifact=history_v2_artifact,
        parameter_result_literal=parameter_result_literal,
    )
    input_raw = encode(request)
    kernel_path = resolve_path(args.kernel_path or DEFAULT_KERNEL_PATH, REPO)
    runner_path = resolve_path(args.runner_path or DEFAULT_RUNNER_PATH, REPO)
    code2 = code2_artifacts(kernel_path, runner_path, kernel_sha256=args.kernel_sha256 or EXPECTED_KERNEL_SHA256, runner_sha256=args.runner_sha256)
    trace_artifacts = {
        "history_v1": history_v1_artifact,
        "history_v2": history_v2_artifact,
        "radius_frame_v1": radius_v1_artifact,
        "radius_frame": radius_v2_artifact,
    }
    source = build_source_packet(parent=parent, current_brief=current_brief, request=request, input_raw=input_raw, parameter_handles_artifact=parameter_artifact, code_artifacts=code2, trace_artifacts=trace_artifacts)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    input_written = write_new(INPUT_SPEC, request)
    diagnostic_raw = write_new(DIAGNOSTIC_CHECKPOINT, b"")
    component_brief = {
        "schema": COMPONENT + "_brief_v1",
        "component": COMPONENT,
        "phase": "prepared_unentered",
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_written)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "history_trace_v1": trace_artifacts["history_v1"],
        "history_trace_v2": trace_artifacts["history_v2"],
        "radius_trace": trace_artifacts["radius_frame"],
        "parameter_handles": parameter_artifact,
        "code2": code2,
        "metadata_only": True,
        "source_only": True,
        "no_runner_execution": True,
        "do_not_reenter": True,
    }
    brief_raw = write_new(BRIEF, component_brief)
    pending = {
        "schema": COMPONENT + "_pending_authority_v1",
        "component": COMPONENT,
        "phase": "prepared_unentered",
        "parent_commit": PARENT_COMMIT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_written),
        "actual_counts": None,
        "new_flags": None,
        "metadata_only": True,
        "source_only": True,
        "do_not_reenter": True,
    }
    main_raw, main_value = update_main_brief({
        "phase": "prepared_unentered",
        "current_component": COMPONENT,
        "component": COMPONENT,
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_written)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "pending_authority": pending,
        "expected_counts": EXPECTED_COUNTS,
        "expected_new_flags": EXPECTED_FLAGS,
        "execution_rules": clone(current_brief["execution_rules"]),
        "metadata_only": True,
        "source_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_task": "run_native_history_radius_bundle_once_after_entry",
    }, expected_sha256=MAIN_BRIEF_SHA256)
    receipt = {
        "schema": COMPONENT + "_prepare_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_written),
        "brief_sha256": sha_bytes(brief_raw),
        "diagnostic_checkpoint_sha256": sha_bytes(diagnostic_raw),
        "main_resume_brief_sha256": sha_bytes(main_raw),
        "pending_authority": pending,
        "source_union_count": source["source_union_count"],
        "parent_commit": PARENT_COMMIT,
        "runtime_model_metadata": None,
        "metadata_only": True,
        "source_only": True,
        "no_runner_execution": True,
    }
    receipt_raw = write_new(PREPARE_RECEIPT, receipt)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw), "main_resume_brief": main_value}


def run(args: argparse.Namespace) -> dict[str, Any]:
    source_raw, source = load_source()
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not RUN_RECEIPT.exists(), "Component was already entered or run")
    require(DIAGNOSTIC_CHECKPOINT.exists(), "Diagnostic checkpoint missing")
    source_kernel_path = resolve_path(source["code2"]["kernel"]["path"], REPO)
    source_runner_path = resolve_path(source["code2"]["runner"]["path"], REPO)
    kernel_path = resolve_path(args.kernel_path or source_kernel_path, REPO)
    runner_path = resolve_path(args.runner_path or source_runner_path, REPO)
    require(kernel_path == source_kernel_path, "Run kernel path differs from frozen code2 path")
    require(runner_path == source_runner_path, "Run runner path differs from frozen code2 path")
    require(sha_path(kernel_path) == source["code2"]["kernel"]["sha256"], "Run kernel bytes differ from frozen code2")
    require(sha_path(runner_path) == source["code2"]["runner"]["sha256"], "Run runner bytes differ from frozen code2")
    input_raw, input_value = read_json(INPUT_SPEC, source["input_spec"]["sha256"])
    require(input_value == source["input"], "Run input differs from frozen source payload")
    require(source["input_spec"].get("bytes") == len(input_raw), "Run input byte count differs from frozen source")
    require(len(DIAGNOSTIC_CHECKPOINT.read_bytes()) == 0, "Diagnostic journal must be empty before entry")
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "entered_at": datetime.now(timezone.utc).isoformat(),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "kernel_path": wire_path(kernel_path),
        "kernel_sha256": source["code2"]["kernel"]["sha256"],
        "runner_path": wire_path(runner_path),
        "runner_sha256": source["code2"]["runner"]["sha256"],
        "derive_contract": ROOT_API_CANDIDATE_SIGNATURE,
        "derive_entry_calls": 0,
        "derive_calls": 0,
        "kernel_import_calls": 1,
        "metadata_only": True,
        "source_only": True,
        "runner_executed": False,
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    append_diagnostic({"event": "manager_entry_started", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw)})
    derive_entry_calls = 0
    derive_calls = 0
    kernel_import_calls = 1
    kernel_imported = False
    try:
        module = import_provider(kernel_path)
        kernel_imported = True
        append_diagnostic({"event": "manager_kernel_imported", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "kernel_path": wire_path(kernel_path), "kernel_sha256": source["code2"]["kernel"]["sha256"]})

        def callback(event: dict[str, Any]) -> None:
            require(isinstance(event, dict), "Provider callback payload must be an object")
            require(event.get("event") in CALLBACK_EVENT_COUNTS, "Provider callback event name changed")
            append_diagnostic(event)

        # Count the sole planned entry before calling the provider so an
        # unexpected exception cannot be misreported as a replayable zero.
        derive_entry_calls = 1
        derive_calls = 1
        component_result = module.derive(source["input"], callback)
        require(isinstance(component_result, dict), "Provider derive result must be a plain object")
    except BaseException as error:
        if not kernel_imported:
            append_diagnostic({"event": "manager_kernel_or_derive_failed", "schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT, "error_type": type(error).__name__})
        component_result = failure_result(error, kernel_import_calls=kernel_import_calls, kernel_imported=kernel_imported, derive_entry_calls=derive_entry_calls, derive_calls=derive_calls)
    durable_result = clone(component_result)
    durable_result["component_result"] = clone(component_result)
    durable_result["component_result_exact"] = True
    durable_result["execution_metadata"] = {
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "runner_executed": False,
        "math_reexecuted": False,
        "old_math_reexecuted": False,
    }
    durable_result["kernel_import_calls"] = kernel_import_calls
    durable_result["kernel_imported"] = kernel_imported
    durable_result["derive_entry_calls"] = derive_entry_calls
    durable_result["derive_calls"] = derive_calls
    durable_result["runner_executed"] = False
    result_raw = write_new(RESULT_CHECKPOINT, durable_result)
    diagnostic_value = diagnostic_summary()
    validation = validate_result(durable_result, source, derive_entry_calls=derive_entry_calls, derive_calls=derive_calls, diagnostic=diagnostic_value)
    receipt = {
        "schema": COMPONENT + "_run_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "diagnostic_summary": diagnostic_value,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "validation": validation,
        "raw_result_saved_before_validation": True,
        "component_result_preserved": True,
        "rerun": False,
        "runner_executed": False,
        "runtime_model_metadata": None,
    }
    receipt_raw = write_new(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}


def _publication_artifact_set(source: Mapping[str, Any], review: Mapping[str, Any]) -> dict[str, Any]:
    current = checkpoint_artifact_entries()
    require(len(current) == 4, "Publication requires all four input checkpoints")
    helper = artifact_entry(Path(__file__).resolve(), embed=True, outside_science_pins=True)
    input_artifact = artifact_entry(INPUT_SPEC, expected_sha256=source["input_spec"]["sha256"], embed=True, outside_science_pins=True)
    _, prepare_receipt = read_json(PREPARE_RECEIPT)
    component_brief = artifact_entry(BRIEF, expected_sha256=prepare_receipt["brief_sha256"], embed=True, outside_science_pins=True)
    code2 = clone(source["code2"])
    trace_artifacts = clone(source["static_trace_artifacts"])
    parameter_artifact = clone(source["saved_parameter_handles_artifact"])
    return {
        "code2": code2,
        "input": input_artifact,
        "checkpoints": current,
        "metadata_helper": helper,
        "component_brief": component_brief,
        "main_brief_preimage": clone(source["current_brief_preimage"]),
        "current_review": clone(dict(review)),
        "traces": trace_artifacts,
        "parameter_handles": parameter_artifact,
    }


def _write_publication_documents(report_pointer: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    milestone_path = REPO / MILESTONE_RELATIVE
    milestone_text = (
        "# Current native history radius bundle (2026-10-27)\n\n"
        f"Report: `{report_pointer}`.\n\n"
        "This is a source-only metadata publication. It serves 90 native FunctionJet definition references and 1 native radius definition across 2 separated branches, with 3 type/context rejections and 98 provider callback rows. It preserves the native ordinary-Z P0 projection and keeps the rebuilt OwnRadius template uninstantiated. No field values, W/P2 values, mathematical identities, numeric bounds, or physical closure are claimed.\n\n"
        "Remaining task: source_native_radius_phase_pullback_for_branch_bundles. Preserve the incoming native Rm owner and radial phase handles while leaving unresolved cross-route identity and all physical caps unresolved.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(milestone_path, milestone_text)
    prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{report_pointer}`.\n"
        "Status: source-only metadata; do not reenter. Definitions: 90 field + 1 radius; branches: 2; rejections: 3; callbacks: 98. Remaining task: source_native_radius_phase_pullback_for_branch_bundles; no numerical closure.\n"
    ).encode("utf-8")
    preimages: list[dict[str, Any]] = []
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(old)[:12]}.md"
        if not preimage.exists():
            write_new(preimage, old)
        preimages.append({"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(old), "raw_bytes_exact": True})
        write_replace(path, prefix + old)
    return artifact_entry(milestone_path, embed=True, outside_science_pins=True), preimages


def publish(args: argparse.Namespace) -> dict[str, Any]:
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT)
    require(PUBLISH_RECEIPT.exists() is False, "Publication receipt already exists; do not republish")
    report_path = REPO / REPORT_RELATIVE
    require(not report_path.exists(), "Publication report already exists; do not republish")
    require(args.review and args.review_sha256, "Publication requires the current saved review path and SHA")
    review = load_current_review(resolve_path(args.review, WORK), args.review_sha256)
    diagnostic = diagnostic_summary()
    validation = validate_result(result, source, derive_entry_calls=result.get("derive_entry_calls", 0), derive_calls=result.get("derive_calls", 0), diagnostic=diagnostic)
    artifacts = _publication_artifact_set(source, review)
    milestone_artifact, handoff_preimages = _write_publication_documents(REPORT_RELATIVE.as_posix())
    current_hashes = {entry["path"]: entry["sha256"] for entry in artifacts["checkpoints"]}
    source_union = merge_hash_maps(source["source_union_hashes"], current_hashes)
    if validation["passed"]:
        failed_hashes = clone(source["failed_artifact_hashes"])
    else:
        failed_additions = {entry["path"]: entry["sha256"] for entry in artifacts["code2"].values()}
        failed_additions[artifacts["input"]["path"]] = artifacts["input"]["sha256"]
        failed_additions.update(current_hashes)
        failed_hashes = merge_hash_maps(source["failed_artifact_hashes"], failed_additions)
    report = {
        "schema": REPORT_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "passed": validation["passed"],
        "metadata_only": True,
        "publication_metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_report": clone(source["parent_report"]),
        "parent_result": clone(source["parent_result"]),
        "parent_report_view": clone(source["parent_report_view"]),
        "parent_false_flags": clone(source["parent_false_flags"]),
        "parent_false_flag_values": clone(source["parent_false_flag_values"]),
        "parent_math_failure_preserved": source["parent_math_failure_preserved"],
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "retained_maps": clone(source["retained_maps"]),
        "preserved_error_cap_metadata": clone(source["retained_maps"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(source["retained_maps"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(source["retained_maps"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(source["retained_maps"]["retained_cap_provenance"]),
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "checkpoint_count": PLANNED_CHECKPOINT_COUNT,
        "checkpoint_hashes": current_hashes,
        "source_pin_count": len(source_union),
        "source_union_count": len(source_union),
        "source_pins": source_union,
        "raw_sha_pins": source_union,
        "preserved_failed_artifact_hashes_union": failed_hashes,
        "preserved_failed_artifact_count": len(failed_hashes),
        "expected_failed_count_if_new_failure": PLANNED_FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_FILE_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_PATH_COUNT,
        "protected_file_hashes": clone(source["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(source["dirty_path_hashes"]),
        "input": clone(source["input"]),
        "input_spec": artifacts["input"],
        "component_brief": artifacts["component_brief"],
        "main_brief_preimage": artifacts["main_brief_preimage"],
        "static_trace_artifacts": artifacts["traces"],
        "saved_parameter_handles_artifact": artifacts["parameter_handles"],
        "code2": artifacts["code2"],
        "provider_api_sha256": EXPECTED_KERNEL_API_SHA256,
        "metadata_helper_artifact": artifacts["metadata_helper"],
        "current_review": artifacts["current_review"],
        "milestone_artifact": milestone_artifact,
        "handoff_preimages": handoff_preimages,
        "read_first": [wire_path(REPO / "AGENTS.md"), wire_path(REPO / MILESTONE_RELATIVE), wire_path(report_path)],
        "current_doc_sha256": milestone_artifact["sha256"],
        "current_document_sha256": milestone_artifact["sha256"],
        "gitattributes_required": {
            "paths": [source["code2"]["kernel"]["path"], source["code2"]["runner"]["path"], wire_path(REPO / REPORT_RELATIVE), wire_path(REPO / MILESTONE_RELATIVE)],
            "rule": "-text",
            "apply_before_freeze": True,
        },
        "full_result": clone(result),
        "full_result_raw_sha256": sha_bytes(result_raw),
        "result_checkpoint": {"path": wire_path(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw), "bytes": len(result_raw)},
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_path(ENTERED_CHECKPOINT),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "actual_counts": validation["actual_counts"],
        "new_flags": validation["new_flags"],
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "callback_event_counts": clone(CALLBACK_EVENT_COUNTS),
        "callback_event_count": CALLBACK_EVENT_COUNT,
        "manager_prefix_event_count": MANAGER_PREFIX_EVENT_COUNT,
        "journal_event_count": JOURNAL_EVENT_COUNT,
        "durable_result_validation": validation,
        "report_scope_attributes": list(REPORT_SCOPE_ATTRIBUTES),
        "lossless_artifacts": [*artifacts["code2"].values(), artifacts["input"], artifacts["component_brief"], artifacts["main_brief_preimage"], *artifacts["checkpoints"], artifacts["metadata_helper"], artifacts["current_review"], *artifacts["traces"].values(), artifacts["parameter_handles"]],
        "raw_bytes_lossless": True,
        "runtime_model_metadata": None,
        "do_not_reenter": True,
        "next_executable_task": "source_native_radius_phase_pullback_for_branch_bundles",
    }
    report_raw = write_new(report_path, report)
    brief_raw, brief = update_main_brief({
        "phase": "published_passed" if report["passed"] else "published_failed",
        "current_component": COMPONENT,
        "report": wire_path(report_path),
        "report_sha256": sha_bytes(report_raw),
        "latest_review": artifacts["current_review"],
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "read_first": [wire_path(REPO / "AGENTS.md"), wire_path(REPO / MILESTONE_RELATIVE), wire_path(report_path)],
        "current_doc_sha256": milestone_artifact["sha256"],
        "current_document_sha256": milestone_artifact["sha256"],
        "gitattributes_required": report["gitattributes_required"],
        "expected_counts": EXPECTED_COUNTS,
        "expected_new_flags": EXPECTED_FLAGS,
        "metadata_only": True,
        "source_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_task": report["next_executable_task"],
    })
    receipt = {
        "schema": COMPONENT + "_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(report_path),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "current_review": artifacts["current_review"],
        "main_resume_brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": brief,
        "runtime_model_metadata": None,
    }
    receipt_raw = write_new(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}


def resume(args: argparse.Namespace) -> dict[str, Any]:
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == REPORT_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    require(report.get("full_result") == result, "Published result changed")
    validation = report.get("durable_result_validation")
    require(isinstance(validation, dict), "Durable validation missing")
    require(validation == validate_result(result, source, derive_entry_calls=result.get("derive_entry_calls", 0), derive_calls=result.get("derive_calls", 0), diagnostic=diagnostic_summary()), "Durable validation changed")
    brief_raw, brief = update_main_brief({
        "phase": "resumed_published",
        "current_component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "latest_review": report.get("current_review"),
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "metadata_only": True,
        "source_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_task": report.get("next_executable_task"),
    })
    receipt = {
        "schema": COMPONENT + "_resume_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "main_resume_brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": brief,
        "do_not_reenter": True,
    }
    receipt_raw = write_replace(RESUME_RECEIPT, receipt)
    return {"resume": receipt, "resume_sha256": sha_bytes(receipt_raw)}


def seal(args: argparse.Namespace) -> dict[str, Any]:
    require(args.review and args.review_sha256, "Seal requires the current saved review path and SHA")
    require(not SEAL_RECEIPT.exists(), "Final seal already exists; do not reseal")
    current_review = load_current_review(resolve_path(args.review, WORK), args.review_sha256)
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == REPORT_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    require(report.get("current_review", {}).get("sha256") == current_review["sha256"], "Seal review SHA differs from published current review")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    require(report.get("full_result") == result, "Published result changed before seal")
    validation = report.get("durable_result_validation")
    require(isinstance(validation, dict), "Durable validation missing before seal")
    require(validation == validate_result(result, source, derive_entry_calls=result.get("derive_entry_calls", 0), derive_calls=result.get("derive_calls", 0), diagnostic=diagnostic_summary()), "Durable validation changed before seal")
    authority = {
        "schema": COMPONENT + "_authoritative_state_v1",
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "passed": report.get("passed") is True,
        "actual_counts": clone(report.get("actual_counts")),
        "new_flags": clone(report.get("new_flags")),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": clone(EXPECTED_FLAGS),
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "current_review": current_review,
        "do_not_reenter": True,
    }
    brief_raw, brief = update_main_brief({
        "phase": "sealed",
        "current_component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "latest_review": current_review,
        "latest_QA": validation,
        "latest_qa": validation,
        "latest_final_qa": validation,
        "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "read_first": clone(report.get("read_first")),
        "current_doc_sha256": report.get("current_doc_sha256"),
        "current_document_sha256": report.get("current_document_sha256"),
        "gitattributes_required": clone(report.get("gitattributes_required")),
        "latest_authoritative_state": authority,
        "metadata_only": True,
        "source_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_task": report.get("next_executable_task"),
    })
    receipt = {
        "schema": SEAL_SCHEMA,
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "current_review": current_review,
        "authority": authority,
        "main_resume_brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": brief,
        "runtime_model_metadata": None,
        "do_not_reenter": True,
    }
    receipt_raw = write_replace(SEAL_RECEIPT, receipt)
    return {"seal": receipt, "seal_sha256": sha_bytes(receipt_raw)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="One-shot native history/radius metadata lifecycle manager")
    sub = parser.add_subparsers(dest="action", required=True)

    prepare_parser = sub.add_parser("prepare", help="freeze source, input, code2, and pending BRIEF before entry")
    prepare_parser.add_argument("--parameter-handles", required=True)
    prepare_parser.add_argument("--parameter-handles-sha256", required=True)
    prepare_parser.add_argument("--kernel-path")
    prepare_parser.add_argument("--runner-path")
    prepare_parser.add_argument("--kernel-sha256")
    prepare_parser.add_argument("--runner-sha256")

    run_parser = sub.add_parser("run", help="enter once, import once, derive once, then save before validation")
    run_parser.add_argument("--kernel-path")
    run_parser.add_argument("--runner-path")

    publish_parser = sub.add_parser("publish", help="publish the durable result with the current saved review")
    publish_parser.add_argument("--review", required=True)
    publish_parser.add_argument("--review-sha256", required=True)

    sub.add_parser("resume", help="resume and revalidate the published durable state")

    seal_parser = sub.add_parser("seal", help="seal the published state against the current saved review")
    seal_parser.add_argument("--review", required=True)
    seal_parser.add_argument("--review-sha256", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> dict[str, Any]:
    args = build_parser().parse_args(argv)
    handlers = {
        "prepare": prepare,
        "run": run,
        "publish": publish,
        "resume": resume,
        "seal": seal,
    }
    result = handlers[args.action](args)
    print(json.dumps(result, ensure_ascii=True, sort_keys=True, separators=(",", ":")))
    return result


def build_one_shot_entry_plan(request: Mapping[str, Any]) -> dict[str, Any]:
    """Return the root-only entry protocol; this function never invokes a provider."""

    require(request.get("contract_status") == ROOT_API_STATUS, "Root API contract is already finalized or malformed")
    return {
        "status": "ready_for_single_entry",
        "request_schema": request.get("schema"),
        "candidate_signature": ROOT_API_CANDIDATE_SIGNATURE,
        "entered_checkpoint_before_provider": True,
        "provider_call_count": 1,
        "expected_callback_event_count": CALLBACK_EVENT_COUNT,
        "kernel_call_only_once": True,
        "complete_result_save_before_validation": True,
        "no_replay_after_failure": True,
        "provider_invocation_required_from_root": True,
        "no_provider_replay_after_entry": True,
    }


def build_schema_inventory() -> dict[str, Any]:
    """Describe the draft contract for the parent agent and later root freeze."""

    return {
        "schema": COMPONENT + "_schema_inventory_v1",
        "component": COMPONENT,
        "status": ROOT_API_STATUS,
        "candidate_signature": ROOT_API_CANDIDATE_SIGNATURE,
        "frozen_static_packet_count": 2,
        "history_trace": {"path": wire_path(HISTORY_TRACE), "sha256": HISTORY_TRACE_SHA256, "schema": TRACE_SCHEMA},
        "radius_trace": {"path": wire_path(RADIUS_TRACE), "sha256": RADIUS_TRACE_SHA256, "schema": RADIUS_TRACE_SCHEMA},
        "saved_parameter_result": {"path": wire_path(PARENT_RESULT), "sha256": SAVED_PARAMETER_RESULT_SHA256, "data_only": True},
        "parent": {
            "component": PARENT_COMPONENT,
            "parent_commit": PARENT_COMMIT,
            "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
            "current_head_commit": CURRENT_HEAD_COMMIT,
            "report_sha256": PARENT_REPORT_SHA256,
            "result_sha256": PARENT_RESULT_SHA256,
            "main_brief_sha256": MAIN_BRIEF_SHA256,
            "source_pin_count": PARENT_SOURCE_PIN_COUNT,
            "checkpoint_count": PARENT_CHECKPOINT_COUNT,
            "failed_artifact_count": PARENT_FAILED_ARTIFACT_COUNT,
            "protected_file_count": PARENT_PROTECTED_FILE_COUNT,
            "dirty_path_count": PARENT_DIRTY_PATH_COUNT,
            "false_flag_count": PARENT_FALSE_FLAG_COUNT,
        },
        "definition_slots": list(DEFINITION_SLOT_KEYS),
        "history_channels": list(HISTORY_KEYS),
        "branches": list(BRANCH_KEYS),
        "parameter_separation": {
            "native": {key: "direct same-chain owner record" for key in PARAMETER_KEYS},
            "rebuilt": {key: "independent source-owner record" for key in PARAMETER_KEYS},
            "same_source_object_claimed": False,
        },
        "report_scope_attributes": list(REPORT_SCOPE_ATTRIBUTES),
        "code_artifact_count": len(CODE2_ARTIFACT_KEYS),
        "input_checkpoint_count": len(INPUT_CHECKPOINT_KEYS),
        "handoff_artifact_count": len(HANDOFF_KEYS),
        "provisional_counts": clone(PROVISIONAL_COUNTS),
        "definition_record_count": EXPECTED_DEFINITION_RECORD_COUNT,
        "radius_definition_count": EXPECTED_RADIUS_RECORD_COUNT,
        "branch_bundle_count": EXPECTED_BRANCH_BUNDLE_COUNT,
        "rejection_count": EXPECTED_REJECTION_COUNT,
        "rejection_keys": list(REJECTION_KEYS),
        "planned_checkpoint_count": PLANNED_CHECKPOINT_COUNT,
        "planned_failed_count_if_new_failure": PLANNED_FAILED_COUNT_IF_NEW_FAILURE,
        "counts_finalized": True,
        "source_guards_finalized": True,
        "requires_current_review_args_for_publish_and_seal": True,
        "lossless_code2_and_input4_checkpoints": True,
        "external_metadata_helper": wire_path(Path(__file__)),
        "importing_module_does_not_execute_provider": True,
        "lifecycle_cli_actions_are_explicit": True,
        "runtime_model_metadata": None,
    }


__all__ = [
    "COMPONENT",
    "INPUT_SCHEMA",
    "TRACE_SCHEMA",
    "RESULT_SCHEMA",
    "REPORT_SCHEMA",
    "SEAL_SCHEMA",
    "build_brief_update",
    "build_brief_preimage_record",
    "build_checkpoint_plan",
    "build_definition_request",
    "build_ledger_preservation",
    "build_one_shot_entry_plan",
    "build_parser",
    "build_report_skeleton",
    "build_schema_inventory",
    "extract_native_definition_bundle",
    "get_native_branch",
    "load_frozen_static_packets",
    "load_current_brief",
    "load_current_review",
    "load_parent_authority",
    "main",
    "prepare",
    "publish",
    "resume",
    "run",
    "seal",
    "reject_duplicate_factorial_request",
    "reject_rebuilt_radius_frame",
    "validate_provider_result",
    "validate_branch_rejection_contract",
    "validate_native_history_trace",
    "validate_radius_trace",
    "validate_saved_parameter_handles",
    "validate_static_trace",
]


if __name__ == "__main__":
    main()

