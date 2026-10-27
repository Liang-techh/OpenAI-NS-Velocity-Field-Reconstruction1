"""Inert metadata lifecycle manager for the native parameter-Z certificate.

This manager is deliberately external to the project tree.  Importing it does
not import a project module, provider, runner, kernel, symbolic package, or
math helper.  It reads JSON and parses Python with :mod:`ast` only.  The
``prepare`` action freezes authenticated source metadata; ``run`` is the
future one-shot entry point that imports only the frozen new kernel after the
durable entry marker and saves its result; ``publish`` consumes the saved
bytes and an independent review.  The run action is not invoked while this
manager is authored or reviewed.

The component's final result schema, literal-wrapper names, counts, callback
events, node visits, and target counts are read from the frozen API contract.
No mathematical count is guessed in this file.
"""

from __future__ import annotations

import argparse
import ast
import base64
from collections.abc import Mapping, Sequence
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO = Path(
    r"C:/Users/z5242/Documents/Codex/2026-09-13/"
    r"before-doing-any-project-work-verify/work/st073-transition-next"
)
WORK = Path(
    r"C:/Users/z5242/Documents/Codex/2026-10-16/"
    r"https-github-com-liang-techh-openai/work"
)

COMPONENT = "native_parameter_Z_independence"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
RESULT_SCHEMA = f"original_{COMPONENT}"
PUBLICATION_SCHEMA = f"{RESULT_SCHEMA}_publication_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
DIAGNOSTIC_SCHEMA = f"{COMPONENT}_diagnostic_v1"

PARENT_COMPONENT = "native_O3_cylindrical_pullback"
PARENT_COMMIT = "0bc7d7eb694467707a954b5a383bbc8dd5fc454d"
CURRENT_HEAD_COMMIT = PARENT_COMMIT
PARENT_REPORT = REPO / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_O3_cylindrical_pullback.json"
)
PARENT_REPORT_SHA256 = "8b4ef714ffad4c7181484d6b29511637fd44c3412bd6783f806b4f7e324de927"
PARENT_RESULT = REPO / "tmp/native_O3_cylindrical_pullback/result.json"
PARENT_RESULT_SHA256 = "6e9a3aa733a080d3203328b91f850ae7ad18d29049df2421a8af409b17249679"

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "aaae8132ab57599c0e766c67c26185aaae65c1c2c0c237b05979cb684a4b597a"
MAIN_BRIEF_BYTES = 80769

KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_parameter_Z_independence.py")
RUNNER_RELATIVE = Path("scripts/run_native_parameter_Z_independence.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_parameter_Z_independence.json"
)
DOC_RELATIVE = Path("docs/CURRENT_NATIVE_PARAMETER_Z_INDEPENDENCE_2026_10_27.md")
HANDOFF_RELATIVES = (
    Path("docs/AGENT_TASKS.md"),
    Path("docs/CROSS_DEVICE_HANDOFF.md"),
    Path("docs/CURRENT_CHECKPOINT.md"),
    Path("docs/NEXT_TASKS.md"),
)
SCOPED_RELATIVES = (
    Path(".gitattributes"),
    KERNEL_RELATIVE,
    RUNNER_RELATIVE,
    REPORT_RELATIVE,
    DOC_RELATIVE,
    *HANDOFF_RELATIVES,
)

CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
JOURNAL_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (
    SOURCE_CHECKPOINT,
    ENTERED_CHECKPOINT,
    RESULT_CHECKPOINT,
    JOURNAL_CHECKPOINT,
)

API_PATH = WORK / "native_parameter_Z_independence_api_2026_10_27.json"
INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"

# These are the two authenticated source packets named by the current brief.
# The API may repeat the pins; it may not silently replace them.
FACT_SPECS = {
    "logC_logP": (
        WORK / "native_logC_logP_scalar_source_literals_2026_10_27.json",
        "a533245e802dda5391601b26e5939da751db1fed9f3f35d5a328f9915b4bfd9d",
        "bounded_static_native_scalar_source_ast_literal_rows_v1",
    ),
    "Ra_phase_origin": (
        WORK / "native_Ra_phase_origin_source_literals_2026_10_27.json",
        "536f693ea89e25066cadda6191f92c57effb48489d74097cb6dc5baacb9beb00",
        "bounded_static_native_phase_parameter_source_ast_literal_rows_v1",
    ),
    "raw_capsule": (
        WORK / "native_parameter_Z_raw_capsule_2026_10_27.json",
        "f5ef9950a04fcf0d33fdfd99dd03b799baaae433c6bc45285a0d7020663fcf7a",
        "native_parameter_Z_raw_graph_records_v1",
    ),
}
EXPECTED_WRAPPER_NAMES = (
    "saved_metric_result",
    "scalar_literals",
    "phase_literals",
    "raw_capsule",
)

PARENT_SOURCE_PIN_COUNT = 2276
PARENT_CHECKPOINT_COUNT = 268
PARENT_FAILED_COUNT = 183
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
PUBLISHED_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_COUNT + 7

REPORT_SCOPE_ATTRIBUTES = (
    "parent_ledgers_and_false_flags",
    "retained_error_caps_and_construction_maps",
    "input_literal_and_static_facts",
    "kernel_runner_api_artifacts",
    "durable_entry_result_and_journal",
    "dynamic_counts_flags_and_certificate_scope",
    "independent_saved_byte_review",
    "parameter_Z_documentation",
    "publication_read_first_and_scope_integrity",
)

EXECUTION_RULES = [
    "Never rerun any entered ancestor mathematical/helper family, including failed leaves, old builders/providers/normalizers/checkers, source-value reducers, quadrature or recurrence. Static source/AST/literal JSON/hash inspection is allowed.",
    "Before a NEW mathematical/source-expression consumer: freeze code and input, persist source checkpoint and pending brief, persist entry before invocation, invoke once, save complete result before review. Preserve a failure and repair in a separate frozen leaf; never reenter.",
    "Authenticate raw gzip SHA BEFORE decoding. Never decode docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json or the huge Rm_loop_density_identity.json; raw SHA only.",
    "Never restart stopped session 17273. Retain nonzero incoming histories and ordinary Z1. Bounds, midpoints and interval overlap are not function values; remove old numerical errors only with a same-source packet theorem.",
    "Use current tree and GitHub; preserve all dirty and untracked files. Commit/push only exact accepted scoped files. Goal remains ACTIVE / INCOMPLETE.",
    "Children: luna_max_worker, gpt-5.6-luna/max fixed role, fork_turns:none, complete brief. Runtime model metadata unavailable means null; never claim it was verified.",
    "A worker function that checks mathematical source identities counts as a source consumer even when it consumes static receipts or is called a pure self-check. Workers return unexecuted code; root alone freezes code and inputs, records entry, invokes once and saves the complete result before review. AST parsing and metadata-only audits remain permitted. Do not execute or reconstruct the rejected memory draft result.",
]


def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def clone(value: Any) -> Any:
    return deepcopy(value)


def encode(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
    )


def wire(path: str | Path) -> str:
    return Path(path).as_posix()


def resolve_path(value: str | Path, base: Path = REPO) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    for candidate in (base / path, REPO / path, WORK / path):
        if candidate.exists():
            return candidate
    return base / path


def read_bytes(path: str | Path, expected: str | None = None) -> bytes:
    path = Path(path)
    require(path.exists(), f"Missing artifact: {wire(path)}")
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), f"Changed artifact: {wire(path)}")
    return raw


def read_json(path: str | Path, expected: str | None = None) -> tuple[bytes, dict[str, Any]]:
    raw = read_bytes(path, expected)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"JSON object required: {wire(path)}")
    return raw, value


def _write(path: str | Path, value: bytes | Any, mode: str) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open(mode) as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_new(path: str | Path, value: bytes | Any) -> bytes:
    return _write(path, value, "xb")


def write_replace(path: str | Path, value: bytes | Any) -> bytes:
    return _write(path, value, "wb")


def artifact_from_bytes(path: str | Path, raw: bytes, *, embed: bool = True) -> dict[str, Any]:
    result = {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
    }
    if embed:
        result["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
    return result


def artifact(path: str | Path, expected: str | None = None, *, embed: bool = True) -> dict[str, Any]:
    return artifact_from_bytes(path, read_bytes(path, expected), embed=embed)


def artifact_pointer(value: Mapping[str, Any]) -> dict[str, Any]:
    """Reference an already embedded artifact without duplicating raw bytes."""
    return {key: clone(item) for key, item in value.items() if key != "raw_bytes_base64"}


def merge_hashes(*maps: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for mapping in maps:
        require(isinstance(mapping, Mapping), "Hash map must be a mapping")
        for name, digest in mapping.items():
            require(valid_sha(digest), f"Malformed SHA for {name}")
            key = wire(name)
            digest = digest.lower()
            require(key not in result or result[key] == digest, f"Conflicting SHA for {key}")
            result[key] = digest
    return result


def brief_safe(value: Any) -> Any:
    """Remove only oversized raw bytes and callback rows from a brief."""
    if isinstance(value, dict):
        return {
            key: brief_safe(child)
            for key, child in value.items()
            if key not in {"raw_bytes_base64", "callback_event_raw_rows"}
        }
    if isinstance(value, list):
        return [brief_safe(child) for child in value]
    return value


def brief_preimage(raw: bytes) -> dict[str, Any]:
    path = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(raw)[:12]}.json"
    return {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "preserved_exact_bytes": True,
        "raw_bytes_base64": base64.b64encode(raw).decode("ascii"),
    }


def load_current_brief(*, exact: bool = False) -> tuple[bytes, dict[str, Any]]:
    raw, brief = read_json(MAIN_BRIEF, MAIN_BRIEF_SHA256 if exact else None)
    if exact:
        require(len(raw) == MAIN_BRIEF_BYTES, "Current main brief byte count changed")
    require(brief.get("execution_rules") == EXECUTION_RULES, "The seven execution rules changed")
    require(
        brief.get("goal_active") is True and brief.get("full_goal_complete") is False,
        "Goal is not ACTIVE / INCOMPLETE",
    )
    require(
        brief.get("metadata_only") is True and brief.get("source_only") is True,
        "Accepted metadata/source-only baseline changed",
    )
    return raw, brief


def update_main_brief(patch: Mapping[str, Any], *, expected_rules: Sequence[str] = EXECUTION_RULES) -> dict[str, Any]:
    raw, brief = read_json(MAIN_BRIEF)
    require(brief.get("execution_rules") == list(expected_rules), "Main brief execution_rules changed")
    preimage = brief_preimage(raw)
    preimage_path = Path(preimage["path"])
    if not preimage_path.exists():
        write_new(preimage_path, raw)
    for key, value in patch.items():
        if key == "in_progress" and isinstance(brief.get(key), dict) and isinstance(value, Mapping):
            merged = clone(brief[key])
            merged.update(clone(value))
            brief[key] = merged
        else:
            brief[key] = clone(value)
    brief["execution_rules"] = list(expected_rules)
    brief["execution_rules_preserved_exact"] = True
    new_raw = encode(brief_safe(brief))
    require(len(new_raw) < 100_000, "Updated main brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, new_raw)
    return {
        "path": wire(MAIN_BRIEF),
        "sha256": sha_bytes(new_raw),
        "bytes": len(new_raw),
        "preimage": {key: value for key, value in preimage.items() if key != "raw_bytes_base64"},
    }


def parent_pointer(path: Path, digest: str, raw: bytes) -> dict[str, Any]:
    """A hash-only pointer; the 41 MB parent report is never embedded."""
    return {
        "path": wire(path),
        "sha256": digest,
        "bytes": len(raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
        "literal_embedded": False,
        "hash_only": True,
    }


def parent_report_view(report: Mapping[str, Any]) -> dict[str, Any]:
    keys = (
        "schema",
        "component",
        "passed",
        "metadata_only",
        "publication_metadata_only",
        "source_only",
        "conditional_geometry_scope",
        "physical_layer_completed",
        "full_goal_complete",
        "actual_counts",
        "new_flags",
        "expected_counts",
        "expected_new_flags",
        "parent_false_flags",
        "parent_false_flag_values",
        "next_executable_task",
        "runtime_model_metadata",
    )
    return {key: clone(report[key]) for key in keys if key in report}


def load_parent() -> dict[str, Any]:
    report_raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    result_raw, result = read_json(PARENT_RESULT, PARENT_RESULT_SHA256)
    require(report.get("component") == PARENT_COMPONENT and report.get("passed") is True, "Parent report identity changed")
    require(
        report.get("metadata_only") is False
        and report.get("source_only") is True
        and report.get("full_goal_complete") is False,
        "Parent report markers changed",
    )
    require(result.get("component") == PARENT_COMPONENT and result.get("passed") is True, "Parent result identity changed")
    require(report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256, "Parent result pointer changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Parent checkpoint count changed")
    maps = (
        "source_pins",
        "protected_file_hashes",
        "legacy_dirty_path_hashes",
        "preserved_failed_artifact_hashes_union",
        "parent_false_flags",
        "parent_false_flag_values",
        "preserved_error_cap_metadata",
        "retained_construction_data",
        "retained_integral_term_inventory",
        "retained_cap_provenance",
    )
    for key in maps:
        require(key in report, f"Parent ledger missing: {key}")
    require(len(report["source_pins"]) == PARENT_SOURCE_PIN_COUNT, "Parent source pin count changed")
    require(len(report["protected_file_hashes"]) == PARENT_PROTECTED_COUNT, "Parent protected count changed")
    require(len(report["legacy_dirty_path_hashes"]) == PARENT_DIRTY_COUNT, "Parent dirty count changed")
    require(
        len(report["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_COUNT,
        "Parent failed-artifact count changed",
    )
    require(len(report["parent_false_flags"]) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag count changed")
    require(
        report["parent_false_flag_values"] == {name: False for name in report["parent_false_flags"]},
        "Parent false-flag values changed",
    )
    for key in (
        "source_pins",
        "protected_file_hashes",
        "legacy_dirty_path_hashes",
        "preserved_failed_artifact_hashes_union",
    ):
        require(all(valid_sha(value) for value in report[key].values()), f"Malformed parent {key}")
    return {
        "report": report,
        "result": result,
        "report_raw": report_raw,
        "result_raw": result_raw,
        "report_meta": parent_pointer(PARENT_REPORT, PARENT_REPORT_SHA256, report_raw),
        "result_meta": parent_pointer(PARENT_RESULT, PARENT_RESULT_SHA256, result_raw),
        "report_view": parent_report_view(report),
    }


def _contract(api: Mapping[str, Any]) -> dict[str, Any]:
    input_contract = api.get("input_contract")
    result_contract = api.get("result_contract")
    diagnostic = api.get("diagnostic", {})
    require(isinstance(input_contract, Mapping), "API input_contract is missing")
    require(isinstance(result_contract, Mapping), "API result_contract is missing")
    input_schema = input_contract.get("schema")
    result_schema = result_contract.get("schema")
    require(isinstance(input_schema, str) and input_schema, "API input schema is missing")
    require(isinstance(result_schema, str) and result_schema, "API result schema is missing")
    saved = input_contract.get("saved_sources")
    if isinstance(saved, Mapping):
        wrapper_names = list(saved)
        wrapper_specs = clone(dict(saved))
    elif isinstance(saved, list):
        wrapper_names = [str(item) for item in saved]
        wrapper_specs = {}
    else:
        shape = input_contract.get("source_shape", {})
        names = shape.get("saved_sources") if isinstance(shape, Mapping) else None
        require(isinstance(names, list) and all(isinstance(item, str) for item in names), "API literal-wrapper names are missing")
        wrapper_names = list(names)
        wrapper_specs = {}
    require(wrapper_names and len(set(wrapper_names)) == len(wrapper_names), "API wrapper names are not unique")
    wrapper = input_contract.get("literal_wrapper", {})
    required_wrapper_keys = wrapper.get("required_keys") if isinstance(wrapper, Mapping) else None
    if required_wrapper_keys is None:
        required_wrapper_keys = ["literal", "raw_sha256"]
    require(list(required_wrapper_keys) == ["literal", "raw_sha256"], "Input wrapper contract is not exact")
    counts = result_contract.get("success_counts_exact")
    if counts is None:
        counts = result_contract.get("expected_counts")
    require(isinstance(counts, Mapping), "API success counts are missing")
    require(all(isinstance(key, str) and isinstance(value, int) and not isinstance(value, bool) and value >= 0 for key, value in counts.items()), "API counts are malformed")
    flags_true = result_contract.get("success_flags_true", [])
    flags_false = result_contract.get("explicit_false_flags", result_contract.get("success_flags_false", []))
    require(isinstance(flags_true, list) and all(isinstance(item, str) for item in flags_true), "API true flags are malformed")
    require(isinstance(flags_false, list) and all(isinstance(item, str) for item in flags_false), "API false flags are malformed")
    require(not set(flags_true) & set(flags_false), "API flag sets overlap")
    events = diagnostic.get("success_event_sequence", diagnostic.get("event_sequence"))
    multiplicity = diagnostic.get("success_event_multiplicity", diagnostic.get("event_multiplicity"))
    if events is not None:
        require(isinstance(events, list) and all(isinstance(item, str) for item in events), "API event sequence is malformed")
        require(isinstance(multiplicity, Mapping), "API event multiplicity is missing")
        require(set(multiplicity) == set(events), "API event names and multiplicity differ")
        require(all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in multiplicity.values()), "API event multiplicity is malformed")
        target_count = diagnostic.get("success_event_count", diagnostic.get("target_count"))
        if target_count is not None:
            require(sum(multiplicity.values()) == target_count, "API event target count changed")
    else:
        events, multiplicity = [], {}
    return {
        "input": clone(dict(input_contract)),
        "result": clone(dict(result_contract)),
        "diagnostic": clone(dict(diagnostic)) if isinstance(diagnostic, Mapping) else {},
        "input_schema": input_schema,
        "result_schema": result_schema,
        "publication_schema": api.get("publication_schema", result_contract.get("publication_schema", PUBLICATION_SCHEMA)),
        "wrapper_names": wrapper_names,
        "wrapper_specs": wrapper_specs,
        "counts": dict(counts),
        "flags_true": list(flags_true),
        "flags_false": list(flags_false),
        "events": list(events),
        "event_counts": dict(multiplicity),
        "new_flag_policy": clone(result_contract.get("new_flag_policy", api.get("new_flag_policy", {}))),
    }


def load_api(path: Path = API_PATH, expected: str | None = None) -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    raw, api = read_json(path, expected)
    require(api.get("component") == COMPONENT, "API component changed")
    artifact_meta = api.get("artifact", {})
    require(
        isinstance(artifact_meta, Mapping)
        and artifact_meta.get("status") in {"draft_unexecuted", "static_draft_unexecuted"}
        and artifact_meta.get("runtime_model_metadata") is None,
        "API is not an inert draft contract",
    )
    require(artifact_meta.get("import_or_execution_performed_during_authoring") is False, "API authoring execution marker changed")
    contract = _contract(api)
    require(contract["input_schema"] == INPUT_SCHEMA, "API input schema changed")
    return raw, api, contract


def load_facts(api_contract: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    values: dict[str, Any] = {}
    artifacts: dict[str, dict[str, Any]] = {}
    for name, (path, digest, schema) in FACT_SPECS.items():
        raw, value = read_json(path, digest)
        require(value.get("schema") == schema, f"Static fact schema changed: {name}")
        expected_status = "raw_extraction_only" if name == "raw_capsule" else "static_only_uncomputed"
        require(value.get("status") == expected_status, f"Static fact status changed: {name}")
        require(value.get("runtime_model_metadata") is None, f"Static fact runtime metadata changed: {name}")
        audit = value.get("audit", {})
        if isinstance(audit, Mapping):
            for key in (
                "source_modules_imported",
                "source_modules_executed",
                "project_math_providers_executed",
                "checkers_normalizers_or_old_helpers_executed",
                "mathematical_function_evaluation",
                "free_variable_union_or_zero_derivative_analysis",
                "calculus_or_simplification",
                "bounds_promoted_to_values",
                "N_node_18_materialized",
            ):
                if key in audit:
                    require(audit[key] is False, f"Static fact audit marker changed: {name}.{key}")
        if name == "raw_capsule":
            boundary = value.get("extraction_boundary", {})
            require(isinstance(boundary, Mapping), "Raw capsule extraction boundary is missing")
            for key in (
                "project_modules_imported",
                "project_modules_executed",
                "oldhelpers_or_checkers_executed",
                "proof_or_identity_analysis",
                "free_variable_analysis",
                "fractional_part_differentiated",
                "N_arithmetic",
                "N_materialized",
            ):
                require(boundary.get(key) is False, f"Raw capsule boundary changed: {key}")
            require(boundary.get("raw_record_dictionaries_unmodified") is True, "Raw capsule records were normalized")
            require(boundary.get("annotated_fields_added_inside_raw_records") is False, "Raw capsule records were annotated")
            union = value.get("record_union", {})
            require(union.get("raw_record_count") == 49, "Raw capsule record count changed")
            packet_inputs = value.get("source_packet_inputs")
            require(isinstance(packet_inputs, list) and len(packet_inputs) == 2, "Raw capsule source-packet count changed")
            require(
                {item.get("raw_sha256") for item in packet_inputs if isinstance(item, Mapping)}
                == {FACT_SPECS["logC_logP"][1], FACT_SPECS["Ra_phase_origin"][1]},
                "Raw capsule source-packet hashes changed",
            )
        values[name] = value
        artifacts[name] = artifact_from_bytes(path, raw, embed=True)
    return values, artifacts


def ast_artifact(path: Path, raw: bytes) -> dict[str, Any]:
    """Hash and syntax-parse source without importing or executing it."""
    text = raw.decode("utf-8")
    tree = ast.parse(text, filename=str(path), type_comments=True)
    dump = ast.dump(tree, annotate_fields=True, include_attributes=True).encode("utf-8")
    return {
        **artifact_from_bytes(path, raw, embed=True),
        "ast_dump_sha256": sha_bytes(dump),
        "syntax_parse_only": True,
        "runtime_model_metadata": None,
    }


def literal_wrapper(literal: Any, raw: bytes) -> dict[str, Any]:
    """The input contract is deliberately exactly two keys."""
    return {"literal": clone(literal), "raw_sha256": sha_bytes(raw)}


def wrapper_name_map(api_contract: Mapping[str, Any]) -> dict[str, str]:
    names = list(api_contract["wrapper_names"])
    require(set(names) == set(EXPECTED_WRAPPER_NAMES), "API literal-wrapper names changed")
    return {
        "accepted_result": "saved_metric_result",
        "logC_logP": "scalar_literals",
        "Ra_phase_origin": "phase_literals",
        "raw_capsule": "raw_capsule",
    }


def certificate_context(api_contract: Mapping[str, Any]) -> dict[str, Any]:
    requested = api_contract["input"].get("requested_context_keys")
    context = {
        "source_only_definitional_certificate": True,
        "ordinary_Z_independence_certificate": True,
        "closure_checked_once": False,
        "operation_semantics_checked_once": False,
        "binding_variables_checked_once": False,
        "context_is_request_policy_not_completed_proof": True,
        "naked_Z_independent_substitution_not_trusted": True,
        "old_dependency_walker": False,
        "old_metric_coefficients_replayed": False,
        "N_materialized": False,
        "phase_materialized": False,
        "raw_capsule_records_copied": False,
        "numeric_values_computed": False,
    }
    if isinstance(requested, Mapping):
        for key, value in requested.items():
            context.setdefault(str(key), clone(value) if isinstance(value, bool) else True)
    elif isinstance(requested, list):
        for key in requested:
            if isinstance(key, str):
                context.setdefault(key, True)
    return context


def build_input(parent: Mapping[str, Any], facts: Mapping[str, Any], fact_artifacts: Mapping[str, Any], api_contract: Mapping[str, Any]) -> dict[str, Any]:
    names = wrapper_name_map(api_contract)
    saved_sources = {
        names["accepted_result"]: literal_wrapper(parent["result"], parent["result_raw"]),
        names["logC_logP"]: literal_wrapper(facts["logC_logP"], read_bytes(FACT_SPECS["logC_logP"][0], FACT_SPECS["logC_logP"][1])),
        names["Ra_phase_origin"]: literal_wrapper(facts["Ra_phase_origin"], read_bytes(FACT_SPECS["Ra_phase_origin"][0], FACT_SPECS["Ra_phase_origin"][1])),
        names["raw_capsule"]: literal_wrapper(facts["raw_capsule"], read_bytes(FACT_SPECS["raw_capsule"][0], FACT_SPECS["raw_capsule"][1])),
    }
    require(all(set(value) == {"literal", "raw_sha256"} for value in saved_sources.values()), "Input wrappers are not exact two-key records")
    return {
        "schema": api_contract["input_schema"],
        "component": COMPONENT,
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "saved_sources": saved_sources,
        "requested_context": certificate_context(api_contract),
    }


def _fact_handles(facts: Mapping[str, Any]) -> dict[str, Any]:
    """Retain only source metadata needed to describe the future consumer."""
    result: dict[str, Any] = {}
    phase = facts.get("Ra_phase_origin", {})
    bindings = phase.get("native_bindings") if isinstance(phase, Mapping) else None
    graph = phase.get("literal_graph_nodes") if isinstance(phase, Mapping) else None
    if isinstance(bindings, Mapping):
        for key in ("chart", "ordinary_Z_symbol", "Rm", "Ra", "phase_origin", "N", "actual_phase"):
            if key in bindings:
                result[key] = clone(bindings[key])
    if isinstance(graph, Mapping):
        for key in ("N_18", "Rm_endpoint_85743", "native_radius_241787", "actual_phase_242551"):
            if key in graph:
                result[key] = clone(graph[key])
    capsule = facts.get("raw_capsule", {})
    if isinstance(capsule, Mapping):
        union = capsule.get("record_union", {})
        result["raw_capsule"] = {
            "schema": capsule.get("schema"),
            "status": capsule.get("status"),
            "raw_record_count": union.get("raw_record_count"),
            "raw_operation_counts": clone(union.get("raw_operation_counts")),
            "union_record_ids": clone(union.get("union_record_ids")),
            "raw_record_dictionaries_unmodified": capsule.get("extraction_boundary", {}).get("raw_record_dictionaries_unmodified"),
        }
    return result


def build_source(
    parent: Mapping[str, Any],
    input_value: Mapping[str, Any],
    input_raw: bytes,
    api_raw: bytes,
    api: Mapping[str, Any],
    api_contract: Mapping[str, Any],
    fact_artifacts: Mapping[str, Mapping[str, Any]],
    kernel_artifact: Mapping[str, Any],
    runner_artifact: Mapping[str, Any],
    main_brief_preimage: Mapping[str, Any],
) -> dict[str, Any]:
    additions = {
        wire(INPUT_SPEC): sha_bytes(input_raw),
        wire(API_PATH): sha_bytes(api_raw),
        wire(KERNEL_PATH): kernel_artifact["sha256"],
        wire(RUNNER_PATH): runner_artifact["sha256"],
    }
    additions.update({item["path"]: item["sha256"] for item in fact_artifacts.values()})
    source_union = merge_hashes(parent["report"]["source_pins"], additions)
    return {
        "schema": api_contract["input_schema"],
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "source_pins": clone(parent["report"]["source_pins"]),
        "protected_file_hashes": clone(parent["report"]["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(parent["report"]["legacy_dirty_path_hashes"]),
        "preserved_failed_artifact_hashes_union": clone(parent["report"]["preserved_failed_artifact_hashes_union"]),
        "parent_false_flags": clone(parent["report"]["parent_false_flags"]),
        "parent_false_flag_values": clone(parent["report"]["parent_false_flag_values"]),
        "preserved_error_cap_metadata": clone(parent["report"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(parent["report"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(parent["report"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(parent["report"]["retained_cap_provenance"]),
        "parent_new_flags": clone(parent["report"].get("new_flags", {})),
        "input_literal": clone(dict(input_value)),
        "input_raw_sha256": sha_bytes(input_raw),
        "input_raw_bytes_base64": base64.b64encode(input_raw).decode("ascii"),
        "input_spec": {"path": wire(INPUT_SPEC), "sha256": sha_bytes(input_raw), "bytes": len(input_raw)},
        "api": artifact_from_bytes(API_PATH, api_raw, embed=True),
        "api_contract": clone(dict(api_contract)),
        "api_source": clone(dict(api)),
        "code2": {"kernel": clone(dict(kernel_artifact)), "runner": clone(dict(runner_artifact))},
        "static_fact_artifacts": clone(dict(fact_artifacts)),
        "certificate_fact_handles": _fact_handles({name: json.loads(base64.b64decode(item["raw_bytes_base64"]).decode("utf-8")) for name, item in fact_artifacts.items()}),
        "expected_counts": clone(api_contract["counts"]),
        "expected_flags_true": clone(api_contract["flags_true"]),
        "expected_flags_false": clone(api_contract["flags_false"]),
        "callback_event_names": clone(api_contract["events"]),
        "callback_event_counts": clone(api_contract["event_counts"]),
        "dynamic_api_targets": clone(api_contract["diagnostic"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "main_brief_preimage": clone(dict(main_brief_preimage)),
        "no_project_imports_before_run": True,
        "no_provider_calls": True,
        "no_kernel_execution": True,
        "input_wrapper_names": clone(wrapper_name_map(api_contract)),
        "certificate_scope": "source-only definitional certificate; not numeric closure",
        "consumer_guards": {
            "closure_checked_once": False,
            "operation_semantics_checked_once": False,
            "binding_variables_checked_once": False,
            "context_is_request_policy_not_completed_proof": True,
            "naked_Z_independent_substitution_rejected": True,
            "old_dependency_walker": False,
            "old_metric_coefficients_replayed": False,
        },
    }


def load_source(*, require_empty_journal: bool = False) -> tuple[bytes, dict[str, Any]]:
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == source.get("api_contract", {}).get("input_schema"), "Source checkpoint schema changed")
    require(source.get("component") == COMPONENT and source.get("parent_commit") == PARENT_COMMIT, "Source identity changed")
    input_raw = base64.b64decode(source["input_raw_bytes_base64"])
    require(sha_bytes(input_raw) == source["input_raw_sha256"], "Frozen input bytes changed")
    require(json.loads(input_raw.decode("utf-8")) == source["input_literal"], "Frozen input literal changed")
    read_bytes(INPUT_SPEC, source["input_raw_sha256"])
    read_bytes(API_PATH, source["api"]["sha256"])
    for item in source["code2"].values():
        read_bytes(resolve_path(item["path"], REPO), item["sha256"])
    for item in source["static_fact_artifacts"].values():
        read_bytes(resolve_path(item["path"], WORK), item["sha256"])
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    require(all(set(value) == {"literal", "raw_sha256"} for value in source["input_literal"]["saved_sources"].values()), "Frozen input wrappers changed")
    if require_empty_journal:
        require(JOURNAL_CHECKPOINT.exists() and JOURNAL_CHECKPOINT.read_bytes() == b"", "Run requires an empty durable journal")
    return raw, source


def append_journal(event: Mapping[str, Any]) -> bytes:
    raw = encode(dict(event))
    JOURNAL_CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    with JOURNAL_CHECKPOINT.open("ab") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def journal_events() -> list[dict[str, Any]]:
    if not JOURNAL_CHECKPOINT.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in JOURNAL_CHECKPOINT.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "Journal rows must be JSON objects")
            rows.append(value)
    return rows


def callback_summary(events: Sequence[Mapping[str, Any]], source: Mapping[str, Any]) -> dict[str, Any]:
    require(len(events) >= 2, "Journal prefix is incomplete")
    require(events[0].get("event") == "manager_entry_started", "Journal entry prefix changed")
    require(events[1].get("event") in {"manager_kernel_imported", "manager_kernel_or_derive_failed"}, "Journal manager prefix changed")
    names = set(source.get("callback_event_names", []))
    counts: dict[str, int] = {}
    rows: list[dict[str, Any]] = []
    for row in events[2:]:
        event = row.get("event")
        if str(event).startswith("manager_"):
            continue
        require(event in names, f"Unexpected callback event: {event}")
        counts[event] = counts.get(event, 0) + 1
        rows.append(dict(row))
    return {
        "manager_prefix_event_count": 2,
        "callback_event_count": len(rows),
        "callback_event_counts": counts,
        "callback_event_raw_rows": rows,
        "journal_event_count": len(events),
        "journal_event_count_expected": 2 + sum(source.get("callback_event_counts", {}).values()),
    }


def zero_counts(source: Mapping[str, Any]) -> dict[str, int]:
    return {key: 0 for key in source.get("expected_counts", {})}


def recovered_counts(source: Mapping[str, Any]) -> dict[str, int]:
    """Recover the durable callback prefix without resetting observed counts."""
    counts = zero_counts(source)
    diagnostic = source.get("api_contract", {}).get("diagnostic", {})
    mapping = diagnostic.get("event_count_keys", diagnostic.get("count_key_by_event", {})) if isinstance(diagnostic, Mapping) else {}
    if not isinstance(mapping, Mapping):
        mapping = {}
    for row in journal_events()[2:]:
        event = row.get("event")
        key = mapping.get(event)
        if key not in counts:
            normalized = str(event).lower().replace("-", "_")
            candidates = [
                normalized,
                normalized.removesuffix("_row"),
                normalized.removesuffix("_entry"),
                normalized + "s",
            ]
            key = next((candidate for candidate in candidates if candidate in counts), None)
        if key in counts:
            counts[key] += 1
    return counts


def false_flags(source: Mapping[str, Any]) -> dict[str, bool]:
    return {name: False for name in (*source.get("expected_flags_true", []), *source.get("expected_flags_false", []))}


def _source_contract_violations(component: Mapping[str, Any]) -> list[dict[str, Any]]:
    violations: list[dict[str, Any]] = []
    payload = component.get("payload")
    if not isinstance(payload, Mapping):
        return violations
    contract = payload.get("source_contract")
    if isinstance(contract, Mapping):
        for key in ("project_modules_imported", "source_files_opened", "source_values_consumed", "old_math_called"):
            if key in contract and contract[key] is not False:
                violations.append({"field": f"payload.source_contract.{key}", "reason": "source-only marker changed"})
    return violations


def _flag_policy_violations(flags: Mapping[str, Any], source: Mapping[str, Any]) -> list[dict[str, Any]]:
    violations: list[dict[str, Any]] = []
    parent_values = source.get("parent_false_flag_values", {})
    for name in parent_values:
        if name in flags and flags[name] is not False:
            violations.append({"field": f"new_flags.{name}", "reason": "inherited global false flag was promoted"})
    parent_new = source.get("parent_new_flags", {})
    for name, value in parent_new.items():
        if name in flags and flags[name] != value:
            violations.append({"field": f"new_flags.{name}", "reason": "parent historical new flag changed"})
    policy = source.get("api_contract", {}).get("new_flag_policy", {})
    allowed = policy.get("allowed_true_flags") if isinstance(policy, Mapping) else None
    if isinstance(allowed, list):
        for name, value in flags.items():
            if value is True and name not in parent_new and name not in allowed:
                violations.append({"field": f"new_flags.{name}", "reason": "API disallows promotion of this flag"})
    else:
        new_true = [name for name, value in flags.items() if value is True and name not in parent_new]
        if len(new_true) > 1:
            violations.append({"field": "new_flags", "reason": "more than one new flag was promoted"})
        if new_true and not any(
            "z" in name.lower()
            and any(token in name.lower() for token in ("independ", "zero", "cert"))
            for name in new_true
        ):
            violations.append({"field": "new_flags", "reason": "the only promotable new flag must be the Rm/Z certificate"})
    return violations


def validate_result(result: Mapping[str, Any], source: Mapping[str, Any], events: Sequence[Mapping[str, Any]], derive_entry_calls: Any, derive_calls: Any) -> dict[str, Any]:
    component = result.get("component_result") if isinstance(result.get("component_result"), Mapping) else result
    violations: list[dict[str, Any]] = []
    result_contract = source.get("api_contract", {}).get("result", {})
    expected_schema = result_contract.get("schema")
    for key, expected in (
        ("schema", expected_schema),
        ("component", COMPONENT),
        ("metadata_only", result_contract.get("metadata_only")),
        ("runtime_model_metadata", None),
        ("source_only", True),
        ("conditional", result_contract.get("conditional")),
        ("physical_layer_completed", False),
        ("full_goal_complete", False),
    ):
        if expected is not None and component.get(key) != expected:
            violations.append({"field": key, "expected": expected, "actual": component.get(key)})
    required = result_contract.get("required_top_level_fields", [])
    if isinstance(required, list):
        for key in required:
            if key not in component:
                violations.append({"field": key, "reason": "API-required result field is missing"})
    counts = component.get("actual_counts")
    expected_counts = source.get("expected_counts", {})
    if not isinstance(counts, Mapping) or set(counts) != set(expected_counts):
        violations.append({"field": "actual_counts", "reason": "dynamic API count keys are required"})
    else:
        for key, value in counts.items():
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                violations.append({"field": f"actual_counts.{key}", "reason": "count is not a nonnegative integer"})
    flags = component.get("new_flags")
    expected_flags = set(source.get("expected_flags_true", [])) | set(source.get("expected_flags_false", []))
    if not isinstance(flags, Mapping) or set(flags) != expected_flags or not all(isinstance(flags.get(name), bool) for name in expected_flags):
        violations.append({"field": "new_flags", "reason": "dynamic API flag keys are required"})
    else:
        violations.extend(_flag_policy_violations(flags, source))
    violations.extend(_source_contract_violations(component))
    try:
        summary = callback_summary(events, source)
    except AssertionError as exc:
        summary = {
            "manager_prefix_event_count": 0,
            "callback_event_count": 0,
            "callback_event_counts": {},
            "callback_event_raw_rows": [],
            "journal_event_count": len(events),
            "journal_event_count_expected": None,
        }
        violations.append({"field": "journal", "reason": str(exc)})
    passed = component.get("passed") is True
    if passed:
        if expected_counts and counts != expected_counts:
            violations.append({"field": "actual_counts", "reason": "counts do not match the frozen API"})
        if expected_flags and any(flags.get(name) is not True for name in source.get("expected_flags_true", [])):
            violations.append({"field": "new_flags", "reason": "API-declared true flags are not true"})
        if expected_flags and any(flags.get(name) is not False for name in source.get("expected_flags_false", [])):
            violations.append({"field": "new_flags", "reason": "API-declared false flags are not false"})
        if source.get("callback_event_counts") and summary.get("callback_event_counts") != source["callback_event_counts"]:
            violations.append({"field": "callback_events", "reason": "callback multiplicity does not match API"})
        if derive_entry_calls != 1 or derive_calls != 1:
            violations.append({"field": "derive_calls", "reason": "successful kernel entry must record exactly one entry and derive"})
    return {
        "schema": f"{COMPONENT}_result_validation_v1",
        "component": COMPONENT,
        "passed": passed and not violations,
        "explicit_passed": component.get("passed"),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "actual_counts": clone(counts),
        "new_flags": clone(flags),
        "callback": summary,
        "violations": violations,
        "raw_result_preserved": True,
        "runtime_model_metadata": None,
    }


def import_kernel(path: Path) -> Any:
    """Import only the frozen new kernel, and only from the entered run."""
    require(path.resolve() != Path(__file__).resolve(), "Kernel path aliases the metadata manager")
    module_name = f"_{COMPONENT}_entry"
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, "Cannot load the frozen new kernel")
    module = importlib.util.module_from_spec(spec)
    # Registration is durable in the in-process import table before execution.
    sys.modules[module_name] = module
    require(sys.modules.get(module_name) is module, "Kernel sys.modules registration failed")
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "Frozen new kernel has no derive API")
    return module


def _failure_result(source: Mapping[str, Any], error: BaseException, kernel_imported: bool) -> dict[str, Any]:
    return {
        "schema": source.get("api_contract", {}).get("result_schema", RESULT_SCHEMA),
        "component": COMPONENT,
        "passed": False,
        "runtime_model_metadata": None,
        "metadata_only": False,
        "source_only": True,
        "conditional": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": recovered_counts(source),
        "new_flags": false_flags(source),
        "payload": {"recovered_callback_prefix": journal_events()[2:]},
        "diagnostics": {"callback_failure_stops_later_work": True, "event_count": len(journal_events())},
        "failure": {"kind": "manager_import_or_derive_failure", "error_type": type(error).__name__, "error": repr(error), "kernel_imported": kernel_imported, "stopped": True},
    }


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    require(not any(path.exists() for path in (INPUT_SPEC, BRIEF, PREPARE_RECEIPT, RUN_RECEIPT, PUBLISH_RECEIPT, *CHECKPOINTS)), "This component already has a durable attempt")
    parent = load_parent()
    current_raw, _ = load_current_brief(exact=True)
    api_raw, api, api_contract = load_api(API_PATH, args.api_sha256.lower() if args.api_sha256 else None)
    facts, fact_artifacts = load_facts(api_contract)
    kernel_raw = read_bytes(REPO / KERNEL_RELATIVE, args.kernel_sha256.lower() if args.kernel_sha256 else None)
    runner_raw = read_bytes(REPO / RUNNER_RELATIVE, args.runner_sha256.lower() if args.runner_sha256 else None)
    # AST parsing is the only inspection of code2.  It is not an import or run.
    kernel_artifact = ast_artifact(REPO / KERNEL_RELATIVE, kernel_raw)
    runner_artifact = ast_artifact(REPO / RUNNER_RELATIVE, runner_raw)
    input_value = build_input(parent, facts, fact_artifacts, api_contract)
    input_raw = write_new(INPUT_SPEC, input_value)
    preimage = brief_preimage(current_raw)
    preimage_path = Path(preimage["path"])
    require(not preimage_path.exists(), "Main brief preimage already exists; metadata repair is external")
    write_new(preimage_path, current_raw)
    source = build_source(parent, input_value, input_raw, api_raw, api, api_contract, fact_artifacts, kernel_artifact, runner_artifact, preimage)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    journal_raw = write_new(JOURNAL_CHECKPOINT, b"")
    pending = {
        "schema": f"{COMPONENT}_brief_v1",
        "component": COMPONENT,
        "phase": "prepared",
        "source_checkpoint": {"path": wire(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "api": artifact_from_bytes(API_PATH, api_raw, embed=False),
        "code2": {"kernel": {key: value for key, value in kernel_artifact.items() if key != "raw_bytes_base64"}, "runner": {key: value for key, value in runner_artifact.items() if key != "raw_bytes_base64"}},
        "facts": {name: {key: value for key, value in item.items() if key != "raw_bytes_base64"} for name, item in fact_artifacts.items()},
        "journal": {"path": wire(JOURNAL_CHECKPOINT), "sha256": sha_bytes(journal_raw)},
        "metadata_only": True,
        "source_only": True,
        "do_not_reenter": True,
        "no_kernel_execution": True,
    }
    pending_raw = write_new(BRIEF, pending)
    main_brief = update_main_brief({
        "phase": "prepared",
        "current_component": COMPONENT,
        "in_progress": {COMPONENT: {"brief": wire(BRIEF), "source": wire(SOURCE_CHECKPOINT), "source_sha256": sha_bytes(source_raw), "metadata_only": True, "source_only": True}},
        "goal_active": True,
        "full_goal_complete": False,
        "next_task": "Root may run the frozen native parameter-Z kernel exactly once after entry; this manager does not execute it during authoring.",
    })
    receipt = {
        "schema": f"{COMPONENT}_prepare_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_raw),
        "journal_checkpoint_sha256": sha_bytes(journal_raw),
        "brief_sha256": sha_bytes(pending_raw),
        "main_brief": main_brief,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "source_union_count": source["source_union_count"],
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "no_kernel_execution": True,
    }
    receipt_raw = write_new(PREPARE_RECEIPT, receipt)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw)}


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Enter once, import the frozen new kernel, and derive exactly once."""
    source_raw, source = load_source(require_empty_journal=True)
    require(PREPARE_RECEIPT.exists(), "Prepare receipt is missing")
    _, prepare_receipt = read_json(PREPARE_RECEIPT)
    require(prepare_receipt.get("source_checkpoint_sha256") == sha_bytes(source_raw), "Prepare receipt does not authenticate source")
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not RUN_RECEIPT.exists(), "Run has already been entered")
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "kernel_path": wire(REPO / KERNEL_RELATIVE),
        "kernel_sha256": source["code2"]["kernel"]["sha256"],
        "import_attempts": 1,
        "import_attempt_recorded_before_import": True,
        "sys_modules_registration_recorded_before_import": True,
        "derive_entry_calls": 0,
        "derive_entry_reserved": True,
        "derive_calls": 0,
        "kernel_imported": False,
        "runtime_model_metadata": None,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
        "do_not_reenter": True,
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    append_journal({
        "schema": DIAGNOSTIC_SCHEMA,
        "event": "manager_entry_started",
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "import_attempts": 1,
        "derive_entry_calls": 0,
        "derive_calls": 0,
        "derive_entry_reserved": True,
    })
    component_result: dict[str, Any]
    kernel_imported = False
    derive_entry_calls = 0
    derive_calls = 0
    try:
        module = import_kernel(KERNEL_PATH)
        kernel_imported = True
        append_journal({
            "schema": DIAGNOSTIC_SCHEMA,
            "event": "manager_kernel_imported",
            "kernel_sha256": source["code2"]["kernel"]["sha256"],
            "sys_modules_registered": True,
        })
        def callback(row: Mapping[str, Any]) -> None:
            require(isinstance(row, Mapping), "Kernel callback row must be a mapping")
            require(row.get("event") in source["callback_event_names"], "Kernel callback event is not in the frozen API")
            append_journal(dict(row))
        derive_entry_calls = 1
        derive_calls = 1
        component_result = module.derive(clone(source["input_literal"]), callback)
        require(isinstance(component_result, dict), "Kernel derive did not return a JSON mapping")
    except BaseException as exc:
        append_journal({"schema": DIAGNOSTIC_SCHEMA, "event": "manager_kernel_or_derive_failed", "error_type": type(exc).__name__, "error": repr(exc), "kernel_imported": kernel_imported, "derive_calls": derive_calls})
        component_result = _failure_result(source, exc, kernel_imported)
    outer = clone(component_result)
    outer["component_result"] = clone(component_result)
    outer["component_result_exact"] = True
    outer["derive_entry_calls"] = derive_entry_calls
    outer["derive_calls"] = derive_calls
    outer["import_or_execution_performed"] = kernel_imported
    outer["execution_metadata"] = {
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "import_attempts": 1,
        "kernel_imported": kernel_imported,
        "runner_executed": False,
        "old_math_reexecuted": False,
    }
    outer["kernel_imported"] = kernel_imported
    outer["passed"] = component_result.get("passed") is True
    outer["actual_counts"] = clone(component_result.get("actual_counts", recovered_counts(source)))
    outer["new_flags"] = clone(component_result.get("new_flags", false_flags(source))) if outer["passed"] else false_flags(source)
    outer["metadata_only"] = False
    outer["source_only"] = True
    result_raw = write_new(RESULT_CHECKPOINT, outer)
    validation = validate_result(outer, source, journal_events(), derive_entry_calls, derive_calls)
    receipt = {
        "schema": f"{COMPONENT}_run_receipt_v1",
        "component": COMPONENT,
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "journal_checkpoint_sha256": sha_bytes(JOURNAL_CHECKPOINT.read_bytes()),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "import_attempts": 1,
        "kernel_imported": kernel_imported,
        "import_or_execution_performed": kernel_imported,
        "durable_result_saved_before_validation": True,
        "callback_prefix_recovered_from_journal": True,
        "validation": validation,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "do_not_reenter": True,
    }
    receipt_raw = write_new(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}


def _review_scope_ok(review: Mapping[str, Any]) -> bool:
    """Accept the actual reviewer's scoped marker, including dict scopes."""
    prohibited = ("numeric", "physical", "full_goal", "closure", "metric_replay")
    tokens = ("parameter_z", "ordinary_z", "rm_z", "z_independ", "definitional", COMPONENT.lower())
    candidates: list[Any] = [review.get("scope"), review.get("review_scope"), review.get("reviewed_scope"), review.get("accepted_scope"), review.get("scope_marker")]
    accepted = False
    blocked = False
    def visit(value: Any, key: str = "") -> None:
        nonlocal accepted, blocked
        if isinstance(value, Mapping):
            for child_key, child_value in value.items():
                visit(child_value, str(child_key).lower())
        elif isinstance(value, list):
            for child in value:
                visit(child, key)
        elif isinstance(value, str):
            text = (key + " " + value).lower()
            if any(token in text for token in tokens):
                accepted = True
        elif value is True:
            if any(token in key for token in tokens):
                accepted = True
            if any(token in key for token in prohibited):
                blocked = True
    for item in candidates:
        visit(item)
    return accepted and not blocked


def load_review(path: str, digest: str) -> dict[str, Any]:
    require(valid_sha(digest), "Saved review SHA is malformed")
    review_path = resolve_path(path, WORK)
    raw, value = read_json(review_path, digest.lower())
    reviewer = value.get("reviewer", {})
    independent = (
        value.get("independent_review") is True
        or value.get("independent") is True
        or (isinstance(reviewer, Mapping) and reviewer.get("independent") is True)
    )
    passed = value.get("passed") is True or value.get("review_passed") is True
    require(independent and passed, "Saved review is not independent and passed")
    require(_review_scope_ok(value), "Saved review scope does not identify the parameter-Z certificate")
    return {
        "path": wire(review_path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "schema": value.get("schema"),
        "passed": True,
        "independent_review": True,
        "review_scope": clone(value.get("scope", value.get("review_scope", value.get("reviewed_scope")))),
        "raw_bytes_exact": True,
        "raw_readonly": True,
        "raw_bytes_base64": base64.b64encode(raw).decode("ascii"),
        "value": value,
    }


def attrs_candidate() -> bytes:
    path = REPO / ".gitattributes"
    old = path.read_bytes() if path.exists() else b""
    text = old.decode("utf-8")
    lines = [
        f"/{KERNEL_RELATIVE.as_posix()} -text",
        f"/{RUNNER_RELATIVE.as_posix()} -text",
        f"/{REPORT_RELATIVE.as_posix()} -text",
        f"/{DOC_RELATIVE.as_posix()} -text",
    ]
    existing = text.splitlines()
    require(all(existing.count(line) <= 1 for line in lines), "Duplicate scope attribute already exists")
    for line in lines:
        if line not in existing:
            text += ("\n" if text and not text.endswith("\n") else "") + line + "\n"
            existing.append(line)
    return text.encode("utf-8")


def document_bytes(validation: Mapping[str, Any], result: Mapping[str, Any]) -> bytes:
    counts = result.get("actual_counts")
    flags = result.get("new_flags")
    return (
        "# Current native parameter-Z independence certificate (2026-10-27)\n\n"
        f"Component: `{COMPONENT}`.\n\n"
        f"Attempt status: `{'passed' if validation.get('passed') else 'failed'}`. The scope is a source-only definitional certificate, not numeric closure or physical completion.\n\n"
        f"Dynamic API-declared counts: `{counts}`.\n\n"
        f"Dynamic API-declared flags: `{flags}`.\n\n"
        "The consumer records the exact native logC/logP and Ra/phase-origin literal packets, the raw graph capsule, the saved metric result, and the ordinary-Z binding context. N remains the opaque double-dyadic expression; Rm, Ra, phase, moving endpoints, and source rows remain data handles.\n\n"
        "The frozen input requests one closure, operation-semantics, and binding-variable check from the future consumer; those checks are not pre-certified by this metadata manager. A naked Z-independent-substitution marker is not trusted, no old dependency walker or metric/Jacobian/Laplacian family is replayed, and no source value, numeric scalar, phase identity, or physical closure is claimed.\n\n"
        "Inherited source pins, false flags, error caps, construction data, integral inventories, and cap provenance remain retained. The goal remains ACTIVE / INCOMPLETE; do not reenter this attempt.\n"
    ).encode("utf-8")


def publish(args: argparse.Namespace) -> dict[str, Any]:
    source_raw, source = load_source()
    report_path = REPO / REPORT_RELATIVE
    require(RUN_RECEIPT.exists() and not report_path.exists() and not PUBLISH_RECEIPT.exists(), "Publication already exists or run receipt is missing")
    _, run_receipt = read_json(RUN_RECEIPT)
    result_raw, result = read_json(RESULT_CHECKPOINT, run_receipt.get("result_checkpoint_sha256"))
    api_raw, api, api_contract = load_api(API_PATH, source["api"]["sha256"])
    validation = validate_result(result, source, journal_events(), run_receipt.get("derive_entry_calls"), run_receipt.get("derive_calls"))
    review = load_review(args.review, args.review_sha256)
    kernel_artifact = artifact(REPO / KERNEL_RELATIVE, source["code2"]["kernel"]["sha256"], embed=True)
    runner_artifact = artifact(REPO / RUNNER_RELATIVE, source["code2"]["runner"]["sha256"], embed=True)
    input_artifact = artifact(INPUT_SPEC, source["input_spec"]["sha256"], embed=True)
    api_artifact = artifact_from_bytes(API_PATH, api_raw, embed=True)
    fact_artifacts = {name: artifact(path, digest, embed=True) for name, (path, digest, _) in FACT_SPECS.items()}
    checkpoint_artifacts = [artifact(path, embed=True) for path in CHECKPOINTS]
    helper_artifact = artifact(Path(__file__), embed=True)
    require(not (REPO / DOC_RELATIVE).exists(), "Publication document already exists; metadata repair is external")
    handoff_old: dict[str, bytes] = {}
    handoff_preimages: dict[str, dict[str, Any]] = {}
    handoff_artifacts: dict[str, dict[str, Any]] = {}
    prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`.\n"
        f"Status: `{'passed' if validation['passed'] else 'failed'}`; source-only definitional certificate; do not reenter.\n"
    )
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        require(f"## {COMPONENT} publication" not in old.decode("utf-8", errors="ignore"), f"Existing handoff already contains {COMPONENT} publication")
        preimage_path = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(old)[:12]}.md"
        require(not preimage_path.exists(), "Handoff preimage already exists; metadata repair is external")
        handoff_old[wire(relative)] = old
        handoff_preimages[wire(relative)] = artifact_from_bytes(preimage_path, old, embed=False)
        handoff_artifacts[wire(relative)] = artifact_from_bytes(path, prefix.encode("utf-8") + old, embed=False)
    doc_raw = document_bytes(validation, result)
    attrs_raw = attrs_candidate()
    current_brief_raw, current_brief = read_json(MAIN_BRIEF)
    require(current_brief.get("execution_rules") == EXECUTION_RULES and current_brief.get("goal_active") is True and current_brief.get("full_goal_complete") is False, "Main brief state changed before publication")
    publication_brief_preimage = brief_preimage(current_brief_raw)
    checkpoints = {item["path"]: item["sha256"] for item in checkpoint_artifacts}
    published_pins = merge_hashes(source["source_union_hashes"], checkpoints)
    failed = clone(source["preserved_failed_artifact_hashes_union"])
    if not validation["passed"]:
        failed = merge_hashes(failed, {item["path"]: item["sha256"] for item in (kernel_artifact, runner_artifact, input_artifact, *checkpoint_artifacts)})
    next_task = api.get("workflow_handoff", {}).get("next_executable_task") if isinstance(api.get("workflow_handoff"), Mapping) else None
    if not isinstance(next_task, str) or not next_task:
        next_task = "root-supplied next task"
    report = {
        "schema": source["api_contract"].get("publication_schema", PUBLICATION_SCHEMA),
        "component": COMPONENT,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "passed": validation["passed"],
        "metadata_only": False,
        "publication_metadata_only": True,
        "source_only": True,
        "source_only_definitional_certificate": True,
        "numeric_closure": False,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "parent_report": source["parent_report"],
        "parent_result": source["parent_result"],
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "checkpoint_count": PUBLISHED_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "preserved_failed_artifact_count": len(failed),
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "source_pin_count": len(published_pins),
        "source_union_count": len(published_pins),
        "source_pins": published_pins,
        "raw_sha_pins": published_pins,
        "checkpoint_hashes": checkpoints,
        "preserved_failed_artifact_hashes_union": failed,
        "protected_file_hashes": clone(source["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(source["legacy_dirty_path_hashes"]),
        "parent_false_flags": clone(source["parent_false_flags"]),
        "parent_false_flag_values": clone(source["parent_false_flag_values"]),
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "preserved_error_cap_metadata": clone(source["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(source["retained_construction_data"]),
        "retained_integral_term_inventory": clone(source["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(source["retained_cap_provenance"]),
        "input": clone(source["input_literal"]),
        "input_spec": input_artifact,
        "code2": {"kernel": kernel_artifact, "runner": runner_artifact},
        "api_artifact": api_artifact,
        "static_fact_artifacts": fact_artifacts,
        "checkpoint_artifacts": checkpoint_artifacts,
        "current_review": review,
        "component_brief": artifact(BRIEF, embed=True),
        "main_brief_preimage": publication_brief_preimage,
        "metadata_helper_artifact": helper_artifact,
        "document_artifact": artifact_from_bytes(REPO / DOC_RELATIVE, doc_raw, embed=True),
        "handoff_preimages": handoff_preimages,
        "handoff_artifacts": handoff_artifacts,
        "gitattributes_required": {"paths": [wire(REPO / path) for path in (KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, DOC_RELATIVE)], "rule": "-text", "apply_before_freeze": True},
        "full_result": clone(result),
        "full_result_raw_sha256": sha_bytes(result_raw),
        "result_checkpoint": {"path": wire(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw), "bytes": len(result_raw)},
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": checkpoints[wire(ENTERED_CHECKPOINT)],
        "diagnostic_checkpoint_sha256": checkpoints[wire(JOURNAL_CHECKPOINT)],
        "actual_counts": validation["actual_counts"],
        "new_flags": validation["new_flags"],
        "expected_counts": clone(source["expected_counts"]),
        "expected_new_flags": {name: True for name in source["expected_flags_true"]} | {name: False for name in source["expected_flags_false"]},
        "callback_event_counts": clone(source["callback_event_counts"]),
        "callback_event_count": sum(source["callback_event_counts"].values()),
        "manager_prefix_event_count": 2,
        "journal_event_count": len(journal_events()),
        "durable_result_validation": validation,
        "report_scope_attributes": list(REPORT_SCOPE_ATTRIBUTES),
        "read_first": [wire(REPO / "AGENTS.md"), wire(REPO / DOC_RELATIVE), wire(REPO / REPORT_RELATIVE)],
        "current_doc_sha256": sha_bytes(doc_raw),
        "current_document_sha256": sha_bytes(doc_raw),
        "runtime_model_metadata": None,
        "do_not_reenter": True,
        "next_executable_task": next_task,
        # Dedicated fields above carry each raw artifact exactly once.  This
        # index is intentionally pointer-only so raw base64 is not duplicated.
        "lossless_artifact_pointers": [
            artifact_pointer(kernel_artifact),
            artifact_pointer(runner_artifact),
            artifact_pointer(input_artifact),
            artifact_pointer(api_artifact),
            *(artifact_pointer(item) for item in fact_artifacts.values()),
            *(artifact_pointer(item) for item in checkpoint_artifacts),
            artifact_pointer(review),
            artifact_pointer(helper_artifact),
        ],
        "raw_bytes_lossless": True,
    }
    # All validation is complete before any publication file is written.
    report_raw = write_new(report_path, report)
    write_new(REPO / DOC_RELATIVE, doc_raw)
    for relative, old in handoff_old.items():
        path = REPO / relative
        write_new(Path(handoff_preimages[relative]["path"]), old)
        write_replace(path, prefix.encode("utf-8") + old)
    write_replace(REPO / ".gitattributes", attrs_raw)
    brief = update_main_brief({
        "phase": "published_passed" if validation["passed"] else "published_failed",
        "current_component": COMPONENT,
        "report": wire(report_path),
        "report_sha256": sha_bytes(report_raw),
        "current_doc_sha256": sha_bytes(doc_raw),
        "current_document_sha256": sha_bytes(doc_raw),
        "read_first": report["read_first"],
        "goal_active": True,
        "full_goal_complete": False,
        "metadata_only": current_brief.get("metadata_only", True),
        "source_only": True,
        "next_task": next_task,
        "do_not_reenter": True,
    })
    receipt = {
        "schema": f"{COMPONENT}_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire(report_path),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "passed": report["passed"],
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "current_review": review,
        "main_resume_brief": brief,
        "checkpoint_count": PUBLISHED_CHECKPOINT_COUNT,
        "source_pin_count": len(published_pins),
        "preserved_failed_artifact_count": len(failed),
        "runtime_model_metadata": None,
        "do_not_reenter": True,
    }
    receipt_raw = write_new(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}


def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("--api-sha256")
    prepare_parser.add_argument("--kernel-sha256")
    prepare_parser.add_argument("--runner-sha256")
    sub.add_parser("run")
    publish_parser = sub.add_parser("publish")
    publish_parser.add_argument("--review", required=True)
    publish_parser.add_argument("--review-sha256", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> dict[str, Any]:
    args = parser().parse_args(argv)
    value = {"prepare": prepare, "run": run, "publish": publish}[args.action](args)
    print(json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":")))
    return value


if __name__ == "__main__":
    main()
